/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in

# =========================================================================*/ 
   
/** ************************************************************************ 
* @brief     stores the information of a 3D TNSE  and 6DOF system matrix 
* @author    Sashikumaar Ganesan  
* @date      13.5.2017
* @History    
 ************************************************************************  */
#include <string.h>
#include <stdlib.h>
#include <Database.h>
#include <System_6DOF.h>
#include <FEDatabase3D.h>
#include <FEVectFunct3D.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>
#include <HexaAffin.h>
#include <HexaIsoparametric.h>
#include <HexaTrilinear.h>
#include <LinAlg.h>


TSystem_6DOF::TSystem_6DOF(TFESpace3D *gridFESpace, TFEVectFunct3D *gridPosFEVect)
{
  GridPosFEVect = gridPosFEVect;
  gridpos = GridPosFEVect->GetValues(); 
  GridFESpace = gridFESpace;
  N_GridDOFs = GridFESpace->GetN_DegreesOfFreedom();
  N_GridActive = GridFESpace->GetActiveBound();  
     
  gridpos_CM = new double[3*N_GridDOFs]; 
//   RotGridPos = new double[3*N_GridDOFs];  
}

TSystem_6DOF::~TSystem_6DOF()
{
  
  
}



void TSystem_6DOF::Init6DOF( )
{
  int i;
  double MOI_Tensor[6]; // six independent components of MOI tensor
  double P_MOI[3];      //  principle moments of inertia
    
  char compz = 'V';
  char uplo = 'U'; 
  
  TCollection *Coll;
  
 
    memset(P_MOI, 0, 3*SizeOfDouble);
    memset(P_Axis, 0, 9*SizeOfDouble);

    GridPosFEVect->GridToData(); 

    memcpy(gridpos_CM, gridpos, 3*N_GridDOFs*SizeOfDouble);
    
    this->GetMOI(MOI_Tensor);

    /** compute the moment of Inertia, which are the eigen values*/
    /** eigenvector of each eigen value is the principal axis vector */
    FindEigenValues(MOI_Tensor, uplo, 3, compz, P_MOI, P_Axis);
    
    for(i=0;i<9;i++)
      cout <<"P_Axis["<<i<<"] "<<P_Axis[i]<<endl;
    for(i=0;i<3;i++)
      cout <<"P_MOI["<<i<<"] "<<P_MOI[i]<<endl;    
    
    Ixx = P_MOI[0]; 
    Iyy = P_MOI[1]; 
    Izz = P_MOI[2];
  
    // x,y,z of grid WRT CM
    Dadd(N_GridDOFs, -CGx_Body, gridpos_CM);              // Remove CMX from grid x
    Dadd(N_GridDOFs, -CGy_Body, gridpos_CM+N_GridDOFs);   // Remove CMY from grid y
    Dadd(N_GridDOFs, -CGz_Body, gridpos_CM+2*N_GridDOFs); // Remove CMZ from grid z   
    
    /** sol */
    memset(sol, 0, 12*SizeOfDouble);    
    memset(oldsol, 0, 12*SizeOfDouble);
    memset(rhs, 0, 12*SizeOfDouble);
    memset(Rot_Vec, 0, 3*SizeOfDouble);
    memset(TotalBodyForce, 0, 6*SizeOfDouble);
}


// solve the 6DOF system using Crank-Nicolson scheme
void TSystem_6DOF::SolveAndUpdate(double tau)
{
  int i, j, N_MaxItr=100;
  double eps =1.e-12, res1, res2, *newgridx, *newgridy, *newgridz;
  double sw1, sw2, sw3, cw1, cw2, cw3, *refgridposx, *refgridposy, *refgridposz;
  double newPaxes[9], t = TDatabase::TimeDB->CURRENTTIME;
  
  //body force 
  TotalBodyForce[0] = 0; TotalBodyForce[1] = 0; TotalBodyForce[2] = 0; // translation force
  TotalBodyForce[3] =  cos(Pi*t/2); TotalBodyForce[4] =  cos(Pi*t/2); TotalBodyForce[5] = cos(Pi*t/2); //rotational force
  
  cout<< "TotalBodyForce[5]  " <<TotalBodyForce[5] <<endl;
  
  //copy old sol to sol
  memcpy(oldsol, sol, 12*SizeOfDouble);
   
  //assemble rhs and add it to working rhs  
  memcpy(B, oldsol, 12*SizeOfDouble);
  AssembleRhs();
  Daxpy(12, 0.5*tau, rhs, B);
  
  for(i=0; i<N_MaxItr; i++) 
   {  
    res1=0;
    for(j=0;j<12;j++)
     res1 += sol[j]*sol[j];
        
    // assemble current rhs
    this->AssembleRhs();
   
    //update the sol
    for(j=0;j<12;j++)
      sol[j] = B[j] + 0.5*tau*rhs[j];
      
    //===================================================
    // check the convergence
    res2=0;
    for(j=0;j<12;j++)
     res2 += sol[j]*sol[j];
    
    if(res2!=0)
     {
      if( (fabs(res2 - res1)/res2)< eps )
       { break; }
     } // if(res2!=0)
    else
     { break; }     
    //===================================================
   } //for i
  
  if(i==N_MaxItr)
   {
    OutPut("6DOF system not converged!!!! N_MaxItr: "<< N_MaxItr <<endl);
   }
//   for(j=3;j<6;j++)
//     cout << i<<" Sol["<<j<<"] "<<sol[j]<<endl; 
  
   //project rotational vector onto std axis
   P_Axis_RotVect[0] = sol[3] - oldsol[3];
   P_Axis_RotVect[1] = sol[4] - oldsol[4];
   P_Axis_RotVect[2] = sol[5] - oldsol[5];
 
//    for(j=0;j<3;j++)
//      cout<<" Rot_Vec["<<j<<"] "<<P_Axis_RotVect[j]<<endl; 
     
   //project rotation vec into Paxis direction
   this->ProjectVect(P_Axis_RotVect); 

//    for(j=0;j<3;j++)
//      cout<<" Rot_Vec["<<j<<"] "<<P_Axis_RotVect[j]<<endl; 
   
   sw1 = sin(Rot_Vec[0]);
   sw2 = sin(Rot_Vec[1]);
   sw3 = sin(Rot_Vec[2]);
   cw1 = cos(Rot_Vec[0]);
   cw2 = cos(Rot_Vec[1]);
   cw3 = cos(Rot_Vec[2]);   
   
   //update the P_Axis 
   for(j=0;j<3;j++)  // rotation of x, y and z coordinates of grid
    {   
     newPaxes[3*j] = P_Axis[3*j]*cw3*cw2 - P_Axis[3*j+1]*sw3*cw2 + P_Axis[3*j+2]*sw2;
     newPaxes[3*j+1] = P_Axis[3*j]*(sw3*cw1 + sw1*sw2*cw3) + P_Axis[3*j+1]*(cw1*cw3 - sw3*sw2*sw1) - P_Axis[3*j+2]*sw1*cw2;
     newPaxes[3*j+2] = P_Axis[3*j]*(sw1*sw3 - cw1*sw2*cw3) + P_Axis[3*j+1]*(cw3*sw1 + cw1*sw2*sw3) + P_Axis[3*j+2]*cw2*cw1;         
    }   
   
   memcpy(P_Axis, newPaxes, 9*sizeof(double)); 
   
   // first rotate the grid      
   newgridx = gridpos;
   newgridy = gridpos+N_GridDOFs;
   newgridz = gridpos+2*N_GridDOFs;
   refgridposx = gridpos_CM;
   refgridposy = gridpos_CM+N_GridDOFs;   
   refgridposz = gridpos_CM+2*N_GridDOFs;
      
   for(j=0;j<N_GridDOFs;j++)  
    {           
     newgridx[j] = refgridposx[j]*cw3*cw2 - refgridposy[j]*sw3*cw2 + refgridposz[j]*sw2;
     newgridy[j] = refgridposx[j]*(sw3*cw1 + sw1*sw2*cw3) + refgridposy[j]*(cw1*cw3 - sw3*sw2*sw1) - refgridposz[j]*sw1*cw2;
     newgridz[j] = refgridposx[j]*(sw1*sw3 - cw1*sw2*cw3) + refgridposy[j]*(cw3*sw1 + cw1*sw2*sw3) + refgridposz[j]*cw2*cw1;         
    }
  
   //since P_Axis alos moved, keep the rotated grid as the reference grid
   memcpy(gridpos_CM, gridpos, 3*N_GridDOFs*sizeof(double));
  
   // now translate gridpos
   Dadd(N_GridDOFs, sol[0]+CGx_Body, gridpos);
   Dadd(N_GridDOFs, sol[1]+CGy_Body, gridpos+N_GridDOFs);    
   Dadd(N_GridDOFs, sol[2]+CGz_Body, gridpos+2*N_GridDOFs);   

   //finally update the grid
   GridPosFEVect->DataToGrid();  
}

//project rot vector into P_axis direction
void TSystem_6DOF::ProjectVect(double *P_Axis_RotVect)
{
  int i;
  double UPaxes[9], norm1, norm2, norm3;
  
  UPaxes[0] = P_Axis[0];
  UPaxes[1] = P_Axis[3];
  UPaxes[2] = P_Axis[6];
  UPaxes[3] = P_Axis[1];
  UPaxes[4] = P_Axis[4];
  UPaxes[5] = P_Axis[7];
  UPaxes[6] = P_Axis[2];
  UPaxes[7] = P_Axis[5];
  UPaxes[8] = P_Axis[8];
  
  norm1 = pow(UPaxes[0]*UPaxes[0] + UPaxes[1]*UPaxes[1] + UPaxes[2]*UPaxes[2],0.5);
  norm2 = pow(UPaxes[3]*UPaxes[3] + UPaxes[4]*UPaxes[4] + UPaxes[5]*UPaxes[5],0.5);
  norm3 = pow(UPaxes[6]*UPaxes[6] + UPaxes[7]*UPaxes[7] + UPaxes[8]*UPaxes[8],0.5);
  
  for(i=0;i<3;i++)
   {
    UPaxes[i] /=norm1;
    UPaxes[i+3] /=norm2;
    UPaxes[i+6] /=norm3;
   }
 
  Rot_Vec[0] = P_Axis_RotVect[0]*UPaxes[0] + P_Axis_RotVect[1]*UPaxes[1] + P_Axis_RotVect[2]*UPaxes[2];
  Rot_Vec[1] = P_Axis_RotVect[0]*UPaxes[3] + P_Axis_RotVect[1]*UPaxes[4] + P_Axis_RotVect[2]*UPaxes[5];
  Rot_Vec[2] = P_Axis_RotVect[0]*UPaxes[6] + P_Axis_RotVect[1]*UPaxes[7] + P_Axis_RotVect[2]*UPaxes[8];
}

//assemble rhs
void TSystem_6DOF::AssembleRhs()
{ 
  rhs[0] = sol[6];
  rhs[1] = sol[7];  
  rhs[2] = sol[8];
  rhs[3] = sol[9];
  rhs[4] = sol[10];
  rhs[5] = sol[11];
  
  rhs[6] = TotalBodyForce[0]/BodyMass;
  rhs[7] = TotalBodyForce[1]/BodyMass;  
  rhs[8] = TotalBodyForce[2]/BodyMass;
  rhs[9] = (TotalBodyForce[3] - (Izz - Iyy)*sol[10]*sol[11])/Ixx;
  rhs[10] = (TotalBodyForce[4] + (Izz - Ixx)*sol[11]*sol[9])/Iyy;
  rhs[11] = (TotalBodyForce[5] - (Iyy - Ixx)*sol[9]*sol[10])/Izz;  
}

// compute MOI
void TSystem_6DOF::GetMOI(double *MOI_Tensor)
{
  int i,j,k, N_Points;
  int N_Cells, N_Joints, MaxApproxOrder = 2;
  
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  JointType jointtype;
  boolean IsIsoparametric;
  RefTrans3D RefTrans;
  TRefTrans3D *rt;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf2;

  double *weights, *xi, *eta, *zeta, Ixy=0., Iyz=0., Izx=0.;
  double Izzcell, Iyycell, Ixxcell, Ixycell, Iyzcell, Izxcell;
  double absdetjk[MaxN_QuadPoints_3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
  double vol, locvol;
  
  Coll = GridFESpace->GetCollection();
  
  CGx_Body = 0; CGy_Body = 0; CGz_Body = 0;
  
  for(i=0;i<N_GridDOFs;i++)
   {
    CGx_Body +=gridpos[i];
    CGy_Body +=gridpos[i+N_GridDOFs];
    CGz_Body +=gridpos[i+2*N_GridDOFs];
   }
  
  CGx_Body /=(double)N_GridDOFs;
  CGy_Body /=(double)N_GridDOFs; 
  CGz_Body /=(double)N_GridDOFs;  // end of CG 
  
  cout <<"CGx is " <<CGx_Body <<endl;
  cout <<"CGy is " <<CGy_Body <<endl;
  cout <<"CGz is " <<CGz_Body <<endl;
    
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
   {
    cell = Coll->GetCell(i);
    N_Joints = cell->GetN_Joints();
    IsIsoparametric = FALSE;
    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      switch(cell->GetType())
      {
	case Tetrahedron:
	  RefTrans = TetraAffin;
	  QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(3);
	break;
	
	case Brick:
	  RefTrans = HexaAffin;
	  QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(2);
	break;

	case Hexahedron:
	  RefTrans = HexaTrilinear;
	  QuadFormula = TFEDatabase3D::GetQFHexaFromDegree(2);
	break;
	
        default:
           Error("Unknown cell->GetType" << endl);
           exit(-1);
      } // endswitch GetType
    } // endfor j

    if(IsIsoparametric)
    {
      switch(RefTrans)
      {
	case HexaAffin:
	case HexaTrilinear:
	RefTrans = HexaIsoparametric;
	break;

	case TetraAffin:
	RefTrans = TetraIsoparametric;
	break;
	
        default:
           Error("Unknown RefTrans" << endl);
           exit(-1);
      }                                                     // endswitch
    }                                                       // endif IsIsoparametric

    qf2 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf2->GetFormulaData(N_Points, weights, xi, eta, zeta);

    rt = TFEDatabase3D::GetRefTrans3D(RefTrans);
    switch(RefTrans)
    {
      case TetraAffin:
        ((TTetraAffin *)rt)->SetCell(cell);
        ((TTetraAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
//       case TetraIsoparametric:
//         ((TTetraIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
//         ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
//         ((TTetraIsoparametric *)rt)->SetCell(cell);
//         ((TTetraIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
//       break;
      case HexaAffin:
        ((THexaAffin *)rt)->SetCell(cell);
        ((THexaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
      case HexaTrilinear:
        ((THexaTrilinear *)rt)->SetCell(cell);
        ((THexaTrilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
      break;
//       case HexaIsoparametric:
//         ((THexaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
//         ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
//         ((THexaIsoparametric *)rt)->SetCell(cell);
//         ((THexaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta, X, Y, Z, absdetjk);
//       break;
      
      default:
        cout<<"Wrong transformation, Isoparameteric not yet implemented/tested!! "<<endl;
	exit(0);
      break;      
    } // endswitch
    
    locvol = 0., Izzcell = 0., Iyycell = 0., Ixxcell = 0.; Ixycell = 0.; Iyzcell = 0.; Izxcell = 0.;
    for(j=0;j<N_Points;j++)
     {
      locvol += weights[j]*absdetjk[j];
      Izzcell+= weights[j]*absdetjk[j]*(pow((X[j] - CGx_Body),2) + pow((Y[j] - CGy_Body),2));
      Iyycell+= weights[j]*absdetjk[j]*(pow((X[j] - CGx_Body),2) + pow((Z[j] - CGz_Body),2));
      Ixxcell+= weights[j]*absdetjk[j]*(pow((Y[j] - CGy_Body),2) + pow((Z[j] - CGz_Body),2));
      Ixycell+= weights[j]*absdetjk[j]*(X[j] - CGx_Body)*(Y[j] - CGy_Body);
      Iyzcell+= weights[j]*absdetjk[j]*(Y[j] - CGy_Body)*(Z[j] - CGz_Body);
      Izxcell+= weights[j]*absdetjk[j]*(Z[j] - CGz_Body)*(X[j] - CGx_Body);      
     }
   BodyMass +=locvol;
   Izz+= Izzcell;
   Iyy+= Iyycell;
   Ixx+= Ixxcell;
   Ixy+= Ixycell;
   Iyz+= Iyzcell;
   Izx+= Izxcell;   
  } // end of i < N_Cells

   //multiply with the density
  Izz *=TDatabase::ParamDB->P0;
  Iyy *=TDatabase::ParamDB->P0;
  Ixx *=TDatabase::ParamDB->P0;
  Ixy *=TDatabase::ParamDB->P0;
  Iyz *=TDatabase::ParamDB->P0;
  Izx *=TDatabase::ParamDB->P0;
  BodyMass*=TDatabase::ParamDB->P0;
  
//      cout<< "Izz is "<< Izz<<endl;
//      cout<< "Iyy is "<< Iyy<<endl;
//      cout<< "Ixx is "<< Ixx<<endl;
//      cout<< "Ixy is "<< Ixy<<endl;
//      cout<< "Iyz is "<< Iyz<<endl;
//      cout<< "Izx is "<< Izx<<endl;
//      cout<< "BodyMass is "<< BodyMass<<endl;

     MOI_Tensor[0] = Ixx; MOI_Tensor[1] = -Ixy; MOI_Tensor[2] = Iyy; //packed symmetric MOI_Tensor columnwise for LinAlg.C
     MOI_Tensor[3] = -Izx; MOI_Tensor[4] = -Iyz; MOI_Tensor[5] = Izz; 
} // end of GetMOI

 
