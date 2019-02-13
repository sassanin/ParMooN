// Navier-Stokes problem, 3D Channel flow
// 
// u(x,y) = unknown
// p(x,y) = unknown
#include<BoundFace.h>

void ExampleFile()
{
#ifdef _MPI
 int rank;
 MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
  { 
   OutPut("Example: Channel_6DOF_ALE.h " << endl);
  }

 TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION=1;
 TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;

 /** solve 6 dof only */
// #define _6DOFONLY_
 #define _NSEONLY_
}

void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0.0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0.0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0.0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0.0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// Description for moving grid
// ========================================================================
void GridU1(double x, double y,double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void GridU2(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void GridU3(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void GridBoundCondition(int BdComp, double x, double y, double z, BoundCond &cond)
{
	 switch(BdComp)
	  {
	   case 0: case 1: case 2: case 4: // side walls
//		    cond = SLIP_FRICTION_PENETRATION_RESISTANCE;

		    TDatabase::ParamDB->MESH_SLIP_WITH_FRICTION=0;
		    TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY=6; // free slip, no penetration, same for velo
		    cond = NEUMANN;
	    break;
	   case 3: // bottom
	   case 5: // top surface, traction bc
	     cond = DIRICHLET;
	   break;

	   default:
	     cout << "wrong boundary part number" << endl;
	     cond = DIRICHLET;
	    break;
	  }
}

// value of boundary condition
void GridBoundValue(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,double *z,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}
//=============================================================================


void ModifyCoords(double x, double y, double z, 
		  double &X, double &Y, double &Z, double t)
{

 X = x;
 Y = y;
 Z = z;
  }
  
void MoveBoundWithVelo(TFEVectFunct3D *Velocity, TFEVectFunct3D *GridPos, double *GridOld, int N_MovVert, TVertex **MovBoundVert,
		       int N_Movfaces, int * Movfaces, TBaseCell ** MovCells, double dt)
{
 int i, j, k, l, Cell_No, Joint_No, MaxLen, pt, BdId, N_LocalDOFs, N_Vertices; 
 int *VeloBeginIndex, *VeloGlobalNumbers, *DOF; 
 int *GridBeginIndex, *GridGlobalNumbers, N_Inner, N_W;
 const int *TmpFV, *TmpLen;

 double *ValuesVX, *ValuesVY, *ValuesVZ, s, t, u, T[4]={0,0,0,0}, S[4]={0,0,0,0};
 double *GridX, *GridY, *GridZ, *OldGridX, *OldGridY, *OldGridZ;
 double xi[8], eta[8], zeta[8],  X[8], Y[8], Z[8], VX[8], VY[8], VZ[8];
 double FunctValues[8][MaxN_BaseFunctions3D];   
   
 TBaseCell* cell;
 TBoundFace* Bdface;
 TBoundComp3D* Bdcomp;
 TVertex* vert; 
 TFESpace3D *VelocitySpace, *GridSpace;
 FE3D FEId;
 TFE3D *Element;
 BaseFunct3D BF;
 TBaseFunct3D *bf;
  
 bool UPDATE; 
  
  VelocitySpace = Velocity->GetFESpace3D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  ValuesVZ = ValuesVY + Velocity->GetLength(); 
 
  GridSpace = GridPos->GetFESpace3D();
  N_Inner = GridSpace->GetN_Inner();
  N_W = GridPos->GetLength();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridX = GridPos->GetValues();
  GridY = GridX + N_W;
  GridZ = GridY + N_W;

  OldGridX = GridOld;
  OldGridY = OldGridX + N_W;
  OldGridZ = OldGridY + N_W;
  
  //indicator = new int[N_W];
//  memset(indicator, 0, N_W*SizeOfInt);

   // sent the Clipbord/marking
   for(i=0;i<N_MovVert;i++)
     MovBoundVert[i]->SetClipBoard(-5);
   
    for(i=0;i<N_Movfaces;i++)
     {
      cell = MovCells[i];
      Cell_No = cell->GetGlobalCellNo();
      Joint_No = Movfaces[i];
      Bdface = (TBoundFace*)cell->GetJoint(Joint_No);
      Bdcomp = Bdface->GetBoundComp();
      BdId = Bdcomp->GetID();
      
      //move the top boundary at the last
//       if(BdId==5) continue; do the sorting of Moving vert accordingly, so that BdId=5 vert comes at the last
      
      cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);     
      pt = TmpLen[Joint_No];
 
//      UPDATE=FALSE;
//      for(int ii=0;ii<pt;ii++)
//       {
//        vert = cell->GetVertex(TmpFV[ Joint_No*MaxLen + ii]);
//        if(vert->GetClipBoard() == -5)
//         {
//          UPDATE=TRUE;
//          break;
//         }
//       }

//     if(UPDATE==FALSE)  continue;

     //compute the velocity at the vertices
     FEId = VelocitySpace->GetFE3D(Cell_No, cell);
     Element = TFEDatabase3D::GetFE3D(FEId);
     bf = Element->GetBaseFunct3D();
     N_LocalDOFs = Element->GetN_DOF();    
      
      switch(pt)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0; xi[3]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1; eta[3] = 0;
          zeta[0]= 0; zeta[1]= 0; zeta[2]= 0; zeta[3]= 1;
  
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          VZ[0] = VZ[1] = VZ[2] = VZ[3] = 0;
          N_Vertices = 4;
          //   cout << "Tetra cell " <<endl;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1; xi[4]  = -1; xi[5]  =  1; xi[6]  = -1; xi[7]  = 1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1; eta[4] = -1; eta[5] = -1; eta[6] =  1; eta[7] = 1;
          zeta[0]= -1; zeta[1]= -1; zeta[2]= -1; zeta[3]= -1; zeta[4]=  1; zeta[5]=  1; zeta[6]=  1; zeta[7]= 1;
	  
          VX[0] = VX[1] = VX[2] = VX[3] = 0; VX[4] = VX[5] = VX[6] = VX[7] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0; VY[4] = VY[5] = VY[6] = VY[7] = 0;
          VZ[0] = VZ[1] = VZ[2] = VZ[3] = 0; VZ[4] = VZ[5] = VZ[6] = VZ[7] = 0;
          N_Vertices = 8;
        break;

        default:
          Error("only triangles and quadrilateral surfaces are allowed!" << endl);
          exit(-1);
      } // endswitch
      
      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D000, xi[j], eta[j], zeta[j], FunctValues[j]);      

      DOF = VeloGlobalNumbers + VeloBeginIndex[Cell_No];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        u = ValuesVZ[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
          VZ[l] += FunctValues[l][j]*u;  
        } // endfor l
       } // endfor j  
          
      FEId = GridSpace->GetFE3D(Cell_No, cell);
      Element = TFEDatabase3D::GetFE3D(FEId);
      BF = Element->GetBaseFunct3D_ID();
      if( (BF != BF_C_T_P1_3D) && (BF != BF_C_H_Q1_3D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct3D();
      N_LocalDOFs = Element->GetN_DOF();      
      
      DOF = GridGlobalNumbers + GridBeginIndex[Cell_No];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
         {
//          // set the clipboard
//          //vert = cell->GetVertex(TmpFV[Joint_No*MaxLen + j]);
//          if(indicator[l] == 0)
//            { indicator[l] == 1; } // bound vertex already updated
//          else
//            { continue;}

          switch(BdId)
	      {
//            case 0: case 1: // wall on z-plane
////             GridX[l] = OldGridX[l] + dt*VX[j];
////             GridY[l] = OldGridY[l] + dt*VY[j];
//             GridX[l] = OldGridX[l] + 0*VX[j];
//             GridY[l] = OldGridY[l] + 0*VY[j];
//
//            break;
//            case 2: case 4: // wall on x-plane
////             GridY[l] = OldGridY[l] + dt*VY[j];
////             GridZ[l] = OldGridZ[l] + dt*VZ[j];
//             GridY[l] = OldGridY[l] + 0*VY[j];
//             GridZ[l] = OldGridZ[l] + 0*VZ[j];
//
//            break;
	        case 5: // free surface
             GridX[l] = OldGridX[l] + dt*VX[j];
             GridY[l] = OldGridY[l] + dt*VY[j];
             GridZ[l] = OldGridZ[l] + dt*VZ[j];

//	           GridX[l] = OldGridX[l] +0.;
//	           GridY[l] = OldGridY[l] + 0.25;
//	           GridZ[l] = OldGridZ[l] + 0.;

	         //    cout << "X " <<  GridX[l] - OldGridX[l]<< " Y " <<  GridY[l] - OldGridY[l] << " Z " << GridZ[l]  - OldGridZ[l]<< endl;

            break;
	    
	        default:
	         cout << "wrong boundary part number, Check Example file" << endl;
	        break;
	      }// endswitch
         } //if(k
        } // endfor j
     }
        
    for(i=0;i<N_Movfaces;i++)
     {
      cell = MovCells[i];
      Joint_No = Movfaces[i];
      cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      Bdface = (TBoundFace*)cell->GetJoint(Joint_No);
      Bdcomp = Bdface->GetBoundComp();

      pt = TmpLen[Joint_No];
      for(int ii=0;ii<pt;ii++)
      {
	   vert = cell->GetVertex(TmpFV[ Joint_No*MaxLen + ii]);
	   Bdcomp->GetTSofXYZ(vert->GetX(), vert->GetY(), vert->GetZ(), T[ii+1], S[ii+1]);
      }

      Bdface->SetParameters(T,S);
     }

// delete [] indicator;
//  cout << " MoveBoundWithVelo" <<endl;
//  exit(0);
}
  
void ModifyBdCoords(int N_MovVert, TVertex **MovBoundVert,
		      int N_Movfaces, int * Movfaces, 
		      TBaseCell ** MovCells, double t)
{
 int i, j, k, l, m;

 double x, y, z, Fact=0.06;
 double disp, theta, T[4]={0,0,0,0}, S[4]={0,0,0,0};

 const int *TmpFV, *TmpLen;
 int MaxLen,pt,jid;
 TBaseCell* cell;
 TBoundFace* Bdface;
 TBoundComp3D* Bdcomp;
 TVertex* vert;

//
//  disp = 0.05*sin(3.93*t);
//
   disp =0.5;

   for(i=0;i<N_MovVert;i++)
    {
     MovBoundVert[i]->GetCoords(x, y, z);
     y += disp;
     MovBoundVert[i]->SetCoords(x, y, z);
    }


 //   cout << " disp " <<disp << endl;

    for(i=0;i<N_Movfaces;i++)
     {
      cell = MovCells[i];
      jid = Movfaces[i];
      cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      Bdface = (TBoundFace*)cell->GetJoint(jid);
      Bdcomp = Bdface->GetBoundComp();

      pt = TmpLen[jid];
      for(int ii=0;ii<pt;ii++)
      {
	   vert = cell->GetVertex(TmpFV[ jid*MaxLen + ii]);
	   Bdcomp->GetTSofXYZ(vert->GetX(), vert->GetY(), vert->GetZ(), T[ii+1], S[ii+1]);
      }

      Bdface->SetParameters(T,S);
     }
     
}// ModifyBdCoords



void  GetMovingBoundData(TCollection *coll, int &N_MovVert,
			             TVertex ** &MovBoundVert, int &N_Movfaces,
			             int * &Movfaces, TBaseCell ** &MovCells)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int temp;
 const int *TmpFV, *TmpLen;
 int MaxLen, N_Points;
 int bdid,pt;
 
 double  x, y, z, x1, y1, z1;
  
 TBaseCell *Me, *temp_cell;
 TJoint *Joint;
 TBoundComp *BoundComp;  
 TBoundFace* Bdface;
 TVertex *temp_vert, **Vertices;

  N_Cells = coll->GetN_Cells();
  N_MovVert = 0;
  N_Movfaces = 0;     
  
    for(j=0;j<N_Cells;j++)
     {
      Me = coll->GetCell(j);
      Me->SetGlobalCellNo(j);
      
      k = Me->GetN_Vertices();
      for(l=0;l<k;l++)
       {
        temp_vert = Me->GetVertex(l);
        temp_vert->SetClipBoard(-1);
       }
     }
  
     for(j=0;j<N_Cells;j++)
      {
       Me = coll->GetCell(j);
       Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);

       k = Me->GetN_Joints();
       
       for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
         if(Joint->GetType() == BoundaryFace)
          {
           Bdface = (TBoundFace*)Joint;
           BoundComp = Bdface->GetBoundComp();
           bdid = BoundComp->GetID();
           if(bdid ==5)
            {
             pt = TmpLen[l];
             N_Movfaces ++; // no need of clips as they belong to atmost 1 cell
	     for(int ii =0;ii<pt;ii++)
	      {
	       temp_vert = Me->GetVertex(TmpFV[ l*MaxLen + ii]);
	       if(temp_vert->GetClipBoard()==-1)
	        {
                 temp_vert->SetClipBoard(1);
	    	  N_MovVert++;
	    	 }
	    	}
	        }
	      }
	  }// endfor l
    }
      
      MovBoundVert = new TVertex*[N_MovVert];
      Movfaces = new int[N_Movfaces];
      MovCells = new TBaseCell*[N_Movfaces];
      
//       cout << "N_MovVert " << N_MovVert  << " N_Movfaces " << N_Movfaces << endl;
      
      N_MovVert = 0;
      N_Movfaces = 0;
      
      for(j=0;j<N_Cells;j++)
      {
       Me = coll->GetCell(j);
       Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        k = Me->GetN_Joints();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
         if(Joint->GetType() == BoundaryFace)
          {
           Bdface = (TBoundFace*)Joint;
           BoundComp = Bdface->GetBoundComp();
           bdid = BoundComp->GetID();
           if(bdid ==5)
           {
        	pt = TmpLen[l];
        	Movfaces[N_Movfaces] = l;
        	if(l>4)
        	{
                  cout << "garbage " << endl;
        	}
	      
        	MovCells[N_Movfaces] = Me;
	   
	    for(int ii =0;ii<pt;ii++)
	    {
	      temp_vert = Me->GetVertex(TmpFV[ l*MaxLen + ii]);
	      if(temp_vert->GetClipBoard()==1)
	      {
		temp_vert->SetClipBoard(-1);
		MovBoundVert[N_MovVert] = temp_vert;
		N_MovVert++; 
	      }	      
	    }	    
	    N_Movfaces++;
	   }
	  }
	 }// endfor l
    }
      
//   cout << "N_MovVert" << N_MovVert  << " N_Movfaces " << N_Movfaces << endl;
       
}// GetMovingBoundData



// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{ 

 switch(CompID)
  {
   case 0: case 1: case 2: case 4: // side walls
      cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
      TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION=1;
	  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY=6; // free slip, no penetration, same for velo
      break;
   case 3: // bottom
     cond = DIRICHLET;
   break;

   case 5: // top surface, traction bc
     cond = NEUMANN;
   break;

   default:
     cout << "wrong boundary part number" << endl;
     cond = DIRICHLET;
    break;
  } 
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{

   switch(CompID)
   {         
    case 0: case 1: case 2:  case 3:  case 4: 
         value=0;
         break;
    case 5: 
         value=0;
         break;
    default: 
      cout << "wrong boundary part number" << endl;
      value=0;
            break;
    }  
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
   switch(CompID)
   {         
    case 0: case 1: case 2:  case 3:  case 4: 
         value=0;
         break;
    case 5: 
         value=0;
         break;
    default: 
      cout << "wrong boundary part number" << endl;
      value=0;
            break;
    }  
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
   switch(CompID)
   {         
    case 0: case 1: case 2:  case 3:  case 4: 
         value=0;
         break;
    case 5: 
         value=0;
         break;
    default: 
      cout << "wrong boundary part number" << endl;
      value=0;
            break;
    }  
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = -1.; // f2
    coeff[3] = 0; // f3
  }
}
