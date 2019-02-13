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
   
#include <Constants.h>
#include <FreeSurface3D.h>
#include <MGLevel3D.h>
#include <Database.h>
#include <ItMethod.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <JacobiIte.h>
#include <MultiGridScaIte.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <FEDatabase3D.h>
#include <BoundFace.h>
#include <IsoJointEqN.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <HexaIsoparametric.h>
#include <TetraAffin.h>
#include <TetraIsoparametric.h>

#include <FEFunction3D.h>

#include <InterfaceJoint3D.h>

#include <NodalFunctional3D.h>

#include <stdlib.h>
#include <string.h>

int CheckNormal(double n1, double n2, double n3, TBaseCell *Cell, int JointNr);
int IsJointDOF(int DOF, int N_JointDOF, int *JointDOF);

// ========================================================================
// calculate all parameters which are needed for calculation
// ========================================================================
void CalculateAllParameters()
{
  const double Lhalf = 1.796755984723713041136085;
  int Law = TDatabase::ParamDB->FS_MAGNETLAW;
  double eta = TDatabase::ParamDB->FS_ETA;
  double rho = TDatabase::ParamDB->FS_RHO;
  double alpha = TDatabase::ParamDB->FS_ALPHA;
  double g = TDatabase::ParamDB->FS_G;
  double Ms = TDatabase::ParamDB->FS_MS;
  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double T = TDatabase::ParamDB->FS_T;
  double HM = TDatabase::ParamDB->FS_HM;
  double DELTA_H = TDatabase::ParamDB->FS_DELTA_H;
  double HT, H0, L, U;

  L = 2*Pi/3 * sqrt(alpha/(rho*g));
  TDatabase::ParamDB->FS_L = L;
  U = L / T;
  TDatabase::ParamDB->FS_U = U;

  if(fabs(Ms) < 1e-8)
  {
    chi0 = 0;
    TDatabase::ParamDB->FS_CHI0 = chi0;
  }

  if(fabs(chi0) < 1e-8)
  {
    Ms = 0;
    TDatabase::ParamDB->FS_MS = Ms;
  }

  if(fabs(HM) > 1e-8)
  {
    H0 = HM;
  }
  else
    if(fabs(DELTA_H) > 1e-8)
    {
      H0 = DELTA_H;
    }
    else
    {
      H0 = 0;
    }

  TDatabase::ParamDB->FS_H0 = H0;

  if(fabs(H0) < 1e-8 || fabs(Ms) < 1e-8)
  {
    // no outer field or non-magnetic medium
    TDatabase::ParamDB->FS_GAMMA = 0;
    TDatabase::ParamDB->FS_HT = 1;
  }
  else
  {
    TDatabase::ParamDB->FS_GAMMA = 3*chi0*H0/Ms;
    TDatabase::ParamDB->FS_HT = Ms*Lhalf/(3*chi0);
  }

  TDatabase::ParamDB->RE_NR = U * L * rho/eta;
  TDatabase::ParamDB->FS_WE = U*U * L * rho / alpha;

  OutPut("Parameters for free surface calculation" << endl);
  OutPut("FS_MAGNETLAW: " << TDatabase::ParamDB->FS_MAGNETLAW << " ");
  switch(Law)
  {
    case 0:
      OutPut("(Langevin)" << endl);
    break;
    case 1:
      OutPut("(Vislovich)" << endl);
    break;
  } // end switch Law

  OutPut("FS_L: " << TDatabase::ParamDB->FS_L << endl);
  OutPut("FS_U: " << TDatabase::ParamDB->FS_U << endl);
  OutPut("FS_T: " << TDatabase::ParamDB->FS_T << endl);

  OutPut("FS_ETA: " << TDatabase::ParamDB->FS_ETA << endl);
  OutPut("FS_RHO: " << TDatabase::ParamDB->FS_RHO << endl);
  OutPut("FS_ALPHA: " << TDatabase::ParamDB->FS_ALPHA << endl);
  OutPut("FS_G: " << TDatabase::ParamDB->FS_G << endl);

  OutPut("FS_MS: " << TDatabase::ParamDB->FS_MS << endl);
  OutPut("FS_CHI0: " << TDatabase::ParamDB->FS_CHI0 << endl);
  OutPut("FS_HM: " << TDatabase::ParamDB->FS_HM << endl);
  OutPut("FS_DELTA_H: " << TDatabase::ParamDB->FS_DELTA_H << endl);
  OutPut("FS_F: " << TDatabase::ParamDB->FS_F << endl);

  OutPut("FS_LH: " << TDatabase::ParamDB->FS_LH << endl);
  OutPut("FS_GAMMA: " << TDatabase::ParamDB->FS_GAMMA << endl);
  OutPut("FS_HT: " << TDatabase::ParamDB->FS_HT << endl);

  OutPut("FS_WRITE: " << TDatabase::ParamDB->FS_WRITE << endl);
  OutPut("FS_READ: " << TDatabase::ParamDB->FS_READ << endl);
  
  OutPut("FS_INNAME: " << TDatabase::ParamDB->FS_INNAME << endl);
  OutPut("FS_OUTNAME: " << TDatabase::ParamDB->FS_OUTNAME << endl);

  OutPut("FS_WE: " << TDatabase::ParamDB->FS_WE << endl);

  OutPut("****************************************" << endl);
  OutPut("calculated RE_NR: " << TDatabase::ParamDB->RE_NR << endl);
  OutPut("****************************************" << endl);
}

// ========================================================================
// calculate all field parameters
// ========================================================================
void CalculateFields()
{
  const double Lhalf = 1.796755984723713041136085;
  int Law = TDatabase::ParamDB->FS_MAGNETLAW;
  double Ms = TDatabase::ParamDB->FS_MS;
  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double T = TDatabase::ParamDB->FS_T;
  double HM = TDatabase::ParamDB->FS_HM;
  double H0 = TDatabase::ParamDB->FS_H0;
  double DELTA_H = TDatabase::ParamDB->FS_DELTA_H;
  double F = TDatabase::ParamDB->FS_F;
  double HT, H, time;

  double x, t;
  int i;

  time = TDatabase::TimeDB->CURRENTTIME;

  if(fabs(H0) > 1e-8)
  {
    H = HM + DELTA_H * cos( 2*Pi*F * T*time );
    TDatabase::ParamDB->FS_H2 = H/H0;

    switch(Law)
    {
      case 0:
        // Langevin function as magnetisation law
        if(fabs(Ms) > 1e-8)
        {
          t = (Ms+Lhalf*Ms/(3*chi0)-H)*3*chi0/(2*Ms);
          x = -t + sqrt(t*t + H*Lhalf*3*chi0/Ms);

          for(i=0;i<10;i++)
          {
            x -= (Ms*(1/tanh(x) - 1/x) + Ms*x/(3*chi0) - H ) /
                  (Ms*( 1/(x*x) - 1/(sinh(x)*sinh(x)) ) + Ms/(3*chi0));
          }
          TDatabase::ParamDB->FS_H1 = x*Ms/(3*chi0*H0);
        }
        else
        {
          TDatabase::ParamDB->FS_H1 = TDatabase::ParamDB->FS_H2;
        }
      break;

      case 1:
        // Vislovich approximation as magnetisation law
        if(fabs(Ms) > 1e-8)
        {
          HT = TDatabase::ParamDB->FS_HT;
          x =  -0.5*(Ms+HT-H);
          TDatabase::ParamDB->FS_H1 = (x + sqrt(x*x + H*HT)) / H0;
        }
        else
        {
          TDatabase::ParamDB->FS_H1 = TDatabase::ParamDB->FS_H2;
        }
      break;

      default:
        Error("Unknown magnetisation law chosen!" << endl);
        Error("Program terminated!" << endl);
        exit(-1);
    } // end switch Law
  }
  else
  {
    H = 0;
    TDatabase::ParamDB->FS_H1 = 0;
    TDatabase::ParamDB->FS_H2 = 0;
  }

  OutPut(TDatabase::TimeDB->CURRENTTIME);
  OutPut(" H0: " << H);
  OutPut(" H1: " << TDatabase::ParamDB->FS_H1);
  OutPut(" H2: " << TDatabase::ParamDB->FS_H2 << endl);
}

// ========================================================================
// declaration for grid handling
// ========================================================================
void GridBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// ========================================================================
// auxiliary routines
// ========================================================================
int CompareVerticesXYZ(TVertex *v0, TVertex *v1)
{
  double x0, y0, z0, x1, y1, z1;
  int ret = 1;
  
  v0->GetCoords(x0, y0, z0);
  v1->GetCoords(x1, y1, z1);

  if(fabs(x0-x1) < 1e-8)
  {
    if(fabs(y0-y1) < 1e-8)
    {
      ret = z0<z1;
    }
    else
    {
      ret = y0<y1;
    }
  }
  else
  {
    ret = x0<x1;
  }

  return ret;
} // end CompareVerticesXYZ

void SortVerticesXYZ(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr=new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(CompareVerticesXYZ(Array[i], Mid)) i++;

        while(CompareVerticesXYZ(Mid, Array[j])) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete rr;
} // end SortVerticesXYZ

void SortVertices(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr=new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete rr;
} // end SortVertices

// calculate the volume which is ocupied by all cell in Coll
double GetVolume(TCollection *Coll)
{
  int i,j,k;
  int N_Cells, N_Joints;
  TBaseCell *cell;
  TJoint *joint;
  JointType jointtype;
  bool IsIsoparametric;
  RefTrans3D RefTrans;
  TRefTrans3D *rt;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf2;
  int N_Points;
  double *weights, *xi, *eta, *zeta;
  double absdetjk[MaxN_QuadPoints_3D];
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double vol, locvol;
  int MaxApproxOrder = 2;

  vol = 0;
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
          QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(2);
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
	    
	  break;
      } // endswitch GetType

      jointtype = joint->GetType();
      if(jointtype == BoundaryFace)
      {
        if( ((TBoundFace*)joint)->GetBoundComp()->GetType() != Plane)
          IsIsoparametric = TRUE;
      }

      if(jointtype == InterfaceJoint3D)
      {
        if( ((TInterfaceJoint3D*)joint)->GetBoundComp()->GetType() != Plane)
          IsIsoparametric = TRUE;
      }

      if(jointtype == IsoInterfaceJoint3D)
        IsIsoparametric = TRUE;

      if(jointtype == IsoBoundFace)
        IsIsoparametric = TRUE;

      if(jointtype == IsoJointEqN)
        IsIsoparametric = TRUE;
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
	    
	  break;
      } // endswitch
    } // endif IsIsoparametric

    qf2 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf2->GetFormulaData(N_Points, weights, xi, eta, zeta);

    rt = TFEDatabase3D::GetRefTrans3D(RefTrans);
    switch(RefTrans)
    {
      case TetraAffin:
        // cout << "TetraAffin" << endl;
        ((TTetraAffin *)rt)->SetCell(cell);
        ((TTetraAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                           X, Y, Z, absdetjk);
      break;
      case TetraIsoparametric:
        // cout << "TetraIsoparametric" << endl;
        ((TTetraIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
        ((TTetraIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTetraIsoparametric *)rt)->SetCell(cell);
        ((TTetraIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                           X, Y, Z, absdetjk);
      break;
      case HexaAffin:
        // cout << "HexaAffin" << endl;
        ((THexaAffin *)rt)->SetCell(cell);
        ((THexaAffin *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                           X, Y, Z, absdetjk);
      break;
      case HexaTrilinear:
        // cout << "HexaTrilinear" << endl;
        ((THexaTrilinear *)rt)->SetCell(cell);
        ((THexaTrilinear *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                           X, Y, Z, absdetjk);
      break;
      case HexaIsoparametric:
        // cout << "HexaIsoparametric" << endl;
        ((THexaIsoparametric *)rt)->SetApproximationOrder(MaxApproxOrder);
        ((THexaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((THexaIsoparametric *)rt)->SetCell(cell);
        ((THexaIsoparametric *)rt)->GetOrigFromRef(N_Points, xi, eta, zeta,
                                           X, Y, Z, absdetjk);
      break;
    } // endswitch

    locvol = 0;
    for(j=0;j<N_Points;j++)
      locvol += weights[j]*absdetjk[j];

    vol += locvol;
  } // endfor i

  return vol;
} // end GetVolume

// ========================================================================
// find normal and tangential vectors for slip d.o.f.
// ========================================================================
void FindVectorsForSlipDOF(TFESpace3D *fespace,
        int &N_FaceDOF, int &N_EdgeDOF,
        int* &FaceDOF, int* &EdgeDOF,
        double* &FaceVectors, double* &EdgeVectors)
{
  int i,j,k;
  int N_Cells, N_Faces, N_SlipDOF;
  TCollection *coll;
  TBaseCell *cell;
  TJoint *joint;
  JointType jointtype;
  TBoundFace *bdface;
  int comp, comp1, comp2, comp3;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, *JointDOF, **JointDOFs, N_JointDOF, N, M;
  int *Bounds, SlipStart, SlipEnd, N_Slip;
  int *Flags;
  FE3D feID;
  TFEDesc3D *fedesc;
  TShapeDesc *ShapeDesc;
  const int *TmpFaceVertex, *TmpLen;
  int MaxLen;
  double x1, y1, z1, x2, y2, z2;
  double t11, t12, t13, t21, t22, t23;
  double n1, n2, n3, len, N1, N2, N3;

  Bounds = fespace->GetBoundaryNodesBound();
  SlipStart = Bounds[SLIP-2];
  SlipEnd = Bounds[SLIP-1];
  N_Slip = SlipEnd - SlipStart;
  Flags = new int[N_Slip];
  for(i=0;i<N_Slip;i++)
    Flags[i] = -1;

  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    feID = fespace->GetFE3D(i, cell);
    fedesc = TFEDatabase3D::GetFEDesc3DFromFE3D(feID);
    JointDOFs = fedesc->GetJointDOF();
    N_JointDOF = fedesc->GetN_JointDOF();

    DOF = GlobalNumbers + BeginIndex[i];
    
    N_Faces = cell->GetN_Joints();
    for(j=0;j<N_Faces;j++)
    {
      joint = cell->GetJoint(j);
      jointtype = joint->GetType();
      if(jointtype == BoundaryFace)
      {
        bdface = (TBoundFace *)joint;
        comp = bdface->GetBoundComp()->GetID();
        JointDOF = JointDOFs[j];
        for(k=0;k<N_JointDOF;k++)
        {
          N = DOF[JointDOF[k]]-SlipStart;
          if( (0 <= N) && (N < N_Slip) )
          {
            if(Flags[N] == -1)
            {
              // dof not handled yet, setting Flags to comp
              Flags[N] = comp;
            }
            else
            {
              if(Flags[N] != comp)
                Flags[N] = -2;
            } // endif == -1
          } // endif
        } // endfor k
      } // endif jointtype
    } // endfor j
  } // endfor i

  // count number of d.o.f. associated with edges and faces
  // set Flags to local number (for faces) and
  // to local number +N_Slip (for edges)
  N_FaceDOF = 0;
  N_EdgeDOF = 0;
  for(i=0;i<N_Slip;i++)
  {
    if(Flags[i] == -2)
    {
      Flags[i] = N_EdgeDOF + N_Slip;
      N_EdgeDOF++;
    }
    else
    {
      Flags[i] = N_FaceDOF;
      N_FaceDOF++;
    }
  }

  if(N_EdgeDOF+N_FaceDOF != N_Slip)
  {
    Error("Wrong numbers for d.o.f. at line " << __LINE__);
    Error(" of file " << __FILE__ << endl);
    exit(-1);
  }

  // for each dof store three vectors with three components
  FaceVectors = new double[9*N_FaceDOF];
  EdgeVectors = new double[9*N_EdgeDOF];
  memset(FaceVectors, 0, 9*N_FaceDOF*SizeOfDouble);
  memset(EdgeVectors, 0, 9*N_EdgeDOF*SizeOfDouble);

  FaceDOF = new int[N_FaceDOF];
  EdgeDOF = new int[N_EdgeDOF];
  memset(FaceDOF, 0, N_FaceDOF*SizeOfInt);
  memset(EdgeDOF, 0, N_EdgeDOF*SizeOfInt);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    feID = fespace->GetFE3D(i, cell);
    fedesc = TFEDatabase3D::GetFEDesc3DFromFE3D(feID);
    JointDOFs = fedesc->GetJointDOF();
    N_JointDOF = fedesc->GetN_JointDOF();

    DOF = GlobalNumbers + BeginIndex[i];
    
    N_Faces = cell->GetN_Joints();
    for(j=0;j<N_Faces;j++)
    {
      joint = cell->GetJoint(j);
      jointtype = joint->GetType();
      if( jointtype == BoundaryFace )
      {
        bdface = (TBoundFace *)joint;
        comp = bdface->GetBoundComp()->GetID();
        JointDOF = JointDOFs[j];

        ShapeDesc = cell->GetShapeDesc();
        ShapeDesc->GetFaceVertex(TmpFaceVertex, TmpLen, MaxLen);
        for(k=0;k<N_JointDOF;k++)
        {
          cell->GetVertex(TmpFaceVertex[j*MaxLen+0])->GetCoords(x1, y1, z1);
        
          cell->GetVertex(TmpFaceVertex[j*MaxLen+1])->GetCoords(x2, y2, z2);
          t11 = x2-x1; t12 = y2-y1; t13 = z2-z1;
          len = sqrt(t11*t11 + t12*t12 + t13*t13);
          t11 /= len; t12 /= len; t13 /= len;
        
          cell->GetVertex(TmpFaceVertex[j*MaxLen+(TmpLen[j]-1)])->GetCoords(x2, y2, z2);
          t21 = x2-x1; t22 = y2-y1; t23 = z2-z1;
          len = sqrt(t21*t21 + t22*t22 + t23*t23);
          t21 /= len; t22 /= len; t23 /= len;

          N1 = t12*t23 - t13*t22;
          N2 = t13*t21 - t11*t23;
          N3 = t11*t22 - t12*t21;
          len = sqrt(N1*N1 + N2*N2 + N3*N3);
          N1 /= len; N2 /= len; N3 /= len;

          N = DOF[JointDOF[k]]-SlipStart;
          if( (0 <= N) && (N < N_Slip) )
          {
            if(Flags[N] < N_Slip)
            {
              // face d.o.f.
              M = 9*Flags[N];
              FaceDOF[ Flags[N] ] = N;
              if(fabs(FaceVectors[M+0]) + fabs(FaceVectors[M+1]) +
                  fabs(FaceVectors[M+2]) < 1e-8)
              {
                // not handled yet
                FaceVectors[M+0] = t11;
                FaceVectors[M+1] = t12;
                FaceVectors[M+2] = t13;

                FaceVectors[M+3] = t21;
                FaceVectors[M+4] = t22;
                FaceVectors[M+5] = t23;

                FaceVectors[M+6] = N1;
                FaceVectors[M+7] = N2;
                FaceVectors[M+8] = N3;
              } // endif < 1e-8
            } // endif < N_Slip
            else
            {
              // edge d.o.f.
              M = 9 * (Flags[N]-N_Slip);
              EdgeDOF[ Flags[N]-N_Slip ] = N;
              if(fabs(EdgeVectors[M+0]) + fabs(EdgeVectors[M+1]) +
                  fabs(EdgeVectors[M+2]) < 1e-8)
              {
                // not handled yet
                EdgeVectors[M+0] = N1;
                EdgeVectors[M+1] = N2;
                EdgeVectors[M+2] = N3;
              }
              else
              {
                t21 = EdgeVectors[M+0];
                t22 = EdgeVectors[M+1];
                t23 = EdgeVectors[M+2];
                t11 = N1; t12 = N2; t13 = N3;
                if( fabs(t12*t23 - t13*t22) + fabs(t13*t21 - t11*t23)
                   +fabs(t11*t22 - t12*t21) > 1e-8)
                {
                  EdgeVectors[M+3] = N1;
                  EdgeVectors[M+4] = N2;
                  EdgeVectors[M+5] = N3;

                  n1 = -t12*t23 + t13*t22;
                  n2 = -t13*t21 + t11*t23;
                  n3 = -t11*t22 + t12*t21;
                  len = sqrt(n1*n1 + n2*n2 + n3*n3);
                  n1 /= len;
                  n2 /= len;
                  n3 /= len;

                  EdgeVectors[M+6] = n1;
                  EdgeVectors[M+7] = n2;
                  EdgeVectors[M+8] = n3;
                } // endif > 1e-8
              } // endif < 1e-8
            } // >= N_Slip
          } // endif
        } // endfor k
      } // endif BoundaryFace
    } // endfor j
  } // endfor i

  for(i=0;i<N_FaceDOF;i++)
    FaceDOF[i] += SlipStart;

  for(i=0;i<N_EdgeDOF;i++)
    EdgeDOF[i] += SlipStart;

/*
  cout << "FaceDOFs" << endl;
  for(i=0;i<N_FaceDOF;i++)
  {
    cout << "FaceDOF: " << FaceDOF[i] << endl;
    cout << "t1: " << FaceVectors[9*i+0] << " " << FaceVectors[9*i+1] << " " << FaceVectors[9*i+2] << endl;
    cout << "t2: " << FaceVectors[9*i+3] << " " << FaceVectors[9*i+4] << " " << FaceVectors[9*i+5] << endl;
    cout << "n: "  << FaceVectors[9*i+6] << " " << FaceVectors[9*i+7] << " " << FaceVectors[9*i+8] << endl;
    cout << "det: " << FaceVectors[9*i+0]*(FaceVectors[9*i+4]*FaceVectors[9*i+8]-FaceVectors[9*i+5]*FaceVectors[9*i+7])
             -FaceVectors[9*i+1]*(FaceVectors[9*i+3]*FaceVectors[9*i+8]-FaceVectors[9*i+5]*FaceVectors[9*i+6])
             +FaceVectors[9*i+2]*(FaceVectors[9*i+3]*FaceVectors[9*i+7]-FaceVectors[9*i+4]*FaceVectors[9*i+6]) << endl;
  }
  cout << endl;
  
  cout << "EdgeDOFs" << endl;
  for(i=0;i<N_EdgeDOF;i++)
  {
    cout << "EdgeDOF: " << EdgeDOF[i] << endl;
    cout << "n1: " << EdgeVectors[9*i+0] << " " << EdgeVectors[9*i+1] << " " << EdgeVectors[9*i+2] << endl;
    cout << "n2: " << EdgeVectors[9*i+3] << " " << EdgeVectors[9*i+4] << " " << EdgeVectors[9*i+5] << endl;
    cout << "t: "  << EdgeVectors[9*i+6] << " " << EdgeVectors[9*i+7] << " " << EdgeVectors[9*i+8] << endl;
    cout << "det: " << EdgeVectors[9*i+0]*(EdgeVectors[9*i+4]*EdgeVectors[9*i+8]-EdgeVectors[9*i+5]*EdgeVectors[9*i+7])
             -EdgeVectors[9*i+1]*(EdgeVectors[9*i+3]*EdgeVectors[9*i+8]-EdgeVectors[9*i+5]*EdgeVectors[9*i+6])
             +EdgeVectors[9*i+2]*(EdgeVectors[9*i+3]*EdgeVectors[9*i+7]-EdgeVectors[9*i+4]*EdgeVectors[9*i+6]) << endl;
  }
  cout << endl;
*/
  delete Flags;
} // end FindVectorsForSlipDOF

// ========================================================================
// manipulate square matrices due to u.n=0 constraint
// ========================================================================
void ManipulateSquareMatrices(TSquareStructure3D *sqstructure,
        double *a11, double *a12, double *a13,
        double *a21, double *a22, double *a23,
        double *a31, double *a32, double *a33,
        int N_FaceDOF, int N_EdgeDOF,
        int *FaceDOF, int *EdgeDOF,
        double* FaceVectors, double* EdgeVectors)
{
  int i,j,k;
  int DOF, Begin, End;
  int *RowPtr, *KCol;
  double v1, v2, v3;
  double t11, t12, t13, t21, t22, t23;
  double n1, n2, n3;
  double n11, n12, n13, n21, n22, n23;
  double t1, t2, t3;

  RowPtr = sqstructure->GetRowPtr();
  KCol = sqstructure->GetKCol();

  for(i=0;i<N_FaceDOF;i++)
  {
    DOF = FaceDOF[i];
    Begin = RowPtr[DOF];
    End = RowPtr[DOF+1];
    t11 = FaceVectors[9*i+0];
    t12 = FaceVectors[9*i+1];
    t13 = FaceVectors[9*i+2];

    t21 = FaceVectors[9*i+3];
    t22 = FaceVectors[9*i+4];
    t23 = FaceVectors[9*i+5];

    n1  = FaceVectors[9*i+6];
    n2  = FaceVectors[9*i+7];
    n3  = FaceVectors[9*i+8];

    // cout << DOF << endl;
    // cout << "t1: " << t11 << " " << t12 << " " << t13 << endl;
    // cout << "t2: " << t21 << " " << t22 << " " << t23 << endl;
    // cout << "n: " << n1 << " " << n2 << " " << n3 << endl;
    // cout << endl;
    
    for(j=Begin;j<End;j++)
    {
      // first column of blocks
      v1 = a11[j];
      v2 = a21[j];
      v3 = a31[j];

      a11[j] = t11*v1 + t12*v2 + t13*v3;
      a21[j] = t21*v1 + t22*v2 + t23*v3;

      // second column of blocks
      v1 = a12[j];
      v2 = a22[j];
      v3 = a32[j];

      a12[j] = t11*v1 + t12*v2 + t13*v3;
      a22[j] = t21*v1 + t22*v2 + t23*v3;

      // third column of blocks
      v1 = a13[j];
      v2 = a23[j];
      v3 = a33[j];

      a13[j] = t11*v1 + t12*v2 + t13*v3;
      a23[j] = t21*v1 + t22*v2 + t23*v3;

      if(KCol[j] == DOF)
      {
        a31[j] = n1;
        a32[j] = n2;
        a33[j] = n3;
      }
      else
      {
        a31[j] = 0;
        a32[j] = 0;
        a33[j] = 0;
      }
    } // endfor j
  } // endfor i

  for(i=0;i<N_EdgeDOF;i++)
  {
    DOF = EdgeDOF[i];
    Begin = RowPtr[DOF];
    End = RowPtr[DOF+1];
    n11 = EdgeVectors[9*i+0];
    n12 = EdgeVectors[9*i+1];
    n13 = EdgeVectors[9*i+2];

    n21 = EdgeVectors[9*i+3];
    n22 = EdgeVectors[9*i+4];
    n23 = EdgeVectors[9*i+5];

    t1  = EdgeVectors[9*i+6];
    t2  = EdgeVectors[9*i+7];
    t3  = EdgeVectors[9*i+8];
    for(j=Begin;j<End;j++)
    {
      // first column of blocks
      a11[j] = t1*a11[j] + t2*a21[j] + t3*a31[j];

      // second column of blocks
      a12[j] = t1*a12[j] + t2*a22[j] + t3*a32[j];

      // third column of blocks
      a13[j] = t1*a13[j] + t2*a23[j] + t3*a33[j];

      if(KCol[j] == DOF)
      {
        a21[j] = n11;
        a22[j] = n12;
        a23[j] = n13;

        a31[j] = n21;
        a32[j] = n22;
        a33[j] = n23;
      }
      else
      {
        a21[j] = 0;
        a22[j] = 0;
        a23[j] = 0;

        a31[j] = 0;
        a32[j] = 0;
        a33[j] = 0;
      }
    } // endfor j
  } // endfor i
} // end of ManipulateSquareMatrices

// ========================================================================
// manipulate square matrices due to u.n=0 constraint
// ========================================================================
void ManipulateMatricesAndRhs(TStructure3D *structure,
        double *b1t, double *b2t, double *b3t,
        double *f1, double *f2, double *f3,
        int N_FaceDOF, int N_EdgeDOF,
        int *FaceDOF, int *EdgeDOF,
        double* FaceVectors, double* EdgeVectors)
{
  int i,j,k;
  int DOF, Begin, End;
  int *RowPtr, *KCol;
  double v1, v2, v3, w1, w2, w3;
  double t11, t12, t13, t21, t22, t23;
  double t1, t2, t3;

  RowPtr = structure->GetRowPtr();
  KCol = structure->GetKCol();

  for(i=0;i<N_FaceDOF;i++)
  {
    DOF = FaceDOF[i];
    Begin = RowPtr[DOF];
    End = RowPtr[DOF+1];
    t11 = FaceVectors[9*i+0];
    t12 = FaceVectors[9*i+1];
    t13 = FaceVectors[9*i+2];

    t21 = FaceVectors[9*i+3];
    t22 = FaceVectors[9*i+4];
    t23 = FaceVectors[9*i+5];

    for(j=Begin;j<End;j++)
    {
      v1 = b1t[j];
      v2 = b2t[j];
      v3 = b3t[j];

      w1 = t11*v1 + t12*v2 + t13*v3;
      w2 = t21*v1 + t22*v2 + t23*v3;

      b1t[j] = w1;
      b2t[j] = w2;
      b3t[j] = 0;
    } // endfor j

    v1 = f1[DOF];
    v2 = f2[DOF];
    v3 = f3[DOF];

    w1 = t11*v1 + t12*v2 + t13*v3;
    w2 = t21*v1 + t22*v2 + t23*v3;

    f1[DOF] = w1;
    f2[DOF] = w2;
    f3[DOF] = 0;
  } // endfor i

  for(i=0;i<N_EdgeDOF;i++)
  {
    DOF = EdgeDOF[i];
    Begin = RowPtr[DOF];
    End = RowPtr[DOF+1];
    t1  = EdgeVectors[9*i+6];
    t2  = EdgeVectors[9*i+7];
    t3  = EdgeVectors[9*i+8];
    for(j=Begin;j<End;j++)
    {
      b1t[j] = t1*b1t[j] + t2*b2t[j] + t3*b3t[j];
      b2t[j] = 0;
      b3t[j] = 0;
    } // endfor j

    f1[DOF] = t1*f1[DOF] + t2*f2[DOF] + t3*f3[DOF];
    f2[DOF] = 0;
    f3[DOF] = 0;
  } // endfor i
} // end of ManipulateMatricesAndRhs

// find all joints which form the free surface
void FindFreeSurface(TCollection *Coll,
                     int &N_SurfaceJoints, int* &CellNumbers,
                     int* &JointNumbers)
{
  int i,j,k,l;
  int N_Cells, N_Joints, N_Joints1;
  TBaseCell *cell, *neigh, *neigh0, *neigh1;
  TJoint *joint;
  TIsoJointEqN *isojoint;

  k = 0;
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->GetSubGridID() == 1)
    {
      // fluid cell
      N_Joints = cell->GetN_Joints();
      for(j=0;j<N_Joints;j++)
      {
        joint = cell->GetJoint(j);
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          if(neigh->GetSubGridID() == 2)
          {
            // upper cell => joint is on free surface
            k++;
            if(joint->GetType() == JointEqN)
            {
              neigh0 = joint->GetNeighbour(0);
              neigh1 = joint->GetNeighbour(1);
              isojoint = new TIsoJointEqN(neigh0, neigh1);
              cell->SetJoint(j, isojoint);
              N_Joints1 = neigh->GetN_Joints();
              for(l=0;l<N_Joints1;l++)
                if(neigh->GetJoint(l) == joint)
                  neigh->SetJoint(l, isojoint);
//               delete (TJoint *)joint;
              isojoint->SetMapType();
            }
            else
            {
              Error("Wrong joint type, JointEqN = 1 excepted!" << endl);
              exit(-1);
            }
          }
        } // endif neigh
      } // endfor j
    } // end ID == 1
  } // endfor i

  N_SurfaceJoints = k;

  CellNumbers = new int[N_SurfaceJoints];
  JointNumbers = new int[N_SurfaceJoints];

  k = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->GetSubGridID() == 1)
    {
      // fluid cell
      N_Joints = cell->GetN_Joints();
      for(j=0;j<N_Joints;j++)
      {
        joint = cell->GetJoint(j);
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          if(neigh->GetSubGridID() == 2)
          {
            // upper cell => joint is on free surface
            CellNumbers[k] = i;
            JointNumbers[k] = j;
            k++;
          }
        } // endif neigh
      } // endfor j
    } // end ID == 1
  } // endfor i
} // end FindFreeSurface

void FindFreeSurfaceFromJointType(TCollection *Coll, JointType type,
				  int &N_SurfaceJoints, int* &CellNumbers,
				  int* &JointNumbers)
{
  int N_Cells, N_Joints, k;
  
  TBaseCell *Cell;
  TJoint *Joint;
  
  N_Cells = Coll->GetN_Cells();
  
  k = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    N_Joints = Cell->GetN_Joints();
    for (int j=0;j<N_Joints;++j)
    {
      Joint = Cell->GetJoint(j);
      
      if ( Joint->GetType() == type )
	++k;
    }
  }
  
  if ( k == 0 ) 
  {
    CellNumbers = NULL;
    JointNumbers = NULL;
    N_SurfaceJoints = 0;
    return;
  }
  
  N_SurfaceJoints = k;
  CellNumbers = new int [N_SurfaceJoints];
  JointNumbers = new int [N_SurfaceJoints];
  
  k = 0;
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    N_Joints = Cell->GetN_Joints();
    for (int j=0;j<N_Joints;++j)
    {
      Joint = Cell->GetJoint(j);
      
      if ( Joint->GetType() == type )
      {
	CellNumbers[k] = i;
	JointNumbers[k] = j;
	++k;
      }
    }
  }
}

// calculate normal vectors on free surface
void CalculateNormals(TCollection *Coll,
                     int N_SurfaceJoints, int *CellNumbers,
                     int *JointNumbers,
                     TFESpace3D *fespace,
                     double* &n1, double* &n2, double* &n3,
                     double* &len)
{
  int i,j,k,l,m;
  int N_Cells;
  int CellNr, JointNr;
  TCollection *coll;
  TBaseCell *cell;
  int N_DOF;
  FE3D FEId;
  TFE3D *ele;
  BF3DRefElements RefElement;
  RefTrans3D RefTrans;
  TRefTrans3D *F_K;
  QuadFormula2D QF2;
  TQuadFormula2D *qf2;
  QuadFormula3D QF3;
  TQuadFormula3D *qf3;
  int N_Points2, N_Points3;
  double *Weights2, *t1, *t2;
  double *Weights3, *xi, *eta, *zeta;
  int *GlobalNumbers, *BeginIndex, *DOF, *LocalJointDOF;
  int N_UsedElements = 1;
  BaseFunct3D BaseFuncts[1];
  bool SecondDer[1] = { FALSE };
  TFEDesc3D *fedesc;
  int N_JointDOF, JointDOF[MaxN_BaseFunctions3D];
  double **AllDerivX, **AllDerivY, **AllDerivZ, **AllValues;
  double *DerivX, *DerivY, *DerivZ, *Values;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double absdetjk[MaxN_QuadPoints_3D];
  double mult;
  double a1, a2, a3, b1, b2, b3;
  double N1, N2, N3, LEN;

  N_DOF = fespace->GetN_DegreesOfFreedom();
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  if(n1 == NULL)
  {
    n1 = new double[N_DOF];
    n2 = new double[N_DOF];
    n3 = new double[N_DOF];
    len = new double[N_DOF];
  }

  memset(n1, 0, N_DOF*SizeOfDouble);
  memset(n2, 0, N_DOF*SizeOfDouble);
  memset(n3, 0, N_DOF*SizeOfDouble);
  memset(len, 0, N_DOF*SizeOfDouble);

  for(i=0;i<N_SurfaceJoints;i++)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    cell = Coll->GetCell(CellNr);

    FEId = fespace->GetFE3D(CellNr, cell);
    ele = TFEDatabase3D::GetFE3D(FEId);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);

    BaseFuncts[0] = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(FEId);

    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

    switch(RefElement)
    {
      case BFUnitHexahedron:
        RefTrans = HexaIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);

        QF2 = TFEDatabase3D::GetQFQuadFromDegree(2*l);
        qf2 = TFEDatabase3D::GetQuadFormula2D(QF2);
        qf2->GetFormulaData(N_Points2, Weights2, t1, t2);

        QF3 = TFEDatabase3D::GetQFHexaFromDegree(2*l);
        qf3 = TFEDatabase3D::GetQuadFormula3D(QF3);
        qf3->GetFormulaData(N_Points3, Weights3, xi, eta, zeta);

        ((THexaIsoparametric *)F_K)->SetApproximationOrder(l);
        ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((THexaIsoparametric *)F_K)->SetCell(cell);
      break;

      case BFUnitTetrahedron:
        Error("Nothing is implemented for tetrahedra!" << endl);
        exit(-1);
      break;
    } // endswitch
      
    DOF = GlobalNumbers + BeginIndex[CellNr];
    fedesc = ele->GetFEDesc3D();
    LocalJointDOF = fedesc->GetJointDOF(JointNr);
    N_JointDOF = fedesc->GetN_JointDOF();
    for(j=0;j<N_JointDOF;j++)
      JointDOF[j] = DOF[LocalJointDOF[j]];

    ((THexaIsoparametric *)F_K)->GetOrigValues(N_UsedElements,
                             BaseFuncts,
                             N_Points3, xi, eta, zeta,
                             QF3, SecondDer);

    ((THexaIsoparametric *)F_K)->GetOrigFromRef(N_Points3, xi, eta, zeta,
                                         X, Y, Z, absdetjk);

    AllDerivX = TFEDatabase3D::GetOrigElementValues(BaseFuncts[0], D100);
    AllDerivY = TFEDatabase3D::GetOrigElementValues(BaseFuncts[0], D010);
    AllDerivZ = TFEDatabase3D::GetOrigElementValues(BaseFuncts[0], D001);

    for(j=0;j<N_Points3;j++)
    {
      DerivX = AllDerivX[j]; 
      DerivY = AllDerivY[j]; 
      DerivZ = AllDerivZ[j]; 

      mult = absdetjk[j]*Weights3[j];

      for(k=0;k<N_JointDOF;k++)
      {
        l = JointDOF[k];
        m = LocalJointDOF[k];

        n1[l] += mult*DerivX[m];
        n2[l] += mult*DerivY[m];
        n3[l] += mult*DerivZ[m];
      } // endfor k
    } // endfor j

    TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)->MakeRefElementData(QF2);

    AllValues = TFEDatabase3D::GetJointValues3D(BaseFuncts[0],
                                                QF2, JointNr);

    for(j=0;j<N_Points2;j++)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetTangentVectors(
                JointNr, t1[j], t2[j], a1, a2, a3, b1, b2, b3);
        break;

        case BFUnitTetrahedron:
          Error("Nothing is implemented for tetrahedra!" << endl);
          exit(-1);
        break;
      } // endswitch

      N1 = a2*b3 - a3*b2;
      N2 = a3*b1 - a1*b3;
      N3 = a1*b2 - a2*b1;

      LEN = sqrt(N1*N1 + N2*N2 + N3*N3);

      Values = AllValues[j];

      mult = LEN*Weights2[j];

      for(k=0;k<N_JointDOF;k++)
      {
        l = JointDOF[k];
        m = LocalJointDOF[k];

        len[l] += mult*Values[m];
      } // endfor k
    } // endfor j

  } // endfor i

  for(i=0;i<N_DOF;i++)
  {
    if(fabs(len[i])>1e-12)
    {
      n1[i] /= len[i];
      n2[i] /= len[i];
      n3[i] /= len[i];
    }
    else
      len[i] = 0;
  }

  /*
  for(i=0;i<N_DOF;i++)
  {
    if(fabs(n1[i])<1e-12) n1[i] = 0;
    if(fabs(n2[i])<1e-12) n2[i] = 0;
    if(fabs(n3[i])<1e-12) n3[i] = 0;
    OutPut(setw(4) << i << setw(13) << n1[i] << setw(13) << n2[i] << setw(13) << n3[i]);
    OutPut(setw(13) << sqrt(n1[i]*n1[i]+n2[i]*n2[i]+n3[i]*n3[i]));
    OutPut(setw(15) << len[i] << endl);
  }
  */
} // end CalculateNormalVectors

// ========================================================================
// ATTENTION !!!
// the current implementation does not check the boundary condition
// it is assumed that all given faces are at the free boundary !!!
// ========================================================================
void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
                 int *CellNumbers, int *JointNumbers,
                 TFEFunction3D *potential, double dt,
                 TSquareMatrix3D *Aii,
                 double *rhs1, double *rhs2, double *rhs3)
{
  int i,j,k,l,m;
  int CellNr, JointNr;
  TBaseCell *cell;
  TJoint *joint;
  double InvWe;
  int N_BaseFunct, *N_BaseFuncts;
  BaseFunct3D *BaseFuncts;
  TFESpace3D *fespace;
  int *BeginIndex, *GlobalNumbers;
  int *RowPtr, *KCol;
  double *ValuesAii;
  FE3D FEId;
  TFE3D *ele;
  RefTrans3D RefTrans;
  TRefTrans3D *F_K;
  BF3DRefElements RefElement;
  double Param1[4], Param2[4];
  QuadFormula2D QuadFormula;
  QuadFormula3D QF3;
  TQuadFormula2D *qf;
  int N_Points;
  double *Weights, *p1, *p2;
  double **uref, **uxiref, **uetaref, **uzetaref;
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  int *DOF, TestDOF, AnsatzDOF;
  double a1, a2, a3, b1, b2, b3;
  double n1, n2, n3, len, val;
  double d1, d2, d3, ngrad;
  double e1, e2, e3, ngrad2;
  int index1, index2;

  TCollection *CollPot;
  int N_CellsPot, CellNrPot;
  int *BeginIndexPot, *GlobalNumbersPot;
  TFESpace3D *fespacePot;
  FE3D FEIdPot;
  double *valuesPot, value;
  int *DOFPot;
  double H, Hx, Hy, Hz, L;
  double mag1, mag2, Hn, arg;

  double U = TDatabase::ParamDB->FS_U;
  double rho = TDatabase::ParamDB->FS_RHO;
  double Ms = TDatabase::ParamDB->FS_MS;
  double chi0 = TDatabase::ParamDB->FS_CHI0;
  double gamma = TDatabase::ParamDB->FS_GAMMA;
  double Factor = 0.4*Pi*Ms*Ms/(U*U*rho); // for dimensionless equation

  fespacePot = potential->GetFESpace3D();
  BeginIndexPot = fespacePot->GetBeginIndex();
  GlobalNumbersPot = fespacePot->GetGlobalNumbers();

  valuesPot = potential->GetValues();

  CollPot = fespacePot->GetCollection();
  N_CellsPot = CollPot->GetN_Cells();
  for(i=0;i<N_CellsPot;i++)
    CollPot->GetCell(i)->SetClipBoard(i);

  InvWe = 1.0/TDatabase::ParamDB->FS_WE;
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  fespace = Aii->GetFESpace();
  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  ValuesAii = Aii->GetEntries();

  RowPtr = Aii->GetRowPtr();
  KCol = Aii->GetKCol();

  for(i=0;i<N_BoundFaces;i++)
  {
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    cell = Coll->GetCell(CellNr);
    joint = cell->GetJoint(JointNr);
    // OutPut(endl << "cell: " << CellNr << " joint: " << JointNr << endl);

    CellNrPot = cell->GetClipBoard();
    FEIdPot = fespacePot->GetFE3D(CellNrPot, cell);

    FEId = fespace->GetFE3D(CellNr, cell);

    if(FEId != FEIdPot)
    {
      Error("Finite elements do not match!" << endl);
      Error(__FILE__ << " at line " << __LINE__ << endl);
      exit(-1);
    }
    
    ele = TFEDatabase3D::GetFE3D(FEId);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FEId);

    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FEId);

    switch(RefElement)
    {
      case BFUnitHexahedron:
        RefTrans = HexaIsoparametric;
        F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);

        QF3 = TFEDatabase3D::GetQFHexaFromDegree(2*l);
        ((THexaIsoparametric *)F_K)->SetQuadFormula(QF3);
        ((THexaIsoparametric *)F_K)->SetApproximationOrder(l);
        ((THexaIsoparametric *)F_K)->SetCell(cell);
      break;

      case BFUnitTetrahedron:
        Error("Nothing is implemented for tetrahedra!" << endl);
        exit(-1);
      break;
    } // endswitch

    QuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
    qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
    qf->GetFormulaData(N_Points, Weights, p1, p2);

    TFEDatabase3D::GetBaseFunct3DFromFE3D(FEId)
        ->MakeRefElementData(QuadFormula);

    DOF = GlobalNumbers + BeginIndex[CellNr];
    DOFPot = GlobalNumbersPot + BeginIndexPot[CellNrPot];
    N_BaseFunct = N_BaseFuncts[FEId];

    for(k=0;k<N_Points;k++)
    {
      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetTangentVectors(
                JointNr, p1[k], p2[k], a1, a2, a3, b1, b2, b3);
        break;

        case BFUnitTetrahedron:
          Error("Nothing is implemented for tetrahedra!" << endl);
          exit(-1);
        break;
      } // endswitch

      n1 = a2*b3 - a3*b2;
      n2 = a3*b1 - a1*b3;
      n3 = a1*b2 - a2*b1;

      // OutPut("xi: " << p1[k] << " eta: " << p2[k] << endl);
      // OutPut("point: " << k << " t1: " << a1 << " " << a2 << " " << a3 << endl);
      // OutPut("point: " << k << " t2: " << b1 << " " << b2 << " " << b3 << endl);

      len = sqrt(n1*n1 + n2*n2 + n3*n3);
      // OutPut("len: " << len << endl);

      n1 /= len;
      n2 /= len;
      n3 /= len;

      // OutPut("len: " << len << endl);
      // OutPut("point: " << k << " normal: " << n1 << " " << n2 << " " << n3 << endl);

      uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FEId],
                QuadFormula, JointNr);
      uxiref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, JointNr, D100);
      uetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, JointNr, D010);
      uzetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FEId],
                QuadFormula, JointNr, D001);

      switch(RefElement)
      {
        case BFUnitHexahedron:
          ((THexaIsoparametric *)F_K)->GetOrigValues(
                JointNr, p1[k], p2[k], N_BaseFunct,
                uref[k], uxiref[k], uetaref[k], uzetaref[k],
                uorig, uxorig, uyorig, uzorig);
        break;

        case BFUnitTetrahedron:
          Error("Nothing is implemented for tetrahedra!" << endl);
          exit(-1);
        break;
      } // endswitch

      // calculate dimensionless magnetic field
      // mag1: int(M dH)
      // mag2: (M.n)^2
      if(fabs(chi0) > 1e-8)
      {
        Hx = 0; Hy = 0; Hz = 0;
        for(l=0;l<N_BaseFunct;l++)
        {
          value = valuesPot[DOFPot[l]];
          Hx += value*uxorig[l];
          Hy += value*uyorig[l];
          Hz += value*uzorig[l];
        }
        H = sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

        if(fabs(H) > 1e-8)
        {
          Hn = (Hx*n1 + Hy*n2 + Hz*n3) / H;
        }
        else
        {
          Hn = 0;
        }

        arg = gamma*H;
        if(arg>1e-2)
        {
          mag1 = log( sinh(arg)/arg ) / (3*chi0);
          L = ( 1/tanh(arg) - 1/arg );
        }
        else
        {
          L = arg/3;
          mag1 = arg*arg/6;
        }
        mag2 = 0.5 * L*L * Hn*Hn;
      }
      else
      {
        mag1 = 0;
        mag2 = 0;
      }

      for(l=0;l<N_BaseFunct;l++)
      {
        TestDOF = DOF[l];

        ngrad = uxorig[l]*n1 + uyorig[l]*n2 + uzorig[l]*n3;
        d1 = uxorig[l] - ngrad*n1;
        d2 = uyorig[l] - ngrad*n2;
        d3 = uzorig[l] - ngrad*n3;

        // NOTE: dt from time integration
        val = InvWe * ( (1-n1*n1)*d1 + (0-n1*n2)*d2 + (0-n1*n3)*d3 );
        val -= Factor * (uorig[l]*n1) * (mag1 + mag2);
        val *= Weights[k] * len;
        rhs1[TestDOF] -= val;

        val = InvWe * ( (0-n2*n1)*d1 + (1-n2*n2)*d2 + (0-n2*n3)*d3 );
        val -= Factor * (uorig[l]*n2) * (mag1 + mag2);
        val *= Weights[k] * len;
        rhs2[TestDOF] -= val;

        val = InvWe * ( (0-n3*n1)*d1 + (0-n3*n2)*d2 + (1-n3*n3)*d3 );
        val -= Factor * (uorig[l]*n3) * (mag1 + mag2);
        val *= Weights[k] * len;
        rhs3[TestDOF] -= val;

        index2 = RowPtr[TestDOF+1];
        for(m=0;m<N_BaseFunct;m++)
        {
          AnsatzDOF = DOF[m];
          index1 = RowPtr[TestDOF];
          while(KCol[index1] != AnsatzDOF) index1++;

          ngrad2 = uxorig[m]*n1 + uyorig[m]*n2 + uzorig[m]*n3;
          e1 = uxorig[m] - ngrad2*n1;
          e2 = uyorig[m] - ngrad2*n2;
          e3 = uzorig[m] - ngrad2*n3;

          val = d1*e1 + d2*e2 + d3*e3;

          // NOTE: one dt from time integration
          //       the other from implicit treatment
          val *= Weights[k] * len * dt * InvWe;
          ValuesAii[index1] += val;
        } // endfor m
      } // endfor l
    } // endfor k
  } // endfor i
} // end FreeSurfInt

void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
		 int *CellNumbers, int *JointNumbers,
		 double dt, TSquareMatrix3D **Aii,
		 double *rhs1, double *rhs2, double *rhs3)
{
  int CellNr, JointNr, l, m, N_Points, N_JointDOF, *JointDOF;
  int *N_BaseFuncts, N_BaseFunct;
  double *weights, *p1, *p2, xi, eta, zeta, X, Y, Z;
  double a1, a2, a3, b1, b2, b3;
  double n1, n2, n3, len, fact;
  double **uref, **uxiref, **uetaref, **uzetaref;
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double test001, test010, test100, ansatz001, ansatz010, ansatz100;
  double ngrad1, ngrad2, val;
  double **EntriesAii;
  double invWe = TDatabase::ParamDB->WB_NR;
  double valx, valy, valz, surf=0;
  int *KColAii;
  int *RowPtrAii;
  int *BeginIndex, *GlobalNumbers, *DOF;
  int TestDOF, AnsatzDOF, index;
  const int *TmpFV, *TmpLen;
  int MaxLen;
  int addlhs = (int) TDatabase::ParamDB->P2;
  int sphere = (int) TDatabase::ParamDB->P3;
  int normal_change_count=0;
  bool Isoparametric = FALSE;
  
  if ( sphere == 1 ) addlhs = 0;
  
  BaseFunct3D *BaseFuncts;
  TBaseCell *Cell;
  TFESpace3D *fesp;
  BF3DRefElements RefElement;
  FE3D FeID;
  TFE3D *ele;
  RefTrans3D RefTrans;
  TRefTrans3D *F_K;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2d;
  TJoint *Joint;
   
  OutPut("Adding surface term ");
  if ( sphere == 1 )
  {
    OutPut("with exact curvature ... ");
  }
  
  fesp = Aii[0]->GetFESpace();
  KColAii = Aii[0]->GetKCol();
  RowPtrAii = Aii[0]->GetRowPtr();
  
  switch (TDatabase::ParamDB->NSTYPE)
  {
    case 2:
      EntriesAii = new double* [1];
      EntriesAii[0] = Aii[0]->GetEntries();
      
      break;
      
    case 4:
      EntriesAii = new double* [3];
      
      EntriesAii[0] = Aii[0]->GetEntries();
      EntriesAii[1] = Aii[1]->GetEntries();
      EntriesAii[2] = Aii[2]->GetEntries();
      
      break;
  }  
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  GlobalNumbers = fesp->GetGlobalNumbers();
  BeginIndex    = fesp->GetBeginIndex();
  
  if ( invWe == 0.0) 
  {
    OutPut("We = 0, skipping surface integral!" << endl);
    return;
  }
  else invWe = 1.0/invWe; 
  
  for (int i=0;i<N_BoundFaces;++i)
  {   
    CellNr = CellNumbers[i];
    JointNr = JointNumbers[i];
    Cell = Coll->GetCell(CellNr);
    Joint = Cell->GetJoint(JointNr);
    
    if ( Joint->GetType() == IsoBoundFace )
    {
      Isoparametric = TRUE;
//       cout << "iso" << endl;
    }
    
    FeID = fesp->GetFE3D(CellNr, Cell);
    ele = TFEDatabase3D::GetFE3D(FeID);
    DOF = GlobalNumbers + BeginIndex[CellNr];
    
    N_JointDOF = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID)->GetN_JointDOF();
    JointDOF = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID)->GetJointDOF(JointNr);
    
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FeID);
    
    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FeID);
    
    switch (RefElement)
    {
      case BFUnitTetrahedron:
	if ( Isoparametric )
	{
	  RefTrans = TetraIsoparametric;
	  F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
	  
	  
	  ((TTetraIsoparametric*) F_K)->SetQuadFormula(P2Tetra);
	  ((TTetraIsoparametric*) F_K)->SetApproximationOrder(2);
	  ((TTetraIsoparametric*) F_K)->SetCell(Cell);
	}
	else
	{
	  RefTrans = TetraAffin;
	  F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
	  
	  ((TTetraAffin*) F_K)->SetCell(Cell);
	}
// 	QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
	QuadFormula = Degree19Tria;
	break;
	
      case BFUnitHexahedron:
	cerr << __FILE__ << ":" << __LINE__ << ": ";
	cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
	exit(0);
	break;
	
      default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
    }	
    
    qf2d = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
    qf2d->GetFormulaData(N_Points, weights, p1, p2);
    
    TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID)
        ->MakeRefElementData(QuadFormula);
	
    N_BaseFunct = N_BaseFuncts[FeID];

    uref = TFEDatabase3D::GetJointValues3D(BaseFuncts[FeID],
                QuadFormula, JointNr);
    uxiref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FeID],
                QuadFormula, JointNr, D100);
    uetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FeID],
                QuadFormula, JointNr, D010);
    uzetaref = TFEDatabase3D::GetJointDerivatives3D(BaseFuncts[FeID],
                QuadFormula, JointNr, D001);
		
//     Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen); 
//     for (m=0;m<TmpLen[JointNr];++m)
//     {	
//       double n1, n2, n3;
// 	
//       ((TTetraAffin*) F_K)->GetOuterNormal(JointNr, 0, 0, n1, n2, n3);
//       Cell->GetVertex(TmpFV[JointNr*MaxLen+m])->AddNormal(n1, n2, n3);
//     }

    for (int k=0;k<N_Points;++k)
    {
      switch (RefElement)
      {
	case BFUnitTetrahedron:	  
	  if ( Isoparametric )
	  {
	    ((TTetraIsoparametric*) F_K)->GetTangentVectors(JointNr, p1[k], p2[k],
							    a1, a2, a3, b1, b2, b3);
							    
	    ((TTetraIsoparametric*) F_K)->GetOrigBoundFromRef(JointNr, p1[k], p2[k],
							      X, Y, Z);
	  }
	  else
	  {
	    switch (JointNr)
	    {
	      case 0:
		xi = p1[k];
		eta = p2[k];
		zeta = 0;
		break;
	      case 1:
		xi = p2[k];
		eta = 0;
		zeta = p1[k];
		break;
	      case 2:
		xi = p1[k];
		eta = 1-p1[k]-p2[k];
		zeta = p2[k];
		break;
	      case 3:
		xi = 0;
		eta = p1[k];
		zeta = p2[k];
		break;
	    }	    
	    
	    ((TTetraAffin*) F_K)->GetTangentVectors(JointNr, p1[k], p2[k],
						    a1, a2, a3, b1, b2, b3);
	    ((TTetraAffin*) F_K)->GetOrigFromRef(xi, eta, zeta, X, Y, Z);
	  }
						   
	  
	  break;
	  
	case BFUnitHexahedron:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
	  cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
	  exit(0);
	  break;
	  
	default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
      } // endswitch
      
      n1 = (a2*b3 - a3*b2);
      n2 = (a3*b1 - a1*b3);
      n3 = (a1*b2 - a2*b1);
      
      len = sqrt(n1*n1 + n2*n2 + n3*n3);
      
      n1 /= len;
      n2 /= len;
      n3 /= len;
      
//       X /= sqrt(X*X+Y*Y+Z*Z);
//       Y /= sqrt(X*X+Y*Y+Z*Z);
//       Z /= sqrt(X*X+Y*Y+Z*Z);
     
      surf += weights[k] * len /** X * n1*/;

      switch(RefElement)
      {
        case BFUnitTetrahedron:
	  if ( Isoparametric )
	  {
	    ((TTetraIsoparametric*) F_K)->GetOrigValues(JointNr, p1[k], p2[k], N_BaseFunct,
							uref[k], uxiref[k], uetaref[k], uzetaref[k],
							uorig, uxorig, uyorig, uzorig);
	  }
	  else
	  {
	    ((TTetraAffin*) F_K)->GetOrigValues(0.0, 0.0, 0.0, N_BaseFunct,
					      uref[k], uxiref[k], uetaref[k], uzetaref[k],
					      uorig, uxorig, uyorig, uzorig);				      
	  }

        break;

        case BFUnitHexahedron:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
	  exit(0);
        break;
	
	default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
      } // endswitch
      
      fact = weights[k] * len * invWe;
      
      for (l=0;l<N_BaseFunct;++l)
      {	
	TestDOF = DOF[l];
	
// 	if ( IsJointDOF(l, N_JointDOF, JointDOF) == 0 )
// 	{
// 	  continue;
// 	}
	
	valx = uxorig[l];
	valy = uyorig[l];
	valz = uzorig[l];
	
	ngrad1 = valx*n1 + valy*n2 + valz*n3;
	test100 = valx - ngrad1*n1;
	test010 = valy - ngrad1*n2;
        test001 = valz - ngrad1*n3;
	
	if ( sphere == 0 )
	{
	  // rhs (dt due to time integration)
	  val = (1.0-n1*n1)*test100 + (0.0-n1*n2)*test010 + (0.0-n1*n3)*test001;
	  val *= 0.5*fact;
  // 	val *= n1;
	  rhs1[TestDOF] -= dt*val;
	  
	  val = (0.0-n2*n1)*test100 + (1.0-n2*n2)*test010 + (0.0-n2*n3)*test001;
	  val *= 0.5*fact;
  // 	val *= n2;
	  rhs2[TestDOF] -= dt*val;
	  
	  val = (0.0-n3*n1)*test100 + (0.0-n3*n2)*test010 + (1.0-n3*n3)*test001;
	  val *= 0.5*fact;
  // 	val *= n3;
	  rhs3[TestDOF] -= dt*val;
	}
	else
	{	  
	  val = fact*n1*uorig[l];
	  rhs1[TestDOF] += dt*val;
	  
	  val = fact*n2*uorig[l];
	  rhs2[TestDOF] += dt*val;
	  
	  val = fact*n3*uorig[l];
	  rhs3[TestDOF] += dt*val;
	}
	
	for (m=0;m<N_BaseFunct && addlhs == 1;++m)
	{	  
	  AnsatzDOF = DOF[m];
	  for(index=RowPtrAii[TestDOF];index<RowPtrAii[TestDOF+1];++index)
	  {
	    if ( KColAii[index] == AnsatzDOF )
	      break;
	  }
	  
	  if ( index == RowPtrAii[TestDOF+1] )
	  {
	    cerr << __FILE__ << ":" << __LINE__ << ": dof not found. Exit!" << endl;
	    exit(0);
	  }
	  
	  valx = uxorig[m];
	  valy = uyorig[m];
	  valz = uzorig[m];
	  
	  ngrad2 = valx*n1 + valy*n2 + valz*n3;
          ansatz100 = valx - ngrad2*n1;
          ansatz010 = valy - ngrad2*n2;
          ansatz001 = valz - ngrad2*n3;
	  
	  val = test001*ansatz001 + test010*ansatz010 + test100*ansatz100;
	  val *= fact*dt*dt; // one dt due to time integration
	  
	  switch (TDatabase::ParamDB->NSTYPE)
	  {
	    case 4:
	      EntriesAii[1][index] += val;
	      EntriesAii[2][index] += val;
	    case 2:
	      EntriesAii[0][index] += val;
	      
	      break;
	  }
	} // end for(m=0; ...
      } // end for(l=0; ...)
    } // end for (int k=0; ...
  } // end for boundfaces
  
  OutPut("done: " << surf <<endl);
//   cout << "leave FreeSurfInt()" << endl;
}

void FreeSurfInt_new(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers,
		     double dt, TSquareMatrix3D **Aii, double *rhs1, double *rhs2, double *rhs3)
{
  int CellNr, JointNr, N_QuadPoints, N_BaseFunct, N_Cells, N_Points;
  int *GlobalNumbers, *BeginIndex, *DOF, TestDOF;
  double *p1, *p2, *weights, len, fact;
  double a1, a2, a3, b1, b2, b3, n1, n2, n3;
  double **uref, *urefvalues, val;
  double values[MaxN_BaseFunctions3D];
  double x, y, z, *xi, *eta, *zeta;
  TBaseCell *Cell;
  TFESpace3D *fespace;
  FE3D FeID;
  TBaseFunct3D *bf;
  TTetraAffin *F_aff;
  TTetraIsoparametric *F_iso;
  bool iso = FALSE;
  TQuadFormula2D *qf2d;
  TNodalFunctional3D *nf;
  double invWe = TDatabase::ParamDB->WB_NR;
  double t[MaxN_QuadPoints_3D], s[MaxN_QuadPoints_3D];
  double Xi, Eta, Zeta;
  
  invWe = 1. / invWe;
  
//   OutPut("Adding surface integral (new) ... ");
  
  fespace = Aii[0]->GetFESpace();
  
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex    = fespace->GetBeginIndex();
  
  // alloc arrays
  uref = new double* [MaxN_QuadPoints_3D];
  for (int i=0;i<MaxN_QuadPoints_3D;++i)
    uref[i] = new double [MaxN_BaseFunctions3D];
  
  N_Cells = Coll->GetN_Cells();
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    int N_BoundJoints = 0;
    for (int j=0;j<4;++j)
    {
      if ( Cell->GetJoint(j)->GetType() == BoundaryFace )
      {
	N_BoundJoints++;
	JointNr = j;
      }
    }
    
    if ( N_BoundJoints > 1 )
    {
      cerr << "askdjkaksjdkl " << endl;
      exit(0);
    }
    
    if ( N_BoundJoints == 1 )
    {
      F_aff = (TTetraAffin*) TFEDatabase3D::GetRefTrans3D(TetraAffin);
      F_aff->SetCell(Cell);
      
//       OutPut("JointNr: " << JointNr << endl);
      
      FeID = fespace->GetFE3D(i, Cell);
      bf = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
      nf = TFEDatabase3D::GetNodalFunctional3DFromFE3D(FeID);
      N_BaseFunct = bf->GetDimension();
      
      qf2d = TFEDatabase3D::GetQuadFormula2D(Degree8Tria);
      qf2d->GetFormulaData(N_QuadPoints, weights, p1, p2);
      
      bf->GetValues(N_QuadPoints, p1, p2, JointNr, uref);
      
      DOF = GlobalNumbers + BeginIndex[i];
      
      for (int k=0;k<N_QuadPoints;++k)
      {
	F_aff->GetOrigBoundFromRef(JointNr, p1[k], p2[k], x, y, z);
	F_aff->GetTangentVectors(JointNr, p1[k], p2[k],
				 a1, a2, a3, b1, b2, b3);
				 
	len = sqrt(x*x+y*y+z*z);
	x /= len;
	y /= len;
	z /= len;
	
	n1 = a2*b3 - a3*b2;
	n2 = a3*b1 - a1*b3;
	n3 = a1*b2 - a2*b1;
	
	len = sqrt(n1*n1+n2*n2+n3*n3);
	
	n1 /= len;
	n2 /= len;
	n3 /= len;
	
	fact = weights[k]*len*invWe;
	for (int j=0;j<N_BaseFunct;++j)
	{
	  TestDOF = DOF[j];
	  
	  val = dt*fact*x*uref[k][j];
	  rhs1[TestDOF] += val;
	  
	  val = dt*fact*y*uref[k][j];
	  rhs2[TestDOF] += val;
	  
	  val = dt*fact*z*uref[k][j];
	  rhs3[TestDOF] += val;
	}
	
// 	OutPut(weights[k]<<"  "<<p1[k]<<"  "<<p2[k]<<endl);
      }
      
//       OutPut(rhs1[DOF[0]]<<endl);
//       exit(0);
    }    
  }
  
  for (int i=0;i<MaxN_QuadPoints_3D;++i)
    delete [] uref[i];
  
  delete [] uref;
  
  OutPut("done"<<endl);
}
    
void FreeSurfInt_Sphere(TFESpace3D *fespace, double dt,
			double *rhs1, double *rhs2, double *rhs3)
{
  int N_Cells, l, N_Points, *N_BaseFuncts, N_BaseFunct;
  int TestDOF, AnsatzDOF, InnerBound, N_DOF;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double *weights, *xi, *eta, *zeta;
  double **uref, **uxiref, **uetaref, **uzetaref;
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double val, curvature, fact, scale;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D], absdetjk[MaxN_QuadPoints_3D];
  TBaseCell *Cell;
  TCollection *Coll;
  FE3D FeID;
  RefTrans3D RefTrans;
  TTetraAffin *F_aff;
  TTetraIsoparametric *F_iso;
  BF3DRefElements RefElement;
  QuadFormula3D QuadFormula;
  TQuadFormula3D *qf3d;
  BaseFunct3D *BaseFuncts;
  double invWe = TDatabase::ParamDB->WB_NR;
  double vol=0;
  bool iso = FALSE;
  
  if ( invWe == 0.0) 
  {
    OutPut("We = 0, skipping surface integral!" << endl);
    return;
  }
  else invWe = 1.0/invWe; 
  
  OutPut("Adding exact curvature ... ");
  
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  N_DOF = fespace->GetN_DegreesOfFreedom();
  InnerBound = fespace->GetInnerBound();
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex    = fespace->GetBeginIndex();
  
  curvature = 1.0;
  scale = dt*curvature*invWe;
  
  for (int i=0;i<N_Cells;++i)
  {
    Cell = Coll->GetCell(i);
    
    FeID = fespace->GetFE3D(i, Cell);
    RefElement = TFEDatabase3D::GetRefElementFromFE3D(FeID);
    l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(FeID);
    
    DOF = GlobalNumbers + BeginIndex[i];
    
    // check joints
    l = Cell->GetN_Joints();
    for (int j=0;j<l;++j)
    {
      if ( Cell->GetJoint(j)->GetType() == IsoBoundFace )
      {
	iso = TRUE;
	break;
      }
    }
    
    switch (RefElement)
    {
      case BFUnitTetrahedron:
	if (iso)
	{
// 	  OutPut("iso"<<endl);
	  RefTrans = TetraIsoparametric;
	  F_iso = (TTetraIsoparametric*) TFEDatabase3D::GetRefTrans3D(RefTrans);
	  
	  QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(2*l);
	  QuadFormula = P8Tetra;
	  
	  F_iso->SetApproximationOrder(2);
	  F_iso->SetQuadFormula(QuadFormula);
	  F_iso->SetCell(Cell);
	}
	else
	{
	  RefTrans = TetraAffin;
	  F_aff = (TTetraAffin*) TFEDatabase3D::GetRefTrans3D(RefTrans);
	  
	  QuadFormula = TFEDatabase3D::GetQFTetraFromDegree(2*l);
	  QuadFormula = P8Tetra;
	  
	  F_aff->SetCell(Cell);
	}
	
	break;
	
      case BFUnitHexahedron:
	cerr << __FILE__ << ":" << __LINE__ << ": ";
	cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
	exit(0);
	break;
	
      default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
    }
    
//     switch (QuadFormula)
//     {
//       case P2Tetra:
// 	OutPut("P2Tetra" << endl);
// 	break;
// 	
//       case P4Tetra:
// 	OutPut("P4Tetra" << endl);
// 	break;
// 	
//       case P5Tetra:
// 	OutPut("P5Tetra" << endl);
// 	break;
// 	
//       case P8Tetra:
// 	OutPut("P8Tetra" << endl);
// 	break;
// 	
//       default:
// 	OutPut("wrong quad: " << QuadFormula << endl);
//     }
    
    qf3d = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
    qf3d->GetFormulaData(N_Points, weights, xi, eta, zeta);
    
    TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID)->MakeRefElementData(QuadFormula);
    
    uref = TFEDatabase3D::GetRefElementValues(BaseFuncts[FeID], QuadFormula, D000);
    uxiref = TFEDatabase3D::GetRefElementValues(BaseFuncts[FeID], QuadFormula, D100);
    uetaref = TFEDatabase3D::GetRefElementValues(BaseFuncts[FeID], QuadFormula, D010);
    uzetaref = TFEDatabase3D::GetRefElementValues(BaseFuncts[FeID], QuadFormula, D001);
    
    N_BaseFunct = N_BaseFuncts[FeID];
    
    switch (RefElement)
    {
      case BFUnitTetrahedron:
	if (iso)
	{
	  F_iso->GetOrigFromRef(N_Points, xi, eta, zeta,
				X, Y, Z, absdetjk);
	}
	else
	{
	  F_aff->GetOrigFromRef(N_Points, xi, eta, zeta,
				X, Y, Z, absdetjk);
	}
	break;
	
      default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
    }
    
    // quadrature
    for (int k=0;k<N_Points;++k)
    {
      switch(RefElement)
      {
        case BFUnitTetrahedron:
	  if (iso)
	  {
	    F_iso->GetOrigValues(xi[k], eta[k], zeta[k], N_BaseFunct,
				uref[k], uxiref[k], uetaref[k], uzetaref[k],
				uorig, uxorig, uyorig, uzorig);
	  }
	  else
	  {
	    F_aff->GetOrigValues(xi[k], eta[k], zeta[k], N_BaseFunct,
				uref[k], uxiref[k], uetaref[k], uzetaref[k],
				uorig, uxorig, uyorig, uzorig);
	  }
	  break;

        case BFUnitHexahedron:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
	  exit(0);
        break;
	
	default:
	  cerr << __FILE__ << ":" << __LINE__ << ": ";
          cerr << "something went wrong. Exit!" << endl;
	  exit(0);
      } // endswitch 
    
      fact = scale*weights[k]*absdetjk[k];
    
      vol += weights[k]*absdetjk[k];
      
      for(int m=0;m<N_BaseFunct;++m)
      {
	TestDOF = DOF[m];
	
	val = fact*uxorig[m];
	rhs1[TestDOF] += val;
	
	val = fact*uyorig[m];
	rhs2[TestDOF] += val;
	
	val = fact*uzorig[m];
	rhs3[TestDOF] += val;
      }
    } // end for k
  }
  
//   for (int i=0;i<N_DOF;++i)
//   {
//     if ( i<InnerBound )
//     {
//       OutPut("inner: " << rhs1[i] << ", " << rhs2[i] << ", " << rhs3[i] << endl);
//     }
//     else
//     {
//       OutPut("boundary: " << rhs1[i] << ", " << rhs2[i] << ", " << rhs3[i] << endl);
//     }
//   }
  
  OutPut("done: " << vol <<endl);
}

int CheckNormal(double n1, double n2, double n3, TBaseCell *Cell, int JointNr)
{
  double barycenter[3], facebarycenter[3];
  const int *TmpFV, *TmpLen;
  int MaxLen, index;
  double x, y, z, scalar;
  
  barycenter[0] = 0.0;
  barycenter[1] = 0.0;
  barycenter[2] = 0.0;
  
  facebarycenter[0] = 0.0;
  facebarycenter[1] = 0.0;
  facebarycenter[2] = 0.0;
  
  Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
  
  for (int i=0;i<4;++i)
  {
    Cell->GetVertex(i)->GetCoords(x,y,z);
    
    barycenter[0] += x;
    barycenter[1] += y;
    barycenter[2] += z;
  }
  
  for (int i=0;i<3;++i)
  {
    index = TmpFV[MaxLen*JointNr+i];
    
    Cell->GetVertex(index)->GetCoords(x,y,z);
    
    facebarycenter[0] += x;
    facebarycenter[1] += y;
    facebarycenter[2] += z;
  }  
  
  x = barycenter[0] - facebarycenter[0];
  y = barycenter[1] - facebarycenter[1];
  z = barycenter[2] - facebarycenter[2];
  
  scalar =  n1*x + n2*y + n3*z;
  
  if ( scalar > 0 ) return 1;
  else if ( scalar < 0 ) return -1;
  else return 0;
}

int IsJointDOF(int DOF, int N_JointDOF, int *JointDOF)
{
  for (int i=0;i<N_JointDOF;++i)
  {
    if ( DOF == JointDOF[i] )
      return 1;
  }
  
  return 0;
}
