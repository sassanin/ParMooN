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
   
// =======================================================================
// @(#)NSE_MGLevel5.C        1.11 07/03/00
//
// Class:       TNSE_MGLevel5
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              type 5 (A)
//
// Author:      Gunar Matthies 28.10.2003
//
// History:     28.10.2003 start of implementation
//
// =======================================================================
#include <NSE_MGLevel5.h>
#include <Database.h>
#include <Constants.h>

#ifdef __2D__
  #include <FESpace2D.h>
  #include <FEDatabase2D.h>
#endif  
#ifdef __3D__
  #include <FESpace3D.h>
  #include <FEDatabase3D.h>
#endif  

#include <LinAlg.h>
#include <Solver.h>
#include <ItMethod.h>
#include <FgmresIte.h>

#include <Joint.h>

#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

/*
extern "C" {
#include "umfpack.h"
}
*/


/** constructor */
#ifdef __2D__
  TNSE_MGLevel5::TNSE_MGLevel5(int level, TSquareMatrixNSE2D *a, 
                               double *f1, double *u1, 
                               int n_aux, double *al, 
                               int velocity_space, int pressure_space,
                               TFESpace2D *pspace, TCollection *Coll)
#endif    

#ifdef __3D__
  TNSE_MGLevel5::TNSE_MGLevel5(int level, TSquareMatrixNSE3D *a, 
                               double *f1, double *u1, 
                               int n_aux, double *al, 
                               int velocity_space, int pressure_space,
                               TFESpace3D *pspace, TCollection *Coll)
#endif    
  : TNSE_MGLevel(level, f1, u1, n_aux, al,
                 velocity_space, pressure_space, Coll)
{
  int i;
  double *aux;

  A = a;
#ifdef __2D__
  StructureA = (TStructureNSE2D *)(A->GetMatrixStructure());
#endif

#ifdef __3D__
  StructureA = (TStructureNSE3D *)(A->GetMatrixStructure());
#endif
  ARowPtr = A->GetRowPtr();
  AKCol = A->GetKCol();
  AEntries = A->GetEntries();

  USpace = A->GetFESpace();

#ifdef __2D__
  PSpace = pspace;
#endif

#ifdef __3D__
  PSpace = pspace;
#endif

  if(VankaColl == NULL)
    VankaColl = USpace->GetCollection();
  
  N_Cells = VankaColl->GetN_Cells();

  N_Active = USpace->GetActiveBound();
  HangingNodeBound = USpace->GetHangingBound();
  N_Dirichlet = USpace->GetN_Dirichlet();

  N_UDOF = USpace->GetN_DegreesOfFreedom();
  N_PDOF = N_Cells;

  N_DOF = GEO_DIM*N_UDOF+N_PDOF;

  U1 = u1;
  U2 = u1 + N_UDOF;
  P  = u1 + GEO_DIM*N_UDOF;

  Rhs1 = f1;
  Rhs2 = f1 + N_UDOF;
  RhsP = f1 + GEO_DIM*N_UDOF;

#ifdef __3D__
  U3 = u1 + 2*N_UDOF;
  Rhs3 = f1 + 2*N_UDOF;
#endif  

  N_Aux = n_aux;
  Aux = new double* [N_Aux]; 
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();

  Type = 5;
}

/** destructor */
TNSE_MGLevel5::~TNSE_MGLevel5()
{
} // ~TNSE_MGLevel5

/** calculate defect */
void TNSE_MGLevel5::Defect(double *u, double *f, double *d, 
                           double &res)
{
  // compute defect
  CoupledDefect(A,u,f,d);

  // project defect into L20
#ifdef __2D__
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(d+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
#ifdef __3D__
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(d+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  

  memset(d+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(d+N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __3D__
  memset(d+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#endif  
  
  // compute residual
  res = sqrt(Ddot(N_DOF,d,d));
}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel5::CorrectNodes(double *u1)
{
  int i,j,k, index;
  double s, t, u, *u2, *u3;
  
  u2 = u1+N_UDOF;
#ifdef __3D__
  u3 = u2+N_UDOF;
#endif  
  
  // set Dirichlet nodes
  memset(u1+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
  memset(u2+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#ifdef __3D__
  memset(u3+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#endif  
}

#ifdef __2D__
/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel5::CellVanka(double *u1, double *rhs1, 
                              double *aux, 
                              int N_Parameters, double *Parameters, 
                              int smoother,int N_Levels)
{
  const int MaxDim = 5, LDB = 1, N_Rhs = 1;
  int LocalDim;
  int i,j,k,l,m,n;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int N_Joints, N_JointsNeigh;
  int DOFtype[MaxN_BaseFunctions2D];
  FE2D UElement, UElementNeigh;
  TFE2D *ele, *eleneigh;
  TFEDesc2D *FEDesc_Obj, *FEDesc_ObjNeigh;
  int N_Inner, *InnerDOF, N_LocalDOF;
  int *GlobalNumbers, *BeginIndex, *DOF, *DOFNeigh;
  int **JointDOFs, *EdgeDOF;
  int *EdgeDOFNeigh;
  int *LocJb; 
  int CompNumber[MaxDim];
  double *LocAlpha, bentry, aentry;
  double LocalSystem[MaxDim*MaxDim], LocalRhs[MaxDim];
  int LocalDOF[MaxDim];
  int begin, end, N_Eqn;
  double omega;
  double *u2, *p, *rhs2, *rhs3;
  double *uold, *pold;
  double r, s;
  double s11, s12, s21, s22, r1, r2;
  double m11, m12, m21, m22;
  // double tau = TDatabase::ParamDB->P8;
  double tau = 1;

  switch(N_Parameters)
  {
    case 1:
      omega = Parameters[0];
    break;
    case 2:
      if (Level==N_Levels-1)
        omega = Parameters[1];
      else
        omega = Parameters[0];
    break;
    default:
      omega = 0.7;
  }

  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  u2 = u1 + N_UDOF;
  p  = u2 + N_UDOF;

  rhs2 = rhs1 + N_UDOF;
  rhs3 = rhs2 + N_UDOF;

  // set Dirichlet values
  memcpy(u1+N_Active, rhs1+N_Active, N_Dirichlet*SizeOfDouble);
  memcpy(u2+N_Active, rhs2+N_Active, N_Dirichlet*SizeOfDouble);

  // save old values
  uold = aux;
  pold = aux + 2*N_UDOF;
  memcpy(uold, u1, N_DOF*SizeOfDouble);

  for(i=0;i<N_Cells;i++)
  {
    cell = VankaColl->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(UElement);
    FEDesc_Obj = ele->GetFEDesc2D();

    N_Inner = FEDesc_Obj->GetN_InnerDOF();
    InnerDOF = FEDesc_Obj->GetInnerDOF();

    N_LocalDOF = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    JointDOFs = FEDesc_Obj->GetJointDOF();

    // reset array of DOF types:
    //  0: usual node
    // -1: bubble node
    //  1: special edge dof in first component
    //  2: special edge dof in second component
    memset(DOFtype, 0, N_LocalDOF*SizeOfDouble);
    for(j=0;j<N_Inner;j++)
      DOFtype[InnerDOF[j]] = -1;

    LocJb = jb + BeginJb[i];

    LocalDim = N_Joints+1;

    // reset local system
    memset(LocalSystem, 0, LocalDim*LocalDim*SizeOfDouble);
    memset(LocalRhs, 0, LocalDim*SizeOfDouble);

    // number of pressure degree of freedom (dummy use)
    LocalDOF[N_Joints] = i;
    LocalRhs[N_Joints] += rhs3[i];

    // collect information from all edges
    for(j=0;j<N_Joints;j++)
    {
      EdgeDOF = JointDOFs[j];
      m = LocJb[j];

      if(m<N_LocalDOF)
      {
        DOFtype[m] = 1;
        CompNumber[j] = 1;
        l = DOF[m];
        LocalRhs[j] += rhs1[l];
      }
      else
      {
        m -= N_LocalDOF;
        DOFtype[m] = 2;
        CompNumber[j] = 2;
        l = DOF[m];
        LocalRhs[j] += rhs2[l];
      }
      LocalDOF[j] = l;

      // find edge local number of dof m
      for(k=0;k<N_DOFperJoint;k++)
        if(EdgeDOF[k] == m) break;

      LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
      if(CompNumber[j] == 1)
      {
        bentry = -LocAlpha[k];
        LocalRhs[N_Joints] -= bentry*u1[l];
      }
      else
      {
        bentry = -LocAlpha[k+N_DOFperJoint];
        LocalRhs[N_Joints] -= bentry*u2[l];
      }

      LocalSystem[N_Joints+LocalDim*j] = bentry;

      if(l < N_Active)
      {
        LocalSystem[j+LocalDim*N_Joints] = bentry;
        LocalRhs[j] -= bentry*p[i];

        // get B entry from neighbour for right-hand side
        // this the negative value of the own B entry
        joint = cell->GetJoint(j);
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          LocalRhs[j] += bentry*p[neigh->GetClipBoard()];
        }
      }
    } // endfor j (edge loop)

    // collect matrix entries for velocity dofs
    for(j=0;j<N_Joints;j++)
    {
      n = LocalDOF[j];

      if(n<N_Active)
      {
        // Active dof
        begin = ARowPtr[n]; 
        end = ARowPtr[n+1];
        for(k=begin;k<end;k++)
        {
          l = AKCol[k];
          if(CompNumber[j] == 1)
          {
            s11 = AEntries[4*k + 0];
            s12 = AEntries[4*k + 1];
            LocalRhs[j] -= s11*u1[l] + s12*u2[l];
          }
          else
          {
            s21 = AEntries[4*k + 2];
            s22 = AEntries[4*k + 3];
            LocalRhs[j] -= s21*u1[l] + s22*u2[l];
          }

          for(m=0;m<N_Joints;m++)
            if(LocalDOF[m] == l) break;

          if(m<N_Joints)
          {
            aentry = AEntries[4*k + (CompNumber[j]-1)*2 + (CompNumber[m]-1)];
            LocalSystem[j+LocalDim*m] += aentry;
          }
        } // endfor k
      }
      else
      {
        // Dirichlet node
        LocalSystem[j+LocalDim*j] = 1;
        LocalRhs[j] = 0;
      }
    } // endfor j

    N_Eqn = N_Joints + 1;

    /*
    // print local system
    for(j=0;j<N_Eqn;j++)
    {
      for(k=0;k<N_Eqn;k++)
        cout << setw(15) << LocalSystem[j+LocalDim*k];
      cout << setw(20) << LocalRhs[j] << endl;
    }
    cout << endl;
    */

    // solve linear system
    SolveLinearSystemLapack((double *)LocalSystem, (double *)LocalRhs, N_Eqn, LocalDim);
    
    /*
    // print solution
    for(j=0;j<N_Eqn;j++)
      cout << i << " " << LocalRhs[j] << endl;
    */ 

    // update solution
    for(j=0;j<N_Joints;j++)
    {
      if(CompNumber[j] == 1)
        u1[LocalDOF[j]] += tau * LocalRhs[j];
      else
        u2[LocalDOF[j]] += tau * LocalRhs[j];
    }
    p[i] += tau * LocalRhs[N_Joints]; 

    for(j=0;j<N_LocalDOF;j++)
    {
      switch(DOFtype[j])
      {
        case 0:
          // solve for both components at once
          m = DOF[j];
    
          if(m<N_Active)
          {
            // Active dof
            begin = ARowPtr[m]; 
            end = ARowPtr[m+1];
            r1 = rhs1[m];
            r2 = rhs2[m];
            m11 = 0;
            m12 = 0;
            m21 = 0;
            m22 = 0;
            for(k=begin;k<end;k++)
            {
              l = AKCol[k];
              s11 = AEntries[4*k + 0];
              s12 = AEntries[4*k + 1];
              s21 = AEntries[4*k + 2];
              s22 = AEntries[4*k + 3];

              if(m == l)
              {
                m11 += s11;
                m12 += s12;
                m21 += s21;
                m22 += s22;
              }
              r1 -= s11*u1[l] + s12*u2[l];
              r2 -= s21*u1[l] + s22*u2[l];
            }
            s = 1/(m11*m22 - m12*m21);
            // cout << 0 << " " << m11 << " " << m12 << " " << r1 << endl;
            // cout << 0 << " " << m21 << " " << m22 << " " << r2 << endl;
            u1[m] += tau * (m22*r1-r2*m12)*s;
            u2[m] += tau * (m11*r2-r1*m21)*s;
          }
        break;

        case 1:
          // only second component, first is special edge dof
          m = DOF[j];
    
          if(m<N_Active)
          {
            // Active dof
            begin = ARowPtr[m]; 
            end = ARowPtr[m+1];
            r = rhs2[m];
            s = 0;
            for(k=begin;k<end;k++)
            {
              l = AKCol[k];
              s21 = AEntries[4*k + 2];
              s22 = AEntries[4*k + 3];
              r -= s21*u1[l] + s22*u2[l];

              if(m == l)
                s += s22;
            }
            // cout << 1 << " " << s << " " << r << endl;
            u2[m] += tau * r/s;
          }
        break;

        case 2:
          // only first component, second is special edge dof
          m = DOF[j];
    
          if(m<N_Active)
          {
            // Active dof
            begin = ARowPtr[m]; 
            end = ARowPtr[m+1];
            r = rhs1[m];
            s = 0;
            for(k=begin;k<end;k++)
            {
              l = AKCol[k];
              s11 = AEntries[4*k + 0];
              s12 = AEntries[4*k + 1];
              r -= s11*u1[l] + s12*u2[l];

              if(m == l)
                s += s11;
            }
            // cout << 2 << " " << s << " " << r << endl;
            u1[m] += tau * r/s;
          }
        break;

        case -1:
          // do nothing
        break;
      } // end switch DOFtype
    } // endfor j (dof loop)
  } // endfor i (cell loop)

  // apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] = uold[j]+omega*(u1[j]-uold[j]);

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound, 
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound, 
         N_Dirichlet*SizeOfDouble);

} // end Vanka
#endif // __2D__

#ifdef __3D__
/** cellwise Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel5::CellVanka(double *u1, double *rhs1, 
                              double *aux, 
                              int N_Parameters, double *Parameters, 
                              int smoother,int N_Levels)
{
  const int MaxDim = 7, LDB = 1, N_Rhs = 1;
  int LocalDim;
  int i,j,k,l,m,n;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int N_Joints, N_JointsNeigh;
  int DOFtype[MaxN_BaseFunctions3D];
  FE3D UElement, UElementNeigh;
  TFE3D *ele, *eleneigh;
  TFEDesc3D *FEDesc_Obj, *FEDesc_ObjNeigh;
  int N_Inner, *InnerDOF, N_LocalDOF;
  int *GlobalNumbers, *BeginIndex, *DOF, *DOFNeigh;
  int **JointDOFs, *EdgeDOF;
  int *EdgeDOFNeigh;
  int *LocJb; 
  int CompNumber[MaxDim];
  double *LocAlpha, bentry, aentry;
  double LocalSystem[MaxDim*MaxDim], LocalRhs[MaxDim];
  int LocalDOF[MaxDim];
  int begin, end, N_Eqn;
  double omega;
  double *u2, *u3, *p, *rhs2, *rhs3, *rhs4;
  double *uold, *pold;
  double r, s, t;
  double s11, s12, s13, s21, s22, s23, s31, s32, s33, r1, r2, r3;
  double m11, m12, m13, m21, m22, m23, m31, m32, m33;
  // double tau = TDatabase::ParamDB->P8;
  double tau = 1;
  double System3[3][3], System2[2][2], Rhs3[3], Rhs2[2];

  switch(N_Parameters)
  {
    case 1:
      omega = Parameters[0];
    break;
    case 2:
      if (Level==N_Levels-1)
        omega = Parameters[1];
      else
        omega = Parameters[0];
    break;
    default:
      omega = 0.7;
  }

  omega = 0.7;

  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  u2 = u1 + N_UDOF;
  u3 = u2 + N_UDOF;
  p  = u3 + N_UDOF;

  rhs2 = rhs1 + N_UDOF;
  rhs3 = rhs2 + N_UDOF;
  rhs4 = rhs3 + N_UDOF;

  // set Dirichlet values
  memcpy(u1+N_Active, rhs1+N_Active, N_Dirichlet*SizeOfDouble);
  memcpy(u2+N_Active, rhs2+N_Active, N_Dirichlet*SizeOfDouble);
  memcpy(u3+N_Active, rhs3+N_Active, N_Dirichlet*SizeOfDouble);

  // save old values
  uold = aux;
  pold = aux + 3*N_UDOF;
  memcpy(uold, u1, N_DOF*SizeOfDouble);

  for(i=0;i<N_Cells;i++)
  {
    cell = VankaColl->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(UElement);
    FEDesc_Obj = ele->GetFEDesc3D();

    N_Inner = FEDesc_Obj->GetN_InnerDOF();
    InnerDOF = FEDesc_Obj->GetInnerDOF();

    N_LocalDOF = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    JointDOFs = FEDesc_Obj->GetJointDOF();

    // reset array of DOF types:
    //  0: usual node
    // -1: bubble node
    //  1: special edge dof in first component
    //  2: special edge dof in second component
    //  3: special edge dof in second component
    memset(DOFtype, 0, N_LocalDOF*SizeOfDouble);
    for(j=0;j<N_Inner;j++)
      DOFtype[InnerDOF[j]] = -1;

    LocJb = jb + BeginJb[i];

    LocalDim = N_Joints+1;

    // reset local system
    memset(LocalSystem, 0, LocalDim*LocalDim*SizeOfDouble);
    memset(LocalRhs, 0, LocalDim*SizeOfDouble);

    // number of pressure degree of freedom (dummy use)
    LocalDOF[N_Joints] = i;
    LocalRhs[N_Joints] += rhs4[i];

    // collect information from all edges
    for(j=0;j<N_Joints;j++)
    {
      EdgeDOF = JointDOFs[j];
      m = LocJb[j];

      if(m<N_LocalDOF)
      {
        DOFtype[m] = 1;
        CompNumber[j] = 1;
        l = DOF[m];
        LocalRhs[j] += rhs1[l];
      }
      else
      {
        if(m<2*N_LocalDOF)
        {
          m -= N_LocalDOF;
          DOFtype[m] = 2;
          CompNumber[j] = 2;
          l = DOF[m];
          LocalRhs[j] += rhs2[l];
        }
        else
        {
          m -= 2*N_LocalDOF;
          DOFtype[m] = 3;
          CompNumber[j] = 3;
          l = DOF[m];
          LocalRhs[j] += rhs3[l];
        }
      }
      LocalDOF[j] = l;

      // find edge local number of dof m
      for(k=0;k<N_DOFperJoint;k++)
        if(EdgeDOF[k] == m) break;

      LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
      switch(CompNumber[j])
      {
        case 1:
          bentry = -LocAlpha[k];
          LocalRhs[N_Joints] -= bentry*u1[l];
        break;

        case 2:
          bentry = -LocAlpha[k+N_DOFperJoint];
          LocalRhs[N_Joints] -= bentry*u2[l];
        break;

        case 3:
          bentry = -LocAlpha[k+2*N_DOFperJoint];
          LocalRhs[N_Joints] -= bentry*u3[l];
        break;
      }

      LocalSystem[N_Joints+LocalDim*j] = bentry;

      if(l < N_Active)
      {
        LocalSystem[j+LocalDim*N_Joints] = bentry;
        LocalRhs[j] -= bentry*p[i];

        // get B entry from neighbour for right-hand side
        // this the negative value of the own B entry
        joint = cell->GetJoint(j);
        neigh = joint->GetNeighbour(cell);
        if(neigh)
        {
          LocalRhs[j] += bentry*p[neigh->GetClipBoard()];
        }
      }
    } // endfor j (edge loop)

    // collect matrix entries for velocity dofs
    for(j=0;j<N_Joints;j++)
    {
      n = LocalDOF[j];

      if(n<N_Active)
      {
        // Active dof
        begin = ARowPtr[n]; 
        end = ARowPtr[n+1];
        for(k=begin;k<end;k++)
        {
          l = AKCol[k];
          switch(CompNumber[j])
          {
            case 1:
              s11 = AEntries[9*k + 0];
              s12 = AEntries[9*k + 1];
              s13 = AEntries[9*k + 2];
              LocalRhs[j] -= s11*u1[l] + s12*u2[l] + s13*u3[l];
            break;
            
            case 2:
              s21 = AEntries[9*k + 3];
              s22 = AEntries[9*k + 4];
              s23 = AEntries[9*k + 5];
              LocalRhs[j] -= s21*u1[l] + s22*u2[l] + s23*u3[l];
            break;
            
            case 3:
              s31 = AEntries[9*k + 6];
              s32 = AEntries[9*k + 7];
              s33 = AEntries[9*k + 8];
              LocalRhs[j] -= s31*u1[l] + s32*u2[l] + s33*u3[l];
            break;
          }

          for(m=0;m<N_Joints;m++)
            if(LocalDOF[m] == l) break;

          if(m<N_Joints)
          {
            aentry = AEntries[9*k + (CompNumber[j]-1)*3 + (CompNumber[m]-1)];
            LocalSystem[j+LocalDim*m] += aentry;
          }
        } // endfor k
      }
      else
      {
        // Dirichlet node
        LocalSystem[j+LocalDim*j] = 1;
        LocalRhs[j] = 0;
      }
    } // endfor j

    N_Eqn = N_Joints + 1;

    /*
    // print local system
    for(j=0;j<N_Eqn;j++)
    {
      for(k=0;k<N_Eqn;k++)
        cout << setw(15) << LocalSystem[j+LocalDim*k];
      cout << setw(20) << LocalRhs[j] << endl;
    }
    cout << endl;
    */

    // solve linear system
    SolveLinearSystemLapack((double *)LocalSystem, (double *)LocalRhs, N_Eqn, LocalDim);
    
    /*
    // print solution
    for(j=0;j<N_Eqn;j++)
      cout << i << " " << LocalRhs[j] << endl;
    */ 

    // update solution
    for(j=0;j<N_Joints;j++)
    {
      switch(CompNumber[j])
      {
        case 1:
          u1[LocalDOF[j]] += tau * LocalRhs[j];
        break;

        case 2:
          u2[LocalDOF[j]] += tau * LocalRhs[j];
        break;

        case 3:
          u3[LocalDOF[j]] += tau * LocalRhs[j];
        break;
      }
    }
    p[i] += tau * LocalRhs[N_Joints]; 

    for(j=0;j<N_LocalDOF;j++)
    {
      m = DOF[j];
    
      if(m<N_Active)
      {
        // Active dof
        begin = ARowPtr[m]; 
        end = ARowPtr[m+1];
        r1 = rhs1[m];
        r2 = rhs2[m];
        r3 = rhs3[m];
        m11 = 0;
        m12 = 0;
        m13 = 0;
        m21 = 0;
        m22 = 0;
        m23 = 0;
        m31 = 0;
        m32 = 0;
        m33 = 0;
        for(k=begin;k<end;k++)
        {
          l = AKCol[k];
          s11 = AEntries[9*k + 0];
          s12 = AEntries[9*k + 1];
          s13 = AEntries[9*k + 2];
          s21 = AEntries[9*k + 3];
          s22 = AEntries[9*k + 4];
          s23 = AEntries[9*k + 5];
          s31 = AEntries[9*k + 6];
          s32 = AEntries[9*k + 7];
          s33 = AEntries[9*k + 8];

          if(m == l)
          {
            m11 += s11;
            m12 += s12;
            m13 += s13;
            m21 += s21;
            m22 += s22;
            m23 += s23;
            m31 += s31;
            m32 += s32;
            m33 += s33;
          }
          r1 -= s11*u1[l] + s12*u2[l] + s13*u3[l];
          r2 -= s21*u1[l] + s22*u2[l] + s23*u3[l];
          r3 -= s31*u1[l] + s32*u2[l] + s33*u3[l];
        }
        /*
        if(DOFtype[j] != -1)
        {
          if(DOFtype[j] != 1) 
            u1[m] += tau * r1/m11;
          if(DOFtype[j] != 2) 
            u2[m] += tau * r2/m22;
          if(DOFtype[j] != 3) 
            u3[m] += tau * r3/m33;
        }
        */
        switch(DOFtype[j])
        {
          case 0:
            // solve for all three components
            System3[0][0] = m11;
            System3[0][1] = m12;
            System3[0][2] = m13;
            System3[1][0] = m21;
            System3[1][1] = m22;
            System3[1][2] = m23;
            System3[2][0] = m31;
            System3[2][1] = m32;
            System3[2][2] = m33;

            Rhs3[0] = r1;
            Rhs3[1] = r2;
            Rhs3[2] = r3;
    
            N_Eqn = 3;
            LocalDim = 3;
            SolveLinearSystemLapack((double *)System3, (double *)Rhs3, N_Eqn, LocalDim);

            u1[m] += tau * Rhs3[0];
            u2[m] += tau * Rhs3[1];
            u3[m] += tau * Rhs3[2];
          break;

          case 1:
            // solve for second and third component
            System2[0][0] = m22;
            System2[0][1] = m23;
            System2[1][0] = m32;
            System2[1][1] = m33;

            Rhs2[0] = r2;
            Rhs2[1] = r3;

            N_Eqn = 2;
            LocalDim = 2;
            SolveLinearSystemLapack((double *)System2, (double *)Rhs2, N_Eqn, LocalDim);

            u2[m] += tau * Rhs2[0];
            u3[m] += tau * Rhs2[1];
          break;

          case 2:
            // solve for first and third component
            System2[0][0] = m11;
            System2[0][1] = m13;
            System2[1][0] = m31;
            System2[1][1] = m33;

            Rhs2[0] = r1;
            Rhs2[1] = r3;

            N_Eqn = 2;
            LocalDim = 2;
            SolveLinearSystemLapack((double *)System2, (double *)Rhs2, N_Eqn, LocalDim);

            u1[m] += tau * Rhs2[0];
            u3[m] += tau * Rhs2[1];
          break;

          case 3:
            // solve for first and second component
            System2[0][0] = m11;
            System2[0][1] = m12;
            System2[1][0] = m21;
            System2[1][1] = m22;

            Rhs2[0] = r1;
            Rhs2[1] = r2;

            N_Eqn = 2;
            LocalDim = 2;
            SolveLinearSystemLapack((double *)System2, (double *)Rhs2, N_Eqn, LocalDim);

            u1[m] += tau * Rhs2[0];
            u2[m] += tau * Rhs2[1];
          break;
        }
      }
    } // endfor j (dof loop)
  } // endfor i (cell loop)

  // apply damping
  for(j=0;j<N_DOF;j++)
    u1[j] = uold[j]+omega*(u1[j]-uold[j]);

  // set Dirichlet values
  memcpy(u1+HangingNodeBound, rhs1+HangingNodeBound, 
         N_Dirichlet*SizeOfDouble);
  memcpy(u2+HangingNodeBound, rhs2+HangingNodeBound, 
         N_Dirichlet*SizeOfDouble);
  memcpy(u3+HangingNodeBound, rhs3+HangingNodeBound, 
         N_Dirichlet*SizeOfDouble);
}
#endif // __3D__

/** nodal Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel5::NodalVanka(double *u1, double *rhs1, double *aux, 
        int N_Parameters, double *Parameters, int smoother,int N_Levels)
{
} // end Vanka

#ifdef __2D__
/** solve exact on this level, old version */
void TNSE_MGLevel5::SolveExact(double *u1, double *rhs1)
{  
  int i,j,k,l,m;
  int N_Cells, N_Joints;
  int N_Rows, N_Eqn, LDA, LDB, N_Rhs;
  int NP, DirichletBound, N_U;
  double *System, *Rhs;
  double *ConstantP;
  TCollection *Coll;
  TBaseCell *cell;
  int *ColInd, *RowPtr, *DOF, *GlobalNumbers, *BeginIndex, end;
  double *Entries;
  FE2D UElement;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  int **JointDOFs, *EdgeDOF, *InnerDOF;
  int *LocJb;
  double *LocAlpha, val;
  int N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC, *BeginP;
  double *C, *Pcoeff;

  Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  N_Rows = USpace->GetN_DegreesOfFreedom();
  N_Eqn = 2*N_Rows+N_Cells;
  LDA = N_Eqn;
  LDB = 1;

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();

  BeginC = A->GetBeginC();
  C = A->GetC();

  BeginP = A->GetBeginP();
  Pcoeff = A->GetP();

  Entries = A->GetEntries();
  RowPtr = A->GetRowPtr();
  ColInd = A->GetKCol();

  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  System = new double[N_Eqn*N_Eqn];
  Rhs = new double[N_Eqn];

  // fill rhs
  memcpy(Rhs, rhs1, N_Eqn*SizeOfDouble);

  // reset system
  memset(System, 0, N_Eqn*N_Eqn*SizeOfDouble);

  NP = PSpace->GetN_DegreesOfFreedom();


  DirichletBound = USpace->GetHangingBound();
  // fill A-blocks
  for(i=0;i<DirichletBound;i++)
  {
    end=RowPtr[i+1];
    for(j=RowPtr[i];j<end;j++)
    {
      k = ColInd[j];
      System[ i        *LDA + k       ] = Entries[4*j+0]; // A11
      System[ i        *LDA + k+N_Rows] = Entries[4*j+1]; // A12
      System[(i+N_Rows)*LDA + k       ] = Entries[4*j+2]; // A21
      System[(i+N_Rows)*LDA + k+N_Rows] = Entries[4*j+3]; // A22
    }
  }

  // fill B-blocks and add entries for dummy bubbles
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(UElement);
    FEDesc_Obj = ele->GetFEDesc2D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // cout << endl;
    // cout << "B entries" << endl;
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        System[ DOF[l]     *LDA + (2*N_Rows+i)] = -LocAlpha[k];
        System[(2*N_Rows+i)*LDA + DOF[l]      ] = -LocAlpha[k];
      }
      else
      {
        // jb in second component
        l -= N_U;
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        System[(DOF[l]+N_Rows)*LDA + (2*N_Rows+i)   ] =
                -LocAlpha[k+N_DOFperJoint];
        System[(2*N_Rows+i)   *LDA + (DOF[l]+N_Rows)] =
                -LocAlpha[k+N_DOFperJoint];
      }
    } // endfor j

    // internal bubbles
    N_Inner = FEDesc_Obj->GetN_InnerDOF();
    InnerDOF = FEDesc_Obj->GetInnerDOF();

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      System[k*LDA+k] = 1;
      Rhs[k] = 0;

      k += N_Rows;
      System[k*LDA+k] = 1;
      Rhs[k] = 0;
    } // endfor j
  } // endfor i

  // realise Dirichlet coundary conditions
  m=RowPtr[DirichletBound];
  for(i=DirichletBound;i<N_Rows;i++)
  {
    memset(System+i*LDA, 0, LDA*SizeOfDouble);
    System[ i        *LDA + i       ] = 1;

    memset(System+(i+N_Rows)*LDA, 0, LDA*SizeOfDouble);
    System[(i+N_Rows)*LDA + i+N_Rows] = 1;
  }

  /*
  // print whole system
  cout << setprecision(20);
  for(i=0;i<N_Eqn;i++)
  {
    for(j=0;j<N_Eqn;j++)
    {
      cout << setw(5) << i+1 << setw(5) << j+1 << setw(30);
      cout << System[i*LDA+j] << endl;
    } // endfor j
  } // endfor i

  // write rhs
  for(i=0;i<N_Eqn;i++)
  {
    cout << setw(5) << i+1 << setw(30) << Rhs[i] << endl;
  }
  */

  // set pressure constant
  memset(System+(N_Eqn-1)*LDA, 0, LDA*SizeOfDouble);
  for(j=0;j<NP;j++)
    System[(N_Eqn-1)*LDA+2*N_Rows+j] = 1;
  Rhs[N_Eqn-1] = 0;

  // solve
  N_Rhs = 1;
  SolveMultipleSystemsNew((double *)System, (double *)Rhs, N_Eqn,
                          LDA, LDB, N_Rhs);

  /*
  cout << "velocity" << endl;
  for(i=0;i<N_Rows;i++)
  {
    if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
    if(fabs(Rhs[i+N_Rows]) < 1e-10) Rhs[i+N_Rows] = 0;
    cout << setw(5) << i << setw(15) << Rhs[i];
    cout << setw(15) << Rhs[i+N_Rows] << endl;
  }
  cout << endl;
  
  cout << "pressure" << endl;
  for(i=2*N_Rows;i<N_Eqn;i++)
  {
    if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
    cout << setw(5) << Rhs[i] << endl;
  }
  cout << endl;
  */

  memcpy(u1, Rhs, N_Eqn*SizeOfDouble);
  delete System;
  delete Rhs;
}

// /** solve exact on this level */
// void TNSE_MGLevel5::SolveExact(double *u1, double *rhs1)
// {  
//   int i,j,k,l,m;
//   int N_Cells, N_Joints;
//   int N_Rows, N_Eqn, LDA, LDB, N_Rhs;
//   int NP, DirichletBound, N_U;
//   double *Rhs;
//   double *ConstantP;
//   TCollection *Coll;
//   TBaseCell *cell;
//   int *ColInd, *RowPtr, *DOF, *GlobalNumbers, *BeginIndex;
//   double *Entries;
//   FE2D UElement;
//   TFE2D *ele;
//   TFEDesc2D *FEDesc_Obj;
//   int **JointDOFs, *EdgeDOF, *InnerDOF;
//   int *LocJb;
//   double *LocAlpha, val;
//   int N_Inner;
//   int *BeginJb, *jb, N_DOFperJoint;
//   double *Alpha;
//   int *BeginC, *BeginP;
//   double *C, *Pcoeff;
// 
//   // new date
//   int *KCol, *Row, *Current;
//   double *Values, value;
//   double *null = (double *) NULL;
//   int begin, end;
//   double t1, t2, t3, t4;
//   void *Symbolic, *Numeric;
//   double *SOL;
// 
//   Coll = USpace->GetCollection();
//   N_Cells = Coll->GetN_Cells();
//   
//   N_Rows = USpace->GetN_DegreesOfFreedom();
//   N_Eqn = 2*N_Rows+N_Cells;
//   LDA = N_Eqn;
//   LDB = 1;
// 
//   BeginJb = A->GetBeginJb();
//   jb = A->GetJb();
//   N_DOFperJoint = A->GetN_DOFperJoint();
//   Alpha = A->GetAlpha();
// 
//   BeginC = A->GetBeginC();
//   C = A->GetC();
// 
//   BeginP = A->GetBeginP();
//   Pcoeff = A->GetP();
// 
//   Entries = A->GetEntries();
//   RowPtr = A->GetRowPtr();
//   ColInd = A->GetKCol();
// 
//   GlobalNumbers = USpace->GetGlobalNumbers();
//   BeginIndex = USpace->GetBeginIndex();
// 
//   Row = new int[N_Eqn+1];
//   memset(Row, 0, (N_Eqn+1)*SizeOfInt);
// 
//   Rhs = new double[N_Eqn];
// 
//   // fill rhs
//   memcpy(Rhs, rhs1, N_Eqn*SizeOfDouble);
// 
//   NP = PSpace->GetN_DegreesOfFreedom();
// 
//   DirichletBound = USpace->GetHangingBound();
// 
//   // fill A-blocks
//   for(i=0;i<DirichletBound;i++)
//   {
//     Row[i] += 2*(RowPtr[i+1]-RowPtr[i]);
//     Row[i+N_Rows] += 2*(RowPtr[i+1]-RowPtr[i]);
//   }
// 
//   // fill B-blocks and add entries for dummy bubbles
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     N_Joints = cell->GetN_Joints();
// 
//     UElement = USpace->GetFE2D(i, cell);
//     ele = TFEDatabase2D::GetFE2D(UElement);
//     FEDesc_Obj = ele->GetFEDesc2D();
//     JointDOFs = FEDesc_Obj->GetJointDOF();
//     N_U = FEDesc_Obj->GetN_DOF();
// 
//     DOF = GlobalNumbers + BeginIndex[i];
// 
//     //  B entries
//     LocJb = jb+BeginJb[i];
// 
//     // cout << endl;
//     // cout << "B entries" << endl;
//     for(j=0;j<N_Joints;j++)
//     {
//       // get local data
//       EdgeDOF = JointDOFs[j]; 
//         
//       l = LocJb[j];
//       if(l<N_U)
//       {
//         // jb in first component
//         for(k=0;k<N_DOFperJoint;k++)
//           if(EdgeDOF[k] == l)
//             break;
// 
//         if(DOF[l]<DirichletBound) Row[DOF[l]]++;
//         Row[2*N_Rows+i]++;
//       }
//       else
//       {
//         // jb in second component
//         l -= N_U;
//         for(k=0;k<N_DOFperJoint;k++)
//           if(EdgeDOF[k] == l)
//             break;
// 
//         if(DOF[l]<DirichletBound) Row[DOF[l]+N_Rows]++;
//         Row[2*N_Rows+i]++;
//       }
//     } // endfor j
// 
//     // internal bubbles
//     N_Inner = FEDesc_Obj->GetN_InnerDOF();
//     InnerDOF = FEDesc_Obj->GetInnerDOF();
// 
//     for(j=0;j<N_Inner;j++)
//     {
//       k = DOF[InnerDOF[j]];
//       Row[k]++;
// 
//       k += N_Rows;
//       Row[k]++;
//     } // endfor j
//   } // endfor i
// 
//   // realise Dirichlet coundary conditions
//   m=RowPtr[DirichletBound];
//   for(i=DirichletBound;i<N_Rows;i++)
//   {
//     Row[i]=1;
//     Row[i+N_Rows]=1;
//   }
// 
//   if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
//   {
//     // for pressure constant
//     Row[N_Eqn-1] = NP;
//   }
// 
//   // for(i=0;i<N_Eqn;i++)
//   //   cout << setw(3) << i << setw(10) << Row[i] << endl;
// 
//   // sum up Row array
//   m = 0;
//   for(i=0;i<N_Eqn;i++)
//   {
//     l = m;
//     m += Row[i];
//     Row[i] = l;
//   }
//   Row[N_Eqn] = m;
// 
//   Current = new int[N_Eqn+1];
//   memcpy(Current, Row, (N_Eqn+1)*SizeOfInt);
//   Values = new double[m];
//   KCol = new int[m];
// 
//   // fill A-blocks
//   for(i=0;i<DirichletBound;i++)
//   {
//     end=RowPtr[i+1];
//     for(j=RowPtr[i];j<end;j++)
//     {
//       k = ColInd[j];
// 
//       Values[Current[i]] = Entries[4*j+0]; // A11
//       KCol[Current[i]] = k;
//       Current[i]++;
// 
//       Values[Current[i]] = Entries[4*j+1]; // A12
//       KCol[Current[i]] = k+N_Rows;
//       Current[i]++;
// 
//       Values[Current[i+N_Rows]] = Entries[4*j+2]; // A21
//       KCol[Current[i+N_Rows]] = k;
//       Current[i+N_Rows]++;
// 
//       Values[Current[i+N_Rows]] = Entries[4*j+3]; // A22
//       KCol[Current[i+N_Rows]] = k+N_Rows;
//       Current[i+N_Rows]++;
//     }
//   }
// 
//   // fill B-blocks and add entries for dummy bubbles
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     N_Joints = cell->GetN_Joints();
// 
//     UElement = USpace->GetFE2D(i, cell);
//     ele = TFEDatabase2D::GetFE2D(UElement);
//     FEDesc_Obj = ele->GetFEDesc2D();
//     JointDOFs = FEDesc_Obj->GetJointDOF();
//     N_U = FEDesc_Obj->GetN_DOF();
// 
//     DOF = GlobalNumbers + BeginIndex[i];
// 
//     //  B entries
//     LocJb = jb+BeginJb[i];
// 
//     // cout << endl;
//     // cout << "B entries" << endl;
//     for(j=0;j<N_Joints;j++)
//     {
//       // get local data
//       LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
//       EdgeDOF = JointDOFs[j]; 
//         
//       l = LocJb[j];
//       if(l<N_U)
//       {
//         // jb in first component
//         for(k=0;k<N_DOFperJoint;k++)
//           if(EdgeDOF[k] == l)
//             break;
// 
//         if(DOF[l]<DirichletBound)
//         {
//           Values[Current[DOF[l]]] = -LocAlpha[k];
//           KCol[Current[DOF[l]]] = 2*N_Rows+i;
//           Current[DOF[l]]++;
//         }
// 
//         Values[Current[2*N_Rows+i]] = -LocAlpha[k];
//         KCol[Current[2*N_Rows+i]] = DOF[l];
//         Current[2*N_Rows+i]++;
//       }
//       else
//       {
//         // jb in second component
//         l -= N_U;
//         for(k=0;k<N_DOFperJoint;k++)
//           if(EdgeDOF[k] == l)
//             break;
// 
//         if(DOF[l]<DirichletBound)
//         {
//           Values[Current[DOF[l]+N_Rows]] = -LocAlpha[k+N_DOFperJoint];
//           KCol[Current[DOF[l]+N_Rows]] = 2*N_Rows+i;
//           Current[DOF[l]+N_Rows]++;
//         }
// 
//         Values[Current[2*N_Rows+i]] = -LocAlpha[k+N_DOFperJoint];
//         KCol[Current[2*N_Rows+i]] = DOF[l]+N_Rows;
//         Current[2*N_Rows+i]++;
//       }
//     } // endfor j
// 
//     // internal bubbles
//     N_Inner = FEDesc_Obj->GetN_InnerDOF();
//     InnerDOF = FEDesc_Obj->GetInnerDOF();
// 
//     for(j=0;j<N_Inner;j++)
//     {
//       k = DOF[InnerDOF[j]];
//       Rhs[k] = 0;
//       Values[Current[k]] = 1;
//       KCol[Current[k]] = k;
//       Current[k]++;
// 
//       k += N_Rows;
//       Rhs[k] = 0;
//       Values[Current[k]] = 1;
//       KCol[Current[k]] = k;
//       Current[k]++;
//     } // endfor j
//   } // endfor i
// 
//   // realise Dirichlet coundary conditions
//   m=RowPtr[DirichletBound];
//   for(i=DirichletBound;i<N_Rows;i++)
//   {
//     Values[Current[i]] = 1;
//     KCol[Current[i]] = i;
//     Current[i]++;
// 
//     Values[Current[i+N_Rows]] = 1;
//     KCol[Current[i+N_Rows]] = i+N_Rows;
//     Current[i+N_Rows]++;
//   }
// 
//   if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
//   {
//     Current[N_Eqn-1] = Row[N_Eqn-1];
//     // set pressure constant
//     for(j=0;j<NP;j++)
//     {
//       Values[Current[N_Eqn-1]] = 1;
//       KCol[Current[N_Eqn-1]] = 2*N_Rows+j;
//       Current[N_Eqn-1]++;
//     }
//     Rhs[N_Eqn-1] = 0;
//   }
// 
//   SOL = new double[N_Eqn];
// 
//   // sort matrix
//   for(i=0;i<N_Eqn;i++)
//   {
//     begin=Row[i];
//     end=Row[i+1];
// 
//     for(j=begin;j<end;j++)
//     {
//       for(k=j+1;k<end;k++)
//       {
//         if(KCol[j] > KCol[k])
//         {
//                 l = KCol[j];     value = Values[j];
//           KCol[j] = KCol[k]; Values[j] = Values[k];
//           KCol[k] = l;       Values[k] = value;
//         } // endif
//       } // endfor k
//     } // endfor j
//   } // endfor i
// 
//   // t1 = GetTime();
//   i = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
//                           &Symbolic, null, null);
//   // t2 = GetTime();
//   // OutPut("symbolic: " << i << " " << t2-t1 << endl);
// 
// 
//   i = umfpack_di_numeric(Row, KCol, Values, Symbolic,
//                          &Numeric, null, null);
//   umfpack_di_free_symbolic(&Symbolic);
//   // t3 = GetTime();
//   // OutPut("numeric: " << i << " "  << t3-t2 << endl);
// 
//   i = umfpack_di_solve(UMFPACK_At, Row, KCol, Values,
//                        SOL, Rhs, Numeric, null, null);
//   umfpack_di_free_numeric(&Numeric);
//   // t4 = GetTime();
//   // OutPut("solve: " << i << " " << t4-t3 << endl);
// 
//   memcpy(Rhs, SOL, N_Eqn*SizeOfDouble);
// 
//   /*
//   cout << "velocity" << endl;
//   for(i=0;i<N_Rows;i++)
//   {
//     if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
//     if(fabs(Rhs[i+N_Rows]) < 1e-10) Rhs[i+N_Rows] = 0;
//     cout << setw(5) << i << setw(15) << Rhs[i];
//     cout << setw(15) << Rhs[i+N_Rows] << endl;
//   }
//   cout << endl;
//   
//   cout << "pressure" << endl;
//   for(i=2*N_Rows;i<N_Eqn;i++)
//   {
//     if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
//     cout << setw(5) << Rhs[i] << endl;
//   }
//   cout << endl;
//   */
// 
//   memcpy(u1, Rhs, N_Eqn*SizeOfDouble);
//   delete Rhs;
// 
//   delete Row;
//   delete KCol;
//   delete Values;
//   delete Current;
//   delete SOL;
// }
#endif

#ifdef __3D__
/** solve exact on this level, old version */
void TNSE_MGLevel5::SolveExact(double *u1, double *rhs1)
{  
  int i,j,k,l,m;
  int N_Cells, N_Joints;
  int N_Rows, N_Eqn, LDA, LDB, N_Rhs;
  int NP, DirichletBound, N_U;
  double *System, *Rhs;
  double *ConstantP;
  TCollection *Coll;
  TBaseCell *cell;
  int *ColInd, *RowPtr, *DOF, *GlobalNumbers, *BeginIndex, end;
  double *Entries;
  FE3D UElement;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  int **JointDOFs, *EdgeDOF, *InnerDOF;
  int *LocJb;
  double *LocAlpha, val;
  int N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC, *BeginP;
  double *C, *Pcoeff;

  Coll = USpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  
  N_Rows = USpace->GetN_DegreesOfFreedom();
  N_Eqn = 3*N_Rows+N_Cells;
  LDA = N_Eqn;
  LDB = 1;

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();

  BeginC = A->GetBeginC();
  C = A->GetC();

  BeginP = A->GetBeginP();
  Pcoeff = A->GetP();

  Entries = A->GetEntries();
  RowPtr = A->GetRowPtr();
  ColInd = A->GetKCol();

  GlobalNumbers = USpace->GetGlobalNumbers();
  BeginIndex = USpace->GetBeginIndex();

  System = new double[N_Eqn*N_Eqn];
  Rhs = new double[N_Eqn];

  // fill rhs
  memcpy(Rhs, rhs1, N_Eqn*SizeOfDouble);

  // reset system
  memset(System, 0, N_Eqn*N_Eqn*SizeOfDouble);

  NP = PSpace->GetN_DegreesOfFreedom();

  DirichletBound = USpace->GetHangingBound();
  // fill A-blocks
  for(i=0;i<DirichletBound;i++)
  {
    end=RowPtr[i+1];
    for(j=RowPtr[i];j<end;j++)
    {
      k = ColInd[j];
      System[ i          *LDA + k         ] = Entries[9*j+0]; // A11
      System[ i          *LDA + k+  N_Rows] = Entries[9*j+1]; // A12
      System[ i          *LDA + k+2*N_Rows] = Entries[9*j+2]; // A13
      System[(i+  N_Rows)*LDA + k         ] = Entries[9*j+3]; // A21
      System[(i+  N_Rows)*LDA + k+  N_Rows] = Entries[9*j+4]; // A22
      System[(i+  N_Rows)*LDA + k+2*N_Rows] = Entries[9*j+5]; // A23
      System[(i+2*N_Rows)*LDA + k         ] = Entries[9*j+6]; // A31
      System[(i+2*N_Rows)*LDA + k+  N_Rows] = Entries[9*j+7]; // A32
      System[(i+2*N_Rows)*LDA + k+2*N_Rows] = Entries[9*j+8]; // A33
    }
  }

  // fill B-blocks and add entries for dummy bubbles
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Joints = cell->GetN_Joints();

    UElement = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(UElement);
    FEDesc_Obj = ele->GetFEDesc3D();
    JointDOFs = FEDesc_Obj->GetJointDOF();
    N_U = FEDesc_Obj->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    //  B entries
    LocJb = jb+BeginJb[i];

    // cout << endl;
    // cout << "B entries" << endl;
    for(j=0;j<N_Joints;j++)
    {
      // get local data
      LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
      EdgeDOF = JointDOFs[j]; 
        
      l = LocJb[j];
      if(l<N_U)
      {
        // jb in first component
        for(k=0;k<N_DOFperJoint;k++)
          if(EdgeDOF[k] == l)
            break;

        System[ DOF[l]     *LDA + (3*N_Rows+i)] = -LocAlpha[k];
        System[(3*N_Rows+i)*LDA + DOF[l]      ] = -LocAlpha[k];
      }
      else
      {
        if(l<2*N_U)
        {
          // jb in second component
          l -= N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          System[(DOF[l]+N_Rows)*LDA + (3*N_Rows+i)   ] =
                  -LocAlpha[k+N_DOFperJoint];
          System[(3*N_Rows+i)   *LDA + (DOF[l]+N_Rows)] =
                  -LocAlpha[k+N_DOFperJoint];
        }
        else
        {
          // jb in third component
          l -= 2*N_U;
          for(k=0;k<N_DOFperJoint;k++)
            if(EdgeDOF[k] == l)
              break;

          System[(DOF[l]+2*N_Rows)*LDA + (3*N_Rows+i)   ] =
                  -LocAlpha[k+2*N_DOFperJoint];
          System[(3*N_Rows+i)   *LDA + (DOF[l]+2*N_Rows)] =
                  -LocAlpha[k+2*N_DOFperJoint];
        }
      }
    } // endfor j

    // internal bubbles
    N_Inner = FEDesc_Obj->GetN_InnerDOF();
    InnerDOF = FEDesc_Obj->GetInnerDOF();

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      System[k*LDA+k] = 1;
      Rhs[k] = 0;

      k += N_Rows;
      System[k*LDA+k] = 1;
      Rhs[k] = 0;

      k += N_Rows;
      System[k*LDA+k] = 1;
      Rhs[k] = 0;
    } // endfor j
  } // endfor i

  // realise Dirichlet coundary conditions
  m=RowPtr[DirichletBound];
  for(i=DirichletBound;i<N_Rows;i++)
  {
    memset(System+i*LDA, 0, LDA*SizeOfDouble);
    System[ i        *LDA + i       ] = 1;

    memset(System+(i+N_Rows)*LDA, 0, LDA*SizeOfDouble);
    System[(i+N_Rows)*LDA + i+N_Rows] = 1;

    memset(System+(i+2*N_Rows)*LDA, 0, LDA*SizeOfDouble);
    System[(i+2*N_Rows)*LDA + i+2*N_Rows] = 1;
  }

  /*
  // print whole system
  cout << setprecision(20);
  for(i=0;i<N_Eqn;i++)
  {
    for(j=0;j<N_Eqn;j++)
    {
      cout << setw(5) << i+1 << setw(5) << j+1 << setw(30);
      cout << System[i*LDA+j] << endl;
    } // endfor j
  } // endfor i

  // write rhs
  for(i=0;i<N_Eqn;i++)
  {
    cout << setw(5) << i+1 << setw(30) << Rhs[i] << endl;
  }
  */

  // set pressure constant
  memset(System+(N_Eqn-1)*LDA, 0, LDA*SizeOfDouble);
  for(j=0;j<NP;j++)
    System[(N_Eqn-1)*LDA+3*N_Rows+j] = 1;
  Rhs[N_Eqn-1] = 0;

  // solve
  N_Rhs = 1;
  SolveMultipleSystemsNew((double *)System, (double *)Rhs, N_Eqn,
                          LDA, LDB, N_Rhs);

  /*
  cout << "velocity" << endl;
  for(i=0;i<N_Rows;i++)
  {
    if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
    if(fabs(Rhs[i+N_Rows]) < 1e-10) Rhs[i+N_Rows] = 0;
    if(fabs(Rhs[i+2*N_Rows]) < 1e-10) Rhs[i+2*N_Rows] = 0;
    cout << setw(5) << i << setw(15) << Rhs[i];
    cout << setw(15) << Rhs[i+N_Rows];
    cout << setw(15) << Rhs[i+2*N_Rows] << endl;
  }
  cout << endl;
  
  cout << "pressure" << endl;
  for(i=3*N_Rows;i<N_Eqn;i++)
  {
    if(fabs(Rhs[i]) < 1e-10) Rhs[i] = 0;
    cout << setw(5) << Rhs[i] << endl;
  }
  cout << endl;
  */

  memcpy(u1, Rhs, N_Eqn*SizeOfDouble);
  delete System;
  delete Rhs;
}
#endif

/** solve exact on this level */
void TNSE_MGLevel5::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{

  OutPut("TNSE_MGLevel5::SolveExactUMFPACK: Are we here?" << endl);
  OutPut("Not yet implemented!!!" << endl);
  exit(4711);
  
  //DirectSolver(A11, A12, A21, A22, B1T, B2T, B1, B2, rhs1, u1, 3);
}

/** step length control for smoother */
double TNSE_MGLevel5::StepLengthControl (double *u1, double *u1old, 
                                         double *def1, 
                                         int N_Parameters, 
                                         double *Parameters)
{
  double *x,*y,omega,numerator,nominator;
  int i,j;

  // allocate auxiliary array
  x = new double[2*N_DOF];
  y = x+N_DOF;

  // save current update in array x
  for (i=0;i<N_DOF;i++)
    x[i] = u1[i]-u1old[i];
  memset(y,0,N_DOF*SizeOfDouble);
  
  // compute matrix times update
  CoupledMatVect(A, x, y);

  numerator = Ddot(N_DOF,def1,y);
  nominator = Ddot(N_DOF,y,y);

  if (nominator > 0)
    omega = numerator/nominator;
  else
    {
      if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.5;
      OutPut("MESSAGE : Step length control failed. Set omega = " << omega<< endl);
    }
  if (fabs(omega)<0.0001)
    {
       if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.9;
    }     
  delete x;
  if (TDatabase::ParamDB->SC_VERBOSE>=2)
    OutPut("step length control " << omega << endl);
  return(omega); 
}

/** print all matrices and oth right hand sides */
void TNSE_MGLevel5::PrintAll()
{
}

/** Braess--Sarazin smoother  */
void TNSE_MGLevel5::BraessSarazin(double *u1, double *rhs1,
				   double *aux, int N_Parameters, 
				   double *Parameters,int N_Levels)
{  
}

#ifdef __2D__
/** calculate ustar-representation from u-representation */
void TNSE_MGLevel5::GetUstarFromU(double *u, double *ustar)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE2D feid;
  TFE2D *ele;
  TFEDesc2D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  /*
  int N_Outer_Q2 = 8;
  int OuterDOF_Q2[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };

  int N_Outer_N2 = 8;
  int OuterDOF_N2[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

  int N_Outer_N1 = 4;
  int OuterDOF_N1[4] = { 0, 1, 2, 3 };

  int N_Outer_Q3 = 12;
  int OuterDOF_Q3[12] = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 }; 

  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case -1:
      N_Outer = N_Outer_N1;
      OuterDOF = OuterDOF_N1;
    break;

    case -2:
      N_Outer = N_Outer_N2;
      OuterDOF = OuterDOF_N2;
    break;

    case 12:
      N_Outer = N_Outer_Q2;
      OuterDOF = OuterDOF_Q2;
    break;

    case 13:
      N_Outer = N_Outer_Q3;
      OuterDOF = OuterDOF_Q3;
    break;
  } // endswitch
  */

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(ustar, u, 2*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(feid);
    fedesc = ele->GetFEDesc2D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = u[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val += LocAlpha[k]*u[DOF[m]];

            val += LocAlpha[k+N_DOFperJoint]*u[DOF[m]+N_DOF];
          } // endfor k
          ustar[n] = val;
        }
        else
        {
          // jb in second component
          l -= N_U;
          n = DOF[l];
          val = u[n+N_DOF];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val += LocAlpha[k+N_DOFperJoint]*u[DOF[m]+N_DOF];

            val += LocAlpha[k]*u[DOF[m]];
          } // endfor k
          ustar[n+N_DOF] = val;
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      ustar[k] = 0;
      ustar[k+N_DOF] = 0;
    }

  } // endfor i
  // cout << "GetUstarFromU" << endl;
}

/** calculate u-representation from ustar-representation */
void TNSE_MGLevel5::GetUFromUstar(double *ustar, double *u)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE2D feid;
  TFE2D *ele;
  TFEDesc2D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  /*
  int N_Outer_Q2 = 8;
  int OuterDOF_Q2[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };

  int N_Outer_N2 = 8;
  int OuterDOF_N2[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

  int N_Outer_N1 = 4;
  int OuterDOF_N1[4] = { 0, 1, 2, 3 };

  int N_Outer_Q3 = 12;
  int OuterDOF_Q3[12] = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 }; 

  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case -1:
      N_Outer = N_Outer_N1;
      OuterDOF = OuterDOF_N1;
    break;

    case -2:
      N_Outer = N_Outer_N2;
      OuterDOF = OuterDOF_N2;
    break;

    case 12:
      N_Outer = N_Outer_Q2;
      OuterDOF = OuterDOF_Q2;
    break;

    case 13:
      N_Outer = N_Outer_Q3;
      OuterDOF = OuterDOF_Q3;
    break;
  } // endswitch
  */

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(u, ustar, 2*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(feid);
    fedesc = ele->GetFEDesc2D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = ustar[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val -= LocAlpha[k]*ustar[DOF[m]];

            val -= LocAlpha[k+N_DOFperJoint]*ustar[DOF[m]+N_DOF];
          } // endfor k
          u[n] = val;
        }
        else
        {
          // jb in second component
          l -= N_U;
          n = DOF[l];
          val = ustar[n+N_DOF];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val -= LocAlpha[k+N_DOFperJoint]*ustar[DOF[m]+N_DOF];

            val -= LocAlpha[k]*ustar[DOF[m]];
          } // endfor k
          u[n+N_DOF] = val;
        }
      }
    } // endfor j

    // modification due to internal bubbles
    // use already modified values from u !!!
    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();
    fedesc->GetOuterDOF(N_Outer, OuterDOF);

    LocC = C + BeginC[i];

    for(j=0;j<N_Inner;j++)
    {
      val1 = 0; // LocC[ j         *(2*N_Outer+1) + 2*N_Outer];
      val2 = 0; // LocC[(j+N_Inner)*(2*N_Outer+1) + 2*N_Outer];

      for(k=0;k<N_Outer;k++)
      {
        m = DOF[OuterDOF[k]];
        val1 -= LocC[ j         *(2*N_Outer+1) + k        ] * u[m];
        val1 -= LocC[ j         *(2*N_Outer+1) + k+N_Outer] * u[m+N_DOF];
        val2 -= LocC[(j+N_Inner)*(2*N_Outer+1) + k        ] * u[m];
        val2 -= LocC[(j+N_Inner)*(2*N_Outer+1) + k+N_Outer] * u[m+N_DOF];
      } // endfor k

      m = DOF[InnerDOF[j]];
      u[m] = val1;
      u[m+N_DOF] = val2;
    } // endfor j
  } // endfor i
  // cout << "GetUFromUstar" << endl;
} // endof GetUFromUstar

/** calculate dstar-representation from d-representation */
void TNSE_MGLevel5::GetDstarFromD(double *d, double *dstar)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE2D feid;
  TFE2D *ele;
  TFEDesc2D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  /*
  int N_Outer_Q2 = 8;
  int OuterDOF_Q2[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };

  int N_Outer_N2 = 8;
  int OuterDOF_N2[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

  int N_Outer_N1 = 4;
  int OuterDOF_N1[4] = { 0, 1, 2, 3 };

  int N_Outer_Q3 = 12;
  int OuterDOF_Q3[12] = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 }; 

  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case -1:
      N_Outer = N_Outer_N1;
      OuterDOF = OuterDOF_N1;
    break;

    case -2:
      N_Outer = N_Outer_N2;
      OuterDOF = OuterDOF_N2;
    break;

    case 12:
      N_Outer = N_Outer_Q2;
      OuterDOF = OuterDOF_Q2;
    break;

    case 13:
      N_Outer = N_Outer_Q3;
      OuterDOF = OuterDOF_Q3;
    break;
  } // endswitch
  */

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(dstar, d, 2*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(feid);
    fedesc = ele->GetFEDesc2D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();
    fedesc->GetOuterDOF(N_Outer, OuterDOF);

    DOF = GlobalNumbers + BeginIndex[i];

    // mofication due to cell bubbles
    LocC = C + BeginC[i];

    for(j=0;j<N_Inner;j++)
    {
      m = DOF[InnerDOF[j]];
      val1 = dstar[m];
      val2 = dstar[m+N_DOF];
      for(k=0;k<N_Outer;k++)
      {
        m = DOF[OuterDOF[k]];
        dstar[m]       -= LocC[ j         *(2*N_Outer+1) + k        ] * val1;
        dstar[m+N_DOF] -= LocC[ j         *(2*N_Outer+1) + k+N_Outer] * val1;
        dstar[m]       -= LocC[(j+N_Inner)*(2*N_Outer+1) + k        ] * val2;
        dstar[m+N_DOF] -= LocC[(j+N_Inner)*(2*N_Outer+1) + k+N_Outer] * val2;
      } // endfor k
    } // endfor j

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = dstar[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              dstar[DOF[m]] -= LocAlpha[k]*val;

            dstar[DOF[m]+N_DOF] -= LocAlpha[k+N_DOFperJoint]*val;
          } // endfor k
        }
        else
        {
          // jb in second component
          l -= N_U;
          n = DOF[l];
          val = dstar[n+N_DOF];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              dstar[DOF[m]+N_DOF] -= LocAlpha[k+N_DOFperJoint]*val;

            dstar[DOF[m]] -= LocAlpha[k]*val;
          } // endfor k
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      dstar[k] = 0;
      dstar[k+N_DOF] = 0;
    }

  } // endfor i
  // cout << "GetDstarFromD" << endl;
} // endof GetDstarFromD

/** calculate d-representation from dstar-representation */
void TNSE_MGLevel5::GetDFromDstar(double *dstar, double *d)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE2D feid;
  TFE2D *ele;
  TFEDesc2D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  /*
  int N_Outer_Q2 = 8;
  int OuterDOF_Q2[8] = { 0, 1, 2, 3, 5, 6, 7, 8 };

  int N_Outer_N2 = 8;
  int OuterDOF_N2[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };

  int N_Outer_N1 = 4;
  int OuterDOF_N1[4] = { 0, 1, 2, 3 };

  int N_Outer_Q3 = 12;
  int OuterDOF_Q3[12] = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 }; 

  switch(TDatabase::ParamDB->VELOCITY_SPACE)
  {
    case -1:
      N_Outer = N_Outer_N1;
      OuterDOF = OuterDOF_N1;
    break;

    case -2:
      N_Outer = N_Outer_N2;
      OuterDOF = OuterDOF_N2;
    break;

    case 12:
      N_Outer = N_Outer_Q2;
      OuterDOF = OuterDOF_Q2;
    break;

    case 13:
      N_Outer = N_Outer_Q3;
      OuterDOF = OuterDOF_Q3;
    break;
  } // endswitch
  */

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(d, dstar, 2*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE2D(i, cell);
    ele = TFEDatabase2D::GetFE2D(feid);
    fedesc = ele->GetFEDesc2D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();
    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 2*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = d[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              d[DOF[m]] += LocAlpha[k]*val;

            d[DOF[m]+N_DOF] += LocAlpha[k+N_DOFperJoint]*val;
          } // endfor k
        }
        else
        {
          // jb in second component
          l -= N_U;
          n = DOF[l];
          val = d[n+N_DOF];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              d[DOF[m]+N_DOF] += LocAlpha[k+N_DOFperJoint]*val;

            d[DOF[m]] += LocAlpha[k]*val;
          } // endfor k
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      d[k] = 0;
      d[k+N_DOF] = 0;
    }
  } // endfor i
  // cout << "GetDFromDstar" << endl;
} // GetDFromDstar

#endif // __2D__

#ifdef __3D__

/** calculate ustar-representation from u-representation */
void TNSE_MGLevel5::GetUstarFromU(double *u, double *ustar)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE3D feid;
  TFE3D *ele;
  TFEDesc3D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2, val3;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(ustar, u, 3*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(feid);
    fedesc = ele->GetFEDesc3D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = u[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val += LocAlpha[k]*u[DOF[m]];

            val += LocAlpha[k+  N_DOFperJoint]*u[DOF[m]+  N_DOF];
            val += LocAlpha[k+2*N_DOFperJoint]*u[DOF[m]+2*N_DOF];
          } // endfor k
          ustar[n] = val;
        }
        else
        {
          if(l<2*N_U)
          {
            // jb in second component
            l -= N_U;
            n = DOF[l];
            val = u[n+N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                val += LocAlpha[k+N_DOFperJoint]*u[DOF[m]+N_DOF];

              val += LocAlpha[k                ]*u[DOF[m]        ];
              val += LocAlpha[k+2*N_DOFperJoint]*u[DOF[m]+2*N_DOF];
            } // endfor k
            ustar[n+N_DOF] = val;
          }
          else
          {
            // jb in third component
            l -= 2*N_U;
            n = DOF[l];
            val = u[n+2*N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                val += LocAlpha[k+2*N_DOFperJoint]*u[DOF[m]+2*N_DOF];

              val += LocAlpha[k              ]*u[DOF[m]      ];
              val += LocAlpha[k+N_DOFperJoint]*u[DOF[m]+N_DOF];
            } // endfor k
            ustar[n+2*N_DOF] = val;
          }
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      ustar[k] = 0;
      ustar[k+N_DOF] = 0;
      ustar[k+2*N_DOF] = 0;
    }

  } // endfor i
  // cout << "GetUstarFromU" << endl;
}

/** calculate u-representation from ustar-representation */
void TNSE_MGLevel5::GetUFromUstar(double *ustar, double *u)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE3D feid;
  TFE3D *ele;
  TFEDesc3D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2, val3;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(u, ustar, 3*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(feid);
    fedesc = ele->GetFEDesc3D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = ustar[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              val -= LocAlpha[k]*ustar[DOF[m]];

            val -= LocAlpha[k+  N_DOFperJoint]*ustar[DOF[m]+  N_DOF];
            val -= LocAlpha[k+2*N_DOFperJoint]*ustar[DOF[m]+2*N_DOF];
          } // endfor k
          u[n] = val;
        }
        else
        {
          if(l<2*N_U)
          {
            // jb in second component
            l -= N_U;
            n = DOF[l];
            val = ustar[n+N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                val -= LocAlpha[k+N_DOFperJoint]*ustar[DOF[m]+N_DOF];

              val -= LocAlpha[k                ]*ustar[DOF[m]        ];
              val -= LocAlpha[k+2*N_DOFperJoint]*ustar[DOF[m]+2*N_DOF];
            } // endfor k
            u[n+N_DOF] = val;
          }
          else
          {
            // jb in second component
            l -= 2*N_U;
            n = DOF[l];
            val = ustar[n+2*N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                val -= LocAlpha[k+2*N_DOFperJoint]*ustar[DOF[m]+2*N_DOF];

              val -= LocAlpha[k              ]*ustar[DOF[m]      ];
              val -= LocAlpha[k+N_DOFperJoint]*ustar[DOF[m]+N_DOF];
            } // endfor k
            u[n+2*N_DOF] = val;
          }
        }
      }
    } // endfor j

    // modification due to internal bubbles
    // use already modified values from u !!!
    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();
    fedesc->GetOuterDOF(N_Outer, OuterDOF);

    LocC = C + BeginC[i];

    for(j=0;j<N_Inner;j++)
    {
      val1 = 0; // LocC[ j           *(3*N_Outer+1) + 3*N_Outer];
      val2 = 0; // LocC[(j+  N_Inner)*(3*N_Outer+1) + 3*N_Outer];
      val3 = 0; // LocC[(j+2*N_Inner)*(3*N_Outer+1) + 3*N_Outer];

      for(k=0;k<N_Outer;k++)
      {
        m = DOF[OuterDOF[k]];
        val1 -= LocC[ j           *(3*N_Outer+1) + k          ] * u[m        ];
        val1 -= LocC[ j           *(3*N_Outer+1) + k+  N_Outer] * u[m+  N_DOF];
        val1 -= LocC[ j           *(3*N_Outer+1) + k+2*N_Outer] * u[m+2*N_DOF];
        val2 -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k          ] * u[m        ];
        val2 -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k+  N_Outer] * u[m+  N_DOF];
        val2 -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k+2*N_Outer] * u[m+2*N_DOF];
        val3 -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k          ] * u[m        ];
        val3 -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k+  N_Outer] * u[m+  N_DOF];
        val3 -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k+2*N_Outer] * u[m+2*N_DOF];
      } // endfor k

      m = DOF[InnerDOF[j]];
      u[m        ] = val1;
      u[m  +N_DOF] = val2;
      u[m+2*N_DOF] = val3;
    } // endfor j
  } // endfor i
  // cout << "GetUFromUstar" << endl;
} // endof GetUFromUstar

/** calculate dstar-representation from d-representation */
void TNSE_MGLevel5::GetDstarFromD(double *d, double *dstar)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE3D feid;
  TFE3D *ele;
  TFEDesc3D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2, val3;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(dstar, d, 3*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(feid);
    fedesc = ele->GetFEDesc3D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();

    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();
    fedesc->GetOuterDOF(N_Outer, OuterDOF);

    DOF = GlobalNumbers + BeginIndex[i];

    // mofication due to cell bubbles
    LocC = C + BeginC[i];

    for(j=0;j<N_Inner;j++)
    {
      m = DOF[InnerDOF[j]];
      val1 = dstar[m        ];
      val2 = dstar[m+  N_DOF];
      val3 = dstar[m+2*N_DOF];
      for(k=0;k<N_Outer;k++)
      {
        m = DOF[OuterDOF[k]];
        dstar[m        ] -= LocC[ j           *(3*N_Outer+1) + k          ] * val1;
        dstar[m+  N_DOF] -= LocC[ j           *(3*N_Outer+1) + k+  N_Outer] * val1;
        dstar[m+2*N_DOF] -= LocC[ j           *(3*N_Outer+1) + k+2*N_Outer] * val1;
        dstar[m        ] -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k          ] * val2;
        dstar[m+  N_DOF] -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k+  N_Outer] * val2;
        dstar[m+2*N_DOF] -= LocC[(j+  N_Inner)*(3*N_Outer+1) + k+2*N_Outer] * val2;
        dstar[m        ] -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k          ] * val3;
        dstar[m+  N_DOF] -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k+  N_Outer] * val3;
        dstar[m+2*N_DOF] -= LocC[(j+2*N_Inner)*(3*N_Outer+1) + k+2*N_Outer] * val3;
      } // endfor k
    } // endfor j

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = dstar[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              dstar[DOF[m]] -= LocAlpha[k]*val;

            dstar[DOF[m]+  N_DOF] -= LocAlpha[k+  N_DOFperJoint]*val;
            dstar[DOF[m]+2*N_DOF] -= LocAlpha[k+2*N_DOFperJoint]*val;
          } // endfor k
        }
        else
        {
          if(l<2*N_U)
          {
            // jb in second component
            l -= N_U;
            n = DOF[l];
            val = dstar[n+N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                dstar[DOF[m]+N_DOF] -= LocAlpha[k+N_DOFperJoint]*val;

              dstar[DOF[m]] -= LocAlpha[k]*val;
              dstar[DOF[m]+2*N_DOF] -= LocAlpha[k+2*N_DOFperJoint]*val;
            } // endfor k
          }
          else
          {
            // jb in third component
            l -= 2*N_U;
            n = DOF[l];
            val = dstar[n+2*N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                dstar[DOF[m]+2*N_DOF] -= LocAlpha[k+2*N_DOFperJoint]*val;

              dstar[DOF[m]] -= LocAlpha[k]*val;
              dstar[DOF[m]+N_DOF] -= LocAlpha[k+N_DOFperJoint]*val;
            } // endfor k
          }
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      dstar[k] = 0;
      dstar[k+N_DOF] = 0;
      dstar[k+2*N_DOF] = 0;
    }

  } // endfor i
  // cout << "GetDstarFromD" << endl;
} // endof GetDstarFromD

/** calculate d-representation from dstar-representation */
void TNSE_MGLevel5::GetDFromDstar(double *dstar, double *d)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Joints,N_DOF, N_U;
  int *GlobalNumbers, *BeginIndex;
  int *DOF, **JointDOFs, *EdgeDOF;
  FE3D feid;
  TFE3D *ele;
  TFEDesc3D *fedesc;
  TCollection *coll;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  int *LocJb;
  double *LocAlpha, val;
  double *LocC, val1, val2, val3;
  int *OuterDOF, N_Outer;
  int *InnerDOF, N_Inner;
  int *BeginJb, *jb, N_DOFperJoint;
  double *Alpha;
  int *BeginC;
  double *C;

  coll = USpace->GetCollection();
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();
  N_DOF = USpace->GetN_DegreesOfFreedom();

  memcpy(d, dstar, 3*N_DOF*SizeOfDouble);

  BeginJb = A->GetBeginJb();
  jb = A->GetJb();
  N_DOFperJoint = A->GetN_DOFperJoint();
  Alpha = A->GetAlpha();
  BeginC = A->GetBeginC();
  C = A->GetC();

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    coll->GetCell(i)->SetClipBoard(i);

  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);

    N_Joints = cell->GetN_Joints();

    feid = USpace->GetFE3D(i, cell);
    ele = TFEDatabase3D::GetFE3D(feid);
    fedesc = ele->GetFEDesc3D();
    JointDOFs = fedesc->GetJointDOF();
    N_U = fedesc->GetN_DOF();
    N_Inner = fedesc->GetN_InnerDOF();
    InnerDOF = fedesc->GetInnerDOF();

    DOF = GlobalNumbers + BeginIndex[i];

    // modification due to edge bubble (jb)
    LocJb = jb+BeginJb[i];

    for(j=0;j<N_Joints;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);
      if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
      {
        // edge update has to be done by this cell

        // get local data
        LocAlpha = Alpha + 3*N_DOFperJoint*(BeginJb[i]+j);
        EdgeDOF = JointDOFs[j]; 
        
        l = LocJb[j];
        if(l<N_U)
        {
          // jb in first component
          n = DOF[l];
          val = d[n];
          for(k=0;k<N_DOFperJoint;k++)
          {
            m = EdgeDOF[k];
            if(m != l)
              d[DOF[m]] += LocAlpha[k]*val;

            d[DOF[m]+N_DOF] += LocAlpha[k+N_DOFperJoint]*val;
            d[DOF[m]+2*N_DOF] += LocAlpha[k+2*N_DOFperJoint]*val;
          } // endfor k
        }
        else
        {
          if(l<2*N_U)
          {
            // jb in second component
            l -= N_U;
            n = DOF[l];
            val = d[n+N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                d[DOF[m]+N_DOF] += LocAlpha[k+N_DOFperJoint]*val;

              d[DOF[m]] += LocAlpha[k]*val;
              d[DOF[m]+2*N_DOF] += LocAlpha[k+2*N_DOFperJoint]*val;
            } // endfor k
          }
          else
          {
            // jb in third component
            l -= 2*N_U;
            n = DOF[l];
            val = d[n+2*N_DOF];
            for(k=0;k<N_DOFperJoint;k++)
            {
              m = EdgeDOF[k];
              if(m != l)
                d[DOF[m]+2*N_DOF] += LocAlpha[k+2*N_DOFperJoint]*val;

              d[DOF[m]] += LocAlpha[k]*val;
              d[DOF[m]+N_DOF] += LocAlpha[k+N_DOFperJoint]*val;
            } // endfor k
          }
        }
      }
    } // endfor j

    for(j=0;j<N_Inner;j++)
    {
      k = DOF[InnerDOF[j]];
      d[k] = 0;
      d[k+N_DOF] = 0;
      d[k+2*N_DOF] = 0;
    }
  } // endfor i
  // cout << "GetDFromDstar" << endl;
} // GetDFromDstar
#endif // __3D__
