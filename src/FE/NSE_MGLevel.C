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
// @(#)NSE_MGLevel.C        1.16 07/03/00
//
// Class:       TNSE_MGLevel
// Purpose:     store all data for one level in a multi grid method
//              for solving a Stokes-/ Navier-Stokes system
//              abstract super class
//
// Author:      Volker John, 28.07.2000
//
// History:     28.07.2000 start of implementation
//
// =======================================================================

#include <NSE_MGLevel.h>
#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __2D__
   #include <FESpace2D.h>
   #include <FEDatabase2D.h>
#endif  
#ifdef __3D__
   #include <FESpace3D.h>
   #include <FEDatabase3D.h>
#endif  

#include <stdlib.h>
#include <string.h>

#include <LinAlg.h>

/** constructor */
TNSE_MGLevel::TNSE_MGLevel(int level,
             double *f1,
             double *u1,
             int n_aux, double *al, int velocity_space, 
             int pressure_space, TCollection *coll)
{
  int i;
  double *aux;

  Level = level;

  Rhs1 = f1;

  U1 = u1;

  VelocitySpace = velocity_space;
  PressureSpace = pressure_space;

  Additional = NULL;

  alpha = al[0];
  beta = al[1];

  VankaColl = coll;

  Type = 0;
}

/** destructor */
TNSE_MGLevel::~TNSE_MGLevel()
{
  delete Aux[0];
  delete Aux;

  if(Additional) 
    delete Additional;
} // ~TNSE_MGLevel

/** return i-th auxiliary vector */
double *TNSE_MGLevel::GetAuxVector(int i)
{
  double *ret;

  if(i<N_Aux)
    ret = Aux[i];
  else
  {
    cerr << "Not enough aux vectors in NSE_MGLevel!" << endl;
    exit(-1);
    ret = NULL;
  }

  return ret;
} // GetAuxVector

/** calculate defect */
void TNSE_MGLevel::Defect(double *u1, double *f1, double *d1, double &res)
{
}

/** update solution */
void TNSE_MGLevel::Update(double *u1, double *v1)
{
  int i;
  double omega=beta; // value from GMG DAMP FACTOR SADDLE / FINE

  for(i=0;i<N_DOF;i++)
    u1[i] += omega*v1[i];
}

/** correct Dirichlet and hanging nodes */
void TNSE_MGLevel::CorrectNodes(double *u1)
{
  int i,j,k, index;
  double s, t;

  memset(u1+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
  memset(u1+N_UDOF+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);
#ifdef __2D__
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(u1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
#ifdef __3D__
  memset(u1+2*N_UDOF+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(u1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
}

/** set hanging nodes in order to satisfy coupling condition */
void TNSE_MGLevel::SetHangingNodes(double *u1)
{
  int i, j, k;
  double value, value1, value2;
  int N_Hanging, N_Nodes;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling;
  int *DOF;
  int FirstHangingNodeNumber;
  double *u2;

  u2 = u1 + N_UDOF;
 
#ifdef __2D__
  N_Hanging = USpace->GetN_Hanging();
  if(N_Hanging)
  {
    HangingNodes = USpace->GetHangingNodes();
    FirstHangingNodeNumber = USpace->GetActiveBound();
    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      N_Nodes = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      // cout << "HN: " << FirstHangingNodeNumber+i << endl;
      
      value1 = 0;
      value2 = 0;
      for(j=0;j<N_Nodes;j++)
      {
        value = Coupling[j];
        k = DOF[j];
        // cout << j << " " << k << " " << Coupling[j] << endl;
        value1 += value * u1[k];
        value2 += value * u2[k];
      }
      // cout << u1[FirstHangingNodeNumber+i]-value1 << " " << u2[FirstHangingNodeNumber+i]-value2 << endl;
      u1[FirstHangingNodeNumber+i] = value1;
      u2[FirstHangingNodeNumber+i] = value2;
    } // endfor i
  } // endif N_Hanging
#endif  
}

/** correct defect */
void TNSE_MGLevel::CorrectDefect(double *v1)
{
  int i;
  double sum;

  int N_Hanging, N_Nodes;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling;
  int *DOF;
  int FirstHangingNodeNumber;
  double value, value1, value2;
  int j,k;
  double *v2;

  v2 = v1 + N_UDOF;

  /*
  cout << "Level: " << Level << endl;
  cout << "N_UDOF: " << N_UDOF << endl;
  cout << "N_PDOF: " << N_PDOF << endl;
  cout << "FirstHangingNodeNumber: " << USpace->GetActiveBound() << endl;
  cout << "N_Hanging: " << USpace->GetN_Hanging() << endl;
  for(i=0;i<N_UDOF;i++)
  {
    cout << setw(4) << i << setw(20) << v1[i] << setw(20) << v2[i] << endl;
  }
  */

// /*
#ifdef __2D__
  N_Hanging = USpace->GetN_Hanging();
  if(N_Hanging)
  {
    HangingNodes = USpace->GetHangingNodes();
    FirstHangingNodeNumber = USpace->GetActiveBound();
    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      N_Nodes = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      // cout << "HN: " << FirstHangingNodeNumber+i << endl;
      
      value1 = v1[FirstHangingNodeNumber+i];
      value2 = v2[FirstHangingNodeNumber+i];
      for(j=0;j<N_Nodes;j++)
      {
        value = Coupling[j];
        k = DOF[j];
        // cout << j << " " << k << " " << Coupling[j] << endl;
        v1[k] += value1*value;
        v2[k] += value2*value;
      }
    } // endfor i
  } // endif N_Hanging
#endif  
// */

  memset(v1+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
  memset(v2+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));
#ifdef __2D__
  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector2D(v1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
#ifdef __3D__
  memset(v1+2*N_UDOF+N_Active, 0, SizeOfDouble*(N_UDOF-N_Active));

  if (TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE) 
    IntoL20Vector3D(v1+GEO_DIM*N_UDOF, N_PDOF, PressureSpace);
#endif  
}

/** reset vector to zero */
void TNSE_MGLevel::Reset(double *v1)
{
  memset(v1, 0, N_DOF*SizeOfDouble);
}

/** Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel::CellVanka(double *u1, double *rhs1, double *aux,
                              int N_Parameters, double *Parameters, 
                              int smoother, int N_Levels)
{
} // end Vanka

/** nodal Vanka smoother, GAUSS-SEIDEL type */
void TNSE_MGLevel::NodalVanka(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters, int smoother, int N_Levels)
{
} // end Vanka


/** solve exact on this level */
void TNSE_MGLevel::SolveExact(double *u1, double *rhs1)
{
}

/** solve exact on this level */
void TNSE_MGLevel::SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag)
{
}

/** Braess--Sarazin smoother */
void TNSE_MGLevel::BraessSarazin(double *u1, double *rhs1, double *aux,
    int N_Parameters, double *Parameters,int N_Levels)
{
}
/** step length control for Vanka */
double TNSE_MGLevel::StepLengthControl(double *u1, double *u1old, double *def1, 
                                        int N_Parameters, double *Parameters)
{
  return 1;
}

/** print all matrices and oth right hand sides */
void TNSE_MGLevel::PrintAll()
{
}
#ifdef _MPI
 void TNSE_MGLevel::UpdateHaloRhs(double*a, double*b){
 } 
#endif
