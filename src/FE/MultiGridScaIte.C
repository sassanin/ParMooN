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
// @(#)MultiGridScaIte.C        1.24 06/27/00
//
// Class:       TMultiGridScaIte
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <MultiGridScaIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>

#ifdef __2D__
  #include <MultiGrid2D.h>
#else
  #include <MultiGrid3D.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** constructor with initialization */
#ifdef __2D__
TMultiGridScaIte::TMultiGridScaIte(MatVecProc *MatVec, 
                             DefectProc *Defect, 
                             TItMethod *Prec,
                             int n_aux, int n_dof,
                             TMultiGrid2D *MG,
                             int zero_start)
#else
TMultiGridScaIte::TMultiGridScaIte(MatVecProc *MatVec, 
                             DefectProc *Defect, 
                             TItMethod *Prec,
                             int n_aux, int n_dof,
                             TMultiGrid3D *MG,
                             int zero_start)
#endif
  : TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
                                    
{
  // initialize only mg
  // every other thing is initialized in the definition of the
  // multigrid
  mg = MG;
  N_DOF = n_dof;
  N_Aux = 0;
  Zero_Start=zero_start;
}


TMultiGridScaIte::~TMultiGridScaIte()
{

}

double tCyc=0.0;
int TMultiGridScaIte::Iterate (TSquareMatrix **sqmat,
                               TMatrix **mat, double *sol, 
                               double *rhs)
{
  double res, *mgsol, *mgrhs;
  double t1,t2;
  
  // set data for multigrid cycle
  mgsol= mg->GetLevel(mg->GetN_Levels()-1)->GetSolution();
  mgrhs= mg->GetLevel(mg->GetN_Levels()-1)->GetRhs();
  //memcpy(mgsol, sol, N_DOF*SizeOfDouble);
  if (Zero_Start)
  { memset(mgsol, 0, N_DOF*SizeOfDouble);}
  else
  { memcpy(mgsol, sol, N_DOF*SizeOfDouble);}
  memcpy(mgrhs, rhs, N_DOF*SizeOfDouble);
  mg->SetDirichletNodes(mg->GetN_Levels()-1);
  mg->SetRecursion(mg->GetN_Levels()-1);
  // one multigrid cycle
#ifdef _MPI
      //do not know why SetDirichletNodes is not needed for MPI, Check!!!
  //cout << "TMultiGridScaIte::Iterate " <<endl;
  //exit(0);
  t1 = MPI_Wtime();
#else
  t1 = GetTime();
#endif
  
 mg->Cycle(mg->GetN_Levels()-1, res); 
 
#ifdef _MPI
  t2 = MPI_Wtime();
#else
  t2 = GetTime();
#endif
  tCyc += (t2-t1);
   
  // store solution on rhs
  memcpy(sol, mgsol, N_DOF*SizeOfDouble);

  return(0);
}                        
