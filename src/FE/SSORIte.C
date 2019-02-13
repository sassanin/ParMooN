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
// @(#)SSORIte.C        1.24 06/27/00
//
// Class:       TSSORIte
// Purpose:     iteration methods
//
// Author:      Volker John 24.10.2000
//
// History:     24.10.2000 start of implementation
//
// =======================================================================
#include <ItMethod.h>
#include <SSORIte.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#include <stdlib.h>
#include <stdio.h>

/** constructor with initialization */
TSSORIte::TSSORIte(MatVecProc *MatVec,
DefectProc *Defect,
TItMethod *Prec,
int n_aux, int n_dof,
int scalar)
: TItMethod(MatVec, Defect, Prec, n_aux, n_dof)

{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
}


// destructor
TSSORIte::~TSSORIte()
{
  if (N_Aux>0)
  {
    delete AuxArray;
    delete AuxArray[0];
  }
}


// SSOR iteration
int TSSORIte::Iterate (TSquareMatrix **sqmat,
TMatrix **mat, double *sol,
double *rhs)
{
  int i, j, k, l, m, *ARowPtr, *AKCol;
  double *AEntries, om, s, diag;

  ARowPtr = sqmat[0]->GetRowPtr();
  AKCol = sqmat[0]->GetKCol();
  AEntries = sqmat[0]->GetEntries();
  // damping parameter
  om = TDatabase::ParamDB->SC_SOR_OMEGA;

  // forward sweep
  for (i=0; i<N_DOF;i++)
  {
    j = ARowPtr[i];
    l = ARowPtr[i+1];
    s = 0.0;
    for (k=j;k<l;k++)
    {
      m = AKCol[k];
      // diagonal entry
      if (m==i)
        diag = AEntries[k];
      s += AEntries[k]*sol[m];
    }
    sol[i] =  om*(rhs[i]-s)/diag;
  }

  // backward sweep
  for (i=N_DOF-1; i>=0; i--)
  {
    j = ARowPtr[i];
    l = ARowPtr[i+1];
    s = 0.0;
    for (k=j;k<l;k++)
    {
      m = AKCol[k];
      // diagonal entry
      if (m==i)
        diag = AEntries[k];
      s += AEntries[k]*sol[m];
    }
    sol[i] =  om*(rhs[i]-s)/diag;
  }
  return(0);
}
