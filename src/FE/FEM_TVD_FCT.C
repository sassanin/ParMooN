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
// FEM_TVD_FCT.C
//
// Purpose:     fem-tvd and fem-fct methods for scalar cdr equations
//
// Authors:     Volker John, 2008/02/07
//
// =======================================================================

#include <Database.h>
#include <LinAlg.h>
#include <Solver.h>
#ifdef __2D__
#include <SquareMatrix2D.h>
#include <DiscreteForm2D.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <FE2D.h>
#endif
#ifdef __3D__
#include <SquareMatrix3D.h>
#include <DiscreteForm3D.h>
#include <FEFunction3D.h>
#include <FEDatabase3D.h>
#include <FE3D.h>
#endif
#include <FEM_TVD_FCT.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>

#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

//**************************************************************
//compute the lumped matrix
//output is a vector
//**************************************************************
#ifdef __2D__
void LumpMassMatrixToVector(TSquareMatrix2D *M, double *lump_mass)
#endif
#ifdef __3D__
void LumpMassMatrixToVector(TSquareMatrix3D *M, double *lump_mass)
#endif
{
  double *Entries;
  int *RowPtr, i, j, rows, j0, j1;

  RowPtr        = M->GetRowPtr();
  Entries       = M->GetEntries();
  rows          = M->GetN_Rows();

  memset(lump_mass, 0, rows*SizeOfDouble);
  for (i=0; i<rows; i++)
  {
    lump_mass[i]=0.0;
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for (j=j0;j<j1;j++)
      lump_mass[i] += Entries[j];

    if(lump_mass[i]==0)
    {
      OutPut("zero entry in lumped matrix "<< i << " "  << lump_mass[i] << endl);
      exit(4711);
    }
  }
}


/*******************************************************************************/
//
// FEM_TVD_ForConvDiff for steady-state cdr equations
// following D. Kuzmin (2007)
//
// inputs : *sqmatrix         - system matrix
//          N_U               - total number of dof
//          N_Active          - active dof, should be the same as N_U
//          *matrix_D_Entries - entries of matrix D
//          *sol              - current solution
//          N_neum_to_diri    - number of Dirichlet dof which are treated as
//                              Neumann dof, MUST BE ORDERED !!!
//          *neum_to_diri     - array which contains the indices of Dirichlet
//                              dof which are treated as Neumann dof
//          compute_matrix_D  - flag which says if matrix_D_Entries should be
//                              computed
//
// outputs: *rhs              - right hand side of FEM-TVD problem
//          *matrix_D_Entries - entries of matrix D if compute_matrix_D == 1
//
/*******************************************************************************/
#ifdef __2D__
void FEM_TVD_ForConvDiff(TSquareMatrix2D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D)
#endif
#ifdef __3D__
void FEM_TVD_ForConvDiff(TSquareMatrix3D *sqmatrix, int N_U, int N_Active,
double *matrix_D_Entries, double *sol, double *rhs,
int N_neum_to_diri, int *neum_to_diri, int compute_matrix_D)
#endif
{
  int *ColInd, *RowPtr, N_Entries;
  int i,j,j0,j1,j2,j3,jj,index;
  double nenner, zaehler;
  double *Entries, *F;
  double *P_plus, *P_minus, *Q_plus, *Q_minus, *R_plus, *R_minus, val;
 
  if (N_Active < N_U)
  {
    OutPut("N_Active < N_U ("<< N_Active<< "<"<< N_U <<
      ") !!! FOR APPLYING ALGEBRAIC FLUX CORRECTION, THE BOUNDARY CONDITIONS SHOULD BE NEUMANN !!!" << endl);
    exit(4711);
  }
  // get pointers to columns, rows and entries of matrix A
  ColInd = sqmatrix->GetKCol();
  RowPtr = sqmatrix->GetRowPtr();
  Entries = sqmatrix->GetEntries();
  N_Entries = sqmatrix->GetN_Entries();

  // allocate memory for array F
  F = new double[N_Entries+6*N_U];
  memset(F, 0, (N_Entries+6*N_U)*SizeOfDouble);
  P_plus = F + N_Entries;
  P_minus = P_plus + N_U;
  Q_plus = P_minus + N_U;
  Q_minus = Q_plus + N_U;
  R_plus = Q_minus + N_U;
  R_minus = R_plus + N_U;

  // D has to be computed only once
  if (compute_matrix_D)
  {
    // loop over all rows
    for(i=0;i<N_U;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // compute first the matrix D
      for(j=j0;j<j1;j++)
      {
        // column
        index = ColInd[j];
        // only off-diagonals
        if (index!=i)
        {
          if (Entries[j] > 0)
            matrix_D_Entries[j] = -Entries[j];
          // now check the transposed entry
          j2 = RowPtr[index];
          j3 = RowPtr[index+1];
          for (jj=j2;jj<j3;jj++)
          {
            if (ColInd[jj]==i)
            {
              if (-Entries[jj]<matrix_D_Entries[j])
                matrix_D_Entries[j] = -Entries[jj];
              break;
            }
          }
        }
      }
    }

    // compute diagonal entry of D
    // loop over all rows
    for(i=0;i<N_Active;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      val = 0;
      // add all entries of i-th row
      for(j=j0;j<j1;j++)
      {
        val +=  matrix_D_Entries[j];
        index = ColInd[j];
        if (index==i)
          jj = j;
      }
      matrix_D_Entries[jj] = -val;
    }
    // matrix D is computed
  }
  // add this matrix to A giving \tilde A (Entries)
  // this is the matrix with the properties of an M matrix
  Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);

  // compute matrix F
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      // d_ij (u_i - u_j)
      F[j] = matrix_D_Entries[j] * (sol[index]-sol[i]);
    }
  }
  // matrix F is computed

  // compute flux limiters
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for(j=j0;j<j1;j++)
    {
      if (Entries[j] > 0)
        continue;
      // column
      index = ColInd[j];
      // check transposed entry -> jj
      // diagonal
      if (index==i)
        continue;
      j2 = RowPtr[index];
      j3 = RowPtr[index+1];
      for (jj=j2;jj<j3;jj++)
      {
        if (ColInd[jj]==i)
        {
          //OutPut(Entries[j] << " " << Entries[jj] <<endl);
          break;
        }
      }
      // check upwind condition
      // this ensures that the 'link' between i and index is treated only once
      if (Entries[jj] > Entries[j])
        continue;
      // only the active part of the matrix
      if (F[j] > 0)
      {
        P_plus[i] += F[j];
        if (index<N_U)
          Q_plus[index] += F[j];
        Q_minus[i] -= F[j];
      }
      if (F[j] < 0)
      {
        P_minus[i] += F[j];
        Q_plus[i] -= F[j];
        if (index<N_U)
          Q_minus[index] +=  F[j];
      }
    }                                             // end loop j
  }

  // apply the nodal correction factor evaluated at the upwind node i
  // loop over all nodes
  if (TDatabase::ParamDB->P9!=4711)
  { // original but discontinuous proposal
    for(i=0;i<N_U;i++)
    {  
      if (fabs(P_plus[i])>0)
      {
        R_plus[i] = Q_plus[i]/P_plus[i];
        if (R_plus[i] >1)
          R_plus[i] = 1;
      }
      if (fabs(P_minus[i])>0)
      {
        R_minus[i] = Q_minus[i]/P_minus[i];
        if (R_minus[i] >1)
          R_minus[i] = 1;
      }
      //OutPut(" P " << P_plus[i] << " " <<  P_minus[i] << " ");
      //OutPut(R_plus[i] << " " <<  R_minus[i] << endl);
    }
  }
  else
  {
    // continuous proposal 
  
    for(i=0;i<N_U;i++)
    { 
      zaehler =  Q_plus[i];
      if (-Q_minus[i] < zaehler)
        zaehler = -Q_minus[i];
      nenner = 1e-32;
      if (P_plus[i] > nenner)
        nenner = P_plus[i];
      if (-P_minus[i] > nenner)
        nenner = -P_minus[i];
      //OutPut(zaehler << " " << nenner );
      R_plus[i] = zaehler/nenner;
      if (R_plus[i] > 1)
        R_plus[i] = 1;
      R_minus[i] = R_plus[i];
      //OutPut("new " << i << " P " << P_plus[i] << " "  << P_minus[i] <<" Q " << Q_plus[i] << " "  << Q_minus[i] << 
      //" R " << R_plus[i] << " " << R_minus[i] << endl);
    }
  }


  // treat Dirichlet nodes
  for (j=0;j<N_neum_to_diri;j++)
  {
    i=neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }

  // apply the flux limiters
  // loop over all rows
  for(i=0;i<N_Active;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for(j=j0;j<j1;j++)
    {
      // column
      index = ColInd[j];
      if (index==i)
        continue;
       // this should not happen
      if (Entries[j] > 0)
      {
        OutPut("positive entry in FEMTVD " << i << " " << j << " " << Entries[j] << endl);
        exit(4711);
        continue;
      }

      // check transposed entry
      j2 = RowPtr[index];
      j3 = RowPtr[index+1];
      for (jj=j2;jj<j3;jj++)
      {
        if (ColInd[jj]==i)
        {
          //OutPut(Entries[j] << " " << Entries[jj] <<endl);
          break;
        }
      }
      
      if (TDatabase::ParamDB->P8!=4711)
      {
        // original, symmetric application 
        // check upwind condition
        // this ensures that the 'link' between i and index is treated only once
        if (Entries[jj] > Entries[j])
          continue;
        //OutPut(R_plus[i] << " " << R_minus[i] << " : " << R_plus[index] << " " << R_minus[index] << "::");
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row 
        rhs[i] += F[j];
        // update rhs wrt to current column 
        // note that F[j] = -F[jj] and alpha_j = alpha_jj (symmetry of alpha matrix)
        if (index<N_Active)
          rhs[index] -= F[j];
      }
      else
      { // nonsymmetric application 
        // compute contribution to rhs
        if (F[j] > 0)
        {
          F[j] = R_plus[i]*F[j];
        }
        else
          F[j] = R_minus[i]*F[j];
        // update rhs of current row 
        rhs[i] += F[j]; 
      }
    }
  }
  delete [] F;
}

//**************************************************************
//compute system matrix for FEM-FCT
//output is a vector
//**************************************************************
#ifdef __2D__
void FEM_FCT_SystemMatrix(TSquareMatrix2D *M_C, TSquareMatrix2D *A,
double *lump_mass,int N_U)
#endif
#ifdef __3D__
void FEM_FCT_SystemMatrix(TSquareMatrix3D *M_C, TSquareMatrix3D *A,
double *lump_mass,int N_U)
#endif
{
  int *ColInd, *RowPtr, N_Entries, *ColInd_M, *RowPtr_M, N_Entries_M;
  int i,j,j0,j1,index;
  double *Entries,*Entries_M, cfl;
  double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;

  // get pointers to columns, rows and entries of matrix M_C
  ColInd_M = M_C->GetKCol();
  RowPtr_M = M_C->GetRowPtr();
  Entries_M = M_C->GetEntries();
  N_Entries_M = M_C->GetN_Entries();

  // get pointers to columns, rows and entries of matrix A
  ColInd = A->GetKCol();
  RowPtr = A->GetRowPtr();
  Entries = A->GetEntries();
  N_Entries = A->GetN_Entries();

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr_M[i];
    j1 = RowPtr_M[i+1];

    for(j=j0;j<j1;j++)
    {
      index = ColInd_M[j];
      if(i == index)
      {
        //diagonale
        Entries_M[j] = lump_mass[i]+delta_t*theta1*Entries[j];
        cfl = lump_mass[i]/Entries[j];
        if (theta2 >0)
        {
          cfl /= theta2;
          if (delta_t > cfl)
          {
            OutPut(lump_mass[i] << " " << "cfl violated: cfl " << cfl << " delta t:" << delta_t <<endl);
          }
        }
      }
      else
      {
        Entries_M[j]= delta_t*theta1*Entries[j];
        if ( Entries_M[j]>1e-10)
          OutPut( Entries_M[j] << " ");

      }
    }
  }
}


//**************************************************************
// MINMOD prelimiter
//**************************************************************
double MinMod(double a, double b)
{
  if (a*b < 0)
  {
    return 0.0;
  }
  if (fabs(a) <  fabs(b))
  {
    return a;
  }
  else
  {
    return b;
  }
}


/*******************************************************************************/
//
// FCT-FEM algorithm
// following D. Kuzmin, M. M"oller (2005) (nonlinear scheme)
//           D. Kuzmin (2008) (linear scheme)
//
// inputs : M_C              - consistent mass matrix
//          A                - stiffness matrix
//          N_U              - number of all unknowns
//          N_Active         - number of active unknows (should be the same
//                             as N_U !!!)
//          lump_mass        - lumped mass matrix, stored as array
//          matrix_D_Entries - set only once in each discrete time, in this
//                             routine
//          sol              - array for solution, contains the current
//                             approximation on the solution
//          oldsol           - array with solution from previous discrete time
//          rhs              - array with right hand side f from current
//                             discrete time
//          rhs_old          - rhs f from previous discrete time
//          tilde_u          - low order solution at (Delta t)/2
//          N_neum_to_diri   - number of dofs which should become finally
//                             Dirichlet nodes, MUST BE ORDERED !!!
//          neum_to_diri     - array containing the indices of the dofs which
//                             should become Dirichlet nodes
//          neum_to_diri_bdry- array containing the number of the boundary part
//          neum_to_diri_param - array containing the parameter of the boundary part
//          compute_matrix_D - flag which says if to compute matrix_D_entries
//          BoundaryValue    - pointer to function containing the boundary values
//          BoundaryValues   - contains the boundary values if no pointer is specified
//
// output : matrix_D_Entries - if compute_matrix_D is true
//          B                - array for right hand side for the solver
//          tilde_u          - if compute_matrix_D is true
//
/*******************************************************************************/

#ifdef __2D__
void FEM_FCT_ForConvDiff(TSquareMatrix2D *M_C,TSquareMatrix2D *A,
			 int N_U, int N_Active,
			 double *lump_mass, double *matrix_D_Entries,
			 double *sol, double *oldsol,
			 double *B, double *rhs, double *rhs_old,
			 double *tilde_u,
			 int N_neum_to_diri, int *neum_to_diri,
			 int *neum_to_diri_bdry,
			 double *neum_to_diri_param,
			 int compute_matrix_D,
			 BoundValueFunct2D *BoundaryValue,
			 double *BoundaryValues)
#endif
#ifdef __3D__
void FEM_FCT_ForConvDiff(TSquareMatrix3D *M_C,TSquareMatrix3D *A,
			 int N_U, int N_Active,
			 double *lump_mass, double *matrix_D_Entries,
			 double *sol, double *oldsol,
			 double *B, double *rhs, double *rhs_old,
			 double *tilde_u,
			 int N_neum_to_diri, int *neum_to_diri,
			 double *neum_to_diri_x, double *neum_to_diri_y, double *neum_to_diri_z,
			 int compute_matrix_D, 
			 BoundValueFunct3D *BoundaryValue,
			 double *BoundaryValues)
#endif
{
  int *ColInd, *RowPtr, N_Entries, *ColInd_M, *RowPtr_M, N_Entries_M;
  int i,j,j0,j1,j2,j3,jj,index,linear;
  int solver_param[3];
  double *Entries,*Entries_M, *res, *aux_v, *alpha, eps = 1e-10, help;
  double *P_plus, *P_minus, *Q_plus, *Q_minus, *R_plus, *R_minus, val, val1;
  double solver_param_d[1];
  double delta_t = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta2 = TDatabase::TimeDB->THETA2;
  double theta3 = TDatabase::TimeDB->THETA3;
  double theta4 = TDatabase::TimeDB->THETA4;

  /*
  if (N_Active < N_U)
  {
    OutPut("N_Active < N_U ("<< N_Active<< "<"<< N_U <<
      ") !!! FOR APPLYING ALGEBRAIC FLUX CORRECTION, THE BOUNDARY CONDITIONS SHOULD BE NEUMANN !!!" << endl);
      
#ifdef _MPI      
    MPI_Finalize(); 
#endif         
    exit(4711);
  }
  */

  // get pointers to columns, rows and entries of matrix M_C
  ColInd_M = M_C->GetKCol();
  RowPtr_M = M_C->GetRowPtr();
  Entries_M = M_C->GetEntries();
  N_Entries_M = M_C->GetN_Entries();

  // get pointers to columns, rows and entries of matrix A
  ColInd = A->GetKCol();
  RowPtr = A->GetRowPtr();
  Entries = A->GetEntries();
  N_Entries = A->GetN_Entries();

  // allocate memory and set pointers
  res = new double[2*N_Entries+7*N_U];
  memset(res, 0, (2*N_Entries+7*N_U)*SizeOfDouble);
  alpha = res + N_Entries;
  P_plus = alpha + N_Entries;
  P_minus = P_plus + N_U;
  Q_plus = P_minus + N_U;
  Q_minus = Q_plus + N_U;
  R_plus = Q_minus + N_U;
  R_minus = R_plus + N_U;
  aux_v = R_minus + N_U;

  //compute the auxiliary solution
  // D has to be computed only once
  if (compute_matrix_D)
  {
    memset(matrix_D_Entries , 0, N_Entries*SizeOfDouble);
    // loop over all rows
    for(i=0;i<N_U;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      // compute first the matrix D
      for(j=j0;j<j1;j++)
      {
        // column
        index = ColInd[j];
        // only the active part of the matrix
        //if (index>=N_U)
        //	  continue;
        // only off-diagonals
        if (index!=i)
        {
          if (Entries[j] > 0)
            matrix_D_Entries[j] = -Entries[j];
          // now check the transposed entry
          j2 = RowPtr[index];
          j3 = RowPtr[index+1];
          for (jj=j2;jj<j3;jj++)
          {
            if (ColInd[jj]==i)
            {
              if (-Entries[jj]<matrix_D_Entries[j])
                matrix_D_Entries[j] = -Entries[jj];
              break;
            }
          }
        }
      }
    }

    // compute diagonal entry of D
    // loop over all rows
    for(i=0;i<N_U;i++)
    {
      // i-th row of sqmatrix
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      val = 0.0;
      // add all entries of i-th row
      for(j=j0;j<j1;j++)
      {
        val +=  matrix_D_Entries[j];
        index = ColInd[j];
        if (index==i)
          jj = j;
      }
      matrix_D_Entries[jj] = -val;
    }
    // matrix D is computed
    // add this matrix to A giving A+D=L
    Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);
  }

//   OutPut("entries " << Ddot(N_Entries,Entries,Entries) << 
//   	 " lump " << Ddot(N_U,lump_mass,lump_mass) << endl);
//   exit(0);
  
  //===================================================
  // different schemes
  //===================================================
  
  switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
  {
    // simplest linear scheme, from the preprint of Kuzmin (2008)
    // linear schemes only for Crank-Nicolson
    case 1:
      // L u_{k-1}
      MatVect(A,oldsol,aux_v);

      // M_L^(-1)(f_{k-1} - L u_{k-1})
      for(i=0;i<N_U;i++)
      {                                           //OutPut("oldsol"<<oldsol[i] << endl );//auxv ist falsch an dieser Stelle oldsol auch
	  //OutPut(i << " lump " << lump_mass[i] <<" aux "<< aux_v[i]<< " rhs_old: " << rhs_old[i]<< endl );
        aux_v[i] = (rhs_old[i] - aux_v[i])/lump_mass[i];
        //OutPut(i << "  lump " << lump_mass[i] <<" aux "<< aux_v[i]<< " rhs_old: " << rhs_old[i]<< 
	//       " " << rhs[i] << endl );
        //OutPut(i << " lump " << lump_mass[i] << endl);
      }
      //OutPut("auxv0 " << Ddot(N_U,aux_v,aux_v) <<
      // " " << Ddot(N_U,lump_mass,lump_mass) <<
      // " " <<  Ddot(N_U,rhs_old,rhs_old) <<
      // " " <<  Ddot(N_U,oldsol,oldsol) << endl);
      val = delta_t/2.0;
      for(i=0;i<N_U;i++)
      {
        tilde_u[i] = oldsol[i] + val * aux_v[i];
      }

      TDatabase::TimeDB->CURRENTTIME -= val;
      // set correct boundary conditions for intermediate solution
      for (j=0 ;j<N_neum_to_diri;j++)
      {
        i=neum_to_diri[j];
#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
#endif
#ifdef __3D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(0, neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
#endif
        aux_v[i] = 2*(tilde_u[i] - oldsol[i])/delta_t;
      }
      TDatabase::TimeDB->CURRENTTIME += val;
      // compute matrix res
      // loop over all rows
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];
          if (index==i)
            continue;
          val = - matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]);
          res[j]= Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
          // prelimiting with MINMOD
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
            res[j] = delta_t*MinMod(res[j]/delta_t, val);
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
            if (res[j]*(tilde_u[index]-tilde_u[i])>0)
              res[j] = 0;
          }
          res[j] *= delta_t;
        }
      }
      
  
      
      break;
    case 2:
    case 3:
    case 4:
      if (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME > 2)
      {
        theta1 = 0.5;
        theta2 = 1-theta1;
        // compute auxiliary solution, eq. (55) from [Kuz09]
        // rhs comes to aux_v
        // L u_{k-1}, (note: A = -L)
        MatVect(A,oldsol,aux_v);
        // delta_t Lu_{k-1}/2
        val = -theta2 * delta_t;
        Dscal(N_U, val, aux_v);
      }
      val = theta1 * delta_t;
      for(i=0;i<N_U;i++)
      {
        // (M  + delta_t L/2) u_{k-1}
        aux_v[i] += oldsol[i]*lump_mass[i];
        // (M  + delta_t L/2) u_{k-1} + delta_t/2 rhs_old
        aux_v[i] += val * rhs_old[i];
        // (M  + delta_t L/2) u_{k-1} + delta_t/2 rhs_old + delta_t/2 rhs
        aux_v[i] += val * rhs[i];
      }

      // matrix comes to A
      for (i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        // j-th column
        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd[j];
          if (index==i)
          {
            Entries[j] = lump_mass[i] + val * Entries[j];
          }
          else
          {
            Entries[j] *= val;
          }
        }
      }

      // prepare solver
      // use bicgstab with ssor preconditioner
      solver_param[0] = TDatabase::ParamDB->SC_SOLVER_SCALAR;
      TDatabase::ParamDB->SC_SOLVER_SCALAR = 13;
      solver_param[1] = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 3;
      solver_param[2] = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100;
      solver_param_d[0] = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-15;
      Solver(A, aux_v, tilde_u);

      // reset solver
      TDatabase::ParamDB->SC_SOLVER_SCALAR = solver_param[0];
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = solver_param[1];
      TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = solver_param[2];
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = solver_param_d[0];

      // reset matrix
      for (i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        // j-th column
        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd[j];
          if (index==i)
          {
            Entries[j] = (Entries[j] - lump_mass[i])/val;
          }
          else
          {
            Entries[j] /= val;
          }
        }
      }

      val = delta_t;
      TDatabase::TimeDB->CURRENTTIME -= val;
      // set correct boundary conditions for auxiliary solution
      for (j=0 ;j<N_neum_to_diri;j++)
      {
        i=neum_to_diri[j];
#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
#endif
#ifdef __3D__
	if (BoundaryValue!=NULL)
	{
	    BoundaryValue(0, neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
			  tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
#endif
      }

      TDatabase::TimeDB->CURRENTTIME += val;

      // intermediate solution should be nonnegative
      for(i=0;i<N_U;i++)
      {
       if (tilde_u[i] < 0)
        OutPut("int sol neg " << tilde_u[i] << endl);
      }

      switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
      {
        case 2:
          for(i=0;i<N_U;i++)
          {
            aux_v[i] = tilde_u[i] - oldsol[i];
	  }
          break;
        case 3:
          // eq. (60) from [Kuz09]
          // L tilde_u
          MatVect(A,tilde_u,aux_v);
          // approximation of derivative of tilde_u
          for(i=0;i<N_U;i++)
          {
            aux_v[i] = -aux_v[i]/lump_mass[i];
            //if (delta_t * aux_v[i] + tilde_u[i]  < -1e-13)
            //OutPut(delta_t * aux_v[i] + tilde_u[i] << " " << endl);
          }
          break;
        case 4:
          // restore matrix A
          Daxpy(N_Entries, -1.0, matrix_D_Entries, Entries);
          // A tilde_u
          MatVect(A, tilde_u, R_minus);
          // restore matrix D
          Daxpy(N_Entries, 1.0, matrix_D_Entries, Entries);
          // prepare solver
          // use bicgstab with ssor preconditioner
          solver_param[0] = TDatabase::ParamDB->SC_SOLVER_SCALAR;
          TDatabase::ParamDB->SC_SOLVER_SCALAR = 13;
          solver_param[1] = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
          TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = 3;
          solver_param[2] = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
          TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = 100;
          solver_param_d[0] = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
          TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-11;

          Solver(M_C, R_minus, aux_v);

          // reset solver
          TDatabase::ParamDB->SC_SOLVER_SCALAR = solver_param[0];
          TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = solver_param[1];
          TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR = solver_param[2];
          TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = solver_param_d[0];
          memset(R_minus, 0 , N_U*SizeOfDouble);
          break;
      }

      // HAS TO BE CHECKED
      /*for (j=0 ;j<N_neum_to_diri;j++)
      {
          i=neum_to_diri[j];
          aux_v[i] = 0;
	  }*/
      // compute matrix res, eq. (57) in [Kuz09]
      // loop over all rows
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];
          if (index==i)
            continue;
          switch (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME)
          {
            case 2:
              val = - matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]);
	      val *= delta_t;
              res[j]= 2 * Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
	      break;
            case 3:
            case 4:
              theta1 = TDatabase::TimeDB->THETA1;
              theta2 = TDatabase::TimeDB->THETA2;
              val = -theta1*tilde_u[i] - theta2*oldsol[i];
              val += theta1*tilde_u[index] + theta2*oldsol[index];
              val *= matrix_D_Entries[j];
              res[j] = Entries_M[j] * (aux_v[i]-aux_v[index]) + val;
              res[j] *= delta_t;
              /*val = -theta1*tilde_u[i] - theta2*oldsol[i];
              val += theta1*tilde_u[index] + theta2*oldsol[index];
              val *= delta_t * matrix_D_Entries[j];
              res[j] = Entries_M[j] *(tilde_u[i] - oldsol[i]);
              res[j] += Entries_M[j] *(tilde_u[index] - oldsol[index]);
              res[j] += val;*/
              break;
          }
          //prelimiting with MINMOD
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
              (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
             val =  -matrix_D_Entries[j]*(tilde_u[i]-tilde_u[index]);
              res[j] = delta_t*MinMod(res[j]/delta_t, val);
          }
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
              (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
	    if (res[j]*(tilde_u[index]-tilde_u[i])>0)
	      res[j] = 0;
	  }
        }
      }
      break;
    case 0:
      //=======================================================
      //nonlinear scheme
      //=======================================================
      if (compute_matrix_D)
      {
	  // reset parameter if necessary
	  if (TDatabase::ParamDB->FEM_FCT_PRELIMITING==0)
	  {
	      TDatabase::ParamDB->FEM_FCT_PRELIMITING = 2;
	      OutPut("nonlinear FEM-FCT works badly without prelimiting, "
		     << " set ParamDB->FEM_FCT_PRELIMITING to "
		     << TDatabase::ParamDB->FEM_FCT_PRELIMITING << endl);
	  }
        // L u_{k-1}
        MatVect(A,oldsol,aux_v);

        // M_L^(-1)(f_{k-1} - L u_{k-1})
        val = delta_t/2.0;
        for(i=0;i<N_U;i++)
        {
          val1 = (rhs_old[i] - aux_v[i])/lump_mass[i];
          tilde_u[i] = oldsol[i] + val * val1;
        }
        //auxiliary variable is computed
        TDatabase::TimeDB->CURRENTTIME -= val;
        // set correct boundary conditions for intermediate solution
        for (j=0;j<N_neum_to_diri;j++)
        {
          i=neum_to_diri[j];
#ifdef __2D__
	if (BoundaryValue!=NULL)
	{
          BoundaryValue(neum_to_diri_bdry[j], neum_to_diri_param[j],
            tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}
	
          //aux_v[i] = 2*(tilde_u[i] - oldsol[i])/delta_t;
#endif
#ifdef __3D__
          // WRONG
	if (BoundaryValue!=NULL)
	{
          BoundaryValue(0, neum_to_diri_x[j], neum_to_diri_y[j], neum_to_diri_z[j],
            tilde_u[i]);
	}
	else
	{
	    tilde_u[i] =   BoundaryValues[j];
	}	    
#endif
        }
        TDatabase::TimeDB->CURRENTTIME += val;
      }
      // compute matrix res
      // loop over all rows
      val = theta1 * delta_t;
      val1 = theta2 * delta_t;
      for(i=0;i<N_U;i++)
      {
        // i-th row of sqmatrix
        j0 = RowPtr_M[i];
        j1 = RowPtr_M[i+1];

        for(j=j0;j<j1;j++)
        {
          // column
          index = ColInd_M[j];

          res[j]= -Entries_M[j] * ((sol[index]-sol[i]) - (oldsol[index] - oldsol[i]))
            + val * matrix_D_Entries[j] * (sol[index] - sol[i])
            + val1 * matrix_D_Entries[j] * (oldsol[index]- oldsol[i]);
          // prelimiting with MINMOD, THIS HAS TO BE CHECKED
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==1)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
            res[j] = delta_t*MinMod(res[j]/delta_t, -matrix_D_Entries[j] * (tilde_u[i]-tilde_u[index]));
          // more prelimiting
          if ((TDatabase::ParamDB->FEM_FCT_PRELIMITING==2)||
            (TDatabase::ParamDB->FEM_FCT_PRELIMITING==3))
          {
            if (res[j]*(tilde_u[index]-tilde_u[i])>0)
              res[j] = 0;
          }
        }
      }
      break;
    default:
      OutPut("FEM_FCT_LINEAR_TYPE " <<
        TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE << " DOES NOT EXIST !!!"<<endl);
      exit(4711);
  }
  /*
    double t_min = 0, t_max = 0;
    for(i=0;i<N_U;i++)
    {
        //if (tilde_u[i] < 0)
  //	  OutPut("tilde " << i << " "  << tilde_u[i] << endl);
        if (tilde_u[i]  < t_min)
      t_min = tilde_u[i];
        if (tilde_u[i]  > t_max)
      t_max = tilde_u[i];
    }
  OutPut("minmax " << t_min << "  " << t_max <<endl);
  */
// exit(0);

  //=======================================================
  //Zalesaks limiter
  //=======================================================

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr_M[i];
    j1 = RowPtr_M[i+1];
    for(j=j0;j<j1;j++)
    {
      index = ColInd_M[j];
      // only the active part of the matrix
      if(i != index)
      {
        // if (i==1111)
        //  OutPut(i << " " << index << " " << res[j] << " : ");
        if (res[j] > 0.0)
          P_plus[i] += res[j];
        if (res[j] <= 0.0)
          P_minus[i] += res[j];

        help = tilde_u[index]-tilde_u[i];

        // if (i==1111)
        //    OutPut(tilde_u[index] << " " << tilde_u[i] << " " << help << endl);
        if (help > Q_plus[i])
          Q_plus[i]= help;

        if (help < Q_minus[i])
          Q_minus[i]= help;
      }
    }                                             // end loop j
  }
  //for(i=0;i<N_U;i++)
  //    OutPut(i << " "<< P_plus[i] << " " <<  P_minus[i] << " " << Q_plus[i] << " " <<  Q_minus[i] << endl);
  //exit(1);
  /*
    for(i=0;i<N_U;i++)
    {
        if (tilde_u[i] + Q_minus[i] < t_min)
      OutPut("tilde1 " << i << " "  << tilde_u[i] + Q_minus[i] <<
       " " << tilde_u[i]  << " "  <<   Q_minus[i] <<
       " " << t_min << endl);
    }
  */

  for(i=0;i<N_U;i++)
  {
    if(fabs(P_plus[i]) == 0.0)
      R_plus[i] = 1.0;
    else
    {
      //OutPut(Q_plus[i] << " "  << P_plus[i] << " " << lump_mass[i] << " ");
      help = (lump_mass[i] * Q_plus[i])/P_plus[i];
      if(help < 1.0)
        R_plus[i] = help;
      else
        R_plus[i] = 1.0;
    }
    if(fabs(P_minus[i]) == 0.0)
      R_minus[i] = 1.0;
    else
    {
      //OutPut(Q_minus[i] << " "  << P_minus[i] << " " << lump_mass[i] << " ");
      help = (lump_mass[i] * Q_minus[i])/P_minus[i];
      if(help < 1.0)
        R_minus[i] = help;
      else
        R_minus[i] = 1.0;
    }
  }

  // treat Dirichlet nodes
  for (j=0;j<N_neum_to_diri;j++)
  {
    i = neum_to_diri[j];
    R_plus[i] = 1;
    R_minus[i] = 1;
  }

  for(i=0;i<N_U;i++)
  {
    // i-th row of sqmatrix
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];

    for(j=j0;j<j1;j++)
    {
      index = ColInd[j];
      if(res[j] > 0.0)
      {
        //Initialisation
        alpha[j] = R_plus[i];
        if(alpha[j] > R_minus[index])
          alpha[j] = R_minus[index];
      }
      if(res[j]<=0.0)
      {
        //initialisation
        alpha[j] = R_minus[i];
        if(alpha[j] > R_plus[index])
          alpha[j] = R_plus[index];
      }
      // clipping, see Kuzmin (2009), end of Section 5
      //if ((fabs(Q_plus[i])< eps)&&(fabs(Q_minus[i])< eps))
      //  alpha[j] = 0;
    }                                             //end loop j
  }

  //=======================================================
  //correct right hand side
  //=======================================================

  // Crank--Nicolson
  if ((fabs(theta2-0.5)<eps) && (fabs(theta3-0.5)<eps)&&
    (TDatabase::ParamDB->INTERNAL_LINEAR_SCHEME<2))
  {
    for(i=0;i<N_U;i++)
    {
      val=0.0;
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
        index = ColInd[j];
        if(i != index)
        {
          val +=  alpha[j] * res[j];
        }
      }
      B[i] = tilde_u[i] * lump_mass[i] + theta4*delta_t*rhs[i] + val;
    }
  }
  else
  {
    // L u_{k-1}
    MatVect(A,oldsol,aux_v);
    val1 = theta2*delta_t;
    for(i=0;i<N_U;i++)
    {
      val =0.0;
      j0 = RowPtr[i];
      j1 = RowPtr[i+1];
      for(j=j0;j<j1;j++)
      {
        index = ColInd[j];
        if(i != index)
        {
          val +=  alpha[j] * res[j];
        }
      }
      B[i] = lump_mass[i] * oldsol[i] - val1 * aux_v[i] + delta_t * theta3 * rhs_old[i]
        + delta_t * theta4*rhs[i] + val;
      /* if (B[i] < -1e-10)
      {
          OutPut(i << " " << B[i] <<  " " << lump_mass[i] * oldsol[i]  << " " <<
           -val1 * aux_v[i] << " " <<
           delta_t * theta3 * rhs_old[i] << " " <<
           delta_t * theta4*rhs[i] << " " <<  val << endl);
           }*/
    }
  }

  //OutPut(" vecB " << Ddot(N_U,B,B) << endl);
  delete [] res;
}
