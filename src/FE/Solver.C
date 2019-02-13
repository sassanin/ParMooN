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
// @(#)Solver.C        1.16 06/27/00
//
// Purpose:     solve equation system
//
// Author:      Gunar Matthies (17.08.98)
//
// History:     start of implementation 17.08.98 (Gunar Matthies)
//
// =======================================================================
#include <Database.h>
#include <Solver.h>
#include <MooNMD_Io.h>
#include <LinAlg.h>

#include <string.h>
#include <stdlib.h>
#ifdef __MAC64__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

#include <ItMethod.h>
#include <FEDatabase2D.h>
#include <FixedPointIte.h>
#include <FgmresIte.h>
#include <Cg.h>
#include <Bcgs.h>
#include <JacobiIte.h>
#include <SSORIte.h>
#include <MultiGridScaIte.h>
#include <DirectSolver.h>

extern "C"
{
  #include <amg_solve_main.h>
}


/******************************************************************************/
// Solver
// solves linear system
// output in sol
// scalar problems: ns_type == 0
/******************************************************************************/

void Solver(TSquareMatrix **sqmatrices, TMatrix **matrices,
double *rhs, double *sol,
MatVecProc *MatVect, DefectProc *Defect,
TMultiGrid2D *MG,
int N_Unknowns, int ns_type
#ifdef _MPI   
     ,TParFECommunicator3D *ParComm
#endif        
)
{
  int solver_type, prec_type;
  double *itmethod_sol, *itmethod_rhs, t1, t2;
  TItMethod *itmethod, *prec;
  TSquareMatrix2D *A;

  t1 = GetTime();
//  memset(sol, 0, N_Unknowns * SizeOfDouble);

  if (!ns_type)
  {
    prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
    solver_type = TDatabase::ParamDB->SC_SOLVER_SCALAR;
  }
  else
  {
    prec_type = TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE;
    solver_type = TDatabase::ParamDB->SC_SOLVER_SADDLE;
  }

  // extract the matrices
  switch(ns_type)
  {
    case 0:
      A = (TSquareMatrix2D*) sqmatrices[0];
      break;
  }

  // choose the solver
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
    case DIRECT:
      switch(ns_type)
      {
        case 0:
        //  OutPut("DIRECT!");
          DirectSolver(A, rhs, sol);
          break;
      }
      //TDatabase::ParamDB->SOLVER_TYPE = 1;
      break;

    case AMG_SOLVE:
      switch(ns_type)
      {
        case 0:
          Solver(A, rhs, sol);
          break;
      }
      break;

    case GMG:
      switch (prec_type)
      {
        case 1:
          prec = new TJacobiIte(MatVect, Defect, NULL,
            0, N_Unknowns, 1
#ifdef _MPI   
                  , ParComm
#endif    
    
                   );
          break;
        case 3:
          prec = new TSSORIte(MatVect, Defect, NULL,
            0, N_Unknowns, 1);
          break;
        case 5:
#ifdef __2D__
          prec = new TMultiGridScaIte(MatVect, Defect, NULL,
              0, N_Unknowns, MG, 0);
#else
          prec = new TMultiGridScaIte(MatVect, Defect, NULL,
              0, N_Unknowns, NULL, 0);
#endif
          break;
        default:
          OutPut("Unknown preconditioner !!!" << endl);
          exit(4711);
      }
      switch (solver_type)
      {
        case 11:
          itmethod = new TFixedPointIte(MatVect, Defect, prec,
            0, N_Unknowns, 1);
          if (prec_type == 5)
          {
            itmethod_sol = new double[N_Unknowns];
            itmethod_rhs = new double[N_Unknowns];
            memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
            memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
          }
          else
          {
            itmethod_sol = sol;
            itmethod_rhs = rhs;
          }
          break;
        case 12:
    itmethod = new TCg(MatVect, Defect, prec,
           0, N_Unknowns, 1);
          if (prec_type == 5)
          {
            itmethod_sol = new double[N_Unknowns];
            itmethod_rhs = new double[N_Unknowns];
            memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
            memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
          }
          else
          {
            itmethod_sol = sol;
            itmethod_rhs = rhs;
          }
          break;

        case 13:
          itmethod = new TBcgs(MatVect, Defect, prec,
                      0, N_Unknowns, 1);
          if (prec_type == 5)
          {
            itmethod_sol = new double[N_Unknowns];
            itmethod_rhs = new double[N_Unknowns];
            memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
            memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
          }
          else
          {
            itmethod_sol = sol;
            itmethod_rhs = rhs;
          }
          break;

        case 16:
          itmethod = new TFgmresIte(MatVect, Defect, prec,
            0, N_Unknowns, 1);
          if (prec_type == 5)
          {
            itmethod_sol = new double[N_Unknowns];
            itmethod_rhs = new double[N_Unknowns];
            memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
            memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
          }
          else
          {
            itmethod_sol = sol;
            itmethod_rhs = rhs;
          }
          break;
        default:
          OutPut("Unknown solver !!!" << endl);
          exit(4711);
      }
      // solve linear system
      itmethod->Iterate(sqmatrices,NULL,itmethod_sol,itmethod_rhs);

      delete prec;
      delete itmethod;

      switch (solver_type)
      {
        case 11:
          if (prec_type == 5)
          {
            memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
            memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            delete itmethod_sol;
            delete itmethod_rhs;
          }
          break;
        case 12:
        case 13:
        case 16:
          if (prec_type == 5)
          {
            memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
            memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            delete itmethod_sol;
            delete itmethod_rhs;
          }
          break;
      }
      break;
  }
  t2 = GetTime();
  if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut("time for solving: " << t2-t1 << endl);
}


void SetAMGDefaults(AMG_CoarsenContext &cc, AMG_SolverContext &sc)
{
  // THESE ARE THE DEFAULTS, DO NOT CHANGE
  // coarsen context
  cc.alpha = TDatabase::ParamDB->CC_ALPHA;
  cc.beta = TDatabase::ParamDB->CC_BETA;
  cc.mincluster = TDatabase::ParamDB->CC_MINCLUSTER;
  cc.maxcluster = TDatabase::ParamDB->CC_MAXCLUSTER;
  cc.maxdistance = TDatabase::ParamDB->CC_MAXDISTANCE;
  cc.maxconnectivity = TDatabase::ParamDB->CC_MAXCONNECTIVITY;
  cc.depthtarget = TDatabase::ParamDB->CC_DEPTHTARGET;
  cc.coarsentarget = TDatabase::ParamDB->CC_COARSENTARGET;
  cc.coarsenrate = TDatabase::ParamDB->CC_COARSENRATE;
  cc.major = TDatabase::ParamDB->CC_MAJOR;
  cc.dependency = TDatabase::ParamDB->CC_DEPENDENCY;
  cc.verbose = TDatabase::ParamDB->CC_VERBOSE;

  // THESE ARE THE DEFAULTS, DO NOT CHANGE
  // solver context
  sc.solver = TDatabase::ParamDB->SC_SOLVER_SCALAR;
  sc.preconditioner = TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
  sc.maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
  sc.red_factor = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
  sc.res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
  sc.amg_prec_it = TDatabase::ParamDB->SC_AMG_PREC_IT;
  sc.amg_prec_red_factor = TDatabase::ParamDB->SC_AMG_PREC_RED_FACTOR;
  sc.ex_maxit = TDatabase::ParamDB->SC_EX_MAXIT;
  sc.gmres_restart = TDatabase::ParamDB->SC_GMRES_RESTART;
  sc.lcd_start_vector = TDatabase::ParamDB->SC_LCD_START_VECTOR;
  sc.ilu_beta = TDatabase::ParamDB->SC_ILU_BETA;
  sc.sor_omega = TDatabase::ParamDB->SC_SOR_OMEGA;
  sc.gamma = TDatabase::ParamDB->SC_MG_CYCLE_SCALAR;
  sc.smoother = TDatabase::ParamDB->SC_SMOOTHER_SCALAR;
  sc.n1[0] = TDatabase::ParamDB-> SC_PRE_SMOOTH_SCALAR;
  sc.n2[0] = TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;
  sc.omega[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  sc.step_length_control_fine=TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR;
  sc.step_length_control_all = TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR;
  sc.smoother_red_factor = TDatabase::ParamDB->SC_SMOOTHER_RED_FACTOR;
  sc.coarse_smoother = TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR;
  sc.coarse_maxit = TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
  sc.coarse_red_factor = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR;
  sc.omega_coarse[0] = TDatabase::ParamDB->SC_OMEGA_COARSE_0;
  sc.omega_p[0] = TDatabase::ParamDB->SC_OMEGA_P_0;
  sc.schur_inv_of_A = TDatabase::ParamDB->SC_SCHUR_INV_OF_A;
  sc.schur_inv_of_A_maxit = TDatabase::ParamDB->SC_SCHUR_INV_OF_A_MAXIT;
  sc.schur_iteration_damp = TDatabase::ParamDB->SC_SCHUR_ITERATION_DAMP;
  sc.schur_iteration_maxit = TDatabase::ParamDB->SC_SCHUR_ITERATION_MAXIT;
  sc.schur_step_length_control=TDatabase::ParamDB->SC_SCHUR_STEP_LENGTH_CONTROL;
  sc.ilut_tol = TDatabase::ParamDB->SC_ILUT_TOL;
  sc.ilut_absolute_fillin = TDatabase::ParamDB->SC_ILUT_ABSOLUTE_FILLIN;
  sc.ilut_relative_fillin = TDatabase::ParamDB->SC_ILUT_RELATIVE_FILLIN;
  sc.ilut_sort = TDatabase::ParamDB->SC_ILUT_SORT;
  sc.mixed_bcgs_cgs_switch_tol= TDatabase::ParamDB->SC_MIXED_BCGS_CGS_SWITCH_TOL;
  sc.div_factor= TDatabase::ParamDB->SC_DIV_FACTOR;
  sc.smoothing_steps=TDatabase::ParamDB->SC_SMOOTHING_STEPS;
  sc.n1_param=TDatabase::ParamDB->SC_N1_PARAM;
  sc.n2_param=TDatabase::ParamDB->SC_N2_PARAM;
  sc.minit=TDatabase::ParamDB->SC_MINIT;
  sc.vas_laz_delta=TDatabase::ParamDB->SC_VAS_LAZ_DELTA;
  sc.row_equilibration=TDatabase::ParamDB->SC_ROW_EQUILIBRATION;
  sc.braess_sarazin_matrix=TDatabase::ParamDB->SC_BRAESS_SARAZIN_MATRIX;
  sc.braess_sarazin_alpha=TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA;
  sc.verbose=TDatabase::ParamDB->SC_VERBOSE_AMG;
}


/*******************************************************************/
/*        SCALAR PROBLEMS                                          */
/*******************************************************************/
void Solver(TSquareMatrix *matrix, double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n,m,i,j, k, *changes;
  double tmp;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;
 info = mstats();
memory[0]=memory[1]=memory[2]=0.;
#else  
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  
  sc.system_type=SCALAR;
  n = matrix->GetN_Columns();
  m = matrix->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=1;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = matrix->GetRowPtr();
  A0->ja = matrix->GetKCol();
  A0->a = matrix->GetEntries();

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  if (TDatabase::ParamDB->INTERNAL_SORT_AMG)
  {
      // sort matrix, diagonal entry
      if (TDatabase::ParamDB->SC_VERBOSE>1)
      OutPut("AMG solver - reordering matrix structure, diagonal first"<< endl);
      changes = new int[n];
      for (i=0;i<n;i++)
      {                              // first entry in row
	  for (j= A0->ra[i];j<A0->ra[i+1];j++)
	  {
	      if (A0->ja[j]==i)
	      {
		  // save information about the changes
		  changes[i] = j;
		  k = A0->ra[i];
		  // set non diagonal entry of ja
		  A0->ja[j] =A0->ja[k];
		  // diagonal entry gets number of row entries
		  A0->ja[k] =A0->ra[i+1]-A0->ra[i];
		  // change matrix entries
		  tmp = A0->a[j];
		  A0->a[j] = A0->a[k];
		  A0->a[k] =tmp;
		  break;
	      }                          // endif
	  }                            // endfor j
      }                              // endfor i
  }
  // fill vector structures
  b->n=n;
  b->b=1;
  b->x=rhs;

  x->n=n;
  x->b=1;
  x->x=sol;
  /*for(i=0;i<n;i++)
  {
   OutPut("sol("<<i+1<< ")=" << sol[i] << ";" << endl);
    OutPut("rhs("<<i+1<< ")=" << rhs[i] << ";" << endl);
    }*/

  SetAMGDefaults(cc, sc);
  
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else   
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  OutPut("call AMG main"<< endl);
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else    
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

  // redo reordering of the matrix, default
  if (TDatabase::ParamDB->INTERNAL_SORT_AMG)
  {
      if (TDatabase::ParamDB->SC_VERBOSE>1)
      OutPut("AMG solver - reordering matrix structure redone"<< endl);
      for (i=0;i<n;i++)
      {                              // first entry in row
	  j = A0->ra[i];
	  k = changes[i];
	  A0->ja[j] = A0->ja[k];
	  A0->ja[k] = i;
	  tmp = A0->a[j];
	  A0->a[j] = A0->a[k];
	  A0->a[k] = tmp;
      }
  }

  delete A0;
  delete B[0];
  delete b;
  delete x;
  delete changes;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/*        SCALAR PROBLEMS WITH VARIOUS RIGHT HAND SIDES            */
/*******************************************************************/
void Solver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n,m,i,j,k, *changes;
  double tmp;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=0.;
#else    
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  
  switch(N_Rhs)
  {
    case 1:
      sc.system_type=SCALAR;
      break;
    case 2:
      sc.system_type=SCALAR2;
      break;
    case 3:
      sc.system_type=SCALAR3;
      break;
    case 6:
      sc.system_type=SCALAR6;
      break;
    default:
      OutPut("SYSTEM_TYPE NOT IMPLEMENTED !!!"<< endl);
      exit(4711);
  }
  n = matrix->GetN_Columns();
  m = matrix->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;
  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=N_Rhs;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = matrix->GetRowPtr();
  A0->ja = matrix->GetKCol();
  A0->a = matrix->GetEntries();

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  // sort matrix, diagonal entry
  if (TDatabase::ParamDB->SC_VERBOSE>1)
  OutPut("AMG solver - reordering matrix structure, diagonal first"<< endl);
  changes = new int[n];

  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
	  // save information about the changes
	  changes[i] = j;
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  // fill vector structures
  b->n=N_Rhs*n;
  b->b=1;
  b->x=rhs;

  x->n=N_Rhs*n;
  x->b=1;
  x->x=sol;
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  
  AMG(&sc,&cc,A0,B,x,b);
  
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

 // redo reordering of the matrix
  for (i=0;i<n;i++)
  {                              // first entry in row
      j = A0->ra[i];
      k = changes[i];
      A0->ja[j] = A0->ja[k];
      A0->ja[k] = i;
      tmp = A0->a[j];
      A0->a[j] = A0->a[k];
      A0->a[k] = tmp;
  }

  /*
    for (i=0;i<n;i++)
      cout << i << " " << sol[i] << " " <<  rhs[i] << endl;
  */

  delete changes;
  delete A0;
  delete B[0];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


void Solver(int what, TSquareMatrix *matrix, double *rhs, double *sol,
int N_Rhs)
{
  static AMG_CoarsenContext cc;
  static AMG_SolverContext sc;
  static AMG_MATRIX *A0,*B[2];
  static AMG_VECTOR *x,*b;
  static int n,m,i,j;
  static double tmp;

  switch(what)
  {
    case 0:                      // initialize
      switch(N_Rhs)
      {
        case 1:
          sc.system_type=SCALAR;
          break;
        case 2:
          sc.system_type=SCALAR2;
          break;
        case 3:
          sc.system_type=SCALAR3;
          break;
        case 6:
          sc.system_type=SCALAR6;
          break;
        default:
          OutPut("SYSTEM_TYPE NOT IMPLEMENTED !!!"<< endl);
          exit(4711);
      }
      n = matrix->GetN_Columns();
      m = matrix->GetN_Entries();

      A0 = new AMG_MATRIX;       // allocate matrix structure
      B[0] = new AMG_MATRIX;     // allocate matrix structure
      b = new AMG_VECTOR;        // allocate vector structures
      x = new AMG_VECTOR;
      // fill matrix structure
      A0->m = n;
      A0->n = n;
      A0->b = 1;
      A0->bb = 1;
      A0->system_as_scalar=1;
      A0->blocks_in_diag=N_Rhs;
      A0->nonzeros = m;
      A0->bandwidth=-1;
      A0->connections = m;
      A0->ra = matrix->GetRowPtr();
      A0->ja = matrix->GetKCol();
      A0->a = matrix->GetEntries();

      // fill matrix structure
      B[0]->m = n;
      B[0]->n = 0;
      B[0]->b = 1;
      B[0]->bb = 1;
      B[0]->system_as_scalar=1;
      B[0]->blocks_in_diag=1;
      B[0]->nonzeros = 0;
      B[0]->bandwidth=-1;
      B[0]->connections = 0;
      B[0]->ra = NULL;
      B[0]->ja = NULL;
      B[0]->a =  NULL;

      B[1] = B[0];
          
      // sort matrix, diagonal entry
      for (i=0;i<n;i++)
      {                          // first entry in row
	  // save information about the changes
        for (j= A0->ra[i];j<A0->ra[i+1];j++)
        {
          if (A0->ja[j]==i)
          {
            // set non diagonal entry of ja
            A0->ja[j] =A0->ja[A0->ra[i]];
            // diagonal entry gets number of row entries
            A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
            // change matrix entries
            tmp = A0->a[j];
            A0->a[j] = A0->a[A0->ra[i]];
            A0->a[A0->ra[i]] =tmp;
            break;
          }                      // endif
        }                        // endfor j
      }                          // endfor i

      // fill vector structures
      b->n=N_Rhs*n;
      b->b=1;
      b->x=rhs;

      x->n=N_Rhs*n;
      x->b=1;
      x->x=sol;
      SetAMGDefaults(cc, sc);
      AMG_Build(&sc,&cc,A0,B,x,b);
      break;

    case 1:                      // solve
      // fill vector structures
      b->n=N_Rhs*n;
      b->b=1;
      b->x=rhs;

      x->n=N_Rhs*n;
      x->b=1;
      x->x=sol;
      AMG_Solve(&sc,&cc,A0,B,x,b);
      break;

    case 2:                      // free memory
      AMG_Delete(&sc,&cc,A0,B,x,b);
      for (i=0;i<n;i++)          //reset array of columns
        A0->ja[A0->ra[i]] = i;

      delete A0;
      delete B[0];
      delete b;
      delete x;
      break;

    default:
      Error("Unknown type: " << what << endl);
  }                              // endswitch
}


/*******************************************************************/
/*        SCALAR MORTAR PROBLEMS                                   */
/*******************************************************************/

void Solver(TSquareMatrix *matrix, TMatrix *matrixmortar, double *rhs,
double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n, m, i, j;
  double *new_rhs, tmp;
  int memory[4];
  
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=0.;
#else   
  struct mallinfo MALLINFO;
  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type=SADDLE_1_TYPE_1;
  n = matrix->GetN_Columns();
  m = matrix->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0]= new AMG_MATRIX;          // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = n;                     // fill matrix structure
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar = 1;
  A0->blocks_in_diag = 1;
  A0->nonzeros = m;
  A0->bandwidth = -1;
  A0->connections = m;
  A0->ra = matrix->GetRowPtr();
  A0->ja = matrix->GetKCol();
  A0->a = matrix->GetEntries();

                                 // fill matrix structure
  B[0]->m = matrixmortar->GetN_Rows();
  B[0]->n = matrixmortar->GetN_Columns();
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar = 1;
  B[0]->blocks_in_diag = 1;
  B[0]->nonzeros =  matrixmortar->GetN_Entries();
  B[0]->bandwidth = -1;
  B[0]->connections = matrixmortar->GetN_Entries();
  B[0]->ra = matrixmortar->GetRowPtr();
  B[0]->ja = matrixmortar->GetKCol();
  B[0]->a = matrixmortar->GetEntries();

  B[1] = B[0];

  // sort matrix, diagonal entry
  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  new_rhs = new double[n+B[0]->m];
  for (i=0;i<n;i++)
    new_rhs[i] = rhs[i];

  for (i=n;i<n+B[0]->m;i++)
    new_rhs[i] = 0.0;

  b->n = n+B[0]->m;              // fill vector structures
  b->b = 1;
  b->x = new_rhs;
  x->n = n+B[0]->m;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
  double tmp2,tmp1;
  tmp2=tmp1=0;
  for (i=0;i<x->n;i++)
  {
    tmp2+=x->x[i];
    tmp1+=b->x[i];
  }
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc, &cc, A0, B, x, b);
  
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;
  delete A0;
  delete B[0];
  delete b;
  delete x;
  delete new_rhs;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3]- memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2]- memory[1] << endl;
}


/*******************************************************************/
/* CONNECT 2 BY 2 SYSTEM TO A SCALAR SYSTEM                        */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TMatrix *matrixa12,
TMatrix *matrixa21, TSquareMatrix *matrixa22,
double *rhs1, double *rhs2,
double *sol1, double *sol2)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n,m,i,j, *RowPtr, *ColPtr, n1, n2, *ARowPtr, *AKCol, *BRowPtr, *BKCol;
  double tmp, *Entries, *sol, *rhs, *AEntries, *BEntries;
  int memory[4], ja, ja1, ja2, jb, jb1, jb2, k;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  // connect blocks to single matrix
  n1 = matrixa11->GetN_Rows();
  n2 = matrixa21->GetN_Rows();
  n = n1 + n2;

  // compute total number of entries
  m = matrixa11->GetN_Entries()+ matrixa12->GetN_Entries()
    + matrixa21->GetN_Entries()+ matrixa22->GetN_Entries();

  // allocate new arrays
  RowPtr = new int[n+1];
  ColPtr = new int[m];
  Entries = new double[m];
  sol = new double[n];
  rhs = new double[n];

  // copy sol and rhs piece by piece
  for (i=0;i<n1;i++)
    sol[i] = sol1[i];
  for (i=n1;i<n;i++)
    sol[i] = sol2[i-n1];
  for (i=0;i<n1;i++)
    rhs[i] = rhs1[i];
  for (i=n1;i<n;i++)
    rhs[i] = rhs2[i-n1];

  // copy blocks in single matrix (upper part)
  ARowPtr = matrixa11->GetRowPtr();
  AKCol = matrixa11->GetKCol();
  AEntries =  matrixa11->GetEntries();
  BRowPtr = matrixa12->GetRowPtr();
  BKCol = matrixa12->GetKCol();
  BEntries =  matrixa12->GetEntries();

  RowPtr[0] = 0;
  // for the first n1 rows
  for (i=0;i<n1;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    RowPtr[i+1] =  RowPtr[i] + ja1 - ja + jb1 - jb;
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n1;
      Entries[k] = BEntries[jb+j];
    }
  }

  // copy blocks in single matrix (lower part)
  ARowPtr = matrixa21->GetRowPtr();
  AKCol = matrixa21->GetKCol();
  AEntries =  matrixa21->GetEntries();
  BRowPtr = matrixa22->GetRowPtr();
  BKCol = matrixa22->GetKCol();
  BEntries =  matrixa22->GetEntries();

  // for the second n2 rows
  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    //OutPut(ja << " " << ja1 << " " << jb << " " << jb1 << endl );
    RowPtr[i+1+n1] =  RowPtr[i+n1] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a21
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n1] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a22
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n1] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n1;
      Entries[k] = BEntries[jb+j];
    }
  }

  /*  for(i=0;i<n;i++)
    {
      for (j=RowPtr[i];j< RowPtr[i+1]; j++)
         OutPut("a(" << i+1 << ", " << ColPtr[j]+1 << ") =  " << Entries[j] << ";"<< endl);
      OutPut("sol("<<i+1<< ")=" << sol[i] << ";" << endl);
      OutPut("rhs("<<i+1<< ")=" << rhs[i] << ";" << endl);
    }
    exit(1);*/

  sc.system_type=SCALAR;

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=1;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = RowPtr;
  A0->ja = ColPtr;
  A0->a =  Entries;

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  // sort matrix, diagonal entry
  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i
  /*for(i=0;i<n;i++)
  {
    for (j=0;j< A0->ja[A0->ra[i]] ;j++)
    {
      OutPut("a0 " << A0->ra[i] << " " <<  A0->ja[A0->ra[i]+j] << " "
             << A0->a[A0->ra[i]+j] << endl);
    }
    }*/

  // fill vector structures
  b->n=n;
  b->b=1;
  b->x=rhs;

  x->n=n;
  x->b=1;
  x->x=sol;
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

  /*
    for (i=0;i<n;i++)
      cout << i << " " << sol[i] << " " <<  rhs[i] << endl;
  */

  // copy sol and rhs piece by piece
  for (i=0;i<n1;i++)
    sol1[i] = sol[i];
  for (i=n1;i<n;i++)
    sol2[i-n1] = sol[i];

  delete RowPtr;
  delete ColPtr;
  delete Entries;
  delete sol;
  delete rhs;
  delete A0;
  delete B[0];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* CONNECT 2 BY 2 SYSTEM TO A SCALAR SYSTEM                        */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12,
TSquareMatrix *matrixa21, TSquareMatrix *matrixa22,
double *rhs1, double *rhs2,
double *sol1, double *sol2)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n,m,i,j, *RowPtr, *ColPtr, n2, *ARowPtr, *AKCol, *BRowPtr, *BKCol;
  double tmp, *Entries, *sol, *rhs, *AEntries, *BEntries;
  int memory[4], ja, ja1, ja2, jb, jb1, jb2, k;

#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else 
   struct mallinfo MALLINFO;
  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  // connect blocks to single matrix
  // compute dimension (use that all blocks are squares)
  n2 = matrixa11->GetN_Columns();
  n = 2*n2;

  // compute total number of entries
  m = matrixa11->GetN_Entries()+ matrixa12->GetN_Entries()
    + matrixa21->GetN_Entries()+ matrixa22->GetN_Entries();

  // allocate new arrays
  RowPtr = new int[n+1];
  ColPtr = new int[m];
  Entries = new double[m];
  sol = new double[n];
  rhs = new double[n];

  // copy sol and rhs piece by piece
  for (i=0;i<n2;i++)
    sol[i] = sol1[i];
  for (i=n2;i<n;i++)
    sol[i] = sol2[i-n2];
  for (i=0;i<n2;i++)
    rhs[i] = rhs1[i];
  for (i=n2;i<n;i++)
    rhs[i] = rhs2[i-n2];

  // copy blocks in single matrix (upper part)
  ARowPtr = matrixa11->GetRowPtr();
  AKCol = matrixa11->GetKCol();
  AEntries =  matrixa11->GetEntries();
  BRowPtr = matrixa12->GetRowPtr();
  BKCol = matrixa12->GetKCol();
  BEntries =  matrixa12->GetEntries();

  RowPtr[0] = 0;
  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    RowPtr[i+1] =  RowPtr[i] + ja1 - ja + jb1 - jb;
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n2;
      Entries[k] = BEntries[jb+j];
    }
  }

  // copy blocks in single matrix (lower part)
  ARowPtr = matrixa21->GetRowPtr();
  AKCol = matrixa21->GetKCol();
  AEntries =  matrixa21->GetEntries();
  BRowPtr = matrixa22->GetRowPtr();
  BKCol = matrixa22->GetKCol();
  BEntries =  matrixa22->GetEntries();

  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    RowPtr[i+1+n2] =  RowPtr[i+n2] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n2] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n2] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n2;
      Entries[k] = BEntries[jb+j];
    }
  }

  /*for(i=0;i<n;i++)
  {
    for (j=RowPtr[i];j< RowPtr[i+1]; j++)
      OutPut("a " << i << " " << ColPtr[j] << " " << Entries[j] << endl);
  }
  return ;*/
  sc.system_type=SCALAR;

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=1;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = RowPtr;
  A0->ja = ColPtr;
  A0->a =  Entries;

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  // sort matrix, diagonal entry
  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i
  /*for(i=0;i<n;i++)
  {
    for (j=0;j< A0->ja[A0->ra[i]] ;j++)
    {
      OutPut("a0 " << A0->ra[i] << " " <<  A0->ja[A0->ra[i]+j] << " "
             << A0->a[A0->ra[i]+j] << endl);
    }
    }*/

  // fill vector structures
  b->n=n;
  b->b=1;
  b->x=rhs;

  x->n=n;
  x->b=1;
  x->x=sol;
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

  /*
    for (i=0;i<n;i++)
      cout << i << " " << sol[i] << " " <<  rhs[i] << endl;
  */

  // copy sol and rhs piece by piece
  for (i=0;i<n2;i++)
    sol1[i] = sol[i];
  for (i=n2;i<n;i++)
    sol2[i-n2] = sol[i];

  delete RowPtr;
  delete ColPtr;
  delete Entries;
  delete sol;
  delete rhs;
  delete A0;
  delete B[0];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* CONNECT 3 BY 3 SYSTEM TO A SCALAR SYSTEM                        */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TSquareMatrix *matrixa12,
TSquareMatrix *matrixa22, TSquareMatrix *matrixa23,
TSquareMatrix *matrixa32, TSquareMatrix *matrixa33,
double *rhs1, double *rhs2, double *rhs3,
double *sol1, double *sol2, double *sol3)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int n,m,i,j, *RowPtr, *ColPtr, n2, n22, *ARowPtr, *AKCol, *BRowPtr, *BKCol;
  double tmp, *Entries, *sol, *rhs, *AEntries, *BEntries;
  int memory[4], ja, ja1, ja2, jb, jb1, jb2, k;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else 
  
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  // connect blocks to single matrix
  // compute dimension (use that all blocks are squares)
  n2 = matrixa11->GetN_Columns();
  n22 = 2*n2;
  n = 3*n2;

  // compute total number of entries
  m = matrixa11->GetN_Entries()+ 2*matrixa12->GetN_Entries()
    + matrixa22->GetN_Entries()+ matrixa23->GetN_Entries()
    + matrixa32->GetN_Entries()+ matrixa33->GetN_Entries();

  // allocate new arrays
  RowPtr = new int[n+1];
  ColPtr = new int[m];
  Entries = new double[m];
  sol = new double[n];
  rhs = new double[n];

  // copy sol and rhs piece by piece
  for (i=0;i<n2;i++)
    sol[i] = sol1[i];
  for (i=n2;i<n22;i++)
    sol[i] = sol2[i-n2];
  for (i=n22;i<n;i++)
    sol[i] = sol3[i-n22];
  for (i=0;i<n2;i++)
    rhs[i] = rhs1[i];
  for (i=n2;i<n22;i++)
    rhs[i] = rhs2[i-n2];
  for (i=n22;i<n;i++)
    rhs[i] = rhs2[i-n22];

  // copy blocks in single matrix (upper part)
  ARowPtr = matrixa11->GetRowPtr();
  AKCol = matrixa11->GetKCol();
  AEntries =  matrixa11->GetEntries();
  BRowPtr = matrixa12->GetRowPtr();
  BKCol = matrixa12->GetKCol();
  BEntries =  matrixa12->GetEntries();

  RowPtr[0] = 0;
  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    // row pointers of a11
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    // row pointers of a12
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    // compute row pointer of new matrix
    RowPtr[i+1] =  RowPtr[i] + ja1 - ja + 2*(jb1 - jb);
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12 and  -matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n2;
      Entries[k] = BEntries[jb+j];
      k = RowPtr[i] + ja2 + jb2+ j;
      ColPtr[k] = BKCol[jb+j]+n22;
      Entries[k] = -BEntries[jb+j];
    }
  }

  // copy blocks in single matrix (middle part)
  ARowPtr = matrixa22->GetRowPtr();
  AKCol = matrixa22->GetKCol();
  AEntries =  matrixa22->GetEntries();
  BRowPtr = matrixa23->GetRowPtr();
  BKCol = matrixa23->GetKCol();
  BEntries =  matrixa23->GetEntries();

  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    RowPtr[i+1+n2] =  RowPtr[i+n2] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n2] + j;
      ColPtr[k] = AKCol[ja+j]+n2;
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n2] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n22;
      Entries[k] = BEntries[jb+j];
    }
  }

  // copy blocks in single matrix (lower part)
  ARowPtr = matrixa32->GetRowPtr();
  AKCol = matrixa32->GetKCol();
  AEntries =  matrixa32->GetEntries();
  BRowPtr = matrixa33->GetRowPtr();
  BKCol = matrixa33->GetKCol();
  BEntries =  matrixa33->GetEntries();

  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    RowPtr[i+1+n22] =  RowPtr[i+n22] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a32
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n22] + j;
      ColPtr[k] = AKCol[ja+j]+n2;
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a33
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n22] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n22;
      Entries[k] = BEntries[jb+j];
    }
  }

  for(i=0;i<n;i++)
  {
    for (j=RowPtr[i];j< RowPtr[i+1]; j++)
      OutPut("a " << i << " " << ColPtr[j] << " " << Entries[j] << endl);
  }
  exit(1) ;
  sc.system_type=SCALAR;

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=1;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = RowPtr;
  A0->ja = ColPtr;
  A0->a =  Entries;

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  // sort matrix, diagonal entry
  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i
  /*for(i=0;i<n;i++)
  {
    for (j=0;j< A0->ja[A0->ra[i]] ;j++)
    {
      OutPut("a0 " << A0->ra[i] << " " <<  A0->ja[A0->ra[i]+j] << " "
             << A0->a[A0->ra[i]+j] << endl);
    }
    }*/

  // fill vector structures
  b->n=n;
  b->b=1;
  b->x=rhs;

  x->n=n;
  x->b=1;
  x->x=sol;
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

  //for (i=0;i<n;i++)
  //  OutPut(i << " " << sol[i] << " " <<  rhs[i] << endl);

  // copy sol and rhs piece by piece
  for (i=0;i<n2;i++)
    sol1[i] = sol[i];
  for (i=n2;i<n22;i++)
    sol2[i-n2] = sol[i];
  for (i=n22;i<n;i++)
    sol3[i-n22] = sol[i];

  delete RowPtr;
  delete ColPtr;
  delete Entries;
  delete sol;
  delete rhs;
  delete A0;
  delete B[0];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* CONNECT SYSTEM FOR VMM [KL02]                                   */
/*******************************************************************/
void Solver(TSquareMatrix *matrixa11, TMatrix *matrixa12, TMatrix *matrixa13,
TMatrix *matrixa21, TMatrix *matrixa31,  TSquareMatrix *matrixa22,
double *rhs1, double *rhs2, double *rhs3,
double *sol1, double *sol2, double *sol3)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[4];
  AMG_VECTOR *x,*b;
  int n,m,i,j, *RowPtr, *ColPtr, n1, n2, *ARowPtr, *AKCol, *BRowPtr, *BKCol;
  int  *CRowPtr, *CKCol;
  double tmp, *Entries, *sol, *rhs, *AEntries, *BEntries, *CEntries;
  int memory[4], ja, ja1, ja2, jb, jb1, jb2, jc, jc1, jc2, k;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  // connect blocks to single matrix
  n1 = matrixa11->GetN_Rows();
  n2 = matrixa21->GetN_Rows();
  n = n1 + 2*n2;

  // compute total number of entries
  m = matrixa11->GetN_Entries()+ matrixa12->GetN_Entries() +  matrixa13->GetN_Entries()
    + matrixa21->GetN_Entries()+  matrixa31->GetN_Entries() + 2*matrixa22->GetN_Entries();

  // allocate new arrays
  RowPtr = new int[n+1];
  ColPtr = new int[m];
  Entries = new double[m];
  sol = new double[n];
  rhs = new double[n];

  // copy sol and rhs piece by piece
  for (i=0;i<n1;i++)
    sol[i] = sol1[i];
  for (i=n1;i<n1+n2;i++)
    sol[i] = sol2[i-n1];
  for (i=n1+n2;i<n;i++)
    sol[i] = sol3[i-n1-n2];
  for (i=0;i<n1;i++)
    rhs[i] = rhs1[i];
  for (i=n1;i<n1+n2;i++)
    rhs[i] = rhs2[i-n1];
  for (i=n1+n2;i<n;i++)
    rhs[i] = rhs3[i-n1-n2];

  // copy blocks in single matrix (upper part)
  ARowPtr = matrixa11->GetRowPtr();
  AKCol = matrixa11->GetKCol();
  AEntries =  matrixa11->GetEntries();
  BRowPtr = matrixa12->GetRowPtr();
  BKCol = matrixa12->GetKCol();
  BEntries =  matrixa12->GetEntries();
  CRowPtr = matrixa13->GetRowPtr();
  CKCol = matrixa13->GetKCol();
  CEntries =  matrixa13->GetEntries();

  RowPtr[0] = 0;
  // for the first n1 rows
  for (i=0;i<n1;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    jc = CRowPtr[i];
    jc1 = CRowPtr[i+1];
    RowPtr[i+1] =  RowPtr[i] + ja1 - ja + jb1 - jb+ jc1 - jc;
    // copy values from matrix_a11
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a12
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n1;
      Entries[k] = BEntries[jb+j];
    }
    // copy values from matrix_a13
    jc2 = jc1 - jc;
    for (j = 0; j < jc2; j++)
    {
      k = RowPtr[i] + ja2 + jb2 + j;
      ColPtr[k] = CKCol[jc+j]+n1+n2;
      Entries[k] = CEntries[jc+j];
    }
  }

  // copy blocks in single matrix (first lower part)
  ARowPtr = matrixa21->GetRowPtr();
  AKCol = matrixa21->GetKCol();
  AEntries =  matrixa21->GetEntries();
  BRowPtr = matrixa22->GetRowPtr();
  BKCol = matrixa22->GetKCol();
  BEntries =  matrixa22->GetEntries();

  // for the second n2 rows
  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    //OutPut(ja << " " << ja1 << " " << jb << " " << jb1 << endl );
    RowPtr[i+1+n1] =  RowPtr[i+n1] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a21
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n1] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a22
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n1] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n1;
      Entries[k] = BEntries[jb+j];
    }
  }

  // copy blocks in single matrix (second lower part)
  ARowPtr = matrixa31->GetRowPtr();
  AKCol = matrixa31->GetKCol();
  AEntries =  matrixa31->GetEntries();
  BRowPtr = matrixa22->GetRowPtr();
  BKCol = matrixa22->GetKCol();
  BEntries =  matrixa22->GetEntries();

  // for the second n2 rows
  for (i=0;i<n2;i++)
  {
    // compute entries in row i
    ja = ARowPtr[i];
    ja1 = ARowPtr[i+1];
    jb = BRowPtr[i];
    jb1 = BRowPtr[i+1];
    //OutPut(ja << " " << ja1 << " " << jb << " " << jb1 << endl );
    RowPtr[i+1+n1+n2] =  RowPtr[i+n1+n2] +  ja1 - ja + jb1 - jb;
    // copy values from matrix_a31
    ja2 = ja1 - ja;
    for (j = 0; j < ja2; j++)
    {
      k = RowPtr[i+n1+n2] + j;
      ColPtr[k] = AKCol[ja+j];
      Entries[k] = AEntries[ja+j];
    }
    // copy values from matrix_a22
    jb2 = jb1 - jb;
    for (j = 0; j < jb2; j++)
    {
      k = RowPtr[i+n1+n2] + ja2 + j;
      ColPtr[k] = BKCol[jb+j]+n1+n2;
      Entries[k] = BEntries[jb+j];
    }
  }

  /*  for(i=0;i<n;i++)
    {
      for (j=RowPtr[i];j< RowPtr[i+1]; j++)
         OutPut("a(" << i+1 << ", " << ColPtr[j]+1 << ") =  " << Entries[j] << ";"<< endl);
      OutPut("sol("<<i+1<< ")=" << sol[i] << ";" << endl);
      OutPut("rhs("<<i+1<< ")=" << rhs[i] << ";" << endl);
    }
    exit(1);
  */
  sc.system_type=SCALAR;

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  // fill matrix structure
  A0->m = n;
  A0->n = n;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=1;
  A0->nonzeros = m;
  A0->bandwidth=-1;
  A0->connections = m;
  A0->ra = RowPtr;
  A0->ja = ColPtr;
  A0->a =  Entries;

  // fill matrix structure
  B[0]->m = n;
  B[0]->n = 0;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros = 0;
  B[0]->bandwidth=-1;
  B[0]->connections = 0;
  B[0]->ra = NULL;
  B[0]->ja = NULL;
  B[0]->a =  NULL;

  B[1] = B[0];

  // sort matrix, diagonal entry
  for (i=0;i<n;i++)
  {                              // first entry in row
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i
  /*for(i=0;i<n;i++)
  {
    for (j=0;j< A0->ja[A0->ra[i]] ;j++)
    {
      OutPut("a0 " << A0->ra[i] << " " <<  A0->ja[A0->ra[i]+j] << " "
             << A0->a[A0->ra[i]+j] << endl);
    }
    }*/

  // fill vector structures
  b->n=n;
  b->b=1;
  b->x=rhs;

  x->n=n;
  x->b=1;
  x->x=sol;
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<n;i++)              //reset array of columns
    A0->ja[A0->ra[i]] = i;

  /*
    for (i=0;i<n;i++)
      cout << i << " " << sol[i] << " " <<  rhs[i] << endl;
  */

  // copy sol and rhs piece by piece
  for (i=0;i<n1;i++)
    sol1[i] = sol[i];
  for (i=n1;i<n;i++)
    sol2[i-n1] = sol[i];

  TDatabase::ParamDB->INTERNAL_AMG_SOLVES++;
  TDatabase::ParamDB->INTERNAL_AMG_PREPARE_TIME+=cc.time;

  //OutPut("Preparation time in AMG was " << cc.time << endl);

  delete RowPtr;
  delete ColPtr;
  delete Entries;
  delete sol;
  delete rhs;
  delete A0;
  delete B[0];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


#ifdef __2D__
/*******************************************************************/
/*        STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 1)              */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
TMatrix *matrixB2, double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4];
  
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_2_TYPE_1;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/*        STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 2)              */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1, TMatrix *matrixB2,
TMatrix *matrixB3, TMatrix *matrixB4, double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[4];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_2_TYPE_2;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB3->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA->GetActiveBound();
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_U;
  B[0]->n = N_P_regular;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_U;
  B[1]->n = N_P_regular;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  B[3]->m = N_P_regular;
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros =  matrixB4->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections = matrixB4->GetN_Entries();
  B[3]->ra =  matrixB4->GetRowPtr();
  B[3]->ja =  matrixB4->GetKCol();
  B[3]->a =   matrixB4->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/*        STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 3)              */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
TMatrix *matrixB1, TMatrix *matrixB2,
double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[5];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_2_TYPE_3;

  N_U = matrixA11->GetN_Rows();
  N_EntriesA = matrixA11->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  B[4] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA11->GetActiveBound();
  A0->ra = matrixA11->GetRowPtr();
  A0->ja = matrixA11->GetKCol();
  A0->a = matrixA11->GetEntries();

  B[2]->m = N_U;                 // fill matrix structure
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros = matrixA12->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections =  matrixA12->GetN_Entries();
  B[2]->active = matrixA12->GetActiveBound();
  B[2]->ra = matrixA12->GetRowPtr();
  B[2]->ja = matrixA12->GetKCol();
  B[2]->a = matrixA12->GetEntries();

  B[3]->m = N_U;                 // fill matrix structure
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros = matrixA21->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections =  matrixA21->GetN_Entries();
  B[3]->active = matrixA21->GetActiveBound();
  B[3]->ra = matrixA21->GetRowPtr();
  B[3]->ja = matrixA21->GetKCol();
  B[3]->a = matrixA21->GetEntries();

  B[4]->m = N_U;                 // fill matrix structure
  B[4]->n = N_U;
  B[4]->b = 1;
  B[4]->bb = 1;
  B[4]->system_as_scalar=1;
  B[4]->blocks_in_diag=1;
  B[4]->nonzeros = matrixA22->GetN_Entries();
  B[4]->bandwidth=-1;
  B[4]->connections =  matrixA22->GetN_Entries();
  B[4]->active = matrixA22->GetActiveBound();
  B[4]->ra = matrixA22->GetRowPtr();
  B[4]->ja = matrixA22->GetKCol();
  B[4]->a = matrixA22->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries of A0
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        // change matrix entries of B[2]
        tmp = B[2]->a[j];
        B[2]->a[j] = B[2]->a[B[2]->ra[i]];
        B[2]->a[B[2]->ra[i]] =tmp;
        // change matrix entries of B[3]
        tmp = B[3]->a[j];
        B[3]->a[j] = B[3]->a[B[3]->ra[i]];
        B[3]->a[B[3]->ra[i]] =tmp;
        // change matrix entries of B[4]
        tmp = B[4]->a[j];
        B[4]->a[j] = B[4]->a[B[4]->ra[i]];
        B[4]->a[B[4]->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete B[4];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else 
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/*        STOKES TYPE SADDLE POINT PROBLEM (NSTYPE 4)              */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
TMatrix *matrixB1, TMatrix *matrixB2,
TMatrix *matrixB3, TMatrix *matrixB4,
double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[7];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_2_TYPE_4;

  N_U = matrixA11->GetN_Rows();
  N_EntriesA = matrixA11->GetN_Entries();
  N_P = matrixB3->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  B[4] = new AMG_MATRIX;         // allocate matrix structure
  B[5] = new AMG_MATRIX;         // allocate matrix structure
  B[6] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA11->GetActiveBound();
  A0->ra = matrixA11->GetRowPtr();
  A0->ja = matrixA11->GetKCol();
  A0->a = matrixA11->GetEntries();

  B[4]->m = N_U;                 // fill matrix structure
  B[4]->n = N_U;
  B[4]->b = 1;
  B[4]->bb = 1;
  B[4]->system_as_scalar=1;
  B[4]->blocks_in_diag=1;
  B[4]->nonzeros = matrixA12->GetN_Entries();
  B[4]->bandwidth=-1;
  B[4]->connections =  matrixA12->GetN_Entries();
  B[4]->active = matrixA12->GetActiveBound();
  B[4]->ra = matrixA12->GetRowPtr();
  B[4]->ja = matrixA12->GetKCol();
  B[4]->a = matrixA12->GetEntries();

  B[5]->m = N_U;                 // fill matrix structure
  B[5]->n = N_U;
  B[5]->b = 1;
  B[5]->bb = 1;
  B[5]->system_as_scalar=1;
  B[5]->blocks_in_diag=1;
  B[5]->nonzeros = matrixA21->GetN_Entries();
  B[5]->bandwidth=-1;
  B[5]->connections =  matrixA21->GetN_Entries();
  B[5]->active = matrixA21->GetActiveBound();
  B[5]->ra = matrixA21->GetRowPtr();
  B[5]->ja = matrixA21->GetKCol();
  B[5]->a = matrixA21->GetEntries();

  B[6]->m = N_U;                 // fill matrix structure
  B[6]->n = N_U;
  B[6]->b = 1;
  B[6]->bb = 1;
  B[6]->system_as_scalar=1;
  B[6]->blocks_in_diag=1;
  B[6]->nonzeros = matrixA22->GetN_Entries();
  B[6]->bandwidth=-1;
  B[6]->connections =  matrixA22->GetN_Entries();
  B[6]->active = matrixA22->GetActiveBound();
  B[6]->ra = matrixA22->GetRowPtr();
  B[6]->ja = matrixA22->GetKCol();
  B[6]->a = matrixA22->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_U;
  B[0]->n = N_P_regular;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_U;
  B[1]->n = N_P_regular;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  B[3]->m = N_P_regular;
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros =  matrixB4->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections = matrixB4->GetN_Entries();
  B[3]->ra =  matrixB4->GetRowPtr();
  B[3]->ja =  matrixB4->GetKCol();
  B[3]->a =   matrixB4->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries of A0
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        // change matrix entries of B[4]
        tmp = B[4]->a[j];
        B[4]->a[j] = B[4]->a[B[4]->ra[i]];
        B[4]->a[B[4]->ra[i]] =tmp;
        // change matrix entries of B[5]
        tmp = B[5]->a[j];
        B[5]->a[j] = B[5]->a[B[5]->ra[i]];
        B[5]->a[B[5]->ra[i]] =tmp;
        // change matrix entries of B[6]
        tmp = B[6]->a[j];
        B[6]->a[j] = B[6]->a[B[6]->ra[i]];
        B[6]->a[B[6]->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete B[4];
  delete B[5];
  delete B[6];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}
#endif

int CheckBraessSarazinParameters()
{
  // set appropriate parameters for preconditioner
  switch(TDatabase::ParamDB->SC_BRAESS_SARAZIN_MATRIX)
  {
    case 0:                      // original matrix
      break;
    case 1:                      // identity
    case 2:                      // diagonal matrix
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = AMG_SCHUR_GMRES;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A = AMG_DJAC;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A_MAXIT = 1;
      break;
    case 3:                      // ILU
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = AMG_SCHUR_GMRES;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A = AMG_ILU;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A_MAXIT = 1;
      break;
    case 4:                      // ILUT
      TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR = AMG_SCHUR_GMRES;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A = AMG_ILUT;
      TDatabase::ParamDB->SC_SCHUR_INV_OF_A_MAXIT = 1;
      break;
  }
  return (0);
}


/*******************************************************************/
/* BRAESS--SARAZIN SMOOTHER (NSTYPE 1)                             */
/*******************************************************************/
#ifdef __2D__
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
TMatrix *matrixB2, double *rhs, double *sol,
int para0)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = BRAESS_SARAZIN_SADDLE_2_TYPE_1;

  alpha_change = 0;
  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  if (TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA<=0)
  {
    alpha_change = 1;
    tmp = 0.0;
    for (i=0;i<N_U;i++)
    {
      if(A0->a[A0->ra[i]] > tmp)
        tmp = A0->a[A0->ra[i]];
    }
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = 2*tmp;
  }

  CheckBraessSarazinParameters();
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  if (alpha_change)
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = -1;

  delete A0;
  delete B[0];
  delete B[1];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* BRAESS--SARAZIN SMOOTHER (NSTYPE 2)                             */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1, TMatrix *matrixB2,
TMatrix *matrixB3, TMatrix *matrixB4, double *rhs, double *sol,
int para0)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[4];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = BRAESS_SARAZIN_SADDLE_2_TYPE_2;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB3->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA->GetActiveBound();
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_U;
  B[0]->n = N_P_regular;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_U;
  B[1]->n = N_P_regular;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  B[3]->m = N_P_regular;
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros =  matrixB4->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections = matrixB4->GetN_Entries();
  B[3]->ra =  matrixB4->GetRowPtr();
  B[3]->ja =  matrixB4->GetKCol();
  B[3]->a =   matrixB4->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  if (TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA<=0)
  {
    alpha_change = 1;
    tmp = 0.0;
    for (i=0;i<N_U;i++)
    {
      if(A0->a[A0->ra[i]] > tmp)
        tmp = A0->a[A0->ra[i]];
    }
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = 2*tmp;
  }

  CheckBraessSarazinParameters();
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* BRAESS--SARAZIN SMOOTHER (NSTYPE 3)                             */
/*******************************************************************/
void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
TMatrix *matrixB1,
TMatrix *matrixB2, double *rhs, double *sol,
int para0)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[5];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = BRAESS_SARAZIN_SADDLE_2_TYPE_3;

  alpha_change = 0;
  N_U = matrixA11->GetN_Rows();
  N_EntriesA = matrixA11->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  B[4] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->ra = matrixA11->GetRowPtr();
  A0->ja = matrixA11->GetKCol();
  A0->a = matrixA11->GetEntries();

  B[2]->m = N_U;                 // fill matrix structure
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros = matrixA12->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections =  matrixA12->GetN_Entries();
  B[2]->active = matrixA12->GetActiveBound();
  B[2]->ra = matrixA12->GetRowPtr();
  B[2]->ja = matrixA12->GetKCol();
  B[2]->a = matrixA12->GetEntries();

  B[3]->m = N_U;                 // fill matrix structure
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros = matrixA21->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections =  matrixA21->GetN_Entries();
  B[3]->active = matrixA21->GetActiveBound();
  B[3]->ra = matrixA21->GetRowPtr();
  B[3]->ja = matrixA21->GetKCol();
  B[3]->a = matrixA21->GetEntries();

  B[4]->m = N_U;                 // fill matrix structure
  B[4]->n = N_U;
  B[4]->b = 1;
  B[4]->bb = 1;
  B[4]->system_as_scalar=1;
  B[4]->blocks_in_diag=1;
  B[4]->nonzeros = matrixA22->GetN_Entries();
  B[4]->bandwidth=-1;
  B[4]->connections =  matrixA22->GetN_Entries();
  B[4]->active = matrixA22->GetActiveBound();
  B[4]->ra = matrixA22->GetRowPtr();
  B[4]->ja = matrixA22->GetKCol();
  B[4]->a = matrixA22->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries of A0
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        // change matrix entries of B[2]
        tmp = B[2]->a[j];
        B[2]->a[j] = B[2]->a[B[2]->ra[i]];
        B[2]->a[B[2]->ra[i]] =tmp;
        // change matrix entries of B[3]
        tmp = B[3]->a[j];
        B[3]->a[j] = B[3]->a[B[3]->ra[i]];
        B[3]->a[B[3]->ra[i]] =tmp;
        // change matrix entries of B[4]
        tmp = B[4]->a[j];
        B[4]->a[j] = B[4]->a[B[4]->ra[i]];
        B[4]->a[B[4]->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  if (TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA<=0)
  {
    alpha_change = 1;
    tmp = 0.0;
    for (i=0;i<N_U;i++)
    {
      if(A0->a[A0->ra[i]] > tmp)
        tmp = A0->a[A0->ra[i]];
    }
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = 2*tmp;
  }

  CheckBraessSarazinParameters();
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  if (alpha_change)
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = -1;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete B[4];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* BRAESS--SARAZIN SMOOTHER (NSTYPE 4)                             */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA11,TSquareMatrix *matrixA12,
TSquareMatrix *matrixA21,TSquareMatrix *matrixA22,
TMatrix *matrixB1, TMatrix *matrixB2,
TMatrix *matrixB3, TMatrix *matrixB4,
double *rhs, double *sol, int para0)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[7];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = BRAESS_SARAZIN_SADDLE_2_TYPE_4;

  alpha_change = 0;
  N_U = matrixA11->GetN_Rows();
  N_EntriesA = matrixA11->GetN_Entries();
  N_P = matrixB3->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  B[4] = new AMG_MATRIX;         // allocate matrix structure
  B[5] = new AMG_MATRIX;         // allocate matrix structure
  B[6] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA11->GetActiveBound();
  A0->ra = matrixA11->GetRowPtr();
  A0->ja = matrixA11->GetKCol();
  A0->a = matrixA11->GetEntries();

  B[4]->m = N_U;                 // fill matrix structure
  B[4]->n = N_U;
  B[4]->b = 1;
  B[4]->bb = 1;
  B[4]->system_as_scalar=1;
  B[4]->blocks_in_diag=1;
  B[4]->nonzeros = matrixA12->GetN_Entries();
  B[4]->bandwidth=-1;
  B[4]->connections =  matrixA12->GetN_Entries();
  B[4]->active = matrixA12->GetActiveBound();
  B[4]->ra = matrixA12->GetRowPtr();
  B[4]->ja = matrixA12->GetKCol();
  B[4]->a = matrixA12->GetEntries();

  B[5]->m = N_U;                 // fill matrix structure
  B[5]->n = N_U;
  B[5]->b = 1;
  B[5]->bb = 1;
  B[5]->system_as_scalar=1;
  B[5]->blocks_in_diag=1;
  B[5]->nonzeros = matrixA21->GetN_Entries();
  B[5]->bandwidth=-1;
  B[5]->connections =  matrixA21->GetN_Entries();
  B[5]->active = matrixA21->GetActiveBound();
  B[5]->ra = matrixA21->GetRowPtr();
  B[5]->ja = matrixA21->GetKCol();
  B[5]->a = matrixA21->GetEntries();

  B[6]->m = N_U;                 // fill matrix structure
  B[6]->n = N_U;
  B[6]->b = 1;
  B[6]->bb = 1;
  B[6]->system_as_scalar=1;
  B[6]->blocks_in_diag=1;
  B[6]->nonzeros = matrixA22->GetN_Entries();
  B[6]->bandwidth=-1;
  B[6]->connections =  matrixA22->GetN_Entries();
  B[6]->active = matrixA22->GetActiveBound();
  B[6]->ra = matrixA22->GetRowPtr();
  B[6]->ja = matrixA22->GetKCol();
  B[6]->a = matrixA22->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_U;
  B[0]->n = N_P_regular;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_U;
  B[1]->n = N_P_regular;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  B[3]->m = N_P_regular;
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros =  matrixB4->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections = matrixB4->GetN_Entries();
  B[3]->ra =  matrixB4->GetRowPtr();
  B[3]->ja =  matrixB4->GetKCol();
  B[3]->a =   matrixB4->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries of A0
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        // change matrix entries of B[4]
        tmp = B[4]->a[j];
        B[4]->a[j] = B[4]->a[B[4]->ra[i]];
        B[4]->a[B[4]->ra[i]] =tmp;
        // change matrix entries of B[5]
        tmp = B[5]->a[j];
        B[5]->a[j] = B[5]->a[B[5]->ra[i]];
        B[5]->a[B[5]->ra[i]] =tmp;
        // change matrix entries of B[6]
        tmp = B[6]->a[j];
        B[6]->a[j] = B[6]->a[B[6]->ra[i]];
        B[6]->a[B[6]->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  if (TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA<=0)
  {
    alpha_change = 1;
    tmp = 0.0;
    for (i=0;i<N_U;i++)
    {
      if(A0->a[A0->ra[i]] > tmp)
        tmp = A0->a[A0->ra[i]];
    }
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = 2*tmp;
  }

  CheckBraessSarazinParameters();
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  if (alpha_change)
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = -1;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete B[4];
  delete B[5];
  delete B[6];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* STOKES TYPE SADDLE POINT PROBLEM WITH MORTAR (NSTYPE 2)         */
/*******************************************************************/
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1, TMatrix *matrixB2,
TMatrix *matrixB3, TMatrix *matrixB4, TMatrix *matrixmortar,
double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[5];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular, N_Mortar;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_2_TYPE_2_MORTAR;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB3->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();
  N_Mortar = matrixmortar->GetN_Rows();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  B[3] = new AMG_MATRIX;         // allocate matrix structure
  B[4] = new AMG_MATRIX;         // allocate matrix structure

  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA->GetActiveBound();
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_U;
  B[0]->n = N_P_regular;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_U;
  B[1]->n = N_P_regular;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  B[3]->m = N_P_regular;
  B[3]->n = N_U;
  B[3]->b = 1;
  B[3]->bb = 1;
  B[3]->system_as_scalar=1;
  B[3]->blocks_in_diag=1;
  B[3]->nonzeros =  matrixB4->GetN_Entries();
  B[3]->bandwidth=-1;
  B[3]->connections = matrixB4->GetN_Entries();
  B[3]->ra =  matrixB4->GetRowPtr();
  B[3]->ja =  matrixB4->GetKCol();
  B[3]->a =   matrixB4->GetEntries();

                                 // fill matrix structure
  B[4]->m = matrixmortar->GetN_Rows();
  B[4]->n = matrixmortar->GetN_Columns();
  B[4]->b = 1;
  B[4]->bb = 1;
  B[4]->system_as_scalar = 1;
  B[4]->blocks_in_diag = 1;
  B[4]->nonzeros =  matrixmortar->GetN_Entries();
  B[4]->bandwidth = -1;
  B[4]->connections = matrixmortar->GetN_Entries();
  B[4]->ra = matrixmortar->GetRowPtr();
  B[4]->ja = matrixmortar->GetKCol();
  B[4]->a = matrixmortar->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U + N_P_regular + 2* N_Mortar;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U + N_P_regular+ 2* N_Mortar;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete B[3];
  delete B[4];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/*  AUXILIARY PROBLEM FOR VASSILEVSKI/LAZAROV APPROACH             */
/*******************************************************************/

void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1, TMatrix *matrixB2,
double *rhs, double *sol,double delta)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[2];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB;
  int memory[4];
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SCALAR_VAL_LAZ_2;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar=1;
  A0->blocks_in_diag=2;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth=-1;
  A0->connections = N_EntriesA;
  A0->active = matrixA->GetActiveBound();
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  B[0]->m = N_P;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 2*N_U;
  b->b = 1;
  b->x = rhs;

  x->n = 2*N_U;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;

}
#endif

#ifdef __3D__
/*******************************************************************/
/* BRAESS--SARAZIN SMOOTHER (NSTYPE 1)                             */
/*******************************************************************/
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
TMatrix *matrixB2, TMatrix *matrixB3,
double *rhs, double *sol,
int para0)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[3];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = BRAESS_SARAZIN_SADDLE_3_TYPE_1;

  alpha_change = 0;
  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar = 1;
  A0->blocks_in_diag = 3;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth = -1;
  A0->connections = N_EntriesA;
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 3*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 3*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  if (TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA<=0)
  {
    alpha_change = 1;
    tmp = 0.0;
    for (i=0;i<N_U;i++)
    {
      if(A0->a[A0->ra[i]] > tmp)
        tmp = A0->a[A0->ra[i]];
    }
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = 2*tmp;
  }

  CheckBraessSarazinParameters();
  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  if (alpha_change)
    TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA = -1;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}


/*******************************************************************/
/* STOKES-TYPE SYSTEM (NSTYPE 1)                                   */
/*******************************************************************/
void Solver(TSquareMatrix *matrixA, TMatrix *matrixB1,
TMatrix *matrixB2, TMatrix *matrixB3,
double *rhs, double *sol)
{
  AMG_CoarsenContext cc;
  AMG_SolverContext sc;
  AMG_MATRIX *A0,*B[3];
  AMG_VECTOR *x,*b;
  int i,j;
  double tmp;
  int N_U, N_P, N_EntriesA, N_EntriesB, N_P_regular;
  int memory[4], alpha_change;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;  
 info = mstats();
 memory[0]=memory[1]=memory[2]=memory[3]=0.;
#else   
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  sc.system_type = SADDLE_3_TYPE_1;

  N_U = matrixA->GetN_Rows();
  N_EntriesA = matrixA->GetN_Entries();
  N_P = matrixB1->GetN_Rows();
  N_EntriesB = matrixB1->GetN_Entries();

  A0 = new AMG_MATRIX;           // allocate matrix structure
  B[0] = new AMG_MATRIX;         // allocate matrix structure
  B[1] = new AMG_MATRIX;         // allocate matrix structure
  B[2] = new AMG_MATRIX;         // allocate matrix structure
  b = new AMG_VECTOR;            // allocate vector structures
  x = new AMG_VECTOR;

  A0->m = N_U;                   // fill matrix structure
  A0->n = N_U;
  A0->b = 1;
  A0->bb = 1;
  A0->system_as_scalar = 1;
  A0->blocks_in_diag = 3;
  A0->nonzeros = N_EntriesA;
  A0->bandwidth = -1;
  A0->connections = N_EntriesA;
  A0->ra = matrixA->GetRowPtr();
  A0->ja = matrixA->GetKCol();
  A0->a = matrixA->GetEntries();

  N_P_regular = N_P;
  B[0]->m = N_P_regular ;
  B[0]->n = N_U;
  B[0]->b = 1;
  B[0]->bb = 1;
  B[0]->system_as_scalar=1;
  B[0]->blocks_in_diag=1;
  B[0]->nonzeros =  matrixB1->GetN_Entries();
  B[0]->bandwidth=-1;
  B[0]->connections = matrixB1->GetN_Entries();
  B[0]->ra =  matrixB1->GetRowPtr();
  B[0]->ja =  matrixB1->GetKCol();
  B[0]->a =   matrixB1->GetEntries();

  B[1]->m = N_P_regular;
  B[1]->n = N_U;
  B[1]->b = 1;
  B[1]->bb = 1;
  B[1]->system_as_scalar=1;
  B[1]->blocks_in_diag=1;
  B[1]->nonzeros =  matrixB2->GetN_Entries();
  B[1]->bandwidth=-1;
  B[1]->connections = matrixB2->GetN_Entries();
  B[1]->ra =  matrixB2->GetRowPtr();
  B[1]->ja =  matrixB2->GetKCol();
  B[1]->a =   matrixB2->GetEntries();

  B[2]->m = N_P_regular;
  B[2]->n = N_U;
  B[2]->b = 1;
  B[2]->bb = 1;
  B[2]->system_as_scalar=1;
  B[2]->blocks_in_diag=1;
  B[2]->nonzeros =  matrixB3->GetN_Entries();
  B[2]->bandwidth=-1;
  B[2]->connections = matrixB3->GetN_Entries();
  B[2]->ra =  matrixB3->GetRowPtr();
  B[2]->ja =  matrixB3->GetKCol();
  B[2]->a =   matrixB3->GetEntries();

  // sort matrix, diagonal entry first entry in row
  for (i=0;i<N_U;i++)
  {
    for (j= A0->ra[i];j<A0->ra[i+1];j++)
    {
      if (A0->ja[j]==i)
      {
        // set non diagonal entry of ja
        A0->ja[j] =A0->ja[A0->ra[i]];
        // diagonal entry gets number of row entries
        A0->ja[A0->ra[i]] =A0->ra[i+1]-A0->ra[i];
        // change matrix entries
        tmp = A0->a[j];
        A0->a[j] = A0->a[A0->ra[i]];
        A0->a[A0->ra[i]] =tmp;
        break;
      }                          // endif
    }                            // endfor j
  }                              // endfor i

  b->n = 3*N_U + N_P_regular;
  b->b = 1;
  b->x = rhs;

  x->n = 3*N_U + N_P_regular;
  x->b = 1;
  x->x = sol;

  SetAMGDefaults(cc, sc);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  AMG(&sc,&cc,A0,B,x,b);
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  for (i=0;i<N_U;i++)            //reset array of columns
    A0->ja[A0->ra[i]] = i;

  delete A0;
  delete B[0];
  delete B[1];
  delete B[2];
  delete b;
  delete x;
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
#else  
  MALLINFO = mallinfo();
  memory[3] = MALLINFO.usmblks+MALLINFO.uordblks;
#endif
  if (memory[3] - memory[0])
    cout << "WARNING : Solver did not set all memory free !!!" <<
      memory[3] - memory[0] << endl;
  if (memory[2] - memory[1])
    cout << "WARNING : AMG did not set all memory free !!!" <<
      memory[2] - memory[1] << endl;
}
#endif
