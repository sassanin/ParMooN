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
   
/****************************************************************************/
/*                                                                          */
/* File:          amg_solve_main.c                                          */
/*                                                                          */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                          */
/* History:   1998/02/19 start using this library for MooN_MD               */
/*                                                                         */
/* Remarks:   1998/02/24 ILU                                                */
/*              1998/02/27 step length control                              */
/*              1998/03/02 GMRES                                            */
/*              1998/03/12 flexible GMRES                                   */
/*              1998/04/01 exact solver                                     */
/*              1998/06/03 ILUT                                             */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*    system include files                                                    */
/*    application include files                                             */
/*                                                                            */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// #ifdef __MAC64__
// #include <malloc/malloc.h>
// #else
// #include <malloc.h>
// #endif

#include <amg_header.h>
#include <amg_low.h>
#include <amg_sp.h>
#include <amg_blas.h>
#include <amg_iter.h>
#include <amg_coarsen.h>
#include <amg_solve_main.h>
#include <amg_solvers.h>
#include <amg_1d_prec.h>
#include <amg_2d_prec.h>

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Iteration procedure, this may be multigrid, sor, jacobi, ilu, ...        */
/* The aim is to solve Mx=b with the following constraints:                 */
/* M is an auxiliary matrix containing e.g. an ILU decomposition of A       */
/* x contains an iterate on entry and holds the new iterate on exit         */
/* b is the right hand side THAT MAY NOT BE CHANGED BY THE PROCEDURE !      */
/* d contains the defect d=b-Ax on entry, arbitrary on exit                 */
/*                                                                          */
/****************************************************************************/

extern IterProcPtr smoother,coarse_smoother;
extern IterProcPtr preconditioner,schur_preconditioner,preconditioner_trans;

extern MultProcPtr dmatmul,dmatminus,A_dmatmul,A_dmatminus,B_dmatmul,B_Trans_dmatmul,
  MGC_dmatminus, dmattransmul;

extern RestrictProcPtr restriction;
extern InterpolationProcPtr interpolation;

extern double start_residual,end_residual,residuals[AMG_CONV_RATE_BACK];
extern int residual_cnt,iteration_cnt;
extern double elapsed_time,TIME_WRAP;
extern clock_t start,finish,start1,finish1;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

extern int coarse_grid;          /* multigrid is on coarse grid */
extern char buf[128];

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/* global data for solvers */
extern AMG_MATRIX *A[AMG_MAX_LEVELS];
extern AMG_MATRIX *B[AMG_MAX_LEVELS];
extern AMG_MATRIX *B_Diri[AMG_MAX_LEVELS];
extern AMG_MATRIX *schur_matrix[AMG_MAX_LEVELS];
extern AMG_MATRIX *M[AMG_MAX_LEVELS];
extern AMG_GRAPH  *G[AMG_MAX_LEVELS];
extern AMG_GRAPH  *G_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *x[AMG_MAX_LEVELS];
extern AMG_VECTOR *b[AMG_MAX_LEVELS];
extern AMG_VECTOR *d[AMG_MAX_LEVELS];
extern AMG_VECTOR *z[AMG_MAX_LEVELS];
extern AMG_VECTOR *r[AMG_MAX_LEVELS];
extern AMG_VECTOR *q;
extern AMG_VECTOR *p[AMG_MAX_LEVELS];
extern AMG_VECTOR *w;
extern AMG_VECTOR *old_def[AMG_MAX_LEVELS];
extern AMG_VECTOR *A_times_update[AMG_MAX_LEVELS];
extern AMG_VECTOR *s,*cosi,*sn;
extern AMG_VECTOR *H[AMG_MAX_GMRES_RESTART+1];
extern AMG_VECTOR *v[AMG_MAX_GMRES_RESTART+1];
extern AMG_VECTOR *zv[AMG_MAX_GMRES_RESTART+1];
extern AMG_VECTOR *old_def[AMG_MAX_LEVELS];        /* array of vectors for step length control*/
extern AMG_VECTOR *A_times_update[AMG_MAX_LEVELS]; /* array of vectors for step length control*/
extern AMG_VECTOR *row_equilibration;
extern int depth;
extern int mgc_recursion[AMG_MAX_LEVELS];
extern AMG_CoarsenContext *global_cc;
extern AMG_SolverContext *global_sc;

extern AMG_VECTOR *s_schur,*cosi_schur,*sn_schur;
extern AMG_VECTOR *H_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
extern AMG_VECTOR *v_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
extern AMG_VECTOR *schur_velo[3][AMG_MAX_LEVELS], *schur_press[5];

extern AMG_VECTOR *u,*v1,*w,*w1,*t[AMG_MAX_LEVELS],*tilde_r;

extern AMG_VECTOR *velo_prolong[AMG_MAX_LEVELS],*pres_prolong[AMG_MAX_LEVELS];
extern AMG_VECTOR *velo_result,*pres_result;

extern AMG_VECTOR *velo_rhs,*velo_help,*pres_help[3];
extern AMG_VECTOR *d_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *z_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *r_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *p_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *w_bcgs_schur;

/****************************************************************************/
/*                                                                            */
/* PREPARE THE SOLVERS                                                      */
/*                                                                            */
/****************************************************************************/
int prepare_scalar_system(AMG_SolverContext *sc, AMG_CoarsenContext *cc, 
                          AMG_MATRIX *A_in, int *A_length,
                          int *B_length,int *prepare_amg)
{
  int k,blocks;

  switch(sc->system_type)
    {
    case SCALAR1 :  
      dmatmul = AMG_dmatmul;                     /* set matrix vector routines */
      dmatminus =  AMG_dmatminus;
      A_dmatmul =  AMG_dmatmul; 
      A_dmatminus =  AMG_dmatminus;
      B_dmatmul = B_Trans_dmatmul = NULL;
      MGC_dmatminus =  dmatminus;
      dmattransmul = AMG_dmattransmul;                     
      A_length[0] = AMG_MATRIX_N(A_in);          /* length of A-components solution vector */
      for (k=0;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                         /* length of B-components solution vector */
      if  ((sc->preconditioner==AMG_SCHUR_COMPLEMENT)
           ||(sc->preconditioner==AMG_SCHUR_CG)
           ||(sc->preconditioner==AMG_SCHUR_GMRES)
           ||(sc->preconditioner==AMG_SCHUR_GMRES_BCGS)
           ||(sc->preconditioner==AMG_BRAESAR))
        {
          AMG_Print("MESSAGE : Saddle point type preconditioner does not work for\n");
          AMG_Print("MESSAGE : scalar systems.\n");                             
          AMG_Print("MESSAGE : Preconditioner changed to AMG_MGC (algebraic multigrid)\n");     
          sc->preconditioner = AMG_MGC;
        }
      if (sc->row_equilibration)
        ApplyRowEquilibration(sc,A_in);
      blocks =1 ;
      break;
    case SCALAR2 :
      dmatmul = AMG_dmatmul_SCALAR2;                     /* set matrix vector routines */
      dmatminus =  AMG_dmatminus_SCALAR2;
      A_dmatmul =  AMG_dmatmul_SCALAR2; 
      A_dmatminus =  AMG_dmatminus_SCALAR2;
      B_dmatmul = B_Trans_dmatmul = NULL;
      MGC_dmatminus =  dmatminus;
      A_length[0] = 2*AMG_MATRIX_N(A_in);          /* length of A-components solution vector */
      for (k=0;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                         /* length of B-components solution vector */
      if  ((sc->preconditioner==AMG_SCHUR_COMPLEMENT)
           ||(sc->preconditioner==AMG_SCHUR_CG)
           ||(sc->preconditioner==AMG_SCHUR_GMRES)
           ||(sc->preconditioner==AMG_SCHUR_GMRES_BCGS)
           ||(sc->preconditioner==AMG_BRAESAR))
        /* ||(sc->preconditioner==AMG_MGC)*/
        {
          AMG_Print("MESSAGE : Preconditioner does not work for SCALAR2\n");
          AMG_Print("MESSAGE : Preconditioner changed to AMG_SSOR\n");     
          sc->preconditioner = AMG_SSOR;
        }
      blocks = 2;
      sc->row_equilibration=0;
     break;
   case SCALAR3 :
      dmatmul = AMG_dmatmul_SCALAR3;                     /* set matrix vector routines */
      dmatminus =  AMG_dmatminus_SCALAR3;
      A_dmatmul =  AMG_dmatmul_SCALAR3; 
      A_dmatminus =  AMG_dmatminus_SCALAR3;
      B_dmatmul = B_Trans_dmatmul = NULL;
      MGC_dmatminus =  dmatminus;
      A_length[0] = 3*AMG_MATRIX_N(A_in);          /* length of A-components solution vector */
      for (k=0;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                         /* length of B-components solution vector */
      if  ((sc->preconditioner==AMG_SCHUR_COMPLEMENT)
           ||(sc->preconditioner==AMG_SCHUR_CG)
           ||(sc->preconditioner==AMG_SCHUR_GMRES)
           ||(sc->preconditioner==AMG_SCHUR_GMRES_BCGS)
           ||(sc->preconditioner==AMG_BRAESAR))
        /* ||(sc->preconditioner==AMG_MGC)*/
        {
          AMG_Print("MESSAGE : Preconditioner does not work for SCALAR3\n");
          AMG_Print("MESSAGE : Preconditioner changed to AMG_SSOR\n");     
          sc->preconditioner = AMG_SSOR;
        }
      blocks = 3;
      sc->row_equilibration=0;
     break;
    case SCALAR6 :
      dmatmul = AMG_dmatmul_SCALAR6;                     /* set matrix vector routines */
      dmatminus =  AMG_dmatminus_SCALAR6;
      A_dmatmul =  AMG_dmatmul_SCALAR6; 
      A_dmatminus =  AMG_dmatminus_SCALAR6;
      B_dmatmul = B_Trans_dmatmul = NULL;
      MGC_dmatminus =  dmatminus;
      A_length[0] = 6*AMG_MATRIX_N(A_in);          /* length of A-components solution vector */
      for (k=0;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                         /* length of B-components solution vector */
      if  ((sc->preconditioner==AMG_SCHUR_COMPLEMENT)
           ||(sc->preconditioner==AMG_SCHUR_CG)
           ||(sc->preconditioner==AMG_SCHUR_GMRES)
           ||(sc->preconditioner==AMG_SCHUR_GMRES_BCGS)
           ||(sc->preconditioner==AMG_BRAESAR))
        /* ||(sc->preconditioner==AMG_MGC)*/
        {
          AMG_Print("MESSAGE : Preconditioner does not work for SCALAR6\n");
          AMG_Print("MESSAGE : Preconditioner changed to AMG_SSOR\n");     
          sc->preconditioner = AMG_SSOR;
        }
      blocks = 6;
      sc->row_equilibration=0;
     break;
    case SCALAR_VAL_LAZ_2 : 
      dmatmul = AMG_dmatmul_SCALAR_VAL_LAZ_2; 
      dmatminus =  AMG_dmatminus_SCALAR_VAL_LAZ_2;
      A_dmatmul = AMG_dmatmul_SCALAR_VAL_LAZ_2;
      A_dmatminus =  AMG_dmatminus_SCALAR_VAL_LAZ_2;
      B_dmatmul = NULL;
      B_Trans_dmatmul = NULL;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = 0;                             /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      sc->amg_prec_it=0;                          /* no preconditioning possible */ 
      sc->row_equilibration=0;
      break;
    }
  if ((sc->amg_prec_it==0) || 
      (sc->solver==AMG_EXACT) || 
      (sc->preconditioner == AMG_NO_PRECONDITIONER))
    {
      depth=0;
      A[0]=A_in;
      sc->amg_prec_it=0;
      sc->preconditioner = AMG_NO_PRECONDITIONER; 
      return(0);
    }
  switch (sc->preconditioner)
    {
    case AMG_MGC :                            /* prepare amg */
      *prepare_amg=1;
      depth=AMG_BuildHierarchy(cc,A_in,A,G);  /* build amg */
      if (depth<0) 
        {
          AMG_Print("Could not set up coarse grid matrices\n");
          exit(4711);
        }
      for (k=1;k<=depth;k++)                  /* no of dof on coarse grids */
        {
          A_length[k] = blocks*AMG_MATRIX_N(A[k]);
          A[k]->level=k;
          A[k]->m=A[k]->n;
        }
      if ((sc->solver==AMG_GMRES_LEFT)||(sc->solver==AMG_GMRES_RIGHT))  
        {                                     /* no appropriate combinations */
          AMG_Print("MESSAGE : Static preconditioned GMRES with multigrid will\n");
          AMG_Print("MESSAGE : in general not converge !!!\n");                             
          AMG_Print("MESSAGE : Solver changed to flexible GMRES (AMG_GMRES_FLEX)\n");     
          sc->solver = AMG_GMRES_FLEX;
        }
      if (sc->system_type == SCALAR1)
      {
        restriction = pc_restrict;              /* set grid transfer operations */
        interpolation = pc_prolongate_auto;
      }
      if (sc->system_type == SCALAR2)
      {
        restriction = pc_restrict_2d;              /* set grid transfer operations */
        interpolation = pc_prolongate_auto_2d;
      }
      if (sc->system_type == SCALAR3)
      {
        restriction = pc_restrict_3d;              /* set grid transfer operations */
        interpolation = pc_prolongate_auto_3d;
      }
      if (sc->system_type == SCALAR6)
      {
        restriction = pc_restrict_6d;              /* set grid transfer operations */
        interpolation = pc_prolongate_auto_6d;
      }
      break;
    default:
      depth=0;
      A[0]=A_in;
    }  
  return(0);
}
int prepare_saddle_system(AMG_SolverContext *sc, AMG_CoarsenContext *cc, 
                          AMG_MATRIX *A_in, int *A_length,
                          int *B_length,int *prepare_amg)
{
  int k,blocks;
  sc->row_equilibration = 0;
  switch(sc->system_type)
    {
    case SADDLE_1 :
      dmatmul = AMG_dmatmul_SADDLE_1 ; 
      dmatminus =  AMG_dmatminus_SADDLE_1;
      A_dmatmul = AMG_A_dmatmul_SADDLE_1;
      A_dmatminus =  AMG_dmatminus;
      B_dmatmul = AMG_B_dmatmul_SADDLE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_1;
      A_length[0] = AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_N(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks=1;
      break;
    case SADDLE_1_TYPE_1 : 
      dmatmul = AMG_dmatmul_SADDLE_1_TYPE_1 ; 
      dmatminus =  AMG_dmatminus_SADDLE_1_TYPE_1;
      A_dmatmul = AMG_A_dmatmul_SADDLE_1;
      A_dmatminus =  AMG_dmatminus;
      B_dmatmul = AMG_B_dmatmul_SADDLE_1_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_1_TYPE_1;
      A_length[0] = AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 1;
     break;
    case SADDLE_2_TYPE_1 : 
      dmatmul = AMG_dmatmul_SADDLE_2_TYPE_1 ; 
      dmatminus =  AMG_dmatminus_SADDLE_2_TYPE_1;
      A_dmatmul = AMG_A_dmatmul_SADDLE_2_TYPE_1;
      A_dmatminus =  AMG_A_dmatminus_SADDLE_2_TYPE_1;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case SADDLE_2_TYPE_2 : 
      dmatmul = AMG_dmatmul_SADDLE_2_TYPE_2; 
      dmatminus =  AMG_dmatminus_SADDLE_2_TYPE_2; 
      A_dmatmul = AMG_A_dmatmul_SADDLE_2_TYPE_1;
      A_dmatminus =  AMG_A_dmatminus_SADDLE_2_TYPE_1;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_2;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[2]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case SADDLE_2_TYPE_3 : 
      dmatmul = AMG_dmatmul_SADDLE_2_TYPE_3; 
      dmatminus =  AMG_dmatminus_SADDLE_2_TYPE_3; 
      A_dmatmul = AMG_A_dmatmul_SADDLE_2_TYPE_3;   
      A_dmatminus =  AMG_A_dmatminus_SADDLE_2_TYPE_3; 
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1;
      A_length[0] = 2*AMG_MATRIX_N(A_in);         /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case SADDLE_2_TYPE_4 : 
      dmatmul = AMG_dmatmul_SADDLE_2_TYPE_4; 
      dmatminus =  AMG_dmatminus_SADDLE_2_TYPE_4; 
      A_dmatmul = AMG_A_dmatmul_SADDLE_2_TYPE_4;   
      A_dmatminus =  AMG_A_dmatminus_SADDLE_2_TYPE_4; 
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_2;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2;
      A_length[0] = 2*AMG_MATRIX_N(A_in);         /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[2]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
      dmatmul = AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_1 ; 
      dmatminus =  AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_1;
      A_dmatmul = AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_1;
      A_dmatminus =  AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_1;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
      dmatmul = AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_2 ; 
      dmatminus =  AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_2;
      A_dmatmul = AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_1;
      A_dmatminus =  AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_1;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_2;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[2]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
      dmatmul = AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_3 ; 
      dmatminus =  AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_3;
      A_dmatmul = AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_3;
      A_dmatminus =  AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_3;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
     case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
      dmatmul = AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_4; 
      dmatminus =  AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_4;
      A_dmatmul = AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_4;
      A_dmatminus =  AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_4;
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_2;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[2]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
      dmatmul = AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1 ; 
      dmatminus =  AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1;
      A_dmatmul = AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1;
      A_dmatminus =  AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1;
      B_dmatmul = AMG_B_dmatmul_SADDLE_3_TYPE_1;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_3_TYPE_1;
      A_length[0] = 3*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[0]);           /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 3;
      break;
   case SADDLE_2_TYPE_2_MORTAR  : 
      dmatmul = AMG_dmatmul_SADDLE_2_TYPE_2_MORTAR ; 
      dmatminus =  AMG_dmatminus_SADDLE_2_TYPE_2_MORTAR ; 
      A_dmatmul = AMG_A_dmatmul_SADDLE_2_TYPE_1;   
      A_dmatminus =  AMG_A_dmatminus_SADDLE_2_TYPE_1; 
      B_dmatmul = AMG_B_dmatmul_SADDLE_2_TYPE_2_MORTAR ;
      B_Trans_dmatmul = AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2_MORTAR ;
      A_length[0] = 2*AMG_MATRIX_N(A_in);           /* length of A-components solution vector */
      B_length[0] = AMG_MATRIX_M(B[2])+2*AMG_MATRIX_M(B[4]); /* length of B-components solution vector */
      for (k=1;k<AMG_MAX_LEVELS;k++)
        B_length[k] = 0;                          /* length of B-components solution vector */
      blocks = 2;
      break;
    }
  if (sc->solver==AMG_EXACT)
    {
      AMG_Print("MESSAGE : Exact solver for saddle point problems not implemented !!!\n");
      AMG_Print("MESSAGE : Solver changed to flexible GMRES (AMG_GMRES_FLEX)\n");     
      sc->solver = AMG_GMRES_FLEX;
    }
  if ((sc->amg_prec_it==0)||(sc->preconditioner == AMG_NO_PRECONDITIONER))           /* no preconditioner */    
    {
      depth=0;
      A[0]=A_in;
      sc->amg_prec_it=0;
      sc->preconditioner = AMG_NO_PRECONDITIONER; 
      return(0);
    }
  if  ((sc->preconditioner==AMG_DJAC)
       ||(sc->preconditioner==AMG_SOR)
       ||(sc->preconditioner==AMG_SSOR)
       ||(sc->preconditioner==AMG_ILU)
       ||(sc->preconditioner==AMG_ILUT))
    {
      AMG_Print("MESSAGE : Scalar type preconditioner does not work for\n");
      AMG_Print("MESSAGE : saddle point systems.\n");                             
      AMG_Print("MESSAGE : Preconditioner changed to AMG_SCHUR_GMRES\n");
      sc->preconditioner = AMG_SCHUR_GMRES;
    }

      
  switch (sc->preconditioner)
    {
    case AMG_SCHUR_COMPLEMENT :             /* all types of schur complement */
    case AMG_SCHUR_CG :                     /* preconditioners */
    case AMG_SCHUR_GMRES :
    case AMG_SCHUR_GMRES_BCGS :
  
      switch (sc->schur_inv_of_A)           /* how to approximate the inverse of A */
        {
        case AMG_MGC :                      /* with AMG */
          *prepare_amg=1;
          MGC_dmatminus = A_dmatminus;
          depth=AMG_BuildHierarchy(cc,A_in,A,G);
          if (depth<0) 
            {
              AMG_Print("Could not set up coarse grid matrices\n");
              exit(4711);
            }
          for (k=1;k<=depth;k++)
            {
              A_length[k] = blocks*AMG_MATRIX_N(A[k]);
              A[k]->level = k;
              A[k]->m=A[k]->n;
            }
          if (sc->solver==AMG_GMRES_LEFT)
            {
              AMG_Print("MESSAGE : Static preconditioned GMRES with multigrid will\n");
              AMG_Print("MESSAGE : in general not converge !!!\n");                             
              AMG_Print("MESSAGE : Solver changed to flexible GMRES (AMG_GMRES_FLEX)\n");     
              sc->solver = AMG_GMRES_FLEX;
            }
          switch(sc->system_type)
          {
            case SADDLE_2_TYPE_1: 
            case SADDLE_2_TYPE_2: 
            case SADDLE_2_TYPE_3: 
            case SADDLE_2_TYPE_4: 
            case BRAESS_SARAZIN_SADDLE_2_TYPE_1:
            case BRAESS_SARAZIN_SADDLE_2_TYPE_2:
            case BRAESS_SARAZIN_SADDLE_2_TYPE_3:
            case BRAESS_SARAZIN_SADDLE_2_TYPE_4:
              restriction = pc_restrict_2d;
              interpolation = pc_prolongate_auto_2d;
              break;
            case SADDLE_3_TYPE_1: 
            case SADDLE_3_TYPE_2: 
            case SADDLE_3_TYPE_3: 
            case SADDLE_3_TYPE_4: 
            case BRAESS_SARAZIN_SADDLE_3_TYPE_1:
            case BRAESS_SARAZIN_SADDLE_3_TYPE_2:
            case BRAESS_SARAZIN_SADDLE_3_TYPE_3:
            case BRAESS_SARAZIN_SADDLE_3_TYPE_4:
              restriction = pc_restrict_3d;
              interpolation = pc_prolongate_auto_3d;
              break;
            default :
              restriction = pc_restrict;
              interpolation = pc_prolongate_auto;
              break;
          }
          if (sc->preconditioner==AMG_SCHUR_GMRES_BCGS)
            prepare_schur_bcgs_solve(sc,A_length,B_length,depth);
          break;
          
        default:
          depth=0;
          A[0]=A_in;
        }
      break;                           /* schur complement precs finished */
        
    case AMG_MGC :                     /* coupled amg */
      *prepare_amg=1;
      MGC_dmatminus =  AMG_dmatminus_gal_SADDLE_2_TYPE_1;
      Build_B_Block_for_Dirichlet(sc,A_in,B);
      BuildDiagSchurComplement(sc,A_in,B) ;
      depth=AMG_BuildHierarchy_Saddle (cc,A_in,A,G,schur_matrix[0],
                                       schur_matrix,G_schur,B);
      if (depth<0) 
        {
          AMG_Print("Could not set up coarse grid matrices\n");
          exit(4711);
        }
      for (k=1;k<=depth;k++)                  /* no of dof on coarse grids */
        {
          A_length[k] = blocks*AMG_MATRIX_N(A[k]);
          B_length[k] = AMG_MATRIX_N(schur_matrix[k]);
          A[k]->level = k;
          A[k]->m=A[k]->n;
          AMG_MATRIX_B(schur_matrix[k])=1;
          printf("%d %d \n", A_length[k], B_length[k]);
        }
      if ((sc->solver==AMG_GMRES_LEFT)||(sc->solver==AMG_GMRES_RIGHT))  
        {                                     /* no appropriate combinations */
          AMG_Print("MESSAGE : Static preconditioned GMRES with multigrid will\n");
          AMG_Print("MESSAGE : in general not converge !!!\n");                             
          AMG_Print("MESSAGE : Solver changed to flexible GMRES (AMG_GMRES_FLEX)\n");     
          sc->solver = AMG_GMRES_FLEX;
        }
      restriction = pc_restrict_saddle_2d;
      interpolation = pc_prolongate_auto_saddle_2d;
      prepare_coupled_mgc(sc,A_length,B_length,depth);
      if (sc->smoother!=AMG_BRAESAR)
        {
          AMG_Print("MESSAGE : prescribed smoother will not work for coupled amg.\n");
          AMG_Print("MESSAGE : smoother changed to braess-sarazin smoother.\n");
          sc->smoother=AMG_BRAESAR;
        }
      if (sc->coarse_smoother!=AMG_BRAESAR)
        {
          AMG_Print("MESSAGE : prescribed coarse smoother will not work for coupled amg.\n");
          AMG_Print("MESSAGE : coarse smoother changed to braess-sarazin smoother.\n");
          sc->coarse_smoother=AMG_BRAESAR;
        }       
      break;
    default:
      depth=0;
      A[0]=A_in;
      
    }     

  
  return(0);
}

static int solver_build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, 
                         AMG_MATRIX *A_in,int user_solver)
{
  int k,solver,A_length[AMG_MAX_LEVELS],B_length[AMG_MAX_LEVELS];
  int i,prepare_amg;


  prepare_amg=0;
  A_in->level = 0;
  switch(sc->system_type)
    {
    case SCALAR1 : 
    case SCALAR2 : 
    case SCALAR3 : 
    case SCALAR6 : 
    case SCALAR_VAL_LAZ_2 : 
      prepare_scalar_system(sc,cc,A_in,A_length,B_length,&prepare_amg); 
      break;
    case SADDLE_1 :
    case SADDLE_1_TYPE_1 : 
    case SADDLE_2_TYPE_1 : 
    case SADDLE_2_TYPE_2 : 
    case SADDLE_2_TYPE_3 : 
    case SADDLE_2_TYPE_4 : 
    case SADDLE_2_TYPE_2_MORTAR  : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
      prepare_saddle_system(sc,cc,A_in,A_length,B_length,&prepare_amg);
     break;
    default : 
      AMG_Print("Unknown system type !!!\n");
      exit (4711);
      break;
    }

  solver=sc->solver;
  switch (solver)                               /* some consistency checks */                       
    {                                         
    case AMG_GMRES_FLEX :
      if (sc->amg_prec_it==0)
        {                                             /* no appropriate combinations */
          AMG_Print("MESSAGE : flexible GMRES without preconditioning does not work!!!\n");
          AMG_Print("MESSAGE : Solver changed to GMRES (AMG_GMRES_LEFT)\n");     
          solver = AMG_GMRES_LEFT;                    /* change solver */
          sc->solver = AMG_GMRES_LEFT;
        }
      break;
    }

  switch (solver)                             /* allocate matrices and vectors for the */
    {                                         /* different solvers */ 
    case AMG_GMRES_LEFT:
    case AMG_GMRES_RIGHT:
      prepare_gmres_solve(sc,A_length,B_length,depth);
      break;
    case AMG_GMRES_FLEX:
      prepare_gmres_flexible_solve(sc,A_length,B_length,depth);
      break;
    case AMG_CG:
      prepare_cg_solve(sc,A_length,B_length,depth);      
      break;
    case AMG_BCGS:
      prepare_bcgs_solve(sc,A_length,B_length,depth);      
      break;
    case AMG_MIXED_BCGS_CGS:
      prepare_mixed_bcgs_cgs_solve(sc,A_length,B_length,depth); 
      break;
    case AMG_LS:
      prepare_ls_solve(sc,A_length,B_length,depth);
      break;
    case AMG_EXACT:
      prepare_exact_solve(sc,A_length,B_length,depth);
      sc->preconditioner = AMG_DJAC;
      sc->step_length_control_fine = 0;
      sc->step_length_control_all = 0;
      sc->omega[0] = 1.0;
      break;
    case AMG_LCD:
      prepare_lcd_solve(sc,A_length,B_length,depth);
      break;
    default:
      AMG_Print("ERROR: invalid solver\n");
      return(AMG_FATAL);
    }
  if (solver==AMG_EXACT)  return(AMG_OK);
                                             /* allocate vectors for step length control */

  if (sc->step_length_control_all)                              
    prepare_steplength_control(A_length,B_length,depth);
  else
    if (sc->step_length_control_fine) 
      prepare_steplength_control(A_length,B_length,0);
   
                                            /* preconditioner preprocess */

  for (k=0; k<=depth; k++) M[k]=A[k];       /* default, change only for ILU */

  switch (sc->preconditioner)
    {
    case AMG_DJAC:
      preconditioner = jac;
      preconditioner_trans = jac;
      break;
    case AMG_SOR:
      preconditioner = sor;
      preconditioner_trans = sor_trans;
      break;
    case AMG_SSOR:
      /* if no step length control */
      if ((!sc->step_length_control_fine) && (!sc->step_length_control_all)
	  && (sc->omega[0]==1.0) && (sc->omega_coarse[0]==1.0))
      {
	preconditioner = ssor;
	preconditioner_trans = ssor_trans;
      }
      else 
	preconditioner = ssor_slc;
      break;
    case AMG_ILU:
      preconditioner = ilu;
      preconditioner_trans = ilu_trans;
      prepare_ilu(sc,A,M,0,0);
      break;
    case AMG_ILUT:
      preconditioner = ilu;
      M[0] = ILUTDecomposition(sc,A[0]);                                    
      break;
    case AMG_MGC:
      preconditioner = mgc;
      break;
    case AMG_SCHUR_COMPLEMENT:
      preconditioner = schur_complement_fixed;
      prepare_schur_complement_fixed(A_length,B_length,depth);
      break;
    case AMG_SCHUR_CG:
      preconditioner = schur_complement_cg;
      prepare_schur_complement_cg(A_length,B_length,depth);
     break;
    case AMG_SCHUR_GMRES:
      preconditioner = schur_complement_gmres;
      prepare_schur_complement_gmres(sc,A_length,B_length,depth);
      break;
    case AMG_SCHUR_GMRES_BCGS:
      preconditioner = schur_complement_gmres_bcgs;
      prepare_schur_complement_gmres(sc,A_length,B_length,depth);
      prepare_schur_bcgs_solve(sc,A_length,B_length,depth);
    case AMG_NO_PRECONDITIONER:
      preconditioner = NULL;
      break;
    default:
      AMG_Print("invalid preconditioner\n");
      exit(4711);
    }

  if (  !((sc->system_type)==SCALAR1 || (sc->system_type)==SCALAR2 || (sc->system_type)==SCALAR3 || (sc->system_type)==SCALAR6 )       
      && ( (sc->preconditioner)==AMG_SCHUR_COMPLEMENT || (sc->preconditioner)==AMG_SCHUR_CG || (sc->preconditioner)==AMG_SCHUR_GMRES 
       || (sc->preconditioner)==AMG_SCHUR_GMRES_BCGS) )
    switch(sc->schur_inv_of_A)
      {
      case AMG_DJAC :
        schur_preconditioner = jac;
        break;
      case AMG_SOR :
        schur_preconditioner = sor;
        break;
      case AMG_SSOR :
        schur_preconditioner = ssor;
        break;
      case AMG_ILU :
        schur_preconditioner = ilu;
        prepare_ilu(sc,A,M,0,0);
        break;
      case AMG_ILUT :
        schur_preconditioner = ilu;
        M[0] = ILUTDecomposition(sc,A[0]);                                    
        break;
      case AMG_MGC :
        schur_preconditioner = mgc;
        break;
      case AMG_EX :
        M[0]=prepare_ex(A[0]);
        if (M[0]==NULL) {
          AMG_Print("error in prepare_ex\n");
          exit(4711);
        }
        schur_preconditioner = ex;
        break;
      default:
        AMG_Print("invalid schur preconditioner\n");
        exit(4711);
      }
  if (!prepare_amg) return(AMG_OK);
                                              /* smoother preprocess */
  switch (sc->smoother)
    {
    case AMG_DJAC:
      smoother = jac;
      break;
    case AMG_SOR:
      smoother = sor;
      break;
    case AMG_SSOR:
      smoother = ssor;
      break;
    case AMG_ILU:
      smoother = ilu;
      prepare_ilu(sc,A,M,0,depth-1);
      break;
    case AMG_ILUT:
      smoother = ilu;
      for (k=0; k<=depth; k++)
        M[k] = ILUTDecomposition(sc,A[k]);                                    
      break;
    case AMG_BRAESAR:
      smoother =  braess_sarazin_smoother;
      prepare_braess_sarazin_smoother(sc,A_length,B_length,0);
      break;
    default:
      AMG_Print("invalid smoother\n");
      exit(4711);
    }
                                               /* coarse smoother preprocess */
  switch (sc->coarse_smoother)
    {
    case AMG_DJAC:
      coarse_smoother = jac;
      break;
    case AMG_SOR:
      coarse_smoother = sor;
      break;
    case AMG_SSOR:
      coarse_smoother = ssor;
      break;
    case AMG_ILU:
      coarse_smoother = ilu;
      prepare_ilu(sc,A,M,depth,depth);
      break;
    case AMG_ILUT:
      coarse_smoother = ilu;
      M[depth] = ILUTDecomposition(sc,A[depth]);
      break;
    case AMG_EX:
      M[depth]=prepare_ex(A[depth]);
      if (M[depth]==NULL) {
        AMG_Print("error in prepare_ex\n");
        exit(4711);
      }
      coarse_smoother = ex;
      break;
    case AMG_BRAESAR :
      coarse_smoother = braess_sarazin_smoother;
      break;
    default:
      AMG_Print("invalid coarse smoother\n");
      exit(4711);
    }
  switch (sc->smoothing_steps)
    {
    case CONSTANT:
      for (k=1;k<AMG_MAX_LEVELS;k++)
        {
          sc->n1[k]=sc->n1[0];
          sc->n2[k]=sc->n2[0];
        }
      break;
    case PLUS_CONSTANT:
      for (k=1;k<AMG_MAX_LEVELS;k++)
        {
          sc->n1[k]=sc->n1[k-1]+sc->n1_param;
          sc->n2[k]=sc->n2[k-1]+sc->n2_param;
        }
      break;
    case TIMES_CONSTANT:
      for (k=1;k<AMG_MAX_LEVELS;k++)
        {
          sc->n1[k]=sc->n1[k-1]*sc->n1_param;
          sc->n2[k]=sc->n2[k-1]*sc->n2_param;
        }
      break;
    case SQUARED:
      for (k=1;k<AMG_MAX_LEVELS;k++)
        {
          sc->n1[k]=sc->n1[k-1]*sc->n1[k-1];
          sc->n2[k]=sc->n2[k-1]*sc->n2[k-1];
        }
      break;
    default :
      for (k=1;k<AMG_MAX_LEVELS;k++)
        {
          sc->n1[k]=sc->n1[0];
          sc->n2[k]=sc->n2[0];
        }
      break;
    }
  if (sc->gamma>0)
    for (k=0;k<AMG_MAX_LEVELS;k++)        
      mgc_recursion[k] = sc->gamma;
  else                /* F -- cycle */
    for (k=0;k<AMG_MAX_LEVELS;k++)        
      mgc_recursion[k] = 2;

  return(AMG_OK);
}

static int ReportSolverParameters(AMG_SolverContext *sc)
{
  int i;

  AMG_Print("\n");
  AMG_Print("***** SOLVER PARAMETERS *****\n");
  AMG_Print("\n");
  AMG_Print("hostname : ");
  gethostname(buf,80);
  AMG_Print(buf);
  AMG_Print("\n");
  AMG_Print("system type : ");
  switch(sc->system_type)
  {
    case SCALAR1 : 
      sprintf(buf,"SCALAR1 (= %d)\n",SCALAR1);
      AMG_Print(buf);
      break;
    case SCALAR2 : 
      sprintf(buf,"SCALAR2 : 1D problem with 2 right hand sides (= %d)\n",SCALAR2);
      AMG_Print(buf);
      break;
    case SCALAR3 : 
      sprintf(buf,"SCALAR3 : 1D problem with 3 right hand sides (= %d)\n",SCALAR3);
      AMG_Print(buf);
      break;
    case SCALAR6 : 
      sprintf(buf,"SCALAR6 : 1D problem with 6 right hand sides (= %d)\n",SCALAR6);
      AMG_Print(buf);
      break;
    case SADDLE_1_TYPE_1 :
      sprintf(buf,"SADDLE_1_TYPE_1 : 1D saddle point, type 1 (= %d)\n",SADDLE_1_TYPE_1);
      AMG_Print(buf);
      break;
    case SADDLE_2_TYPE_1 :
      sprintf(buf,"SADDLE_2_TYPE_1 : 2D saddle point, type 1 (= %d)\n",SADDLE_2_TYPE_1);
      AMG_Print(buf);
      break;
    case SADDLE_2_TYPE_2 :
      sprintf(buf,"SADDLE_2_TYPE_2 : 2D saddle point, type 2 (= %d)\n",SADDLE_2_TYPE_2);
      AMG_Print(buf);
      break;
    case SADDLE_2_TYPE_2_MORTAR  :
      sprintf(buf,"SADDLE_2_TYPE_2_MORTAR  : 2D saddle point for mortar discretization, type 2 (= %d)\n",
              SADDLE_2_TYPE_2_MORTAR );
      AMG_Print(buf);
      break;
    case SADDLE_2_TYPE_4 :
      sprintf(buf,"SADDLE_2_TYPE_4 : 2D saddle point, type 4 (= %d)\n",SADDLE_2_TYPE_4);
      AMG_Print(buf);
      break;
    case SADDLE_2_TYPE_3 :
      sprintf(buf,"SADDLE_2_TYPE_3 : 2D saddle point, type 3 (= %d)\n",SADDLE_2_TYPE_3);
      AMG_Print(buf);
      break;
    case SCALAR_VAL_LAZ_2 :
      sprintf(buf,"SCALAR_VAL_LAZ_2 : 2D Vassilevski/Lazarov preconditioner (= %d)\n",SCALAR_VAL_LAZ_2);
      AMG_Print(buf);
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 :
      sprintf(buf,"BRAESS_SARAZIN_SADDLE_2_TYPE_1 : Braess-Sarazin 2D saddle point, type 1 (= %d)\n"
              ,BRAESS_SARAZIN_SADDLE_2_TYPE_1);
      AMG_Print(buf);
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 :
      sprintf(buf,"BRAESS_SARAZIN_SADDLE_2_TYPE_2 : Braess-Sarazin 2D saddle point, type 2 (= %d)\n"
              ,BRAESS_SARAZIN_SADDLE_2_TYPE_2);
      AMG_Print(buf);
      break;
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 :
      sprintf(buf,"BRAESS_SARAZIN_SADDLE_2_TYPE_3 : Braess-Sarazin 2D saddle point, type 3 (= %d)\n"
              ,BRAESS_SARAZIN_SADDLE_2_TYPE_3);
      AMG_Print(buf);
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 :
      sprintf(buf,"BRAESS_SARAZIN_SADDLE_2_TYPE_4 : Braess-Sarazin 2D saddle point, type 4 (= %d)\n"
              ,BRAESS_SARAZIN_SADDLE_2_TYPE_4);
      AMG_Print(buf);
      break;
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 :
      sprintf(buf,"BRAESS_SARAZIN_SADDLE_3_TYPE_1 : Braess-Sarazin 3D saddle point, type 1 (= %d)\n"
              ,BRAESS_SARAZIN_SADDLE_3_TYPE_1);
      AMG_Print(buf);
      break;
    default :
       AMG_Print("UNKWOWN !!!\n");
      break;     
    }
  AMG_Print("solver : ");
  switch(sc->solver)
    {
    case AMG_LS :
      sprintf(buf,"AMG_LS - linear iteration (= %d)\n",AMG_LS);
      AMG_Print(buf);
      break;
    case AMG_CG :
      sprintf(buf,"AMG_CG - conjugate gradients (= %d)\n",AMG_CG);
      AMG_Print(buf);
      break;
    case AMG_BCGS :
      sprintf(buf,"AMG_BCGS - biconjugate gradients stabilized (= %d)\n",AMG_BCGS);
      AMG_Print(buf);
      break;
    case AMG_GMRES_LEFT :
      sprintf(buf,"AMG_GMRES_LEFT - left preconditioned GMRES (= %d)\n",AMG_GMRES_LEFT);
      AMG_Print(buf);
      break;
    case AMG_GMRES_RIGHT :
      sprintf(buf,"AMG_GMRES_RIGHT - right preconditioned GMRES (= %d)\n",AMG_GMRES_RIGHT);
      AMG_Print(buf);
      break;
    case AMG_GMRES_FLEX :
      sprintf(buf,"AMG_GMRES_FLEX - flexible GMRES (= %d)\n",AMG_GMRES_FLEX);
      AMG_Print(buf);
      break;
    case AMG_EXACT :
      sprintf(buf,"AMG_EXACT - exact (= %d)\n",AMG_EXACT);
      AMG_Print(buf);
      break;
    case AMG_MIXED_BCGS_CGS :
      sprintf(buf,"AMG_MIXED_BCGS_CGS - mixed cgs/bcgs method (= %d)\n",AMG_MIXED_BCGS_CGS);
      AMG_Print(buf);
      break;
    case AMG_LCD :
      sprintf(buf,"AMG_LCD - left conjugate direction (= %d)\n",AMG_LCD);
      AMG_Print(buf);
      break;
    default :
       AMG_Print("UNKWOWN !!!\n");
      break;     
    }       

  sprintf(buf,"maximal number of iterations : %d\n",sc->maxit);
  AMG_Print(buf);
  AMG_Print("execute exactly maximal number of iterations : ");
  if (sc->ex_maxit)
    AMG_Print("yes\n");
  else
    AMG_Print("no\n");
   sprintf(buf,"minimal number of iterations : %d\n",sc->minit);
  AMG_Print(buf);
  
  sprintf(buf,"reduction factor for norm of initial residual : %g\n",sc->red_factor);
  AMG_Print(buf);
  sprintf(buf,"stopping criterion for norm of residual : %g\n",sc->res_norm_min);
  AMG_Print(buf);
  AMG_Print("row equilibration : ");
  if (sc->row_equilibration)
    AMG_Print("yes\n");
  else
    AMG_Print("no\n");

  AMG_Print("preconditioner : ");
  switch(sc->preconditioner)
    {
    case AMG_DJAC :
      sprintf(buf,"AMG_DJAC - Jacobi (= %d)\n",AMG_DJAC);
      AMG_Print(buf);
      break;
    case AMG_SOR :
      sprintf(buf,"AMG_SOR - succesive over relaxation (= %d)\n",AMG_SOR);
      AMG_Print(buf);
      break;
    case AMG_SSOR :
      sprintf(buf,"AMG_SSOR - symmetric succesive over relaxation (= %d)\n",AMG_SSOR);
      AMG_Print(buf);
      break;
    case AMG_ILU :
      sprintf(buf,"AMG_ILU - incomplete lu factorization (= %d)\n",AMG_ILU);
      AMG_Print(buf);
      break;
    case AMG_MGC :
      sprintf(buf,"AMG_MGC - algebraic multigrid cycle (= %d)\n",AMG_MGC);
      AMG_Print(buf);
      break;
    case AMG_ILUT :
      sprintf(buf,"AMG_ILUT - incomplete lu factorization with fillin and threshhold (= %d)\n",AMG_ILUT);
      AMG_Print(buf);
      break;
    case AMG_SCHUR_COMPLEMENT :
      sprintf(buf,"AMG_SCHUR_COMPLEMENT - fixed point schur complement iteration (= %d)\n",AMG_SCHUR_COMPLEMENT);
      AMG_Print(buf);
      break;
    case AMG_SCHUR_CG :
      sprintf(buf,"AMG_SCHUR_CG - schur complement iteration with cg (= %d)\n",AMG_SCHUR_CG);
      AMG_Print(buf);
      break;
    case AMG_SCHUR_GMRES :
      sprintf(buf,"AMG_SCHUR_GMRES - schur complement iteration with gmres (= %d)\n",AMG_SCHUR_GMRES);
      AMG_Print(buf);
      break;
   case AMG_EXACT:
      sprintf(buf,"AMG_EXACT - exact (= %d)\n",AMG_EXACT);
      AMG_Print(buf);
      break;
    default :
      AMG_Print("UNKWOWN !!!\n");
      break;     
    }       
  sprintf(buf,"number of preconditioner iterations : %d\n",sc->amg_prec_it);
  AMG_Print(buf);
  sprintf(buf,"reduction factor for preconditioner iteration : %g\n",sc->amg_prec_red_factor);
  AMG_Print(buf);

  AMG_Print("multigrid cycle : ");
  switch(sc->gamma)
    {
    case 1 :
      sprintf(buf,"V--cycle (gamma = %d)\n",sc->gamma);
      AMG_Print(buf);
      break;
    case 2 :
      sprintf(buf,"W--cycle (gamma = %d)\n",sc->gamma);
      AMG_Print(buf);
      break;
    default :
      if (sc->gamma<1)
        {
          sprintf(buf,"F--cycle (gamma = %d)\n",sc->gamma);
          AMG_Print(buf);
        }
      else
        {
          sprintf(buf,"no--name--cycle (gamma = %d)\n",sc->gamma);
          AMG_Print(buf);
        }        
      break;     
    }       

  AMG_Print("smoother : ");
  switch(sc->smoother)
    {
    case AMG_DJAC :
      sprintf(buf,"AMG_DJAC - Jacobi (= %d)\n",AMG_DJAC);
      AMG_Print(buf);
      break;
    case AMG_SOR :
      sprintf(buf,"AMG_SOR - succesive over relaxation (= %d)\n",AMG_SOR);
      AMG_Print(buf);
      break;
    case AMG_SSOR :
      sprintf(buf,"AMG_SSOR - symmetric succesive over relaxation (= %d)\n",AMG_SSOR);
      AMG_Print(buf);
      break;
    case AMG_ILU :
      sprintf(buf,"AMG_ILU - incomplete lu factorization (= %d)\n",AMG_ILU);
      AMG_Print(buf);
      break;
    case AMG_ILUT :
      sprintf(buf,"AMG_ILUT - incomplete lu factorization with fillin and threshhold (= %d)\n",AMG_ILUT);
      AMG_Print(buf);
      break;
    case AMG_BRAESAR :  
      sprintf(buf,"AMG_BRAESAR - braess-sarazin smoother for saddle point problems (= %d)\n",AMG_BRAESAR);
      AMG_Print(buf);
      break;
     default :
      AMG_Print("UNKWOWN !!!\n");
      break;     
    }
       
  sprintf(buf,"number of pre smoothing steps on finest level : %d\n",sc->n1[0]);
  AMG_Print(buf);
  if (sc->smoothing_steps==CONSTANT )
    {
      sprintf(buf,"number of pre smoothing steps on all other levels : %d\n",sc->n1[0]);
      AMG_Print(buf);
    }
  else
    {
      for (i=1;i<AMG_MAX_LEVELS;i++)
        {
          sprintf(buf,"l %d %d  ",i,sc->n1[i]);
          AMG_Print(buf);
        }
      AMG_Print("\n");
    }
  sprintf(buf,"number of post smoothing steps on finest level : %d\n",sc->n2[0]);
  AMG_Print(buf);
  if (sc->smoothing_steps==CONSTANT )
    {
      sprintf(buf,"number of post smoothing steps on all other levels : %d\n",sc->n2[0]);
      AMG_Print(buf);
    }
  else
    {
      for (i=1;i<AMG_MAX_LEVELS;i++)
        {
          sprintf(buf,"l %d %d  ",i,sc->n2[i]);
          AMG_Print(buf);
        }
      AMG_Print("\n");
    }
  sprintf(buf,"reduction factor for smoothing iteration : %g\n",sc->smoother_red_factor);
  AMG_Print(buf);
  sprintf(buf,"damping factor for smoothing iteration : %g\n",sc->omega[0]);
  AMG_Print(buf);

  AMG_Print("coarse smoother : ");
  switch(sc->coarse_smoother)
    {
    case AMG_DJAC :
      sprintf(buf,"AMG_DJAC - Jacobi (= %d)\n",AMG_DJAC);
      AMG_Print(buf);
      break;
    case AMG_SOR :
      sprintf(buf,"AMG_SOR - succesive over relaxation (= %d)\n",AMG_SOR);
      AMG_Print(buf);
      break;
    case AMG_SSOR :
      sprintf(buf,"AMG_SSOR - symmetric succesive over relaxation (= %d)\n",AMG_SSOR);
      AMG_Print(buf);
      break;
    case AMG_ILU :
      sprintf(buf,"AMG_ILU - incomplete lu factorization (= %d)\n",AMG_ILU);
      AMG_Print(buf);
      break;
    case AMG_ILUT :
      sprintf(buf,"AMG_ILUT - incomplete lu factorization with fillin and threshhold (= %d)\n",AMG_ILUT);
      AMG_Print(buf);
      break;
   case AMG_BRAESAR :  
      sprintf(buf,"AMG_BRAESAR - braess-sarazin smoother for saddle point problems (= %d)\n",AMG_BRAESAR);
      AMG_Print(buf);
      break;
     case AMG_EXACT:
      sprintf(buf,"AMG_EXACT : exact (= %d)\n",AMG_EXACT);
      AMG_Print(buf);
      break;
    default :
      AMG_Print("UNKWOWN !!!\n");
      break;     
    }
       
  sprintf(buf,"maximal number of coarse smoothing steps : %d\n",sc->coarse_maxit);
  AMG_Print(buf);
  sprintf(buf,"reduction factor for coarse smoothing iteration : %g\n",sc->coarse_red_factor);
  AMG_Print(buf);
  sprintf(buf,"damping factor for smoothing iteration : %g\n",sc->omega_coarse[0]);
  AMG_Print(buf);

  sprintf(buf,"damping factor for interpolation : %g\n",sc->omega_p[0]);
  AMG_Print(buf);

  AMG_Print("step length control on the finest level (all preconditioners) : ");
  if (sc->step_length_control_fine)
    AMG_Print("yes\n");
  else
    AMG_Print("no\n");
  AMG_Print("step length control on all levels (multigrid) : ");
  if (sc->step_length_control_all)
    AMG_Print("yes\n");
  else
    AMG_Print("no\n");

  sprintf(buf,"restart in gmres : %d\n",sc->gmres_restart);
  AMG_Print(buf);
  sprintf(buf,"omega in jacobi/sor/ssor : %g\n",sc->sor_omega);
  AMG_Print(buf);
  sprintf(buf,"beta in ilu : %g\n",sc->ilu_beta);
  AMG_Print(buf);
  sprintf(buf,"absolute fill in in ilut : %d\n",sc->ilut_absolute_fillin);
  AMG_Print(buf);
  sprintf(buf,"relative fill in in ilut : %g \n",sc->ilut_relative_fillin);
  AMG_Print(buf);
  sprintf(buf,"dropping tolerance in ilut : %g\n",sc->ilut_tol);
  AMG_Print(buf);
  AMG_Print("sort algorithm in ilut : ");
  switch(sc->ilut_sort)
    {
    case  ILUT_QUICK_SPLIT_0 :      
      sprintf(buf,"quick split, version 0 (= %d)\n",ILUT_QUICK_SPLIT_0);
      AMG_Print(buf);
      break;
    case  ILUT_QUICK_SPLIT_1 :      
      sprintf(buf,"quick split, version 1 (= %d)\n",ILUT_QUICK_SPLIT_1);
      AMG_Print(buf);
      break;
    case  ILUT_QUICK_SPLIT_2 :      
      sprintf(buf,"quick split, version 2 (= %d)\n",ILUT_QUICK_SPLIT_2);
      AMG_Print(buf);
      break;
     default :
       AMG_Print("UNKWOWN !!!\n");
      break;     
    }

  AMG_Print("approximate inverse of A in Schur iteration by : ");
  switch(sc->schur_inv_of_A)
    {
    case AMG_DJAC :
      sprintf(buf,"AMG_DJAC - Jacobi (= %d)\n",AMG_DJAC);
      AMG_Print(buf);
      break;
    case AMG_SOR :
      sprintf(buf,"AMG_SOR - SOR (= %d)\n",AMG_SOR);
      AMG_Print(buf);
      break;
    case AMG_SSOR :
      sprintf(buf,"AMG_SSOR - SSOR (= %d)\n",AMG_SSOR);
      AMG_Print(buf);
      break;
    case AMG_ILU :
      sprintf(buf,"AMG_ILU - ILU (= %d)\n",AMG_ILU);
      AMG_Print(buf);
      break;
    case AMG_MGC :
      sprintf(buf,"AMG_MGC - AMG (= %d)\n",AMG_MGC);
      AMG_Print(buf);
      break;
    case AMG_ILUT :
      sprintf(buf,"AMG_ILUT - ILUT (= %d)\n",AMG_ILUT);
      AMG_Print(buf);
      break;
   case AMG_EXACT :
      sprintf(buf,"AMG_EXACT - exact (= %d)\n",AMG_EXACT);
      AMG_Print(buf);
      break;
    default :
      AMG_Print("UNKWOWN !!!\n");
      break;     
    }       

  sprintf(buf,"number of iterations in approximating inverse of A in Schur complement : %d\n",sc->schur_inv_of_A_maxit);
  AMG_Print(buf);
  sprintf(buf,"damping factor for Schur complement fixed point iteration : %g\n",sc->schur_iteration_damp);
  AMG_Print(buf);
  sprintf(buf,"maximal number of iterations for Schur complement solver : %d\n",sc->schur_iteration_maxit);
  AMG_Print(buf);
  AMG_Print("step length control in Schur complement fixed point iteration : ");
  if (sc->schur_step_length_control)
    AMG_Print("yes\n");
  else
    AMG_Print("no\n");

  AMG_Print("matrix in Braess-Sarazin smoother : ");
  switch(sc->braess_sarazin_matrix)
    {
    case 0 :
      sprintf(buf,"original matrix (= 0)\n");
      AMG_Print(buf);
      break;
    case 1 :
      sprintf(buf,"idendity (= 1)\n");
      AMG_Print(buf);
      break;
    case 2 :
      sprintf(buf,"diagonal of original matrix (= 2)\n");
      AMG_Print(buf);
      break;
    case 3 :
      sprintf(buf,"ilu decomposition of original matrix (= 3)\n");
      AMG_Print(buf);
      break;
    case 4 :
      sprintf(buf,"ilut decomposition of original matrix (= 4)\n");
      AMG_Print(buf);
      break;
    default :
      AMG_Print("UNKWOWN !!!\n");
      break;     
    }       
  sprintf(buf,"damping factor in Braess-Sarazin smoother : %g\n",
          sc->braess_sarazin_alpha);
  AMG_Print(buf);

  AMG_Print("\n");
  AMG_Print("*****************************\n");
  AMG_Print("\n");
  return(0);
}
                
static int AMG_ClearSolverOverhead (AMG_SolverContext *sc, AMG_CoarsenContext *cc)
{
  int k,solver,i,clear_amg=0;
  
  solver=sc->solver;
                                  /* free coarse grid matrices */
  switch (sc->preconditioner)
    {
    case AMG_MGC:
      
      if ((sc->system_type==SCALAR1)||(sc->system_type==SCALAR3)||(sc->system_type==SCALAR2)
          ||(sc->system_type==SCALAR6))
        clear_mgc(A,G,depth);
      else
        clear_coupled_mgc(A,G,depth);
      clear_amg=1;
      break;
    default:
      /*free(G[0]->ca);
      free(G[0]->na);
      free(G[0]->la);
      free(G[0]->da);
      free(G[0]);*/          
      break;
    }

  switch (solver)                             /* free matrices and vectors for the */
    {                                         /* different solvers */ 
    case AMG_GMRES_LEFT:
    case AMG_GMRES_RIGHT:
      clear_gmres_solve(sc,depth);
      break;
    case AMG_GMRES_FLEX:
      clear_gmres_flexible_solve(sc,depth);
      break;
    case AMG_CG:
      clear_cg_solve(sc,depth);
      break;
    case AMG_BCGS:
      clear_bcgs_solve(sc,depth);
      break;
    case AMG_MIXED_BCGS_CGS:
      clear_mixed_bcgs_cgs_solve(sc,depth);
      break;
    case AMG_LS:
      clear_ls_solve(sc,depth);
      break;
    case AMG_EXACT:
      clear_exact_solve(sc);
      break;
    case AMG_LCD:
      clear_lcd_solve(sc,depth);
      break;
    default:
      AMG_Print("ERROR: invalid solver\n");
      return(AMG_FATAL);
    }
                                             /* free vectors for step length control */
  if (sc->step_length_control_all)                              
    clear_steplength_control(depth);
  else
    if (sc->step_length_control_fine) 
      clear_steplength_control(0);
 
  
/* preconditioner postprocess */

  switch (sc->preconditioner)
    {
    case AMG_ILU :
      clear_ilu(M,0,0);
      break;
    case AMG_ILUT :
      clear_ilut(M,0,0);
      break;
    case AMG_SCHUR_COMPLEMENT :
      clear_schur_complement_fixed(depth);
      break;
    case AMG_SCHUR_CG :
      clear_schur_complement_cg(depth);
      break;
    case AMG_SCHUR_GMRES :
      clear_schur_complement_gmres(sc,depth);
      break;
    }
  if ((!((sc->system_type==SCALAR1)||(sc->system_type==SCALAR2)||(sc->system_type==SCALAR3)
         ||(sc->system_type==SCALAR6)))&&
      ((sc->preconditioner==AMG_SCHUR_COMPLEMENT)
       || (sc->preconditioner==AMG_SCHUR_CG)
       || (sc->preconditioner==AMG_SCHUR_GMRES)))
    switch(sc->schur_inv_of_A)
      {
      case AMG_ILU :
        clear_ilu(M,0,0);
        break;
      case AMG_ILUT :
        clear_ilut(M,0,0);
        break;
      case AMG_MGC :
        clear_mgc(A,G,depth);
        clear_amg=1;
        break; 
      case AMG_EX :
        clear_ex(M,0);
        break;
      }
  
  if (!clear_amg)   return(AMG_OK);
                                              /* smoother postprocess */
  switch (sc->smoother)
    {
    case AMG_ILU:
      clear_ilu(M,0,depth-1);
      break;
    case AMG_ILUT:
      clear_ilut(M,0,depth-1);
      break;
    case AMG_BRAESAR:
      clear_braess_sarazin_smoother(0);
      break;
    }
                                               /* coarse smoother postprocess */
  switch (sc->coarse_smoother)
    {
    case AMG_ILU:
      clear_ilu(M,depth,depth);
      break;
    case AMG_ILUT:
      clear_ilut(M,depth,depth);
      break;
    case AMG_EX:
      clear_ex(M,depth);
      break;
    }
 
  return(AMG_OK);
}
                        
/****************************************************************************/
/*
   amg - solve system of linear equations with algebraic multigrid
   
   SYNOPSIS:
   
   int AMG (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b)
   
   PARAMETERS:
   
   DESCRIPTION:
   
   RETURN VALUE:
                                                                            */
/****************************************************************************/

int AMG(AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b0)
{
  int rv,i,j,memory[3];
 
 
//   struct mallinfo MALLINFO;
  printf("AMG main started\n"); 
  /*sprintf(buf,"elapsed_time: %g sec\n", elapsed_time);
    AMG_Print(buf);*/      

//   MALLINFO = mallinfo();
//   memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
  elapsed_time=0.0;
  start=clock();  
  /* store context */
  global_cc=cc; 
  global_sc=sc;
  
  B[0] = B0[0];
  B[1] = B0[1];
  B[2] = B0[2];
  B[3] = B0[3];
  B[4] = B0[4];
  B[5] = B0[5];
  B[6] = B0[6];
  b[0] = b0;

  A->ratr = NULL;
  A->jatr = NULL;
  A->postr = NULL;

  switch (sc->solver)
  {
    case AMG_CG:
      rv=solver_build(sc,cc,A,AMG_CG);
      break;
      
    case AMG_BCGS:
      rv=solver_build(sc,cc,A,AMG_BCGS);
      break;
      
    case AMG_LS:
      rv=solver_build(sc,cc,A,AMG_LS);
      break;
      
    case AMG_GMRES_LEFT:
      rv=solver_build(sc,cc,A,AMG_GMRES_LEFT);
      break;
    case AMG_GMRES_RIGHT:
      rv=solver_build(sc,cc,A,AMG_GMRES_LEFT);
      break;
    case AMG_GMRES_FLEX:
      rv=solver_build(sc,cc,A,AMG_GMRES_FLEX);
      break;
    case AMG_EXACT:
      rv=solver_build(sc,cc,A,AMG_EXACT);
      break;
    case AMG_MIXED_BCGS_CGS:
      rv=solver_build(sc,cc,A,AMG_MIXED_BCGS_CGS);
      break;
    case AMG_LCD:
      rv=solver_build(sc,cc,A,AMG_LCD);
      break;
    default:
      AMG_Print("solver not implemented\n");
      exit(4711);
  }
  if (sc->row_equilibration)
    EquilibrateRHS(b0,row_equilibration);

  finish = clock();
  if (finish>=start)
    elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
  else
    elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
  start = clock();

  cc->time = elapsed_time;

  if (sc->verbose>0)
  {
//     MALLINFO = mallinfo();
//     memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
//     sprintf(buf,"AMG: PREPARATION TIME %g sec, MEMORY ALLOCATED %7.3g MB\n",
//             elapsed_time,(memory[1]-memory[0])/(1048576.0));
//     AMG_Print(buf);      
  }

  if (sc->verbose>1)
    ReportSolverParameters(sc);

  residual_cnt=0;  
  for(i=0;i<10;i++) residuals[i]=0.0;
  if (sc->verbose>0)
  {
//     MALLINFO = mallinfo();
//     memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
//     sprintf(buf,"AMG: TOTAL MEMORY IN USE %7.3g MB\n", memory[2]/(1048576.0));
//     AMG_Print(buf);      
  }
  switch (global_sc->solver)
  {
    case AMG_CG:
      rv=cg_solve(x,b0);
      break;
      
    case AMG_BCGS:
      rv=bcgs_solve(x,b0);
      break;
      
    case AMG_LS:
      rv=ls_solve(x,b0);
      break;
      
    case AMG_GMRES_LEFT:
      rv=gmres_left_solve(x,b0);
      break;
      
    case AMG_GMRES_RIGHT:
      rv=gmres_right_solve(x,b0);
      break;

    case AMG_GMRES_FLEX:
      rv=gmres_flexible_solve(x,b0);
      break;

    case AMG_EXACT:
      rv=exact_solve(x,b0);
      break;
 
   case AMG_MIXED_BCGS_CGS:
      rv=mixed_bcgs_cgs_solve(x,b0);
      break;

   case AMG_LCD:
      rv=lcd_solve(x,b0);
      break;

    default:
      AMG_Print("solver not implemented\n");
      exit(4711);
  }

  if (global_sc->verbose>=1)
  {
    if(iteration_cnt>0)
    {
      sprintf(buf,"convergence rate (from first iteration) : %g\n",
              (double)pow(end_residual/start_residual,1.0/((double)(iteration_cnt))));
      AMG_Print(buf);
    }
    if (residuals[(residual_cnt)%AMG_CONV_RATE_BACK]>0.0)
    {
      sprintf(buf,"convergence rate (last 10 iterations)   : %g\n",
              (double)pow(end_residual/residuals[(residual_cnt)%AMG_CONV_RATE_BACK],((double)(1.0/AMG_CONV_RATE_BACK))));
      AMG_Print(buf);
    }
  }
  if (sc->row_equilibration)
    DisApplyRowEquilibration(sc,A,b0);
  AMG_ClearSolverOverhead(global_sc,global_cc);
  finish = clock();
  if (finish>=start)
    elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
  else
    elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
  start = clock();

  if (sc->verbose>0)
  {
    sprintf(buf,"SOLUTION TIME %g sec\n",elapsed_time);
    AMG_Print(buf);
  }

  return(rv);                                
}

int AMG_Build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b)
{
  int rv,i,j,memory[2];
//   struct mallinfo MALLINFO;

  /*sprintf(buf,"elapsed_time: %g sec\n", elapsed_time);
  AMG_Print(buf);      
  */

//   MALLINFO = mallinfo();
//   memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
  elapsed_time=0.0;
  start=clock();  
  /* store context */
  global_cc=cc; global_sc=sc;
  
  B[0] = B0[0];
  B[1] = B0[1];
  B[2] = B0[2];
  B[3] = B0[3];
  B[4] = B0[4];
  B[5] = B0[5];
  B[6] = B0[6];
 
  switch (sc->solver)
  {
    case AMG_CG:
      rv=solver_build(sc,cc,A,AMG_CG);
      break;
      
    case AMG_BCGS:
      rv=solver_build(sc,cc,A,AMG_BCGS);
      break;
      
    case AMG_LS:
      rv=solver_build(sc,cc,A,AMG_LS);
      break;
      
    case AMG_GMRES_LEFT:
      rv=solver_build(sc,cc,A,AMG_GMRES_LEFT);
      break;
    case AMG_GMRES_RIGHT:
      rv=solver_build(sc,cc,A,AMG_GMRES_LEFT);
      break;
    case AMG_GMRES_FLEX:
      rv=solver_build(sc,cc,A,AMG_GMRES_FLEX);
      break;
    case AMG_EXACT:
      rv=solver_build(sc,cc,A,AMG_EXACT);
      break;
    case AMG_MIXED_BCGS_CGS:
      rv=solver_build(sc,cc,A,AMG_MIXED_BCGS_CGS);
      break;
    case AMG_LCD:
      rv=solver_build(sc,cc,A,AMG_LCD);
      break;
    default:
      AMG_Print("solver not implemented\n");
      exit(4711);
  }
  if (sc->row_equilibration)
    EquilibrateRHS(b,row_equilibration);

  finish = clock();
  if (finish>=start)
    elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
  else
    elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
  start = clock();

  if (sc->verbose>0)
  {
//     MALLINFO = mallinfo();
//     memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
//     sprintf(buf,"PREPARATION TIME %g sec, MEMORY ALLOCATED %7.3g MB\n",
//             elapsed_time,(memory[1]-memory[0])/(1048576.0));
//     AMG_Print(buf);      
  }

  if (sc->verbose>1)
    ReportSolverParameters(sc);

  return(rv);                                
}

int AMG_Solve (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b)
{
  int rv,i,j,memory[2];
//   struct mallinfo MALLINFO;

  elapsed_time=0.0;
  start=clock();  

  residual_cnt=0;  
  for(i=0;i<10;i++) residuals[i]=0.0;
  switch (global_sc->solver)
  {
    case AMG_CG:
      rv=cg_solve(x,b);
      break;
      
    case AMG_BCGS:
      rv=bcgs_solve(x,b);
      break;
      
    case AMG_LS:
      rv=ls_solve(x,b);
      break;
      
    case AMG_GMRES_LEFT:
      rv=gmres_left_solve(x,b);
      break;
      
    case AMG_GMRES_RIGHT:
      rv=gmres_right_solve(x,b);
      break;

    case AMG_GMRES_FLEX:
      rv=gmres_flexible_solve(x,b);
      break;

    case AMG_EXACT:
      rv=exact_solve(x,b);
      break;
 
   case AMG_MIXED_BCGS_CGS:
      rv=mixed_bcgs_cgs_solve(x,b);
      break;

   case AMG_LCD:
      rv=lcd_solve(x,b);
      break;


    default:
      AMG_Print("solver not implemented\n");
      exit(4711);
  }

  if (global_sc->verbose>=1)
  {
    if(iteration_cnt>0)
    {
      sprintf(buf,"convergence rate (from first iteration) : %g\n",
              (double)pow(end_residual/start_residual,1.0/((double)(iteration_cnt))));
      AMG_Print(buf);
    }
    if (residuals[(residual_cnt)%AMG_CONV_RATE_BACK]>0.0)
    {
      sprintf(buf,"convergence rate (last 10 iterations)   : %g\n",
              (double)pow(end_residual/residuals[(residual_cnt)%AMG_CONV_RATE_BACK],((double)(1.0/AMG_CONV_RATE_BACK))));
      AMG_Print(buf);
    }
  }

  finish = clock();
  if (finish>=start)
    elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
  else
    elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
  start = clock();

  if (sc->verbose>0)
  {
    sprintf(buf,"SOLUTION TIME %g sec\n",elapsed_time);
    AMG_Print(buf);
  }

  return(rv);                                
}


int AMG_Delete (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b)
{
  if (sc->row_equilibration)
    DisApplyRowEquilibration(sc,A,b);
  AMG_ClearSolverOverhead(global_sc,global_cc);

  return(AMG_OK);                                
}
