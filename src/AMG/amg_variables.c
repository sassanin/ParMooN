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
   
#include <time.h>

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

IterProcPtr smoother,coarse_smoother;
IterProcPtr preconditioner,schur_preconditioner,preconditioner_trans;

MultProcPtr dmatmul,dmatminus,A_dmatmul,A_dmatminus,B_dmatmul,B_Trans_dmatmul,
  MGC_dmatminus, dmattransmul;

RestrictProcPtr restriction;
InterpolationProcPtr interpolation;

double start_residual,end_residual,residuals[AMG_CONV_RATE_BACK];
int residual_cnt,iteration_cnt;
double elapsed_time,TIME_WRAP=4295.0;
clock_t start,finish,start1,finish1;

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

int coarse_grid    = 0;          /* multigrid is on coarse grid */
char buf[128];

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/* global data for solvers */
AMG_MATRIX *A[AMG_MAX_LEVELS];
AMG_MATRIX *B[AMG_MAX_LEVELS];
AMG_MATRIX *B_Diri[AMG_MAX_LEVELS];
AMG_MATRIX *schur_matrix[AMG_MAX_LEVELS];
AMG_MATRIX *M[AMG_MAX_LEVELS];
AMG_GRAPH  *G[AMG_MAX_LEVELS];
AMG_GRAPH  *G_schur[AMG_MAX_LEVELS];
AMG_VECTOR *x[AMG_MAX_LEVELS];
AMG_VECTOR *b[AMG_MAX_LEVELS];
AMG_VECTOR *d[AMG_MAX_LEVELS];
AMG_VECTOR *z[AMG_MAX_LEVELS];
AMG_VECTOR *r[AMG_MAX_LEVELS];
AMG_VECTOR *q;
AMG_VECTOR *p[AMG_MAX_LEVELS];
AMG_VECTOR *w;
AMG_VECTOR *old_def[AMG_MAX_LEVELS];
AMG_VECTOR *A_times_update[AMG_MAX_LEVELS];
AMG_VECTOR *s,*cosi,*sn;
AMG_VECTOR *H[AMG_MAX_GMRES_RESTART+1];
AMG_VECTOR *v[AMG_MAX_GMRES_RESTART+1];
AMG_VECTOR *zv[AMG_MAX_GMRES_RESTART+1];
AMG_VECTOR *old_def[AMG_MAX_LEVELS];        /* array of vectors for step length control*/
AMG_VECTOR *A_times_update[AMG_MAX_LEVELS]; /* array of vectors for step length control*/
AMG_VECTOR *row_equilibration;
int depth;
int mgc_recursion[AMG_MAX_LEVELS];
AMG_CoarsenContext *global_cc;
AMG_SolverContext *global_sc;

AMG_VECTOR *s_schur,*cosi_schur,*sn_schur;
AMG_VECTOR *H_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
AMG_VECTOR *v_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
AMG_VECTOR *schur_velo[3][AMG_MAX_LEVELS], *schur_press[5];

AMG_VECTOR *u,*v1,*w,*w1,*t[AMG_MAX_LEVELS],*tilde_r;

AMG_VECTOR *velo_prolong[AMG_MAX_LEVELS],*pres_prolong[AMG_MAX_LEVELS];
AMG_VECTOR *velo_result,*pres_result;

AMG_VECTOR *velo_rhs,*velo_help,*pres_help[3];
AMG_VECTOR *d_bcgs_schur[AMG_MAX_LEVELS];
AMG_VECTOR *z_bcgs_schur[AMG_MAX_LEVELS];
AMG_VECTOR *r_bcgs_schur[AMG_MAX_LEVELS];
AMG_VECTOR *p_bcgs_schur[AMG_MAX_LEVELS];
AMG_VECTOR *w_bcgs_schur;

int *ratr, *jatr, *postr;
