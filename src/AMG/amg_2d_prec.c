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
/*                                                                            */
/* File:          amg_2d_prec.c                                                   */
/*                                                                            */
/* Purpose:   preconditioners for 2D saddle point problems                   */
/*                                                                            */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                            */
/* History:   1998/02/19 start using this library for MooN_MD                    */
/*                                                                            */
/* Remarks:   1998/05/27 Schur complement iteration                                */
/*                                                                            */
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

#include <amg_header.h>
#include <amg_low.h>
#include <amg_sp.h>
#include <amg_blas.h>
#include <amg_iter.h>
#include <amg_coarsen.h>
#include <amg_solve_main.h>
#include <amg_1d_prec.h>
#include <amg_2d_prec.h>


void AMG_ApplyPlaneRotation(double *dx, double *dy, double cs, double sn);
void AMG_GeneratePlaneRotation(double dx, double dy, double *cs, double *sn);
void AMG_Update(AMG_VECTOR *x, int Len_x, int k, AMG_VECTOR **h, AMG_VECTOR *s, AMG_VECTOR **v);

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/*    compile time constants defining static data size (i.e. arrays)            */
/*    other constants                                                            */
/*    macros                                                                    */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                            */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* definition of exported global variables                                    */
/*                                                                            */
/****************************************************************************/

extern IterProcPtr smoother,coarse_smoother;
extern IterProcPtr preconditioner,schur_preconditioner;

extern MultProcPtr dmatmul,dmatminus,A_dmatmul,B_dmatmul,B_Trans_dmatmul,A_dmatminus;

extern double start_residual,end_residual,residuals[AMG_CONV_RATE_BACK];
extern int residual_cnt,iteration_cnt;
extern double elapsed_time,TIME_WRAP;
extern clock_t start,finish,start1,finish1;

extern int coarse_grid;          /* multigrid is on coarse grid */

extern AMG_VECTOR *old_def[AMG_MAX_LEVELS];        /* array of vectors for step length control*/
extern AMG_VECTOR *A_times_update[AMG_MAX_LEVELS]; /* array of vectors for step length control*/
extern AMG_MATRIX *A[AMG_MAX_LEVELS],*B[AMG_MAX_LEVELS];
extern AMG_MATRIX *schur_matrix[AMG_MAX_LEVELS];
extern AMG_MATRIX *B_Diri[AMG_MAX_LEVELS];

extern AMG_VECTOR *s_schur,*cosi_schur,*sn_schur;
extern AMG_VECTOR *H_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
extern AMG_VECTOR *v_schur[AMG_MAX_SCHUR_GMRES_RESTART+1];
extern AMG_VECTOR *schur_velo[3][AMG_MAX_LEVELS], *schur_press[5];
extern AMG_VECTOR *pres_prolong[AMG_MAX_LEVELS];

extern AMG_GRAPH  *G_schur[AMG_MAX_LEVELS];
extern AMG_SolverContext *global_sc;

extern AMG_VECTOR *velo_rhs,*velo_help,*pres_help[3];
extern AMG_VECTOR *d_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *z_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *r_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *p_bcgs_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *w_bcgs_schur;


extern  char buf[128];

/****************************************************************************/
/*                                                                            */
/* definition of variables global to this source file only (static!)            */
/*                                                                            */
/****************************************************************************/
/****************************************************************************/
/*                                                                            */
/* forward declarations of functions used before they are defined            */
/*                                                                            */
/****************************************************************************/

int SymmetrizeSaddleProblem (AMG_SolverContext *sc,AMG_MATRIX *A,
                             AMG_MATRIX *B[AMG_MAX_LEVELS],AMG_VECTOR *b_)
{
  int *ra,*ja,*ara,*aja,i,k,m_b0,n_a,start,end,current,next_start,columns;
  double *a,*b,eps=1e-16;
  
  if (sc->verbose>1)
    AMG_Print("symmetrizing problem\n");

  switch(sc->system_type)
    {
    case SCALAR1 : 
      return(0);
      break;
    case SADDLE_1_TYPE_1 :
      m_b0 = B[0]->m;
      n_a = A->n;
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      ara = AMG_MATRIX_RA(A);
      aja = AMG_MATRIX_JA(A);
      b = b_->x;
      for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                b[i+n_a] -= b[ja[k]]*a[k]; /* update pressure rhs with Dirichlet velo */
                a[k] = 0.0;              /* set matrix entry to 0 */
              }
        }
      break;
    case SADDLE_2_TYPE_1 :     
    case SADDLE_2_TYPE_3 :     
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 :     
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 :     
      m_b0 = B[0]->m;
      n_a = A->n;
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      ara = AMG_MATRIX_RA(A);
      aja = AMG_MATRIX_JA(A);
      b = b_->x;
      for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                b[i+2*n_a] -= b[ja[k]]*a[k]; /* update pressure rhs with Dirichlet velo */
                a[k] = 0.0;              /* set matrix entry to 0 */
              }
        }
      
      a = AMG_MATRIX_A(B[1]);
      ra = AMG_MATRIX_RA(B[1]);
      ja = AMG_MATRIX_JA(B[1]);
      
      for (i=0; i<m_b0; i++)              /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                b[i+2*n_a] -= b[ja[k]+n_a]*a[k]; /* update pressure rhs with Dirichlet velo */
                a[k] = 0.0;              /* set matrix entry to 0 */
              }
        }
      break;
    }

  sc->condense = 1;
  
  if (!sc->condense) return(0);

  switch(sc->system_type)
    {
    case SADDLE_1_TYPE_1 :
      m_b0 = B[0]->m;
      n_a = A->n;
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      ara = AMG_MATRIX_RA(A);
      aja = AMG_MATRIX_JA(A);
      b = b_->x;
      for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                b[i+n_a] -= b[ja[k]]*a[k]; /* update pressure rhs with Dirichlet velo */
                a[k] = 0.0;              /* set matrix entry to 0 */
              }
        }
      break;
    case SADDLE_2_TYPE_1 :     
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 :     
      m_b0 = B[0]->m;
      current = 0;
      next_start = 0;
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          columns = 0;
          start=next_start; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (fabs(a[k])>eps)          /* nonzero entry */
              {
                a[current]=a[k];  /* copy entry on first free place */
                ja[current] = ja[k];
                current++;
                columns++;
              }
          next_start = ra[i+1];
          ra[i+1]=ra[i]+columns;
        }
      if (sc->verbose>1)
        {
          sprintf(buf,"B[0] condensed from %d nonzeros to %d nonzeros !\n",
                  B[0]->nonzeros,current);
          AMG_Print(buf);      
        }
      B[0]->nonzeros = B[0]->connections = current;
      a = AMG_MATRIX_A(B[1]);
      ra = AMG_MATRIX_RA(B[1]);
      ja = AMG_MATRIX_JA(B[1]);
      current = 0;
      next_start = 0;
       for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          columns = 0;
          start=next_start; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (fabs(a[k])>eps)          /* nonzero entry */
              {
                a[current]=a[k];  /* copy entry on first free place */
                ja[current] = ja[k];
                current++;
                columns++;
              }
          next_start = ra[i+1];
          ra[i+1]=ra[i]+columns;
        }
      if (sc->verbose>1)
        {
          sprintf(buf,"B[1] condensed from %d nonzeros to %d nonzeros !\n",
                  B[1]->nonzeros,current);
          AMG_Print(buf);      
        }
      B[1]->nonzeros =  B[1]->connections = current;
     
      break;
    }
  
  return(0);
}

int Build_B_Block_for_Dirichlet(AMG_SolverContext *sc,AMG_MATRIX *A, 
                                AMG_MATRIX *B[AMG_MAX_LEVELS])
{
  int *ra,*ja,*ara,*aja,i,k,m_b0,n_a,start,end,current,next_start,columns;
  double *a;
  
  if (sc->verbose>1)
    AMG_Print("build B_Diri\n");

  switch(sc->system_type)
    {
    case SCALAR1 : 
      return(0);
      break;
    case SADDLE_1_TYPE_1 :
      m_b0 = B[0]->m;
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      ara = AMG_MATRIX_RA(A);
      aja = AMG_MATRIX_JA(A);
      B_Diri[0] = AMG_NewMatrix(B[0]->n,1,B[0]->nonzeros,1,1,1,1,"B_Diri");
      B_Diri[0]->m = m_b0;
      for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                B_Diri[0]->a[k] = 0;     /* set matrix entry to 0 */
              }
            else 
              B_Diri[0]->a[k] = a[k];
        }
      break;
    case SADDLE_2_TYPE_1 :     
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 :     
      m_b0 = B[0]->m;
      B_Diri[0]=AMG_Malloc(sizeof(AMG_MATRIX));
      B_Diri[0]->m = B[0]->m;
      B_Diri[0]->n = B[0]->n;
      B_Diri[0]->b = B[0]->b;
      B_Diri[0]->bb = B[0]->bb;
      B_Diri[0]->system_as_scalar = B[0]->system_as_scalar;
      B_Diri[0]->blocks_in_diag = B[0]->blocks_in_diag;
      B_Diri[0]->bandwidth = B[0]->bandwidth;
      B_Diri[0]->nonzeros = B[0]->nonzeros;
      B_Diri[0]->connections = B[0]->connections;
      B_Diri[0]->level = B[0]->level;
      B_Diri[0]->ra = B[0]->ra;
      B_Diri[0]->ja = B[0]->ja;
      B_Diri[0]->a = AMG_Malloc(B_Diri[0]->nonzeros*sizeof(double));
      B[2]=B_Diri[0];
      a = AMG_MATRIX_A(B[0]);
      ra = AMG_MATRIX_RA(B[0]);
      ja = AMG_MATRIX_JA(B[0]);
      ara = AMG_MATRIX_RA(A);
      aja = AMG_MATRIX_JA(A);
     for (i=0; i<m_b0; i++)     /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                 B_Diri[0]->a[k] = 0.0;              /* set matrix entry to 0 */
              }
            else
              B_Diri[0]->a[k] = a[k];
        }
      B_Diri[1]=AMG_Malloc(sizeof(AMG_MATRIX));
      B_Diri[1]->m = B[1]->m;
      B_Diri[1]->n = B[1]->n;
      B_Diri[1]->b = B[1]->b;
      B_Diri[1]->bb = B[1]->bb;
      B_Diri[1]->system_as_scalar = B[1]->system_as_scalar;
      B_Diri[1]->blocks_in_diag = B[1]->blocks_in_diag;
      B_Diri[1]->bandwidth = B[1]->bandwidth;
      B_Diri[1]->nonzeros = B[1]->nonzeros;
      B_Diri[1]->connections = B[1]->connections;
      B_Diri[1]->level = B[1]->level;
      B_Diri[1]->ra = B[1]->ra;
      B_Diri[1]->ja = B[1]->ja;
      B_Diri[1]->a = AMG_Malloc(B_Diri[1]->nonzeros*sizeof(double));
      B[3]=B_Diri[1];
      a = AMG_MATRIX_A(B[1]);
      ra = AMG_MATRIX_RA(B[1]);
      ja = AMG_MATRIX_JA(B[1]);
      
      for (i=0; i<m_b0; i++)              /* columnwise */ 
        {
          start=ra[i]; 
          end = ra[i+1]-1;
          for (k=start; k<=end; k++)
            if (aja[ara[ja[k]]] == 1)    /* row ja[k] in A has only one entry -> Diriclet */
              {
                B_Diri[1]->a[k] = 0.0;              /* set matrix entry to 0 */
              }
            else
              B_Diri[1]->a[k] = a[k];
        }
      break;
    }

  return(0);
}

/****************************************************************************/
/*                                                                          */
/* schur_complement_fixed                                                   */
/*                                                                          */
/* preconditioner for saddle point problems via the Schur complement        */
/* approach, schur complement equation is solved with a simple fixed point  */
/* iteration                                                                */
/*                                                                          */
/* schur_press[0] : current pressure iterate p                              */
/* schur_press[1] : B^TA^-1Bp and  B^TA^-1Bp+s-BA^-1r                       */
/* schur_press[2] : s-BA^-1r                                                */
/* schur_press[3] : B^TA^-1B*update in step length control                  */
/* schur_velo[0][0] : x=A^-1b                                               */
/* schur_velo[1][0] : rhs b of Ax=b                                         */
/* schur_velo[2][0] : defect b-Ax                                           */
/*                                                                          */
/****************************************************************************/
int prepare_schur_complement_fixed(int *A_length,int *B_length, int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      schur_velo[k][i]= AMG_NewVector(A_length[i],AMG_MATRIX_B(A[i]),"schur_velo"); 
  for (k=0;k<4;k++)
    {
      schur_press[k]= AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"schur_press"); 
    }
  return(0);
}

int clear_schur_complement_fixed(int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      {
        free(schur_velo[k][i]->x);
        free(schur_velo[k][i]);
      }
  for (k=0;k<4;k++)
    {     
      free(schur_press[k]->x);
      free(schur_press[k]);
    }
  return(0);
}

int schur_complement_fixed (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha,schur_alpha,numerator,nominator,eps=1e-12,alpha_eps=0.05;
  int m_b,n_b,iter,max_iter,i,schur_step_length_control,j;
  int  schur_iteration_maxit,schur_inv_of_A_iterations;
  
  schur_alpha = sc->schur_iteration_damp;          /* set parameters */   
  max_iter=sc->schur_iteration_maxit;
  schur_step_length_control=sc-> schur_step_length_control; 
  schur_inv_of_A_iterations=sc->schur_inv_of_A_maxit;
  
  switch(sc->system_type)
    {
    case SADDLE_1 :
      m_b = AMG_MATRIX_M(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_1_TYPE_1 :
      m_b = AMG_MATRIX_N(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_1 : 
    case SADDLE_2_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2 : 
    case SADDLE_2_TYPE_4 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2_MORTAR  : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0])+2*AMG_MATRIX_M(B[4]);/* # columns of B */
       break;
    case SADDLE_3_TYPE_1 : 
    case SADDLE_3_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]) + AMG_MATRIX_N(B[2]); 
                                                    /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    default :
      AMG_Print("system type not defined in schur_complement_gmres\n");
      exit(4711);
      break;
    }
  for (i=0;i<n_b;i++)  schur_press[2]->x[i] = -d[0]->x[i+m_b];  /* -s: rhs in pressure space */
  for (i=0;i<m_b;i++)  
    schur_velo[2][0]->x[i] = schur_velo[1][0]->x[i] = d[0]->x[i];/* r: rhs and defect in velocity space */
   AMG_dset(schur_velo[0][0],0.0);                              /* initial iterate */        
                                                               /* apply approx of A^{-1} */
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]);   
  for (j=1;j<schur_inv_of_A_iterations;j++)                    /* more than one iteration */
    {                                                          /* in approx inverse of A */
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);            /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);   /* compute new defect */
                                                               /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    }  
  B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
  AMG_daxpy(schur_press[2],1,schur_press[0]);                  /* compute rhs of schur equ. -s+B^TA^-1r */      
  AMG_dset(schur_press[0],0.0);                              /* initial iterate of fixed point iteration */        
  for (iter=0;iter<max_iter;iter++)                           /* start fixed point iteration */
    {
      if (iter>0)                                             /* compute rhs for next iteration */ 
        {
          B_dmatmul(schur_velo[2][0],A[0],B,schur_press[0]);  /* multiply last iterate with B */
          AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);       /* copy result to defect vector since */
          AMG_dset(schur_velo[0][0],0.0);                     /* initial iterate is set to zero */
                                                              /* apply approx of A^{-1} */
          schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
          for (j=1;j<schur_inv_of_A_iterations;j++)            /* more than one iteration */
            {                                                  /* in approx inverse of A */
              AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);          /* copy rhs to defect */ 
              A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]); /* compute new defect */
                                                               /* apply approx of A^{-1} */
              schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
            } 
          B_Trans_dmatmul(schur_press[1],A[0],B,schur_velo[0][0]);   /*  multiply with B^T */
        }
      else                                                     /* first iteration */
        AMG_dset(schur_press[1],0.0);
                                                               /* end if */
      AMG_daxpy(schur_press[1],-1.0,schur_press[2]);             /* compute B^TA^-1Bp^m - B^TA^-1r+s */
 
      if (sc->verbose>1)
        {
          sprintf(buf,"schur %d %g\n",iter,sqrt(AMG_ddot(schur_press[1],schur_press[1])));
          AMG_Print(buf);
        }

      if (sqrt(AMG_ddot(schur_press[1],schur_press[1]))<eps)  /* stopping criterion */
        break;

      if (schur_step_length_control)                         /* apply step length control */
        {
          B_dmatmul(schur_velo[2][0],A[0],B,schur_press[1]);          
          AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);       /* copy result to defect vector since */
          AMG_dset(schur_velo[0][0],0.0);                     /* initial iterate is set to zero */
                                                              /* apply approx of A^{-1} */
          schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
          for (j=1;j<schur_inv_of_A_iterations;j++)            /* more than one iteration */
            {
              AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);          /* copy rhs to defect */ 
              A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]); /* compute new defect */
                                                               /* apply approx of A^{-1} */
              schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
            } 
          B_Trans_dmatmul(schur_press[3],A[0],B,schur_velo[0][0]);   /*  multiply with B^T */
          
          numerator =  AMG_ddot(schur_press[1],schur_press[3]);
          nominator =  AMG_ddot(schur_press[3],schur_press[3]);
          
          schur_alpha= numerator/nominator;
          if (sc->verbose>2)
            {
              sprintf(buf,"schur alpha %g %g %g\n",schur_alpha,numerator,nominator);
              AMG_Print(buf);
            }
        }
      AMG_daxpy(schur_press[0],-schur_alpha,schur_press[1]);           /* new iterate */
    }                                                         /* end fixed point iteration */ 
  
  for (i=0;i<n_b;i++)  d[0]->x[i+m_b] = schur_press[0]->x[i]; /* copy last iterate to solution vector */
  
  B_dmatmul(schur_velo[0][0],A[0],B,schur_press[0]);          /* multiply last iterate with B */  
  for (i=0;i<m_b;i++)  schur_velo[2][0]->x[i] = d[0]->x[i];
  AMG_daxpy(schur_velo[2][0],-1,schur_velo[0][0]);
  AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);            
  AMG_dset(schur_velo[0][0],0.0); 
                                                              /* apply approx of A^{-1} */
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
  for (j=1;j<schur_inv_of_A_iterations;j++)                   /* more than one iteration */
    {
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);           /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);  /* compute new defect */
                                                              /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    } 
  
  for (i=0;i<m_b;i++)  d[0]->x[i] =schur_velo[0][0]->x[i] ;  /* copy velo result to solution vector */
  
                                                  /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    {                                            /* old defect in b[k] */
      dmatmul(A_times_update[k],A[k],B,d[k]);    /* current update in d[k] */
      
      numerator =  AMG_ddot(A_times_update[k],b[k]);
      nominator =  AMG_ddot(A_times_update[k],A_times_update[k]);
      if (nominator > 0)
        alpha = numerator/nominator;
      else
        {
          alpha = 0.05;
          sprintf(buf,"MESSAGE : Step length control failed on level %d, set alpha = %g !!\n",k,alpha);
          AMG_Print(buf);
        }
      if (fabs(alpha)<alpha_eps)
        {
          /* sprintf(buf,"MESSAGE : Step length factor %g, too small, set alpha = %g\n",
             alpha,alpha_eps);
             AMG_Print(buf);*/
          alpha = alpha_eps;
        }
          if (sc->verbose>2)
            {
              sprintf(buf,"alpha %g\n",alpha);
              AMG_Print(buf);
            }
    }
  else 
    alpha =sc->omega[0];                                 /* no step length control */

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */
  
  return(AMG_OK);
}

/****************************************************************************/
/*                                                                            */
/* schur_complement_cg                                                      */
/*                                                                            */
/* preconditioner for saddle point problems via the Schur complement             */
/* approach, schur complement equation is solved with a cg-method           */
/*                                                                            */
/* schur_press[0] : x (in cg as solver)                                           */
/* schur_press[1] : z[0] (in cg as solver)                                         */
/* schur_press[2] : r[0] (in cg as solver)                                    */
/* schur_press[3] : d[0] (in cg as solver)                                     */
/* schur_press[4] : q    (in cg as solver)                                     */
/* schur_velo[0][0] : x=A^-1b                                                    */
/* schur_velo[1][0] : rhs b of Ax=b                                             */
/* schur_velo[2][0] : defect b-Ax                                            */
/*                                                                            */
/****************************************************************************/
int prepare_schur_complement_cg(int *A_length,int *B_length, int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      schur_velo[k][i]= AMG_NewVector(A_length[i],AMG_MATRIX_B(A[i]),"schur_velo"); 
  for (k=0;k<5;k++)
    schur_press[k]= AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"schur_press");   
  return(0);
}

int clear_schur_complement_cg(int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      {
        free(schur_velo[k][i]->x);
        free(schur_velo[k][i]);
      }
  for (k=0;k<5;k++)
    {
      free(schur_press[k]->x);
      free(schur_press[k]);
    }
  return(0);
}

int schur_complement_cg (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha,schur_alpha,numerator,nominator,eps=1e-12,alpha_eps=0.05;
  int m_b,n_b,iter,max_iter,i,schur_step_length_control,j;
  int  schur_iteration_maxit,schur_inv_of_A_iterations;
  double dnorm,dnorm0,dnormlast,rho,rho_last;

  schur_alpha = sc->schur_iteration_damp;          /* set parameters */   
  max_iter=sc->schur_iteration_maxit;
  schur_step_length_control=sc-> schur_step_length_control; 
  schur_inv_of_A_iterations=sc->schur_inv_of_A_maxit;
  switch(sc->system_type)
    {
    case SADDLE_1 :
      m_b = AMG_MATRIX_M(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_1_TYPE_1 :
      m_b = AMG_MATRIX_N(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_1 : 
    case SADDLE_2_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2 : 
    case SADDLE_2_TYPE_4 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2_MORTAR  : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0])+2*AMG_MATRIX_M(B[4]);/* # columns of B */
       break;
    case SADDLE_3_TYPE_1 : 
    case SADDLE_3_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1])+ AMG_MATRIX_N(B[2]); 
                                                    /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
    default :
      AMG_Print("system type not defined in schur_complement_gmres\n");
      exit(4711);
      break;
    }
  for (i=0;i<n_b;i++)  schur_press[2]->x[i] = -d[0]->x[i+m_b];  /* -s: rhs in pressure space */
  for (i=0;i<m_b;i++)  
    schur_velo[2][0]->x[i] = schur_velo[1][0]->x[i] = d[0]->x[i];/* r: rhs and defect in velocity space */
  AMG_dset(schur_velo[0][0],0.0);                              /* initial iterate */        
                                                               /* apply approx of A^{-1} */                
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
  for (j=1;j<schur_inv_of_A_iterations;j++)                    /* more than one iteration */
    {                                                          /* in approx inverse of A */
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);            /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);   /* compute new defect */
                                                               /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    }  
  B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
  AMG_daxpy(schur_press[2],1,schur_press[0]);                  /* compute rhs of schur equ. -s+B^TA^-1r */      
  AMG_dset(schur_press[0],0.0);                                 /* initial iterate of cg iteration */
  rho_last = 1.0;                                                /* prepare cg */
  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(schur_press[2],schur_press[2]));
  if (sc->verbose>1)
    {
      sprintf(buf,"Entering Schur_Complement_CG: %4d %12.4lE \n",0,dnorm);
      AMG_Print(buf);
    } 
  for (iter=0;iter<max_iter;iter++)                           /* start cg iteration */
    {
      AMG_dcopy(schur_press[1],schur_press[2]);
      rho = AMG_ddot(schur_press[2],schur_press[1]);
      if (i>0) AMG_daxpy(schur_press[1],rho/rho_last,schur_press[4]);
      rho_last=rho;
      AMG_dcopy(schur_press[4],schur_press[1]);
      B_dmatmul(schur_velo[2][0],A[0],B,schur_press[4]);  /* multiply last iterate with B */
      AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);       /* copy result to defect vector since */
      AMG_dset(schur_velo[0][0],0.0);                          /* initial iterate */        
                                                               /* apply approx of A^{-1} */                
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
      for (j=1;j<schur_inv_of_A_iterations;j++)                    /* more than one iteration */
        {                                                          /* in approx inverse of A */
          AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);            /* copy rhs to defect */ 
          A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);   /* compute new defect */
                                                                        /* apply approx of A^{-1} */
          schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
        }  
      B_Trans_dmatmul(schur_press[3],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
      alpha=rho/AMG_ddot(schur_press[4],schur_press[3]);
      AMG_daxpy(schur_press[0],alpha,schur_press[4]);
      AMG_daxpy(schur_press[2],-alpha,schur_press[3]);
      dnorm=sqrt(AMG_ddot(schur_press[2],schur_press[2]));
      if (sc->verbose>0)
        {
          /*sprintf(buf,"Iteration0 %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);*/
          sprintf(buf,"CG Iteration %4d %g %g %12.4lg\n",iter+1,dnorm,sc->res_norm_min,dnorm/dnormlast);
          AMG_Print(buf);
        }
      dnormlast=dnorm;
      /*if (dnorm<dnorm0*sc->red_factor) break;
        if (dnorm<sc->res_norm_min) break;*/
    }                                                         /* end cg iteration */ 
  
  for (i=0;i<n_b;i++)  d[0]->x[i+m_b] = schur_press[0]->x[i]; /* copy last iterate to solution vector */
  
  B_dmatmul(schur_velo[0][0],A[0],B,schur_press[0]);          /* multiply last iterate with B */  
  for (i=0;i<m_b;i++)  schur_velo[2][0]->x[i] = d[0]->x[i];
  AMG_daxpy(schur_velo[2][0],-1,schur_velo[0][0]);
  AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);            
  AMG_dset(schur_velo[0][0],0.0); 
                                                              /* apply approx of A^{-1} */
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
  for (j=1;j<schur_inv_of_A_iterations;j++)                   /* more than one iteration */
    {
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);           /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);  /* compute new defect */
                                                              /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    } 
  
  for (i=0;i<m_b;i++)  d[0]->x[i] =schur_velo[0][0]->x[i] ;  /* copy velo result to solution vector */
  
                                                  /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    {                                            /* old defect in b[k] */
      dmatmul(A_times_update[k],A[k],B,d[k]);    /* current update in d[k] */
      
      numerator =  AMG_ddot(A_times_update[k],b[k]);
      nominator =  AMG_ddot(A_times_update[k],A_times_update[k]);
      if (nominator > 0)
        alpha = numerator/nominator;
      else
        {
          alpha = 0.05;
          sprintf(buf,"MESSAGE : Step length control failed on level %d, set alpha = %g !!\n",k,alpha);
          AMG_Print(buf);
        }
      if (fabs(alpha)<alpha_eps)
        {
          /* sprintf(buf,"MESSAGE : Step length factor %g, too small, set alpha = %g\n",
             alpha,alpha_eps);
             AMG_Print(buf);*/
          alpha = alpha_eps;
        }
          if (sc->verbose>2)
            {
              sprintf(buf,"alpha %g\n",alpha);
              AMG_Print(buf);
            }
    }
  else 
    alpha =sc->omega[0];                                 /* no step length control */

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */
  
  return(AMG_OK);
}

/****************************************************************************/
/*                                                                            */
/* schur_complement_gmres                                                   */
/*                                                                            */
/* preconditioner for saddle point problems via the Schur complement             */
/* approach, schur complement equation is solved with a gmres-method        */
/*                                                                            */
/* schur_press[0] : x,d[0] (in gmres solver)                                      */
/* schur_press[1] : z[0],r[0] (in gmres solver)                                    */
/* schur_velo[0][0] : x=A^-1b                                                    */
/* schur_velo[1][0] : rhs b of Ax=b                                             */
/* schur_velo[2][0] : defect b-Ax                                            */
/*                                                                            */
/****************************************************************************/
int prepare_schur_complement_gmres(AMG_SolverContext *sc, int *A_length,int *B_length, int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      schur_velo[k][i]= AMG_NewVector(A_length[i],AMG_MATRIX_B(A[i]),"schur_velo"); 
  for (k=0;k<1;k++)
    schur_press[k]= AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"schur_press");   

  if (sc->schur_iteration_maxit>AMG_MAX_SCHUR_GMRES_RESTART)
    {
      sprintf(buf,"MESSAGE : schur_iteration_maxit %d greater than AMG_MAX_SCHUR_GMRES_RESTART %d!\n",
              sc->schur_iteration_maxit,AMG_MAX_SCHUR_GMRES_RESTART);
      AMG_Print(buf);
      sprintf(buf,"MESSAGE : set schur_iteration_maxit to %d!\n",AMG_MAX_SCHUR_GMRES_RESTART);
      AMG_Print(buf);
      sc->gmres_restart=AMG_MAX_SCHUR_GMRES_RESTART;
    }

  s_schur = AMG_NewVector(sc->schur_iteration_maxit+1,AMG_MATRIX_B(A[0]),"s");   
  cosi_schur = AMG_NewVector(sc->schur_iteration_maxit+1,AMG_MATRIX_B(A[0]),"cs"); 
  sn_schur = AMG_NewVector(sc->schur_iteration_maxit+1,AMG_MATRIX_B(A[0]),"sn");
  for (k=0;k<sc->schur_iteration_maxit+1;k++)    /* allocate matrices depending on schur_iteration_maxit*/
    {
      H_schur[k] = AMG_NewVector(sc->schur_iteration_maxit+1,AMG_MATRIX_B(A[0]),"H");
      AMG_dset(H_schur[k],0.0);
      v_schur[k] = AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"v");
      AMG_dset(v_schur[k],0.0);      
    }
  return(0);
}

int clear_schur_complement_gmres(AMG_SolverContext *sc,int depth)
{
  int i,k;
  for (k=0;k<3;k++)
    for (i=0;i<=depth;i++)
      {
        free(schur_velo[k][i]->x);
        free(schur_velo[k][i]);
      }
  for (k=0;k<1;k++)
    {
      free(schur_press[k]->x);
      free(schur_press[k]);
    }
     
  free(s_schur->x);
  free(s_schur);
  free(cosi_schur->x);
  free(cosi_schur);
  free(sn_schur->x);
  free(sn_schur);
  for (k=0;k<sc->schur_iteration_maxit+1;k++)
    {
      free(H_schur[k]->x);
      free(H_schur[k]);
      free(v_schur[k]->x);
      free(v_schur[k]);
    }
  return(0);
}

int schur_complement_gmres (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha,schur_alpha,numerator,nominator,eps=1e-12,alpha_eps=0.05;
  int m_b,n_b,iter,max_iter,i,schur_step_length_control,j,ite=0,k1;
  int  schur_iteration_maxit,schur_inv_of_A_iterations;
  double resid,beta,residlast;

  schur_alpha = sc->schur_iteration_damp;          /* set parameters */   
  schur_iteration_maxit=sc->schur_iteration_maxit;
  schur_step_length_control=sc-> schur_step_length_control; 
  schur_inv_of_A_iterations=sc->schur_inv_of_A_maxit;
  switch(sc->system_type)
    {
    case SADDLE_1 :
      m_b = AMG_MATRIX_M(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_1_TYPE_1 :
      m_b = AMG_MATRIX_N(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_1 : 
    case SADDLE_2_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2 : 
    case SADDLE_2_TYPE_4 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2_MORTAR  : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0])+2*AMG_MATRIX_M(B[4]);/* # columns of B */
       break;
    case SADDLE_3_TYPE_1 : 
    case SADDLE_3_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1])+ AMG_MATRIX_N(B[2]); 
                                                    /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    default :
      AMG_Print("system type not defined in schur_complement_gmres\n");
      exit(4711);
      break;
    }
  for (i=0;i<n_b;i++)  v_schur[0]->x[i] = -d[0]->x[i+m_b];  /* -s: rhs in pressure space */
  for (i=0;i<m_b;i++)  
    schur_velo[2][0]->x[i] = schur_velo[1][0]->x[i] = d[0]->x[i];/* r: rhs and defect in velocity space */
  AMG_dset(schur_velo[0][0],0.0);                              /* initial iterate */        
                                                               /* apply approx of A^{-1} */                
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
  for (j=1;j<schur_inv_of_A_iterations;j++)                    /* more than one iteration */
    {                                                          /* in approx inverse of A */
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);            /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);   /* compute new defect */
                                                               /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    }  
  B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
  AMG_daxpy(v_schur[0],1,schur_press[0]);                  /* compute rhs of schur equ. -s+B^TA^-1r */      

  resid = beta=sqrt(AMG_ddot(v_schur[0],v_schur[0]));           /* norm of residual */
  AMG_dscale(v_schur[0],1.0/beta);
  AMG_dset(s_schur,0.0);
  s_schur->x[0]=beta;
  for (i=0;i<schur_iteration_maxit;i++)
    {
      ite++;
      B_dmatmul(schur_velo[2][0],A[0],B,v_schur[i]);  /* multiply last iterate with B */
      AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);       /* copy result to defect vector since */
      AMG_dset(schur_velo[0][0],0.0);                          /* initial iterate */        
                                                               /* apply approx of A^{-1} */                
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
      for (j=1;j<schur_inv_of_A_iterations;j++)                    /* more than one iteration */
        {                                                          /* in approx inverse of A */
          AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);            /* copy rhs to defect */ 
          A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);   /* compute new defect */
                                                                        /* apply approx of A^{-1} */
          schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
        }  
      B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
      /*AMG_dcopy(schur_press[0],schur_press[1]);*/
            
    for(k1=0; k1<=i; k1++)
        {
          H_schur[k1]->x[i]=AMG_ddot(schur_press[0], v_schur[k1]);
          AMG_daxpy(schur_press[0],-H_schur[k1]->x[i], v_schur[k1]);
        }                                          
      
     H_schur[i+1]->x[i]=sqrt(AMG_ddot(schur_press[0],schur_press[0]));
     AMG_dcopy(v_schur[i+1],schur_press[0]);
     AMG_dscale(v_schur[i+1], 1.0/(H_schur[i+1]->x[i]));
      
     for(k1=0;k1<i;k1++)
       AMG_ApplyPlaneRotation(&H_schur[k1]->x[i], &H_schur[k1+1]->x[i], cosi_schur->x[k1], sn_schur->x[k1]);
     
      AMG_GeneratePlaneRotation(H_schur[i]->x[i], H_schur[i+1]->x[i], &cosi_schur->x[i], &sn_schur->x[i]);
      AMG_ApplyPlaneRotation(&H_schur[i]->x[i], &H_schur[i+1]->x[i], cosi_schur->x[i], sn_schur->x[i]);
      AMG_ApplyPlaneRotation(&s_schur->x[i], &s_schur->x[i+1], cosi_schur->x[i], sn_schur->x[i]);
     
      residlast=resid;
      resid=fabs(s_schur->x[i+1]);
      if (sc->verbose>1)
        {
          sprintf(buf,"Schur GMRES Iteration %4d %4d %12.4lE %12.4lg\n",ite,i,resid,resid/residlast);
          AMG_Print(buf);
        }
     if (resid<1.0e-12) break;
 
    }                                                   /* endfor of gmres iteration */
 
  AMG_dset(schur_press[0],0.0);                         /* transform to obtain solution */
  AMG_Update(schur_press[0],AMG_VECTOR_N(schur_press[0]),ite-1, H_schur, s_schur, v_schur);
      
  if (sc->verbose>1)
    {
      sprintf(buf,"Schur GMRES : iterations %d residual %g \n",ite,resid);
      AMG_Print(buf);
    }
  for (i=0;i<n_b;i++)  d[0]->x[i+m_b] = schur_press[0]->x[i]; /* copy last iterate to solution vector */
  
  B_dmatmul(schur_velo[0][0],A[0],B,schur_press[0]);          /* multiply last iterate with B */  
  for (i=0;i<m_b;i++)  schur_velo[2][0]->x[i] = d[0]->x[i];
  AMG_daxpy(schur_velo[2][0],-1,schur_velo[0][0]);
  AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);            
  AMG_dset(schur_velo[0][0],0.0); 
                                                             /* apply approx of A^{-1} */
  schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
  for (j=1;j<schur_inv_of_A_iterations;j++)                   /* more than one iteration */
    {
      AMG_dcopy(schur_velo[2][0],schur_velo[1][0]);           /* copy rhs to defect */ 
      A_dmatminus(schur_velo[2][0],A[0],B,schur_velo[0][0]);  /* compute new defect */
                                                              /* apply approx of A^{-1} */
      schur_preconditioner(sc,0,depth,A,G,M,B,schur_velo[0],schur_velo[1],schur_velo[2]); 
    } 

  for (i=0;i<m_b;i++)  d[0]->x[i] =schur_velo[0][0]->x[i] ;  /* copy velo result to solution vector */

                                                  /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    {                                            /* old defect in b[k] */
      dmatmul(A_times_update[k],A[k],B,d[k]);    /* current update in d[k] */
      
      numerator =  AMG_ddot(A_times_update[k],b[k]);
      nominator =  AMG_ddot(A_times_update[k],A_times_update[k]);
      if (nominator > 0)
        alpha = numerator/nominator;
      else
        {
          alpha = 0.05;
          sprintf(buf,"MESSAGE : Step length control failed on level %d, set alpha = %g !!\n",k,alpha);
          AMG_Print(buf);
        }
      if (fabs(alpha)<alpha_eps)
        {
          /* sprintf(buf,"MESSAGE : Step length factor %g, too small, set alpha = %g\n",
             alpha,alpha_eps);
             AMG_Print(buf);*/
          alpha = alpha_eps;
        }
          if (sc->verbose>2)
            {
              sprintf(buf,"alpha %g\n",alpha);
              AMG_Print(buf);
            }
    }
  else 
    alpha =sc->omega[0];                                 /* no step length control */

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                            */
/* BiConjugate Gradient Stabilized                                          */
/*                                                                            */
/****************************************************************************/
int prepare_schur_bcgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;

  w_bcgs_schur = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"w");
  z_bcgs_schur[0] = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"z");
  r_bcgs_schur[0] = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"r");
  p_bcgs_schur[0] = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"p");
  d_bcgs_schur[0] = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"d");
  if (sc->schur_inv_of_A!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      z_bcgs_schur[k] = AMG_NewVector(A_length[k],AMG_MATRIX_B(A[k]),"z");
      r_bcgs_schur[k] = AMG_NewVector(A_length[k],AMG_MATRIX_B(A[k]),"r");
      p_bcgs_schur[k] = AMG_NewVector(A_length[k],AMG_MATRIX_B(A[k]),"p");
      d_bcgs_schur[k] = AMG_NewVector(A_length[k],AMG_MATRIX_B(A[k]),"d");
    }
  return(0);
}
int clear_schur_bcgs_solve(AMG_SolverContext *sc,int depth)
{
  int k;

  free(w_bcgs_schur->x);
  free(w_bcgs_schur);
  free(z_bcgs_schur[0]->x);
  free(z_bcgs_schur[0]);
  free(d_bcgs_schur[0]->x);
  free(d_bcgs_schur[0]);
  free(r_bcgs_schur[0]->x);
  free(r_bcgs_schur[0]);
  free(p_bcgs_schur[0]->x);
  free(p_bcgs_schur[0]);
  if (sc->schur_inv_of_A!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(z_bcgs_schur[k]->x);
      free(z_bcgs_schur[k]);
      free(d_bcgs_schur[k]->x);
      free(d_bcgs_schur[k]);
      free(r_bcgs_schur[k]->x);
      free(r_bcgs_schur[k]);
      free(p_bcgs_schur[k]->x);
      free(p_bcgs_schur[k]);
    }
  return(0);
}

int schur_bcgs_solve (AMG_SolverContext *sc,int depth,        AMG_MATRIX *AA[AMG_MAX_LEVELS],
                      AMG_GRAPH *GG[AMG_MAX_LEVELS],AMG_MATRIX *MM[AMG_MAX_LEVELS],
                      AMG_MATRIX *BB[AMG_MAX_LEVELS], 
                      AMG_VECTOR *x,AMG_VECTOR *def,
                      AMG_VECTOR *b)
{
  int i,j;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last,alpha,beta,omega,delta;
  
                                                 /* iterate */
  AMG_dcopy(r_bcgs_schur[0],b);
  A_dmatminus(r_bcgs_schur[0],AA[0],BB,x);                    /* r[0] is initial residual */
  AMG_dcopy(b,r_bcgs_schur[0]);                                 /* b is initial residual, r^tilde in algorithm */

  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r_bcgs_schur[0],r_bcgs_schur[0]));
  start_residual=dnorm;
  residuals[0]=dnorm;
  residual_cnt++;
  
  if (sc->verbose>0)
    {
      sprintf(buf,"Entering SCHUR_BCGS: %4d %12.4lE \n",0,dnorm);
      AMG_Print(buf);
    }
  for (i=0; i<sc->schur_inv_of_A_maxit; i++)
    {
      rho=AMG_ddot(r_bcgs_schur[0],b);                       /* b is r[0]^tilde */
      if (fabs(rho)<1.0E-50)
        {
          AMG_Print("BCGS break down\n");
          iteration_cnt=i;
          end_residual=dnorm;                
          exit(4711);
        }
      
      if (i==0)                                   /* compute p */
        AMG_dcopy(p_bcgs_schur[0],r_bcgs_schur[0]);
      else
        {
          beta=(rho/rho_last)*(alpha/omega);
          AMG_dscale(p_bcgs_schur[0],beta);
          AMG_daxpy(p_bcgs_schur[0],1.0,r_bcgs_schur[0]);
          AMG_daxpy(p_bcgs_schur[0],-beta*omega,w_bcgs_schur);
        }
      rho_last=rho;
      
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dset(z_bcgs_schur[0],0.0);
          AMG_dcopy(d_bcgs_schur[0],p_bcgs_schur[0]);
          for (j=0;j<sc->amg_prec_it;j++)
            {
              schur_preconditioner(sc,0,depth,AA,GG,MM,BB,z_bcgs_schur,p_bcgs_schur,d_bcgs_schur);
              AMG_dcopy(d_bcgs_schur[0],p_bcgs_schur[0]);
              A_dmatminus(d_bcgs_schur[0],AA[0],BB,z_bcgs_schur[0]);
            }
        }
      else
          AMG_dcopy(z_bcgs_schur[0],p_bcgs_schur[0]);        
      A_dmatmul(w_bcgs_schur,AA[0],BB,z_bcgs_schur[0]);
      delta=AMG_ddot(b,w_bcgs_schur);
      if (fabs(delta)<1.0E-50) delta=1.0E-50;
      alpha=rho/delta;
      AMG_daxpy(r_bcgs_schur[0],-alpha,w_bcgs_schur);
      AMG_daxpy(x,alpha,z_bcgs_schur[0]);
      dnorm=sqrt(AMG_ddot(r_bcgs_schur[0],r_bcgs_schur[0]));
      
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dset(z_bcgs_schur[0],0.0);
          AMG_dcopy(d_bcgs_schur[0],r_bcgs_schur[0]);
          for (j=0;j<sc->amg_prec_it;j++)
            {
              schur_preconditioner(sc,0,depth,AA,GG,MM,BB,z_bcgs_schur,r_bcgs_schur,d_bcgs_schur);
              AMG_dcopy(d_bcgs_schur[0],r_bcgs_schur[0]);
              A_dmatminus(d_bcgs_schur[0],AA[0],BB,z_bcgs_schur[0]);
            }
        }
      else
          AMG_dcopy(z_bcgs_schur[0],r_bcgs_schur[0]);
        
      A_dmatmul(d_bcgs_schur[0],AA[0],BB,z_bcgs_schur[0]);
      delta=AMG_ddot(d_bcgs_schur[0],d_bcgs_schur[0]);
      if (fabs(delta)<1.0E-50) delta=1.0E-50;
      omega=AMG_ddot(d_bcgs_schur[0],r_bcgs_schur[0])/delta;
      AMG_daxpy(r_bcgs_schur[0],-omega,d_bcgs_schur[0]);
      AMG_daxpy(x,omega,z_bcgs_schur[0]);
      dnorm=sqrt(AMG_ddot(r_bcgs_schur[0],r_bcgs_schur[0]));
      if (sc->verbose>0)
        {
          sprintf(buf,"SCHUR_BCGS Iteration %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
          AMG_Print(buf);
        }
      dnormlast=dnorm;
    }

  if (sc->verbose>=0)
    {
      sprintf(buf,"SCHUR_BCGS : iterations %d residual %g \n",i,dnorm);
      AMG_Print(buf);
    }
  return(0);
}                        

int schur_complement_gmres_bcgs (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha,schur_alpha,numerator,nominator,eps=1e-12,alpha_eps=0.05;
  int m_b,n_b,iter,max_iter,i,schur_step_length_control,j,ite=0,k1;
  int  schur_iteration_maxit,schur_inv_of_A_iterations;
  double resid,beta,residlast;

  printf("schur -1\n");
  schur_alpha = sc->schur_iteration_damp;          /* set parameters */   
  schur_iteration_maxit=sc->schur_iteration_maxit;
  schur_step_length_control=sc-> schur_step_length_control; 
  schur_inv_of_A_iterations=sc->schur_inv_of_A_maxit;
  switch(sc->system_type)
  {
    case SADDLE_1 :
      m_b = AMG_MATRIX_M(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_1_TYPE_1 :
      m_b = AMG_MATRIX_N(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_1 : 
    case SADDLE_2_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2 : 
    case SADDLE_2_TYPE_4 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_2 : 
    case BRAESS_SARAZIN_SADDLE_2_TYPE_4 : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_2_MORTAR  : 
      m_b = AMG_MATRIX_M(B[0])+ AMG_MATRIX_M(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0])+2*AMG_MATRIX_M(B[4]);/* # columns of B */
      break;
    case SADDLE_3_TYPE_1 : 
    case SADDLE_3_TYPE_3 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_1 : 
    case BRAESS_SARAZIN_SADDLE_3_TYPE_3 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1])+ AMG_MATRIX_N(B[2]); 
                                                    /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    default :
      AMG_Print("system type not defined in schur_complement_gmres\n");
      exit(4711);
      break;
   }
  for (i=0;i<n_b;i++)  v_schur[0]->x[i] = -d[0]->x[i+m_b];  /* -s: rhs in pressure space */
  for (i=0;i<m_b;i++)  
    schur_velo[2][0]->x[i] = schur_velo[1][0]->x[i] = d[0]->x[i];/* r: rhs and defect in velocity space */
  AMG_dset(schur_velo[0][0],0.0);                              /* initial iterate */        
  printf("schur 0\n");
                                                               /* apply approx of A^{-1} */                
  schur_bcgs_solve(sc,depth,A,G,M,B,schur_velo[0][0],schur_velo[1][0],schur_velo[2][0]); 
  printf("schur 1\n");
  B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
  AMG_daxpy(v_schur[0],1,schur_press[0]);                  /* compute rhs of schur equ. -s+B^TA^-1r */      

  resid = beta=sqrt(AMG_ddot(v_schur[0],v_schur[0]));           /* norm of residual */
  AMG_dscale(v_schur[0],1.0/beta);
  AMG_dset(s_schur,0.0);
  s_schur->x[0]=beta;
  printf("schur 2\n");
  for (i=0;i<schur_iteration_maxit;i++)
    {
      ite++;
      B_dmatmul(schur_velo[2][0],A[0],B,v_schur[i]);  /* multiply last iterate with B */
      AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);       /* copy result to defect vector since */
      AMG_dset(schur_velo[0][0],0.0);                          /* initial iterate */        
                                                               /* apply approx of A^{-1} */                
      schur_bcgs_solve(sc,depth,A,G,M,B,schur_velo[0][0],schur_velo[1][0],schur_velo[2][0]); 
      B_Trans_dmatmul(schur_press[0],A[0],B,schur_velo[0][0]);      /* multiply with B^T */
      /*AMG_dcopy(schur_press[0],schur_press[1]);*/
            
    for(k1=0; k1<=i; k1++)
        {
          H_schur[k1]->x[i]=AMG_ddot(schur_press[0], v_schur[k1]);
          AMG_daxpy(schur_press[0],-H_schur[k1]->x[i], v_schur[k1]);
        }                                          
      
     H_schur[i+1]->x[i]=sqrt(AMG_ddot(schur_press[0],schur_press[0]));
     AMG_dcopy(v_schur[i+1],schur_press[0]);
     AMG_dscale(v_schur[i+1], 1.0/(H_schur[i+1]->x[i]));
      
     for(k1=0;k1<i;k1++)
       AMG_ApplyPlaneRotation(&H_schur[k1]->x[i], &H_schur[k1+1]->x[i], cosi_schur->x[k1], sn_schur->x[k1]);
     
      AMG_GeneratePlaneRotation(H_schur[i]->x[i], H_schur[i+1]->x[i], &cosi_schur->x[i], &sn_schur->x[i]);
      AMG_ApplyPlaneRotation(&H_schur[i]->x[i], &H_schur[i+1]->x[i], cosi_schur->x[i], sn_schur->x[i]);
      AMG_ApplyPlaneRotation(&s_schur->x[i], &s_schur->x[i+1], cosi_schur->x[i], sn_schur->x[i]);
     
      residlast=resid;
      resid=fabs(s_schur->x[i+1]);
      if (sc->verbose>1)
        {
          sprintf(buf,"Schur GMRES Iteration %4d %4d %12.4lE %12.4lg\n",ite,i,resid,resid/residlast);
          AMG_Print(buf);
        }
     if (resid<1.0e-12) break;
 
    }                                                   /* endfor of gmres iteration */
 
  AMG_dset(schur_press[0],0.0);                         /* transform to obtain solution */
  AMG_Update(schur_press[0],AMG_VECTOR_N(schur_press[0]),ite-1, H_schur, s_schur, v_schur);
      
  if (sc->verbose>1)
    {
      sprintf(buf,"Schur GMRES : iterations %d residual %g \n",ite,resid);
      AMG_Print(buf);
    }
  for (i=0;i<n_b;i++)  d[0]->x[i+m_b] = schur_press[0]->x[i]; /* copy last iterate to solution vector */
  
  B_dmatmul(schur_velo[0][0],A[0],B,schur_press[0]);          /* multiply last iterate with B */  
  for (i=0;i<m_b;i++)  schur_velo[2][0]->x[i] = d[0]->x[i];
  AMG_daxpy(schur_velo[2][0],-1,schur_velo[0][0]);
  AMG_dcopy(schur_velo[1][0],schur_velo[2][0]);            
  AMG_dset(schur_velo[0][0],0.0); 
                                                             /* apply approx of A^{-1} */
  schur_bcgs_solve(sc,depth,A,G,M,B,schur_velo[0][0],schur_velo[1][0],schur_velo[2][0]); 

  for (i=0;i<m_b;i++)  d[0]->x[i] =schur_velo[0][0]->x[i] ;  /* copy velo result to solution vector */

                                                  /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    {                                            /* old defect in b[k] */
      dmatmul(A_times_update[k],A[k],B,d[k]);    /* current update in d[k] */
      
      numerator =  AMG_ddot(A_times_update[k],b[k]);
      nominator =  AMG_ddot(A_times_update[k],A_times_update[k]);
      if (nominator > 0)
        alpha = numerator/nominator;
      else
        {
          alpha = 0.05;
          sprintf(buf,"MESSAGE : Step length control failed on level %d, set alpha = %g !!\n",k,alpha);
          AMG_Print(buf);
        }
      if (fabs(alpha)<alpha_eps)
        {
          /* sprintf(buf,"MESSAGE : Step length factor %g, too small, set alpha = %g\n",
             alpha,alpha_eps);
             AMG_Print(buf);*/
          alpha = alpha_eps;
        }
          if (sc->verbose>2)
            {
              sprintf(buf,"alpha %g\n",alpha);
              AMG_Print(buf);
            }
    }
  else 
    alpha =sc->omega[0];                                 /* no step length control */

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */
  
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                            */
/* braess_sarazin_smoother                                                  */
/*                                                                            */
/* preconditioner for saddle point problems via a coupled multigrid              */
/* method with braess--sarazin smoother                                     */
/*                                                                            */
/* ( A B ) (u) = (r)                                                            */
/* ( B'0 ) (p) = (s)                                                            */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/* schur_press[0] : x,d[0] (in gmres solver)                                      */
/* schur_press[1] : z[0],r[0] (in gmres solver)                                    */
/* schur_velo[0][0] : x=A^-1b                                                    */
/* schur_velo[1][0] : rhs b of Ax=b                                             */
/* schur_velo[2][0] : defect b-Ax                                            */
/*                                                                            */
/****************************************************************************/
int prepare_braess_sarazin_smoother(AMG_SolverContext *sc,int *A_length,
                                    int *B_length, int depth)
{
  int i,k;

  for (k=0;k<2;k++)
    for (i=0;i<=depth;i++)
      schur_velo[k][i]= AMG_NewVector(A_length[i],AMG_MATRIX_B(A[i]),"schur_velo"); 
  for (k=0;k<4;k++)
    schur_press[k]= AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"schur_velo"); 
  for (k=0;k<3;k++)
    pres_help[k] = AMG_NewVector(B_length[0],AMG_MATRIX_B(A[0]),"pres_help"); 
  velo_rhs = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"velo_rhs"); 
  velo_help = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"velo_help"); 
  return(0);
}

int clear_braess_sarazin_smoother(int depth)
{
  int i,k;

  for (k=0;k<2;k++)
    for (i=0;i<=depth;i++)
      {
        free(schur_velo[k][i]->x);
        free(schur_velo[k][i]);
      }
  for (k=0;k<4;k++)
    {
      free(schur_press[k]->x);
      free(schur_press[k]);
    }
  for (k=0;k<3;k++)
    {
      free(pres_help[k]->x);
      free(pres_help[k]);
    }
  free(velo_rhs->x);
  free(velo_rhs);
  free(velo_help);
  
  free (B_Diri[0]->a);
  free (B_Diri[0]);
  free (B_Diri[1]->a);
  free (B_Diri[1]);
  
  return(0);
}


int braess_sarazin_cg (AMG_SolverContext *sc,
        AMG_MATRIX *A, AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,max_iter,n_b,m_b,m_bb;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last,alpha;
  
  max_iter=sc->schur_iteration_maxit;
  max_iter=200;
  switch(sc->system_type)
    {
    case SADDLE_1 :
      m_b = AMG_MATRIX_M(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_N(B[0]);                     /* # columns of B */
      break;
    case SADDLE_1_TYPE_1 :
      m_b = AMG_MATRIX_N(B[0]);                     /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);                     /* # columns of B */
      break;
    case SADDLE_2_TYPE_1 : 
      m_b = AMG_MATRIX_N(B[0])+ AMG_MATRIX_N(B[1]); /* # rows of B = rows, columns of A */
      n_b = AMG_MATRIX_M(B[0]);
      break;
    }

  n_b = x->n;  
  m_b = 2 * A->n;
  m_bb = A->n;
  schur_press[0]->n = n_b;             /* set current length for vectors */
  schur_press[1]->n = n_b;
  schur_press[2]->n = n_b;
  schur_press[3]->n = n_b;
  schur_velo[0][0]->n = m_b;  
  
  AMG_dset(x,0.0);                               /* initial iterate */
  for (i=0;i<n_b;i++) 
    { 
      schur_press[0]->x[i] = b->x[i];  /* rhs in pressure space */
      schur_press[2]->x[i] = b->x[i];  /* initial defect */
    }       
                                  
  rho_last = 1.0;        
  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(schur_press[2],schur_press[2]));
  if (sc->verbose>0)
    {
      sprintf(buf,"Entering Schur_CG_: %4d %12.4lE n %d na %d\n",0,dnorm,n_b,m_b);
      AMG_Print(buf);
    }
 
  if (dnorm<1e-10)
    return(0);
  for (i=0; i<max_iter; i++)
    {
      AMG_dset(schur_press[3],0.0);
      AMG_dcopy(schur_press[0],schur_press[2]);
      AMG_dcopy(schur_press[3],schur_press[2]);
      rho = AMG_ddot(schur_press[2],schur_press[3]);
      if (i>0) AMG_daxpy(schur_press[3],rho/rho_last,schur_press[1]);
      rho_last=rho;
      AMG_dcopy(schur_press[1],schur_press[3]);

      AMG_B_dmatmul_gal_SADDLE_2_TYPE_1 (schur_velo[0][0],A,B,schur_press[1]);
      /*for (j=0;j<n_b;j++)
        printf("schur_press[1](%d)= %g\n",j+1,schur_press[1]->x[j]);*/
      /*for (j=0;j<m_b;j++)
        printf("rhs_r(%d)= %g\n",j+1,rhs_r[0]->x[j]);*/
      for (j=0;j<m_bb;j++)
        {
          schur_velo[0][0]->x[j] =schur_velo[0][0]->x[j] / A->a[A->ra[j]];
          schur_velo[0][0]->x[j+m_b/2] =schur_velo[0][0]->x[j+m_b/2] / A->a[A->ra[j]];
        }
      AMG_B_Trans_dmatmul_gal_SADDLE_2_TYPE_1 (schur_press[0], A,B,schur_velo[0][0]);
      /*for (j=0;j<n_b;j++)
        printf("schur_press[0] %g\n",schur_press[0]->x[j]);*/
    
      alpha=rho/AMG_ddot(schur_press[1],schur_press[0]);
      /*printf("alpha %g\n",alpha);*/
      AMG_daxpy(x,alpha,schur_press[1]);
      AMG_daxpy(schur_press[2],-alpha,schur_press[0]);
      dnorm=sqrt(AMG_ddot(schur_press[2],schur_press[2]));
      if (sc->verbose>1)
        {
          /*sprintf(buf,"Schur_CG Iteration %4d %g %g %12.4lg\n",i+1,dnorm,sc->res_norm_min,dnorm/dnormlast);
          AMG_Print(buf)*/;
          /* for (j=0;j<n_b;j++)
            {
              printf("x(%d) = %g\n",j,x->x[j]);
            }*/
        }
      dnormlast=dnorm;
      if (sc->ex_maxit) continue;
      if (dnorm<dnorm0*1e-10) break;
    }
 
      if (sc->verbose>1)
        {
          sprintf(buf,"Schur_CG_Stop %4d %g\n",i+1,dnorm);
          AMG_Print(buf);
        }
 
  return(0);
}                        

int braess_sarazin_smoother (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{                                                    /* rhs of system to solve is d[k] */
                                                     /* do not change b */ 
  double alpha,sum,*help3;
  int i,m,n_a,n_b,l,start,end;

  printf ("entering braess--sarazin\n");
  alpha=1.0;
   
  m = A[k]->n;                       /* compute size of matrices */
  switch(sc->system_type)
    {
    case SADDLE_1_TYPE_1 :
      n_a = m;                       /* size of A  */
      break;
    case SADDLE_2_TYPE_1 : 
      n_a = 2*m;                     /* size of A  */
      break;
    }
  n_b =d[k]->n-n_a;                  /* size of B */  
  
  velo_rhs->n=n_a;                   /* length of vectors */
  velo_help->n=n_a;
  pres_help[0]->n=n_b;
  pres_help[1]->n=n_b;
  pres_help[2]->n=n_b;
  
  for (i=0;i<n_b;i++)                /* store pressure defect s */
    pres_help[0]->x[i]=d[k]->x[n_a+i];

  help3= AMG_Malloc(sizeof(double)*d[k]->n);

  m = A[k]->n;
  schur_velo[1][0]->n = m;        /* store diag of A, current length */
  for (i=0;i<m;i++)
    schur_velo[1][0]->x[i] = A[k]->a[A[k]->ra[i]];
    
  for(i=0;i<d[k]->n;i++)                  /* store defect */              
    help3[i] = d[k]->x[i];
    
  sum = 0.0;                 /* compute average of pressure defect */
  for (i=0;i<n_b;i++) 
    sum+=d[k]->x[i+n_a];
  sum = sum/n_b;
  printf("sum pressure defect %g \n",sum);
  for (i=0;i<n_b;i++)               /* average pressure defect (for solvability) */ 
    help3[i+n_a]= pres_help[2]->x[i] = d[k]->x[i+n_a] = d[k]->x[n_a+i]-sum; 
  
  printf("dnorm defect %g p %g\n",sqrt(AMG_ddot(d[k],d[k])),sqrt(AMG_ddot(pres_help[2],pres_help[2])));
  
  switch(sc->system_type)           /* divide velo defect by diag of A */ 
    {                               /* -C^-1(f-Au^k-Bp^k) */ 
    case SADDLE_1_TYPE_1 :
      for (i=0;i<m;i++)
        velo_rhs->x[i] = d[k]->x[i] / schur_velo[1][0]->x[i];
      break;
    case SADDLE_2_TYPE_1 :
      for (i=0;i<m;i++)
        {
          velo_rhs->x[i] = d[k]->x[i] / schur_velo[1][0]->x[i];
          velo_rhs->x[i+m] = d[k]->x[i+m] / schur_velo[1][0]->x[i];
        }
      break;
    }

  AMG_B_Trans_dmatmul_gal_SADDLE_2_TYPE_1 (pres_help[1], A[k],B,velo_rhs); /* -B^TC^-1 r */
  
  /*  for (i=0;i<pres_help[1]->n;i++)
    {printf("p0(%d) = %g %g\n",i,pres_help[1]->x[i],pres_help[0]->x[i]);}*/
  AMG_daxpy(pres_help[1],-alpha,pres_help[0]);      /* -B^TC^-1 r - alpha s */
  AMG_dset(pres_help[0],0.0);                       /* initialize solution of cg */
  AMG_dcopy(pres_help[2],pres_help[1]);
  
  /*sum = 0.0;                      
  for (i=0;i<n_b;i++)
    sum+=pres_help[2]->x[i];
  sum = sum/n_b;
  printf("sum %g \n",sum);*/
  /*for (i=0;i<n_b;i++)
    pres_help[2]->x[i]-=sum;*/
  
  /*for (i=0;i<pres_help[0]->n;i++)
    {
      printf("before cg def_p(%d) = %g\n",i+1,pres_help[2]->x[i]);
    }*/
  printf("schur %d 0 %d 1 %d 2 %d\n",schur_matrix[k]->n,pres_help[0]->n,pres_help[1]->n,pres_help[2]->n);
  for (i=0;i<0;i++)
    {
      printf("%d defect %g\n",i,AMG_ddot(pres_help[2],pres_help[2]));
      jac(sc,0,depth,&schur_matrix[k],G_schur,M,B,pres_help,&pres_help[1],&pres_help[2]);
      AMG_dcopy(pres_help[2],pres_help[1]);                       
      AMG_dmatminus(pres_help[2],schur_matrix[k],B,pres_help[0]);                /* compute new defect */      
    }
  
  braess_sarazin_cg(sc,A[k],B,pres_help[0],pres_help[2]); /* q^{k+1} */
  
  for (i=0;i<pres_help[0]->n;i++)
    { 
      x[k]->x[n_a+i]+=pres_help[0]->x[i];
      /*printf("p(%d) = %g\n",i,pres_help[0]->x[i]);*/
    }
  AMG_B_dmatmul_gal_SADDLE_2_TYPE_1(velo_rhs,A[k],B,pres_help[0]); /* B q^{k+1} */
  velo_help->x=d[k]->x;                        /* store -r */
  AMG_daxpy(velo_help,-1.0,velo_rhs);          /* -r-B q^{k+1} */        
  
  m=A[k]->n;
  
  switch(sc->system_type)                            /* divide velo-rhs by diag of A */ 
    {
    case SADDLE_1_TYPE_1 :
      for (i=0;i<m;i++)
        x[k]->x[i] += velo_help->x[i] / (schur_velo[1][0]->x[i]*alpha);
      break;
    case SADDLE_2_TYPE_1 :
      for (i=0;i<m;i++)
        {
          x[k]->x[i] += velo_help->x[i] / (schur_velo[1][0]->x[i]*alpha);
          x[k]->x[i+m] += velo_help->x[i+m] / (schur_velo[1][0]->x[i]*alpha);
          /*printf("u[%d] = %g u[%d] = %g\n",i,x[k]->x[i],i+m,x[k]->x[i+m]);*/
        }
      break;
    }

  /* for(i=0;i<n_a+n_b;i++)
    printf("x(%d) = %g rhs(%d) = %g  def(%d) =  %g \n",i,x[k]->x[i],i,help3[i],i,d[k]->x[i]);*/
  
  
  for(i=0;i<m;i++)
    {
      velo_rhs->x[i] = velo_help->x[i] / (schur_velo[1][0]->x[i]*alpha);
      velo_rhs->x[i+m] = velo_help->x[i+m] / (schur_velo[1][0]->x[i]*alpha);
    }
  AMG_B_Trans_dmatmul_gal_SADDLE_2_TYPE_1 (pres_help[1], A[k],B,velo_rhs);
  for(i=0;i<n_b;i++)
    pres_help[2]->x[i] = pres_help[0]->x[i];
  for(i=0;i<n_b;i++)
    d[k]->x[i+n_a] = pres_help[1]->x[i];
  /*for(i=0;i<n_b;i++)
    pres_help[0]->x[i] = x[k]->x[i+n_a];*/
  AMG_B_dmatmul_gal_SADDLE_2_TYPE_1 (velo_rhs,A[k],B,pres_help[2]);
  for(i=0;i<m;i++)
    {
      d[k]->x[i] = velo_rhs->x[i] + 
        velo_help->x[i]*alpha * schur_velo[1][0]->x[i]/(schur_velo[1][0]->x[i]*alpha);
      d[k]->x[i+m] = velo_rhs->x[i+m] + 
        velo_help->x[i+m]*alpha  * schur_velo[1][0]->x[i]/(schur_velo[1][0]->x[i]*alpha);
    }
  /*for(i=0;i<n_a+n_b;i++)
    printf("x(%d) = %g rhs(%d) = %g  def(%d) =  %g \n",i,x[k]->x[i],i,help3[i],i,d[k]->x[i]);*/
  sum=0.0;
  for(i=0;i<d[k]->n;i++)
    sum+=(help3[i]-d[k]->x[i])*(help3[i]-d[k]->x[i]); 
  printf("defectdnorm %g\n",sum);
  free(help3);
  
 return(AMG_OK);
}

int BuildDiagSchurComplement(AMG_SolverContext *sc,AMG_MATRIX *A, AMG_MATRIX **B) 
{
  double *schur_value,sum,*A_diag;
  int *schur_row,*schur_column,schur_n,i,j,k,kj,end,endj,start,startj,nonzeros,found;
  int row_nonzeros,tmp,startk,endk,row_k,mem;
 

  schur_n = B[0]->m;                               /* dimension of schur complement */
  schur_row = AMG_Malloc((schur_n+1)*sizeof(int));     /* allocate vector for row starts */
  mem=tmp = AMG_MAX(1000,(schur_n*schur_n)/4);
  schur_column = AMG_Malloc(tmp*sizeof(int)); /* allocate vector for column indices */
  for (i=0;i<tmp;i++)
    schur_column[i]=-1;
  schur_row[0]=0;
  nonzeros=0;
  A_diag = AMG_Malloc(A->m*sizeof(double));
  for (i=0;i<A->m;i++)
    A_diag[i] = A->a[A->ra[i]];
  /*for (i=0;i<A->m;i++)
    printf("diag[%d] = %g\n",i, A_diag[i]);*/

  switch(sc->system_type)
    {
    case SADDLE_1_TYPE_1 :                         /* scalar mortar type system */
      for (i=0;i<schur_n;i++)
        {                                          /* i-th row in B^T */
          row_nonzeros=0;
          start = B[0]->ra[i];     
          end = B[0]->ra[i+1]-1;
          for (k=start; k<=end; k++)
            {
              tmp = B[0]->ja[k];                           /* index to compare with */
              for (j=0;j<schur_n;j++)
                {
                  found=0;
                  /*printf("row %d, col %d, colc %d \n",i,tmp,j);*/
                  startj=schur_row[i];
                  for (kj=startj;kj<startj+row_nonzeros;kj++)
                    if (schur_column[kj]==j)        /* index already in list */
                      {/*printf("%d %d \n",schur_column[kj],j);*/
                      found=2;}
                  if (found==2) continue;
                  startj = B[0]->ra[j];             /* j-th column of B */    
                  endj = B[0]->ra[j+1]-1;
                  found=0;
                  for (kj=startj; kj<=endj; kj++)
                    if (B[0]->ja[kj]==tmp)          /* indices equal -> nonzero in schur complement */
                      {
                        
                        found=1;
                        break;
                      }
                    
                  if (found)
                    {
                      printf("found row %d, col %d, colc %d \n",i,tmp,j);
                      schur_column[nonzeros]=j;
                      nonzeros++;
                      row_nonzeros++;
                    }
                }
            }
          schur_row[i+1]=schur_row[i]+ row_nonzeros;
        }
                        
      break;
    case SADDLE_2_TYPE_1 :                         /* 2D Stokes system */
      for (i=0;i<schur_n;i++)
        {                                          /* i-th row in B^T */
          row_nonzeros=0;
          start = B[0]->ra[i];     
          end = B[0]->ra[i+1]-1;
          for (k=start; k<=end; k++)
            {
              tmp = B[0]->ja[k];                           /* index to compare with */
              for (j=0;j<schur_n;j++)
                {
                  found=0;
                  /*printf("row %d, col %d, colc %d row %d col %d \n",i,tmp,j,schur_row[i],schur_column[schur_row[i]]);*/
                  startj=schur_row[i];                     /* start of row i */
                  for (kj=startj;kj<startj+row_nonzeros;kj++)
                    { 
                      /*printf("%d %d \n",schur_column[kj],j);*/
                      if (schur_column[kj]==j)        /* index already in list */
                        {
                          /*printf("break %d %d \n",schur_column[kj],j);*/
                          found=2;
                        }
                    }
                  if (found==2) continue;
                  startj = B[0]->ra[j];             /* j-th column of B */    
                  endj = B[0]->ra[j+1]-1;
                  found=0;
                  for (kj=startj; kj<=endj; kj++)
                    if (B[0]->ja[kj]==tmp)          /* indices equal -> nonzero in schur complement */
                      {
                        
                        found=1;
                        break;
                      }
                    
                  if (found)
                    {
                      /*printf("found row %d, col %d, colc %d \n",i,tmp,j);*/
                      schur_column[nonzeros]=j;
                      nonzeros++;
                      row_nonzeros++;
                    }
                }
            }
          start = B[1]->ra[i];     
          end = B[1]->ra[i+1]-1;
          for (k=start; k<=end; k++)
            {
              tmp = B[1]->ja[k];                           /* index to compare with */
              for (j=0;j<schur_n;j++)
                {
                  found=0;
                  /*printf("11 row %d, col %d, colc %d \n",i,tmp,j);*/
                  startj=schur_row[i];
                  for (kj=startj;kj<startj+row_nonzeros;kj++)
                    if (schur_column[kj]==j)        /* index already in list */
                      {/*printf("11 %d %d \n",schur_column[kj],j);*/found=2;}
                  if (found==2) continue;
                  startj = B[1]->ra[j];             /* j-th column of B */    
                  endj = B[1]->ra[j+1]-1;
                  found=0;
                  for (kj=startj; kj<=endj; kj++)
                    if (B[1]->ja[kj]==tmp)          /* indices equal -> nonzero in schur complement */
                      {
                        
                        found=1;
                        break;
                      }
                    
                  if (found)
                    {
                      /*printf("11 found row %d, col %d, colc %d \n",i,tmp,j);*/
                      schur_column[nonzeros]=j;
                      nonzeros++;
                      row_nonzeros++;
                    }
                }
            }
          schur_row[i+1]=schur_row[i]+ row_nonzeros;
        }
      break;
    default :
       AMG_Print("SYSTEM TYPE UNKWOWN !!!\n");
      break;     
    }
  /*for (i=0;i<nonzeros;i++)
    printf("s[%d] = %d \n",i,schur_column[i]);
  for (i=0;i<schur_n;i++)
  printf("rs[%d] = %d %d\n",i,schur_row[i],schur_n);*/
  if (nonzeros>mem)
    {
      sprintf(buf,"not enough memory for schur matrix %d %d!!\n",mem,nonzeros);
      AMG_Print(buf);
      exit(4711);
    }
  if (sc->verbose>1)
    {
      sprintf(buf,"schur complement matrix nonzeros : %d \n",nonzeros);
      AMG_Print(buf);
    } 
  
  for (i=0;i<schur_n;i++)                      /* reorder columns */
    {                                          /* i-th row schur matrix*/
      start = schur_row[i];     
      end = schur_row[i+1];
      for (k=start; k<end; k++)
        {
          if (schur_column[k]==i)
            {
              schur_column[k] = schur_column[start]; /* diagonal entry first */
              schur_column[start] = end-start;       /* nnz in row */
              break;
            }
        }
    }

  /*for (i=0;i<nonzeros;i++)
    printf("ja[%d] = %d \n",i,schur_column[i]);*/
  /*for (i=0;i<schur_n;i++)
  printf("rs[%d] = %d %d\n",i,schur_row[i],schur_n);*/
                                                 
  if (realloc(schur_column,nonzeros*sizeof(int))==NULL)     /* reduce size of schur_column */ 
  {
      exit(4711);
  }
  /*for (i=0;i<nonzeros;i++)
    printf("ja[%d] = %d \n",i,schur_column[i]);*/

  schur_matrix[0]=AMG_Malloc(sizeof(AMG_MATRIX)); /* build schur matrix */
  schur_matrix[0]->m=schur_n;
  schur_matrix[0]->n=schur_n;
  schur_matrix[0]->b=1;
  schur_matrix[0]->bb=1;
  schur_matrix[0]->system_as_scalar=1;
  schur_matrix[0]->blocks_in_diag=1;
  schur_matrix[0]->bandwidth=-1;
  schur_matrix[0]->nonzeros=nonzeros;
  schur_matrix[0]->connections=nonzeros;
  schur_matrix[0]->ra=schur_row;
  schur_matrix[0]->ja=schur_column;
  schur_value = AMG_Malloc(nonzeros*sizeof(double)); /* allocate vector for values */
  schur_matrix[0]->a=schur_value;

  /*  A_diag = AMG_Malloc(A->m*sizeof(double));
  for (i=0;i<A->m;i++)
  A_diag[i] = A->a[A->ra[i]];*/
    

  /*for (i=0;i<nonzeros;i++)
    printf("s[%d] = %d \n",i,schur_matrix[0]->ja[i]);
  for (i=0;i<schur_n;i++)
    printf("rs[%d] = %d %d\n",i,schur_matrix[0]->ra[i],schur_n);
  for (i=0;i<A->m;i++)
  printf("diag[%d] = %g\n",i, A_diag[i]);*/

  switch(sc->system_type)
    {
    case SADDLE_1_TYPE_1 :                         /* scalar mortar type system */
      for (i=0;i<schur_n;i++)                      /* compute matrix values */
        {                                          /* i-th row B^T matrix*/
          start = B[0]->ra[i];                   
          end = B[0]->ra[i+1];
          sum = 0.0;                               /* diagonal entry */
          for (k=start; k<end; k++)
            sum+=B[0]->a[k]*B[0]->a[k]/A_diag[B[0]->ja[k]];
          schur_matrix[0]->a[schur_matrix[0]->ra[i]] = sum;
          
          startk= schur_matrix[0]->ra[i];
          endk= schur_matrix[0]->ra[i+1];
          for (k=startk+1; k<endk; k++)      /* multiply row i and schur_matrix->ja[k] */
            {
              sum=0.0;          
              row_k = schur_matrix[0]->ja[k];
              startj= B[0]->ra[row_k];
              endj= B[0]->ra[row_k+1];
              for (j=startj; j<endj; j++)  /* j-th column in row_k */ 
                {
                  for (kj=start; kj<end; kj++)  /* i-th column in row_k */ 
                    if (B[0]->ja[j]==B[0]->ja[kj])
                      {
                        sum+= B[0]->a[j]*B[0]->a[kj]/A_diag[B[0]->ja[kj]];
                        /*printf("a %d %d %g\n",k,kj,sum);*/
                      }
                }
              schur_matrix[0]->a[k] = sum;
              
            }
        }
      break;
    case SADDLE_2_TYPE_1 :                         /* 2D Stokes system */
       for (i=0;i<schur_n;i++)                     /* compute matrix values */
        {                                          /* i-th row schur  matrix */
          start = B[0]->ra[i];                     /* start of i-th row of B^T */                 
          end = B[0]->ra[i+1];                     /* end of i-th row of B^T */ 
          sum = 0.0;                               /* diagonal entry in schur matrix*/
          for (k=start; k<end; k++)                /* go trough i-th row of B^T */
            sum+=B[0]->a[k]*B[0]->a[k]/A_diag[B[0]->ja[k]];
            /*sum+=B[0]->a[k]*B[0]->a[k];*/
          schur_matrix[0]->a[schur_matrix[0]->ra[i]] = sum; /* diag at start of row */
          
          startk= schur_matrix[0]->ra[i];    /* go through i-th row schur  matrix, start */
          endk= schur_matrix[0]->ra[i+1];    /* end of i-th row schur  matrix*/
          for (k=startk+1; k<endk; k++)      /* multiply row i of B^T and row no. schur_matrix->ja[k] */
            {
              sum=0.0;          
              row_k = schur_matrix[0]->ja[k];/* row to multiply with i-th row */
              startj= B[0]->ra[row_k];       /* find row in B^T */
              endj= B[0]->ra[row_k+1];
              for (j=startj; j<endj; j++)    /* j-th column in row_k */ 
                {
                  for (kj=start; kj<end; kj++)     /* i-th column in row_k */ 
                    if (B[0]->ja[j]==B[0]->ja[kj]) /* column entries equal */
                      {
                        sum+= B[0]->a[j]*B[0]->a[kj]/A_diag[B[0]->ja[kj]];
                        /*sum+= B[0]->a[j]*B[0]->a[kj];*/
                        /*printf("a %d %d %g %g %g\n",k,kj,B[0]->a[j],B[0]->a[kj],sum);*/
                      }
                }
              schur_matrix[0]->a[k] = sum;
              
            }
        }
       for (i=0;i<schur_n;i++)                      /* compute matrix values */
        {                                          /* i-th row B^T matrix*/
          start = B[1]->ra[i];                   
          end = B[1]->ra[i+1];
          sum = 0.0;                               /* diagonal entry */
          for (k=start; k<end; k++)
            sum+=B[1]->a[k]*B[1]->a[k]/A_diag[B[0]->ja[k]];
          schur_matrix[0]->a[schur_matrix[0]->ra[i]] += sum;
          
          startk= schur_matrix[0]->ra[i];
          endk= schur_matrix[0]->ra[i+1];
          for (k=startk+1; k<endk; k++)      /* multiply row i and schur_matrix->ja[k] */
            {
              sum=0.0;          
              row_k = schur_matrix[0]->ja[k];
              startj= B[1]->ra[row_k];
              endj= B[1]->ra[row_k+1];
              for (j=startj; j<endj; j++)  /* j-th column in row_k */ 
                {
                  for (kj=start; kj<end; kj++)  /* i-th column in row_k */ 
                    if (B[1]->ja[j]==B[1]->ja[kj]){
                      sum+= B[1]->a[j]*B[1]->a[kj]/A_diag[B[1]->ja[kj]];
                      /*sum+= B[1]->a[j]*B[1]->a[kj];*/
                      /*printf("b %d %d %g\n",k,kj,sum);*/}
                }
              schur_matrix[0]->a[k] += sum;
              
            }
        }
     break;
    default :
      AMG_Print("SYSTEM TYPE UNKWOWN !!!\n");
      break;     
    }
  
  /*for (i=0;i<nonzeros;i++)
    printf("s[%d] = %d %g\n",i,schur_matrix[0]->ja[i],schur_matrix[0]->a[i]);
  for (i=0;i<schur_n;i++)
  printf("rs[%d] = %d %d\n",i,schur_matrix[0]->ra[i],schur_n);*/

  return(0);
}
int NewBuildDiagSchurComplement(AMG_SolverContext *sc,AMG_MATRIX *A, AMG_MATRIX **B) 
{
  double *A_diag;
  int *schur_row,*schur_column,schur_n,i,j,k,kj,end,endj,start,startj,nonzeros,found;
  int row_nonzeros,tmp,startk,endk,row_k,schur_m,*columns[2],n,*start_col,*row_ind;
 

  schur_n = B[0]->m;                                /* dimension of schur complement */
  schur_row = AMG_Malloc((schur_n+1)*sizeof(int));     /* allocate vector for row starts */
  tmp = AMG_MAX(1000,(schur_n*schur_n)/10);
  schur_column = AMG_Malloc(tmp*sizeof(int)); /* allocate vector for column indices */
  for (i=0;i<tmp;i++)
    schur_column[i]=-1;
  schur_row[0]=0;
  nonzeros=0;
  A_diag = AMG_Malloc(A->m*sizeof(double));
  for (i=0;i<A->m;i++)
    A_diag[i] = A->a[A->ra[i]];
  for (i=0;i<A->m;i++)
    printf("diag[%d] = %g\n",i, A_diag[i]);

  switch(sc->system_type)
    {
    case SADDLE_1_TYPE_1 :                         /* scalar mortar type system */

      break;
    case SADDLE_2_TYPE_1 :                         /* 2D Stokes system */
      schur_m = B[0]->n+B[1]->n;                        /* columns of B^T */
      start_col = AMG_Malloc((schur_m+1)*sizeof(int));  /* allocate vector for column starts of B */
      row_ind = AMG_Malloc((B[0]->nonzeros+B[1]->nonzeros)* sizeof(int));
      for (i=0;i<=schur_m;i++) start_col[i] = 0;        /* initialize */
      for (i=0;i<=B[0]->nonzeros+B[1]->nonzeros;i++) row_ind[i] = -1;        /* initialize */
      columns[0] = B[0]->ja;
      columns[1] = B[1]->ja;
      n = B[0]->nonzeros;
      for (i=0;i<n;i++)                           /* i-th row in B^T */
        {                                         /* count column i (in  start_col[i+1])*/
          start_col[columns[0][i]+1]++;        
        }
      n = B[1]->nonzeros;
      for (i=0;i<n;i++)                           /* i-th row in B^T */
        {                                               /* column start vector for rows */
          start_col[B[0]->n+columns[1][i]+1]++;        
        }
      for (i=0;i<=schur_m;i++) printf("%d %d\n",i,start_col[i]);        /* initialize */
      for (i=1;i<=schur_m;i++)  
        {
          start_col[i]+=start_col[i-1];                                 /* indices for column starts */
        }            
      for (i=0;i<=schur_m;i++) printf("%d %d\n",i,start_col[i]);        /* initialize */
      for (i=0;i<schur_n;i++)                              /* go through rows of B^T */
        {
          start = B[0]->ra[i];                   
          end = B[0]->ra[i+1];
          for (k=start;k<end;k++)                          /* column B[0]->ja[k] */
            {
              startk = start_col[B[0]->ja[k]];             /* start of column */ 
              while(row_ind[startk]!= -1) startk++;        /* look for next free entry */
              row_ind[startk] = i ;                        /* put row index i in this entry */
            }
        }
       for (i=0;i<schur_n;i++)                              /* go through rows of B^T */
        {
          start = B[1]->ra[i];                   
          end = B[1]->ra[i+1];
          for (k=start;k<end;k++)                          /* column B[0]->ja[k] */
            {
              startk = start_col[B[1]->ja[k]+B[0]->n];     /* start of column */ 
              while(row_ind[startk]!= -1) startk++;        /* look for next free entry */
              row_ind[startk] = i ;                        /* put row index i in this entry */
            }
        }
      for (i=0;i<B[0]->nonzeros+B[1]->nonzeros;i++) printf("row %d %d\n",i,row_ind[i]);        /* initialize */

      for (i=0;i<schur_n;i++)
        {                                          /* i-th row in B^T */
          row_nonzeros=0;
          start = B[0]->ra[i];     
          end = B[0]->ra[i+1]-1;
          for (k=start; k<=end; k++)
            {
              tmp = B[0]->ja[k];                           /* index to compare with */
              startk = start_col[tmp];                     /* look for rows  tmp-th column */
              endk = start_col[tmp+1];
              for (j=startk;j<endk;j++)                    /* indices of first block */
                {
                  kj = row_ind[j];
                  schur_column[nonzeros]=kj;
                  nonzeros++;
                  row_nonzeros++;
                }
              for (j=0;j<schur_n;j++)
                {
                  found=0;
                  printf("row %d, col %d, colc %d row %d col %d \n",i,tmp,j,schur_row[i],schur_column[schur_row[i]]);
                  startj=schur_row[i];                     /* start of row i */
                  for (kj=startj;kj<startj+row_nonzeros;kj++)
                    { printf("%d %d \n",schur_column[kj],j);
                    if (schur_column[kj]==j)        /* index already in list */
                      {printf("break %d %d \n",schur_column[kj],j);found=2;}}
                  if (found==2) continue;
                  startj = B[0]->ra[j];             /* j-th column of B */    
                  endj = B[0]->ra[j+1]-1;
                  found=0;
                  for (kj=startj; kj<=endj; kj++)
                    if (B[0]->ja[kj]==tmp)          /* indices equal -> nonzero in schur complement */
                      {
                        
                        found=1;
                        break;
                      }
                    
                  if (found)
                    {
                      printf("found row %d, col %d, colc %d \n",i,tmp,j);
                      schur_column[nonzeros]=j;
                      nonzeros++;
                      row_nonzeros++;
                    }
                }
            }
          start = B[1]->ra[i];     
          end = B[1]->ra[i+1]-1;
          for (k=start; k<=end; k++)
            {
              tmp = B[1]->ja[k];                           /* index to compare with */
              for (j=0;j<schur_n;j++)
                {
                  found=0;
                  printf("11 row %d, col %d, colc %d \n",i,tmp,j);
                  startj=schur_row[i];
                  for (kj=startj;kj<startj+row_nonzeros;kj++)
                    if (schur_column[kj]==j)        /* index already in list */
                      {printf("11 %d %d \n",schur_column[kj],j);found=2;}
                  if (found==2) continue;
                  startj = B[1]->ra[j];             /* j-th column of B */    
                  endj = B[1]->ra[j+1]-1;
                  found=0;
                  for (kj=startj; kj<=endj; kj++)
                    if (B[1]->ja[kj]==tmp)          /* indices equal -> nonzero in schur complement */
                      {
                        
                        found=1;
                        break;
                      }
                    
                  if (found)
                    {
                      printf("11 found row %d, col %d, colc %d \n",i,tmp,j);
                      schur_column[nonzeros]=j;
                      nonzeros++;
                      row_nonzeros++;
                    }
                }
            }
          schur_row[i+1]=schur_row[i]+ row_nonzeros;
        }

     exit(1);
    default :
       AMG_Print("SYSTEM TYPE UNKWOWN !!!\n");
      break;     
    }
  return(0);
}
