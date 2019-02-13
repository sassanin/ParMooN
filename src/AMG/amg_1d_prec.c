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
/* File:          amg_1d_prec.c                                             */
/*                                                                          */
/* Purpose:   preconditioners for scalar problems                           */
/*                                                                          */
/* Author:    Volker John                                                   */
/*            Otto--von--Guericke Universitaet Magdeburg                    */
/*            Institut fuer Analysis und Numerik                            */
/*            Postfach 4120                                                 */
/*            39016 Magdeburg                                               */
/*            email : volker.john@mathematik.uni-magdeburg.de               */
/*                                                                          */
/* History:   1998/02/19 start using this library for MooN_MD               */
/*                                                                          */
/* Remarks:   1998/02/24 ILU                                                */
/*              1998/02/27 step length control                              */
/*              1998/06/03 ILUT                                             */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*    system include files                                                  */
/*    application include files                                             */
/*                                                                          */
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



int AMG_sorfb_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_, double *omega);
int AMG_ilub_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_);
int AMG_iluf_trans (AMG_MATRIX *A, AMG_VECTOR *v_, AMG_VECTOR *d_);

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
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern IterProcPtr smoother,coarse_smoother;
extern IterProcPtr preconditioner,schur_preconditioner,preconditioner_trans;

extern MultProcPtr dmatmul,dmatminus,A_dmatmul,B_dmatmul,B_Trans_dmatmul,A_dmatminus,
  MGC_dmatminus;

extern RestrictProcPtr restriction;
extern InterpolationProcPtr interpolation;

extern double start_residual,end_residual,residuals[AMG_CONV_RATE_BACK];
extern int residual_cnt,iteration_cnt;
extern double elapsed_time,TIME_WRAP;
extern clock_t start,finish,start1,finish1;

extern int coarse_grid;          /* multigrid is on coarse grid */
extern int mgc_recursion[AMG_MAX_LEVELS];

extern AMG_VECTOR *old_def[AMG_MAX_LEVELS];        /* array of vectors for step length control*/
extern AMG_VECTOR *A_times_update[AMG_MAX_LEVELS]; /* array of vectors for step length control*/
extern AMG_MATRIX *A[AMG_MAX_LEVELS],*B[AMG_MAX_LEVELS];
extern AMG_MATRIX *schur_matrix[AMG_MAX_LEVELS];
extern AMG_GRAPH  *G_schur[AMG_MAX_LEVELS];
extern AMG_VECTOR *velo_prolong[AMG_MAX_LEVELS],*pres_prolong[AMG_MAX_LEVELS];
extern AMG_VECTOR *velo_result,*pres_result;
extern AMG_VECTOR *row_equilibration;

extern  char buf[128];

extern int *ratr, *jatr, *postr;

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

/****************************************************************************/
/*                                                                            */
/* GRID TRANSFER OPERATIONS                                                   */
/*                                                                            */
/****************************************************************************/

int pc_restrict (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse)
{

  int *ca=AMG_GRAPH_CA(g);                   /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;

  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse))
    {
      printf("pc_restrict : number of blocks mismatched !!!\n");
      exit(4711);
    }
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g))
     {
      printf("pc_restrict : vector length mismatched %d %d !!!\n",AMG_VECTOR_N(fine),AMG_GRAPH_N(g));
      exit(4711);
    }
 
  if (b==1)
    {
      for (i=0; i<nc; i++) c[i] = 0.0;
      for (i=0; i<nf; i++) c[ca[i]] += f[i];
    }
  else
    {
      for (i=0; i<nc; i++) c[i] = 0.0;
      for (i=0; i<nf; i++) c[ca[i/b]*b+(i%b)] += f[i];
    }
  return(AMG_OK);
}        

int pc_restrict_2d (AMG_GRAPH *g, AMG_GRAPH *g1,  AMG_VECTOR *fine, AMG_VECTOR *coarse)
{

  int *ca=AMG_GRAPH_CA(g);                   /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,j,lc,lf;
  
  lc=nc/2;
  lf=nf/2;
  
  /* initialize coarse grid vector */
  for (i=0; i<nc; i++) 
    c[i] = 0.0;

  /* compute coarse grid vector */
  for (i=0; i<lf; i++) 
  {
    j = ca[i];
    c[j] += f[i];
    c[j+lc] += f[i+lf];
  }
    
  return(AMG_OK);
}
        
int pc_restrict_3d (AMG_GRAPH *g, AMG_GRAPH *g1,  AMG_VECTOR *fine, AMG_VECTOR *coarse)
{

  int *ca=AMG_GRAPH_CA(g);                   /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,j,lc,lf, lc2, lf2;

  /* initialize coarse grid vector */
  for (i=0;i<nc;i++)
    c[i] = 0.0;

  /* compute coarse grid vector */
  lc=nc/3;
  lf=nf/3;  
  lc2 = 2*lc;
  lf2 = 2*lf;
  for (i=0; i<lf; i++) 
  {
    j=ca[i];
    c[j] += f[i];
    c[j+lc] += f[i+lf];
    c[j+lc2] += f[i+lf2];
  }
    
  return(AMG_OK);
}        

int pc_restrict_6d (AMG_GRAPH *g, AMG_GRAPH *g1,  AMG_VECTOR *fine, AMG_VECTOR *coarse)
{

  int *ca=AMG_GRAPH_CA(g);                   /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,j,lc,lf, lc2, lf2, lc3, lf3, lc4, lf4, lc5, lf5;

  /* initialize coarse grid vector */
  for (i=0; i<nc; i++) 
    c[i] = 0.0;

  lc=nc/6;
  lf=nf/6;  
  lf2 = 2*lf;
  lc2 = 2*lc;
  lf3 = 3*lf;
  lc3 = 3*lc;
  lf4 = 4*lf;
  lc4 = 4*lc;
  lf5 = 5*lf;
  lc5 = 5*lc;
  for (i=0; i<lf; i++) 
    {
      j=ca[i];
      c[j] += f[i];
      c[j+lc] += f[i+lf];
      c[j+lc2] += f[i+lf2];
      c[j+lc3] += f[i+lf3];
      c[j+lc4] += f[i+lf3];
      c[j+lc5] += f[i+lf5];
    }
    
  return(AMG_OK);
}        

int pc_restrict_saddle_2d (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse)
{

  int *ca_velo=AMG_GRAPH_CA(g);          /* cluster mapping information */
  int *ca_pres=AMG_GRAPH_CA(g1);          /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int clust_velo = g->clusters;          /* # coarse velo dof (one component) */
  int clust_pres = g1->clusters;         /* # coarse pressure dof */
  int i,n_velo,n_pres;
  
  n_velo = AMG_GRAPH_N(g)/2;               /* velocity dof on fine grid */
  n_pres = AMG_GRAPH_N(g1);              /* pressure dof on fine grid */
  /* consistency */
  if (AMG_VECTOR_N(fine)!=2*AMG_GRAPH_N(g)+AMG_GRAPH_N(g1)) 
    {
      AMG_Print("Error in  pc_restrict_saddle_2d!\n");
      exit(4711);
    }
  
  for (i=0; i<nc; i++) c[i] = 0.0;                /* initialize coarse vector */
  for (i=0; i<n_velo; i++) 
    {
      c[ca_velo[i]] += f[i];                      /* restrict velocity, 1. component */
      c[ca_velo[i]+clust_velo] += f[i+n_velo];            /* restrict velocity, 2. component */
    }
  n_velo = 2 * n_velo;
  clust_velo = 2 * clust_velo;
  for (i=0; i<n_pres; i++) c[ca_pres[i]+clust_velo] += f[i+n_velo]; /* restrict pressure */

   return(AMG_OK);
}        

int pc_prolongate (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{

  int *ca=AMG_GRAPH_CA(g);                   /* cluster mapping information */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;
  double om;
  
  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse)) return(AMG_FATAL);
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g)) return(AMG_FATAL);
  
  if (b==1)
    {
      om=damp[0];
      for (i=0; i<nf; i++)  f[i] += c[ca[i]]*om;
    }
  else
    {
      for (i=0; i<nf; i++) f[i] += c[ca[i/b]*b+(i%b)]*damp[i%b];
    }
  return(AMG_OK);
}        


int pc_prolongate_auto (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);               /* cluster mapping information */
  float *da=AMG_GRAPH_DA(g);             /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i;
  double w1,w2;
  
  /* consistency */
  if (AMG_VECTOR_B(fine)!=AMG_VECTOR_B(coarse)) 
    {
      printf("pc_prolongate_auto : number of blocks mismatched !!!\n");
      exit(4711);
    }
  if (AMG_VECTOR_N(fine)!=AMG_GRAPH_N(g)) 
    {
      printf("pc_prolongate_auto : vector sizes mismatched !!!\n");
      exit(4711);
    }
  
  if (b==1)
    {
      w1=2.0-damp[0]; w2=damp[0]-1;
      for (i=0; i<nf; i++)  f[i] += c[ca[i]]*(w1+w2*da[i]);
    }
  else
    {
      for (i=0; i<nf; i++) f[i] += c[ca[i/b]*b+(i%b)]*damp[i%b];
    }
  return(AMG_OK);
}        
int pc_prolongate_auto_2d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);               /* cluster mapping information */
  float *da=AMG_GRAPH_DA(g);             /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,lc,lf;
  double w1,w2;
  
  lf=nf/2;lc=nc/2;
  w1=2.0-damp[0]; 
  w2=damp[0]-1;
  for (i=0; i<lf; i++) 
    { 
      f[i] += c[ca[i]]*(w1+w2*da[i]);
      f[i+lf] += c[ca[i]+nc/2]*(w1+w2*da[i]);
    }
   return(AMG_OK);
}        

int pc_prolongate_auto_3d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, 
                           AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);               /* cluster mapping information */
  float *da=AMG_GRAPH_DA(g);             /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,j,lc,lf, lf2, lc2;
  double w1,w2,w3;
  
  lf=nf/3;
  lc=nc/3;
  w1=2.0-damp[0]; 
  w2=damp[0]-1;
  lf2 = 2*lf;
  lc2 = 2*lc;
  for (i=0; i<lf; i++) 
  { 
    j = ca[i];
    w3 = w1+w2*da[i];
    f[i] += c[j]*w3;
    f[i+lf] += c[j+lc]*w3;
    f[i+lf2] += c[j+lc2]*w3;
  }
   return(AMG_OK);
}        

int pc_prolongate_auto_6d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, 
                           AMG_VECTOR *coarse, double *damp)
{
  int *ca=AMG_GRAPH_CA(g);               /* cluster mapping information */
  float *da=AMG_GRAPH_DA(g);             /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int i,j,lc,lf, lf2, lc2, lf3, lc3, lf4, lc4, lf5, lc5;
  double w1,w2,w3;
  
  w1=2.0-damp[0]; 
  w2=damp[0]-1;
  lf=nf/6;
  lc=nc/6;
  lf2 = 2*lf;
  lc2 = 2*lc;
  lf3 = 3*lf;
  lc3 = 3*lc;
  lf4 = 4*lf;
  lc4 = 4*lc;
  lf5 = 5*lf;
  lc5 = 5*lc;
  for (i=0; i<lf; i++) 
  { 
    j = ca[i];
    w3 = w1+w2*da[i];
    f[i] += c[j]*w3;
    f[i+lf] += c[j+lc]*w3;
    f[i+lf2] += c[j+lc2]*w3;
    f[i+lf3] += c[j+lc3]*w3;
    f[i+lf4] += c[j+lc4]*w3;
    f[i+lf5] += c[j+lc5]*w3;
  }
   return(AMG_OK);
}        

int pc_prolongate_auto_saddle_2d (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp)
{
  int *ca_velo=AMG_GRAPH_CA(g);               /* cluster mapping information */
  float *da_velo=AMG_GRAPH_DA(g);             /* auto damping factor */
  int *ca_pres=AMG_GRAPH_CA(g1);               /* cluster mapping information */
  float *da_pres=AMG_GRAPH_DA(g);             /* auto damping factor */
  double *f=AMG_VECTOR_X(fine);
  double *c=AMG_VECTOR_X(coarse);
  int b=AMG_VECTOR_B(fine);
  int nf=b*AMG_VECTOR_N(fine);
  int nc=b*AMG_VECTOR_N(coarse);
  int clust_velo = g->clusters;          /* # coarse velo dof (one component) */
  int clust_pres = g1->clusters;         /* # coarse pressure dof */
  int i,n_velo,n_pres;
  double w1,w2;
  
  n_velo = AMG_GRAPH_N(g);
  n_pres = AMG_GRAPH_N(g1);
  
  /* consistency */
  if (AMG_VECTOR_N(fine)!=2*AMG_GRAPH_N(g)+AMG_GRAPH_N(g1)) 
    {
      AMG_Print("Error in pc_prolongate_auto_saddle_2d!\n");
      exit(4711);
    }
  
  w1=2.0-damp[0]; 
  w2=damp[0]-1;
  for (i=0; i<n_velo; i++) 
    { 
      f[i] += c[ca_velo[i]]*(w1+w2*da_velo[i]);
      f[i+ n_velo] += c[ca_velo[i]+clust_velo]*(w1+w2*da_velo[i]);
    }
  n_velo = 2* n_velo; 
  clust_velo = 2*clust_velo; 
  for (i=0; i<n_pres; i++)  f[i+n_velo] += c[ca_pres[i]+clust_velo]*(w1+w2*da_pres[i]);
     
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                            */
/* PRECONDITIONERS (SMOOTHERS)                                                */
/*                                                                            */
/****************************************************************************/
int prepare_steplength_control(int *A_length, int *B_length, int depth)
{
  int k;
 
  for (k=0; k<=depth; k++)
    {
      old_def[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"old_def");
      A_times_update[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"A_times_update");
    }
  return(0);
}

int clear_steplength_control(int depth)
{
  int k;
  for (k=0; k<=depth; k++)
    {
      free(old_def[k]->x);
      free(old_def[k]);
      free(A_times_update[k]->x);
      free(A_times_update[k]);
    }
  return(0);
}
int initialize_steplength_control(AMG_SolverContext *sc,AMG_VECTOR *d, int k,int *old_def_size)
{
  if (old_def[k]->n>d->n)
    {
      *old_def_size=old_def[k]->n;
      old_def[k]->n=A_times_update[k]->n=d->n;
    }
  if (old_def[k]->n<d->n)
    {
      AMG_Print("Not enough memory for step length control !!!\n");
      AMG_Print("Step length control switced off !!!\n");
      sc->step_length_control_all=sc->step_length_control_fine=0;
    }
  AMG_dcopy(old_def[k],d);
  return(0);
}


double steplength_control(AMG_MATRIX *A,AMG_VECTOR *d,int k,int old_def_size)
{
  double alpha,nominator,numerator,alpha_eps=1.0e-4;

  A_dmatmul(A_times_update[k],A,B,d);
  numerator =  AMG_ddot(A_times_update[k],old_def[k]);
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
  if (old_def_size>-1)
    {
      A_times_update[k]->n=old_def[k]->n= old_def_size;
    }
  return(alpha);
}

/* SSOR without possibility of step lenght control or damping */
int ssor  (AMG_SolverContext *sc, int k, int depth,
          AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
          AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
          AMG_VECTOR *x[AMG_MAX_LEVELS],
          AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
                                               
  AMG_sorfb(A[k],x[k],b[k],&sc->sor_omega);           /* compute new iterate */        
  return(AMG_OK);
}

/* SSOR with possibility of damping or step lenght control */
int ssor_slc (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                /* defect d=b-Ax is valid on entry */
                                                /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);
  
  AMG_dcopy(d[k],x[k]);                           /* save old iterate */
  AMG_sorfb(A[k],x[k],b[k],&sc->sor_omega);       /* compute new iterate */        
  AMG_daxpby(d[k],1.0,x[k],-1.0,d[k]);            /* compute update */

                                                  /* if step length control */          
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
        alpha = sc->omega[0];                    /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }

  AMG_daxpy(x[k],alpha-1.0,d[k]);                        /* update solution x */

  return(AMG_OK);
}
            
int sor (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                /* defect d=b-Ax is valid on entry */  
                                                /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);

  AMG_sorf(A[k],d[k],d[k],&sc->sor_omega);           /* compute correction, overwrite d */        
 
                                                /* if step length control */      
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
        alpha = sc->omega[0];             /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */

  return(AMG_OK);
}
                
int jac (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                /* defect d=b-Ax is valid on entry */  
                                                /* if step length control should be applied */ 
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);

  AMG_jac(A[k],d[k],d[k],&sc->sor_omega);                 /* compute correction, overwrite d */        
  
                                                /* if step length control */          
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
  {
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  }
  else
  {
    if (!coarse_grid) 
      alpha = sc->omega[0];             /* set fixed damping factor */
    else
      alpha = sc->omega_coarse[0];
  }
  
  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */
  return(AMG_OK);
}
                
int ex (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int i,old_def_size=-1;
                                                 /* defect d=b-Ax is valid on entry */
                                                 /* if step length control should be applied */
 
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
   initialize_steplength_control(sc,d[k],k,&old_def_size);    

  AMG_EXApplyLU(AMG_MATRIX_A(M[k]),AMG_MATRIX_BW(M[k]),AMG_MATRIX_N(M[k]),AMG_MATRIX_BLDI(M[k]),AMG_VECTOR_X(d[k]));
                                              /* if step length control */          
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
        alpha = sc->omega[0];             /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */

  return(AMG_OK);
}

int prepare_ilu(AMG_SolverContext *sc,AMG_MATRIX *A[AMG_MAX_LEVELS],
                 AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth)
{
  int k;

  for (k=fine;k<=depth;k++)
    M[k] = ILUDecomposition(sc,A[k]);
  return(0);
}
int clear_ilu(AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth)
{
  int k;
  for (k=fine;k<=depth;k++)
  {
    free(M[k]->a);
    free(M[k]);
  }
  return(0);
}
int ilu (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                       /* defect d=b-Ax is valid on entry */
                                                /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);

  AMG_iluf(M[k],d[k],d[k]);                     /* forward step, overwrite d */        
  AMG_ilub(M[k],d[k],d[k]);                     /* backward step overwrite d */

                                               /* if step length control */          
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
            alpha = sc->omega[0];             /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }     

  AMG_daxpy(x[k],alpha,d[k]);
        
  return(AMG_OK);
}
/****************************************************************************/
/*                                                                            */
/* PREPARATION OF PRECONDITIONERS (SMOOTHERS)                               */
/*                                                                            */
/****************************************************************************/
AMG_MATRIX *prepare_ex (AMG_MATRIX *A)
{
  int bw,rl;
  int i,n=A->n,k,start,end;
  int *ra=A->ra, *ja=A->ja;
  double *a=A->a,*lu;
  AMG_MATRIX *new;
  char buf[128];

  /* compute bandwith */
  bw=0;
  for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      for (k=start+1; k<end; k++)
        {bw=AMG_MAX(bw,AMG_ABS(i-ja[k]));}
    }
  
  /*if (sc->verbose>1)*/
  {
    sprintf(buf,"EXACT SOLVER : band width %d\n",bw);
    AMG_Print(buf);
  }
  
  /* allocate new matrix */
  new = AMG_NewMatrix(n,1,n*(2*bw+1),AMG_MATRIX_SAS(A),AMG_MATRIX_BLDI(A),1,1,"ex matrix");
  if (new==NULL) return(new);
  lu=new->a;
  AMG_MATRIX_BW(new)=bw;
  
  /* insert & copy entries */
  for (i=0; i<n*(2*bw+1); i++) lu[i]=0.0;
  for (i=0; i<n; i++)
    {
      start = ra[i]; end = start+ja[start];
      AMG_EX_MAT(lu,bw,i,i) = a[start];
      for (k=start+1; k<end; k++)
        AMG_EX_MAT(lu,bw,i,ja[k]) = a[k];
    }
  
  /* decompose */
  if (AMG_EXDecomposeMatrix(lu,bw,n)) return(AMG_NULL);
  
  /* return matrix */
  return(new);
}
int clear_ex(AMG_MATRIX *M[AMG_MAX_LEVELS],int depth)
{
  free(M[depth]->ra);
  free(M[depth]->ja);
  free(M[depth]->a);
  free(M[depth]);
  return(0);
}

AMG_MATRIX *ILUDecomposition(AMG_SolverContext *sc, AMG_MATRIX *A)
{
  
  int i,j,k,l,row_start,row_end,row_start_j,row_end_j,found;
  int n=A->n,b=A->b,bb=A->bb,nonzeros=A->nonzeros;
  int *ra=A->ra, *ja=A->ja;
  double ilu_eps =1e-10,diag_ii,update,beta_ilu=sc->ilu_beta,pivot;
  AMG_MATRIX *M;
  char buf[128];

  if (sc->verbose>1)
    {
      sprintf(buf,"Entering ILU Decomposition\n");
      AMG_Print(buf);
    }
  M = AMG_NewMatrix(n,b,nonzeros,A->system_as_scalar,AMG_MATRIX_BLDI(A),0,0,"ilu");   /* allocate matrix for ilu deco */
  if (M==AMG_NULL) 
    {
      AMG_Print("ILU decomposition failed, cannot allocate ilu matrix !\n");
      exit(4711);
    }

  M->ra = ra;
  M->ja = ja;
  for (k=0; k<nonzeros*bb; k++) M->a[k] = A->a[k];
 
  for (i=0;i<n;i++)                                             /* loop over all rows */
    {
      diag_ii = M->a[ra[i]];                                    /* diagonal entry in row i */ 
      if (fabs(diag_ii)<ilu_eps)                                /* diagonal entry is zero */
        {
          printf("i %d %g %g\n",i,M->a[ra[i]],A->a[ra[i]]);
          AMG_Print("ILU decomposition failed by zero diagonal entry !\n");
          exit(4711);
        }
      row_start = ra[i];                                        /* start of row i in array M->a */          
      row_end = ra[i]+ja[row_start];                            /* end of row i in array M->a */
      for (j=row_start+1;j<row_end;j++)                         /* loop over columns of row i */
        {
          if (ja[j]>i)                                          /* upper triangular entry */
            {
              row_start_j = ra[ja[j]];                          /* find row ja[j] */
              row_end_j = ra[ja[j]] + ja[row_start_j];
              found = 0;
              for (k=row_start_j+1;k<row_end_j;k++)             /* search a(ja[j],i) */
                {
                  if (ja[k]!=i) continue;
                  pivot = M->a[k]/diag_ii;                      /* pivot element, store in L */
                  M->a[k] = pivot;
                  found = 1;
                  break;
                }
              if (!found) continue;                             /* row of Dirichlet dof */ 
              M->a[row_start_j] -= pivot * M->a[j];             /* diagonal in row ja[j], is in P(A) */
              for (k=row_start+1;k<row_end;k++)                 /* update row ja[j] */
                {
                  if (ja[k]<=ja[j]) continue;                   /* left part already done */
                  update = M->a[k]*pivot;
                  for (l=row_start_j+1;l<row_end_j;l++)         /* check if (k,l) in P(A) */
                    {
                      if (ja[k]==ja[l])                         /* same column, compute coeff of U */
                        {
                          M->a[l] -= update;
                          update = 0.0;
                          break;
                        }
                    }
                  
                  M->a[row_start_j] += beta_ilu * fabs(update); /* if not in pattern, add to diagonal */ 
                }
            }
        }
    }
  if (sc->verbose>1)
    {
      sprintf(buf,"Leaving ILU Decomposition\n");
      AMG_Print(buf);
    }
  return(M);
}
int quick_split_1(int begin_index,int length, double compare,double *w, int *w_ind)
    {
      int tmp,i,j,right;
      double tmp_1;
    
      right=begin_index-1;
     
      for (i=begin_index;i<begin_index+length;i++)
        {
          if (fabs(w[i])<compare)
            {
              for (j=right;j<i;j++)
                {
                  tmp_1=w[j+1];
                  w[j+1]=w[right];
                  w[right]=tmp_1;
                  tmp=w_ind[j+1];
                  w_ind[j+1]=w_ind[right];
                  w_ind[right]=tmp;
                }
              right++;
            }
        }
       
      return(right);
    }

int quick_split_0(int begin_index,int length, double compare,double *w, int *w_ind)
    {
      int tmp,i,j,right;
      double tmp_1;
    
      right=begin_index-1;
     
      for (i=begin_index;i<begin_index+length;i++)
        {
          if (fabs(w[i])<compare)
            {
              tmp_1=w[i];
              w[i]=w[right];
              w[right]=tmp_1;
              tmp=w_ind[i];
              w_ind[i]=w_ind[right];
              w_ind[right]=tmp;
              right++;
            }
        }
      return(right);
    }

int quick_split_2(int first, int last, double compare, double *w, int *w_ind)
{
     int tmp,j,right;
     double itmp;

     right=first;
     
     for (j=right;j>=last;j--)
       {
         if(fabs(w[j]) > compare)
           {
             right--;
             tmp=w[right];
             itmp=w_ind[right];
             w[right]=w[j];
             w_ind[right]=w_ind[j];
             w[j]=tmp;
             w_ind[j]=itmp;
           }
       }

     if (first!=right)
       {
         tmp=w[right];
         itmp=w_ind[right];
         w[right]=w[first];
         w_ind[right]=w_ind[first];
         w[first]=tmp;
         w_ind[first]=itmp;
       }

     return(right);
}
int clear_ilut(AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth)
{
  int k;
  for (k=fine;k<=depth;k++)
    {
      free(M[k]->a);
      free(M[k]->ra);
      free(M[k]->ja);
      free(M[k]);
    }
  return(0);
}

AMG_MATRIX *ILUTDecomposition(AMG_SolverContext *sc, AMG_MATRIX *A)
{
  
  int i,j,k,l,row_start,row_end,row_start_j,row_end_j,found,p,*w_ind,ij;
  int n=A->n,b=A->b,bb=A->bb,nonzeros=A->nonzeros;
  int *ra=A->ra, *ja=A->ja;
  int A_lower, A_upper, M_lower, M_upper, M_index,ilut_sort;
  int mid,left,begin_index,length,first,last,itmp;
  double ilu_eps =1e-10,diag_ii,update,beta_ilu=sc->ilu_beta,pivot,*w, min_w,t;
  double compare,norm,tol_norm,drop_1,tmp;
  AMG_MATRIX *M;
  char buf[128];
  
  ilut_sort = sc->ilut_sort;
  p = sc->ilut_absolute_fillin;
  drop_1 = sc->ilut_relative_fillin;
  if (drop_1*n/100<p)
    p = (int)floor(drop_1*n/100);
  if (p==0) p=1;

  tol_norm = sc->ilut_tol;
  
  if (sc->verbose>1)
    {
      sprintf(buf,"Entering ILUT Decomposition\n");
      AMG_Print(buf);
    }
  w_ind=malloc(n*sizeof(int));
  w=malloc(n*sizeof(double));                                          /* allocate vector w */
  
  M = AMG_NewMatrix(n,b,nonzeros+2*p*n,A->system_as_scalar,
                    A->blocks_in_diag,1,1,"ilut");/* allocate matrix for ilu deco */
  if (M==AMG_NULL) 
    {
      AMG_Print("ILUT decomposition failed, cannot allocate ilu matrix !\n");
      return(AMG_NULL);
    }
  M_index=0;
  for (i=0;i<n;i++)
    {
      w[i]=0;                                                     /* Fuelle w mit Nullen */
      w_ind[i]=i;
   }
  
  for (i=0;i<n;i++)                                             /* loop over all rows */
    {
      
      finish = clock();
      if (finish>=start)
        elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
      else
        elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
      start = clock();
      row_start = ra[i];                                        /* start of i-th row of A */
      row_end=row_start + ja[row_start];                        /* end of i-th row of A */
      w[i]=A->a[row_start];                                     /* a_ii */
      norm=w[i]*w[i];
      A_lower=A_upper=0;                                        /* initialize counters */
      for (j=row_start+1;j<row_end;j++)                         /* go through row i */
        {
          w[ja[j]]=A->a[j];                                     /* copy row i into w */
          norm+=w[ja[j]]*w[ja[j]];                              /* compute norm */
          if (ja[j]>i)
            A_upper++;                                          /* nonzero in U */
          else
            A_lower++;                                          /* nonzero in L */
        }
      drop_1=tol_norm*sqrt(norm);                               /* dropping tolerance */
      
      for (k=0;k<i;k++)                                         /* go throu all lines above i */
        {
          if (w[k]==0) continue;                                /* if w[k] == 0 continue */
          if (M->a[M->ra[k]]==0) 
            {                        
              AMG_Print("ILUT : Division durch Null\n");
              exit(4711);
            }
          else
            {
              w[k]=w[k]/M->a[M->ra[k]];                            /* compute entry in L */
            }                               
          
          if (fabs(w[k])<drop_1)                                /* apply a dropping rule */          
            { 
              w[k]=0;
              continue;
            }
          
          row_start = M->ra[k];                                 /* start of k-th row in M */
          row_end = row_start + M->ja[row_start];               /* end of k-th row in M */
          for (j=row_start+1;j<row_end;j++)                     /* go through k-th row in M */
            if (M->ja[j]>k)                                     /* if column greater than k */
              {                                                 /* i.e. if element in U */
                 w[M->ja[j]]=w[M->ja[j]]-w[k]*M->a[j];           /* update w=w-w[k]*u_k */ 
              }
        
        } /* end k */
      
      /* apply a dropping rule */
      switch(ilut_sort)
        {
        case ILUT_QUICK_SPLIT_0 :        
        case ILUT_QUICK_SPLIT_1 :        
          
          mid=left = 0;                                         /* first "split point" */
          for (k=0;k<i;k++)                                     /* check the lower part of w */
            {
              if (fabs(w[k])==0)                                /* look for zero entries in w */
                {
                  if (w[mid]==0)                                /* if split point zero */
                    mid++;                                      /* define new split point */
                  else
                    {
                      w[k]=w[mid];
                      w_ind[k]=w_ind[mid];
                      w[mid]=0;
                      mid++;
                    }
                }
            }                                                   /* all zero entries at the beginning of w */
          length = i-mid-1;
          begin_index = mid+1;
          compare = fabs(w[mid]);
          while(1)
            { 
              if (mid >= i-A_lower-p) break;                    /* found the p largest entries in w */
              length = i-1-mid;                                 /* length of the vector which we want to check */
              begin_index = mid+1;                              /* starting point */
              compare = fabs(w[mid]);                           /* entry used for comparison */
              if (ilut_sort ==ILUT_QUICK_SPLIT_1)               /* which routine ? */
                mid=quick_split_1(begin_index,length,compare,w,w_ind);  /* starting quick_split_1 */
              else 
                mid=quick_split_0(begin_index,length,compare,w,w_ind); /* starting quick_split_0 */
              
              if (mid == i-A_lower-p) break;                    /* found the p largest entries */
              
              if (ilut_sort ==ILUT_QUICK_SPLIT_1)               /* criterion when using quick_split */
                {
                  if ((mid < i-A_lower-p) && (mid==begin_index-1))   /* compare too small */
                    mid++;
                }
              else                                              /* criterion when using quick_split_old */
                if (mid < i-A_lower-p)                          /* compare too small */
                  mid++;
              
              if (mid > i-A_lower-p)                            /* compare too large */
                mid=begin_index-1;
              
            }
          
          mid=left = i+1;                                       /* check the upper part */
          for (k=i+1;k<n;k++)                                   /* same procedure like for the lower part */
            {
              if (fabs(w[k])==0)
                {
                  if (w[mid]==0)
                    mid++;
                  else
                    {
                      w[k]=w[mid];
                      w_ind[k]=w_ind[mid];
                      w[mid]=0;
                      mid++;
                    }
                }
            }
          length = n-mid-1;
          begin_index = mid+1;
          compare = fabs(w[mid]);
          while(1)
            { 
              if (mid >= n-A_upper-p) break;
              if (A_upper+p==n-mid) break;
              
              length = n-mid-1;
              begin_index = mid+1;
              compare = fabs(w[mid]);
              
              if (ilut_sort ==ILUT_QUICK_SPLIT_1)
                mid=quick_split_1(begin_index,length,compare,w,w_ind);
              else 
                mid=quick_split_0(begin_index,length,compare,w,w_ind);
              
              
               if (mid == n-A_upper-p) break;
              
              if (ilut_sort ==ILUT_QUICK_SPLIT_1)
                {
                  if ((mid < n-A_upper-p) && (mid==begin_index-1))   /* compare too small */
                    mid++;
                }
              else
                if (mid < n-A_upper-p)                              /* compare too small */
                  mid++;
              
              if (mid > n-A_upper-p)                                /* compare too large */
                mid=begin_index-1;
              
            }
          
          M_upper=M_lower=0;
          M->ra[i]=M_index;
          M->a[M_index]=w[i];
          found=M_index;
          M_index++;
          for (k=i-1;(k>i-1-A_lower-p) && (k>=0) ;k--)              /* loop over p largest entries of the lower part*/
            { 
              if (fabs(w[k])==0) continue;
              M->a[M_index]=w[k];                                   /* update M */
              M->ja[M_index]=w_ind[k];                              /* update M */
              M_index++;
              M_lower++;
            }
          for (k=n-1;(k>n-1-A_upper-p)  && (k > i) ;k--)            /* loop over p largest entries of the upper part */
            {
              if (fabs(w[k])==0) continue;
              M->a[M_index]=w[k];                                   /* update M */
              M->ja[M_index]=w_ind[k];                              /* update M */
              M_index++;
              M_upper++;
            }
          M->ja[found]=M_lower+M_upper+1;
          for (k=0;k<n;k++)                                         /* loop over all entries */
            {
              w[k]=0;                                               /* fill with zeros*/
              w_ind[k]=k;
            }
          break;

        case ILUT_QUICK_SPLIT_2 :
          first=i-1;                /* end of lower triangle */
          last=0;                   /* begin of lower triangle */ 


          while(1)                  /* skip leading zeros */
            {
              if (w[last]==0)
                last++;
              else
                break;
            }

          for (j=last;j<i;j++)      /* put all zeros to begin of vector */       
            {
              if (fabs(w[j])==0)
                {                     
                  w[j]=w[last];
                  w_ind[j]=w_ind[last];
                  w[last]=0;
                  last++;                    
                }
            }

          if (last<i-A_lower-p)           /* too much nonzeros */ 
            while(1)
              {                            
                mid=first;                /* set split point */
                compare=fabs(w[mid]);
                
                for (j=first-1;j>=last;j--)  /* go backward */
                  {
                    if (fabs(w[j])>compare)  /* larger than split point */
                      {
                        mid--;               /* change smaller value next to split point */ 
                        tmp=w[mid];          /* with larger value */
                        w[mid]=w[j];         /* larger value stands now left of split point */
                        w[j]=tmp;
                        itmp=w_ind[mid];
                        w_ind[mid]=w_ind[j];
                        w_ind[j]=itmp;
                      }
                  }
                
                tmp=w[mid];                    /* change most left of larger values */
                w[mid]=w[first];               /* with split point */  
                w[first]=tmp;                  /* values right of former split point are */
                itmp=w_ind[mid];               /* larger, former split point is on position mid */
                w_ind[mid]=w_ind[first];       /* values left of mid are smaller than former */
                w_ind[first]=itmp;             /* split point (or equal) */
                
                if (mid==i-A_lower-p) break;   /* i-A_lower-p largest values found */
                                               /* w[mid],...,w[i-1] */ 
                if (mid>i-A_lower-p)           /* not enough largest values found */
                  first=mid-1;                 /* search in lower sub vector again */
                else                           /* too much largest values found */
                  last=mid+1;                       /* search in upper sub vector again */
            }

          first=n-1;    /* end of upper triangle */
          last=i+1;     /* begin of upper triangle */

          while(1)      /* skip leading zeros */     
            {
              if (w[last]==0)
                last++;
              else
                break;
            }
          
          for (j=last;j<n;j++)     /* put zeros at begin of the vector */
            {
              if (fabs(w[j])==0)
                {
                  w[j]=w[last];
                  w_ind[j]=w_ind[last];
                  w[last]=0;
                  last++;
                }
            }
            

          if (last<n-A_upper-p)     /* too much nonzeros */
            while(1)
              {
                for (j=0;j<n;j++)
        
                mid=first;          /* split point */
                compare=fabs(w[mid]);
              
                for (j=first-1;j>last-1;j--)
                  {
                    if (fabs(w[j])>compare) /* put larger values than */
                      {                     /* of split point next to the */
                        mid--;              /* left of split point */
                        tmp=w[mid];
                        w[mid]=w[j];
                        w[j]=tmp;
                        itmp=w_ind[mid];
                        w_ind[mid]=w_ind[j];
                        w_ind[j]=itmp;
                      }
                  }
              
                tmp=w[mid];                 /* put split point left to all */
                w[mid]=w[first];            /* values larger than split point */
                w[first]=tmp;
                itmp=w_ind[mid];
                w_ind[mid]=w_ind[first];
                w_ind[first]=itmp; 
                            
                if (mid==n-A_upper-p) break; /* the A_upper+p largest values found */
        
              
                if (mid>n-A_upper-p)   /* not enough largest values found */     
                  first=mid-1;         /* search in lower sub vector again */
                else                   /* too much largest values found */
                  last=mid+1;          /* search in upper sub vector again */
              }
          

          M_upper=M_lower=0;
          M->ra[i]=M_index;
          M->a[M_index]=w[i];
          found=M_index;
          M_index++;
          for (k=i-1;(k>i-1-A_lower-p) && (k>=0) ;k--)              /* loop over p largest entries of the lower part*/
            { 
              if (fabs(w[k])==0) continue;
              M->a[M_index]=w[k];                                   /* update M */
              M->ja[M_index]=w_ind[k];                              /* update M */
              M_index++;
              M_lower++;
            }
          for (k=n-1;(k>n-1-A_upper-p)  && (k > i) ;k--)            /* loop over p largest entries of the upper part */
            {
              if (fabs(w[k])==0) continue;
              M->a[M_index]=w[k];                                   /* update M */
              M->ja[M_index]=w_ind[k];                              /* update M */
              M_index++;
              M_upper++;
            }
          M->ja[found]=M_lower+M_upper+1;
          for (k=0;k<n;k++)                                         /* loop over all entries */
            {
              w[k]=0;                                               /* fill with zeros*/
              w_ind[k]=k;
            }

          break;
        default :
          AMG_Print("No selection routine for ILUT found\n");
          exit(4711);
        }
    } /* end i */

  free(w);
  free(w_ind);

  if (sc->verbose>1)
    {
      sprintf(buf,"Leaving ILUT Decomposition\n");
      AMG_Print(buf);
    }
  return(M);
}

int AverageRHS(AMG_MATRIX *A,AMG_VECTOR *b,AMG_VECTOR *d)
{
  int n_a,n_b,i;
  double sum;

  n_a = 2*A->n;
  n_b = b->n-n_a;
  sum = 0.0;
  for (i=0;i<n_b;i++)
    sum+=b->x[i+n_a];
  sum = sum/n_b;
  for (i=0;i<n_b;i++)
    b->x[i+n_a] = d->x[i+n_a] -= sum;

  return(0);
}

/****************************************************************************/
/*                                                                            */
/* MULTIGRID CYCLE                                                          */
/*                                                                            */
/****************************************************************************/
int clear_mgc(AMG_MATRIX *A[AMG_MAX_LEVELS],AMG_GRAPH *G[AMG_MAX_LEVELS],int depth)
{
  int k;
 
  for (k=1;k<depth;k++)                        
    {
      free(A[k]->ra);
      free(A[k]->ja);
      free(A[k]->a);
      free(A[k]);
      free(G[k]->na);
      free(G[k]->la);
      free(G[k]->da);
      free(G[k]->ca);
      free(G[k]);
    }
 
  if (depth!=0)
    {
      free(G[0]->ca);
      free(G[0]->na);
      free(G[0]->la);
      free(G[0]->da);
      free(G[0]);           
      free(A[depth]->ra);
      free(A[depth]->ja);
      free(A[depth]->a);
      free(A[depth]);
    }

  return(0);
}
int prepare_coupled_mgc(AMG_SolverContext *sc,int *A_length,int *B_length, int depth)
{
  int k;
  
  velo_result = AMG_NewVector(A_length[0],AMG_MATRIX_B(A[0]),"velo_result");
  pres_result = AMG_NewVector(B_length[0],AMG_MATRIX_B(schur_matrix[0]),"pres_result");
  for (k=0;k<=depth;k++)                  
    {
      velo_prolong[k] = AMG_NewVector(A_length[k],AMG_MATRIX_B(A[k]),"velo_prolong");
      pres_prolong[k] = AMG_NewVector(B_length[k],AMG_MATRIX_B(schur_matrix[k]),"pres_prolong");
    }
  
  return(0);
}
int clear_coupled_mgc(AMG_MATRIX *A[AMG_MAX_LEVELS],AMG_GRAPH *G[AMG_MAX_LEVELS],int depth)
{
  int k;

   for (k=1;k<depth;k++)                        
    {
      free(A[k]->ra);
      free(A[k]->ja);
      free(A[k]->a);
      free(A[k]);
      free(G[k]->na);
      free(G[k]->la);
      free(G[k]->da);
      free(G[k]->ca);
      free(G[k]);
      free(schur_matrix[k]->ra);
      free(schur_matrix[k]->ja);
      free(schur_matrix[k]->a);
      free(schur_matrix[k]);
      free(G_schur[k]->na);
      free(G_schur[k]->la);
      free(G_schur[k]->da);
      free(G_schur[k]->ca);
      free(G_schur[k]);
    }
 
  if (depth!=0)
    {
      free(G[0]->ca);
      free(G[0]->na);
      free(G[0]->la);
      free(G[0]->da);
      free(G[0]);           
      free(G_schur[0]->ca);
      free(G_schur[0]->na);
      free(G_schur[0]->la);
      free(G_schur[0]->da);
      free(G_schur[0]);           
      free(A[depth]->ra);
      free(A[depth]->ja);
      free(A[depth]->a);
      free(A[depth]);
      free(schur_matrix[depth]->ra);
      free(schur_matrix[depth]->ja);
      free(schur_matrix[depth]->a);
      free(schur_matrix[depth]);
    }
  free(velo_result->x);
  free(velo_result);
  free(pres_result->x);
  free(pres_result);
  for (k=0;k<=depth;k++)                  
    {
      free(velo_prolong[k]->x);
      free(velo_prolong[k]);
      free(pres_prolong[k]->x);
      free(pres_prolong[k]);
    }
  return(0);
}

int mgc (AMG_SolverContext *sc, int k, int depth,                                 
                AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],                 
                AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
                AMG_VECTOR *x[AMG_MAX_LEVELS],          
                AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])

{
  int i;
  double dnorm, dnorm0;
  
                                                  /* defect d=b-Ax is valid on entry */
  if (k==depth)                               /* coarse grid solve */
    {
      if ((sc->system_type!=SCALAR1)||(sc->system_type!=SCALAR3)) 
        if (sc->preconditioner==AMG_MGC)
          AverageRHS(A[k],b[k],d[k]);
      dnorm0=sqrt(AMG_ddot(d[k],d[k]))*sc->coarse_red_factor;
      coarse_grid = 1;
      for (i=0; i<sc->coarse_maxit; i++)
        {
          coarse_smoother(sc,k,depth,A,G,M,B,x,b,d);
          AMG_dcopy(d[k],b[k]);                                /* recompute d=b-Ax */
          MGC_dmatminus(d[k],A[k],B,x[k]);
          dnorm=sqrt(AMG_ddot(d[k],d[k]));
          /*sprintf(buf,"coarse dnorm %g\n",dnorm);
            AMG_Print(buf);*/          
          if (dnorm<=dnorm0) 
            {
              coarse_grid = 0;
	      i++;
              break;
            }
        }
      if (sc->verbose>2) 
      {
        sprintf(buf,"coarse grid solver, iterations %d residual %g\n",i,dnorm);
	AMG_Print(buf);
      }
      coarse_grid = 0;
      return(0);                              /* to avoid the final 'if' in this routine */
    }
  else
    {                                         /* pre smoothing */
      if (sc->verbose>2) 
      {
        sprintf(buf,"restriction before, level %d residual %g\n",k,sqrt(AMG_ddot(d[k],d[k])));
	AMG_Print(buf);
      }
      if (sc->n1[k]>1)
        dnorm0=sqrt(AMG_ddot(d[k],d[k]))*sc->smoother_red_factor;
      for (i=0; i<sc->n1[k]; i++) 
        {
          smoother(sc,k,depth,A,G,M,B,x,b,d);
          AMG_dcopy(d[k],b[k]);                                /* recompute d=b-Ax */
          MGC_dmatminus(d[k],A[k],B,x[k]);
          if (i<sc->n1[k])
            if (sqrt(AMG_ddot(d[k],d[k])) <= dnorm0)
	    {
		i++;
		break;
	    }
        }
      if (sc->verbose>2) 
      {
        sprintf(buf,"restriction, level %d iterations %d residual %g\n",k,i,sqrt(AMG_ddot(d[k],d[k])));
	AMG_Print(buf);
      }
      restriction(G[k],G_schur[k],d[k],b[k+1]);
      AMG_dcopy(d[k+1],b[k+1]);
      AMG_dset(x[k+1],0.0);
      for (i=0; i<AMG_MIN(mgc_recursion[k],depth-k); i++) /* coarsest grid only once */
        {
          mgc(sc,k+1,depth,A,G,M,B,x,b,d);
          if (i+1==AMG_MIN(mgc_recursion[k],depth-k)) break; /* after coarsest grid */
          AMG_dcopy(d[k+1],b[k+1]);            /* because d must cont defect on entry */
          MGC_dmatminus(d[k+1],A[k+1],B,x[k+1]);
        }
      if (sc->gamma<1) mgc_recursion[k] = 1;              /* F--cycle */
      interpolation(G[k],G_schur[k],x[k],x[k+1],sc->omega_p); /* post smoothing */ 
      for (i=0; i<sc->n2[k]; i++) 
        {
          AMG_dcopy(d[k],b[k]);                                /* recompute d=b-Ax */
          MGC_dmatminus(d[k],A[k],B,x[k]);
          if ((i==0)&&(sc->n2[k]>1)) 
            dnorm0=sqrt(AMG_ddot(d[k],d[k]))*sc->smoother_red_factor;
          else
            {
              if (i<sc->n2[k]-1)
                if (sqrt(AMG_ddot(d[k],d[k])) <= dnorm0)
		{ 
		    i--;
		    break;
		}
            }
          smoother(sc,k,depth,A,G,M,B,x,b,d);
        }
      if (sc->verbose>2) 
      {
        sprintf(buf,"prolongation, level %d iterations %d residual %g\n",k,i,sqrt(AMG_ddot(d[k],d[k])));
	AMG_Print(buf);
      }
    }

  if ((sc->gamma<1)&&(k==0))            /* set parameters for next F--cycle */
    for (i=0;i<AMG_MAX_LEVELS;i++)        
      mgc_recursion[i] = 2;

  return(AMG_OK);
}

int ApplyRowEquilibration(AMG_SolverContext *sc,AMG_MATRIX *A)
{
  int i,start,end,*ra,*ja,A_length,k;
  double s,*a;

  A_length =  AMG_MATRIX_N(A);
  row_equilibration = AMG_NewVector(A_length,AMG_MATRIX_B(A),"row_equi");

  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  
  for (i=0; i<A_length; i++)
    {
      row_equilibration->x[i] = 0;         /* compute equilibration factor */
      start = ra[i]; end = start+ja[start];
      s = fabs(a[start]);
      for (k=start+1; k<end; k++) 
        s += fabs(a[k]);
      row_equilibration->x[i] = s; 
      a[start]/=s;                         /* apply equilibration */
      for (k=start+1; k<end; k++) 
        a[k]/=s;
    }
      
  return(0);
}
int EquilibrateRHS(AMG_VECTOR *b,AMG_VECTOR *equi)
{
  int i,n;
  n = AMG_VECTOR_N(b);
  for (i=0;i<n;i++)
    b->x[i]/=equi->x[i];      
  return(0);
}
int DisApplyRowEquilibration(AMG_SolverContext *sc,AMG_MATRIX *A,AMG_VECTOR *b)
{
  int i,start,end,*ra,*ja,A_length,k;
  double s,*a;

  a = AMG_MATRIX_A(A);
  ra = AMG_MATRIX_RA(A);
  ja = AMG_MATRIX_JA(A);
  A_length =  AMG_MATRIX_N(A);
  
  for (i=0; i<A_length; i++)
    {
      s = row_equilibration->x[i];
      start = ra[i]; end = start+ja[start];
      a[start]*=s;                         /* apply equilibration */
      for (k=start+1; k<end; k++) 
        a[k]*=s;
      b->x[i]*=s;
    }
      
  free(row_equilibration->x);
  free(row_equilibration);
 
  return(0);
}

int sor_trans (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                /* defect d=b-Ax is valid on entry */  
                                                /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);

  if (A[k]->ratr==NULL)
    ComputeArraysForTransposedMatrix(A[k]);

  AMG_sorf_trans(A[k],d[k],d[k],&sc->sor_omega);           /* compute correction, overwrite d */        
 
                                                /* if step length control */      
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
        alpha = sc->omega[0];             /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }

  AMG_daxpy(x[k],alpha,d[k]);                        /* update solution x */

  return(AMG_OK);
}
/* SSOR with transposed matrix */
int ssor_trans  (AMG_SolverContext *sc, int k, int depth,
          AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
          AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
          AMG_VECTOR *x[AMG_MAX_LEVELS],
          AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
                                               
  AMG_sorfb_trans(A[k],x[k],b[k],&sc->sor_omega);           /* compute new iterate */        
/*  AMG_sorf_trans(A[k],x[k],b[k],&sc->sor_omega); */          /* compute new iterate */        
  return(AMG_OK);
}

int ilu_trans (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS])
{
  double alpha;
  int old_def_size=-1;
                                                       /* defect d=b-Ax is valid on entry */
                                                /* if step length control should be applied */
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    initialize_steplength_control(sc,d[k],k,&old_def_size);
  
  if (M[k]->ratr==NULL)
    ComputeArraysForTransposedMatrix(M[k]);

  AMG_ilub_trans(M[k],d[k],d[k]);                     /* backward step, overwrite d */        
  AMG_iluf_trans(M[k],d[k],d[k]);                     /* forward step overwrite d */

                                               /* if step length control */          
  if ((sc->step_length_control_all)||
      (((sc->step_length_control_fine)&&(k==0))))
    alpha = steplength_control(A[k],d[k],k,old_def_size);
  else
    {
      if (!coarse_grid) 
            alpha = sc->omega[0];             /* set fixed damping factor */
      else
        alpha = sc->omega_coarse[0];
    }     

  AMG_daxpy(x[k],alpha,d[k]);
        
  return(AMG_OK);
}

int ComputeArraysForTransposedMatrix(AMG_MATRIX *A)
{
  int i,start,end,*ra,*ja,A_length,k,entries;
  int k1,start1,end1,col;

  ra = AMG_MATRIX_RA(A);             /* row pointer */
  ja = AMG_MATRIX_JA(A);             /* column indices */
  A_length = AMG_MATRIX_N(A);        /* dimension of A */
  entries = AMG_MATRIX_NONZEROS(A);  /* number of nonzeros of A */

  ratr = AMG_Malloc((A_length)*sizeof(int));
  for (i=0; i<A_length; i++)
      ratr[i] = 0;

  jatr = AMG_Malloc(entries*sizeof(int));
  for (i=0; i<entries; i++)
      jatr[i] = -1;

  postr = AMG_Malloc(entries*sizeof(int));
  

  for (i=0; i<A_length; i++)         /* count the entries in each column */
  {
    start = ra[i]; 
    end = start+ja[start];
    ratr[i]++;                        /* diagonal element */
    for (k=start+1; k<end; k++) 
        ratr[ja[k]]++;
  }

  k = 0;
  for (i=0; i<A_length;i++)          /* store info in accumulated form */
  {
      jatr[k] = ratr[i];
      k += ratr[i];
      ratr[i] = k-ratr[i];
  }

  for (i=0; i<A_length; i++)         /* fill arrays */
  {
    start = ra[i]; 
    end = start+ja[start];
    postr[ratr[i]] = start;          /* diagonal element */
    for (k=start+1; k<end; k++) 
    {
	col = ja[k];
	start1 = ratr[col];
	end1 = start1 + jatr[start1];
	for (k1=start1+1; k1<end1; k1++)
	{
	    if (jatr[k1]==-1)
	    {
		jatr[k1] = i;
		postr[k1] = k;
		break;
	    }
	}
    }
  }
/* for (i=0;i<A_length;i++)
      printf("%d ",ra[i]);
  printf("\n");
 for (i=0;i<entries;i++)
      printf("%d ",ja[i]);
  printf("\n");
  for (i=0;i<A_length;i++)
      printf("%d ",ratr[i]);
  printf("\n");
 for (i=0;i<entries;i++)
      printf("%d ",jatr[i]);
  printf("\n");
 for (i=0;i<entries;i++)
      printf("%d ",postr[i]);
  printf("\n");
*/

  A->ratr = ratr;
  A->jatr = jatr;
  A->postr = postr;

  return(0);
}
