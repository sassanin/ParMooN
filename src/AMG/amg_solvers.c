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
/* File:          amg_solve_solver.c                                        */
/*                                                                          */
/* Purpose:   contains the solvers                                          */
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
/* Remarks:   1998/03/02 GMRES                                              */
/*              1998/03/12 flexible GMRES                                   */
/*              1998/04/01 exact solver                                     */
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
#include <amg_2d_prec.h>
#include <amg_solvers.h>

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
/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern IterProcPtr smoother,coarse_smoother;
extern IterProcPtr preconditioner,schur_preconditioner,preconditioner_trans;


extern MultProcPtr dmatmul,dmatminus,A_dmatmul,B_dmatmul,B_Trans_dmatmul,dmattransmul;

extern double start_residual,end_residual,residuals[AMG_CONV_RATE_BACK];
extern int residual_cnt,iteration_cnt;
extern double elapsed_time,TIME_WRAP;
extern clock_t start,finish,start1,finish1;

extern int coarse_grid;                            /* multigrid is on coarse grid */

extern AMG_VECTOR *old_def[AMG_MAX_LEVELS];        /* array of vectors for step length control*/
extern AMG_VECTOR *A_times_update[AMG_MAX_LEVELS]; /* array of vectors for step length control*/
extern AMG_MATRIX *B[AMG_MAX_LEVELS];

char buffer[128];

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* SOLVERS                                                                  */
/*                                                                          */
/****************************************************************************/

/* global data for solvers */
extern        AMG_MATRIX *A[AMG_MAX_LEVELS];
extern        AMG_MATRIX *M[AMG_MAX_LEVELS];
extern        AMG_GRAPH  *G[AMG_MAX_LEVELS];
extern        AMG_VECTOR *x[AMG_MAX_LEVELS];
extern        AMG_VECTOR *b[AMG_MAX_LEVELS];
extern        AMG_VECTOR *d[AMG_MAX_LEVELS];
extern        AMG_VECTOR *z[AMG_MAX_LEVELS];
extern        AMG_VECTOR *r[AMG_MAX_LEVELS];
extern        AMG_VECTOR *q;
extern        AMG_VECTOR *p[AMG_MAX_LEVELS];
extern        AMG_VECTOR *w;
extern        AMG_VECTOR *old_def[AMG_MAX_LEVELS];
extern        AMG_VECTOR *A_times_update[AMG_MAX_LEVELS];
extern        AMG_VECTOR *s,*cosi,*sn;
extern        AMG_VECTOR *H[AMG_MAX_GMRES_RESTART+1];
extern        AMG_VECTOR *v[AMG_MAX_GMRES_RESTART+1];
extern        AMG_VECTOR *zv[AMG_MAX_GMRES_RESTART+1];
extern        int depth;
extern  AMG_CoarsenContext *global_cc;
extern  AMG_SolverContext *global_sc;

extern AMG_VECTOR *u,*v1,*w,*w1,*t[AMG_MAX_LEVELS],*tilde_r;
/****************************************************************************/
/*                                                                          */
/* Linear Solver                                                            */
/*                                                                          */
/****************************************************************************/
int prepare_ls_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */

  for (k=1; k<=depth; k++)
    {
      x[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"x");
      b[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"b");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
    }
  return(0);
}

int clear_ls_solve(AMG_SolverContext *sc,int depth)
{
  int k;
   
  free(d[0]->x);
  free(d[0]);
 
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(x[k]->x);
      free(x[k]);
      free(d[k]->x);
      free(d[k]);
      free(b[k]->x);
      free(b[k]);
    }
  return(0);
}

int ls_solve (AMG_VECTOR *x_in, AMG_VECTOR *b_in)
{
  int i,maxit;
  double dnorm, dnorm0, dnormlast;
  AMG_SolverContext *sc=global_sc;

  
  x[0] = x_in; 
  b[0] = b_in;
  /* iterate */
  AMG_dcopy(d[0],b[0]);
  dmatminus(d[0],A[0],B,x[0]);
  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(d[0],d[0]));
  start_residual=dnorm;
  residuals[0]=dnorm;
  residual_cnt++;
  
  if (sc->verbose>0)
    {
      sprintf(buffer,"FixedPoint Iteration %4d %12.4lE \n",0,dnorm);
      AMG_Print(buffer);
    }
  if (dnorm>sc->res_norm_min)
    maxit=sc->maxit;
  else
    maxit=1;
    
  for (i=1; i<=maxit; i++)
    {
      preconditioner(sc,0,depth,A,G,M,B,x,b,d);
      AMG_dcopy(d[0],b[0]);
      dmatminus(d[0],A[0],B,x[0]);
      dnorm=sqrt(AMG_ddot(d[0],d[0]));
      if (sc->verbose>0)
        {
          sprintf(buffer,"FixedPoint Iteration %4d %12.4lE %12.4lg\n",i,dnorm,dnorm/dnormlast);
          AMG_Print(buffer);
        }
      dnormlast=dnorm;
      if (i<sc->maxit)
        {
          residuals[residual_cnt%AMG_CONV_RATE_BACK]=dnorm;
          residual_cnt++;
        }                
      finish = clock();
      if (finish>=start)
        elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
      else
        elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
      start = clock();
      if (sc->ex_maxit) continue;
      if (dnorm<dnorm0*sc->red_factor) break;
      if (dnorm<sc->res_norm_min) break;
      if (dnorm>sc->div_factor*dnorm0)
        {
          AMG_Print("FixedPoint iteration diverges !!!\n");
          exit(4711);
        }
    }
  if (i>=sc->maxit && !sc->ex_maxit && dnorm>sc->res_norm_min &&dnorm>dnorm0*sc->red_factor )
    {
      AMG_Print("solver not converged\n");
      sprintf(buffer,"FixedPoint Iteration: (maximal) iterations %d residual %g \n",i-1,dnorm);
      AMG_Print(buffer);
      iteration_cnt=i;
      end_residual=dnorm;                
      return(-1);
    }
  sprintf(buffer,"FixedPoint iteration : iterations %d residual %g \n",AMG_MIN(i,maxit),dnorm);
  AMG_Print(buffer);
  iteration_cnt=i;
  end_residual=dnorm;                
  return(i);
}                        

/****************************************************************************/
/*                                                                            */
/* Conjugate Gradient                                                       */
/*                                                                            */
/****************************************************************************/
int prepare_cg_solve(AMG_SolverContext *sc, int *A_length, int *B_length, 
		     int depth)
{
  int k, n_active, n, *ra, *ja, start, end, i, l;
  double *a, *bx;

  z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  q = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"q");

  if (sc->preconditioner==AMG_MGC) /* to avoid unnecessary allocation for schur complement methods */
  {
      for (k=1; k<=depth; k++)
      {
	  z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
	  r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
	  d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
      }
  }
  return(0);
  
  // modification of A to get a symmetric matrix
  // tried: 2005-10-07, there were no differences in the behavior of the method
  // insert Dirichlet values and make matrix symmetric
  n_active = 0;

  a = AMG_MATRIX_A(A[0]);
  ra = AMG_MATRIX_RA(A[0]);
  ja = AMG_MATRIX_JA(A[0]);
  n = AMG_MATRIX_N(A[0]);
  for (i=0; i<n; i++)
  {
      start = ra[i]; end = start+ja[start];
      if (fabs(a[start]-1)>1e-6)
      {
	  n_active++;
	  continue;
      }
      // main diagonal element is one, it might be a Dirichlet node
      for (k=start+1; k<end; k++) 
      {
	  if (fabs(a[k])>1e-9)
	      continue;
      }
      // first Dirichlet node found
      break;
  }
  sprintf(buffer,"CG: active d.o.f. %d\n",n_active);
  AMG_Print(buffer);

  // insert Dirichlet nodes
  bx = AMG_VECTOR_X(b[0]);
  printf("%d\n",AMG_VECTOR_N(b[0]));
  return(0);
  for (i=0; i<n_active; i++)
  {
      start = ra[i]; 
      end = start+ja[start];
     
      for  (k=start+1; k<end; k++) 
      {
	  // column index
	  l = ja[k];
	  printf("b %d %d\n",i,l);
	  if (l>n_active)
	  {
	      bx[i] -= a[k]*bx[l];
	      a[k] = 0.0;
	  }
      }
  }
  return(0);
}

int clear_cg_solve(AMG_SolverContext *sc,int depth)
{
  int k;

  free(z[0]->x);
  free(z[0]);
  free(d[0]->x);
  free(d[0]);
  free(q->x);
  free(q);
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(z[k]->x);
      free(z[k]);
      free(d[k]->x);
      free(d[k]);
      free(r[k]->x);
      free(r[k]);
    }
  return(0);
}

int cg_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last=1.0,alpha;
  AMG_SolverContext *sc=global_sc;
  
                                              /* iterate */
  r[0] = b;                                   /* overwrite r with residual */
  dmatminus(r[0],A[0],B,x);                   /* compute residual */        
  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r[0],r[0]));

  if (sc->verbose>0)
  {
    sprintf(buffer,"Entering CG: %4d %12.4lE \n",0,dnorm);
    AMG_Print(buffer);
  }
  if (dnorm<=sc->res_norm_min)/* stopping criterion fulfilled */
  {
    if ((sc->minit==0)||(dnorm<1e-20))
    {
      sprintf(buffer,"CG : iterations 0 residual %g \n",dnorm);
      AMG_Print(buffer);
      return(0);
    }
    else
      sc->maxit = sc->minit;
  }
  
  AMG_dset(q,0.0);

  for (i=0; i<sc->maxit; i++)
  {
    if (sc->amg_prec_it)                      /* preconditioner should be applied */
    {
      AMG_dset(z[0],0.0);
      AMG_dcopy(d[0],r[0]);
      for (j=0;j<sc->amg_prec_it;j++)
      {
	preconditioner(sc,0,depth,A,G,M,B,z,r,d);      
	if (j<sc->amg_prec_it-1)
	{
	  AMG_dcopy(d[0],r[0]);
	  dmatminus(d[0],A[0],B,z[0]);
	}
      }
    }
    else
      AMG_dcopy(z[0],r[0]);
    rho = AMG_ddot(r[0],z[0]);
    AMG_daxpby(q,1.0,z[0],rho/rho_last,q);
    rho_last=rho;
    dmatmul(d[0],A[0],B,q);
    alpha=rho/AMG_ddot(q,d[0]);
    AMG_daxpy(x,alpha,q);
    AMG_daxpy(r[0],-alpha,d[0]);
    dnorm=sqrt(AMG_ddot(r[0],r[0]));
    if (sc->verbose>0)
    {
      sprintf(buffer,"CG Iteration %4d %g %g %12.4lg\n",i+1,dnorm,sc->res_norm_min,dnorm/dnormlast);
      AMG_Print(buffer);
    }
    dnormlast=dnorm;
    if (i<sc->minit) continue;
    if (sc->ex_maxit) continue;
    if (dnorm<dnorm0*sc->red_factor) break;
    if (dnorm<sc->res_norm_min) break;
    if (dnorm>sc->div_factor*dnorm0)
    {
      AMG_Print("CG iteration diverges !!!\n");
      exit(4711);
    }
  }
  if (i==sc->maxit && !sc->ex_maxit)
  {
    AMG_Print("solver not converged\n");
    sprintf(buffer,"CG : (maximal) iterations %d residual %g reduction %g\n",i,dnorm,dnorm/dnorm0);
    AMG_Print(buffer);
    return(-1);
  }
  sprintf(buffer,"CG : iterations %d residual %g \n",i+1,dnorm);
  AMG_Print(buffer);
  return(i+1);
}                  
/****************************************************************************/
/*                                                                            */
/* BiConjugate Gradient Stabilized                                          */
/*                                                                            */
/****************************************************************************/
int prepare_bcgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;

  w = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
  r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"r");
  p[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"p");
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
      r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
      p[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"p");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
    }
  return(0);
}
int clear_bcgs_solve(AMG_SolverContext *sc,int depth)
{
  int k;

  free(w->x);
  free(w);
  free(z[0]->x);
  free(z[0]);
  free(d[0]->x);
  free(d[0]);
  free(r[0]->x);
  free(r[0]);
  free(p[0]->x);
  free(p[0]);
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(z[k]->x);
      free(z[k]);
      free(d[k]->x);
      free(d[k]);
      free(r[k]->x);
      free(r[k]);
      free(p[k]->x);
      free(p[k]);
    }
  return(0);
}

int bcgs_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j;
  double dnorm, dnorm0, dnormlast;
  double rho,rho_last,alpha,beta,omega,delta;
  AMG_SolverContext *sc=global_sc;
  
                                                 /* iterate */
  AMG_dcopy(r[0],b);
  dmatminus(r[0],A[0],B,x);                    /* r[0] is initial residual */
  AMG_dcopy(b,r[0]);                                 /* b is initial residual, r^tilde in algorithm */

  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r[0],r[0]));
  start_residual=dnorm;
  residuals[0]=dnorm;
  residual_cnt++;
  
  if (sc->verbose>0)
    {
      sprintf(buffer,"Entering BCGS: %4d %12.4lE \n",0,dnorm);
      AMG_Print(buffer);
    }
  if (dnorm<sc->res_norm_min) 
    {
       iteration_cnt=0;
       end_residual=dnorm;
       return(0);
    }
  for (i=0; i<sc->maxit; i++)
    {
      rho=AMG_ddot(r[0],b);                       /* b is r[0]^tilde */
      if (fabs(rho)<1.0E-50)
        {
          AMG_Print("BCGS break down\n");
          iteration_cnt=i;
          end_residual=dnorm;                
          return(-1);
        }
      
      if (i==0)                                   /* compute p */
        AMG_dcopy(p[0],r[0]);
      else
        {
          beta=(rho/rho_last)*(alpha/omega);
          AMG_dscale(p[0],beta);
          AMG_daxpy(p[0],1.0,r[0]);
          AMG_daxpy(p[0],-beta*omega,w);
        }
      rho_last=rho;
      
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dset(z[0],0.0);
          AMG_dcopy(d[0],p[0]);
          for (j=0;j<sc->amg_prec_it;j++)
            {
              preconditioner(sc,0,depth,A,G,M,B,z,p,d);
              if (j<sc->amg_prec_it-1)
              {
                 AMG_dcopy(d[0],p[0]);
                 dmatminus(d[0],A[0],B,z[0]);
              }
            }
        }
      else
          AMG_dcopy(z[0],p[0]);        
      dmatmul(w,A[0],B,z[0]);
      delta=AMG_ddot(b,w);
      if (fabs(delta)<1.0E-50) delta=1.0E-50;
      alpha=rho/delta;
      AMG_daxpy(r[0],-alpha,w);
      AMG_daxpy(x,alpha,z[0]);
      dnorm=sqrt(AMG_ddot(r[0],r[0]));
      if ( (dnorm<dnorm0*sc->red_factor || dnorm<sc->res_norm_min) && !sc->ex_maxit )
        {
          if (sc->verbose>0)
            {
              sprintf(buffer,"BCGS Iteration %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
              AMG_Print(buffer);
            }
          break;
        }
      
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dset(z[0],0.0);
          AMG_dcopy(d[0],r[0]);
          for (j=0;j<sc->amg_prec_it;j++)
            {
              preconditioner(sc,0,depth,A,G,M,B,z,r,d);
              if (j<sc->amg_prec_it-1)
              {
                 AMG_dcopy(d[0],r[0]);
                 dmatminus(d[0],A[0],B,z[0]);
              }
            }
        }
      else
          AMG_dcopy(z[0],r[0]);
        
      dmatmul(d[0],A[0],B,z[0]);
      delta=AMG_ddot(d[0],d[0]);
      if (fabs(delta)<1.0E-50) delta=1.0E-50;
      omega=AMG_ddot(d[0],r[0])/delta;
      AMG_daxpy(r[0],-omega,d[0]);
      AMG_daxpy(x,omega,z[0]);
      dnorm=sqrt(AMG_ddot(r[0],r[0]));
      if (sc->verbose>0)
        {
          sprintf(buffer,"BCGS Iteration %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
          AMG_Print(buffer);
        }
      if (i<sc->maxit-1)
        {
          residuals[residual_cnt%AMG_CONV_RATE_BACK]=dnorm;
          residual_cnt++;
        }                
      finish = clock();
      if (finish>=start)
        elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
      else
        elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
      start = clock();
      dnormlast=dnorm;
      if (sc->ex_maxit) continue;
      if (dnorm<dnorm0*sc->red_factor) break;
      if (dnorm<sc->res_norm_min) break;
      if (dnorm>sc->div_factor*dnorm0)
        {
          AMG_Print("BCGS iteration diverges !!!\n");
          exit(4711);
        }
    }
  if (i==sc->maxit && !sc->ex_maxit && sc->verbose>=0)
    {
      sprintf(buffer,"BCGS : (maximal) iterations %d residual %g \n",i,dnorm);
      AMG_Print(buffer);
      iteration_cnt=i+1;
      end_residual=dnorm;              
      return(-1);
    }
  if (sc->verbose>=0)
    {
      sprintf(buffer,"BCGS : iterations %d residual %g \n",i+1,dnorm);
      AMG_Print(buffer);
    }
  iteration_cnt=i+1;
  end_residual=dnorm;                
  return(i+1);
}                        

/****************************************************************************/
/*                                                                            */
/* Mixed-BiCGStab-CGS                                                       */
/*                                                                            */
/****************************************************************************/
int prepare_mixed_bcgs_cgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;

  w = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  w1 = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  u = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  v1 = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  q = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  tilde_r = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"w");
  z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
  r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"r");
  p[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"p");
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  t[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"t");
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
      r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
      p[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"p");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
      t[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"t");
    }
  return(0);
}

int clear_mixed_bcgs_cgs_solve(AMG_SolverContext *sc, int depth)
{
  int k;

  free(w->x);
  free(w);
  free(w1->x);
  free(w1);
  free(u->x);
  free(u);
  free(v1->x);
  free(v1);
  free(q->x);
  free(q);
  free(tilde_r->x);
  free(tilde_r);
  free(z[0]->x);
  free(z[0]);
  free(d[0]->x);
  free(d[0]);
  free(r[0]->x);
  free(r[0]);
  free(p[0]->x);
  free(p[0]);
  free(t[0]->x);
  free(t[0]);
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(z[k]->x);
      free(z[k]);
      free(d[k]->x);
      free(d[k]);
      free(r[k]->x);
      free(r[k]);
      free(p[k]->x);
      free(p[k]);
      free(t[k]->x);
      free(t[k]);
    }
  return(0);
}

int mixed_bcgs_cgs_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,k,cgs;
  double dnorm, dnorm0, dnormlast,tol;
  double rho,rho_new,alpha[max_switch],beta[max_switch],omega,delta,tmp;
  AMG_SolverContext *sc=global_sc;

  tol = sc->mixed_bcgs_cgs_switch_tol;
 
  AMG_dcopy(r[0],b);                           /* copy rhs to r[0] */
  dmatminus(r[0],A[0],B,x);                    /* r[0] is initial residual */
  AMG_dcopy(tilde_r,r[0]);                               /* b is initial residual, r^tilde in algorithm */
  /*  AMG_dscale(tilde_r,-1.0);        */                       
  AMG_dcopy(u,r[0]);                               /* initialize */
  AMG_dcopy(v1,r[0]);                               
  AMG_dcopy(p[0],r[0]);                               
  rho = AMG_ddot(tilde_r,r[0]);
  k=0;

  for (i=0;i<max_switch;i++) alpha[i]=beta[i]=0.0;

  dnorm=dnorm0=dnormlast=sqrt(AMG_ddot(r[0],r[0]));
  start_residual=dnorm;
  residuals[0]=dnorm;
  residual_cnt++;
  
  if (sc->verbose>0)
    {
      sprintf(buffer,"Entering Mixed-BiCGStab-CGS: %4d %12.4lE \n",0,dnorm);
      AMG_Print(buffer);
    }
  if (dnorm<sc->res_norm_min) 
    {
       iteration_cnt=0;
       end_residual=dnorm;
       return(0);
    }
  cgs=1;
  for (i=0; i<sc->maxit; i++)
    {
      if (cgs)
        {
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(d[0],p[0]);                     /* rhs is p, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,p,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],p[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],p[0]);
          dmatmul(w,A[0],B,z[0]);                   /* matrix times preconditioner */
          tmp =  AMG_ddot(tilde_r,w);
          for (j=1;j<=k;j++)
            alpha[j]=alpha[j-1];
          alpha[0] = rho/tmp;
          AMG_daxpby(q,1.0,v1,-alpha[0],w);
          AMG_daxpby(t[0],alpha[0],u,alpha[k],q);
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(d[0],t[0]);                     /* rhs is t, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,t,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],t[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],t[0]);            
          AMG_daxpy(x,1.0,z[0]);
          dmatmul(w1,A[0],B,z[0]);                   /* matrix times preconditioner */
          AMG_daxpy(r[0],-1.0,w1);
          dnorm=sqrt(AMG_ddot(r[0],r[0]));
          if (sc->verbose>0)
            {
              sprintf(buffer,"CGS Iteration %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
              AMG_Print(buffer);
            }
          if (i<sc->maxit-1)
            {
              residuals[residual_cnt%AMG_CONV_RATE_BACK]=dnorm;
              residual_cnt++;
            }                
          finish = clock();
          if (finish>=start)
            elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
          else
            elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
          start = clock();
          if ((dnorm/dnormlast > tol)&& (k<max_switch)) cgs=0;      /* switch to bcgs */
          dnormlast=dnorm;
          if (dnorm<dnorm0*sc->red_factor) break;
          if (dnorm<sc->res_norm_min) break;
          if (dnorm>sc->div_factor*dnorm0)
            {
              AMG_Print("CGS/BCGS iteration diverges !!!\n");
              exit(4711);
            }
          rho_new = AMG_ddot(tilde_r,r[0]);
          for (j=1;j<=k;j++)
            beta[j]=beta[j-1];
          beta[0] = alpha[0]*rho_new/(alpha[k]*rho);
          rho = rho_new;
          AMG_daxpby(u,1.0,u,-alpha[k],w);
          AMG_daxpby(u,1.0,r[0],beta[0],u);
          AMG_daxpby(v1,1.0,r[0],beta[k],q);
          AMG_daxpby(p[0],1.0,q,beta[0],p[0]);
          AMG_daxpby(p[0],1.0,u,beta[k],p[0]);
        }
      else
        {                                           /* start like in cgs */
          k++;
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(t[0],u);                     /* rhs is p, defect is d */   
              AMG_dcopy(d[0],t[0]);                     /* rhs is p, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,t,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],t[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],u);             
          dmatmul(w,A[0],B,z[0]);                   /* matrix times preconditioner */
          tmp =  AMG_ddot(tilde_r,w);
          alpha[0] = rho/tmp;
          AMG_daxpby(t[0],1.0,r[0],-alpha[0],w);
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(d[0],t[0]);                     /* rhs is t, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,t,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],t[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],t[0]);            
          dmatmul(w1,A[0],B,z[0]);                   /* matrix times preconditioner */
          omega = AMG_ddot(t[0],w1)/AMG_ddot(w1,w1);
          AMG_daxpby(r[0],1.0,t[0],-omega,w1);
          AMG_daxpby(t[0],alpha[0],u,omega,t[0]);
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(d[0],t[0]);                     /* rhs is t, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,t,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],t[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],t[0]);            
          AMG_daxpy(x,1.0,z[0]);
          dnorm=sqrt(AMG_ddot(r[0],r[0]));
          if (sc->verbose>0)
            {
              sprintf(buffer,"BCGS Iteration %4d %12.4lE %12.4lg\n",i+1,dnorm,dnorm/dnormlast);
              AMG_Print(buffer);
            }
          if (i<sc->maxit-1)
            {
              residuals[residual_cnt%AMG_CONV_RATE_BACK]=dnorm;
              residual_cnt++;        
            }        
          finish = clock();
          if (finish>=start)
            elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
          else
            elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
          start = clock();
          dnormlast=dnorm;
          if (sc->ex_maxit) continue;
          if (dnorm<dnorm0*sc->red_factor) break;
          if (dnorm<sc->res_norm_min) break;
          if (dnorm>sc->div_factor*dnorm0)
            {
              AMG_Print("CGS/BCGS iteration diverges !!!\n");
              exit(4711);
            }
          rho_new = AMG_ddot(tilde_r,r[0]);
          beta[0] = alpha[0]*rho_new/(omega*rho);
          rho = rho_new;
          AMG_daxpy(u,-omega,w);
          AMG_daxpby(u,1.0,r[0],beta[0],u);
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(d[0],p[0]);                     /* rhs is p, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,p,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],p[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],p[0]); 
          dmatmul(w,A[0],B,z[0]);                   /* matrix times preconditioner */
          AMG_daxpy(v1,-alpha[0],w);
          if (sc->amg_prec_it)                      /* preconditioner should be applied */
            {
              AMG_dset(z[0],0.0);                       /* prepare preconditioner */
              AMG_dcopy(t[0],v1);                     /* rhs is p, defect is d */   
              AMG_dcopy(d[0],t[0]);                     /* rhs is p, defect is d */   
              for (j=0;j<sc->amg_prec_it;j++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,z,t,d);   /* solution is z */
                  if (j<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],t[0]);                       
                     dmatminus(d[0],A[0],B,z[0]);                /* compute new defect */
                  }
                }
            }
          else
            AMG_dcopy(z[0],v1); 
          dmatmul(w1,A[0],B,z[0]);                   /* matrix times preconditioner */
          AMG_daxpy(v1,-omega,w1);
          AMG_daxpy(p[0],-omega,w);
          AMG_daxpby(p[0],1.0,v1,beta[0],p[0]);          
          cgs = 1;
               }
    }
   if (i==sc->maxit && !sc->ex_maxit)
    {
      AMG_Print("solver not converged\n");
      sprintf(buffer,"MIXED_BCGS_CGS : (maximal) iterations %d residual %g \n",i+1,dnorm);
      AMG_Print(buffer);
      iteration_cnt=i;
      end_residual=dnorm;              
      return(-1);
    }
  sprintf(buffer,"MIXED_BCGS_CGS : iterations %d residual %g \n",i+1,dnorm);
  AMG_Print(buffer);
  iteration_cnt=i;
  end_residual=dnorm;                
  
  return(i+1);
}
/****************************************************************************/
/*                                                                          */
/* Generalized Minimum RESidual, left preconditioned                        */
/*                                                                          */
/****************************************************************************/

void AMG_GeneratePlaneRotation(double dx, double dy, double *cs, double *sn)
{
  if (dy == 0.0) 
    {
      cs[0] = 1.0;
      sn[0] = 0.0;
    } 
  else if (fabs(dy) > fabs(dx)) 
    {
      double temp = dx / dy;
      sn[0] = 1.0 / sqrt( 1.0 + temp*temp );
      cs[0] = temp * sn[0];
    } 
  else 
    {
      double temp = dy / dx;
      cs[0] = 1.0 / sqrt( 1.0 + temp*temp );
      sn[0] = temp * cs[0];
    }
}

void AMG_ApplyPlaneRotation(double *dx, double *dy, double cs, double sn)
{
  double temp  =  cs * dx[0] + sn * dy[0];
  dy[0] = -sn * dx[0] + cs * dy[0];
  dx[0] = temp;
}

void AMG_Update(AMG_VECTOR *x, int Len_x, int k, AMG_VECTOR **h, 
                AMG_VECTOR *s, AMG_VECTOR **v)
{
  int i,j;
  /* Backsolve: */ 
  
  for (i = k; i >= 0; i--) 
    {
      s->x[i] /= h[i]->x[i];
      for (j = i - 1; j >= 0; j--)
        s->x[j] -= h[j]->x[i] * s->x[i];
    }
  
  for (j = 0; j <= k; j++)
    for(i=0; i<Len_x; i++)
      x->x[i] += v[j]->x[i] * s->x[j];
}


int prepare_gmres_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;

  if (sc->gmres_restart>sc->maxit) sc->gmres_restart=sc->maxit; 
  
  if (sc->gmres_restart>AMG_MAX_GMRES_RESTART)
    {
      sprintf(buffer,"WARNING: restart %d greater than AMG_MAX_GMRES_RESTART %d!\n",
              sc->gmres_restart,AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sprintf(buffer,"WARNING: set restart to %d!\n",AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sc->gmres_restart=AMG_MAX_GMRES_RESTART;
    }
                                             /* allocate vectors depending on restart */
  s = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"s");   
  cosi = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"cs"); 
  sn = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"sn");
        
  for (k=0;k<=sc->gmres_restart;k++)    /* allocate matrices depending on restart */
    {
      H[k] = AMG_NewVector(sc->gmres_restart,AMG_MATRIX_B(A[0]),"H");
      AMG_dset(H[k],0.0);
      v[k] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"v");
      AMG_dset(v[k],0.0);      
    }
   if ((sc->step_length_control_fine)||(sc->step_length_control_all)) 
    {
      AMG_Print("MESSAGE : Static preconditioned GMRES with step length control\n");
      AMG_Print("MESSAGE :    for the preconditioner will in general not converge !\n");
      AMG_Print("MESSAGE :    Use instead flexible GMRES (AMG_GMRES_FLEX) !!!\n");
      AMG_Print("MESSAGE :    Step length control switched off !\n");
      sprintf(buffer,"MESSAGE :    Static damping with %g used instead !\n",sc->omega[0]);
      AMG_Print(buffer);
      sc->step_length_control_fine=0;
      sc->step_length_control_all=0;
    }
   r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"r");
   d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
   z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
   if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
   for (k=1; k<=depth; k++)
    {
      r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
      z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
    }
  return(0);
}

int clear_gmres_solve(AMG_SolverContext *sc,int depth)
{
  int k;

  free(s->x);
  free(s);
  free(cosi->x);
  free(cosi);
  free(sn->x);
  free(sn);        
  for (k=0;k<=sc->gmres_restart;k++)    /* free matrices depending on restart */
    {
      free(H[k]->x);
      free(H[k]);
      free(v[k]->x);
      free(v[k]);
    }
  free(r[0]->x);
  free(r[0]);
  free(d[0]->x);
  free(d[0]);
  free(z[0]->x);
  free(z[0]);
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(r[k]->x);
      free(r[k]);
      free(d[k]->x);
      free(d[k]);
      free(z[k]->x);
      free(z[k]);
    }
  return(0);
}

int gmres_left_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,k,l;
  double beta,resid,residlast;
  AMG_SolverContext *sc=global_sc;
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering GMRES\n");
      AMG_Print(buffer);
    }
  AMG_dcopy(z[0],b);                     /* copy rhs (b) into r */
  

  dmatminus(z[0],A[0],B,x);               /* r = r - A*x, r is overwritten */

  if (sc->amg_prec_it)                      /* preconditioner should be applied */
    {
      AMG_dcopy(d[0],z[0]);                 /* starting defect stored in d */
      AMG_dset(r[0],0.0);                   /* set solution to zero */
      for (l=0;l<sc->amg_prec_it;l++)       /* number of preconditioner iterations */
        {                                   /* apply preconditioner */
          preconditioner(sc,0,depth,A,G,M,B,r,z,d);
	  if (l<sc->amg_prec_it-1)
	  {
	    AMG_dcopy(d[0],z[0]);             /* copy rhs to array d (defect) */
	    dmatminus(d[0],A[0],B,r[0]);    /* compute new defect */
	  }
        }
    } 
  else                                      /* no preconditioning */
    AMG_dcopy(r[0],z[0]);
  
  beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
  resid = beta;
  start_residual=resid;
  residuals[0]=resid;
  residual_cnt++;
  
  if(beta  <= sc->res_norm_min )           /* stopping criterion fulfilled */
    {
      if ((sc->verbose>0)&&(sc->minit==0)) 
        {
          sprintf(buffer,"GMRES Iteration %4d %12.4lE\n",0,beta);
          AMG_Print(buffer);
        }
      iteration_cnt=0;
      end_residual=resid;
      return(0);
    }
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering GMRES iteration cycle\n");
      AMG_Print(buffer);
    }
  j=1;
  while (j<=sc->maxit)
    {
      AMG_dcopy(v[0],r[0]);
      AMG_dscale(v[0],1.0/beta);
      AMG_dset(s,0.0);
      s->x[0]=beta;
      
      for(i=0; i<sc->gmres_restart && j<=sc->maxit; i++, j++)
        {
          dmatmul(z[0],A[0],B,v[i]);               /* d=A*v[i]; */
          
          if (sc->amg_prec_it)
            {
              AMG_dcopy(r[0],z[0]);
              AMG_dset(d[0],0.0);                       /* d = sol, z=rhs ,r=def */
              for (k=0;k<sc->amg_prec_it;k++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,d,z,r);
		  if (k<sc->amg_prec_it-1)
		  {
		    AMG_dcopy(r[0],z[0]);
		    dmatminus(r[0],A[0],B,d[0]);
		  }
                }
            }
          else
              AMG_dcopy(d[0],z[0]);
            
          for(k=0; k<=i; k++)
            {
              H[k]->x[i]=AMG_ddot(d[0], v[k]);
              AMG_daxpy(d[0],-H[k]->x[i], v[k]);
            }                                          
          
          H[i+1]->x[i]=sqrt(AMG_ddot(d[0],d[0]));
          AMG_dcopy(v[i+1],d[0]);
          AMG_dscale(v[i+1], 1.0/(H[i+1]->x[i]));
          
          for(k=0;k<i;k++)
            AMG_ApplyPlaneRotation(&H[k]->x[i], &H[k+1]->x[i], cosi->x[k], sn->x[k]);
          
          AMG_GeneratePlaneRotation(H[i]->x[i], H[i+1]->x[i], &cosi->x[i], &sn->x[i]);
          AMG_ApplyPlaneRotation(&H[i]->x[i], &H[i+1]->x[i], cosi->x[i], sn->x[i]);
          AMG_ApplyPlaneRotation(&s->x[i], &s->x[i+1], cosi->x[i], sn->x[i]);
          
          residlast=resid;
          if (((resid=fabs(s->x[i+1])) < sc->res_norm_min)
            || ((resid=fabs(s->x[i+1]))<sc->red_factor*start_residual))
            {
              AMG_Update(x,AMG_VECTOR_N(x), i, H, s, v);
              if (sc->verbose>0)
                {
                  sprintf(buffer,"GMRES (left) : iterations %d residual %g \n",j,resid);
                  AMG_Print(buffer);
                }
              iteration_cnt=j;
              end_residual=resid;              
              return(j);
            }
          if (sc->verbose>0)
            {
              sprintf(buffer,"GMRES (left) Iteration %4d %12.4lE %12.4lg\n",j,resid,resid/residlast);
                  AMG_Print(buffer);
            }
          if (resid>sc->div_factor*start_residual)
            {
              AMG_Print("GMRES (left) iteration diverges !!!\n");
              exit(4711);
            }
          if (j<sc->maxit)
            {
              residuals[residual_cnt%AMG_CONV_RATE_BACK]=resid;
              residual_cnt++;
            }          
          finish = clock();
          if (finish>=start)
            elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
          else
            elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
          start = clock();
        }                                                   /* endfor i */
      
      AMG_Update(x,AMG_VECTOR_N(x), sc->gmres_restart-1, H, s, v);
      AMG_dcopy(z[0],b);                        /* copy rhs (b) into r */
      dmatminus(z[0],A[0],B,x);               /* r = r - A*x, r is overwritten */
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dcopy(d[0],z[0]);                 /* starting defect stored in d */
          AMG_dset(r[0],0.0);                   /* set solution to zero */
          for (l=0;l<sc->amg_prec_it;l++)       /* number of preconditioner iterations */
            {                                   /* apply preconditioner */
              preconditioner(sc,0,depth,A,G,M,B,r,z,d);
	      if (l<sc->amg_prec_it-1)
	      {
		AMG_dcopy(d[0],z[0]);             /* copy rhs to array d (defect) */
		dmatminus(d[0],A[0],B,r[0]);    /* compute new defect */
	      }
            }
        }
      else
          AMG_dcopy(r[0],z[0]);                 /* no preconditioning */
        
      beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
      
      if (sc->verbose>1)
        {    
          if(fabs(beta-resid)>0.01*beta) 
            printf("restart residual changed  %g %g\n",sqrt(AMG_ddot(r[0],r[0])),resid);
        }
      if((resid=beta) < sc->res_norm_min)
        {
          if (sc->verbose>0)
            {
              sprintf(buffer,"GMRES (left) : iterations %d residual %g \n",j,resid);
              AMG_Print(buffer);
            }
          return(j);
        }
    }
  /* endwhile */
  if (sc->verbose>0)
    {
      sprintf(buffer,"GMRES (left) : (maximal) iterations %d residual %g \n",sc->maxit,resid);
      AMG_Print(buffer);
    }
  iteration_cnt=j-1;
  end_residual=resid;              
  return(j);
}
/****************************************************************************/
/*                                                                            */
/* Generalized Minimum RESidual, right preconditioned                       */
/*                                                                            */
/****************************************************************************/

int gmres_right_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,k,l,count_restarts=0;
  double beta,resid,residlast;
  AMG_SolverContext *sc=global_sc;
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering GMRES\n");
      AMG_Print(buffer);
    }
  AMG_dcopy(r[0],b);                     /* copy rhs (b) into r */
  dmatminus(r[0],A[0],B,x);               /* r = r - A*x, r is overwritten */
  
  beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
  resid = beta;
  start_residual=resid;
  residuals[0]=resid;
  residual_cnt++;
  
  if ((beta  <= sc->res_norm_min )&&(sc->minit==0))           /* stopping criterion fulfilled */
    {
      if (sc->verbose>0)
        {
          sprintf(buffer,"GMRES Iteration %4d %12.4lE\n",0,beta);
          AMG_Print(buffer);
          iteration_cnt=0;
          end_residual=resid;
         }
      return(0);
    }
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering GMRES iteration cycle\n");
      AMG_Print(buffer);
    }
  j=1;
  while (j<=sc->maxit)
    {
      AMG_dcopy(v[0],r[0]);
      AMG_dscale(v[0],1.0/beta);
      AMG_dset(s,0.0);
      s->x[0]=beta;
      
      for(i=0; i<sc->gmres_restart && j<=sc->maxit; i++, j++)
        {
          if (sc->amg_prec_it)
            {
              AMG_dcopy(z[0],v[i]);
              AMG_dcopy(d[0],v[i]);
              AMG_dset(r[0],0.0);                       /* r = sol, z=rhs ,d=def */
              for (k=0;k<sc->amg_prec_it;k++)
                {
                  preconditioner(sc,0,depth,A,G,M,B,r,z,d);
                  if (k<sc->amg_prec_it-1)
                  {
                     AMG_dcopy(d[0],z[0]);
                     dmatminus(d[0],A[0],B,r[0]);
                  }        
                }    
              dmatmul(d[0],A[0],B,r[0]);               /* d=A*v[i]; */
            }
          else
            dmatmul(d[0],A[0],B,v[i]);                 /* d=A*v[i]; */
            
          for(k=0; k<=i; k++)
            {
              H[k]->x[i]=AMG_ddot(d[0], v[k]);
              AMG_daxpy(d[0],-H[k]->x[i], v[k]);
            }                                          
          
          H[i+1]->x[i]=sqrt(AMG_ddot(d[0],d[0]));
          AMG_dcopy(v[i+1],d[0]);
          AMG_dscale(v[i+1], 1.0/(H[i+1]->x[i]));
          
          for(k=0;k<i;k++)
            AMG_ApplyPlaneRotation(&H[k]->x[i], &H[k+1]->x[i], cosi->x[k], sn->x[k]);
          
          AMG_GeneratePlaneRotation(H[i]->x[i], H[i+1]->x[i], &cosi->x[i], &sn->x[i]);
          AMG_ApplyPlaneRotation(&H[i]->x[i], &H[i+1]->x[i], cosi->x[i], sn->x[i]);
          AMG_ApplyPlaneRotation(&s->x[i], &s->x[i+1], cosi->x[i], sn->x[i]);
          
          residlast=resid;
          if (((resid=fabs(s->x[i+1])) < sc->res_norm_min)
            || ((resid=fabs(s->x[i+1]))<sc->red_factor*start_residual))
            {
              AMG_dset(z[0],0.0);                       /* set solution to zero */      
              AMG_Update(z[0],AMG_VECTOR_N(z[0]), i, H, s, v);
              if (sc->amg_prec_it)                      /* preconditioner should be applied */
                {
                  AMG_dcopy(d[0],z[0]);                 /* starting defect stored in d */
                  AMG_dset(r[0],0.0);                   /* set solution to zero */
                  for (l=0;l<sc->amg_prec_it;l++)       /* number of preconditioner iterations */
                    {                                   /* apply preconditioner */
                      preconditioner(sc,0,depth,A,G,M,B,r,z,d);
                      if (l<sc->amg_prec_it-1)
                      {
                         AMG_dcopy(d[0],z[0]);             /* copy rhs to array d (defect) */
                         dmatminus(d[0],A[0],B,r[0]);    /* compute new defect */
                      }
                    }
                  AMG_daxpy(x,1.0,r[0]);                    /* new iterate */
                }
              else
                AMG_daxpy(x,1.0,z[0]); 
              if (sc->verbose>0)
                {
                  sprintf(buffer,"GMRES (right) : iterations %d residual %g \n",j,resid);
                  AMG_Print(buffer);
                }
              iteration_cnt=j;
              end_residual=resid;              
              return(j);
            }
          if (sc->verbose>0)
            {
              sprintf(buffer,"GMRES (right) Iteration %4d %12.4lE %12.4lg\n",j,resid,resid/residlast);
                  AMG_Print(buffer);
            }
          if (resid>sc->div_factor*start_residual)
            {
              AMG_Print("GMRES (right) iteration diverges !!!\n");
              exit(4711);
            }
          if (j<sc->maxit)
            {
              residuals[residual_cnt%AMG_CONV_RATE_BACK]=resid;
              residual_cnt++;
            }
          finish = clock();
          if (finish>=start)
            elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
          else
            elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
          start = clock();
        }                                                   /* endfor i */
                                                /* compute the update=rhs for prec. */
      AMG_dset(z[0],0.0);                       /* set solution to zero */      
      AMG_Update(z[0],AMG_VECTOR_N(z[0]), sc->gmres_restart-1, H, s, v);
      if (sc->amg_prec_it)                      /* preconditioner should be applied */
        {
          AMG_dcopy(d[0],z[0]);                 /* starting defect stored in d */
          AMG_dset(r[0],0.0);                   /* set solution to zero */
          for (l=0;l<sc->amg_prec_it;l++)       /* number of preconditioner iterations */
            {                                   /* apply preconditioner */
              preconditioner(sc,0,depth,A,G,M,B,r,z,d);
              if (l<sc->amg_prec_it-1)
              {
                 AMG_dcopy(d[0],z[0]);             /* copy rhs to array d (defect) */
                 dmatminus(d[0],A[0],B,r[0]);    /* compute new defect */
              }
            }
        }
      else
        AMG_dcopy(r[0],z[0]);
      count_restarts++;
      AMG_daxpy(x,1.0,r[0]);                    /* new iterate */
      AMG_dcopy(r[0],b);                        /* copy rhs (b) into r */
      dmatminus(r[0],A[0],B,x);               /* new defect r = r - A*x, r is overwritten */

          

      beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
      
      if (sc->verbose>1)
        {    
          if(fabs(beta-resid)>0.01*beta) 
            printf("restart residual changed  %g %g\n",sqrt(AMG_ddot(r[0],r[0])),resid);
        }
    }
  /* endwhile */
  if (sc->verbose>0)
    {
      sprintf(buffer,"GMRES (right) : (maximal) iterations %d residual %g \n",sc->maxit,resid);
      AMG_Print(buffer);
    }
  iteration_cnt=j-1;
  end_residual=resid;              
  return(j);
}

/****************************************************************************/
/*                                                                          */
/* Generalized Minimum RESidual, flexible preconditioner                    */
/*                                                                          */
/****************************************************************************/

int prepare_gmres_flexible_solve(AMG_SolverContext *sc, int *A_length, 
                                 int *B_length, int depth)
{
  int k;

  if (sc->gmres_restart>sc->maxit) sc->gmres_restart=sc->maxit; 
  if (sc->gmres_restart>AMG_MAX_GMRES_RESTART)
    {
      sprintf(buffer,"MESSAGE : restart %d greater than AMG_MAX_GMRES_RESTART %d!\n",
              sc->gmres_restart,AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sprintf(buffer,"MESSAGE : set restart to %d!\n",AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sc->gmres_restart=AMG_MAX_GMRES_RESTART;
    }
  /* allocate vectors depending on restart */
  s = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"s");   
  cosi = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"cs"); 
  sn = AMG_NewVector(sc->gmres_restart+1,AMG_MATRIX_B(A[0]),"sn");
  
  for (k=0;k<=sc->gmres_restart;k++)    /* allocate matrices depending on restart */
    {
      H[k] = AMG_NewVector(sc->gmres_restart,AMG_MATRIX_B(A[0]),"H");
      AMG_dset(H[k],0.0);
      v[k] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"v");
      AMG_dset(v[k],0.0);  
      zv[k] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"zv");
      AMG_dset(zv[k],0.0);  
    }
  r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"r");
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  
  for (k=1; k<=depth; k++)
    {
      r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
      z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
    }
  return(0);
}

int clear_gmres_flexible_solve(AMG_SolverContext *sc, int depth)
{
  int k;

  free(s->x);
  free(s);
  free(cosi->x);
  free(cosi);
  free(sn->x);
  free(sn);      
  for (k=0;k<=sc->gmres_restart;k++)    /* free matrices depending on restart */
    {
      free(H[k]->x);
      free(H[k]);
      free(v[k]->x);
      free(v[k]);
      free(zv[k]->x);
      free(zv[k]);
    }
  free(r[0]->x);
  free(r[0]);
  free(d[0]->x);
  free(d[0]);
  free(z[0]->x);
  free(z[0]);
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(r[k]->x);
      free(r[k]);
      free(d[k]->x);
      free(d[k]);
      free(z[k]->x);
      free(z[k]);
    }
  return(0);
}

int gmres_flexible_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,k,l;
  double beta,resid,residlast,dnorm0;
  AMG_SolverContext *sc=global_sc;
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering FGMRES\n");
      AMG_Print(buffer);
    }
  if (sc->amg_prec_it<1)
    {
      sprintf(buffer,"MESSAGE : Number of preconditioner iterations too small: %d\n",sc->amg_prec_it<1);
      AMG_Print(buffer);
      sprintf(buffer,"          Set number of preconditioner iterations to 1 !!!");
      AMG_Print(buffer);
    }
  AMG_dcopy(r[0],b);                     /* copy rhs (b) into r */
  dmatminus(r[0],A[0],B,x);               /* r = r - A*x, r is overwritten */
  
  beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
  resid = beta;
  start_residual=resid;
  residuals[0]=resid;
  residual_cnt++;

  if (((beta  <= sc->res_norm_min ) &&(sc->minit==0))  
      || (beta<1e-20))                 /* stopping criterion fulfilled */
    {
      if (sc->verbose>0)
        {
          sprintf(buffer,"FGMRES Iteration %4d %12.4lE (no iteration)\n",0,beta);
          AMG_Print(buffer);
        }
      iteration_cnt=0;
      end_residual=resid;
      return(0);
    }
  
  if (sc->verbose>0)
    {
      sprintf(buffer,"FGMRES Iteration %4d %12.4lE\n",0,beta);
      AMG_Print(buffer);
        }
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering FGMRES iteration cycle\n");
      AMG_Print(buffer);
    }
  j=1;
  while (j<=sc->maxit)
    {
      AMG_dcopy(v[0],r[0]);
      AMG_dscale(v[0],1.0/beta);
      AMG_dset(s,0.0);
      s->x[0]=beta;
      
      for(i=0; i<sc->gmres_restart && j<=sc->maxit; i++, j++)
        {
          AMG_dcopy(z[0],v[i]);
          AMG_dcopy(d[0],v[i]);
          if (sc->amg_prec_it>1)
            dnorm0 = sqrt(AMG_ddot(d[0],d[0]))*sc->amg_prec_red_factor;
          AMG_dset(r[0],0.0);                       /* r = sol, z=rhs ,d=def */
          for (k=0;k<sc->amg_prec_it;k++)
            {
              preconditioner(sc,0,depth,A,G,M,B,r,z,d);
              if (k<sc->amg_prec_it-1)
              {
                 AMG_dcopy(d[0],z[0]);
                 dmatminus(d[0],A[0],B,r[0]);
                 if (sqrt(AMG_ddot(d[0],d[0]))<= dnorm0)
                    break;
              }
            }
         
          AMG_dcopy(zv[i],r[0]);
          dmatmul(d[0],A[0],B,r[0]);               /* d=A*v[i]; */

          for(k=0; k<=i; k++)
            {
              H[k]->x[i]=AMG_ddot(d[0], v[k]);
              AMG_daxpy(d[0],-H[k]->x[i], v[k]);
            }                                          
          
          H[i+1]->x[i]=sqrt(AMG_ddot(d[0],d[0]));
          AMG_dcopy(v[i+1],d[0]);
          AMG_dscale(v[i+1], 1.0/(H[i+1]->x[i]));
          
          for(k=0;k<i;k++)
            AMG_ApplyPlaneRotation(&H[k]->x[i], &H[k+1]->x[i], cosi->x[k], sn->x[k]);
          
          AMG_GeneratePlaneRotation(H[i]->x[i], H[i+1]->x[i], &cosi->x[i], &sn->x[i]);
          AMG_ApplyPlaneRotation(&H[i]->x[i], &H[i+1]->x[i], cosi->x[i], sn->x[i]);
          AMG_ApplyPlaneRotation(&s->x[i], &s->x[i+1], cosi->x[i], sn->x[i]);
          
          residlast=resid;
          if (((resid=fabs(s->x[i+1])) < sc->res_norm_min)
            || ((resid=fabs(s->x[i+1]))<sc->red_factor*start_residual))
            {
              AMG_dset(r[0],0.0);                       /* set solution to zero */      
              AMG_Update(r[0],AMG_VECTOR_N(r[0]), i, H, s, zv);

              AMG_daxpy(x,1.0,r[0]);                    /* new iterate */
              if (sc->verbose>=1)
                {
                  sprintf(buffer,"FGMRES : iterations %d residual %g \n",j,resid);
                  AMG_Print(buffer);
                }
              iteration_cnt=j;
              end_residual=resid;              
              return(j);
            }
          if (sc->verbose>0)
            {
              sprintf(buffer,"FGMRES Iteration %4d %12.4lE %12.4lg\n",j,resid,resid/residlast);
              AMG_Print(buffer);
            }
          if (resid>sc->div_factor*start_residual)
            {
              AMG_Print("FGMRES iteration diverges !!!\n");
              exit(4711);
            }
          if (j<sc->maxit)
            {
              residuals[residual_cnt%AMG_CONV_RATE_BACK]=resid;
              residual_cnt++;
            }          
          finish = clock();
          if (finish>=start)
            elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
          else
            elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
          start = clock();
        }                                                   /* endfor i */
                                                /* compute the update=rhs for prec. */
      AMG_dset(r[0],0.0);                       /* set solution to zero */      
      AMG_Update(r[0],AMG_VECTOR_N(r[0]), sc->gmres_restart-1, H, s, zv);

      AMG_daxpy(x,1.0,r[0]);                    /* new iterate */
      AMG_dcopy(r[0],b);                        /* copy rhs (b) into r */
      dmatminus(r[0],A[0],B,x);               /* new defect r = r - A*x, r is overwritten */
      beta=sqrt(AMG_ddot(r[0],r[0]));           /* norm of residual */
      
      if (sc->verbose>1)
        {    
          if(fabs(beta-resid)>0.01*beta) 
            printf("restart residual changed  %g %g\n",sqrt(AMG_ddot(r[0],r[0])),resid);
        }
    }
  /* endwhile */
  if (sc->verbose>=0)
    {
      sprintf(buffer,"FGMRES : (maximal) iterations %d residual %g \n",sc->maxit,resid);
      AMG_Print(buffer);
    }
  iteration_cnt=j-1;
  end_residual=resid;              
  return(j);
}

/****************************************************************************/
/*                                                                            */
/* exact solver with LU decomposition                                       */
/*                                                                            */
/****************************************************************************/
int prepare_exact_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{  
  d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
  M[0]=prepare_ex(A[0]);
  if (M[0]==NULL) 
    {
      AMG_Print("Error in preparing exact solver !!!\n");
      AMG_Print("Use a different solver !!!\n");
      exit(1);
    }
  return(0);
}
int clear_exact_solve(AMG_SolverContext *sc)
{
  free(d[0]->x);
  free(d[0]);     
  free(r[0]->x);
  free(r[0]);     
  free(z[0]->x);
  free(z[0]);     
  free(M[0]->ra);
  free(M[0]->ja);
  free(M[0]->a);
  free(M[0]);
  return(0);
}

int exact_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  AMG_SolverContext *sc=global_sc;
  int i;

  AMG_dset(r[0],0.0);                        /* copy rhs (b) into r */
  AMG_dcopy(z[0],b);                        /* copy rhs (b) into r */
  ex(sc,0,0,A,G,M,B,r,d,z);
  AMG_dcopy(x,r[0]);                        /* copy rhs (b) into r */
  AMG_dcopy(z[0],b);                        /* copy rhs (b) into r */
  dmatminus(z[0],A[0],B,x);               /* new defect r = r - A*x, r is overwritten */
  sprintf(buffer,"EXACT SOLVER  residual %12.4lg\n",sqrt(AMG_ddot(z[0],z[0])));
  AMG_Print(buffer);

  return(1);
}

/****************************************************************************/
/*                                                                          */
/* LCD with restart                                                         */
/* Yuan, Golub, Plemmons, Cecilio, BIT 44, 189-207, 2004                    */
/* Catabriga, Coutinho, Franca, IJNME 60, 1513-1534, 2004                   */
/*                                                                          */
/****************************************************************************/

int prepare_lcd_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth)
{
  int k;

  if (sc->gmres_restart>sc->maxit) sc->gmres_restart=sc->maxit; 
  
  if (sc->gmres_restart>AMG_MAX_GMRES_RESTART)
    {
      sprintf(buffer,"WARNING: restart %d greater than AMG_MAX_GMRES_RESTART %d!\n",
              sc->gmres_restart,AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sprintf(buffer,"WARNING: set restart to %d!\n",AMG_MAX_GMRES_RESTART);
      AMG_Print(buffer);
      sc->gmres_restart=AMG_MAX_GMRES_RESTART;
    }
        
  for (k=0;k<=sc->gmres_restart;k++)    /* allocate matrices depending on restart */
    {
      v[k] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"v");
      AMG_dset(v[k],0.0);      
      zv[k] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"zv");
      AMG_dset(zv[k],0.0);      
    }

   r[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"r");
   d[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"d");
   z[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"z");
   p[0] = AMG_NewVector(A_length[0]+B_length[0],AMG_MATRIX_B(A[0]),"p");
   if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
   for (k=1; k<=depth; k++)
    {
      r[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"r");
      d[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"d");
      z[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"z");
      p[k] = AMG_NewVector(A_length[k]+B_length[k],AMG_MATRIX_B(A[k]),"p");
    }
  return(0);
}

int clear_lcd_solve(AMG_SolverContext *sc,int depth)
{
  int k;

  for (k=0;k<=sc->gmres_restart;k++)    /* free matrices depending on restart */
    {
      free(v[k]->x);
      free(v[k]);
      free(zv[k]->x);
      free(zv[k]);
    }
  free(r[0]->x);
  free(r[0]);
  free(d[0]->x);
  free(d[0]);
  free(z[0]->x);
  free(z[0]);
  free(p[0]->x);
  free(p[0]);
  if (A[0]->ratr != NULL)
  {
      free(A[0]->ratr);
      free(A[0]->jatr);
      free(A[0]->postr);
      A[0]->ratr = NULL;
  }
  if (M[0]->ratr != NULL)
  {
      free(M[0]->ratr);
      free(M[0]->jatr);
      free(M[0]->postr);
  }
  if (sc->preconditioner!=AMG_MGC) return(0); /* to avoid unnecessary allocation for schur complement methods */
  for (k=1; k<=depth; k++)
    {
      free(r[k]->x);
      free(r[k]);
      free(d[k]->x);
      free(d[k]);
      free(z[k]->x);
      free(z[k]);
      free(p[k]->x);
      free(p[k]);
    }
  return(0);
}

int lcd_solve (AMG_VECTOR *x, AMG_VECTOR *b)
{
  int i,j,k,l,n, *ra,row;
  double resid,residlast,alpha,alpha_deno, *a,*t,*s;
  AMG_SolverContext *sc=global_sc;
  double *beta_deno;
  
  beta_deno = malloc((sc->gmres_restart+1) *sizeof(double));  
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering LCD\n");
      AMG_Print(buffer);
    }
  n = AMG_MATRIX_N(A[0]);
  ra = AMG_MATRIX_RA(A[0]);
  a = AMG_MATRIX_A(A[0]);

  AMG_dcopy(p[0],b);                        /* copy rhs (b) into r */

  dmatminus(p[0],A[0],B,x);                 /* p = p - A*x, p is overwritten */

  if (sc->amg_prec_it)
  {
    AMG_dcopy(r[0],p[0]);          
    AMG_dset(z[0],0.0);                       /* z = sol, p=rhs ,r=def */
    for (k=0;k<sc->amg_prec_it;k++)
    {
      preconditioner(sc,0,depth,A,G,M,B,z,p,r);
      if (k<sc->amg_prec_it-1)
      {
	AMG_dcopy(r[0],p[0]);
	dmatminus(r[0],A[0],B,z[0]);
      }
    }
  }
  else
    AMG_dcopy(z[0],p[0]);
  
  resid = sqrt(AMG_ddot(z[0],z[0]));           /* norm of residual */
  start_residual=resid;
  residuals[0]=resid;
  residual_cnt++;

  if (sc->verbose>0)
  {
    sprintf(buffer,"LCD Iteration %4d %12.4lE\n",0,resid);
    AMG_Print(buffer);
  }
  
  if(resid  <= sc->res_norm_min )           /* stopping criterion fulfilled */
    {
      if (sc->verbose>0)
        {
          sprintf(buffer,"LCD Iteration %4d %12.4lE\n",0,resid);
          AMG_Print(buffer);
        }
      iteration_cnt=0;
      end_residual=resid;
      free(beta_deno);
      return(0);
    }
  
  if (sc->verbose>1)
    {
      sprintf(buffer,"Entering LCD iteration cycle\n");
      AMG_Print(buffer);
    }
  
  switch(sc->lcd_start_vector)
  {
     case 0:
     case 10:
        /* (b - Ax_0) / diag(A) */
        for (k=0;k<n;k++)
        {
           row =  ra[k];
           v[0]->x[k] = p[0]->x[k]/a[row];
           /* printf("%d %d %d %g \n ",k,row, AMG_MATRIX_NONZEROS(A[0]), a[row]);*/
        }
        break;
     default:
        /* M^{-1} rhs */
        if (sc->amg_prec_it)
        {
           AMG_dcopy(p[0],b);
           AMG_dcopy(r[0],p[0]);
           AMG_dset(d[0],0.0);                       /* d = sol, p=rhs ,r=def */
           for (k=0;k<sc->amg_prec_it;k++)
           {
              preconditioner(sc,0,depth,A,G,M,B,d,p,r);
              if (k<sc->amg_prec_it-1)
              {
                 AMG_dcopy(r[0],p[0]);
                 dmatminus(r[0],A[0],B,d[0]);
              }
           }
           AMG_dcopy(v[0],d[0]);	
        }
        else
           AMG_dcopy(v[0],b);                     /* copy rhs (b) into starting vector v[0] */
        break;
  }


  j=1; 
       
  while (j<=sc->maxit)
  {
    for(i=0; i<sc->gmres_restart && j<=sc->maxit; i++, j++)
    {
      if (sc->amg_prec_it)
      {
	AMG_dcopy(p[0],v[i]);
	AMG_dcopy(r[0],p[0]);
	AMG_dset(d[0],0.0);                       /* d = sol, p=rhs ,r=def */
	for (k=0;k<sc->amg_prec_it;k++)
	{
	  preconditioner_trans(sc,0,depth,A,G,M,B,d,p,r);
	  if (k<sc->amg_prec_it-1)
	  {
	    AMG_dcopy(r[0],p[0]);
	    dmatminus(r[0],A[0],B,d[0]);
	  }
	}
      }
      else
	AMG_dcopy(d[0],v[i]);

      dmattransmul(zv[i],A[0],B,d[0]);               /* zv[i]=A*d; */
      
      alpha = AMG_ddot(v[i],z[0]);
      alpha_deno = AMG_ddot(v[i],zv[i]);
      
      if (fabs(alpha_deno)<1e-120)
      {
	sprintf(buffer,"LCD : denominator in alpha too small\n");
	AMG_Print(buffer);
	exit(4711);
      }
      alpha = alpha/alpha_deno;

      AMG_daxpy(x,alpha,v[i]);                   /* update solution */
                                                 /* starting updating residual */
      dmatmul(p[0],A[0],B,v[i]);                 /* p = A*v[i]; */
          
      if (sc->amg_prec_it)
      {
	AMG_dcopy(r[0],p[0]);
	AMG_dset(d[0],0.0);                       /* d = sol, p=rhs ,r=def */
	for (k=0;k<sc->amg_prec_it;k++)
	{
	  preconditioner(sc,0,depth,A,G,M,B,d,p,r);
          if (k<sc->amg_prec_it-1)
	  {
	    AMG_dcopy(r[0],p[0]);
	    dmatminus(r[0],A[0],B,d[0]);
	  }
	}
      }
      else
	AMG_dcopy(d[0],p[0]);
                            	  
      AMG_daxpy(z[0],-alpha,d[0]);              /* update residual */
      
      residlast=resid;
      resid = sqrt(AMG_ddot(z[0],z[0]));

      if ((resid < sc->res_norm_min) || (resid<sc->red_factor*start_residual))
      {
      /* stopping criterion satisfied */	  
        if (sc->verbose>0)
	{
          sprintf(buffer,"LCD : iterations %d residual %g \n",j,resid);
          AMG_Print(buffer);
	}
	iteration_cnt=j;
	end_residual=resid;
        free(beta_deno);
	return(j);
      }

      if (sc->verbose>0)
      {
	sprintf(buffer,"LCD Iteration %4d %12.4lE %12.4lg\n",j,resid,resid/residlast);
	AMG_Print(buffer);
      }
      if (resid>sc->div_factor*start_residual)
      {
	AMG_Print("LCD iteration diverges !!!\n");
	exit(4711);
      }
      if (j<sc->maxit)
      {
        residuals[residual_cnt%AMG_CONV_RATE_BACK]=resid;
	residual_cnt++;
      }
      AMG_dcopy(v[i+1],z[0]);
      /* compute betas and update v[i+1] */
      for (k=0;k<=i;k++)
      {
	alpha = AMG_ddot(zv[k],v[i+1]);
	if (k==i)
	{
	    beta_deno[k] = AMG_ddot(zv[k],v[k]);
	    if (fabs(beta_deno[k])<1e-120)
	    {
		sprintf(buffer,"LCD : denominator in beta too small\n");
		AMG_Print(buffer);
		exit(4711);
	    }
	}
	alpha = -alpha/beta_deno[k];
	AMG_daxpy(v[i+1],alpha,v[k]);
      }
      
      finish = clock();
      if (finish>=start)
	elapsed_time+=(double)(finish-start)/CLOCKS_PER_SEC;
      else
	elapsed_time+=(double)(finish-start+TIME_WRAP)/CLOCKS_PER_SEC;    
      start = clock();
    }                                                   /* endfor i */
    /* restart */
    switch(sc->lcd_start_vector)
    {
	case 0:
	  AMG_dcopy(p[0],b);                        /* copy rhs (b) into r */
	  dmatminus(p[0],A[0],B,x);                 /* p = p - A*x, p is overwritten */
	  /* (b - Ax_0) / diag(A) */
	  for (k=0;k<n;k++)
	  {
	      row =  ra[k];
	      v[0]->x[k] = p[0]->x[k]/a[row];
	      /* printf("%d %d %d %g \n ",k,row, AMG_MATRIX_NONZEROS(A[0]), a[row]);*/
	  }	    
	  /* residual  */
       case 10:
       case 11:
          for (k=0;k<n;k++)
          {
             row =  ra[k];
             v[0]->x[k] = z[0]->x[k];///a[row];
             /* printf("%d %d %d %g \n ",k,row, AMG_MATRIX_NONZEROS(A[0]), a[row]);*/
          }  
          break;
       default:
          /* last conjugate direction before restart */
          AMG_dcopy(v[0],v[i]);                        /* copy last iterate into starting vector v[0] */
          break;
    }

  }
  /* endwhile */
  if (sc->verbose>0)
  {
    sprintf(buffer,"LCD : (maximal) iterations %d residual %g \n",sc->maxit,resid);
    AMG_Print(buffer);
  }
  iteration_cnt=j-1;
  end_residual=resid;              
  free(beta_deno);
  return(j);
}
