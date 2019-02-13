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
/* File:          amg_solve.h                                                    */
/*                                                                            */
/* Purpose:   solvers for AMG                                                    */
/*                                                                            */
/* Author:          Peter Bastian                                                     */
/*                          Institut fuer Computeranwendungen III             */
/*                          Universitaet Stuttgart                            */
/*                          Pfaffenwaldring 27                                    */
/*                          70550 Stuttgart                                    */
/*                          email: peter@ica3.uni-stuttgart.de                    */
/*                          phone: 0049-(0)711-685-7003                            */
/*                          fax  : 0049-(0)711-685-7000                            */
/*                                                                            */
/* History:   05 FEB 1996 Begin                                                    */
/*              02 OKT 1997 redesign                                            */
/*                                                                            */
/* Remarks:                                                                     */
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
/* Remarks:                                                                 */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                            */
/*                                                                            */
/****************************************************************************/

#ifndef __AMG_SOLVE__
#define __AMG_SOLVE__

#include "amg_coarsen.h"

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/* compile time constants defining static data size (i.e. arrays)            */
/* other constants                                                            */
/* macros                                                                    */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* data structures exported by the corresponding source file                    */
/*                                                                            */
/****************************************************************************/

#define AMG_DJAC                        1                /* smoother types */
#define AMG_SOR                                2
#define AMG_SSOR                        3
#define AMG_ILU                                4
#define AMG_MGC                                5
#define AMG_ILUT                        6
#define AMG_SCHUR_COMPLEMENT                7
#define AMG_EX                                17              /* same as solver */

#define AMG_LS                                11                /* solver types        */
#define AMG_CG                                12
#define AMG_BCGS                        13
#define AMG_GMRES_LEFT                        14
#define AMG_GMRES_RIGHT                        15
#define AMG_GMRES_FLEX                        16
#define AMG_EXACT                        17 

#define SCALAR                          0
#define SADDLE_1                        1
#define SADDLE_2_TYPE_1                 2
#define SADDLE_3                        3

#define ILUT_SELECTION_SORT             0
#define ILUT_BUBBLE_SORT                1

#define AMG_MAX_GMRES_RESTART                1000

#define AMG_CONV_RATE_BACK              10

typedef struct 
{                                        /* parameters for solver */
  int verbose;                                /* be verbose */
  int system_type;                      /* type of the system */

  /* fine grid solver */
  int solver;                                /* type of solver to be used */
  int preconditioner;                        /* type of preconditioner */
  int maxit;                                /* max number of iterations */
  int ex_maxit;                         /* 1 to execute exactly maxit iterations */
  double red_factor;                        /* required reduction in residual norm */
  double res_norm_min;                        /* convergence limit */
  int amg_prec_it;                      /* number preconditioner of iterations */ 
  double amg_prec_red_factor;           /* required reduction in residual norm for preconditioner*/ 
  int gmres_restart;                    /* restart parameter in gmres, also used in lcd */
  int lcd_start_vector;                        /* choose starting lcd vector */
  
  /* coarse grid solver */
  int coarse_smoother;                        /* type of smoother on coarse grid */
  int coarse_maxit;                        /* iteration number for coarse grid sol */
  double coarse_red_factor;                /* required reduction in residual norm */
  double omega_coarse[AMG_MAX_COMP];        /* damping factor per component */
  
  /* multigrid cycle */
  int n1,n2;                                /* pre and post smoothing */
  int gamma;                                /* cycle form */
  double omega_p[AMG_MAX_COMP];                /* damping factor for interpolation */
  
  /* smoother */
  int smoother;                                /* type of smoother to be used */
  double smoother_red_factor;           /* required reduction in residual norm */
  double omega[AMG_MAX_COMP];                /* damping factor per component        */
  int step_length_control_fine;         /* step length control on the finest level */ 
  int step_length_control_all;          /* step length control on all levels */ 
  double ilu_beta;                      /* parameter for ilu(beta) */

  /*ilut*/
  double ilut_tol;                      /* tolerance for first dropping criterion */
  int ilut_fillin;                      /* number of additional fill in's */
  int ilut_sort;                        /* sorting strategy for ilut */

  /* saddle point problems parameter */
  int schur_inv_of_A;                   /* method to approximate inv(A) */
  double schur_iteration_damp;          /* damping factor in schur complement iteration */
  int schur_iteration_maxit;            /* maximal number of it. in schur complement iteration */

} AMG_SolverContext;


/****************************************************************************/
/*                                                                            */
/* functions                                                                    */
/*                                                                            */
/****************************************************************************/

int AMG (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,  AMG_MATRIX **B, 
         AMG_VECTOR *x, AMG_VECTOR *b,char *OutFileName);

#endif
