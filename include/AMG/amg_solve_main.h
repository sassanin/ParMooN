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
/* File:          amg_solve_main.h                                          */
/*                                                                          */
/* Purpose:   header file for amg_solve_main.c                              */
/*                                                                          */
/* Author:          Peter Bastian                                           */
/*                          Institut fuer Computeranwendungen III           */
/*                          Universitaet Stuttgart                          */
/*                          Pfaffenwaldring 27                              */
/*                          70550 Stuttgart                                 */
/*                          email: peter@ica3.uni-stuttgart.de              */
/*                          phone: 0049-(0)711-685-7003                     */
/*                          fax  : 0049-(0)711-685-7000                     */
/*                                                                          */
/* History:   05 FEB 1996 Begin                                             */
/*              02 OKT 1997 redesign                                        */
/*                                                                          */
/* Remarks:                                                                 */
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
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __AMG_SOLVE__
#define __AMG_SOLVE__

#include "amg_coarsen.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/* compile time constants defining static data size (i.e. arrays)           */
/* other constants                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

#define AMG_DJAC                        1                /* smoother types */
#define AMG_SOR                         2
#define AMG_SSOR                        3
#define AMG_ILU                         4
#define AMG_MGC                         5
#define AMG_ILUT                        6
#define AMG_SCHUR_COMPLEMENT            7
#define AMG_SCHUR_CG                    8
#define AMG_SCHUR_GMRES                 9
#define AMG_SCHUR_GMRES_BCGS           10
#define AMG_BRAESAR                    11
#define AMG_EX                         17                /* same as solver */
#define AMG_NO_PRECONDITIONER       -4711

#define AMG_LS                         11                /* solver types        */
#define AMG_CG                         12
#define AMG_BCGS                       13
#define AMG_GMRES_LEFT                 14
#define AMG_GMRES_RIGHT                15
#define AMG_GMRES_FLEX                 16
#define AMG_EXACT                      17 
#define AMG_MIXED_BCGS_CGS             18
#define AMG_LCD                        19

#define SCALAR1                         0         /* scalar system with 1 rhs */
#define SCALAR2                         1         /* scalar system with 2 rhs */
#define SCALAR3                         2         /* scalar system with 3 rhs */
#define SCALAR6                         5         /* scalar system with 6 rhs */

#define SADDLE_1                        10
#define SADDLE_1_TYPE_1                 11         /* scalar mortar systems */
#define SADDLE_2_TYPE_1                 12         /* 2D Stokes type systems */
#define SADDLE_2_TYPE_2                 13         /* 2D Oseen type systems mit SDFEM */
#define SADDLE_2_TYPE_3                 14         /* 2D Oseen type systems mit SDFEM */
#define SADDLE_2_TYPE_4                 15         /* 2D Oseen type systems with all A blocks */
#define BRAESS_SARAZIN_SADDLE_2_TYPE_1  16         /* Braess--Sarzin smoother for 2D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_2_TYPE_2  17         /* Braess--Sarzin smoother for 2D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_2_TYPE_3  18         /* Braess--Sarzin smoother for 2D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_2_TYPE_4  19         /* Braess--Sarzin smoother for 2D Stokes type systems */

#define SADDLE_2_TYPE_2_MORTAR          22         /* 2D Oseen type systems with Mortar */

#define SADDLE_3                        100
#define SADDLE_3_TYPE_1                 112         /* 3D Stokes type systems */
#define SADDLE_3_TYPE_2                 113         /* 3D Oseen type systems mit SDFEM */
#define SADDLE_3_TYPE_3                 114         /* 3D Oseen type systems mit SDFEM */
#define SADDLE_3_TYPE_4                 115         /* 3D Oseen type systems with all A blocks */
#define BRAESS_SARAZIN_SADDLE_3_TYPE_1  116         /* Braess--Sarzin smoother for 3D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_3_TYPE_2  117         /* Braess--Sarzin smoother for 3D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_3_TYPE_3  118         /* Braess--Sarzin smoother for 3D Stokes type systems */
#define BRAESS_SARAZIN_SADDLE_3_TYPE_4  119         /* Braess--Sarzin smoother for 3D Stokes type systems */

#define SCALAR_VAL_LAZ_2                101         /* Vassilevski/Lazarov preconditioner */        

#define ILUT_QUICK_SPLIT_0              0
#define ILUT_QUICK_SPLIT_1              1
#define ILUT_QUICK_SPLIT_2              2

#define AMG_MAX_GMRES_RESTART                1000
#define AMG_MAX_SCHUR_GMRES_RESTART        100
#define max_switch                      100  /* cgs-bcgstab */

#define AMG_CONV_RATE_BACK              10

#define CONSTANT                        0
#define PLUS_CONSTANT                   1
#define TIMES_CONSTANT                  2
#define SQUARED                         3   

typedef struct 
{                                        /* parameters for solver */
  int verbose;                                /* be verbose */
  int system_type;                      /* type of the system */

  /* fine grid solver */
  int solver;                                /* type of solver to be used */
  int preconditioner;                        /* type of preconditioner */
  int maxit;                                /* max number of iterations */
  int minit;                                /* min number of iterations */
  int ex_maxit;                         /* 1 to execute exactly maxit iterations */
  double red_factor;                        /* required reduction in residual norm */
  double res_norm_min;                        /* convergence limit */
  int amg_prec_it;                      /* number preconditioner of iterations */ 
  double amg_prec_red_factor;           /* required reduction in residual norm for preconditioner */ 
  int gmres_restart;                    /* restart parameter in gmres, also used in LCD */
  int lcd_start_vector;                        /* choose starting lcd vector */
  double mixed_bcgs_cgs_switch_tol;     /* switching tolerance in mixed_bcgs_cgs */ 
  double div_factor;                        /* factor for increase in residual norm for divergence */
  int row_equilibration;                /* apply row equilibration */ 
  
  /* coarse grid solver */
  int coarse_smoother;                        /* type of smoother on coarse grid */
  int coarse_maxit;                        /* iteration number for coarse grid sol */
  double coarse_red_factor;                /* required reduction in residual norm */
  double omega_coarse[AMG_MAX_COMP];        /* damping factor per component */
  
  /* multigrid cycle */
  int n1[AMG_MAX_LEVELS];                /* pre smoothing */
  int n2[AMG_MAX_LEVELS];                /* post smoothing */
  int gamma;                                /* cycle form */
  double omega_p[AMG_MAX_COMP];                /* damping factor for interpolation */
  int smoothing_steps;                  /* number of smoothing steps on coarser grids */
  int n1_param;                         /* parameter for changing pre smoothing steps */
  int n2_param;                         /* parameter for changing post smoothing steps */
 
  /* smoother */
  int smoother;                                /* type of smoother to be used */
  double smoother_red_factor;           /* required reduction in residual norm */
  double omega[AMG_MAX_COMP];                /* damping factor per component        */
  int step_length_control_fine;         /* step length control on the finest level */ 
  int step_length_control_all;          /* step length control on all levels */ 
  double ilu_beta;                      /* parameter for ilu(beta) */

  /*ilut*/
  double ilut_tol;                      /* tolerance for first dropping criterion */
  int ilut_absolute_fillin;             /* number of additional fill in's */
  double ilut_relative_fillin;             /* relative fill in to matrix dimension in percent */
  int ilut_sort;                        /* sorting strategy for ilut */
  double sor_omega;                     /* damping factor in Jacobi, SOR, SSOR */

  /* saddle point problems parameter */
  int schur_inv_of_A;                   /* method to approximate inv(A) */
  int schur_inv_of_A_maxit;             /* max number of iterations if inv(A) is approximated by an iterative method */
  double schur_iteration_damp;          /* damping factor in schur complement iteration */
  int schur_iteration_maxit;            /* maximal number of it. in schur complement iteration */
  int schur_step_length_control;        /* step length control in schur complement iteration */
  int condense;                         /* condense after symmetrizing the problem */

  double vas_laz_delta;                    /* delta for Vassilevski/Lazaroc preconditioner */                   
  int braess_sarazin_matrix;            /* matrix for Braess-Sarazin smoother */
  double braess_sarazin_alpha;             /* damping factor for Braess-Sarazin smoother */


} AMG_SolverContext;

typedef int (*IterProcPtr) (AMG_SolverContext *sc, int level, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS],AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],AMG_VECTOR *b[AMG_MAX_LEVELS], 
        AMG_VECTOR *d[AMG_MAX_LEVELS]);

typedef int (*MultProcPtr) (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B, 
                            AMG_VECTOR *y);

typedef int (*RestrictProcPtr) (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);

typedef int (*InterpolationProcPtr) (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);

/****************************************************************************/
/*                                                                            */
/* functions                                                                    */
/*                                                                            */
/****************************************************************************/

int AMG (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,  AMG_MATRIX **B, 
         AMG_VECTOR *x, AMG_VECTOR *b);

int AMG_Build (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b);

int AMG_Solve (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b);

int AMG_Delete (AMG_SolverContext *sc, AMG_CoarsenContext *cc, AMG_MATRIX *A,
         AMG_MATRIX **B0, AMG_VECTOR *x, AMG_VECTOR *b);

#endif
