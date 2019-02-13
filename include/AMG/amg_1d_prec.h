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
/* File:          amg_1d_prec.h                                                    */
/*                                                                            */
/* Purpose:   header file for amg_1d_prec.c                                 */
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
/* Remarks:   1998/06/08 start                                              */
/*                                                                            */
/****************************************************************************/
int prepare_steplength_control(int *A_length, int *B_length, int depth);
int clear_steplength_control(int depth);

int ssor (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int ssor_slc (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int sor (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int jac (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int ex (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int prepare_ilu(AMG_SolverContext *sc,AMG_MATRIX *A[AMG_MAX_LEVELS],
                 AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth); 
int clear_ilu(AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth);
int clear_ilut(AMG_MATRIX *M[AMG_MAX_LEVELS],int fine,int depth);
int ilu (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int clear_mgc(AMG_MATRIX *A[AMG_MAX_LEVELS],AMG_GRAPH *G[AMG_MAX_LEVELS],int depth);
int prepare_coupled_mgc(AMG_SolverContext *sc,int *A_length,int *B_length, int depth);
int clear_coupled_mgc(AMG_MATRIX *A[AMG_MAX_LEVELS],AMG_GRAPH *G[AMG_MAX_LEVELS],int depth);
int mgc (AMG_SolverContext *sc, int k, int depth,                                 
         AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],                 
         AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
         AMG_VECTOR *x[AMG_MAX_LEVELS],          
         AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int clear_ex(AMG_MATRIX *M[AMG_MAX_LEVELS],int depth);
AMG_MATRIX *prepare_ex (AMG_MATRIX *A);
AMG_MATRIX *ILUDecomposition(AMG_SolverContext *sc, AMG_MATRIX *A);
AMG_MATRIX *ILUTDecomposition(AMG_SolverContext *sc, AMG_MATRIX *A);


int pc_restrict (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);
int pc_prolongate_auto (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);
int pc_restrict_2d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);
int pc_prolongate_auto_2d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);
int pc_restrict_saddle_2d (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);
int pc_prolongate_auto_saddle_2d (AMG_GRAPH *g, AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);
int pc_restrict_3d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);
int pc_prolongate_auto_3d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);
int pc_restrict_6d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse);
int pc_prolongate_auto_6d (AMG_GRAPH *g,AMG_GRAPH *g1, AMG_VECTOR *fine, AMG_VECTOR *coarse, double *damp);

int ApplyRowEquilibration(AMG_SolverContext *sc,AMG_MATRIX *A);
int EquilibrateRHS(AMG_VECTOR *b,AMG_VECTOR *equi);
int DisApplyRowEquilibration(AMG_SolverContext *sc,AMG_MATRIX *A,AMG_VECTOR *b);


int sor_trans (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int ssor_trans (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int ilu_trans (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int ComputeArraysForTransposedMatrix(AMG_MATRIX *A);
