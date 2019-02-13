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
/* File:      amg_2d_prec.h                                                       */
/*                                                                            */
/* Purpose:   header file for amg_2d_prec.c                                 */
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
/* Remarks:   1998/02/24 ILU                                                 */
/*              1998/02/27 step length control                                */
/*              1998/06/03 ILUT                                               */
/*                                                                            */
/****************************************************************************/
int prepare_schur_complement_fixed(int *A_length,int *B_length, int depth);
int clear_schur_complement_fixed(int depth);
int schur_complement_fixed (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int prepare_schur_complement_cg(int *A_length,int *B_length, int depth);
int clear_schur_complement_cg(int depth);
int schur_complement_cg(AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int prepare_schur_complement_gmres(AMG_SolverContext *sc,int *A_length,int *B_length, int depth);
int clear_schur_complement_gmres(AMG_SolverContext *sc,int depth);
int schur_complement_gmres(AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);
int prepare_schur_bcgs_solve(AMG_SolverContext *sc,int *A_length, int *B_length, int depth);
int schur_complement_gmres_bcgs (AMG_SolverContext *sc, int k, int depth,
        AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],
        AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS],
        AMG_VECTOR *x[AMG_MAX_LEVELS],
        AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int prepare_braess_sarazin_smoother(AMG_SolverContext *sc,int *A_length,int *B_length, int depth);
int clear_braess_sarazin_smoother(int depth);
int braess_sarazin_smoother (AMG_SolverContext *sc, int k, int depth,                                 
         AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS],                 
         AMG_MATRIX *M[AMG_MAX_LEVELS], AMG_MATRIX *B[AMG_MAX_LEVELS], 
         AMG_VECTOR *x[AMG_MAX_LEVELS],          
         AMG_VECTOR *b[AMG_MAX_LEVELS], AMG_VECTOR *d[AMG_MAX_LEVELS]);

int BuildDiagSchurComplement(AMG_SolverContext *sc,AMG_MATRIX *A, AMG_MATRIX **B); 
int SymmetrizeSaddleProblem (AMG_SolverContext *sc,AMG_MATRIX *A,AMG_MATRIX *B[AMG_MAX_LEVELS],
                       AMG_VECTOR *b_);
int Build_B_Block_for_Dirichlet(AMG_SolverContext *sc,AMG_MATRIX *A, 
                                AMG_MATRIX *B[AMG_MAX_LEVELS]);
