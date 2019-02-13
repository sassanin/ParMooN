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
/* File:          amg_blas.h                                                */
/*                                                                          */
/* Purpose:   BLAS on amg data structures                                   */
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
/* History:   01 FEB 1996 Begin                                             */
/*                          30 SEP 1997 redesign                            */
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
/* History:   1998/04/28 start adding routine for solving saddle            */
/*                         point problems                                   */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __AMG_BLAS__
#define __AMG_BLAS__

#include "amg_sp.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/* compile time constants defining static data size (i.e. arrays)           */
/* other constants                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/* BLAS Level 1 */
int    AMG_dset      (AMG_VECTOR *x, double a);
int    AMG_randomize (AMG_VECTOR *x);
int    AMG_dcopy     (AMG_VECTOR *x, AMG_VECTOR *y);
int    AMG_dscale    (AMG_VECTOR *x, double a);
int    AMG_daxpy     (AMG_VECTOR *x, double a, AMG_VECTOR *y);
double AMG_ddot      (AMG_VECTOR *x, AMG_VECTOR *y);
int AMG_daxpby (AMG_VECTOR *x, double a, AMG_VECTOR *y, double b, AMG_VECTOR *z);

/* BLAS Level 2 */
int AMG_dmatset   (AMG_MATRIX *A, double a);
int AMG_dmatcopy  (AMG_MATRIX *A, AMG_MATRIX *B);

int AMG_dmatmul   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmatminus (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmattransmul (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SCALAR2   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmatminus_SCALAR2 (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);

int AMG_dmatmul_SCALAR3   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmatminus_SCALAR3 (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);

int AMG_dmatmul_SCALAR6   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmatminus_SCALAR6 (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);

int AMG_dmatmul_SCALAR_VAL_LAZ_2   (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);
int AMG_dmatminus_SCALAR_VAL_LAZ_2 (AMG_VECTOR *x, AMG_MATRIX *A, AMG_MATRIX **B,AMG_VECTOR *y);

int AMG_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_SADDLE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_1_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_1_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

/* 2d saddle point problems */
int AMG_dmatmul_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_2_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_2_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SADDLE_2_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_2_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_2_TYPE_2 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_2_TYPE_2_MORTAR  (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatminus_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_gal_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_2 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_3 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_2_TYPE_4 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_);


/* 3d saddle point problems */
int AMG_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_B_Trans_dmatmul_SADDLE_3_TYPE_1 (AMG_VECTOR *x_,  AMG_MATRIX *A, AMG_MATRIX **B, AMG_VECTOR *y_);

int AMG_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatmul_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                  AMG_MATRIX **B, AMG_VECTOR *y_);
int AMG_A_dmatminus_BRAESS_SARAZIN_SADDLE_3_TYPE_1 (AMG_VECTOR *x_, AMG_MATRIX *A, 
                                                    AMG_MATRIX **B, AMG_VECTOR *y_);
#endif
