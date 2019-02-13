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
/* File:          amg_iter.h                                                     */
/*                                                                            */
/* Purpose:   header file for amg_iter.c                                    */
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

#ifndef __AMG_ITER__
#define __AMG_ITER__

#include "amg_sp.h"

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/* compile time constants defining static data size (i.e. arrays)            */
/* other constants                                                            */
/* macros                                                                    */
/*                                                                            */
/****************************************************************************/

#define AMG_EX_MAT(m,b,i,j)                        ((m)[2*(b)*(i) + (j)])


/****************************************************************************/
/*                                                                            */
/* data structures exported by the corresponding source file                    */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* functions                                                                    */
/*                                                                            */
/****************************************************************************/

/* linear iteration kernels */
int         AMG_jac                 (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int         AMG_sorf                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int         AMG_sorb                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int         AMG_sorfb               (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);
int         AMG_iluf                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d);
int         AMG_ilub                (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d);

int         AMG_EXDecomposeMatrix (double *Mat, int bw, int n);
int         AMG_EXApplyLU (double *Mat, int bw, int n,  int blocks,double *Vec);

int         AMG_sorf_trans          (AMG_MATRIX *A, AMG_VECTOR *v, AMG_VECTOR *d, double *omega);

#endif
