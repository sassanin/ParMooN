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
/*                                                                                 */
/* File:          amg_sp.h                                                    */
/*                                                                            */
/* Purpose:   interface to sparse matrix/vector data structure for amg             */
/*                                                                            */
/* Author:          Peter Bastian                                              */
/*                          Institut fuer Computeranwendungen III             */
/*                          Universitaet Stuttgart                            */
/*                          Pfaffenwaldring 27                                    */
/*                          70550 Stuttgart                                    */
/*                          email: peter@ica3.uni-stuttgart.de                    */
/*                          phone: 0049-(0)711-685-7003                            */
/*                          fax  : 0049-(0)711-685-7000                            */
/*                                                                            */
/* History:   01 FEB 1996 Begin                                                    */
/*                          30 SEP 1997 redesign                                    */
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
/* History:   1998/04/28 start changing structures for solving saddle             */
/*                         point problems                                            */
/*                                                                            */
/* Remarks:                                                                 */
/*                                                                            */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                            */
/*                                                                            */
/****************************************************************************/

#ifndef __AMG_SP__
#define __AMG_SP__

#include "amg_header.h"

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/* compile time constants defining static data size (i.e. arrays)            */
/* other constants                                                            */
/* macros                                                                    */
/*                                                                            */
/****************************************************************************/

#define AMG_LINES_PER_PAGE                60        /* repeat legend every this line */
#define AMG_COLS_PER_LINE                3        /* matrix entries per line */

/****************************************************************************/
/*                                                                            */
/* data structures exported by the corresponding source file                    */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/******** a block vector, i.e. n*blocksize entries    ***********************/
/****************************************************************************/

struct amg_vector 
{                                /* data + ptr to nodes structure */
  char name[AMG_NAME_SIZE];        /* name of this vector                 */
  int n;                        /* dimension of vector in blocks */
  int b;                        /* dimension of block                 */
  double *x;                        /* n*b doubles                         */
} ;
typedef 
struct amg_vector AMG_VECTOR;        /* a whole vector                 */

/****************************************************************************/
/******** functions for vectors                             *****************/
/****************************************************************************/

#define AMG_VECTOR_NAME(p)                ((p)->name)
#define AMG_VECTOR_N(p)                        ((p)->n)
#define AMG_VECTOR_B(p)                        ((p)->b)
#define AMG_VECTOR_X(p)                        ((p)->x)
#define AMG_VECTOR_ENTRY(p,i,ii)        ((p)->x[(i)*(p)->b+(ii)])

AMG_VECTOR *AMG_NewVector (int n, int b, char *name);

/****************************************************************************/
/******** pattern of a general m x n sparse block matrix    *****************/
/****************************************************************************/

struct amg_matrix 
{                                /* combines all links */
  char name[AMG_NAME_SIZE];        /* name of this matrix */
  int m;                        /* number of lines */
  int n;                        /* number of columns */
  int b;                        /* dimension of blocks: b x b */
  int bb;                        /* block size in doubles: ie bs=b*b */
  int system_as_scalar;                /* system treated as scalars if = 1 */
  int blocks_in_diag;                /* matrix is a block of block diag matrix with blocks_in_diag blocks*/
  int bandwidth;                /* bandwidth of the matrix */
  int nonzeros;                        /* number of nonzero blocks allocated */
  int connections;                /* nonzeros actually used */
  int active;                   /* active dof (inner + Neumann) */
  int level;                    /* number of level */
  int *ra;                        /* ra[i]: index of first entry of row i        */
  int *ja;                        /* aj[k]: col index of entry k */
  double *a;                        /* the matrix */
  int *ratr;                    /* index of first entry of column i for transposed matrix*/
  int *jatr;                        /* aj[k]: row index of transposed matris */
  int *postr;                        /* position index of transposed matris */
} ;
typedef 
struct amg_matrix  AMG_MATRIX;        /* a matrix */

/****************************************************************************/
/******** functions for vectors                             *****************/
/****************************************************************************/

#define AMG_MATRIX_NAME(p)                ((p)->name)
#define AMG_MATRIX_N(p)                        ((p)->n)
#define AMG_MATRIX_M(p)                        ((p)->m)
#define AMG_MATRIX_B(p)                        ((p)->b)
#define AMG_MATRIX_BB(p)                ((p)->bb)
#define AMG_MATRIX_SAS(p)                ((p)->system_as_scalar)
#define AMG_MATRIX_BLDI(p)                ((p)->blocks_in_diag)
#define AMG_MATRIX_NONZEROS(p)                ((p)->nonzeros)
#define AMG_MATRIX_CONNECTIONS(p)        ((p)->connections)
#define AMG_MATRIX_BW(p)                ((p)->bandwidth)
#define AMG_MATRIX_RA(p)                ((p)->ra)
#define AMG_MATRIX_JA(p)                ((p)->ja)
#define AMG_MATRIX_A(p)                        ((p)->a)

/* Construction */
AMG_MATRIX *AMG_NewMatrix (int n, int b, int nonzeros, int system_as_scalar, 
                           int blocks_in_diag,int allocate_row,int allocate_column,char *name);
AMG_MATRIX *AMG_CopyMatrix (AMG_MATRIX *A, char *name);
int AMG_SetRowLength (AMG_MATRIX *A, int i, int l);
int AMG_FindEntry (AMG_MATRIX *A, int i, int j);
int AMG_InsertEntry (AMG_MATRIX *A, int i, int j);
int AMG_InsertValues (AMG_MATRIX *A, int i, int j, double *aij);
int AMG_AddValues (AMG_MATRIX *A, int i, int j, double *aij);

/* input / output */
int AMG_PrintVector   (int k, AMG_VECTOR **vlist, char *text);
int AMG_PrintMatrix   (AMG_MATRIX *A, char *text);

#endif
