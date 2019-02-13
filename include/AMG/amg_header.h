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
/*                                                                                                                                                        */
/* File:          amg_header.c                                                                                                        */
/*                                                                                                                                                        */
/* Purpose:   general header for common things (return values, misc..)                */
/*                                                                                                                                                        */
/* Author:          Peter Bastian                                                                                                         */
/*                          Institut fuer Computeranwendungen III                                                 */
/*                          Universitaet Stuttgart                                                                                */
/*                          Pfaffenwaldring 27                                                                                        */
/*                          70550 Stuttgart                                                                                                */
/*                          email: peter@ica3.uni-stuttgart.de                                                        */
/*                          phone: 0049-(0)711-685-7003                                                                        */
/*                          fax  : 0049-(0)711-685-7000                                                                        */
/*                                                                                                                                                        */
/* History:   28 Jan 1996 Begin                                                                                                */
/*            02 Apr 1996 new memory allocation strategy                                        */
/*            30 Sep 1997 redesign                                                                                        */
/*                                                                                                                                                        */
/* Remarks:                                                                                                                                 */
/*                                                                                                                                                        */
/****************************************************************************/

/* RCS_ID
$Header: /homes/matthies/ARCHIVE_MooNMD/MooNMD/include/AMG/amg_header.h,v 1.1.1.1 2001/01/12 12:55:18 matthies Exp $
*/

/****************************************************************************/
/*                                                                                                                                                        */
/* auto include mechanism and other include files                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#ifndef __AMGHEADER__
#define __AMGHEADER__

/****************************************************************************/
/*                                                                                                                                                        */
/* defines in the following order                                                                                        */
/*                                                                                                                                                        */
/*                  compile time constants defining static data size (i.e. arrays)        */
/*                  other constants                                                                                                        */
/*                  macros                                                                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

/* general sizes */
#define AMG_NAME_SIZE                        32  /* for names of objects                                        */

/* general return values*/
#define AMG_OK                                        0        /* operation succeded                                        */
#define AMG_NULL                                NULL/* null pointer                                                        */
#define AMG_FATAL                                9999/* fatal error                                                        */

/* misc macros */
#define AMG_MIN(x,y)            (((x)<(y)) ? (x) : (y))
#define AMG_MAX(x,y)            (((x)>(y)) ? (x) : (y))
/* #define AMG_ABS(x)                    (((x)>=0) ?  (x) : (-x))*/
#define AMG_ABS(i) (((i)<0) ? (-(i)) : (i))

#endif
