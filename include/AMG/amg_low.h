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
/* File:          amg_low.c                                                                                                                */
/*                                                                                                                                                        */
/* Purpose:   handlers, memory management for amg                                                        */
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
$Header: /homes/matthies/ARCHIVE_MooNMD/MooNMD/include/AMG/amg_low.h,v 1.2 2003/03/13 10:09:56 matthies Exp $
*/

/****************************************************************************/
/*                                                                                                                                                        */
/* auto include mechanism and other include files                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#ifndef __AMGLOW__
#define __AMGLOW__

#include "amg_header.h"

/****************************************************************************/
/*                                                                                                                                                        */
/* defines in the following order                                                                                        */
/*                                                                                                                                                        */
/*                  compile time constants defining static data size (i.e. arrays)        */
/*                  other constants                                                                                                        */
/*                  macros                                                                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

/****************************************************************************/
/******** typedefs     ******************************************************/
/****************************************************************************/

typedef void (*AMG_PrintFuncPtr) (char *);
typedef void * (*AMG_MallocFuncPtr) (size_t);

/****************************************************************************/
/******** management functions  *********************************************/
/****************************************************************************/

/* string i/o */
int       AMG_InstallPrintHandler (AMG_PrintFuncPtr print);
int       AMG_Print         (char *s);
int           AMG_RedirectToFile (char *name);
int           AMG_RedirectToScreen (void);

/* sp internal memory handling */
int       AMG_InstallMallocHandler (AMG_MallocFuncPtr mall);
void      *AMG_Malloc  (size_t);

#endif
