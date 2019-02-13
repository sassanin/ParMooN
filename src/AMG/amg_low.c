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


/****************************************************************************/
/*                                                                                                                                                        */
/* include files                                                                                                                        */
/*                          system include files                                                                                        */
/*                          application include files                                                                         */
/*                                                                                                                                                        */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <amg_header.h>
#include <amg_low.h>

void Cpp_Print(char *s);

/****************************************************************************/
/*                                                                                                                                                        */
/* defines in the following order                                                                                        */
/*                                                                                                                                                        */
/*                  compile time constants defining static data size (i.e. arrays)        */
/*                  other constants                                                                                                        */
/*                  macros                                                                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#undef DEBUG

/****************************************************************************/
/*                                                                                                                                                        */
/* data structures used in this source file (exported data structures are        */
/*                  in the corresponding include file!)                                                                */
/*                                                                                                                                                        */
/****************************************************************************/



/****************************************************************************/
/*                                                                                                                                                        */
/* definition of exported global variables                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/



/****************************************************************************/
/*                                                                                                                                                        */
/* definition of variables global to this source file only (static!)                */
/*                                                                                                                                                        */
/****************************************************************************/

/* RCS_ID
$Header: /homes/matthies/ARCHIVE_MooNMD/MooNMD/src/AMG/amg_low.c,v 1.3 2002/10/17 08:20:36 matthies Exp $
*/

/* user installable print handler */
static AMG_PrintFuncPtr AMG_UserPrint=NULL;
static AMG_MallocFuncPtr AMG_UserMalloc=NULL;

/****************************************************************************/
/*                                                                                                                                                        */
/* forward declarations of functions used before they are defined                        */
/*                                                                                                                                                        */
/****************************************************************************/


/****************************************************************************/
/*D
   AMG_InstallPrintHandler - install user print function for amg package
   
   SYNOPSIS:
   typedef void (*AMG_PrintFuncPtr) (char *);
   int AMG_InstallPrintHandler (AMG_PrintFuncPtr print);
   
   PARAMETERS:
.  print - Function of type AMG_PrintFuncPtr.
   
   DESCRIPTION:
   This function allows you to define your own print function, i.e. a function
   that writes a string to the console. This is necessary since printf
   is usually not available in windowed environments. If no handler is ever
   installed then AMG will use the standard fputs function to write out strings.
   
   RETURN VALUE:
.n AMG_OK
   
D*/
/****************************************************************************/

int AMG_InstallPrintHandler (AMG_PrintFuncPtr print)
{
        AMG_UserPrint = print;
        return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_Print - function to write a string
   
   SYNOPSIS:
   int AMG_Print (char *s);
   
   PARAMETERS:
.  s - string to print.
   
   DESCRIPTION:
   In order to enable proper working of amg within other codes
   all diagnostic output should be done through these functions.
   The usual way is to use sprintf to provide
   formatted output into a buffer which is then written with AMG_Print
   to the screen.
   
   RETURN VALUE:
.n AMG_OK
   
D*/
/****************************************************************************/

int AMG_Print (char *s)
{
  Cpp_Print(s);
  if (AMG_UserPrint==NULL)
    fputs(s,stdout);
  else
    (*AMG_UserPrint)(s);
  
  return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_InstallMallocHandler - install user malloc function
   
   SYNOPSIS:
   typedef void (*AMG_MallocFuncPtr) (size_t n);
   int AMG_InstallMallocHandler (AMG_MallocFuncPtr mall);
   
   PARAMETERS:
.  mall - pointer to user definable malloc function.
   
   DESCRIPTION:
   This allows to redefine the memory allocation function used by amg
   in order to be able to use amg within other programs having
   their own memory management. All memory allocation is usually static
   within amg, i.e. memory is allocated but never freed again.
   
   RETURN VALUE:
.n AMG_OK
.n AMG_FATAL   
D*/
/****************************************************************************/

int AMG_InstallMallocHandler (AMG_MallocFuncPtr mall)
{
        AMG_UserMalloc = mall;
        return(AMG_OK);
}

/****************************************************************************/
/*D
   AMG_Malloc - malloc function for amg
   
   SYNOPSIS:
   void *AMG_Malloc (size_t nbytes);
   
   PARAMETERS:
.  size_t - nbytes.
   
   DESCRIPTION:
   All dynamic memory allocation in amg is done thru this function.
   It allows redefinition.
   
   RETURN VALUE:
.n pointer to memory or NULL
   
D*/
/****************************************************************************/

void *AMG_Malloc (size_t nbytes)
{
  void *ptr;
  if (AMG_UserMalloc==NULL)
    {
      ptr = malloc(nbytes);
      if (ptr==NULL)
        {
          AMG_Print("ERROR : Cannot allocate memory in amg library !!!\n");
          AMG_Print("ERROR : Program teminated !!!\n");
          exit(4711);
        }
      else
        return(ptr);
    }
  else
    return(AMG_UserMalloc(nbytes));
  return(ptr);
}
