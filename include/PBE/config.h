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
   
#ifndef __CONFIG_H
#define __CONFIG_H
/**
 **
 ** provides basic types and configuration settings
 **
 ** arch-tag: 3c370550-398d-496c-ac9e-b4b7cb79c45f
 **/

/*
 * included header
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "seconds.h"

/*
 * basic types
 */

typedef enum { 
	FEPC_FALSE, 
	FEPC_TRUE 
} bool_t;

typedef double  fepc_real_t;

/*
 * consistency check
 */

#if !defined(NDEBUG)
#  define ASSERT( cond )   assert( cond )
#else
#  define ASSERT( cond )
#endif

/*
 * enables debugging/output
 */

#define DEBUG  if ( false )
#define PRINT  if ( false )

/*
 * inline function decl.
 */

#if defined(__USE_ISOC99) && ( defined(__GNUC__) || defined(__ICC) || defined(__ECC) )
#define _INLINE_  inline
#else
#define _INLINE_  static
#endif

#endif  /* __CONFIG_H */
