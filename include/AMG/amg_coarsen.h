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
/* File:          amg_coarsen.h                                                                                                        */
/*                                                                                                                                                        */
/* Purpose:   algebraic multigrid coarse grid setup                                                        */
/*                                                                                                                                                        */
/* Author:          Peter Bastian                                                                                  */
/*                          Institut fuer Computeranwendungen III                                                 */
/*                          Universitaet Stuttgart                                                                                */
/*                          Pfaffenwaldring 27                                                                                        */
/*                          70550 Stuttgart                                                                                                */
/*                          email: peter@ica3.uni-stuttgart.de                                                        */
/*                          phone: 0049-(0)711-685-7003                                                                        */
/*                          fax  : 0049-(0)711-685-7000                                                                        */
/*                                                                                                                                                        */
/* History:   29 FEB 1996 Begin                                                                                                */
/*                          01 OKT 1997 redesign                                                                                        */
/*                                                                                                                                                        */
/* Remarks:                                                                                                                                 */
/*                                                                                                                                                        */
/****************************************************************************/

/* RCS_ID
$Header: /homes/matthies/ARCHIVE_MooNMD/MooNMD/include/AMG/amg_coarsen.h,v 1.3 2004/12/10 08:33:37 john Exp $
*/

/****************************************************************************/
/*                                                                                                                                                        */
/* auto include mechanism and other include files                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#ifndef __AMG_COARSEN__
#define __AMG_COARSEN__

#include "amg_sp.h"

/****************************************************************************/
/*                                                                                                                                                        */
/* defines in the following order                                                                                        */
/*                                                                                                                                                        */
/*                  compile time constants defining static data size (i.e. arrays)        */
/*                  other constants                                                                                                        */
/*                  macros                                                                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#define AMG_MAX_COMP                        5                /* max number of components in syste*/

#define AMG_MAX_LEVELS                    32                /* max no of levels in hierarchy        */
#define AMG_MAX_FRONT                        4096        /* max number of front nodes                */
#define AMG_MAX_CLUSTER                        256                /* max size of a cluster                        */
#define AMG_MAX_STACK                        256                /* size of seed stack                                */
#define AMG_MAX_ROW                                2048                /* max number of nonzeros in row        */


#define AMG_UNSYM                                1                /* unsymmetric dependency strategy        */
#define AMG_SYM                                        2                /* symmetric dependency strategy        */

/****************************************************************************/
/*                                                                                                                                                        */
/* data structures exported by the corresponding source file                                */
/*                                                                                                                                                        */
/****************************************************************************/

typedef struct {                                        /* parameters of coarsening strategy        */
  int verbose;                                        /* be verbose                                                        */
  double alpha;                                        /* "relative" strongness                                 */
  double beta;                                        /* "absolute" strongness                                */
  int mincluster;                                        /* minimal cluster size                                        */
  int maxcluster;                                        /* maximum cluster size                                        */
  int maxdistance;                                /* maximum distance in one cluster                */
  int maxconnectivity;                        /* limit for connectivity                                 */
  int coarsentarget;                                /* quit coarsening if nodes reached                */
  int depthtarget;                                /* create at most so many levels                */
  double coarsenrate;                                /* quit if coarsening is too slow                */
  int major;                                                /* use major component strategy if >=0  */
  int dependency;                                        /* selects dependency strategy                          */
    double time;                                    /* time for coarsening */
} AMG_CoarsenContext;

struct  amg_graph{                              /* graph data structure */
  int n;                                        /* number of nodes (this is n from A) */
  int e;                                        /* number of edges (nonzeros from A)  */
  int *ra,*ja;                                        /* connectivity from A */
  int *ca;                                        /* cluster array (fine to coarse map) */
  char *na;                                        /* node array (flags per node) */
  char *la;                                        /* link array */
  float *da;                                        /* damping array with automatic damping */
  int clusters;                                        /* total number of clusters build */
  int conclusters;                                /* number of nonisolated clusters */
  int system_as_scalar;                                /* copied from matrix */
};  

typedef  struct amg_graph AMG_GRAPH;

#define AMG_GRAPH_N(p)                                (p)->n
#define AMG_GRAPH_E(p)                                (p)->e
#define AMG_GRAPH_SAS(p)                        (p)->system_as_scalar
#define AMG_GRAPH_RA(p)                                (p)->ra
#define AMG_GRAPH_JA(p)                                (p)->ja
#define AMG_GRAPH_CA(p)                                (p)->ca
#define AMG_GRAPH_NA(p)                                (p)->na
#define AMG_GRAPH_LA(p)                                (p)->la
#define AMG_GRAPH_DA(p)                                (p)->da


/****************************************************************************/
/*                                                                                                                                                        */
/* functions                                                                                                                                */
/*                                                                                                                                                        */
/****************************************************************************/

int AMG_BuildHierarchy (AMG_CoarsenContext *cc, AMG_MATRIX *A, 
        AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS]);

int AMG_BuildHierarchy_Saddle (AMG_CoarsenContext *cc, AMG_MATRIX *A, 
                               AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS], 
                               AMG_MATRIX *schur, AMG_MATRIX *H_schur[AMG_MAX_LEVELS], 
                               AMG_GRAPH *G_schur[AMG_MAX_LEVELS],
                               AMG_MATRIX *B[AMG_MAX_LEVELS]);
#endif
