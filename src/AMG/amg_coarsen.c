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
/* File:          amg_coarsen.c                                                    */
/*                                                                            */
/* Purpose:   allocate and compute coarse grid matrices                            */
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
/* History:   29 FEB 1996 Begin                                                    */
/*              01 OKT 1997 redesign                                            */
/*                                                                            */
/* Remarks:                                                                     */
/*                                                                            */
/****************************************************************************/


/****************************************************************************/
/*                                                                            */
/* include files                                                            */
/* system include files                                                            */
/* application include files                                                     */
/*                                                                            */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include <amg_header.h>
#include <amg_low.h>
#include <amg_sp.h>
#include <amg_blas.h>
#include <amg_iter.h>
#include <amg_coarsen.h>

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/*   compile time constants defining static data size (i.e. arrays)            */
/*   other constants                                                            */
/*   macros                                                                    */
/*                                                                            */
/****************************************************************************/

#undef         DEBUG
#undef         DEBUG_CLUSTERING
#undef        DEBUG_STACK

/* flags for the node array */
/* maximum number of components is currently 7 ! */
#define ISOLATED(c)                ((c)&(0x1))            /* isolated node */
#define SET_ISOLATED(c)                (c) = (c)|(0x1)
#define RESET_ISOLATED(c)             (c) = (c)&(~(0x1))

#define VISITED(c)                ((c)&(0x80))           /* visited node */
#define SET_VISITED(c)                (c) = (c)|(0x80)
#define RESET_VISITED(c)             (c) = (c)&(0x7F)

#define FRONT(c)                ((c)&(0x40))           /* front node */
#define SET_FRONT(c)                (c) = (c)|(0x40)
#define RESET_FRONT(c)             (c) = (c)&(0xBF)

/* flags for link array */
#define DEPENDS(c)                ((c)&(0x01))          
#define SET_DEPENDS(c)                (c) = (c)|(0x01)
#define RESET_DEPENDS(c)             (c) = (c)&(0xFE)

#define INFLUENCES(c)                ((c)&(0x02))
#define SET_INFLUENCES(c)        (c) = (c)|(0x02)
#define RESET_INFLUENCES(c)     (c) = (c)&(0xFD)

#define ONEWAY(c)                (((c)&0x03)==1)      /* depends */
#define TWOWAY(c)                (((c)&0x03)==3)      /* depends & influences */
#define STRONG(c)                ((c)|0x03)           /* depends | influences */

/****************************************************************************/
/*                                                                            */
/* data structures used in this source file (exported data structures are   */
/*                  in the corresponding include file!)                            */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* definition of exported global variables                                    */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* definition of variables global to this source file only (static!)            */
/*                                                                            */
/****************************************************************************/

/* data for CVS */
static int clustersize=0;             /* number of nodes in cluster */
static int cluster[AMG_MAX_CLUSTER]; /* current cluster        */
static int frontsize=0;                     /* number of nodes in front */
static int front[AMG_MAX_FRONT];     /* current front */
static int stacksize=0;                     /* number of nodes in seed stack */
static int stackhead=0;                     /* push/pop entry */
static int stack[AMG_MAX_STACK];     /* current stack */
static int connectsize=0;             /* number of neighbors */
static int connect[AMG_MAX_ROW];     /* connectivity list */
static char buf[256];

/****************************************************************************/
/*                                                                            */
/* forward declarations of functions used before they are defined            */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*D
   NewGraph - allocate and initialize graph structure
   
   SYNOPSIS:
   static AMG_GRAPH *NewGraph (AMG_MATRIX *A)
   
   PARAMETERS:
.  A - matrix
   
   DESCRIPTION:
   This function allocates and initializes a graph data structure from
   a matrix.
      
   RETURN VALUE:
.n NULL memory overflow
.n else valid pointer
   
D*/
/****************************************************************************/

static AMG_GRAPH *NewGraph (AMG_MATRIX *A)
{
  AMG_GRAPH *new;
  int i,n,e,*ca;
  char *na,*la;
  float *da;
        
  
  /* allocate matrix structure */
  new = AMG_Malloc(sizeof(AMG_GRAPH));
  if (new==NULL) return(AMG_NULL);
  
  n = AMG_MATRIX_N(A);
  e = AMG_MATRIX_NONZEROS(A);
  
  /* allocate cluster array */
  ca = AMG_Malloc(n*sizeof(int));
  if (ca==NULL) return(AMG_NULL);
  
  /* allocate node array */
  na = AMG_Malloc(n*sizeof(char));
  if (na==NULL) return(AMG_NULL);
  
  /* allocate link array */
  la = AMG_Malloc(e*sizeof(char));
  if (la==NULL) return(AMG_NULL);
  
  /* allocate damping array */
  da = AMG_Malloc(n*sizeof(float));
  if (da==NULL) return(AMG_NULL);

  /* fill data structure */
  new->n = n;
  new->e = e;
  new->ra = AMG_MATRIX_RA(A);
  new->ja = AMG_MATRIX_JA(A);
  new->ca = ca;
  new->na = na;
  new->la = la;
  new->da = da;
  new->system_as_scalar = AMG_MATRIX_SAS(A);
  
  /* initialize data structure */
  for (i=0; i<n; i++) ca[i] = -1;
  for (i=0; i<n; i++) na[i] = 0;
  for (i=0; i<n; i++) da[i] = 1.0;
  for (i=0; i<e; i++) la[i] = 0;
  
  return(new);
}

/****************************************************************************/
/*
   Dependency - compute dependency graph and isolated nodes
   
   SYNOPSIS:
   static int Dependency (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc);
   
   PARAMETERS:
.  A - matrix
.  g - graph associated with matrix
.  cc - coarsen context parameters
   
   DESCRIPTION:
   This function computes the strong dependence graph and the isolated
   nodes. Corresponding flags of this process are left in the na and la arrays:
   
   ISOLATED(na[i],c) is true if component c in node i is isolated
   DEPENDS(la[k]) is true if node i depends strongly on node j and k is the
   index of this link.
   INFLUENCES(la[k]) is true if node i influences node j strongly.
   
   RETURN VALUE:
.n -1   fatal error
.n 0    no strong connections, do not coarsen
.n >0   number of strong connections
   
*/
/****************************************************************************/

static int Dependency (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
  int n = AMG_GRAPH_N(g);          /* number of nodes */
  int e = AMG_GRAPH_E(g);          /* number of edges/nonzeros */
  int *ra = AMG_GRAPH_RA(g);     
  int *ja = AMG_GRAPH_JA(g);
  int *ca = AMG_GRAPH_CA(g);       /* connectivity array */
  char *na = AMG_GRAPH_NA(g);
  char *la = AMG_GRAPH_LA(g);
  double *a = AMG_MATRIX_A(A);
  char buffer[128];
  int comp,sas=AMG_GRAPH_SAS(g);   /* sas - system as scalar */
  
  int i,k,kk,start,end;
  double maxval;
  int strong=0;
                  
  for (i=0; i<n; i++) na[i] = 0;   /* initialize flags, graph node array */
  for (i=0; i<e; i++) la[i] = 0;   /* graph edge array */
          
  for (i=0; i<n; i++)              /* for all rows */
    {
      start = ra[i];               /* i-th row start */
      end = start+ja[start];       /* i-th row end */
      comp=i%sas;                  /* important only for block systems */

      maxval = -1.0E20;            /* compute maximum offdiagonal strength */
      for (k=start+1; k<end; k++)  /* go through row i */
        if ((ja[k]%sas)==comp) maxval = AMG_MAX(maxval,-a[k]); /* this is one possibility */
      
      if (maxval<cc->beta*a[start])/* check for isolated condition */
        {
          SET_ISOLATED(na[i]);     /* component 0, the scalar case */
          continue;                /* next row */
        }      
                                   /* non-isolated node, set flags */
      maxval = maxval*cc->alpha;   /* reduce by alpha */
      for (k=start+1; k<end; k++)  /* go through row i */
        if ( ((ja[k]%sas)==comp) && (-a[k]>=maxval) )
          {
            SET_DEPENDS(la[k]);    /* i depends on ja[k] */
            /* original : kk = AMG_FindEntry(A,ja[k],i);*/
            /*kk = AMG_FindEntry(A,i,ja[k]);*/ /* changed 98/02/22, find a(i,ja[k]) */
            /*if (kk<0) 
              {
                sprintf(buffer,"%12s: no entry (%d,%d) 00\n","Dependency",ja[k],i);
                AMG_Print(buffer);
                return(-1);
              }
            printf("k %d kk %d\n",k,kk);*/
            SET_INFLUENCES(la[k]); /* ja[k] is influenced by i */
            strong++;
          }
    }

  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: %d rows %d nonzeros %d strong connections\n","Dependency",
              n,AMG_MATRIX_CONNECTIONS(A),strong);
      AMG_Print(buffer);
    }
  
  if (cc->verbose==1)
    {
      AMG_Print("U");
    }
  
  return(strong);
}

/****************************************************************************/
/*
   Dependency_sym - compute dependency graph and isolated nodes for matrices
                    with symmetric pattern
                    preferably the matrix should be symmetric itself
   
   SYNOPSIS:
   static int Dependency_sym (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc);
   
   PARAMETERS:
.  A - matrix
.  g - graph associated with matrix
.  cc - coarsen context parameters
   
   DESCRIPTION:
   This function computes the strong dependence graph and the isolated
   nodes. Corresponding flags of this process are left in the na and la arrays:
   
   ISOLATED(na[i],c) is true if component c in node i is isolated
   DEPENDS(la[k]) is true if node i depends strongly on node j and k is the
   index of this link.
   INFLUENCES(la[k]) is true if node i influences node j strongly.
   
   RETURN VALUE:
.n -1   fatal error
.n 0    no strong connections, do not coarsen
.n >0   number of strong connections
   
*/
/****************************************************************************/
static int Dependency_sym (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
  int n = AMG_GRAPH_N(g);
  int e = AMG_GRAPH_E(g);
  int *ra = AMG_GRAPH_RA(g);
  int *ja = AMG_GRAPH_JA(g);
  int *ca = AMG_GRAPH_CA(g);
  char *na = AMG_GRAPH_NA(g);
  char *la = AMG_GRAPH_LA(g);
  double *a = AMG_MATRIX_A(A);
  char buffer[128];
  int comp,sas=AMG_GRAPH_SAS(g);
  
  int i,k,kk,start,end;
  double maxval;
  int strong=0;

  for (i=0; i<n; i++) na[i] = 0;  /* initialize flags, graph node array */
  for (i=0; i<e; i++) la[i] = 0;  /* graph edge array */

  for (i=0; i<n; i++)             /* for all rows */
    {
      start = ra[i];              /* i-th row start */
      end = start+ja[start];      /* i-th row end */ 
      comp=i%sas;                 /* important only for block systems */
                
      maxval = -1.0E20;                  /* compute symmetric maxval */
      for (k=start+1; k<end; k++) /* go through row i */
        if ((ja[k]%sas)==comp)
          {
            kk = AMG_FindEntry(A,ja[k],i); /* find a(ja[k],i) */
            if (kk<0) 
              {
                sprintf(buffer,"%12s: no entry (%d,%d) 11\n","Dependency",ja[k],i);
                AMG_Print(buffer);
                sprintf(buffer,"WARNING : SET cc.dependency=AMG_UNSYM !!!");
                AMG_Print(buffer);
                return(-1);
              }                   /* (a_ij * a_ji)/(a_ii * a_jj) */
            maxval=AMG_MAX(maxval,(a[k]*a[kk])/(a[start]*a[ra[ja[k]]]));
          }
      
      if (maxval<cc->beta)        /* check for isolated condition */
        {
          SET_ISOLATED(na[i]);    /* component 0, the scalar case */
          continue;               /* next row */
        }
      
      maxval = maxval*cc->alpha;  /* reduce by alpha */
      for (k=start+1; k<end; k++) /* go through row i */
        if ((ja[k]%sas)==comp)
          {
            kk = AMG_FindEntry(A,ja[k],i);/* find a(ja[k],i) */
            if (kk<0) 
              {
                sprintf(buffer,"%12s: no entry (%d,%d) 22\n","Dependency",ja[k],i);
                AMG_Print(buffer);
                return(-1);
              }
            if ((a[k]*a[kk])/(a[start]*a[ra[ja[k]]]) < maxval) continue; /* not strong */
            SET_DEPENDS(la[k]);  SET_INFLUENCES(la[kk]); /* i depends on k */ 
            SET_DEPENDS(la[kk]); SET_INFLUENCES(la[k]);  /* k depends on i */ 
            strong++;
          }
    }
  
  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: %d rows %d nonzeros %d strong connections\n","Dependency",
              n,AMG_MATRIX_CONNECTIONS(A),strong);
      AMG_Print(buffer);
    }
  
  if (cc->verbose==1)
    {
      AMG_Print("S");
    }
  
  return(strong);
}

/****************************************************************************/
/*D
   Clustering - cluster nodes
   
   SYNOPSIS:
   static int Clustering (AMG_GRAPH *g, AMG_CoarsenContext *cc);
   
   PARAMETERS:
.  g - the graph structure
.  cc - coarsen context structure
   
   DESCRIPTION:
   This function uses the flags set by dependency to partition the
   nodes into clusters.
   
   RETURN VALUE:
.n AMG_OK if no error
   
D*/
/****************************************************************************/

/* add node to cluster */
static int AddCluster (AMG_GRAPH *g, int s, int clusternumber)
{
  int k,m,n,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  int sas=AMG_GRAPH_SAS(g);
  
  if (clustersize>=AMG_MAX_CLUSTER)
    {
      AMG_Print("AddCluster : clustersize greater than AMG_MAX_CLUSTER!!!\n");
      exit(4711);
    }
  cluster[clustersize++] = s;      /* add node to cluster */
  ca[s] = clusternumber;           /* asign to node the current cluster number */
  if (s%sas != cluster[0]%sas)
        {
          AMG_Print("components of cluster inconsistent!\n");
          return(AMG_FATAL);
        }
  
                                   /* extend connectivity list of cluster, include -1 */
                                   /* check consistency of components */
  start = ra[s]; end = start+ja[start]; /* go through row of new node */
  for (k=start+1; k<end; k++)          
    {
      n=ca[ja[k]];                      /* clusternumber of node in this row */
      for (m=0; m<connectsize; m++)     /* check if there is already a connection */ 
        if (connect[m]==n) break;       /* to this cluster: got it, break */
      if (m==connectsize)               /* no connection found */
        {                               
          if (connectsize>=AMG_MAX_ROW) /* too much connections */
            {
              AMG_Print("AddCluster : connectsize greater than AMG_MAX_ROW!!!\n");
              exit(4711);
            }
          connect[connectsize++]=n;     /* add new connection */
        }
    }
  
  return(AMG_OK);
}

/* initialize cluster with seed node */
static int SeedCluster (AMG_GRAPH *g, int s, int clusternumber)
{
        int *ca=g->ca;

        clustersize = 0;
        connectsize=1;
        connect[0] = clusternumber;
        return(AddCluster(g,s,clusternumber));
}

/* compute "front" nodes */
static int SetFront (AMG_GRAPH *g)
{
  int i,k,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *na=g->na, *la=g->la;
  int comp,sas=AMG_GRAPH_SAS(g);
  
  frontsize = 0;
  for (i=0; i<clustersize; i++)
    {
      start = ra[cluster[i]]; end = start+ja[start]; comp=cluster[i]%sas;
      for (k=start+1; k<end; k++) 
        if ((ja[k]%sas==comp) && (ca[ja[k]]<0) && (!FRONT(na[ja[k]]))) /* not in any cluster and not in front */
          {
            if (frontsize>=AMG_MAX_FRONT) return(AMG_FATAL);
            front[frontsize++] = ja[k];
            SET_FRONT(na[ja[k]]);
          }
    }
  return(AMG_OK);
}

/* clear the front */
static int ResetFront (AMG_GRAPH *g)
{
        int i;
        char *na=g->na;
        
        for (i=0; i<frontsize; i++)
                RESET_FRONT(na[front[i]]);
        frontsize=0;
        return(AMG_OK);
}

/* count number of two-way connections to the cluster */
static int TwoWay (AMG_GRAPH *g, int f, int clusternumber)
{
  int k,start,end,t;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *la=g->la;
  
  t = 0;
  start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++)
    {
      if (ca[ja[k]]!=clusternumber) continue;
      if (((la[k])&0x03)!=3) continue;
      t++;
    }
  return(t);
}

/* count number of one-way connections to the cluster */
static int OneWay (AMG_GRAPH *g, int f, int clusternumber)
{
  int k,start,end,o;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *la=g->la;
  
  o = 0;
  start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++)
    if (ca[ja[k]]==clusternumber && 
        (DEPENDS(la[k]) || INFLUENCES(la[k])) && (!TWOWAY(la[k]))) o++;
  return(o);
}

/* count number of connections to the cluster either way */
static int Connected (AMG_GRAPH *g, int f, int clusternumber)
{
  int k,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *la=g->la;
  
  start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++)
    if ( ca[ja[k]]==clusternumber  ) return(1);
  return(0);
}

/* count number of neighbors in front */
static int FrontNeighbors (AMG_GRAPH *g, int f)
{
  int k,start,end,w;
  int *ra=g->ra, *ja=g->ja;
  char *la=g->la, *na=g->na;
  
  w = 0;
  start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++)
    if (FRONT(na[ja[k]])) w++;
  return(w);
}

/* count number of strong connections to free nodes in same component */
static int UnusedNeighbors (AMG_GRAPH *g, int f)
{
  int k,start,end,u,l;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char  *la=g->la;
  /*int comp,sas=AMG_GRAPH_SAS(g);*/
  
  u = 0;
  start = ra[f]; 
  end = start+ja[start]; 
  /*comp=f%sas;*/
  for (k=start+1; k<end; k++)
    {
      /*l = ja[k];*/
      if (ca[ja[k]]>=0) continue; /* skip nodes already in clusters */
      /*if (l%sas!=comp) continue;*/ /* skip other components */
      l =la[k];
      if (DEPENDS(l)) u++;
      if (INFLUENCES(l)) u++;
    }
  return(u);
}

/* count the connectivity of the new cluster if f is added */
static int Connectivity (AMG_GRAPH *g, int f)
{
  int k,m,n,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  int c;
  
  c=0; start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++) /* look also at other components, this is important! */
    {
      n=ca[ja[k]];
      for (m=0; m<connectsize; m++)
        if (connect[m]==n) break; /* got it, connection already there */
      if (m==connectsize || n<0)  /* ja[k] not yet in a cluster */
        {
          c++; 
          continue;
        } 
      c+=2;                       /* new connection */
    }
  
  return(c);
}

/* count number of strong connections to cluster nodes */
/* since f is a front node it is from the same component */
static int ClusterNeighbors (AMG_GRAPH *g, int f, int clusternumber)
{
  int k,start,end,u;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char  *la=g->la;
  
  u = 0;
  start = ra[f]; end = start+ja[start];
  for (k=start+1; k<end; k++)
    {
      if (ca[ja[k]]!=clusternumber) continue; /* skip nodes already in clusters */
      if (DEPENDS(la[k])) u++;
      if (INFLUENCES(la[k])) u++;
    }
  return(u);
}

/* return a cluster number of a strongly connected non-isolated node of the same component */
static int InnerNeighborCluster (AMG_GRAPH *g, int f)
{
  int k,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char  *na=g->na, *la=g->la;
  int comp,sas=AMG_GRAPH_SAS(g);
  
  start = ra[f]; end = start+ja[start]; comp=f%sas;
  for (k=start+1; k<end; k++)
    if ( (ja[k]%sas==comp) && ca[ja[k]]>=0 && (!ISOLATED(na[ja[k]])) )
      return(ca[ja[k]]);
  
  return(-1);
}


/****************************************************************************/
/*                                                                            */
/* compute the ratio of                                                     */
/*      connections of the type "depend" between the nodes of the cluster   */
/* and                                                                            */
/*        connections of the type "strong" between the nodes of the cluster   */
/*                                                                            */
/****************************************************************************/
static float AutoDamp (AMG_GRAPH *g, int clusternumber)
{
  int f,k,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char  *na=g->na, *la=g->la;
  int strong=0,depends=0;
  float d;
  
  for (f=0; f<clustersize; f++)            /* for all nodes in the cluster */
    {
      start = ra[cluster[f]];              /* row of the node in fine matrix */
      end = start+ja[start]; 
      for (k=start+1; k<end; k++)
        {
          if ( ca[ja[k]]!=clusternumber ) continue;  /* neighbour not in cluster */
          if (DEPENDS(la[k])) depends++;             /* i depends on k */
          if (DEPENDS(la[k]) || INFLUENCES(la[k])) strong++; /* neighbour depends strongly */
        }
    }
  strong = strong/2;                        /* always counted twice */
  if (strong>0)
    d = ((float)depends) / ((float)strong); /* is between 1.0 and 2.0 */
  else
    d=1.0;
  
  return(d);
}

/* compute max distance to any cluster node */
static int Distance (AMG_GRAPH *g, int f, int clusternumber)
{
  int i,k,start,end,first,last,l;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca; 
  char *la=g->la, *na=g->na;
  int visitedsize=0;
  int visited[AMG_MAX_CLUSTER];
  int dist=0;
  
  if (clustersize>=AMG_MAX_CLUSTER) return(100000); /* cluster is full */
  visitedsize=1; visited[0]=f; first=0; last=1;     /* breadth first search */
  while (1)
    {
      for (i=first; i<last; i++)              
        {
          start = ra[visited[i]]; end = start+ja[start];
          for (k=start+1; k<end; k++)
            {
              l = ja[k];
              if ( (ca[l]==clusternumber) && (!VISITED(na[l])) )
                {                                   /* not yet visited */
                  SET_VISITED(na[l]);           /* set visited */
                  visited[visitedsize++] = l;   /* store number of visited node */
                }
            }
        }
      first=last; last=visitedsize;
      if (last>first) dist++;
      if (last==first) break;
      if (dist>1000) break;
    }
  for (i=1; i<visitedsize; i++) RESET_VISITED(na[visited[i]]);
  return(dist);
}


static int Admissible (AMG_GRAPH *g, int f, int clusternumber)
{
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *la=g->la, *na=g->na;
  int i,j,ki,kj,kij,starti,startj,startij,endi,endj,endij,found;
  
  /* situation 1, front node depends on two nodes in cluster */        
  starti = ra[f]; endi = starti+ja[starti];
  for (ki=starti+1; ki<endi; ki++)
    if (DEPENDS(la[ki]) && !INFLUENCES(la[ki]) && ca[ja[ki]]==clusternumber)
      {
        i = ja[ki];                /* f depends one-way on i */
        startj = ra[f]; endj = startj+ja[startj];
        for (kj=startj+1; kj<endj; kj++)
          if ((ja[kj]>i) && DEPENDS(la[kj]) && !INFLUENCES(la[kj]) && ca[ja[kj]]==clusternumber)
            {
              j = ja[kj];          /* f depends one-way on j neq i */
              
                                   /* check if i and j are strongly connected */
              startij = ra[i]; endij = startij+ja[startij]; found=0;
              for (kij=startij+1; kij<endij; kij++)
                if (ja[kij]==j && STRONG(la[kij]))
                  {
                    found=1; break;
                  }
              if (!found) return(0);
            }
      }
  
  /* situation 2, cluster node depends on front node and other cluster node */
  starti = ra[f]; endi = starti+ja[starti];
  for (ki=starti+1; ki<endi; ki++)
    if (!DEPENDS(la[ki]) && INFLUENCES(la[ki]) && ca[ja[ki]]==clusternumber)
      {
        i = ja[ki];               /* i depends one-way on f */
        startj = ra[i]; endj = startj+ja[startj]; /* look at neighbors of i */
        for (kj=startj+1; kj<endj; kj++)
          if (ja[kj]!=f && DEPENDS(la[kj]) && !INFLUENCES(la[kj]) && ca[ja[kj]]==clusternumber)
            {
              j = ja[kj];         /* i depends one-way on j also in cluster (f not in cluster) */
              
                                  /* check if f and j are strongly connected */
              startij = ra[f]; endij = startij+ja[startij]; found=0;
              for (kij=startij+1; kij<endij; kij++)
                if (ja[kij]==j && STRONG(la[kij]))
                  {
                    found=1; break;
                  }
              if (!found) return(0);
            }
      }
  
  return(1);
}

static int ClearStack (void)
{
        stacksize=0;
        stackhead=0;
        return(AMG_OK);
}

static int PopMajor (AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
        int *ca=g->ca;
        char *na=g->na;
        int i,umin,isoumin,majisoumin,majumin;
        int n = AMG_GRAPH_N(g);
        int iso,con;
    int isomaj,conmaj,comp,sas=AMG_GRAPH_SAS(g),major=cc->major;

        if (major<0) return(-1);
        
        /* get from stack */
        while (stacksize>0)
        {
                stackhead = (stackhead+AMG_MAX_STACK-1)%AMG_MAX_STACK;
                stacksize--;
                i = stack[stackhead];
                if (ca[i]<0) return(i);
        }
        
        /* stack is empty, are there any nodes left to fill the stack ? */
        stackhead=0; stacksize=0;
        majisoumin=majumin=isoumin=umin=1000000;
        isomaj=conmaj=iso=con=0;
        for (i=0; i<n; i++)
        {
                if ( ca[i]>=0 ) continue;
                if ( ISOLATED(na[i]) )
                {
                        isoumin =  AMG_MIN(isoumin,UnusedNeighbors(g,i));
                        iso++;
                        if (major>=0 && sas>1 && (i%sas==major) ) 
                        {
                                isomaj++;
                                majisoumin = AMG_MIN(majisoumin,UnusedNeighbors(g,i));
                        }
                }
                else
                {
                        umin = AMG_MIN(umin,UnusedNeighbors(g,i));
                        con++;
                        if (major>=0 && sas>1 && (i%sas==major) ) 
                        {
                                conmaj++;
                                majumin = AMG_MIN(majumin,UnusedNeighbors(g,i));
                        }
                }
        }
        if (iso+con==0) return(-1); /* early exit */

        #ifdef DEBUG_STACK
        sprintf(buf,"RESTACKING n=%d iso=%d con=%d\n",n,iso,con);
        AMG_Print(buf);
        for (i=0; i<n; i++)
        {
                if (ca[i]>=0) continue;
            sprintf(buf,"i=%d, ca[i]=%d, isolated=%d, comp=%d\n",
                        i,ca[i],ISOLATED(na[i]),i%sas);
            AMG_Print(buf);
        }
        #endif

        if (conmaj+isomaj==0)
        {
                AMG_Print("PopMajor: only nonmajor nodes left?! \n");
                return(-1);
        }

        if (conmaj>0)
        {
                /* there are major components in major mode ... */
        for (i=0; i<n; i++)
            if ( (i%sas==major) && ca[i]<0 && (!ISOLATED(na[i])) && UnusedNeighbors(g,i)==majumin )
            {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
            }
        #ifdef DEBUG_STACK
        sprintf(buf,"restacking connected size=%d\n",stacksize);
        AMG_Print(buf);
        #endif
        }
        else
        {
                /* there are major components in major mode ... */
        for (i=0; i<n; i++)
            if ( (i%sas==major) && ca[i]<0 && (ISOLATED(na[i])) && UnusedNeighbors(g,i)==majisoumin )
            {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
            }
        #ifdef DEBUG_STACK
        sprintf(buf,"restacking connected size=%d\n",stacksize);
        AMG_Print(buf);
        #endif
        }
        
        /* now try again */
        while (stacksize>0)
        {
                stackhead = (stackhead+AMG_MAX_STACK-1)%AMG_MAX_STACK;
                stacksize--;
                i = stack[stackhead];
                if (ca[i]<0) return(i);
        }

        /* no, all nodes have been processed return invalid node number */
        return(-1);
}

static int Pop (AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
  int *ca=g->ca;
  char *na=g->na;
  int i,umin,isoumin,majisoumin,majumin;
  int n = AMG_GRAPH_N(g);
  int iso,con;
  int isomaj,conmaj,comp,sas=AMG_GRAPH_SAS(g),major=cc->major;
  
  /* get from stack */
  while (stacksize>0)
    {
      stackhead = (stackhead+AMG_MAX_STACK-1)%AMG_MAX_STACK;
      stacksize--;
      i = stack[stackhead];
      if (ca[i]<0) return(i);
    }
        
  /* stack is empty, are there any nodes left to fill the stack ? */
  stackhead=0; stacksize=0;
  majisoumin=majumin=isoumin=umin=1000000;
  isomaj=conmaj=iso=con=0;
  for (i=0; i<n; i++)
    {
      if ( ca[i]>=0 ) continue;
      if ( ISOLATED(na[i]) )
        {
          if (con) continue;
          isoumin =  AMG_MIN(isoumin,UnusedNeighbors(g,i));
          iso++;
          /*if (major>=0 && sas>1 && (i%sas==major) ) 
            {
              isomaj++;
              majisoumin = AMG_MIN(majisoumin,UnusedNeighbors(g,i));
              }*/
        }
      else
        {
          umin = AMG_MIN(umin,UnusedNeighbors(g,i));
          con++;
          /*if (major>=0 && sas>1 && (i%sas==major) ) 
            {
              conmaj++;
              majumin = AMG_MIN(majumin,UnusedNeighbors(g,i));
              }*/
        }
    }
  if (iso+con==0) return(-1); /* early exit */
#ifdef DEBUG_STACK
  sprintf(buf,"RESTACKING n=%d iso=%d con=%d\n",n,iso,con);
  AMG_Print(buf);
  for (i=0; i<n; i++)
    {
      if (ca[i]>=0) continue;
      sprintf(buf,"i=%d, ca[i]=%d, isolated=%d, comp=%d\n",
              i,ca[i],ISOLATED(na[i]),i%sas);
      AMG_Print(buf);
    }
#endif
  if (con>0)
    {
      if (conmaj==0)
        {
          /*if (major>=0) AMG_Print("Pop: pushing connected nonmajor nodes?!\n");*/
          /* push only connected */
          for (i=0; i<n; i++)
            if ( ca[i]<0 && (!ISOLATED(na[i])) && UnusedNeighbors(g,i)==umin )
              {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
              }
#ifdef DEBUG_STACK
          sprintf(buf,"restacking connected size=%d\n",stacksize);
          AMG_Print(buf);
#endif
        }
      else
        {
          /* there are major components in major mode ... */
          for (i=0; i<n; i++)
            /*if ( (i%sas==major) && ca[i]<0 && (!ISOLATED(na[i])) && UnusedNeighbors(g,i)==majumin )*/
            if ( ca[i]<0 && (!ISOLATED(na[i])) && UnusedNeighbors(g,i)==majumin )
              {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
              }
#ifdef DEBUG_STACK
          sprintf(buf,"restacking connected size=%d\n",stacksize);
          AMG_Print(buf);
#endif
        }
    }
  else
    {
      if (isomaj==0)
        {
          /*if (major>=0) AMG_Print("Pop: pushing isolated nonmajor nodes?!\n");*/
          /* push only isolated */
          for (i=0; i<n; i++)
            if ( ca[i]<0 && (ISOLATED(na[i])) && UnusedNeighbors(g,i)==isoumin )
              {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
              }
#ifdef DEBUG_STACK
          sprintf(buf,"restacking isolated size=%d\n",stacksize);
          AMG_Print(buf);
#endif
        }
      else
        {
          /* there are major components in major mode ... */
          for (i=0; i<n; i++)
            /*if ( (i%sas==major) && ca[i]<0 && (ISOLATED(na[i])) && UnusedNeighbors(g,i)==majisoumin )*/
            if ( ca[i]<0 && (ISOLATED(na[i])) && UnusedNeighbors(g,i)==majisoumin )
              {
                stack[stackhead] = i;
                stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
                stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
              }
#ifdef DEBUG_STACK
          sprintf(buf,"restacking connected size=%d\n",stacksize);
          AMG_Print(buf);
#endif
        }
    }
  
  /* now try again */
  while (stacksize>0)
    {
      stackhead = (stackhead+AMG_MAX_STACK-1)%AMG_MAX_STACK;
      stacksize--;
      i = stack[stackhead];
      if (ca[i]<0) return(i);
    }
  
  /* no, all nodes have been processed return invalid node number */
  return(-1);
}

static int Push (AMG_GRAPH *g, int i, AMG_CoarsenContext *cc)
{
        int *ca=g->ca;
        int sas=g->system_as_scalar, major=cc->major;
        char buffer[128];

        if ( ca[i]>=0 ) return(0);

        if (major>=0 && i%sas!=major)
        {
                sprintf(buffer,"pushing nonmajor component %d\n",i);
                AMG_Print(buffer);
        }

        stack[stackhead] = i;
        stacksize=AMG_MIN(stacksize+1,AMG_MAX_STACK);
        stackhead=(stackhead+AMG_MAX_STACK+1)%AMG_MAX_STACK;
        
        return(1);
}

/* find a neighbor with some cluster, same ISO, and strong connection */
static int MergeNeighbor (AMG_GRAPH *g, int f)
{
  int k,start,end,j,startj,endj;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char  *la=g->la, *na=g->na;
  int comp,sas=AMG_GRAPH_SAS(g);
  
  start = ra[f]; end = start+ja[start]; comp=f%sas;
  for (k=start+1; k<end; k++)
    {
      if (ja[k]%sas!=comp) continue;        /* skip other components */
      if (ca[ja[k]]<0) continue;            /* must be a valid cluster */
      if ((ISOLATED(na[ja[k]])) != (ISOLATED(na[f]))) continue;
      if (!ISOLATED(na[f]))                 /* node not isolated */
        {          
          if ( (DEPENDS(la[k]) || INFLUENCES(la[k])) && Admissible(g,f,ca[ja[k]]) )
            {
              startj=ra[ja[k]];
              endj=startj+ja[startj];
              for (j=startj+1;j<endj;j++)
                if (ja[j]==f)
                  return(ja[k]);
            }
        }
      else                                  /* node f isolated */ 
        {
          return(ja[k]);                    /* node ja[k] isolated too */
        }
    }
  return(-1);
}

/****************************************************************************/
/*                                                                            */
/* build small clusters of nodes which are as far as possible strongly      */
/* connected                                                                    */
/*                                                                            */
/****************************************************************************/
static int Clustering (AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
  int m,s,f,i,n=g->n;
  int *ca=g->ca;
  char *na=g->na;
  int C,D,T,O,W,N,Cmax,Nmax,Tmax,Omax,Wmax,candidate,c;
  int clusternumber,conclusters,isoclusters;
  char buffer[128];
  int comp,sas=AMG_GRAPH_SAS(g),major=cc->major;
  int noniso,iso;
  clock_t start_coarsen,finish_coarsen;  
                        
  /* initialization */
  for (i=0; i<n; i++) {
    ca[i]=-1; 
    RESET_FRONT(na[i]); 
    RESET_VISITED(na[i]);
  }
  ClearStack();
  clusternumber=0;
  conclusters=isoclusters=0;
  
  /* greedy search */
  while ( 1 )
    {
      if (major>=0)             /* find seed for a new cluster */
        s=PopMajor(g,cc); 
      else 
        s=Pop(g,cc);            /* scalar case */
      if (s<0) break;           /* all nodes are in clusters */
      
      if (cc->verbose==1 && clusternumber>0 && (clusternumber%10000==0) )
        {
          AMG_Print("c");
        }
        
      /* grow the cluster */
      if (ISOLATED(na[s]))                      /* the isolated case */
        {                                       /* s is the seed of our cluster */
          SeedCluster(g,s,clusternumber);        /* initialize new cluster and add seed to it */
          c=InnerNeighborCluster(g,s);          /* get number of cluster where a non-iso neighbor */
                                                /* of s belongs to or -1 */
          while (clustersize<cc->mincluster && connectsize<cc->maxconnectivity)
            {
              SetFront(g); /* init front, find node which are not in any cluster and not in front */
              Wmax=-1; 
              Cmax=-100000; 
              candidate=-1;                     /* no candidate yet */
              for (i=0; i<frontsize; i++)       /* go through front */
                {
                  f = front[i];                 /* examine node f from front */
                  if (!ISOLATED(na[f])) continue; /* only isolated nodes here */
                  D=Distance(g,f,clusternumber);
                  if (D>cc->maxdistance) continue; /* too far away */
                  
                  if (c>=0)                     /* if there is already a neighbor cluster */
                    {
                                                /* front node  must be connected to same inner cluster */
                                                /* if it is connected to any inner cluster at all      */
                      if ( (InnerNeighborCluster(g,f)>=0) && (!Connected(g,f,c)) ) 
                        continue;  
                    }
                  C=Connectivity(g,f);

                  if (C==Cmax)
                    {
                      W=FrontNeighbors(g,f);
                      if (W>Wmax)
                        {
                          Wmax=W;  
                          candidate=f;
                        }
                    }
                  if (C>Cmax)
                    {
                      Cmax=C; 
                      Wmax=FrontNeighbors(g,f); 
                      candidate=f;
                    }
                }                                 /* all front nodes examined */
              ResetFront(g);                          /* we do not need the front anymore */
              if (candidate<0) break;             /* no more candidates found */                                
              if (AddCluster(g,candidate,clusternumber)!=AMG_OK) return(AMG_FATAL); /* extend cluster */

                                                  /* fill seed stack with front */
              SetFront(g); 
              if (major>=0 && (s%sas==major))
                {
                  /* push only major components */
                  for (i=0; i<frontsize; i++)
                    if ( front[i]%sas==major && (ISOLATED(na[front[i]])) )
                      Push(g,front[i],cc);
                }
              else
                {
                  for (i=0; i<frontsize; i++)
                    if (ISOLATED(na[front[i]])) Push(g,front[i],cc);
                }
              ResetFront(g);
              
                                /* let's see if we can find a neighboring non-iso cluster */
              if (c<0) c=InnerNeighborCluster(g,candidate);
            }
        }       /* isolated case done */
      else        /* the interior node case */
        {
                                                /* s is the seed of our cluster */
          SeedCluster(g,s,clusternumber);        /* initialize new cluster and add seed to it */
                
          while (clustersize<cc->mincluster)    /* up to mincluster size */
            {
                                                /* now comes the big one step logic */
              SetFront(g);  /* init front, find node which are not in any cluster and not in front */
                                
                                                /* find a candidate to add to the cluster */
              Tmax=Omax=0; 
              Wmax=-1;                          /* maximal number of neighbours in front */
              Cmax=-100000;                     /* maximal connectivity */
              candidate=-1;
              for (i=0; i<frontsize; i++)       /* go through front */
                {
                  f = front[i];                 /* examine node f from front */
                  if (ISOLATED(na[f])) continue;/* no isolated nodes here */
                  D=Distance(g,f,clusternumber);/* maximal distance of f to cluster nodes */
                  if (D>cc->maxdistance) continue; /* too far away */
                                        
                  /* two-way case, front node depends on cluster nodes and influences cluster nodes */
                  T=TwoWay(g,f,clusternumber);   /* number of two-way connections */ 
                  if (T==Tmax && T>0)            /* equal to current maximal number of two-way connections */ 
                    {
                      C=Connectivity(g,f);       /* connectivity to other clusters if f is added */
                      if (C==Cmax)               /* equal to current maximal connectivity */
                        {
                          W=FrontNeighbors(g,f); /* number of neighbours in front */
                          if (W>Wmax)            /* new maximal number of neighbours in front */
                            {
                              Wmax=W;  
                              candidate=f;       /* f is the current candidate */
                            }
                        }
                      if (C>Cmax)                /* new maximal connectivity */
                        {
                          Cmax=C; 
                          Wmax=FrontNeighbors(g,f); 
                          candidate=f;           /* f is the current candidate */
                        }
                    }
                  if (T>Tmax)                    /* new maximal number of two-way connections */ 
                    {
                      Tmax=T; 
                      Wmax=FrontNeighbors(g,f); 
                      Cmax=Connectivity(g,f);
                      candidate=f;               /* f is the current candidate */
                      Omax=100000;               /* Two way connections preceed one way conn. */
                    }
                  if (T>0) continue;             /* this was a two-way node, skip the rest */
                  if (Omax==100000) continue;    /* there was already a two-way node */
                                                 /* previous line added 98/12/11 V.J. */     
                  /* one-way case */
                  O=OneWay(g,f,clusternumber);   /* number of one way connections */
                  if (O>0 && !Admissible(g,f,clusternumber)) continue; /* must be admissible */ 
                  if (O==Omax && O>0)
                    {
                      C=Connectivity(g,f);      /* connectivity to other clusters if f is added */
                      if (C==Cmax)              /* equal to current maximal connectivity */
                        {
                          W=FrontNeighbors(g,f);/* number of neighbours in front */
                          if (W>Wmax)           /* new maximal number of neighbours in front */
                            {
                              Wmax=W;  
                              candidate=f;      /* f is the current candidate */
                            }
                        }
                      if (C>Cmax)               /* new maximal connectivity */
                        {
                          Cmax=C; 
                          Wmax=FrontNeighbors(g,f); 
                          candidate=f;          /* f is the current candidate */
                        }
                    }
                  if (O>Omax)                   /* new maximal number of one-way connections */ 
                    {
                      Omax=O; 
                      Wmax=FrontNeighbors(g,f); 
                      Cmax=Connectivity(g,f);
                      candidate=f;              /* f is the current candidate */
                    }
                }                               /* all front nodes examined */
                                
              ResetFront(g);                        /* we do not need the front anymore */
              if (candidate<0) break;           /* no more candidates found */
              if (AddCluster(g,candidate,clusternumber)!=AMG_OK) return(AMG_FATAL); /* extend cluster */
                                                
                                                /* fill seed stack with front */
              SetFront(g);                      /* compute new front nodes */
              if (major>=0 && (s%sas==major))
                {
                  /* push only major components */
                  for (i=0; i<frontsize; i++)
                    if ( front[i]%sas==major && (!ISOLATED(na[front[i]])) )
                      Push(g,front[i],cc);
                }
              else
                {
                  for (i=0; i<frontsize; i++)   /* put non-iso node in stack if it doesn't belong */
                    if (!ISOLATED(na[front[i]])) Push(g,front[i],cc); /* to any cluster */
                }
              ResetFront(g);
            }    /* now clustersize==mincluster or no candidates found */

          while (clustersize<cc->maxcluster)    /* "rounding step" */
            {                                        /* now comes the big one step logic */
              SetFront(g);                      /* init front */
              candidate=-1;                     /* find a candidate to add to the cluster */
              for (i=0; i<frontsize; i++)
                {
                  f = front[i];                 /* examine node f from front */
                  if (ISOLATED(na[f])) continue;/* no isolated nodes here */
                  D=Distance(g,f,clusternumber);
                  if (D>cc->maxdistance) continue; /* too far away */
                  
                  T=TwoWay(g,f,clusternumber);  /* number of two-way connections */ 
                  if (T==0)                     /* no two-way connections */ 
                    {
                      O=OneWay(g,f,clusternumber);/* number of one-way connections */ 
                      if (O==0) continue;         /* no one-way connections */ 
                      if (!Admissible(g,f,clusternumber)) continue;/* not admissible */
                    }
                  /* number of strong connections to nodes in the cluster must be larger than */
                  /* number of strong connections to free nodes */
                  if (ClusterNeighbors(g,f,clusternumber)<=UnusedNeighbors(g,f))
                    continue;       
                  candidate=f;                  /* candidate found */
                  break;
                }
              ResetFront(g);                        /* we do not need the front anymore */
              if (candidate<0) break;           /* no more candidates found */
              if (AddCluster(g,candidate,clusternumber)!=AMG_OK) return(AMG_FATAL); /* extend cluster */
              
                                                /* fill seed stack with front */
              SetFront(g); 
              if (major>=0 && (s%sas==major))
                {
                  /* push only major components */
                  for (i=0; i<frontsize; i++)
                    if ( front[i]%sas==major && (!ISOLATED(na[front[i]])) )
                      Push(g,front[i],cc);
                }
              else
                {
                  for (i=0; i<frontsize; i++)
                    if (!ISOLATED(na[front[i]])) Push(g,front[i],cc);
                }
              ResetFront(g);
            }   /* clustersize==maxcluster or no more candidates */

        }       /* interior node case finished */

      if ((clustersize==1)&&(!ISOLATED(na[s]))) /* second condtion added 98/12/11 V.J. */ 
        m=MergeNeighbor(g,s);                   /* only one node in cluster, find neighbour */
      else                                      /* with strong connections to that node, same iso ... */
        m=-1;
      if ((clustersize==1) && (m>=0) && (!ISOLATED(na[s])) )
        {
          /* merge */
          ca[s]=ca[m];                          /* put s to cluster m */
#ifdef DEBUG_CLUSTERING
          sprintf(buf,"Merge   node %d: %d\n",s,ca[m]); AMG_Print(buf);
#endif
          
          /* major mode for systems: use same clustering for all components */
          if (major>=0)
            {
              if (s%sas!=major)
                {
                  AMG_Print("clustered a nonmajor isolated component !?\n");
                  sprintf(buffer,"clusternumber=%d, s=%d, major=%d, sas=%d\n",
                          clusternumber,s,major,sas);
                  AMG_Print(buffer);
                  for (i=0; i<n; i++)
                    {
                      if (ca[i]>=0) continue;
                      sprintf(buffer,"%d: %d\n",i,ca[i]);
                      AMG_Print(buffer);
                    }
                  AMG_Print(buffer);
                  return(AMG_FATAL);
                }
              
                                /* major strategy: use same clustering for all components */
              for (comp=0; comp<sas; comp++)
                {
                  if (comp==major) continue;
                  ca[(s/sas)*sas+comp]=ca[(m/sas)*sas+comp];
#ifdef DEBUG_CLUSTERING
                  sprintf(buf,"MMerge  node %d: %d\n",(cluster[0]/sas)*sas+comp,ca[(m/sas)*sas+comp]); AMG_Print(buf);
#endif
                }
            }
          
        }
      else    /* clustersize > 1 || s is isolated || no neighbour for s to merge */
        {                
#ifdef DEBUG_CLUSTERING
          sprintf(buf,"Cluster %4d: ",clusternumber); AMG_Print(buf);
          for (i=0; i<clustersize; i++)
            {
              sprintf(buf,"%4d ",cluster[i]); AMG_Print(buf);
            }
          AMG_Print("\n");
#endif
          
          clusternumber++;              /* increment cluster number */
          if (ISOLATED(na[s])) 
            isoclusters++; 
          else 
            conclusters++;
          
          /* major mode for systems: use same clustering for all components */
          if (major>=0)
            {
              if (s%sas!=major)
                {
                  AMG_Print("clustered a nonmajor nonisolated component !?\n");
                  sprintf(buffer,"clusternumber=%d, s=%d, major=%d, sas=%d\n",
                          clusternumber,s,major,sas);
                  for (i=0; i<n; i++)
                    {
                      if (ca[i]>=0) continue;
                      sprintf(buffer,"%d: %d\n",i,ca[i]);
                      AMG_Print(buffer);
                    }
                  AMG_Print(buffer);
                  return(AMG_FATAL);
                }
              
              for (comp=0; comp<sas; comp++)
                {
                  if (comp==major) continue;
                  
                  /* make a new cluster */
                  for (i=0; i<clustersize; i++) {
                    f=(cluster[i]/sas)*sas+comp;
                    ca[f]=clusternumber;
                  }
#ifdef DEBUG_CLUSTERING
                  sprintf(buf,"Major   %4d: ",clusternumber); AMG_Print(buf);
                  for (i=0; i<clustersize; i++)
                    {
                      sprintf(buf,"%4d ",(cluster[i]/sas)*sas+comp); AMG_Print(buf);
                    }
                  AMG_Print("\n");
#endif
                  clusternumber++;
                  if (ISOLATED(na[s])) isoclusters++; else conclusters++; /* this is some approx... */
                }
            }
        }
    }
  
  if (cc->verbose>1)
        {
          sprintf(buffer,"%12s: %d clusters (%d:%d) out of %d nodes (%.4lg)\n","Clustering",
                  clusternumber,conclusters,isoclusters,n,((double)n)/((double)clusternumber));
          AMG_Print(buffer);
        }
  
  if (cc->verbose==1)
    {
      sprintf(buffer,"%dC",conclusters);
      AMG_Print(buffer);
    }
  
  /* store cluster count */
  g->clusters=clusternumber;
  g->conclusters=conclusters;
  
  return(AMG_OK);
}

/****************************************************************************/
/*D
   Coarsen - build coarse grid matrix
   
   SYNOPSIS:
   static AMG_MATRIX *Coarsen (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc)
   
   PARAMETERS:
.  A - fine grid matrix
.  g - corresponding graph with clustering information
.  cc - data structure containing parameters.
   
   DESCRIPTION:
   This function builds the coarse grid matrix from the fine grid matrix
   and the clustering information.
   
   RETURN VALUE:
.n NULL error or memory overflow
   
D*/
/****************************************************************************/

/* reconstruct cluster from seed */
static int ReconstructCluster (AMG_GRAPH *g, int f)
{
  int i,k,start,end,first,last,dist;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *na=g->na, *la=g->la;
  int clusternumber;
  
  /* first, reconstruct the cluster */
  clusternumber=ca[f];
  clustersize=1; cluster[0]=f; first=0; last=1; dist=0; /* breadth first search */
  SET_VISITED(na[f]);
  while (1)
    {
      for (i=first; i<last; i++)
        {
          start = ra[cluster[i]]; end = start+ja[start];
          for (k=start+1; k<end; k++){/*printf("i %d  k %d c %d",i,ja[k],ca[ja[k]]);*/
            if ( (ca[ja[k]]==clusternumber) && (!VISITED(na[ja[k]])) )
              {
                SET_VISITED(na[ja[k]]);
                if (clustersize>=AMG_MAX_CLUSTER) return(AMG_FATAL);
                cluster[clustersize++] = ja[k];
              }
            /*printf("\n");*/
          }
        }
      first=last; last=clustersize;
      if (last>first) dist++;
      if (last==first) break;
      if (dist>1000) return(AMG_FATAL);
    }
  
  return(AMG_OK);
}

/* construct coarse grid connectivity for the cluster with seed f */
static int ConstructConnectivity (AMG_GRAPH *g, int f)
{
  int i,k,m,n,start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *na=g->na, *la=g->la;
  int clusternumber;
  
  clusternumber=ca[f];
  
  /* reconstruct cluster */
  if (ReconstructCluster(g,f)!=AMG_OK) return(AMG_FATAL);
  
  /* construct the "front", ie a list of nodes connected to cluster */
  frontsize = 0;
  for (i=0; i<clustersize; i++)
    {
      start = ra[cluster[i]]; end = start+ja[start];
      for (k=start+1; k<end; k++) 
        if ((ca[ja[k]]!=clusternumber) && (!FRONT(na[ja[k]]))) /* not in our cluster and not in front */ 
          {
            if (frontsize>=AMG_MAX_FRONT) {printf("AMG_MAX_FRONT %d \n",frontsize);return(AMG_FATAL);}
            front[frontsize++] = ja[k];
            SET_FRONT(na[ja[k]]);  
          }
    }
  
  /* construct connectivity list from front nodes */
  connectsize = 1; connect[0]=clusternumber;
  for (i=0; i<frontsize; i++)
    {
      n=ca[front[i]];
      if (n<0) return(AMG_FATAL);
      
      /* check if n is in connectivity list */
      for (m=0; m<connectsize; m++)
        if (connect[m]==n) break; /* got it */
      if (m==connectsize) 
        {
          if (connectsize>=AMG_MAX_ROW) {printf("AMG_MAX_ROW %d %d\n",connectsize,AMG_MAX_ROW);return(AMG_FATAL);}
          connect[connectsize++]=n;
        }
    }
  
  /* reset front flags */
  for (i=0; i<frontsize; i++) RESET_FRONT(na[front[i]]);
  frontsize=0;
  
  return(AMG_OK);
}


/****************************************************************************/
/*                                                                            */
/* compute coarse matrix                                                    */
/*                                                                            */
/****************************************************************************/
static AMG_MATRIX *Coarsen (AMG_MATRIX *A, AMG_GRAPH *g, AMG_CoarsenContext *cc)
{
  int i,k,n=g->n;
  int nonzeros,clusters;
  AMG_MATRIX *new;
  int start,end;
  int *ra=g->ra, *ja=g->ja, *ca=g->ca;
  char *na=g->na;
  double *a=AMG_MATRIX_A(A),rescale;
  int bb=AMG_MATRIX_BB(A);
  char buffer[128];
  int min_con=1000000, max_con=-1;
  float *da=g->da;
  float d;
  
  /* pass 1: renumber clusters, we need to process them in ascending order */
  for (i=0; i<n; i++)             /* initialize */ 
    { 
      RESET_FRONT(na[i]); 
      RESET_VISITED(na[i]); 
    }
  clusters=0;
  for (i=0; i<n; i++)             /* this loop relies on the VISITED flag ! */
    {
      if (VISITED(na[i])) continue; /* we have this cluster already */
      if (ReconstructCluster(g,i)!=AMG_OK) 
        {
          sprintf(buffer,"coarsen; could not reconstruct cluster %d, i=%d, cls=%d\n",
                  ca[i],i,clustersize);
          AMG_Print(buffer);
          exit(4711);
        }
      for (k=0; k<clustersize; k++) /* for all nodes in the cluster */
        ca[cluster[k]]=clusters;    /* asign cluster number in array ca */
      
#ifdef DEBUG
      sprintf(buf,"Renumber %4d: ",clusters); AMG_Print(buf);
      for (k=0; k<clustersize; k++)
        {
          sprintf(buf,"%4d ",cluster[k]); AMG_Print(buf);
        }
      AMG_Print("\n");
#endif
      
      clusters++;                   /* increase current cluster number */
    }
  if (clusters!=g->clusters) /* this should never happen */
    {
      printf("coarsen: clusters %d g->clusters %d\n",clusters,g->clusters);
      AMG_Print("number of clusters inconsistent\n");
      exit(4711);      
    }
  
  /* pass 2: count nonzeros to allocate coarse matrix */
  for (i=0; i<n; i++)                /* initialize */ 
    { 
      RESET_FRONT(na[i]); 
      RESET_VISITED(na[i]); 
    }
  nonzeros=0;
  for (i=0; i<n; i++) 
    {
      if (VISITED(na[i])) continue;  /* we have this cluster already */
      if (ConstructConnectivity(g,i)!=AMG_OK) return(AMG_NULL); /* get row */
      nonzeros+=connectsize;         /* connection of this cluster to other clusters */

      min_con=AMG_MIN(min_con,connectsize);
      max_con=AMG_MAX(max_con,connectsize);
      
      d=AutoDamp(g,ca[cluster[0]]);  /* set automatic damping factor */
      for (k=0; k<clustersize; k++) 
        da[cluster[k]]=d;            /* asign damping factor in array da */
      
    }
  
  new = AMG_NewMatrix(clusters,AMG_MATRIX_B(A),nonzeros,AMG_MATRIX_SAS(A),
                      AMG_MATRIX_BLDI(A),1,1,"auto coarsen"); /* allocate the matrix */
  if (new==NULL) 
    {
      AMG_Print("coarse: could not allocate matrix\n");
      exit(4711);
    }
  
  /* pass 3: set up structure of matrix */
  for (i=0; i<n; i++) 
    { 
      RESET_FRONT(na[i]); 
      RESET_VISITED(na[i]); 
    }
  clusters=0;
  for (i=0; i<n; i++)                 /* for all nodes */  
    {
      if (VISITED(na[i])) continue;   /* we have this cluster already */
      if (ConstructConnectivity(g,i)!=AMG_OK) return(AMG_NULL); /* get row */
      if (AMG_SetRowLength(new,clusters,connectsize)!=AMG_OK) /* set length of row */
        {                                                     /* in first entry of */
          AMG_Print("coarse: could not set row length\n");    /* column vector and */
          exit(4711);                                         /* update row and nonzero */
        }                                                     /* vector */
      
      for (k=0; k<connectsize; k++)                           /* for all connections */
        if (AMG_InsertEntry(new,clusters,connect[k])<0)       /* fill column vector of */
          {                                                   /* the coarse matrix */
            AMG_Print("coarse: could not insert entry\n");
            exit(4711);  
          }
      
      clusters++;
    }
  
  /* pass 4: accumulate matrix entries */        
  for (i=0; i<n; i++)                 /* for all nodes */ 
    {
      start = ra[i]; end = start+ja[start]; /* row of the node */
      if (AMG_AddValues(new,ca[i],ca[i],a+(bb*start))<0) /* accumulate diagonal */
        {                                                /* += a[bb*start] */
          AMG_Print("coarse: could not add diagonal value\n");
          exit(4711);
        }
      for (k=start+1; k<end; k++)                        /* accumulate off diagonal */
        if (AMG_AddValues(new,ca[i],ca[ja[k]],a+(bb*k))<0)/* += a[bb*k] */
          {
            AMG_Print("coarse: could not add off diagonal value\n");
            exit(4711);
          }
    }
  
  /* pass %: rescale matrix */
  rescale = 1/1.8;
  for (i=0; i<nonzeros; i++)                 /* for all matrix entries */ 
    new->a[i]*=rescale;

  /* print statistics */
  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: %d rows %d nonzeros %d min_con %d max_con %d avg\n",
              "Coarsen",clusters,nonzeros,min_con,max_con,nonzeros/clusters);
      AMG_Print(buffer);
    }
  
  if (cc->verbose==1)
    {
      AMG_Print("M");
    }
  
  return(new);
}

/****************************************************************************/
/*
   AMG_BuildHierarchy - build hierarchy of coarse grid equations
   
   SYNOPSIS:
   int AMG_BuildHierarchy (AMG_CoarsenContext *cc, AMG_MATRIX *A, 
      AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS]);
   
   PARAMETERS:
.  cc - data structure containing parameters.
.  A - fine grid matrix
.  H -  array of matrix pointers to be filled
.  G - array of corresponding graph structures to be filled
   
   DESCRIPTION:
   This is the function to be called from outside to create a coarse
   grid hierarchy.
   This function creates an AMG hierarchy. It starts from the
   currently coarsest level untel a prescribed depth or coarseness is reached.
   
   RETURN VALUE:
.n -1 error
.n >=0 depth of the hierarchy
*/   
/****************************************************************************/

static void LexPrint (AMG_MATRIX *A[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS], int level)
{
        int i,j,k,n,d;
        char ch,b;
        
        n = (int) sqrt( (double) AMG_MATRIX_N(A[0]) );
        b='0';
        
        for (j=n-1; j>=0; j--)
        {
                for (i=0; i<n; i++)
                {
                        k = j*n+i;
                        for (d=0; d<level; d++) k = G[d]->ca[k];
                        ch = b+k%16;
                        AMG_Print(&ch);
                }
                AMG_Print("\n");
        }
}
                
static int Statistic (AMG_MATRIX *A, int depth)
{
        int i,n=A->n;
        int *ra=A->ra, *ja=A->ja;
        char buffer[196];
        int min_con=1000000, max_con=-1;

        /* compute min / max connections */
        for (i=0; i<n; i++) /* this loop relies on the VISITED flag ! */
        {
                min_con=AMG_MIN(min_con,ja[ra[i]]);
                max_con=AMG_MAX(max_con,ja[ra[i]]);
        }
        sprintf(buffer,"Level %2d: %7d=rows %7d=nonzeros %12.4lg=avg %4d=min %4d=max\n",
                depth,n,A->connections,(double)(((double)A->connections)/((double)n)),min_con,max_con);
        AMG_Print(buffer);
        
        return(AMG_OK);
}


int AMG_BuildHierarchy (AMG_CoarsenContext *cc, AMG_MATRIX *A, 
                        AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS])
{
  AMG_GRAPH *g;
  AMG_MATRIX *C;
  char buffer[128];
  int i,depth;
  clock_t start_coarsen,finish_coarsen;  

  depth=0;                                     /* make level 0 */
  H[0] = A;
                                               /* maximal numbers of levels to build */
  cc->depthtarget = AMG_MIN(cc->depthtarget,AMG_MAX_LEVELS-1); 
  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: start coarsening process\n","Hierarchy");
      AMG_Print(buffer);
    }

  while (1)                                        /* loop */
    {
      start_coarsen = clock();
                                              /* check targets */
      if (depth>=cc->depthtarget) break;        /* finest level reached */
      if (AMG_MATRIX_N(H[depth])<=cc->coarsentarget) break; /* level coarse enough */      
      if (cc->verbose>=1)
        {
          sprintf(buffer,"[%d:",depth);
          AMG_Print(buffer);
        }
                
                                                /* build clusters */
      g = NewGraph(H[depth]);                   /* allocate new graph */
      if (g==NULL) return(-1);
      G[depth] = g;                             /* graph structure on fine grid */
      if (cc->dependency==AMG_SYM)              /* compute dependencies */
        {
          if (Dependency_sym(H[depth],G[depth],cc)<0) 
            {
              printf("2 \n");
              AMG_Print("Symmetric Dependency failed\n");
              cc->dependency=AMG_UNSYM;        /* try unsymmetric dependency */
              if (Dependency(H[depth],G[depth],cc)<0) 
                {
                  printf("3 \n");
                  AMG_Print("Dependency failed\n");
                  return(-1);
                }
            }
        }
      else
        {
          if (Dependency(H[depth],G[depth],cc)<0) 
            {
              printf("3 \n");
              AMG_Print("Dependency failed\n");
              return(-1);
            }
        }
      
      finish_coarsen = clock();
      if (cc->verbose>1)
        {
          sprintf(buffer," %d dependency building %g sec\n ",depth,
                  (double)(finish_coarsen-start_coarsen)/CLOCKS_PER_SEC);
          AMG_Print(buffer);
        }
      
      start_coarsen = clock();
      if (Clustering(G[depth],cc)!=AMG_OK)    /* cluster nodes in small clusters */
        {
          printf("4 \n");
          AMG_Print("Clustering failed\n");
          return(-1);
        }
       finish_coarsen = clock();
      if (cc->verbose>1)
        {
          sprintf(buffer," %d clustering %g sec\n ",depth,
                  (double)(finish_coarsen-start_coarsen)/CLOCKS_PER_SEC);
          AMG_Print(buffer);
        }
     
      /* check targets */
      if ((((double)AMG_MATRIX_N(H[depth]))/((double)G[depth]->clusters) < cc->coarsenrate )
          || (G[depth]->conclusters==0))
        {                                          /* coarsening too slow */ 
          if (cc->verbose>0)
            AMG_Print("coarsening too slow\n");      /* or no unknowns on coarse grid */
          free(G[depth]->ca);
          free(G[depth]->na);
          free(G[depth]->la);
          free(G[depth]->da);
          free(G[depth]);                     
          break;
        }
      
      /* build coarse matrix */
      start_coarsen = clock();
      C=Coarsen(H[depth],G[depth],cc);
      finish_coarsen = clock();
      if (cc->verbose>1)
        {
          sprintf(buffer," %d coarsening %g sec\n ",depth,
                  (double)(finish_coarsen-start_coarsen)/CLOCKS_PER_SEC);
          AMG_Print(buffer);
        }
      if (C==NULL)
        {
          printf("5 \n");
          AMG_Print("Coarsen failed\n");
          return(-1);
        }
      
      if (cc->verbose==1)
        {
          AMG_Print("]");
        }
      
      /* initialize new level */
      depth++;
      H[depth]=C;
      /*LexPrint(H,G,depth);*/
      /*if (C->n<25) AMG_PrintMatrix(C,"mist");*/
    }            /* end while */
  
  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: %d levels in total\n","Hierarchy",depth+1);
      AMG_Print(buffer);
    }
  
  if (cc->verbose==1)
    {
      AMG_Print("\n");
      for (i=0; i<=depth; i++) Statistic(H[i],i);
    }
  
  return(depth);                                
}

int AMG_BuildHierarchy_Saddle (AMG_CoarsenContext *cc, AMG_MATRIX *A, 
                               AMG_MATRIX *H[AMG_MAX_LEVELS], AMG_GRAPH *G[AMG_MAX_LEVELS], 
                               AMG_MATRIX *schur, AMG_MATRIX *H_schur[AMG_MAX_LEVELS], 
                               AMG_GRAPH *G_schur[AMG_MAX_LEVELS],AMG_MATRIX *B[AMG_MAX_LEVELS])

{
  AMG_GRAPH *g,*g_schur;
  AMG_MATRIX *C,*C_schur;
  char buffer[128];
  int i,depth;

  depth=0;                                     /* make level 0 */
  H[0] = A;
  H_schur[0] = schur;
                                               /* maximal numbers of levels to build */
  cc->depthtarget = AMG_MIN(cc->depthtarget,AMG_MAX_LEVELS-1); 
  if (cc->verbose>1)
    {
      sprintf(buffer,"Entering coarsening for saddle point problems.\n");
      AMG_Print(buffer);
      sprintf(buffer,"%12s: start coarsening process\n","Hierarchy");
      AMG_Print(buffer);
    }

  while (1)                                                    /* loop */
    {
                                                            /* check targets */
      if (depth>=cc->depthtarget) break;                    /* finest level reached */
      if (AMG_MATRIX_N(H[depth])<=cc->coarsentarget) break; /* level coarse enough */      
      if (cc->verbose>=1)
        {
          sprintf(buffer,"[%d:",depth);
          AMG_Print(buffer);
        }
                
                                                            /* build clusters for matrix A */
      g = NewGraph(H[depth]);                               /* allocate new graph */
      if (g==NULL) return(-1);
      G[depth] = g;                                         /* graph structure on fine grid */
      if (cc->dependency==AMG_SYM)                          /* compute dependencies */
        {
          if (Dependency_sym(H[depth],G[depth],cc)<0) 
            {
              AMG_Print("Symmetric Dependency failed\n");
              cc->dependency=AMG_UNSYM;                    /* try unsymmetric dependency */
              if (Dependency(H[depth],G[depth],cc)<0) 
                {
                  printf("3 \n");
                  AMG_Print("Dependency failed\n");
                  return(-1);
                }
            }
        }
      else
        {
          if (Dependency(H[depth],G[depth],cc)<0) 
            {
              printf("3 \n");
              AMG_Print("Dependency failed\n");
              return(-1);
            }
        }
                                                           /* dependency for matrix A o.k. */                 
      if (Clustering(G[depth],cc)!=AMG_OK)                 /* cluster nodes in small clusters */
        {
          printf("4 \n");
          AMG_Print("Clustering failed\n");
          return(-1);
        }
      
                                                           /* check targets */
      if ((((double)AMG_MATRIX_N(H[depth]))/((double)G[depth]->clusters) < cc->coarsenrate)
          || (G[depth]->conclusters==0))
        {                                                   /* coarsening too slow */ 
          AMG_Print("coarsening too slow\n");               /* or no unknowns on coarse grid */
          free(G[depth]->ca);
          free(G[depth]->na);
          free(G[depth]->la);
          free(G[depth]->da);
          free(G[depth]);                     
          break;
        }
                                                            /* build coarse matrix */
      C=Coarsen(H[depth],G[depth],cc);
      if (C==NULL)
        {
          printf("5 \n");
          AMG_Print("Coarsen failed\n");
          return(-1);
        }
      
      if (cc->verbose==1)
        {
          AMG_Print("]");
        }

                                                            /* build clusters for matrix schur */
      g_schur = NewGraph(H_schur[depth]);                   /* allocate new graph */
      if (g==NULL) return(-1);
      G_schur[depth] = g_schur;                             /* graph structure on fine grid */
      if (cc->dependency==AMG_SYM)                          /* compute dependencies */
        {
          if (Dependency_sym(H_schur[depth],G_schur[depth],cc)<0) 
            {
              AMG_Print("Symmetric Dependency failed\n");
              cc->dependency=AMG_UNSYM;                    /* try unsymmetric dependency */
              if (Dependency(H_schur[depth],G_schur[depth],cc)<0) 
                {
                  printf("3 \n");
                  AMG_Print("Dependency failed\n");
                  return(-1);
                }
            }
        }
      else
        {
          if (Dependency(H_schur[depth],G_schur[depth],cc)<0) 
            {
              printf("3 \n");
              AMG_Print("Dependency failed for schur matrix\n");
              return(-1);
            }
        }
                                                           /* dependency for matrix A o.k. */                 
      if (Clustering(G_schur[depth],cc)!=AMG_OK)           /* cluster nodes in small clusters */
        {
          printf("4 \n");
          AMG_Print("Clustering failed for schur matrix\n");
          return(-1);
        }
      
                                                           /* check targets */
      if (G_schur[depth]->conclusters==0)                  /* no unknowns on coarse grid */
        {
          printf("4 \n");
          AMG_Print("No unknowns for schur matrix\n");
          return(-1);
        }
                                                           /* build coarse matrix */
      C_schur=Coarsen(H_schur[depth],G_schur[depth],cc);
      /*C_schur=Coarsen(B[depth],G_schur[depth],cc);*/
      if (C_schur==NULL)
        {
          printf("5 \n");
          AMG_Print("Coarsen failed for schur matrix\n");
          return(-1);
        }
      
      if (cc->verbose==1)
        {
          AMG_Print("]");
        }

      
                                                           /* initialize new level */
      depth++;
      H[depth]=C;
      H_schur[depth]=C_schur;
      /*LexPrint(H,G,depth);*/
    }                                                      /* end while */
  
  if (cc->verbose>1)
    {
      sprintf(buffer,"%12s: %d levels in total\n","Hierarchy",depth+1);
      AMG_Print(buffer);
    }
  
  if (cc->verbose==1)
    {
      AMG_Print("\n");
      for (i=0; i<=depth; i++) Statistic(H[i],i);
    }
  if (cc->verbose==1)
    {
      AMG_Print("schur complement matrix\n");
      for (i=0; i<=depth; i++) 
          Statistic(H_schur[i],i);
    }
  
  return(depth);                                
}
