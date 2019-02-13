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
   
// =======================================================================
// @(#)FESpace.h        1.1 10/30/98
// 
// Class:       TFESpace
// Purpose:     general super class for all finite element spaces
//              special spaces are implemented in subclasses
//
// Author:      Gunar Matthies (03.11.97)
//
// History:     start of implementation 03.11.97 (Gunar Matthies)
//
//              add data for storing information on the used elements
//              29.11.1997 (Gunar Matthies)
//
//              split FESpace into TFESpacexD
//              15.04.1998 (Volker Behns)
//
//              start of reimplementation
//              30.07.1998 (Gunar Matthies)
//
// =======================================================================

#ifndef __FESPACE__
#define __FESPACE__

#include <Collection.h>

/** general super class for all finite element spaces, special spaces are
    implemented in subclasses */
class TFESpace
{
  protected:
// =======================================================================
// administration
// =======================================================================
    /** name of the space */
    char *Name;

    /** some more words describing the space */
    char *Description;

// =======================================================================
// information of cell collection
// =======================================================================
    /** Collection containing the cells used for building this space */
    TCollection *Collection;

    /** number of cells in the triangulation used for building this space */
    int N_Cells;

// =======================================================================
// general information on degrees of freedom and their numbers
// =======================================================================
    /** number of all degrees of freedom */
    int N_DegreesOfFreedom;

    /** array containing all global numbers of local degrees of freedom
        for all elements */
    int *GlobalNumbers;

    /** array containing the begin index in GlobalNumbers array for each
        element */
    int *BeginIndex;

    /** number of used elements */
    int N_UsedElements;

// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** number of different boundary node type, except Dirichlet */
    int N_DiffBoundNodeTypes;

    /** type for each of the different types */
    BoundCond *BoundaryNodeTypes;

    /** number of Dirichlet nodes */
    int N_Dirichlet;

    /** number of nodes for each boundary node type */
    int *N_BoundaryNodes;

    /** number of inner nodes */
    int N_Inner;

    /** 0 <= i < InnerBound for all inner degrees of freedom i */
    int InnerBound;

    /** InnerBound <= i < NeumannBound for all Neumann nodes i */
    int *BoundaryNodesBound;

    /** NeumannBound <= i < DirichletBound for all Dirichlet node i */
    int DirichletBound;

    /** number of inner and non-Dirichlet boundary nodes are less than */
    int ActiveBound;

    /** 0 space for Galerkin disc, 1 - space for DG disc */
    int DGSpace;

#ifdef  _MPI
     /** Maximum number of subdomains associated with any dof */
    int MaxSubDomainPerDof;
#endif

  private:
    /** copying given parameters into inner storage places */
    int InitData(TCollection *coll, char *name, char *description);

  public:
    /** constructor */
    TFESpace(TCollection *coll, char *name, char *description);

    /** destrcutor */
    ~TFESpace();

    /** return name */
    char *GetName() const
    { return Name; }

    /** return description */
    char *GetDescription() const
    { return Description; }

    /** return number of cells in the triangulation used for building 
        this space */
    int GetN_Cells() const
    { return N_Cells; }

    /** return the collection of this space */
    TCollection *GetCollection() const
    { return Collection; }

    /** return global numbers of local degrees of freedom */
    int *GetGlobalNumbers() const
    { return GlobalNumbers; }
    
    void SetGlobalNumbers(int* NewGN)
    { GlobalNumbers=NewGN; }

    /** return begin index for each element */
    int *GetBeginIndex() const
    { return BeginIndex; }
    
    /** @brief return correspondence map from local to global degrees of freedom
     * 
     * set int * DOF=feSpace->GetGlobalDOF(i); then DOF[j]-th global degree of 
     * freedom corresponds to the j-th local degree of freedom.
    */
    int* GetGlobalDOF(int i) const
    { return GlobalNumbers+BeginIndex[i];}

    /** return number of used elements */
    int GetN_UsedElements() const
    { return N_UsedElements; }

    /** return number of all degrees of freedom */
    int GetN_DegreesOfFreedom() const
    { return N_DegreesOfFreedom; }

// =======================================================================
// counts and bounds for different types of degrees of freedom
// =======================================================================
    /** get number of different boundary node types, except Dirichlet */
    int GetN_DiffBoundaryNodeTypes() const
    { return N_DiffBoundNodeTypes; }

    /** return type for each of the different types */
    BoundCond *GetBoundaryNodeTypes() const
    { return BoundaryNodeTypes; }

    /** return number of nodes for each boundary node type */
    int *GetN_BoundaryNodes() const
    { return N_BoundaryNodes; }

    /** return N_Dirichlet */
    int GetN_Dirichlet() const
    { return N_Dirichlet; }

    /** return N_Inner */
    int GetN_Inner() const
    { return N_Inner; }

    /** return InnerBound */
    int GetInnerBound() const
    { return InnerBound; }

    /** return BoundaryNodesBound */
    int *GetBoundaryNodesBound() const
    { return BoundaryNodesBound; }

    /** return DirichletBound */
    int GetDirichletBound() const
    { return DirichletBound; }

    /** return ActiveBound */
    int GetActiveBound() const
    { return ActiveBound; }

    /** write info on fespace into file */
    int Write(const char *filename);

    void SetAsDGSpace()
    { DGSpace = 1; }

    int IsDGSpace() const
    { return DGSpace; }

#ifdef  _MPI
    /** return  MaxSubDomainPerDof */
    void SetMaxSubDomainPerDof(int maxSubDomainPerDof)
    { MaxSubDomainPerDof = maxSubDomainPerDof; }

    /** return  MaxSubDomainPerDof */
    int GetMaxSubDomainPerDof()
    { return MaxSubDomainPerDof; }
#endif

};

#endif
