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
   
/** ************************************************************************ 
*
* @class     TAssemble
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      30.05.15
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLE__
#define __ASSEMBLE__

#include <AllClasses.h>
#include <Constants.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 3D scalar system matrix */
class TAssemble
{
  protected:  
    
    TCollection *Coll;
      
    int N_Cells;
    
    /** Number of Matrices */
    int N_SqMatrices, N_Matrices, N_AllMatrices;
    
    /** Number of FE spaces */
    int N_FeSpaces;    
    
    /** Number of RHS */    
    int N_Rhs;
    
    /** Variable for methods */
    int N_Parameters, **GlobalNumbers, **BeginIndex, **RhsGlobalNumbers, **RhsBeginIndex;
    int **TestGlobalNumbers, **TestBeginIndex, **AnsatzGlobalNumbers, **AnsatzBeginIndex;
    double **HangingEntries, **HangingRhs;
    double **Rhs, **Matrix, ***LocMatrices, **LocRhs, *aux, *auxarray, *rhsaux, *paramaux;
    double *auxmat, **Matrices;
    
    bool *SecondDer;
    
  public:
    /** Constructor*/
    TAssemble(int n_fespaces, int n_sqmatrices, int n_matrices, int n_rhs, double **rhs);
    
    /** destrcutor */
    ~TAssemble();
    
};

#endif
