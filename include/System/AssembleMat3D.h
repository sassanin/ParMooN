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
* @class     TAssembleMat3D
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      30.05.15
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLEMAT3D__
#define __ASSEMBLEMAT3D__

#include <AllClasses.h>
#include <Assemble.h>
#include <Enumerations.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 3D scalar assemble */
class TAssembleMat3D: public TAssemble
{
  protected:  

  TFESpace3D **FeSpaces,  **FeRhs;
  TSquareMatrix3D **SqMatrices; 
  TMatrix3D **RecMatrices;
  TDiscreteForm3D *DiscreteForm;
  BoundCondFunct3D **BoundaryConditions;
  BoundValueFunct3D **BoundaryValues;
  TAuxParam3D *AuxParam;                     
  double *Param[MaxN_QuadPoints_3D], *AuxArray[MaxN_QuadPoints_3D];  

 
  public:
    /** Constructor*/
    TAssembleMat3D(int n_fespaces, TFESpace3D **fespaces,
                      int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                      int n_matrices, TMatrix3D **matrices,
                      int n_rhs, double **rhs, TFESpace3D **ferhs,
                      TDiscreteForm3D *discreteform,
                      BoundCondFunct3D **boundarybonditions,
                      BoundValueFunct3D **boundaryvalues,
                      TAuxParam3D *parameters);
    
    /** destrcutor */
    ~TAssembleMat3D();
    
  /** allocate memorey for all aux array */
  void Init();
    
   /** deallocate memorey for all aux array */
  void DeAllocate();

  /** reset all values of Mat and rhs to zero */
  void Reset();
  
   /** assemble the matrices */
  void Assemble3D(); 
  
  /** set No penetration BC, that is free slip */
  void AssembleNavierSlip();
  
private:  
  /** Add Local Square matrix to Global Matrix */
  void AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);
 
  /** Add Local matrix to Global Matrix */
  void AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);

  /** Add Local rhs to Global rhs */
  void AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct3D *BaseFuncts, RefTrans3D reftrans);
  
  /**  ModifyMatHang */
  void ModifyMatHang();

  /**  print all matrices */  
  void PrintAllMat();
  
};

#endif
