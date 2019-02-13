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
* @class     TAssembleMat2D
* @brief     base class for assembling matrices 
* @author    Sashikumaar Ganesan 
* @date      16.04.16
* @History    
 ************************************************************************  */


#ifndef __ASSEMBLEMAT2D__
#define __ASSEMBLEMAT2D__

#include <AllClasses.h>
#include <Assemble.h>
#include <Enumerations.h>

#ifdef _MPI

#endif

#ifdef _OMPONLY

#endif

/** class for 2D scalar assemble */
class TAssembleMat2D: public TAssemble
{
  protected:  

  TFESpace2D **FeSpaces,  **FeRhs;
  TSquareMatrix2D **SqMatrices; 
  TMatrix2D **RecMatrices;
  TDiscreteForm2D *DiscreteForm;
  BoundCondFunct2D **BoundaryConditions;
  BoundValueFunct2D **BoundaryValues;
  TAuxParam2D *AuxParam;                     
  double *Param[MaxN_QuadPoints_2D], *AuxArray[MaxN_QuadPoints_2D];  

 
  public:
    /** Constructor*/
    TAssembleMat2D(int n_fespaces, TFESpace2D **fespaces,
                      int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                      int n_matrices, TMatrix2D **matrices,
                      int n_rhs, double **rhs, TFESpace2D **ferhs,
                      TDiscreteForm2D *discreteform,
                      BoundCondFunct2D **boundarybonditions,
                      BoundValueFunct2D **boundaryvalues,
                      TAuxParam2D *parameters);
    
    /** destrcutor */
    ~TAssembleMat2D();
    
  /** allocate memorey for all aux array */
  void Init();
    
   /** deallocate memorey for all aux array */
  void DeAllocate();

  /** reset all values of Mat and rhs to zero */
  void Reset();
  
   /** assemble the matrices */
  void Assemble2D(); 
  
  /** asssemble slip with friction BC */
  void AssembleNavierSlip();
  
private:  
  /** Add Local Square matrix to Global Matrix */
  void AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);
 
  /** Add Local matrix to Global Matrix */
  void AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct);

  /** Add Local rhs to Global rhs */
  void AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct2D *BaseFuncts, RefTrans2D reftrans);
  
  /**  ModifyMatHang */
  void ModifyMatHang();

  /**  print all matrices */  
  void PrintAllMat();
  
};

#endif
