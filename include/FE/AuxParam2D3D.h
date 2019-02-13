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
// 
// Class:       TAuxParam2D3D
// Purpose:     store parameter functions and FE functions
//
// Author:      Andreas Hahn (19.08.2010)
//
// History:     start of implementation 19.08.2010 (Andreas Hahn)
//
// ======================================================================= 

#ifndef __AUXPARAM2D3D__
#define __AUXPARAM2D3D__

#include <FEVectFunct3D.h>
#include <Enumerations.h>
#include <BaseCell.h>

class TAuxParam2D3D
{
  protected:
    /** stored FEFunctions3D */
    int mN_FEFunctions;
    TFEFunction3D **mFEFunctions;
    
    /** stored FEVectFuncts3D */
    int mN_FEVectFuncts;
    TFEVectFunct3D **mFEVectFuncts;
    
    /** total number of parameters */
    int mN_Parameters;
    int mN_TotalParameters;
    
    int *mFct_Index;
    
    MultiIndex3D *mDerivatives;
    
    int *mIsGlobal;
    
    /** internal storage */
    double n1[MaxN_QuadPoints_3D];
    double n2[MaxN_QuadPoints_3D];
    double n3[MaxN_QuadPoints_3D];
    
    TBaseCell *mCell;
    int mCellNr;
    int mGlobCellNr;
    int mN_Points;
    bool mProject;
    
  protected:
    void GetFunctionsParams(int Fct, MultiIndex3D Der, int index, double **parameters);
    void GetVectFunctParams(int Fct, int index, double **parameters);
    
  public:
    TAuxParam2D3D (int n_fefunctions, TFEFunction3D **fefunctions,
		   int n_fevectfuncts, TFEVectFunct3D **fevectfuncts,
		   int n_parameters, int *fct_index, MultiIndex3D *derivatives,
		   int *isglobal);
		   
    ~TAuxParam2D3D ();
    
    void GetParameters(TBaseCell *cell, int cellnr, int jointnr, int globnr,
		      RefTrans3D RefTrans,
		      int n_points, double *xi, double *eta,
		      double **parameters);
		      
    int GetN_Parameters()
    { return mN_TotalParameters; }
};

#endif
