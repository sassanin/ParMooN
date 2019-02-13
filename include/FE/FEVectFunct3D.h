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
// @(#)FEVectFunct3D.h        1.2 07/20/99
// 
// Class:       TFEVectFunct3D
// Purpose:     a function from a finite element space in 3D
//
// Author:      Gunar Matthies (13.07.2000)
//
// History:     start of implementation 13.07.2000 (Gunar Matthies)
//
//              WriteSol/ReadSol    13.12.10 (Sashikumaar Ganesan)
// =======================================================================

#ifndef __FEVECTFUNCT3D__
#define __FEVECTFUNCT3D__

#include <FEFunction3D.h>

/** a function from a finite element space */
class TFEVectFunct3D : public TFEFunction3D
{
  protected:
    /** number of components */
    int N_Components;

  public:
    /** constructor with vector initialization */
    TFEVectFunct3D(TFESpace3D *fespace3D, char *name, char *description,
                  double *values, int length, int n_components);

    /** return number of components */
    int GetN_Components()
    { return N_Components; }

    /** return i-th component as FEFunction3D */
    TFEFunction3D *GetComponent(int i)
    {
      return new TFEFunction3D(FESpace3D, Name, Description,
                               Values+i*Length, Length);
    }

    /** convert current grid to vector-values FE function */
    void GridToData();

    /** use current data for grid replacement */
    void DataToGrid();

    /** calculate errors to given vector function */
    void GetDeformationTensorErrors( 
        DoubleFunct3D *Exact, DoubleFunct3D *Exact1,
        DoubleFunct3D *Exact2,
        int N_Derivatives,
        MultiIndex3D *NeededDerivatives,
        int N_Errors, ErrorMethod3D *ErrorMeth, 
        CoeffFct3D *Coeff, TAuxParam3D *Aux,
        int n_fespaces, TFESpace3D **fespaces,
        double *errors);
        
    /** write the solution into a data file **/
    void WriteSol(double t);

    /** Read the solution from a given data file **/
    void ReadSol(char *BaseName);
   
};

#endif
