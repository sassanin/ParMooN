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
   
// ======================================================================
// TCD3D.h
//
// common declaration for time dependent convection diffusion problems
// ====================================================================== 

#include <Constants.h>
#include <Enumerations.h>

namespace TCD3D_Surf
{
  static int N_Terms = 4;
  static MultiIndex3D Derivatives[4] = { D100, D010, D001, D000 };
  static int FESpaceNumbers[4] = { 0, 0, 0, 0 };
  static int N_Matrices = 2;
  static int N_Rhs = 0;
  static int RowSpaces[2] = { 0, 0 };
  static int ColumnSpaces[2] = { 0, 0 };
  static int *RhsSpaces = NULL;
  
  void MatricesAssemble (double Mult, double *coeff, double *param, double hK,
			 double **OrigValues, int *N_BaseFuncts,
			 double ***LocMatrices, double **LocRhs);
			 
  // parameters
  static int N_Parameters = 4;
  static int N_FEFct = 3;
  static int N_FEVectFct = 1;
  static int FctIndex[4] = { 0, 1, 2, 3};
  static MultiIndex3D ParamDerivatives[4] = { D100, D010, D001, D000 };
  static int IsGlobal[4] = { 1, 1, 1, 1 };
}
