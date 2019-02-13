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
// @(#)Upwind.h        1.4 10/18/99
//
// Purpose:     upwind stabilization
//              do upwind assembling for first order nonconforming elements
//
// Authors:     Volker Behns / Gunar Matthies  18.10.99
//
// =======================================================================


#ifndef __UPWIND__
#define __UPWIND__

void UpwindForNavierStokes(CoeffFct2D *Coeff, TSquareMatrix2D *sqmatrix, 
			   TFEFunction2D *u1, TFEFunction2D *u2);

void UpwindForConvDiff(CoeffFct2D *Coeff, 
		       TSquareMatrix2D *sqmatrix, double *RHS,
                       TFESpace2D *fespace, TDiscreteForm2D *DiscreteForm,
		       TFEFunction2D *u1, TFEFunction2D *u2,
		       int ConvIsVelo);

/******************************************************************************/
//
// IMPROVED MIZUKAMI-HUGHES METHOD (Knobloch, CMAME 2007)
//
/******************************************************************************/

void ComputeParametersMizukamiHughes(TBaseCell *cell, int cell_no, 
				     TFEFunction2D *u, CoeffFct2D *Coeffs,
				     BoundCondFunct2D *BoundaryCondition,
                                     int *dof, int ActiveBound,
				     double *c_mh);

void MizukamiHughes(TSquareMatrix2D *sqmatrix, 
		    double *RHS,
		    TFESpace2D *fespace, 
		    TFEFunction2D *u, 
		    CoeffFct2D *Coeffs,
		    BoundCondFunct2D *BoundaryCondition);


#endif
