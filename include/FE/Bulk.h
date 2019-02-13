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
// Buld.h
//
// Purpose:   common routines for bulk precipitation in 2d/3d and 3d/4d
//
// Author:    Volker John  
//
// =======================================================================

#ifndef __BULK__
#define __BULK__

// ======================================================================
// lump matrix to diagonal matrix
// the sparsity pattern of the matrix is not condensed
// ======================================================================

#ifdef __2D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix2D *M);
#endif    
#ifdef __3D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix3D *M);
#endif   

double calculate_dp_50(int N, double *size, double *number);

#endif
