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
// VMS.h
//
// Purpose:     routines for projection-based VMS
//
// Author:       Volker John  2006/05/18
//
// =======================================================================

#ifndef __VMS__
#define __VMS__

#include <SquareMatrix2D.h>
#include <Matrix2D.h>

#ifdef __3D__
  #include <SquareMatrix3D.h>
  #include <Matrix3D.h>
#endif

#ifdef __2D__
void VMSProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
				 TSquareMatrix2D **SQMATRICES, 
				 TMatrix2D **MATRICES);

void LumpMassMatrixToDiag(TSquareMatrix2D *M);

#endif // __2D__

#ifdef __3D__
void VMS_ProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
                             TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES);

void VMS_ProjectionExplUpdateRhs(int N_U, int N_Active, int N_L, TFEVectFunct3D *u,
				 TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES, 
				 double *rhs_vms_expl);

void LumpMassMatrixToDiag(TSquareMatrix3D *M);

void ComputeVMSProjection(TMatrix3D *G11, TMatrix3D *G22, TMatrix3D *G33,
			  TSquareMatrix3D *MatrixL, TFEFunction3D *u_1, 
			  TFEFunction3D *u_2, TFEFunction3D *u_3,
			  TFEVectFunct3D *vms_projection_fe);

void ComputeSizeOfSmallScales(TMatrix3D *matG11, TMatrix3D *matG22, TMatrix3D *matG33,
			  TSquareMatrix3D *MatrixL, TFEFunction3D *u1, 
			  TFEFunction3D *u2, TFEFunction3D *u3,
			  TFEVectFunct3D *vms_projection_fe, double *size_small_scales);

void MeanAndLargestSize( TFESpace3D *projection_space, double *size_small_scales, 
                         double *mean, double *largest_size);

void AdaptProjectionSpace(TFESpace3D *projection_space, 
			  double *size_small_scales, 
			  FE3D  *fes, 
			  double mean, 
			  double mean_time_average, 
			  double largest_size, 
			  double max_time_average, 
			  double *label_space);


#endif // __3D__ 

#endif
