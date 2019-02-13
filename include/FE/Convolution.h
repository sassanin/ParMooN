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
// @(#)Convolution.h        1.2 04/13/00
//
// Purpose:     convolute velocity and tensors
//
// Authors:     Volker John, Gunar Matthies
// =======================================================================

#ifndef __CONVOLUTION__
#define __CONVOLUTION__

#ifdef __2D__
#include <FEVectFunct2D.h>
#endif

#ifdef __3D__
#include <FEFunction3D.h>
#include <FEVectFunct3D.h>
#endif


double CharacteristicFilterWidth(double h);

#ifdef __2D__
double GaussianFilter(double delta, double dist_sq);

void  ConvoluteVelocity(TFEVectFunct2D *u, TFEVectFunct2D *uConv);

void  ConvoluteVelocityFull(TFEVectFunct2D *u, TFEVectFunct2D *uConv);


// ========================================================================
// convolute (grad w grad w)
// ========================================================================
void  ConvoluteDuTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);

void  ConvoluteSymmetricTensor(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);

void  ConvoluteSymmetricTensorFull(TFEVectFunct2D *u, TFEVectFunct2D *duTensor);
#endif

#ifdef __3D__
double GaussianFilter3D(double delta, double dist_sq);

void  ConvoluteVelocity3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv);
void  ConvoluteVelocityFull3D(TFEVectFunct3D *u, TFEVectFunct3D *uConv);
void  ConvoluteSymmetricTensor3D(TFEVectFunct3D *u, TFEVectFunct3D *duTensor);
void  ConvoluteSymmetricTensorFull3D(TFEVectFunct3D *u, 
                                     TFEVectFunct3D *duTensor);
void  ConvolutePressure3D(TFEFunction3D *u, TFEFunction3D *uConv);
#endif

#endif
