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
   
//***************************************************** 
//  Superconvergence header file for Laplace, 
//     Stokes and Navier-Strokes 2D problems
//*****************************************************
//
#include <FEFunction2D.h>

void Transform_Q2Q3_2D(double *fine_values, double *coarse_values);

void Transform_Q2Q4_2D(double *fine_values, double *coarse_values);

void Transform_P1P2_1_2D(double *fine_values, double *coarse_values); 

void Transform_P1P2_2_2D(double *fine_values, double *coarse_values); 

void Superconvergence_Q1Q2_2D(TFEFunction2D *q1_function, 
			    TFEFunction2D *q2_function);

void Superconvergence_Q2Q3_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q3_function);
 
void Superconvergence_Q2Q4_2D(TFEFunction2D *q2_function, 
			    TFEFunction2D *q4_function); 

void Superconvergence_P1P2_2D(int version, TFEFunction2D *p1_function, 
			    TFEFunction2D *p2_function); 

void Superconvergence_NQ1P2_2D(TFEFunction2D *qn1_function, 
			    TFEFunction2D *p2_function); 
