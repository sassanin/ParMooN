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
   
#ifdef __2D__
void ApproximateRFBSolutionQuadNSE2D(TCollection *Coll, TFEFunction2D *u1, 
				     TFEFunction2D *u2, CoeffFct2D *coeff,
				     double *rhs);

void ApproximateRFBSolutionQuad_Q2_NSE2D(TCollection *Coll, TFEFunction2D *u1,
					 TFEFunction2D *u2, CoeffFct2D *Coeffs,
					 double *rhs);
#endif

#ifdef __3D__
void Compute_Q2_Value(double *coeff, double x, double y, double z, double *val);
void Compute_Q2_Value_Gradient(double *coeff, double x, double y, double z, double *val);

void ApproximateRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, CoeffFct3D *Coeffs,
				     double *rhs);



void ApproximateTimeRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *rhs);

void ApproximateTimeRFBSolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs, double *rhs);

void ApproximateTimeRFB_coupled_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
     TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
     double *rhs);

void ApproximateTimeRFB_coupled_cn_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
     TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
     double *old_small_scales, double *rhs);

void ApproximateTimeRFBSolutionQuad_cn_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
					     double *old_small_scales, double *rhs);
#endif
