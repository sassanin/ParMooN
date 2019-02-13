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
   
void ComputeCoordinates(int i, int N, int Nz, double *x, double *y, double *z);

void Build_3D_FDM_Matrix(double h, TFEFunction2D *velocity1, TFEFunction2D *velocity2, 
                         TFEFunction2D *concent_C, double *f_old, 
			 double *velo1, double *velo2);

void Evalute_f_at_outflow(int n_dof_FDM, int Nx, int Nz, double *f);

void Integral_For_Particle_Increase_Term(TFESpace2D *fespace, TFEFunction2D *fefct,
					 int n_dof_FDM, int Nx, int Nz, double *f);
