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
   
/****************************************************************************************
 *                                                                                      *
 *                         Urea_3d4d.h                                                  *
 *                        -------------                                                 *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __UREA_2D3D__
#define __UREA_2D3D__
void grid_generator_3d_urea(TCollection *coll,int &N_x,
int &N_y, double &x_min, double &x_max, double &y_min, double &y_max,
double a_min, double a_max, int N_a,
double* &x_coord, double* &y_coord, double* &a_coord, double* a_layers_coord);
double growth_rate(double c, double temp);
double b_nuc(double c, double temp);

int PSD_bound_cond_from_velo_inflow_urea_2D(double x, double y);
void Urea_FWE_FDM_Upwind_3D(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2,
TFEFunction2D *concent_C, TFEFunction2D *Temp,
double *f_old,  double *rhs_psd,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
double x_min, double x_max, double y_min, double y_max,
double z_min, double z_max,
double *velo1, double *velo2, double *concent_C_array,
			    int *correspond_2dgrid
);

 
void Urea_BWE_FDM_Upwind_3D(TCollection *coll,
			    TFEFunction2D *velocity1, TFEFunction2D *velocity2,
			    TFEFunction2D *concent_C, TFEFunction2D *Temp,
			    double *sol, double *rhs_psd,
			    int *correspond_2dgrid,
			    int N_x, int N_y, int N_z,
			    double *x_coord, double *y_coord, double *z_coord,
			    TSquareMatrix2D *mat);



void Urea_RKV_FDM_3D(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2, TFEFunction2D *concent_C,
TFEFunction2D *Temp,
double *f_old, double *rhs_psd,  double **stages,
int N_x, int N_y, int N_z, double *x_coord, double *y_coord, double *z_coord, 
double *velo1, double *velo2, double *concent_C_array, int *correspond_2dgrid);


void FEM_FCT_Matrix_Q1_GroupFEM_3D_Urea(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2,
TFEFunction2D *concent_C, TFEFunction2D *Temp,
double *sol, double *oldsol, double *rhs_psd, double *rhs_psd_old,
double *lump_mass_PSD, double *matrix_D_Entries_PSD,
int *correspond_2dgrid,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
TSquareMatrix2D *mat,
TSquareMatrix2D *matM_cons,
TSquareMatrix2D *matM,
TSquareMatrix2D *matU1,
TSquareMatrix2D *matU2,
TSquareMatrix2D *matG,
double *psd_coeff,
int N_neum_to_diri,
int *neum_to_diri,
double *neum_to_diri_x,
double *neum_to_diri_y,
double *neum_to_diri_z);

void Compute_Neum_To_Diri_FEM_FCT_urea(int N_x, int N_y, int N_z,
					 double *x_coord, double *y_coord, 
					 double *z_coord,
					 int &N_neum_to_diri, 
					 int* &neum_to_diri,
					 double* &neum_to_diri_x,
					 double* &neum_to_diri_y,
					 double* &neum_to_diri_z);


void Build_3D_FEM_FCT_Matrices_Q1_GroupFEM_Urea(TCollection *coll,
				    int N_x, int N_y, int N_z,
				    double *x_coord, double *y_coord, double *z_coord,
				    TSquareMatrix2D *matM,TSquareMatrix2D *matU1, 
						      TSquareMatrix2D *matU2, TSquareMatrix2D *matG,
				    double *lump_mass_PSD);

double psd_boundary_urea_2D(double x, double y, double a,
double c, double temp);

double InletPSD(double a);
void PrepareAgglomerationUrea2D(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2, TFEFunction2D *temperature, int N_x,
int N_y, int N_a,
double *x_coord, double *y_coord, double *a_layers_coord, double *f, double *rhs_new, int 
*correspond_2dgrid);

void Calculate_PSD_on_node(int N_x, int N_y, int N_a, 
double *x_coord, double *y_coord, double *a_layers_coord, double *sol_psd, double x, double y);

void Calculate_PSD_outflow_2D(TCollection *Coll,int N_x, int N_y,  int N_a,
double *x_coord, double *y_coord,
double *a_layers_coord, double *sol_psd, int step,double x_end);
#endif
