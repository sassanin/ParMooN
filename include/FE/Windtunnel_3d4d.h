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
 *                         windtunnel_3d4d.h                                            *
 *                        -------------------                                           *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __WINDTUNNEL_3D4D__
#define __WINDTUNNEL_3D4D__
double DROPS_bound_cound_from_velo_inflow(int coord_x, int coord_y, int coord_z, double r,
                                          int N_x, int N_y, int N_z, int N_r);
                                           


void Windtunnel_FWE_FDM_Upwind_4D(TCollection *coll,
			    TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            double *f_old,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            double x_min, double x_max, double y_min, double y_max,
                            double z_min, double z_max, double a_min, double a_max,
                            double *velo1, double *velo2, double *velo3,
			    int *correspond_3dgrid, double ***diff_velo_air_drops);


void Windtunnel_BWE_FDM_Upwind_4D(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            double *sol,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a, double x_min, double x_max,
 			    double y_min, double y_max, double z_min,
                            double z_max, double a_min, double a_max,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord, double ***diff_velo_air_drops, TSquareMatrix3D *mat);

void FEM_FCT_Matrix_Q1_4D_Windtunnel(TCollection *coll,
				TFEFunction3D *velocity1, TFEFunction3D *velocity2, 
				TFEFunction3D *velocity3, double ***diff_velo_air_drops,
				double *sol, double *oldsol,
				double *lump_mass_PSD, double *matrix_D_Entries_PSD,
				int *correspond_3dgrid,
				int N_x, int N_y, int N_z, int N_a,
				double *x_coord, double *y_coord, double *z_coord, double *a_coord,
				TSquareMatrix3D *mat,
				TSquareMatrix3D *matM_cons,
				TSquareMatrix3D *matM,
				int *index_test_ansatz,
				int N_neum_to_diri, 
				int *neum_to_diri,
				double *neum_to_diri_x,
				double *neum_to_diri_y,
				double *neum_to_diri_z,
				double *neum_to_diri_a);


void Compute_Neum_To_Diri_FEM_FCT_Windtunnel(int N_x, int N_y, int N_z, int N_a,
				  double *x_coord, double *y_coord, 
				  double *z_coord, double *a_coord,
				  int &N_neum_to_diri, 
				  int* &neum_to_diri,
				  double* &neum_to_diri_x,
				  double* &neum_to_diri_y,
				  double* &neum_to_diri_z,
				  double* &neum_to_diri_a);

void Build_4D_FEM_FCT_MassMatrix_Q1_Windtunnel(TCollection *coll,
                                    int N_x, int N_y, int N_z, int N_a,
                                    double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                                    int* &index_test_ansatz, 
				    TSquareMatrix3D *matM,
                                    double *lump_mass_PSD);

void Build_4D_FEM_FCT_Matrices_Q1_GroupFEM_Windtunnel(TCollection *coll,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
TSquareMatrix3D *matM, TSquareMatrix3D *matU1,
TSquareMatrix3D *matU2, TSquareMatrix3D *matU3,
TSquareMatrix3D *matG,
double *lump_mass_PSD,double *psd_coeff);

void FEM_FCT_Matrix_Q1_GroupFEM_4D_Windtunnel(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3, double ***diff_velo_air_drops,
double *sol, double *oldsol, double *rhs_psd, double *rhs_psd_old,
double *lump_mass_PSD, double *matrix_D_Entries_PSD,
int *correspond_3dgrid,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord, double *a_layers_coord,
TSquareMatrix3D *mat,
TSquareMatrix3D *matM_cons,
TSquareMatrix3D *matM,
TSquareMatrix3D *matU1,
TSquareMatrix3D *matU2,
TSquareMatrix3D *matU3,
TSquareMatrix3D *matG,
double *psd_coeff,
int N_neum_to_diri,
int *neum_to_diri,
double *neum_to_diri_x,
double *neum_to_diri_y,
double *neum_to_diri_z,
double *neum_to_diri_a);


double normal_rand();

double log_normal_rand(double mu, double sigma, double scal);
double calc_log_normal(double x, double sigma, double mu, double scal);

void write_vtk_file_yzlayer( int N_x, int N_y, int N_z, int N_a, double x_fix_coord,
                     double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                     double *f_old, const char *name);
void write_data_file_meanvalue(double ***mean_value_outflow,  double *a_layers_coord, int N_x, int N_y ,int N_z, int N_a, const char *name);
void compute_coordinate(int index, int *coord_x, int *coord_y, 
                   int *coord_z, int *coord_a, int N_x, int N_y, 
                     int N_z, int N_a);

void midpointflow(double *x_coord, double *y_coord, 
                   double *z_coord, double *a_coord, int N_x, int N_y, 
		  int N_z, int N_a, double *sol );
void compute_mean_value_outflow(int*** mean_value_outflow_indices, double*** mean_value_outflow, double*sol, double* x_coord,
                      int N_x, int N_y, int N_z, int N_a, int *only_first_time); 
double calc_velo_u1(double *** diff_velo_air_drops, int coord_int_x, int coord_int_y, int coord_int_z,
       int  a_coord_int, int N_x, int N_y, int N_z, int N_a);
void alloc_cubix(double ****cubix, int dim_x, int dim_y, int dim_z);
void disalloc_cubix(double ****cubix, int dim_x, int dim_y, int dim_z);
void alloc_cubix_int(int ****cubix, int dim_x, int dim_y, int dim_z);
void disalloc_cubix_int(int ****cubix, int dim_x, int dim_y, int dim_z);
void alloc_matrix(double ***matrix, int dim_x, int dim_y);
void disalloc_matrix(double **matrix, int dim_x, int dim_y);
void alloc_matrix_int(int ***matrix, int dim_x, int dim_y);
void disalloc_matrix_int(int **matrix, int dim_x, int dim_y);

void Q_criterion(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *Qcrit);

void calc_inflow_outflow_mass(double *a_layers_coord, double *f, double *x_coord, int N_x,
int N_y, int N_z, int N_r, double *mass_inflow, double  *mass_outflow );

void calc_inflow_outflow_mass_m3(double *a_layers_coord, double *f, double *x_coord, double *y_coord, int N_x,
int N_y, int N_z, int N_r, double *mass_inflow, double  *mass_outflow );

void write_vtk_q_crit( int N_x, int N_y, int N_z, int N_a, double *x_coord, double *y_coord, double *z_coord, double *a_coord,
double *q_crit, const char *name);

void ComputeIntegralLengthScale(TCollection *coll,
TFEFunction3D *velocity1,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
int *correspond_3dgrid);

#endif
