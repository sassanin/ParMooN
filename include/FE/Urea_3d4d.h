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

#ifndef __UREA_3D4D__
#define __UREA_3D4D__
double growth_rate(double c, double temp);
double b_nuc(double c, double temp);
double InletPSD(double a);
double psd_boundary_urea(double x, double y, double z, double a,
                    double c, double temp);
int PSD_bound_cond_from_velo_inflow_urea(double x, double y, double z);
int PSD_bound_cond_from_velo_inflow_urea_pipe(double x, double y, double z);
void Urea_FWE_FDM_Upwind_4D(TCollection *coll,
			    TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            TFEFunction3D *Temp,
                            double *f_old,  double *f_new,
                            double *rhs_psd,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            double x_min, double x_max, double y_min, double y_max,
                            double z_min, double z_max, double a_min, double a_max,
                            double *velo1, double *velo2, double *velo3,
                            double *concent_C_array, double *Temp_array,
			    int *correspond_3dgrid);

void Urea_RKV_FDM_4D(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *f_old, double *rhs_psd,  double **stages, 
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
double *velo1, double *velo2, double *velo3, double *concent_C_array, double *Temp_array,
         int *correspond_3dgrid);


void Urea_BWE_FDM_Upwind_4D(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            TFEFunction3D *Temp,
                            double *sol,
                            double *rhs_psd,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
                            TSquareMatrix3D *mat);


void Urea_Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord,
double *z_coord, double *a_coord,
int &N_neum_to_diri,
int* &neum_to_diri,
double* &neum_to_diri_x,
double* &neum_to_diri_y,
double* &neum_to_diri_z,
				       double* &neum_to_diri_a);

void Urea_Build_4D_FEM_FCT_Matrix_Q1(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *sol, double *oldsol,double *rhs_psd, double *rhs_psd_old,
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

void FEM_FCT_Matrix_Q1_GroupFEM_4D_Urea(TCollection *coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *concent_C,
TFEFunction3D *Temp,
double *sol, double *oldsol, double *rhs_psd, double *rhs_psd_old,
double *lump_mass_PSD, double *matrix_D_Entries_PSD,
int *correspond_3dgrid,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord, double *a_coord,
double *x_coord_pipe, double *y_coord_pipe, double *z_coord_pipe,
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

void Build_4D_FEM_FCT_Matrices_Q1_GroupFEM_Urea(TCollection *coll,
            int N_x, int N_y, int N_z, int N_a,
            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
            TSquareMatrix3D *matM, TSquareMatrix3D *matU1, 
            TSquareMatrix3D *matU2, TSquareMatrix3D *matU3, 
            TSquareMatrix3D *matG,
            double *lump_mass_PSD);

void Build_4D_FEM_FCT_Matrices_Q1_GroupFEM_Urea_Pipe(TCollection *coll, int *correspond_3dgrid,
            int N_x, int N_y, int N_z, int N_a,
            double *x_coord, double *y_coord, double *z_coord, double *a_coord,
            double *x_coord_pipe, double *y_coord_pipe, double *z_coord_pipe, 
            TSquareMatrix3D *matM, TSquareMatrix3D *matU1, 
            TSquareMatrix3D *matU2, TSquareMatrix3D *matU3, 
            TSquareMatrix3D *matG,
            double *lump_mass_PSD);


double calculate_q_3(int N, double a_layers_coord, double sol_psd);

void Evaluate_f_at_outflow1(int N_x, int N_y, int N_z, int N_a,
                           double *x_coord, double *y_coord,double *z_coord,
                           double *a_coord, double *f);


void Calculate_Volume_Distribution(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, 
			   double *x_coord, double *y_coord, double *z_coord,
				   double *a_layers_coord, double *f);

void Calculate_PSD_on_node(int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
			   double *a_layers_coord, double *sol_psd, double x, double y, double z);

void Calculate_PSD_outflow(TCollection *Coll,int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
         double *a_layers_coord, double *sol_psd,int *step, double x_end);

void PrepareAgglomerationBreakage(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3,
TFEFunction3D *temperature,
int N_x, int N_y, int N_z, int N_a,
double *x_coord, double *y_coord, double *z_coord,
				  double *a_layers_coord, double *f, double *rhs_new);



void call_apply_integral_operator(int nx, int ny, int nz, int na, double* input, double* v, 
                                  double* grad_v, double* temp, double* output, double* grid,
                                  double* params);
  
double kernel_urea_old(double x, double x_bar);
double gauss_source_urea_old(int j, int index, double *a_layers_coordinate, 
        double a_layers_coordinate_outerLoop, double *f, int N3, int Na);
double gauss_sink_urea_old(int j, int index, double *a_layers_coordinate, 
      double a_layers_coordinate_outerLoop, double *f, int N3); 
void aggregation_conv_urea_old(TCollection *Coll,int N_x, int N_y, int N_z, int N_a, double *x_coord, double *y_coord, double *z_coord,
          double *f, double *rhs_new, double *a_layers_coord, TFEFunction3D *temperature);      
double kernel_urea(double x, double x_bar, double shear_rate);
double psd_value(double *f, double *a_layers_coord, double val, int index, int N_a, int N3);
double gauss_source_urea(int j, int index, double *a_layers_coordinate, double a_layers_coordinate_outerLoop, double *f, int N3, int Na, double shear_rate);
double gauss_sink_urea(int j, int index, double *a_layers_coordinate, double a_layers_coordinate_outerLoop, double *f, int N3, double shear_rate);
void aggregation_conv_urea(TCollection *Coll,double *x_coord, double *y_coord, double *z_coord, 
         int N_x, int N_y, int N_z, int N_a, 
               double *f, double *rhs_new, double *a_layers_coord, TFEFunction3D *temperature);

#endif
