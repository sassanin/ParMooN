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
 *                         Bulk_3d4d.h                                                  *
 *                        -------------                                                 *
 *                                                                                      *
 ***************************************************************************************/

#ifndef __BULK_3D4D__
#define __BULK_3D4D__

void grid_generator_4d(TCollection *coll, int &N_x, int &N_y, int &N_z,
                       double &x_min, double &x_max, double &y_min, 
                       double &y_max,  double &z_min, double &z_max, 
                       double a_min, double a_max, int N_a,
                       double* &x_coord, double* &y_coord, double* &z_coord, 
                       double* &a_coord, double *a_layers_coord);

void grid_generator_4d_pipe(int N_x, int N_y, int N_z, int N_a,
                            double* &x_coord, double* &y_coord, 
                            double* &z_coord,
                            double* &x_coord_pipe, double* &y_coord_pipe, 
                            double* &z_coord_pipe);

void grid_generator_5d(TCollection *coll, int &N_x, int &N_y, int &N_z,
                       double &x_min, double &x_max, double &y_min, 
                       double &y_max, double &z_min, double &z_max, 
                       double a_min, double a_max, double b_min, double b_max, 
                       int N_a,int N_b, double* &x_coord, double* &y_coord, 
                       double* &z_coord, double* &a_coord, double* &b_coord,
                       double *a_layers_coord,double *b_layers_coord);


int PSD_bound_cound_from_velo_inflow(double x, double y, double z);

void Bulk_FWE_FDM_Upwind_4D(TCollection *coll, TFEFunction3D *velocity1, 
                            TFEFunction3D *velocity2, TFEFunction3D *velocity3,
                            TFEFunction3D *concent_C,
                            double *f_old,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord,
                            double *a_coord,
                            double x_min, double x_max, 
                            double y_min, double y_max,
                            double z_min, double z_max,
                            double a_min, double a_max,
                            double *velo1, double *velo2, double *velo3,
                            double *concent_C_array, int *correspond_3dgrid);


void Bulk_BWE_FDM_Upwind_4D(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2,
                            TFEFunction3D *velocity3, TFEFunction3D *concent_C,
                            double *sol, int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, 
                            double *a_coord, TSquareMatrix3D *mat);

void Compute_Q1_Value_Gradient_4D(double *coeff, double x, double y, double z,
                                  double a, double *val);
void Compute_Q1_Value_4D(double *coeff, double x, double y, double z, double a,
                         double *val);
void Compute_Q1_Gradient_4D(double *coeff, double x, double y, double z,
                            double a, double *val);

void generate_correspond_3d_grid(int N_x, int N_y, int N_z,
                                 double *x_coord, double *y_coord, 
                                 double *z_coord,
                                 TCollection *coll, int *correspond_3dgrid);

void generate_correspond_3d_grid_urea_pipe(int N_x, int N_y, int N_z,
                                           double *x_coord, double *y_coord, 
                                           double *z_coord,
                                           TCollection *coll,
                                           int *correspond_3dgrid);

void CoordsTrafo_SquareToCircle(double x, double y, double& nx, double& ny);


void filling_row_and_col_ptr(int *N_Entries, int Nodes,
                             int N_x, int N_y, int N_z,
                             double x_max, double x_min, 
                             double y_max, double y_min,
                             double z_max, double z_min,
                             double a_max, double a_min,
                             double *x_coord, double *y_coord, double *z_coord,
                             double *a_coord, int *row_ptr, int *col_ptr);


double Bulk_mesh_size_in_convection_direction(double hK, double b1, double b2,
                                              double b3, double b4, double *x, 
                                              double *y, double *z, 
                                              double *a_4d);


void Build_4D_FEM_Matrix_Q1(TCollection *coll,
                            TFEFunction3D *velocity1, TFEFunction3D *velocity2,
                            TFEFunction3D *velocity3, TFEFunction3D *concent_C,
                            double *sol, double *oldsol,
                            double *lump_mass_PSD, double *matrix_D_Entries_PSD,
                            int *correspond_3dgrid,
                            int N_x, int N_y, int N_z, int N_a,
                            double *x_coord, double *y_coord, double *z_coord, 
                            double *a_coord,
                            TSquareMatrix3D *mat, TSquareMatrix3D *matM);

void Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z, int N_a,
                                  double *x_coord, double *y_coord, 
                                  double *z_coord, double *a_coord,
                                  int &N_neum_to_diri, 
                                  int* &neum_to_diri,
                                  double* &neum_to_diri_x,
                                  double* &neum_to_diri_y,
                                  double* &neum_to_diri_z,
                                  double* &neum_to_diri_a);

void Build_4D_FEM_FCT_MassMatrix_Q1(TCollection *coll,
                                    int N_x, int N_y, int N_z, int N_a,
                                    double *x_coord, double *y_coord, 
                                    double *z_coord, double *a_coord,
                                    int* &index_test_ansatz, 
                                    TSquareMatrix3D *matM,
                                    double *lump_mass_PSD);

void Build_4D_FEM_FCT_Matrix_Q1(TCollection *coll,
                                TFEFunction3D *velocity1,
                                TFEFunction3D *velocity2, 
                                TFEFunction3D *velocity3,
                                TFEFunction3D *concent_C,
                                double *sol, double *oldsol,
                                double *lump_mass_PSD,
                                double *matrix_D_Entries_PSD,
                                int *correspond_3dgrid,
                                int N_x, int N_y, int N_z, int N_a,
                                double *x_coord, double *y_coord, 
                                double *z_coord, double *a_coord,
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

void Integral_For_Particle_Increase_Term(TFESpace3D *fespace, 
                                         TFEFunction3D *fefct,
                                         int N_x, int N_y, int N_z, int N_a,
                                         double *x_coord, double *y_coord, 
                                         double *z_coord, double *a_coord,
                                         double *concent_C_array, double *f);

void Integral_For_Particle_Increase_Term_2D(TFESpace3D *fespace, 
                                            TFEFunction3D *fefct,
                                            TFEFunction3D *fefct_2,
                                            int N_x, int N_y, int N_z,
                                            int N_a, int N_b,
                                            double *x_coord, double *y_coord,
                                            double *z_coord,
                                            double *a_coord, double *b_coord,
                                            double *f);           


void Evaluate_f_at_outflow(int N_x, int N_y, int N_z, int N_a, 
                           double *x_coord, double *y_coord, 
                           double *a_coord, double *f,
                           double *average_median, int *average_step);

void write_vtk_file( int N_x, int N_y, int N_z, int N_a, double cut_coord,
                     double *x_coord, double *y_coord, double *z_coord,
                     double *a_coord, double *f_old, const char *name);

void save_f_old_in_txt_file( int Nodes, double *f_old, const char *name);

void read_in_f_old_from_txt_file(int Nodes, double *f_old, const char *name);

void Compute_psd_level(TCollection *coll, TFESpace3D *psd_level_space,
                       TFEFunction3D *psd_level_fe_fct,
                       int N_x, int N_y, int N_z, int N_a, 
                       double *x_coord, double *y_coord, 
                       double *z_coord, double *a_coord, double *f);
#endif
