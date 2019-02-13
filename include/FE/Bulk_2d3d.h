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
// Buld_2d3d.h
//
// Purpose:   routines for bulk precipitation in 2d/3d
//
// Author:    Volker John  
//
// =======================================================================

#ifndef __BULK_2D3D__
#define __BULK_2D3D__

#ifdef __2D__
void Bulk_FWE_FDM_Upwind_3D(TCollection *coll,
TFEFunction2D *velocity1, TFEFunction2D *velocity2,
TFEFunction2D *concent_C,
double *f_old, double *rhs_psd,
int N_x, int N_y, int N_z,
double *x_coord, double *y_coord, double *z_coord,
double x_min, double x_max, double y_min, double y_max,
double z_min, double z_max,
double *velo1, double *velo2, double *concent_C_array,
			    int *correspond_2dgrid
    );

void Bulk_RKV_FDM_3D(TCollection *coll,
         TFEFunction2D *velocity1, TFEFunction2D *velocity2,
         TFEFunction2D *concent_C,
         double *f_old, double **stages, 
         int N_x, int N_y, int N_z,
         double *x_coord, double *y_coord, double *z_coord,
         double *velo1, double *velo2, double *concent_C_array,
         int *correspond_2dgrid);


void Bulk_BWE_FDM_Upwind_3D(TCollection *coll,
			    TFEFunction2D *velocity1, TFEFunction2D *velocity2,
			    TFEFunction2D *concent_C,
			    double *sol,
			    int *correspond_2dgrid,
			    int N_x, int N_y, int N_z,
			    double *x_coord, double *y_coord, double *z_coord,
			    TSquareMatrix2D *mat);

void Evalute_f_at_outflow(int N_x, int N_y, int N_z, double *x_coord, double *z_layers_coord, 
			  double *f, double *average_median, int *average_step);

void Integral_For_Particle_Increase_Term(TFESpace2D *fespace, TFEFunction2D *fefct,
                                         int N_x, int N_y, int N_z,
                                         double *x_coord, double *y_coord, double *z_layers_coord,
                                         double *concent_C_array, double *f);

double Bulk_mesh_size_in_convection_direction(double hK, double b1, double b2, 
					      double b3, double *x, double *y,
					      double *z);

// void Build_3D_FEM_Matrix_Q1(double h, TFEFunction2D *velocity1, TFEFunction2D *velocity2,
//                             TFEFunction2D *concent_C, double *f_old,
//                             double *velo1, double *velo2, TCollection *coll);

void grid_generator_3d(double x_min, double x_max, int N_x,
                       double y_min, double y_max, int N_y,
                       double z_min, double z_max, int N_z,
                       double *x_coord, double *y_coord, double *z_coord,
                       double *z_layers_coord);

void generate_correspond_2d_grid(int N_x, int N_y, double *x_coord, double *y_coord,
                                 TCollection *coll, int *correspond_2dgrid);

void filling_row_and_col_ptr(int *N_Entries, int Nodes, int N_x, int N_y, double x_max, double x_min,
                             double y_max, double y_min, double z_max, double z_min,
                             double *x_coord, double *y_coord, double *z_coord,
                             int *row_ptr, int *col_ptr);

void Build_3D_FEM_Matrix_Q1(TCollection *coll,
                            TFEFunction2D *velocity1, TFEFunction2D *velocity2,
                            TFEFunction2D *concent_C,
                            double *sol, double *oldsol,
                            double *lump_mass_PSD, double *matrix_D_Entries_PSD,
                            int *correspond_2dgrid,
                            int N_x, int N_y, int N_z,
                            double *x_coord, double *y_coord, double *z_coord,
                            TSquareMatrix2D *mat, TSquareMatrix2D *matM);

void Compute_Neum_To_Diri_FEM_FCT(int N_x, int N_y, int N_z,
					 double *x_coord, double *y_coord, 
					 double *z_coord,
					 int &N_neum_to_diri, 
					 int* &neum_to_diri,
					 double* &neum_to_diri_x,
					 double* &neum_to_diri_y,
					 double* &neum_to_diri_z);

void Build_3D_FEM_FCT_MassMatrix_Q1(TCollection *coll,
				    int N_x, int N_y, int N_z,
				    double *x_coord, double *y_coord, double *z_coord,
				    int* &index_test_ansatz, 
//				    int* &index_test_ansatz_diag, 
				    TSquareMatrix2D *matM,
				    double *lump_mass_PSD);

void Build_3D_FEM_FCT_Matrix_Q1(TCollection *coll,
                            TFEFunction2D *velocity1, TFEFunction2D *velocity2,
                            TFEFunction2D *concent_C,
                            double *sol, double *oldsol,
                            double *lump_mass_PSD, double *matrix_D_Entries_PSD,
                            int *correspond_2dgrid,
                            int N_x, int N_y, int N_z,
                            double *x_coord, double *y_coord, double *z_coord,
				TSquareMatrix2D *mat, TSquareMatrix2D *matM_cons,
				TSquareMatrix2D *matM,
				int *index_test_ansatz,
//				int *index_test_ansatz_diag,
				int N_neum_to_diri, 
				int *neum_to_diri,
				double *neum_to_diri_x,
				double *neum_to_diri_y,
				double *neum_to_diri_z);

void Build_3D_FEM_FCT_Matrices_Q1_GroupFEM_Bulk(TCollection *coll, int N_x, 
                                                int N_y, int N_z, 
                                                double *x_coord, 
                                                double *y_coord, 
                                                double *z_coord,
                                                TSquareMatrix2D *matM,
                                                TSquareMatrix2D *matU1, 
                                                TSquareMatrix2D *matU2,
                                                TSquareMatrix2D *matG,
                                                double *lump_mass_PSD);
            
void FEM_FCT_Matrix_Q1_GroupFEM_3D_Bulk(TCollection *coll, 
                                        TFEFunction2D *velocity1, 
                                        TFEFunction2D *velocity2,
                                        TFEFunction2D *concent_C,
                                        double *sol, double *oldsol, 
                                        double *lump_mass_PSD, 
                                        double *matrix_D_Entries_PSD,
                                        int *correspond_2dgrid,
                                        int N_x, int N_y, int N_z,
                                        double *x_coord, double *y_coord, 
                                        double *z_coord, 
                                        TSquareMatrix2D *mat,
                                        TSquareMatrix2D *matM_cons,
                                        TSquareMatrix2D *matM,
                                        TSquareMatrix2D *matU1,
                                        TSquareMatrix2D *matU2,
                                        TSquareMatrix2D *matG,
                                        double *psd_coeff, int N_neum_to_diri,
                                        int *neum_to_diri, 
                                        double *neum_to_diri_x,
                                        double *neum_to_diri_y,
                                        double *neum_to_diri_z);


void write_vtk_psd(int N_x, int N_y, int N_z,
		   double *x_coord, double *y_coord, double *z_coord,
		   double *f_old, const char *name);


void write_psd(int N_x, int N_y, int N_z,
	     double *x_coord, double *y_coord, double *z_layers_coord,
	       double *f_old, const char *name);

void PrepareAgglomerationBULK(TCollection *Coll, TFEFunction2D *velocity1, 
                              TFEFunction2D *velocity2, int N_x, int N_y, int N_a,
                              double *x_coord, double *y_coord, 
                              double *a_layers_coord, double *f, 
                              double *rhs_new, int *correspond_2dgrid);

#endif /// __BULK_2D3D__

#endif // __2D__
