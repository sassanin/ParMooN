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
// @(#)Database.h        1.23 06/27/00
//
// Class:       TDatabase
// Purpose:     database of needed refinement, mapping and
//              shape descriptors as well as iterators
//              parameter database
//
// Author:      Volker Behns  24.07.97
//
// =======================================================================
#if defined(_MPI) || defined(_SMPI)
#  include "mpi.h"
#endif

#ifdef _HYBRID
  #include <omp.h>
#endif

#ifndef __DATABASE__
#define __DATABASE__

#include <Iterator.h>
#include <RefDesc.h>
#include <Constants.h>

#include <Domain.h>

struct TParaDB
{
  int VERSION;
  // indicate what kind of problem is solved (T means time dependent)
  //  0: not set
  //  1: CD
  //  2: TCD
  //  3: Stokes
  //  4: TStokes
  //  5: NSE
  //  6: TNSE
  int PROBLEM_TYPE;

  //======================================================================
  /** parameters data output and input files                            */
  //======================================================================
  TDomain *Domain;
  
  //======================================================================
  /** parameters data output and input files                            */
  //======================================================================
  char *GEOFILE;
  char *BNDFILE;
  char *GEOFILE_INTL;
  char *BNDFILE_INTL;
  char *MAPFILE;
  char *OUTFILE;
  char *PODFILE;
   
  int SAVESOL;
  
  char *BASENAME;
  char *VTKBASENAME;
  char *PSBASENAME;
  char *GRAPEBASENAME;
  char *GNUBASENAME;
  char *READGRAPEBASENAME;
  char *GMVBASENAME;
  char *MATLABBASENAME;
  char *OUTPUTDIR;
  char *SAVE_DATA_FILENAME;
  char *READ_DATA_FILENAME;
  char *MATLAB_MATRIX;
  char *SMESHFILE;
  
  char *POD_FILENAME;
  char *SNAP_FILENAME;
  
  //======================================================================
  /** parameters for controling the program */
  //======================================================================
  int PRECOND_LS;
  int SOLVER_TYPE;
  int WRITE_PS;
  int WRITE_GRAPE;
  int WRITE_GNU;
  int WRITE_GMV;
  int WRITE_AMIRA;
  int WRITE_VTK;
  int WRITE_MATLAB;
  int WRITE_MATLAB_MATRIX;
  int SAVE_DATA;
  int READ_DATA;
  int READ_GRAPE_FILE;
  int MEASURE_ERRORS;
  int ESTIMATE_ERRORS;
  int COMPUTE_VORTICITY_DIVERGENCE;
  int MESH_TYPE;
  int USE_PRM;
 
  int timeprofiling;  //(time profiling)
  
    //======================================================================
  /** Time parameters for multigrid                    */
  //======================================================================
  double time_system_assemble;
  
  double time_solve;
  double time_solve_start;
  double time_solve_end;
  
  double time_vanka;
  double time_vanka_start;
  double time_vanka_end;
  
  double time_vanka_solve;
  double time_vanka_solve_start;
  double time_vanka_solve_end;
  
  double time_communication;
  double time_communication_start;
  double time_communication_end;
  
  double time_GMRES;
  double time_GMRES_start;
  double time_GMRES_end;
  
  double time_MG;
  double time_MG_start;
  double time_MG_end;
  
  double time_projection;
  double time_projection_start;
  double time_projection_end;
  
  double time_restriction;
  double time_restriction_start;
  double time_restriction_end;
  
  //======================================================================
  /** parameters for setting finite element spaces                      */
  //======================================================================
  int ANSATZ_ORDER;
  int TEST_ORDER;

  int VELOCITY_SPACE;
  int PRESSURE_SPACE;
  int PRESSURE_SEPARATION;

  int EXAMPLE; // used (and explained) in derived classes of Example2D
  //======================================================================
  /** parameters for grid generation                                    */
  //======================================================================

  int OMPNUMTHREADS;

  int LEVELS;
  int UNIFORM_STEPS;
  double DRIFT_X;
  double DRIFT_Y;
  double DRIFT_Z;
  int NONLINEARIT_TYPE_NEWTON;
  int REFINEMENT;
  int GRID_TYPE;
  int GRID_TYPE_1;
  int GRID_TYPE_2;
  double CHANNEL_GRID_STRETCH;
  int MESHGEN_ALLOW_EDGE_REF;
  int MESHGEN_REF_QUALITY;

  //======================================================================
  /** parameters for adaptive grid refinement                           */
  //======================================================================
  int ADAPTIVE_REFINEMENT_CRITERION;
  int ERROR_CONTROL;
  int REFINE_STRATEGY;
  int MAX_CELL_LEVEL;
  double REFTOL;
  double COARSETOL;
  double MIN_FRACTION_TO_CHANGE;
  double DECREASE_REFTOL_FACTOR;
  double INCREASE_COARSETOL_FACTOR;
  double FRACTION_OF_ERROR;
  int CONVERT_QUAD_TO_TRI;
  int N_CELL_LAYERS;

  //======================================================================
  /** parameters for setting the discretization                         */
  //======================================================================
  // general
  int DISCTYPE;
  int USE_ISOPARAMETRIC;
  int CELL_MEASURE;
  // parameter for non-conforming elements
  int NC_TYPE;
  
  // DISCTYPE for internal space, PBE
  int  INTL_DISCTYPE;
  
  /** upwind methods */
  int UPWIND_ORDER;
  double UPWIND_FLUX_DAMP;
  int UPWIND_APPLICATION;

  // Shishkin meshes
  int SHISHKIN_MESH;
  double SHISHKIN_DIAM;

  /** stabilization with face integrals */
  int WEAK_BC;
  double FACE_SIGMA;
  double WEAK_BC_SIGMA;
  double TAU;
  double TAU2;
  double TAU3;

  /** SUPG or SDFEM method */
  int    SDFEM_TYPE;
  double DELTA0;
  double DELTA1;
  double DELTA2;
  double SDFEM_POWER0;
  int    SDFEM_NORM_B;
  double ADJOINT_FACTOR_4_OMEGA_EQ_0;
  int    CIP_TYPE;
  /** parameters for SOLD methods */
  int SOLD_TYPE;
  int SOLD_PARAMETER_TYPE;
  double SOLD_CONST;
  double SOLD_POWER;
  double SOLD_S;
  double SOLD_U0;
  int SOLD_PARAMETER_SCALING;
  double SOLD_PARAMETER_SCALING_FACTOR;

  //======================================================================
  /** parameters for vectorial FE (Raviart-Thomas, Brezzi-Douglas-Marini) */
  //======================================================================
  int NORMAL_ORIENTATION_QUAD[4];
  int NORMAL_ORIENTATION_TRIA[3];
  int NORMAL_ORIENTATION_TETRA[4];
  int NORMAL_ORIENTATION_HEXA[6];

  //======================================================================
  /** parameter for local projection stabilization */
  //======================================================================
  int LP_FULL_GRADIENT;
  int LP_STREAMLINE;
  int LP_DIVERGENCE;
  int LP_PRESSURE;
  int LP_CROSSWIND;
  int LP_COEFF_TYPE;

  double LP_FULL_GRADIENT_COEFF;
  double LP_STREAMLINE_COEFF;
  double LP_DIVERGENCE_COEFF;
  double LP_PRESSURE_COEFF;

  double LP_FULL_GRADIENT_EXPONENT;
  double LP_STREAMLINE_EXPONENT;
  double LP_DIVERGENCE_EXPONENT;
  double LP_PRESSURE_EXPONENT;

  int LP_ORDER_DIFFERENCE;
  int LP_FULL_GRADIENT_ORDER_DIFFERENCE;
  int LP_STREAMLINE_ORDER_DIFFERENCE;
  int LP_DIVERGENCE_ORDER_DIFFERENCE;
  int LP_PRESSURE_ORDER_DIFFERENCE;
  
  int LP_CROSSWIND_COEFF_TYPE;
  double LP_CROSSWIND_COEFF;
  double LP_CROSSWIND_EXPONENT; 

  //======================================================================
  /** parameter for a posteriori parameter computation with adjoint problem */
  //======================================================================
  int SOLVE_ADJOINT_PROBLEM;
  int SOLD_ADJOINT;
  int N_STAGES_ADJOINT;
  int SC_NONLIN_ITE_ADJOINT;
  int OPTIMIZATION_ITE_TYPE_ADJOINT;
  int BFGS_VECTORS_ADJOINT;
  int PENALTY_ADJOINT;
  double RELATIVE_DECREASE_ADJOINT;
  double PENALTY_VALUE_AT_ZERO_ADJOINT;
  double PENALTY_SMALLEST_PARAM_FAC_ADJOINT;
  double PENALTY_LARGEST_PARAM_FAC_ADJOINT;
  double WEIGHT_RESIDUAL_L1_ADJOINT;
  double WEIGHT_RESIDUAL_L2_ADJOINT;
  double WEIGHT_GRADIENT_L1_ADJOINT;
  double WEIGHT_GRADIENT_L2_ADJOINT;
  double WEIGHT_STREAM_DER_L1_ADJOINT;
  double WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT;
  double WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT;
  double REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT;
  double MIN_VAL_ADJOINT;
  double MAX_VAL_ADJOINT;
  double MIN_MAX_EXPONENT_ONE_ADJOINT;
  double MIN_MAX_EXPONENT_TWO_ADJOINT;
  double MIN_MAX_FACTOR_ONE_ADJOINT;
  double MIN_MAX_FACTOR_TWO_ADJOINT;
  double WEIGHT_RESIDUAL_LP_ADJOINT;
  double WEIGHT_RESIDUAL_EXP_LP_ADJOINT;
  double WEIGHT_RESIDUAL_CW_ADJOINT;
  double RESIDUAL_LP_ADJOINT;
  int MIN_MAX_ADJOINT;
  int INITIAL_STEEPEST_DESCENT_ADJOINT;


  /** parameter for superconvergence */
  int SUPERCONVERGENCE_ORDER;

  /** the following parameters are for the membrane REACTOR */
  double REACTOR_P0;
  double REACTOR_P1;
  double REACTOR_P2;
  double REACTOR_P3;
  double REACTOR_P4;
  double REACTOR_P5;
  double REACTOR_P6;
  double REACTOR_P7;
  double REACTOR_P8;
  double REACTOR_P9;
  double REACTOR_P10;
  double REACTOR_P11;
  double REACTOR_P12;
  double REACTOR_P13;
  double REACTOR_P14;
  double REACTOR_P15;
  double REACTOR_P16;
  double REACTOR_P17;
  double REACTOR_P18;
  double REACTOR_P19;
  double REACTOR_P20;
  double REACTOR_P21;
  double REACTOR_P22;
  double REACTOR_P23;
  double REACTOR_P24;
  double REACTOR_P25;
  double REACTOR_P26;
  double REACTOR_P27;
  double REACTOR_P28;
  double REACTOR_P29;
  double REACTOR_P30;

  /** parameter for linear FEM-FCT scheme */
  int FEM_FCT_LINEAR_TYPE;
  int FEM_FCT_PRELIMITING;

  int FEM_FCT_GROUP_FEM; 
  int GROUP_FEM; 
   /** parameter for WENO scheme */
  int WENO_TYPE;
  
  //======================================================================
  /** PARAMETERS FOR STOKES AND NAVIER-STOKES PROBLEMS                  */
  //======================================================================

  /** general parameters */
  double RE_NR;
  double RA_NR;
  double ROSSBY_NR;
  double START_RE_NR;
  double RE_NR_INCREMENT;
  int FLOW_PROBLEM_TYPE;
  int NSE_NONLINEAR_FORM;
  int NSTYPE;
  int LAPLACETYPE;
  int DEFECT_CORRECTION_TYPE;
  double OSEEN_ZERO_ORDER_COEFF;

  //======================================================================
  /** PARAMETERS FOR DARCY PROBLEM                  */
  //======================================================================
  int DARCYTYPE; 
  double SIGMA_PERM;
  //======================================================================
  
  double FR_NR;
  double WB_NR;
  double PR_NR;
  double PE_NR;
  double BI_NR;
  double WEI_NR;
  int Axial3D;
  int Axial3DAxis;
  
  /** parameters for LES */
  double FILTER_WIDTH_CONSTANT;
  double FILTER_WIDTH_POWER;
  double GAUSSIAN_GAMMA;
  int CONVOLUTE_SOLUTION;

  /** parameters for turbulent viscosity */
  int TURBULENT_VISCOSITY_TYPE;
  int TURBULENT_VISCOSITY_TENSOR;
  double TURBULENT_VISCOSITY_CONSTANT;
  double TURBULENT_VISCOSITY_POWER;
  double TURBULENT_VISCOSITY_SIGMA;
  int TURBULENT_MOD_TYPE;
  
  double viscosity_max;
  double viscosity_min;

  /** parameters for VMS */
  int VMS_LARGE_VELOCITY_SPACE;
  int VMS_COARSE_MG_SMAGO;

  // constants in AdaptProjectionSpace
  double VMS_ADAPT_LOWER;
  double VMS_ADAPT_MIDDLE;
  double VMS_ADAPT_UPPER;
  int VMS_ADAPT_STEPS;
  int VMS_ADAPT_COMP;

  double ARTIFICIAL_VISCOSITY_CONSTANT;
  double ARTIFICIAL_VISCOSITY_POWER;
  int VMM_COARSE_LEVEL;
  int VMM_COARSE_SPACE_ORDER;
  int RFB_SUBMESH_LAYERS;

  //======================================================================
  /** parameters for slip with friction and penetration with resistance
      boundary conditions                                               */
  //======================================================================
  double FRICTION_CONSTANT;
  double FRICTION_POWER;
  int    FRICTION_TYPE;
  double FRICTION_U0;
  double PENETRATION_CONSTANT;
  double PENETRATION_POWER;

  //======================================================================
  /** parameters for div-div stabilization */
  //======================================================================
  int    DIV_DIV_STAB_TYPE;
  double DIV_DIV_STAB_C1;
  double DIV_DIV_STAB_C2;

  //======================================================================
  // ******** parameters for scalar system *********//
  //======================================================================
  // parameters for nonlinear iteration
  int    SC_NONLIN_ITE_TYPE_SCALAR;
  int    SC_NONLIN_MAXIT_SCALAR;
  double SC_NONLIN_RES_NORM_MIN_SCALAR;
  double SC_NONLIN_DAMP_FACTOR_SCALAR;

  // parameters for linear iteration
  int    SC_SOLVER_SCALAR;
  int    SC_PRECONDITIONER_SCALAR;
  int    SC_LIN_MAXIT_SCALAR;
  double SC_LIN_RED_FACTOR_SCALAR;
  double SC_LIN_RES_NORM_MIN_SCALAR;
  int    SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;

  int    SC_LIN_MAXIT_SCALAR_SOLD;
  double SC_LIN_RED_FACTOR_SCALAR_SOLD;
  double SC_LIN_RES_NORM_MIN_SCALAR_SOLD;

  // parameters which are used in multigrid for scalar problems
  int    SC_MG_TYPE_SCALAR;
  int    SC_MG_CYCLE_SCALAR;
  int    SC_SMOOTHER_SCALAR;
  int    SC_PRE_SMOOTH_SCALAR;
  int    SC_POST_SMOOTH_SCALAR;
  double SC_SMOOTH_DAMP_FACTOR_SCALAR;
  double SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
  double SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR;
  int    SC_COARSE_SMOOTHER_SCALAR;
  int    SC_COARSE_MAXIT_SCALAR;
  double SC_COARSE_RED_FACTOR_SCALAR;
  double SC_GMG_DAMP_FACTOR_SCALAR;
  double SC_GMG_DAMP_FACTOR_FINE_SCALAR;
  int    SC_FIRST_SOLUTION_LEVEL_SCALAR;
  int    SC_COARSEST_LEVEL_SCALAR;

  int    SC_STEP_LENGTH_CONTROL_FINE_SCALAR;
  int    SC_STEP_LENGTH_CONTROL_ALL_SCALAR;

  //======================================================================
  // ******** parameters for saddle point system *********//
  //======================================================================
  // parameters for nonlinear iteration
  int    SC_NONLIN_ITE_TYPE_SADDLE;
  int    SC_NONLIN_MAXIT_SADDLE;
  double SC_NONLIN_RES_NORM_MIN_SADDLE;
  double SC_NONLIN_DAMP_FACTOR_SADDLE;
  int    SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE;

  // parameters for linear iteration
  int    SC_SOLVER_SADDLE;
  int    SC_PRECONDITIONER_SADDLE;
  int    SC_LIN_MAXIT_SADDLE;
  double SC_LIN_RED_FACTOR_SADDLE;
  double SC_LIN_RES_NORM_MIN_SADDLE;

  // parameters which are used in multigrid for saddle point problems
  int    SC_MG_TYPE_SADDLE;
  int    SC_MG_CYCLE_SADDLE;
  int    SC_SMOOTHER_SADDLE;
  int    SC_PRE_SMOOTH_SADDLE;
  int    SC_POST_SMOOTH_SADDLE;
  double SC_SMOOTH_DAMP_FACTOR_SADDLE;
  double SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
  double SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
  int    SC_COARSE_SMOOTHER_SADDLE;
  int    SC_COARSE_MAXIT_SADDLE;
  double SC_COARSE_RED_FACTOR_SADDLE;
  double SC_GMG_DAMP_FACTOR_SADDLE;
  double SC_GMG_DAMP_FACTOR_FINE_SADDLE;

  int    SC_FIRST_SOLUTION_LEVEL_SADDLE;
  int    SC_COARSEST_LEVEL_SADDLE;

  int    SC_STEP_LENGTH_CONTROL_FINE_SADDLE;
  int    SC_STEP_LENGTH_CONTROL_ALL_SADDLE;
  int    SC_LARGEST_DIRECT_SOLVE;
  int    SC_DOWNWIND_TYPE;
  
  //======================================================================
  /** AMG solver parameters */
  //======================================================================

  // coarsen context
  double CC_ALPHA;
  double CC_BETA;
  int    CC_MINCLUSTER;
  int    CC_MAXCLUSTER;
  int    CC_MAXDISTANCE;
  int    CC_MAXCONNECTIVITY;
  int    CC_DEPTHTARGET;
  int    CC_COARSENTARGET;
  double CC_COARSENRATE;
  int    CC_MAJOR;
  int    CC_DEPENDENCY;
  double CC_RESCALE;
  int    CC_VERBOSE;

  // solver context
  int    SC_SYSTEM_TYPE;
  int    SC_AMG_PREC_IT;
  double SC_AMG_PREC_RED_FACTOR;
  int    SC_EX_MAXIT;
  int    SC_GMRES_RESTART;
  int    SC_LCD_START_VECTOR;
  double SC_ILU_BETA;
  double SC_SOR_OMEGA;
  double SC_SMOOTHER_RED_FACTOR;
  double SC_OMEGA_COARSE_0;
  double SC_OMEGA_P_0;
  double SC_ILUT_TOL;
  int    SC_ILUT_ABSOLUTE_FILLIN;
  double SC_ILUT_RELATIVE_FILLIN;
  int    SC_ILUT_SORT;
  int    SC_SCHUR_INV_OF_A;
  int    SC_SCHUR_INV_OF_A_MAXIT;
  double SC_SCHUR_ITERATION_DAMP;
  int    SC_SCHUR_ITERATION_MAXIT;
  int    SC_SCHUR_STEP_LENGTH_CONTROL;
  int    SC_MIXED_BCGS_CGS_SWITCH_TOL;
  double SC_DIV_FACTOR;
  double SC_NONLIN_DIV_FACTOR;
  int    SC_SMOOTHING_STEPS;
  int    SC_N1_PARAM;
  int    SC_N2_PARAM;
  int    SC_MINIT;
  double SC_VAS_LAZ_DELTA;
  int    SC_ROW_EQUILIBRATION;
  int    SC_VERBOSE;
  int    SC_VERBOSE_AMG;

  int    SC_BRAESS_SARAZIN_MATRIX;
  double SC_BRAESS_SARAZIN_ALPHA;

  double CHAR_L0;
  double D_VISCOSITY;
  double SURF_TENSION;
  double IMPACT_ANGLE;
  double Area;
  //======================================================================
  /** the following parameters are for individual use */
  //======================================================================
  double P0;
  double P1;
  double P2;
  double P3;
  double P4;
  double P5;
  double P6;
  double P7;
  double P8;
  double P9;
  double P10;
  double P11;
  double P12;
  double P13;
  double P14;
  double P15;
  
  int MG_DEBUG;
  int DOF_Average;
  int DOF_Reorder;
  int SC_LOCAL_SMOOTH;
  //======================================================================
  /** PARAMETERS FOR APPLICATIONS */
  //======================================================================

  int TENSOR_TYPE;                        // conformation stress tensor or deformation tensor in DEVSS
  int FREE_SURFACE_FLOW;                 // Impinging droplet (free surface flow)
  int TWO_PHASE_FLOW;                    // Two Phase flow
  int PHASE1_TYPE;                      // 1 - Newtonian, 2 - Oldroyd, 3- Giesekus
  int PHASE2_TYPE;                      // 1 - Newtonian, 2 - Oldroyd, 3- Giesekus
  
  //======================================================================
  /** parameters for free surface calculation */
  //======================================================================
  bool INTERFACE_FLOW;                             // free surface or Interface flow or not
  int FS_MAGNETLAW;                               // 0 = Langevin, 1 = Vislovich
  double FS_L;                                    // characteristic length
  double FS_U;                                    // characteristic velocity
  double FS_T;                                    // characteristic time: T = L/U

  double FS_ETA;                                  // dynamic viscosity
  double FS_RHO;                                  // fluid density
  double FS_ALPHA;                                // coefficient of surface tension
  double FS_G;                                    // acceleration due to gravity

  double FS_MS;                                   // saturation magnetisation
  double FS_CHI0;                                 // initial susceptibility

  double FS_HM;                                   // mean field strength
  double FS_DELTA_H;                              // amplitude for field oscillations
  double FS_F;                                    // oscillation frequency

  double FS_H0;                                   // characteristic field strength
  double FS_H1;                                   // field inside for plane interface
  double FS_H2;                                   // field outside for plane interface

  double FS_LH;                                   // layer height
  double FS_GAMMA;                                // Langevin parameter
  double FS_HT;                                   // parameter for Vislovich approximation

  double FS_WE;                                   // Weber number = Re*Ca

  char *FS_INNAME;                                // file name for reading surface position
  char *FS_OUTNAME;                               // file name for writing surface position

  int FS_WRITE;                                   // != 0 => write free surface
  int FS_READ;                                    // != 0 => read free surface
 
  double HEAT_TANGENTIAL_STRESS_FACTOR;           // C_1/\sigma_sa
  double HEAT_SOLID_SURFACE_FACTOR;               // eps = (1/PE_NR)HEAT_SOLID_SURFACE_FACTOR
  double EQ_CONTACT_ANGLE;                        // equilibrium contact angle
  double AD_CONTACT_ANGLE;                        // advancing contact angle 
  double RE_CONTACT_ANGLE;                        // receding contact angle 
  double DY_CONTACT_ANGLE;                        // dynamic contact angle 
  int CONTACT_ANGLE_TYPE;                         // type of contact angle   
  
  //======================================================================
  /** parameters for turbulent channel flow */
  //======================================================================
  int CHANNEL_STATISTICS2_WITH_MODEL;

  //======================================================================
  /** parameters for turbulent flow around a squared cylinder */
  //======================================================================
  double CYLINDER_22000_YPLUS_SIDES;
  double CYLINDER_22000_YPLUS_FRONT;
  double CYLINDER_22000_YPLUS_BACK;
  
  //======================================================================
  /** parameters for BULK computations */
  //======================================================================
  int BULK_REACTION_DISC;
  int BULK_PB_DISC;
  int BULK_PB_DISC_STAB;
  int BULK_PB_DISC_FCT_GROUP;
  int BULK_COUPLING;
  int BULK_GROWTH_RATE;
  int BULK_REACTION_MASS_LUMPING;
  int BULK_METHODS_OF_MOMENTS;
  int BULK_MOM_DISC;
  int BULK_SOLD_PARAMETER_TYPE;
  int N_CELL_LAYERS_PSD;
  int N_CELL_LAYERS_PSD_2;
  int OUTPUT_NODE_LAYER_PSD;
  double BULK_REACTION_C_CUT;

  double BULK_l_infty;
  double BULK_u_infty;
  double BULK_c_infty;
  double BULK_c_C_infty_sat;
  double BULK_c_C_infty;
  double BULK_f_infty;

  double BULK_density;
  double BULK_dynamic_viscosity;

  double BULK_C_g;
  double BULK_C_nuc;
  double BULK_C_sat;
  double BULK_C_2;
  double BULK_D_A;
  double BULK_D_P_0;
  double BULK_D_P_MAX;
  double BULK_k_g;
  double BULK_k_r;
  double BULK_k_nuc;
  double BULK_D_P_MIN;

  //======================================================================
  /** parameters for shear slip mesh update method computations */
  //======================================================================
  double SSMUM_MP_X;
  double SSMUM_MP_Y;
  double SSMUM_INNER_RADIUS;
  double SSMUM_OUTER_RADIUS;
  double SSMUM_ROT_PER_SECOND;
  double SSMUM_ANGLE;
  int SSMUM_MAX_CELLS_LAYERS;  
  int SSMUM_INTERPOLATION;

  //======================================================================
  /** parameters for WINDTUNNEL computations */
  //======================================================================
  double WINDTUNNEL_SHIFT; 
  int WINDTUNNEL_CONFIGURATION; 
  int WINDTUNNEL_INTERPOLATION;
  int WINDTUNNEL_STEADY;
  int WINDTUNNEL_SPATIAL;
  double WINDTUNNEL_BROWNIAN ;
  int WINDTUNNEL_POL_ORDER;
  int WINDTUNNEL_SHEAR_FACTOR_TYPE;
  double  WINDTUNNEL_SHEAR_FACTOR;
  int WINDTUNNEL_QUAD_METHOD;
  int WINDTUNNEL_MEASURE_MASS;
  int WINDTUNNEL_LAYER_NUMBER_X;
  int WINDTUNNEL_DIM_Y;
  int WINDTUNNEL_DIM_Z;
  int WINDTUNNEL_DIM_R;
  //double WINDTUNNEL_Y[WINDTUNNEL_DIM_Y_CONST];
  //double WINDTUNNEL_Z[WINDTUNNEL_DIM_Z_CONST];
  double WINDTUNNEL_BOUND_VAL[WINDTUNNEL_DIM_Y_CONST][WINDTUNNEL_DIM_Z_CONST][2];
  double WINDTUNNEL_BOUND_KOEFF[WINDTUNNEL_DIM_Y_CONST][WINDTUNNEL_DIM_Z_CONST];
  double WINDTUNNEL_DROP_VELO[WINDTUNNEL_LAYER_NUMBER_X_CONST][WINDTUNNEL_DIM_Y_CONST][WINDTUNNEL_DIM_Z_CONST];
  double WINDTUNNEL_BOUND_VAL_DROPS[WINDTUNNEL_DIM_Y_CONST][WINDTUNNEL_DIM_Z_CONST][WINDTUNNEL_DIM_R_CONST][2];
  double WINDTUNNEL_ENVIR_COND;
  double WINDTUNNEL_SUPERSAT;
  double WINDTUNNEL_U_INFTY;
  double WINDTUNNEL_L_INFTY;
  double WINDTUNNEL_R_MIN;
  double WINDTUNNEL_R_INFTY;
  double WINDTUNNEL_R_INFTY_EXPERIMENT;
  double WINDTUNNEL_F_INFTY;
  double WINDTUNNEL_kinematic_viscosity;
  double WINDTUNNEL_dynamic_viscosity;
  double WINDTUNNEL_density; 
  //double WINDTUNNEL_BOUND_KOEFF[WINDTUNNEL_DIM_Y_CONST][WINDTUNNEL_DIM_Z_CONST];

  //======================================================================
  /** parameters for urea synthesis computations */
  //======================================================================
  int UREA_REACTION_DISC;
  int UREA_PB_DISC;
  int UREA_MODEL;
  int UREA_PB_DISC_STAB;
  int UREA_SOLD_PARAMETER_TYPE;
  int UREA_PIPE;
  int UREA_CONC_MAXIT;
  double UREA_inflow_time;

  double UREA_l_infty;
  double UREA_u_infty;
  double UREA_c_infty;
  double UREA_temp_infty;
  double UREA_nu;
  double UREA_rho;
  double UREA_c_p;
  double UREA_lambda;
  double UREA_D_P_0;
  double UREA_D_P_MAX;
  double UREA_k_v;
  double UREA_f_infty;
  double UREA_m_mol;
  double UREA_D_J;
  double UREA_rho_d;
  double UREA_k_g;
  double UREA_delta_h_cryst;
  double UREA_g;
  double UREA_rho_sat_1;
  double UREA_rho_sat_2;
  double UREA_alfa_nuc;
  double UREA_beta_nuc;
  double UREA_INFLOW_SCALE;
  double UREA_CONC_TOL;

  double UREA_AGGR_SPATIAL;
  double UREA_AGGR_BROWNIAN;
  double UREA_AGGR_POL_ORDER;
  double UREA_AGGR_SHEAR_FACTOR_TYPE;
  double UREA_AGGR_SHEAR_FACTOR;
  double UREA_AGGR_BROWNIAN_TEMP;
  double UREA_AGGR_BROWNIAN_SCAL;
  double UREA_PIPE_RADIUS;
  int PB_DISC_TYPE;
  int PB_TIME_DISC;
  
  //======================================================================
  /** parameters for kdp synthesis computations */
  //======================================================================
  
  int KDP_MODEL;
  double KDP_D_P_0_2;
  double KDP_D_P_MAX_2;
  double KDP_l_infty;
  double KDP_u_infty;
  double KDP_c_infty;
  double KDP_temp_infty;
  double KDP_nu;
  double KDP_rho;
  double KDP_rho_water;
  double KDP_c_p;
  double KDP_lambda;
  double KDP_D_P_0;
  double KDP_D_P_MAX;
  double KDP_f_infty;
  double KDP_m_mol;
  double KDP_D_J;
  double KDP_rho_d;
  double KDP_k_g_1;
  double KDP_k_g_2;
  double KDP_k_b;
  double KDP_delta_h_cryst;
  double KDP_g_1;
  double KDP_g_2;
  double KDP_b;
  double KDP_w_sat_1;
  double KDP_w_sat_2;
  double KDP_w_sat_3;
  double KDP_w_sat_1_Ma;
  double KDP_w_sat_2_Ma;
  double KDP_w_sat_3_Ma;
  double KDP_INFLOW_SCALE;
  double KDP_CONC_TOL;
  double KDP_INTERNAL_NUC_A;
  double KDP_INTERNAL_NUC_B;

  int KDP_CONC_MAXIT;
  double KDP_inflow_time;
  
  //======================================================================
  /** parameters for Stokes--Darcy (StoDa) coupling */
  //======================================================================
  int StoDa_interfaceType; //Beavers-Joseph-Saffman or u.t=0
  double StoDa_alpha; // from Beavers-Joseph-Saffman condition on interface
  int StoDa_problemType; // Neumann--Neumann, Robin--Robin, ...
  int StoDa_updatingStrategy; // update of the etas
  double StoDa_theta_f; //damping in Stokes (flow) part
  double StoDa_theta_p; //damping in Darcy (porous) part
  double StoDa_gamma_f; // parameter for Robin condition on interface
  double StoDa_gamma_p; // parameter for Robin condition on interface
  double StoDa_weakGamma; // parameter for enforcing weak boundary conditions
  double StoDa_solutionStrategy; // iterative (0), one big matrix (2), both (1)
  int StoDa_algorithm; // Gauss--Seidel, Jacobi, ...
  int StoDa_StokesFirst; // for Gauss--Seidel type method.
  int StoDa_nIterations; // maximum number of iterations
  // convergence criteria: 
  // interface error e = ( ||uS.n-uD.n||^2_L2  +  ||nTn+pD||^2_L2 )^{1/2}
  double StoDa_relDiff_interfaceError; // (e_k - e_{k+1})/e_k < this number
  // E_k^2 = ( a1 * (||uS_{h,k}-uS_{h,k+1})/uS_{h,k} )^2
  //        +( a2 * (||pS_{h,k}-pS_{h,k+1})/pS_{h,k} )^2
  //        +( a3 * (||pD_{h,k}-pD_{h,k+1})/pD_{h,k} )^2
  // a1,a2,a3 are the following
  double StoDa_relDiff_factor1;
  double StoDa_relDiff_factor2;
  double StoDa_relDiff_factor3;
  double StoDa_relDiff_solution; // E_k < this number
  double StoDa_bigResidual; // residual of big System < this number
  int StoDa_periodicBoundary; // true if there is a periodic boundary
  // a prescribed pressure drop at the periodic boundary (to have a flow at all)
  double StoDa_periodicBoundaryPressureDrop; 
  

  //======================================================================
  /** internal parameters
  cannot be set in the readin file
  are used as global variables */
  //======================================================================
  int    INTERNAL_PROBLEM_LINEAR;
  int    INTERNAL_PROJECT_PRESSURE;
  int    INTERNAL_PRESSURE_SPACE;
  int    INTERNAL_SLIP_WITH_FRICTION;
  int    INTERNAL_SLIP_WITH_FRICTION_IDENTITY;
  int    INPUT_QUAD_RULE;
  int    INTERNAL_QUAD_HEXA;
  int    INTERNAL_QUAD_TETRA;
  int    INTERNAL_QUAD_QUAD;
  int    INTERNAL_QUAD_TRIA;
  int    INTERNAL_QUAD_RULE;
  int    INTERNAL_LOCAL_DOF;
  int    INTERNAL_PERIODIC_IDENTITY;
  int    INTERNAL_PROBLEM_IDENTITY;
  int    INTERNAL_LEVEL;
  int    INTERNAL_CONVECTION_EQ_VELOFIELD;
  int    INTERNAL_STEADY_STATE_MATRICES_OR_RHS;
  int    INTERNAL_AMG_SOLVES;
  double INTERNAL_AMG_PREPARE_TIME;
  int    INTERNAL_GMRES_INFO;
  int    INTERNAL_POLYNOMIAL_DEGREE;
  int    INTERNAL_MESH_CELL_TYPE;
  double INTERNAL_BULK_MEAN;
  double INTERNAL_BULK_SIMULATION;
  double INTERNAL_VERTEX_X[8];
  double INTERNAL_VERTEX_Y[8];
  double INTERNAL_VERTEX_Z[8];
  double INTERNAL_HK_CONVECTION;
  int    INTERNAL_MEAN_COMPUTATION;
  int    INTERNAL_MOMENT;
  int    INTERNAL_LINEAR_SCHEME;
  int    INTERNAL_SOLD_ACTIVE;
  int    INTERNAL_UMFPACK_FLAG;
  int    INTERNAL_SORT_AMG;
  int    INTERNAL_FESPACE_CONSTRUCT; 
  int    INTERNAL_DO_NOT_RESPECT_DIRICHLET_BC; 
  
  double INTERNAL_COERCIVITY;
  double *INTERNAL_P1_Array;
  double *INTERNAL_WEIGHT_Array;
  int    *INTERNAL_INDICATOR_Array;
  int    INTERNAL_FACE_INTEGRALS;
  int    INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS;
  int    INTERNAL_WRONG_NEUMANN_CHECKED;
  double INTERNAL_BFGS_RESTART_ADJOINT;
  int    INTERNAL_ARRAY_LENGTH;
  int    INTERNAL_CELL;
  int    INTERNAL_OUTFLOW_BOUNDARY[10];
  double INTERNAL_WEIGHT_SUPG_ADJOINT;
  double INTERNAL_WEIGHT_SOLD_ADJOINT;
  int    INTERNAL_NEW_MATRICES_B;
  int    INTERNAL_FULL_MATRIX_STRUCTURE;
  int    INTERNAL_DISC_FLAG;
  int    INTERNAL_START_PARAM;
  int    MESH_SLIP_WITH_FRICTION;

  // parameter for tetgen
  double TETGEN_QUALITY;
  double TETGEN_VOLUMEN;
  int 	 TETGEN_STEINER;

  //======================================================================
  /** parameters for individual use in parallel computations */
  //======================================================================
  int Par_P0; //out rank
  int Par_P1; // 1 - root takes part in computation; 0 - not
  int Par_P2; // mesh partition type: 1 - dual; 0 - nodal
  int Par_P3; // 1 - use halocells; 0 - dont
  int Par_P4; // 1-redistribution of masters 0-otherwise
  int Par_P5;
  int Par_P6;
  int Par_P7;
  int Par_P8;
  int Par_P9;
  int MapperType;
  int DSType;
  double Par_P10;
  double Par_P11;
  double Par_P12;
  double Par_P13;
  double Par_P14;
  double Par_P15;
  double Par_P16;
  double Par_P17;
  double Par_P18;
  double Par_P19;
  double Par_P20;

  //======================================================================
  /** parameters for ROM */
  //======================================================================
  int WRITE_SNAPSHOTS;
  int DO_ROM;
  int DO_ROM_P;
  int RANK_OF_BASIS;
  int RANK_OF_BASIS_P;
  int POD_INNER_PRODUCT;
  int POD_INNER_PRODUCT_P;
  int BUILD_PODFILE;
  int POD_FLUCT_FIELD;
  int POD_FLUCT_FIELD_P;
  int P_ROM_METHOD;
  
  //======================================================================
  /** parameters for projection methods (NSE) */
  //======================================================================
  int PROJECTION_METHOD;
  
  //======================================================================
  /** parameters for population balance computations */
  //======================================================================
  int PBE_P0;
  int PBE_P1;
  int PBE_P2;
  int PBE_P3;
  int PBE_P4;
  int PBE_P5;
  int PBE_P6;
  int PBE_P7;
  int PBE_P8;
  int PBE_P9;

  //======================================================================
  /** parameters for DG computations */
  //======================================================================
  double DG_P0;
  double DG_P1;
  double DG_P2;
  double DG_P3;
  double DG_P4;
  double DG_P5;
  double DG_P6;
  double DG_P7;
  double DG_P8;
  double DG_P9;
  
  //======================================================================
  /** parameters for moving domains */
  //======================================================================
  
  int MOVING_BOUNDARY;
  double LameC;
  int DEPENDENT_BASIS;
  int DEPENDENT_BASIS_Q1;
  int DEPENDENT_BASIS_Q2;  
  bool ASSEMBLEMESHMAT;
  #if defined(_MPI) || defined(_SMPI)
  /** MPI_Comm for which the computation is started (should not be changed during coomputation)*/
  MPI_Comm Comm;
 #endif
};

typedef struct TParaDB TParamDB;

struct TTimDB
{
  double CURRENTTIME;
  double CURRENTTIMESTEPLENGTH;
  double TIMESTEPLENGTH;
  double INTERNAL_STARTTIME;
  double MIN_TIMESTEPLENGTH;
  double MAX_TIMESTEPLENGTH;
  double TIMESTEPLENGTH_TOL;
  int TIMESTEPLENGTH_CONTROL;
  int TIMESTEPLENGTH_CONTROLLER;  // mlh
  double TIMESTEPLENGTH_PARA_KK_I;
  double TIMESTEPLENGTH_PARA_KK_P;
  double TIMESTEPLENGTH_PARA_KK_E;
  double TIMESTEPLENGTH_PARA_KK_R;
  double TIMESTEPLENGTH_PARA_KK_D;
  double TIMESTEPLENGTH_PARA_FAC;
  double TIMESTEPLENGTH_PARA_FAC_MAX;
  double TIMESTEPLENGTH_PARA_FAC_MIN;
  double TIMESTEPLENGTH_PARA_TOL;
  double TIMESTEPLENGTH_PARA_ATOL;
  double TIMESTEPLENGTH_PARA_RTOL;
  int FIRST_SSC_STEP;
  int RESET_CURRENTTIME;
  double RESET_CURRENTTIME_STARTTIME;
  double STEADY_STATE_TOL;
  double SCALE_DIVERGENCE_CONSTRAINT;

  int CONTROL;
  double CONTROL_ALPHA;
  double CONTROL_BETA;
  double CONTROL_GAMMA;
  double CONTROL_SAFTY;
  double CONTROL_MAXSCALE;
  double CONTROL_MINSCALE;
  
  // control parameter
  double THETA1;
  double THETA2;
  double THETA3;
  double THETA4;

  int TIME_DISC;
  int TIME_DISC2;

  // start and end time
  double STARTTIME;
  double ENDTIME;
  double EXTRAPOLATE_WEIGHT;
  double EXTRAPOLATE_STEPS;
  int    EXTRAPOLATE_PRESSURE;
  // between time steps
  int    EXTRAPOLATE_VELOCITY;

  // parameter for individual use
  double T0;
  double T1;
  double T2;
  double T3;
  double T4;
  double T5;
  double T6;
  double T7;
  double T8;
  double T9;

  // write only every n-th time step
  int STEPS_PER_IMAGE;
  int STEPS_PER_SNAP;

  // parameters for Rosenbrock methods
  int RB_TYPE;
  int RB_TYPE2;
  int RB_SSC;
  double RB_SSC_TOL;
  double RB_SSC_ALPHA;
  double RB_SSC_ALPHA_MAX;
  double RB_SSC_ALPHA_MIN;
  double RB_SSC_MAX_ERROR;

  int RB_APPROX_J;
  int RB_APPROX_C;
  int RB_APPROX_STEPS;

  double RB_GAMMA_I;
  double RB_GAMMA_II;
  double RB_ALPHA_I;
  double RB_SIGMA_I;
  double RB_A_IJ[10];
  double RB_C_IJ[10];
  double RB_S_IJ[10];
  double RB_M_I;
  double RB_MS_I;
  
  // parameters for higher order Galerkin-type methods
  int INTERNAL_SYSTEMSIZE;  
  double *INTERNAL_ALPHA;
  double *INTERNAL_BETA;
  double *VALUE_AT_ONE; 
  double *VAL_AT_QUAD_POINTS;
  double *DER_AT_QUAD_POINTS;
  double *CORR_AT_QUAD_POINTS;  
  double *DER_CORR_AT_QUAD_POINTS;
  double *DER_AT_START;
  double *DER_AT_ONE;
  double *DER_COR_AT_ONE;
  double NORMW;
  double *ZETA;
  double *WEIGHTS;
  int N_QUADPOINTS;
  int N_DEGREES;
  double *ALPHA0;
  double *BETA0;
  double *GAMMA0;
  double *CORRECTION;
  double *POINTS;

  double RK_A[5][5];
  double RK_c[5];
  double RK_b[5];
  double RK_e[5];
  int RK_ord;   // Ordnung des RK-Verfahrens    mlh
  int RK_ord_e;   // Ordnung des eingebetteten RK-Verfahrens  mlh
  
  //dG time steppings
  int DG_TimeDisc;
  int DG_Order;
  
};

typedef struct TTimDB TTimeDB;

/** database of needed refinement, mapping and
    shape descriptors as well as iterators */
class TDatabase
{
  public:
    /** database of shape descriptors */
    static TShapeDesc **ShapeDB;

    /** database of refinement descriptors */
    static TRefDesc **RefDescDB;

    /** database of mapper */
    static TMapper **MapperDB;

    /** database of iterators */
    static TIterator **IteratorDB;

    /** general parameters */
    static TParamDB *ParamDB;

    /** parameter for time discretization */
    static TTimeDB *TimeDB;

  public:
    // Constructors
    /** initialize the database */
    TDatabase();

    // Methods
#ifdef __MORTAR__
    /** add descriptor for mortar refinement with base edge 0 */
    void AddMortar0(int Mortar_Ni, int N);
    /** add descriptor for mortar refinement with base edge 1 */
    void AddMortar1(int Mortar_Ni, int N);
#endif

    // set default parameters

    static void SetDefaultParameters();

    static void WriteParamDB(char *ExecutedFile);

    static void WriteTimeDB();

    static void CheckParameterConsistencyNSE();
};
#endif
