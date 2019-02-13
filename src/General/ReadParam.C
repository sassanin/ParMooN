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
// @(#)ReadParam.C        1.28 06/27/00
//
// Purpose:     read parameter file and
//              read boundary paramerization
//.
// Author:      Volker Behns  22.07.97
// Version:     1.0
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <BdCircle.h>
#include <BdLine.h>
#include <BdSpline.h>
#include <BdPolygon.h>
#include <BdNonUniformSpline.h>
#include <Domain.h>
#include <Database.h>
#include <Joint.h>
#include <MacroCell.h>
#include <MortarBaseJoint.h>
#include <MooNMD_Io.h>
#include <fstream>
#include <string.h>

#ifdef __3D__
#include <BdNoPRM.h>
#include <BdPlane.h>
#include <BdWall.h>
#include <BdSphere.h>
#endif

// #include <amg_solve_main.h>

#include <stdlib.h>

int TDomain::ReadParam(char *ParamFile)
{
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  char line[100], *aux_char;
  int N_Param = 0, flag[7];
  std::ifstream dat(ParamFile);

  memset(flag,0,7*SizeOfInt);

  TDatabase::SetDefaultParameters();

  if (!dat)
  {
#ifdef _MPI
  if(rank==0)
#endif
    cerr << "cannot open '" << ParamFile << "' for input" << endl;
    exit(-1);
  }

  while (!dat.eof())
  {
    dat >> line;

    if (!strcmp(line, "VERSION:"))
    {
      dat >> TDatabase::ParamDB->VERSION;
      N_Param++;
    }
    if (!strcmp(line, "PROBLEM_TYPE:"))
    {
      dat >> TDatabase::ParamDB->PROBLEM_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "EXAMPLE:"))
    {
      dat >> TDatabase::ParamDB->EXAMPLE;
      N_Param++;
    }
      if (!strcmp(line, "MG_DEBUG:"))
    {
      dat >> TDatabase::ParamDB->MG_DEBUG;
      N_Param++;
    }
    
          if (!strcmp(line, "SC_LOCAL_SMOOTH:"))
    {
      dat >> TDatabase::ParamDB->SC_LOCAL_SMOOTH;
      N_Param++;
    }
    
    
    
    if (!strcmp(line, "DOF_Reorder:"))
    {
      dat >> TDatabase::ParamDB->DOF_Reorder;
      N_Param++;
    }
    
    if (!strcmp(line, "DOF_Average:"))
    {
      dat >> TDatabase::ParamDB->DOF_Average;
      N_Param++;
    }
    
    if (!strcmp(line, "GEOFILE:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->GEOFILE = aux_char;
      N_Param++;
    }
        if (!strcmp(line, "GEOFILE_INTL:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->GEOFILE_INTL = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "BNDFILE:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->BNDFILE = aux_char;
      N_Param++;
    }

        if (!strcmp(line, "BNDFILE_INTL:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->BNDFILE_INTL = aux_char;
      N_Param++;
    }
    
    if (!strcmp(line, "MAPFILE:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->MAPFILE = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "OUTFILE:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->OUTFILE = aux_char;
      N_Param++;
    }
    if (!strcmp(line, "SAVESOL:"))
    {
      dat >> TDatabase::ParamDB->SAVESOL;
      N_Param++;
    }
    
    if (!strcmp(line, "TETGEN_QUALITY:"))
    {
      dat >> TDatabase::ParamDB->TETGEN_QUALITY;
      N_Param++;
    }
    
    if (!strcmp(line, "TETGEN_VOLUMEN:"))
    {
      dat >> TDatabase::ParamDB->TETGEN_VOLUMEN;
      N_Param++;
    }

    if (!strcmp(line, "TETGEN_STEINER:"))
    {
      dat >> TDatabase::ParamDB->TETGEN_STEINER;
      N_Param++;
    }

    if (!strcmp(line, "MESHGEN_ALLOW_EDGE_REF:"))
    {
      dat >> TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
      N_Param++;
    }

    if (!strcmp(line, "MESHGEN_REF_QUALITY:"))
    {
      dat >> TDatabase::ParamDB->MESHGEN_REF_QUALITY;
      N_Param++;
    }

    if (!strcmp(line, "ANSATZ_ORDER:"))
    {
      dat >> TDatabase::ParamDB->ANSATZ_ORDER;
      N_Param++;
    }

    if (!strcmp(line, "TEST_ORDER:"))
    {
      dat >> TDatabase::ParamDB->TEST_ORDER;
      N_Param++;
    }

    if (!strcmp(line, "VELOCITY_SPACE:"))
    {
      dat >> TDatabase::ParamDB->VELOCITY_SPACE;
      N_Param++;
    }

    if (!strcmp(line, "PRESSURE_SPACE:"))
    {
      dat >> TDatabase::ParamDB->PRESSURE_SPACE;
      N_Param++;
    }

    if (!strcmp(line, "PRESSURE_SEPARATION:"))
    {
      dat >> TDatabase::ParamDB->PRESSURE_SEPARATION;
      N_Param++;
    }

    if (!strcmp(line, "LEVELS:"))
    {
      dat >> TDatabase::ParamDB->LEVELS;
      N_Param++;
    }

    if (!strcmp(line, "UNIFORM_STEPS:"))
    {
      dat >> TDatabase::ParamDB->UNIFORM_STEPS;
      N_Param++;
    }
    if (!strcmp(line, "DRIFT_X:"))
    {
      dat >> TDatabase::ParamDB->DRIFT_X;
      N_Param++;
    }
    if (!strcmp(line, "DRIFT_Y:"))
    {
      dat >> TDatabase::ParamDB->DRIFT_Y;
      N_Param++;
    }
    if (!strcmp(line, "DRIFT_Z:"))
    {
      dat >> TDatabase::ParamDB->DRIFT_Z;
      N_Param++;
    }
    
    if (!strcmp(line, "NONLINEARIT_TYPE_NEWTON:"))
    {
      dat >> TDatabase::ParamDB->NONLINEARIT_TYPE_NEWTON;
      N_Param++;
    }

    if (!strcmp(line, "REFINEMENT:"))
    {
      dat >> TDatabase::ParamDB->REFINEMENT;
      N_Param++;
    }

    if (!strcmp(line, "GRID_TYPE:"))
    {
      dat >> TDatabase::ParamDB->GRID_TYPE;
      N_Param++;
    }

   if (!strcmp(line, "ADAPTIVE_REFINEMENT_CRITERION:"))
    {
      dat >> TDatabase::ParamDB->ADAPTIVE_REFINEMENT_CRITERION;
      N_Param++;
    }

    if (!strcmp(line, "ERROR_CONTROL:"))
    {
      dat >> TDatabase::ParamDB->ERROR_CONTROL;
      N_Param++;
    }

    if (!strcmp(line, "REFINE_STRATEGY:"))
    {
      dat >> TDatabase::ParamDB->REFINE_STRATEGY;
      N_Param++;
    }
    if (!strcmp(line, "REFTOL:"))
    {
      dat >> TDatabase::ParamDB->REFTOL;
      N_Param++;
    }
    if (!strcmp(line, "COARSETOL:"))
    {
      dat >> TDatabase::ParamDB->COARSETOL;
      N_Param++;
    }
    if (!strcmp(line, "MIN_FRACTION_TO_CHANGE:"))
    {
      dat >> TDatabase::ParamDB->MIN_FRACTION_TO_CHANGE;
      N_Param++;
    }
    if (!strcmp(line, "DECREASE_REFTOL_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->DECREASE_REFTOL_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "INCREASE_COARSETOL_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->INCREASE_COARSETOL_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "FRACTION_OF_ERROR:"))
    {
      dat >> TDatabase::ParamDB->FRACTION_OF_ERROR;
      N_Param++;
    }
    if (!strcmp(line, "MAX_CELL_LEVEL:"))
    {
      dat >> TDatabase::ParamDB->MAX_CELL_LEVEL;
      N_Param++;
    }

    if (!strcmp(line, "CONVERT_QUAD_TO_TRI:"))
    {
      dat >> TDatabase::ParamDB->CONVERT_QUAD_TO_TRI;
      N_Param++;
    }

    if (!strcmp(line, "N_CELL_LAYERS:"))
    {
      dat >> TDatabase::ParamDB->N_CELL_LAYERS;
      N_Param++;
    }
   if (!strcmp(line, "CHANNEL_GRID_STRETCH:"))
    {
      dat >> TDatabase::ParamDB->CHANNEL_GRID_STRETCH;
      N_Param++;
    }

    if (!strcmp(line, "DISCTYPE:"))
    {
      dat >> TDatabase::ParamDB->DISCTYPE;
      N_Param++;
    }
    if (!strcmp(line, "INTL_DISCTYPE:"))
    {
      dat >> TDatabase::ParamDB->INTL_DISCTYPE;
      N_Param++;
    }   
    
    if (!strcmp(line, "NSTYPE:"))
    {
      dat >> TDatabase::ParamDB->NSTYPE;
      N_Param++;
    }
    if (!strcmp(line, "LAPLACETYPE:"))
    {
      dat >> TDatabase::ParamDB->LAPLACETYPE;
      N_Param++;
    }
     if (!strcmp(line, "DARCYTYPE:"))
     {
       dat >> TDatabase::ParamDB->DARCYTYPE;
       N_Param++;
    }
      if (!strcmp(line, "TENSOR_TYPE:"))
     {
       dat >> TDatabase::ParamDB->TENSOR_TYPE;
       N_Param++;
    }
      if (!strcmp(line, "TWO_PHASE_FLOW:"))
     {
       dat >> TDatabase::ParamDB->TWO_PHASE_FLOW;
       N_Param++;
    }
          if (!strcmp(line, "FREE_SURFACE_FLOW:"))
     {
       dat >> TDatabase::ParamDB->FREE_SURFACE_FLOW;
       N_Param++;
    }
          if (!strcmp(line, "PHASE1_TYPE:"))
     {
       dat >> TDatabase::ParamDB->PHASE1_TYPE;
       N_Param++;
    }
          if (!strcmp(line, "PHASE2_TYPE:"))
     {
       dat >> TDatabase::ParamDB->PHASE2_TYPE;
       N_Param++;
    }
    if (!strcmp(line, "SIGMA_PERM:"))
    {
      dat >> TDatabase::ParamDB->SIGMA_PERM;
      N_Param++;  
    }
    if (!strcmp(line, "USE_ISOPARAMETRIC:"))
    {
      dat >> TDatabase::ParamDB->USE_ISOPARAMETRIC;
      N_Param++;
    }
    if (!strcmp(line, "VMM_COARSE_LEVEL:"))
    {
      dat >> TDatabase::ParamDB->VMM_COARSE_LEVEL;
      N_Param++;
    }
    if (!strcmp(line, "VMM_COARSE_SPACE_ORDER:"))
    {
      dat >> TDatabase::ParamDB->VMM_COARSE_SPACE_ORDER;
      N_Param++;
    }

    if (!strcmp(line, "UPWIND_ORDER:"))
    {
      dat >> TDatabase::ParamDB->UPWIND_ORDER;
      N_Param++;
    }

    if (!strcmp(line, "RFB_SUBMESH_LAYERS:"))
    {
      dat >> TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
      N_Param++;
    }
    if (!strcmp(line, "DEFECT_CORRECTION_TYPE:"))
    {
      dat >> TDatabase::ParamDB->DEFECT_CORRECTION_TYPE;
      N_Param++;
    }

    if (!strcmp(line, "CELL_MEASURE:"))
    {
	dat >> TDatabase::ParamDB->CELL_MEASURE;
      N_Param++;
    }

    if (!strcmp(line, "SAMARSKI_DAMP:"))
    {
      dat >> TDatabase::ParamDB->UPWIND_FLUX_DAMP; 
      N_Param++;
    }

    if (!strcmp(line, "UPWIND_APPLICATION:"))
    {
      dat >> TDatabase::ParamDB->UPWIND_APPLICATION;
      N_Param++;
    }

    if (!strcmp(line, "SHISHKIN_MESH:"))
    {
      dat >> TDatabase::ParamDB->SHISHKIN_MESH;
      N_Param++;
    }

    if (!strcmp(line, "SHISHKIN_DIAM:"))
    {
      dat >> TDatabase::ParamDB->SHISHKIN_DIAM;
      N_Param++;
    }
    if (!strcmp(line, "RE_NR:"))
    {
      dat >> TDatabase::ParamDB->RE_NR;
      N_Param++;
    }
    if (!strcmp(line, "RA_NR:"))
    {
      dat >> TDatabase::ParamDB->RA_NR;
      N_Param++;
    }
    if (!strcmp(line, "ROSSBY_NR:"))
    {
      dat >> TDatabase::ParamDB->ROSSBY_NR;
      N_Param++;
    }
    if (!strcmp(line, "START_RE_NR:"))
    {
      dat >> TDatabase::ParamDB->START_RE_NR;
      N_Param++;
    }
    if (!strcmp(line, "RE_NR_INCREMENT:"))
    {
      dat >> TDatabase::ParamDB->RE_NR_INCREMENT;
      N_Param++;
    }
    if (!strcmp(line, "FLOW_PROBLEM_TYPE:"))
    {
      dat >> TDatabase::ParamDB->FLOW_PROBLEM_TYPE;
      N_Param++;
    }

    if (!strcmp(line, "FACE_SIGMA:"))
    {
      dat >> TDatabase::ParamDB->FACE_SIGMA;
      N_Param++;
    }
      if (!strcmp(line, "WEAK_BC_SIGMA:"))
    {
      dat >> TDatabase::ParamDB->WEAK_BC_SIGMA;
      N_Param++;
    }

      if (!strcmp(line, "WEAK_BC:"))
    {
      dat >> TDatabase::ParamDB->WEAK_BC;
      N_Param++;
    }

     if (!strcmp(line, "FR_NR:"))
    {
      dat >> TDatabase::ParamDB->FR_NR;
      N_Param++;
    }
     if (!strcmp(line, "WB_NR:"))
    {
      dat >> TDatabase::ParamDB->WB_NR;
      N_Param++;
    }

     if (!strcmp(line, "PR_NR:"))
    {
      dat >> TDatabase::ParamDB->PR_NR;
      N_Param++;
    }
      if (!strcmp(line, "PE_NR:"))
    {
      dat >> TDatabase::ParamDB->PE_NR;
      N_Param++;
    }   
    
     if (!strcmp(line, "BI_NR:"))
    {
      dat >> TDatabase::ParamDB->BI_NR;
      N_Param++;
    }
         if (!strcmp(line, "WEI_NR:"))
    {
      dat >> TDatabase::ParamDB->WEI_NR;
      N_Param++;
    }
        if (!strcmp(line, "LameC:"))
    {
      dat >> TDatabase::ParamDB->LameC;
      N_Param++;
    }
     if (!strcmp(line, "Axial3D:"))
    {
      dat >> TDatabase::ParamDB->Axial3D;
      N_Param++;
    }
     if (!strcmp(line, "Axial3DAxis:"))
    {
      dat >> TDatabase::ParamDB->Axial3DAxis;
      N_Param++;
    }


      if (!strcmp(line, "TAU:"))
    {
      dat >> TDatabase::ParamDB->TAU;
      N_Param++;
    }

    if (!strcmp(line, "DELTA1:"))
    {
      dat >> TDatabase::ParamDB->DELTA1;
      N_Param++;
    }

    if (!strcmp(line, "DELTA0:"))
    {
      dat >> TDatabase::ParamDB->DELTA0;
      N_Param++;
    }
    
        if (!strcmp(line, "DELTA2:"))
    {
      dat >> TDatabase::ParamDB->DELTA2;
      N_Param++;
    }
    
    
    if (!strcmp(line, "SDFEM_POWER0:"))
    {
      dat >> TDatabase::ParamDB->SDFEM_POWER0;
      N_Param++;
    }

    if (!strcmp(line, "SDFEM_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SDFEM_TYPE;
      N_Param++;
    }

    if (!strcmp(line, "SDFEM_NORM_B:"))
    {
      dat >> TDatabase::ParamDB->SDFEM_NORM_B;
      N_Param++;
    }

    if (!strcmp(line, "FILTER_WIDTH_CONSTANT:"))
    {
      dat >> TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
      N_Param++;
    }
    if (!strcmp(line, "FILTER_WIDTH_POWER:"))
    {
      dat >> TDatabase::ParamDB->FILTER_WIDTH_POWER;
      N_Param++;
    }
    if (!strcmp(line, "GAUSSIAN_GAMMA:"))
    {
      dat >> TDatabase::ParamDB->GAUSSIAN_GAMMA;
      N_Param++;
    }
    if (!strcmp(line, "CONVOLUTE_SOLUTION:"))
    {
      dat >> TDatabase::ParamDB->CONVOLUTE_SOLUTION;
      N_Param++;
    }
    if (!strcmp(line, "TURBULENT_VISCOSITY_TYPE:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "TURBULENT_VISCOSITY_TENSOR:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
      N_Param++;
    }
    if (!strcmp(line, "TURBULENT_VISCOSITY_CONSTANT:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
      N_Param++;
    }
    if (!strcmp(line, "TURBULENT_VISCOSITY_POWER:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      N_Param++;
    }
    if (!strcmp(line, "TURBULENT_VISCOSITY_SIGMA:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
      N_Param++;
    }
    
    if (!strcmp(line, "TURBULENT_MOD_TYPE:"))
    {
      dat >> TDatabase::ParamDB->TURBULENT_MOD_TYPE;
      N_Param++;
    }
    
    
    if (!strcmp(line, "ARTIFICIAL_VISCOSITY_CONSTANT:"))
    {
      dat >> TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT;
      N_Param++;
    }
    if (!strcmp(line, "ARTIFICIAL_VISCOSITY_POWER:"))
    {
      dat >> TDatabase::ParamDB->ARTIFICIAL_VISCOSITY_POWER;
      N_Param++;
    }
    if (!strcmp(line, "FRICTION_CONSTANT:"))
    {
      dat >> TDatabase::ParamDB->FRICTION_CONSTANT;
      N_Param++;
    }
    if (!strcmp(line, "FRICTION_POWER:"))
    {
      dat >> TDatabase::ParamDB->FRICTION_POWER;
      N_Param++;
    }
    if (!strcmp(line, "FRICTION_U0:"))
    {
      dat >> TDatabase::ParamDB->FRICTION_U0;
      N_Param++;
    }
    if (!strcmp(line, "FRICTION_TYPE:"))
    {
      dat >> TDatabase::ParamDB->FRICTION_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "PENETRATION_CONSTANT:"))
    {
      dat >> TDatabase::ParamDB->PENETRATION_CONSTANT;
      N_Param++;
    }
    if (!strcmp(line, "PENETRATION_POWER:"))
    {
      dat >> TDatabase::ParamDB->PENETRATION_POWER;
      N_Param++;
    }

    if (!strcmp(line, "NSE_NONLINEAR_FORM:"))
    {
      dat >> TDatabase::ParamDB->NSE_NONLINEAR_FORM;
      N_Param++;
    }
    if (!strcmp(line, "OSEEN_ZERO_ORDER_COEFF:"))
    {
      dat >> TDatabase::ParamDB->OSEEN_ZERO_ORDER_COEFF;
      N_Param++;
    }

    if (!strcmp(line, "DIV_DIV_STAB_TYPE:"))
    {
      dat >> TDatabase::ParamDB->DIV_DIV_STAB_TYPE;
      N_Param++;
    }
    
    if (!strcmp(line, "DIV_DIV_STAB_C1:"))
    {
      dat >> TDatabase::ParamDB->DIV_DIV_STAB_C1;
      N_Param++;
    }
    if (!strcmp(line, "DIV_DIV_STAB_C2:"))
    {
      dat >> TDatabase::ParamDB->DIV_DIV_STAB_C2;
      N_Param++;
    }

    if (!strcmp(line, "LP_FULL_GRADIENT:"))
    {
      dat >> TDatabase::ParamDB->LP_FULL_GRADIENT;
      N_Param++;
    }

    if (!strcmp(line, "LP_STREAMLINE:"))
    {
      dat >> TDatabase::ParamDB->LP_STREAMLINE;
      N_Param++;
    }

    if (!strcmp(line, "LP_DIVERGENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_DIVERGENCE;
      N_Param++;
    }

    if (!strcmp(line, "LP_PRESSURE:"))
    {
      dat >> TDatabase::ParamDB->LP_PRESSURE;
      N_Param++;
    }

    if (!strcmp(line, "LP_COEFF_TYPE:"))
    {
      dat >> TDatabase::ParamDB->LP_COEFF_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "LP_FULL_GRADIENT_COEFF:"))
    {
      dat >> TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF;
      N_Param++;
    }

    if (!strcmp(line, "LP_STREAMLINE_COEFF:"))
    {
      dat >> TDatabase::ParamDB->LP_STREAMLINE_COEFF;
      N_Param++;
    }

    if (!strcmp(line, "LP_DIVERGENCE_COEFF:"))
    {
      dat >> TDatabase::ParamDB->LP_DIVERGENCE_COEFF;
      N_Param++;
    }

    if (!strcmp(line, "LP_PRESSURE_COEFF:"))
    {
      dat >> TDatabase::ParamDB->LP_PRESSURE_COEFF;
      N_Param++;
    }

    if (!strcmp(line, "LP_FULL_GRADIENT_EXPONENT:"))
    {
      dat >> TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT;
      N_Param++;
    }

    if (!strcmp(line, "LP_STREAMLINE_EXPONENT:"))
    {
      dat >> TDatabase::ParamDB->LP_STREAMLINE_EXPONENT;
      N_Param++;
    }

    if (!strcmp(line, "LP_DIVERGENCE_EXPONENT:"))
    {
      dat >> TDatabase::ParamDB->LP_DIVERGENCE_EXPONENT;
      N_Param++;
    }

    if (!strcmp(line, "LP_PRESSURE_EXPONENT:"))
    {
      dat >> TDatabase::ParamDB->LP_PRESSURE_EXPONENT;
      N_Param++;
    }

    if (!strcmp(line, "LP_ORDER_DIFFERENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_ORDER_DIFFERENCE;
      N_Param++;
    }

    if (!strcmp(line, "BASENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->BASENAME = aux_char;
      N_Param++;
    }
    if (!strcmp(line, "VTKBASENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->VTKBASENAME = aux_char;
      N_Param++;
    }    
    if (!strcmp(line, "OUTPUTDIR:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->OUTPUTDIR = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "LP_FULL_GRADIENT_ORDER_DIFFERENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE;
      N_Param++;
    }

    if (!strcmp(line, "LP_STREAMLINE_ORDER_DIFFERENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE;
      N_Param++;
    }

    if (!strcmp(line, "LP_DIVERGENCE_ORDER_DIFFERENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE;
      N_Param++;
    }

    if (!strcmp(line, "LP_PRESSURE_ORDER_DIFFERENCE:"))
    {
      dat >> TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE;
      N_Param++;
    }
    if (!strcmp(line, "LP_CROSSWIND:"))
    {
      dat >> TDatabase::ParamDB->LP_CROSSWIND;
      N_Param++;
    }

    if (!strcmp(line, "LP_CROSSWIND_COEFF_TYPE:"))
    {
      dat >> TDatabase::ParamDB->LP_CROSSWIND_COEFF_TYPE;
      N_Param++;
    }

    if (!strcmp(line, "LP_CROSSWIND_COEFF:"))
    {
      dat >> TDatabase::ParamDB->LP_CROSSWIND_COEFF;
      N_Param++;
    }
    
    if (!strcmp(line, "LP_CROSSWIND_EXPONENT:"))
    {
      dat >> TDatabase::ParamDB->LP_CROSSWIND_EXPONENT;
      N_Param++;
    }


   if (!strcmp(line, "SMESHFILE:"))                                                                            
    {                                                                                                           
      dat >> line;                                                                                              
      aux_char = new char[strlen(line) + 1];                                                                    
      strcpy(aux_char, line);                                                                                   
      TDatabase::ParamDB->SMESHFILE = aux_char;                                                                 
      N_Param++;                                                                                                
    }

    if (!strcmp(line, "SAVE_DATA_FILENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->SAVE_DATA_FILENAME = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "READ_DATA_FILENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->READ_DATA_FILENAME = aux_char;
      N_Param++;
    }
    if (!strcmp(line, "POD_FILENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->POD_FILENAME = aux_char;
      N_Param++;
    }
    
    if (!strcmp(line, "SNAP_FILENAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->SNAP_FILENAME = aux_char;
      N_Param++;
    }


    if (!strcmp(line, "SOLD_TYPE:"))
    {
	dat >> TDatabase::ParamDB->SOLD_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_PARAMETER_TYPE:"))
    {
	dat >> TDatabase::ParamDB->SOLD_PARAMETER_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_CONST:"))
    {
	dat >> TDatabase::ParamDB->SOLD_CONST;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_POWER:"))
    {
	dat >> TDatabase::ParamDB->SOLD_POWER;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_S:"))
    {
	dat >> TDatabase::ParamDB->SOLD_S;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_U0:"))
    {
	dat >> TDatabase::ParamDB->SOLD_U0;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_PARAMETER_SCALING:"))
    {
	dat >> TDatabase::ParamDB->SOLD_PARAMETER_SCALING;
      N_Param++;
    }
        if (!strcmp(line, "SOLD_PARAMETER_SCALING_FACTOR:"))
    {
	dat >> TDatabase::ParamDB->SOLD_PARAMETER_SCALING_FACTOR;
      N_Param++;
    }

    if (!strcmp(line, "PRECOND_LS:"))
    {
      dat >> TDatabase::ParamDB->PRECOND_LS;
      N_Param++;
    }

    if (!strcmp(line, "SOLVER_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SOLVER_TYPE;
      N_Param++;
    }

    if (!strcmp(line, "WRITE_PS:"))
    {
      dat >> TDatabase::ParamDB->WRITE_PS;
      N_Param++;
    }

    if (!strcmp(line, "WRITE_GRAPE:"))
    {
      dat >> TDatabase::ParamDB->WRITE_GRAPE;
      N_Param++;
    }

    if (!strcmp(line, "WRITE_GMV:"))
    {
      dat >> TDatabase::ParamDB->WRITE_GMV;
      N_Param++;
    }

    if (!strcmp(line, "WRITE_VTK:"))
    {
      dat >> TDatabase::ParamDB->WRITE_VTK;
      N_Param++;
    }
    
    if (!strcmp(line, "MESH_TYPE:"))
    {
      dat >> TDatabase::ParamDB->MESH_TYPE;
      N_Param++;
    }
    
    if (!strcmp(line, "USE_PRM:"))
    {
      dat >> TDatabase::ParamDB->USE_PRM;
      N_Param++;
      
      if(TDatabase::ParamDB->MESH_TYPE==0)
       {
        //cout<<"USE_PRM=1 for MESH_TYPE=0 !!!!!!!!!!!!" <endl;
        TDatabase::ParamDB->USE_PRM=1;
       }
    }
    
    if (!strcmp(line, "WRITE_AMIRA:"))
    {
      dat >> TDatabase::ParamDB->WRITE_AMIRA;
      N_Param++;
    }
    
    if (!strcmp(line, "WRITE_SNAPSHOTS:"))
    {
      dat >> TDatabase::ParamDB->WRITE_SNAPSHOTS;
      N_Param++;
    }
    
    if (!strcmp(line, "DO_ROM:"))
   {
    dat >> TDatabase::ParamDB->DO_ROM;
      N_Param++;
   }
    
    if (!strcmp(line, "DO_ROM_P:"))
    {
     dat >> TDatabase::ParamDB->DO_ROM_P;
      N_Param++;
   }
    
    if (!strcmp(line, "RANK_OF_BASIS:"))
    {
      dat >> TDatabase::ParamDB->RANK_OF_BASIS;
      N_Param++;
    }
    
    if (!strcmp(line, "RANK_OF_BASIS_P:"))
    {
      dat >> TDatabase::ParamDB->RANK_OF_BASIS_P;
      N_Param++;
    }
    
    if (!strcmp(line, "POD_INNER_PRODUCT:"))
    {
      dat >> TDatabase::ParamDB->POD_INNER_PRODUCT;
      N_Param++;
    }
    
    if (!strcmp(line, "POD_INNER_PRODUCT_P:"))
    {
      dat >> TDatabase::ParamDB->POD_INNER_PRODUCT_P;
      N_Param++;
    }

    if( !strcmp(line, "PODFILE:" ))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->PODFILE = aux_char;

      N_Param++;
    }
    if( !strcmp(line, "BUILD_PODFILE:" ))
    {
      dat >> TDatabase::ParamDB->BUILD_PODFILE;
      N_Param++;
    }
    
    if( !strcmp(line, "POD_FLUCT_FIELD:" ))
    {
      dat >> TDatabase::ParamDB->POD_FLUCT_FIELD;
      N_Param++;
    }
    
    if( !strcmp(line, "POD_FLUCT_FIELD_P:" ))
    {
      dat >> TDatabase::ParamDB->POD_FLUCT_FIELD_P;
      N_Param++;
    }
    
    if( !strcmp(line, "P_ROM_METHOD:" ))
    {
      dat >> TDatabase::ParamDB->P_ROM_METHOD;
      N_Param++;
    }
    
    if (!strcmp(line, "PROJECTION_METHOD:"))
    {
      dat >> TDatabase::ParamDB->PROJECTION_METHOD;
      N_Param++;
    }

    if (!strcmp(line, "SAVE_DATA:"))
    {
	dat >> TDatabase::ParamDB->SAVE_DATA;
      N_Param++;
    }
    if (!strcmp(line, "READ_DATA:"))
    {
	dat >> TDatabase::ParamDB->READ_DATA;
      N_Param++;
    }

    if (!strcmp(line, "READ_GRAPE_FILE:"))
    {
      dat >> TDatabase::ParamDB->READ_GRAPE_FILE;
      N_Param++;
    }
    if (!strcmp(line, "WRITE_GNU:"))
    {
      dat >> TDatabase::ParamDB->WRITE_GNU;
      N_Param++;
    }

    if (!strcmp(line, "MEASURE_ERRORS:"))
    {
      dat >> TDatabase::ParamDB->MEASURE_ERRORS;
      N_Param++;
    }

    if (!strcmp(line, "ESTIMATE_ERRORS:"))
    {
      dat >> TDatabase::ParamDB->ESTIMATE_ERRORS;
      N_Param++;
    }
    if (!strcmp(line, "SOLVE_ADJOINT_PROBLEM:"))
    {
      dat >> TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM;
      N_Param++;
    }

    if (!strcmp(line, "COMPUTE_VORTICITY_DIVERGENCE:"))
    {
      dat >> TDatabase::ParamDB->COMPUTE_VORTICITY_DIVERGENCE;
      N_Param++;
    }
    if (!strcmp(line, "P0:"))
    {
      dat >> TDatabase::ParamDB->P0;
      N_Param++;
    }

    if (!strcmp(line, "P1:"))
    {
      dat >> TDatabase::ParamDB->P1;
      N_Param++;
    }

    if (!strcmp(line, "P2:"))
    {
      dat >> TDatabase::ParamDB->P2;
      N_Param++;
    }

    if (!strcmp(line, "P3:"))
    {
      dat >> TDatabase::ParamDB->P3;
      N_Param++;
    }

    if (!strcmp(line, "P4:"))
    {
      dat >> TDatabase::ParamDB->P4;
      N_Param++;
    }

    if (!strcmp(line, "P5:"))
    {
      dat >> TDatabase::ParamDB->P5;
      N_Param++;
    }

    if (!strcmp(line, "P6:"))
    {
      dat >> TDatabase::ParamDB->P6;
      N_Param++;
    }

    if (!strcmp(line, "P7:"))
    {
      dat >> TDatabase::ParamDB->P7;
      N_Param++;
    }

    if (!strcmp(line, "P8:"))
    {
      dat >> TDatabase::ParamDB->P8;
      N_Param++;
    }

    if (!strcmp(line, "P9:"))
    {
      dat >> TDatabase::ParamDB->P9;
      N_Param++;
    }
    if (!strcmp(line, "P10:"))
    {
      dat >> TDatabase::ParamDB->P10;
      N_Param++;
    }

    if (!strcmp(line, "P11:"))
    {
      dat >> TDatabase::ParamDB->P11;
      N_Param++;
    }

    if (!strcmp(line, "P12:"))
    {
      dat >> TDatabase::ParamDB->P12;
      N_Param++;
    }

    if (!strcmp(line, "P13:"))
    {
      dat >> TDatabase::ParamDB->P13;
      N_Param++;
    }

    if (!strcmp(line, "P14:"))
    {
      dat >> TDatabase::ParamDB->P14;
      N_Param++;
    }

    if (!strcmp(line, "P15:"))
    {
      dat >> TDatabase::ParamDB->P15;
      N_Param++;
    }

    if (!strcmp(line, "CC_ALPHA:"))
    {
      dat >> TDatabase::ParamDB->CC_ALPHA;
      N_Param++;
    }

    if (!strcmp(line, "CC_BETA:"))
    {
      dat >> TDatabase::ParamDB->CC_BETA;
      N_Param++;
    }

    if (!strcmp(line, "CC_MINCLUSTER:"))
    {
      dat >> TDatabase::ParamDB->CC_MINCLUSTER;
      N_Param++;
    }

    if (!strcmp(line, "CC_MAXCLUSTER:"))
    {
      dat >> TDatabase::ParamDB->CC_MAXCLUSTER;
      N_Param++;
    }

    if (!strcmp(line, "CC_MAXDISTANCE:"))
    {
      dat >> TDatabase::ParamDB->CC_MAXDISTANCE;
      N_Param++;
    }

    if (!strcmp(line, "CC_MAXCONNECTIVITY:"))
    {
      dat >> TDatabase::ParamDB->CC_MAXCONNECTIVITY;
      N_Param++;
    }
    if (!strcmp(line, "CC_DEPTHTARGET:"))
    {
      dat >> TDatabase::ParamDB->CC_DEPTHTARGET;
      N_Param++;
    }
    if (!strcmp(line, "CC_COARSENTARGET:"))
    {
      dat >> TDatabase::ParamDB->CC_COARSENTARGET;
      N_Param++;
    }
    if (!strcmp(line, "CC_COARSENRATE:"))
    {
      dat >> TDatabase::ParamDB->CC_COARSENRATE;
      N_Param++;
    }
    if (!strcmp(line, "CC_MAJOR:"))
    {
      dat >> TDatabase::ParamDB->CC_MAJOR;
      N_Param++;
    }
    if (!strcmp(line, "CC_DEPENDENCY:"))
    {
      dat >> TDatabase::ParamDB->CC_DEPENDENCY;
      N_Param++;
    }
    if (!strcmp(line, "CC_RESCALE:"))
    {
      dat >> TDatabase::ParamDB->CC_RESCALE;
      N_Param++;
    }
    if (!strcmp(line, "CC_VERBOSE:"))
    {
      dat >> TDatabase::ParamDB->CC_VERBOSE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SYSTEM_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SC_SYSTEM_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SOLVER:"))
    {
      dat >> TDatabase::ParamDB->SC_SOLVER_SCALAR;
      TDatabase::ParamDB->SC_SOLVER_SADDLE =
        TDatabase::ParamDB->SC_SOLVER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRECONDITIONER:"))
    {
      dat >> TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE =
         TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_AMG_PREC_IT:"))
    {
      dat >> TDatabase::ParamDB->SC_AMG_PREC_IT;
      N_Param++;
    }
    if (!strcmp(line, "SC_AMG_PREC_RED_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_AMG_PREC_RED_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "SC_EX_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_EX_MAXIT;
      N_Param++;
    }
    if (!strcmp(line, "SC_LCD_START_VECTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_LCD_START_VECTOR;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMRES_RESTART:"))
    {
      dat >> TDatabase::ParamDB->SC_GMRES_RESTART;
      N_Param++;
    }
    if (!strcmp(line, "SC_ILU_BETA:"))
    {
      dat >> TDatabase::ParamDB->SC_ILU_BETA;
      N_Param++;
    }
    if (!strcmp(line, "SC_SOR_OMEGA:"))
    {
      dat >> TDatabase::ParamDB->SC_SOR_OMEGA;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTHER_RED_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTHER_RED_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "SC_OMEGA_COARSE_0:"))
    {
      dat >> TDatabase::ParamDB->SC_OMEGA_COARSE_0;
      N_Param++;
    }
    if (!strcmp(line, "SC_OMEGA_P_0:"))
    {
      dat >> TDatabase::ParamDB->SC_OMEGA_P_0;
      N_Param++;
    }
    if (!strcmp(line, "SC_ILUT_TOL:"))
    {
      dat >> TDatabase::ParamDB->SC_ILUT_TOL;
      N_Param++;
    }
    if (!strcmp(line, "SC_ILUT_ABSOLUTE_FILLIN:"))
    {
      dat >> TDatabase::ParamDB->SC_ILUT_ABSOLUTE_FILLIN;
      N_Param++;
    }
    if (!strcmp(line, "SC_ILUT_RELATIVE_FILLIN:"))
    {
      dat >> TDatabase::ParamDB->SC_ILUT_RELATIVE_FILLIN;
      N_Param++;
    }
    if (!strcmp(line, "SC_ILUT_SORT:"))
    {
      dat >> TDatabase::ParamDB->SC_ILUT_SORT;
      N_Param++;
    }
    if (!strcmp(line, "SC_SCHUR_INV_OF_A:"))
    {
      dat >> TDatabase::ParamDB->SC_SCHUR_INV_OF_A;
      N_Param++;
    }
    if (!strcmp(line, "SC_SCHUR_INV_OF_A_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_SCHUR_INV_OF_A_MAXIT;
      N_Param++;
    }
    if (!strcmp(line, "SC_SCHUR_ITERATION_DAMP:"))
    {
      dat >> TDatabase::ParamDB->SC_SCHUR_ITERATION_DAMP;
      N_Param++;
    }
    if (!strcmp(line, "SC_SCHUR_ITERATION_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_SCHUR_ITERATION_MAXIT;
      N_Param++;
    }
    if (!strcmp(line, "SC_SCHUR_STEP_LENGTH_CONTROL:"))
    {
      dat >> TDatabase::ParamDB->SC_SCHUR_STEP_LENGTH_CONTROL;
      N_Param++;
    }
    if (!strcmp(line, "SC_MIXED_BCGS_CGS_SWITCH_TOL:"))
    {
      dat >> TDatabase::ParamDB->SC_MIXED_BCGS_CGS_SWITCH_TOL;
      N_Param++;
    }
    if (!strcmp(line, "SC_DIV_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_DIV_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_DIV_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_DIV_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTHING_STEPS:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTHING_STEPS;
      N_Param++;
    }
    if (!strcmp(line, "SC_N1_PARAM:"))
    {
      dat >> TDatabase::ParamDB->SC_N1_PARAM;
      N_Param++;
    }
    if (!strcmp(line, "SC_N2_PARAM:"))
    {
      dat >> TDatabase::ParamDB->SC_N2_PARAM;
      N_Param++;
    }
    if (!strcmp(line, "SC_MINIT:"))
    {
      dat >> TDatabase::ParamDB->SC_MINIT;
      N_Param++;
    }
    if (!strcmp(line, "SC_VAS_LAZ_DELTA:"))
    {
      dat >> TDatabase::ParamDB->SC_VAS_LAZ_DELTA;
      N_Param++;
    }
    if (!strcmp(line, "SC_ROW_EQUILIBRATION:"))
    {
      dat >> TDatabase::ParamDB->SC_ROW_EQUILIBRATION;
      N_Param++;
    }
    if (!strcmp(line, "SC_VERBOSE:"))
    {
      dat >> TDatabase::ParamDB->SC_VERBOSE;
      N_Param++;
    }
    if (!strcmp(line, "SC_VERBOSE_AMG:"))
    {
      dat >> TDatabase::ParamDB->SC_VERBOSE_AMG;
      N_Param++;
    }

    if (!strcmp(line, "SC_BRAESS_SARAZIN_MATRIX:"))
    {
      dat >> TDatabase::ParamDB->SC_BRAESS_SARAZIN_MATRIX;
      N_Param++;
    }

    if (!strcmp(line, "SC_BRAESS_SARAZIN_ALPHA:"))
    {
      dat >> TDatabase::ParamDB->SC_BRAESS_SARAZIN_ALPHA;
      N_Param++;
    }

    if (!strcmp(line, "CHAR_L0:"))
    {
      dat >> TDatabase::ParamDB->CHAR_L0;
      N_Param++;
    }
    if (!strcmp(line, "D_VISCOSITY:"))
    {
      dat >> TDatabase::ParamDB->D_VISCOSITY;
      N_Param++;
    }
    if (!strcmp(line, "SURF_TENSION:"))
    {
      dat >> TDatabase::ParamDB->SURF_TENSION;
      N_Param++;
    }

    if (!strcmp(line, "IMPACT_ANGLE:"))
    {
      dat >> TDatabase::ParamDB->IMPACT_ANGLE;
      N_Param++;
    }

    if (!strcmp(line, "Area:"))
    {
      dat >> TDatabase::ParamDB->Area;
      N_Param++;
    }


    // *********** PARAMETERS FOR SCALAR SOLVER ***********

    if (!strcmp(line, "SC_NONLIN_ITE_TYPE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_MAXIT_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_RES_NORM_MIN_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_DAMP_FACTOR_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SOLVER_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_SOLVER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRECONDITIONER_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_LIN_MAXIT_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      N_Param++;
      if (!flag[6])
         TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR_SOLD = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
    }
    if (!strcmp(line, "SC_LIN_RED_FACTOR_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
      N_Param++;
      if (!flag[4])
         TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    }
     if (!strcmp(line, "SC_LIN_RES_NORM_MIN_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      N_Param++;
      if (!flag[5])
         TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    }

    if (!strcmp(line, "SC_LIN_RED_FACTOR_SCALAR_SOLD:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD;
      N_Param++;
      flag[4] = 1;
    }
     if (!strcmp(line, "SC_LIN_RES_NORM_MIN_SCALAR_SOLD:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD;
      N_Param++;
      flag[5] = 1;
    }
    if (!strcmp(line, "SC_LIN_MAXIT_SCALAR_SOLD:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR_SOLD;
      flag[6] = 1;
    }
    if (!strcmp(line, "SC_FLEXIBLE_KRYLOV_SPACE_SOLVER:"))
    {
      dat >> TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;
      N_Param++;
    }

    if (!strcmp(line, "CIP_TYPE:"))
    {
      dat >> TDatabase::ParamDB->CIP_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "SOLD_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->SOLD_ADJOINT;
      N_Param++;
    }

    if (!strcmp(line, "N_STAGES_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->N_STAGES_ADJOINT;
      N_Param++;
    }


    if (!strcmp(line, "SC_NONLIN_ITE_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_ITE_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "ADJOINT_FACTOR_4_OMEGA_EQ_0:"))
    {
      dat >> TDatabase::ParamDB->ADJOINT_FACTOR_4_OMEGA_EQ_0;
      N_Param++;
    }
    if (!strcmp(line, "OPTIMIZATION_ITE_TYPE_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->OPTIMIZATION_ITE_TYPE_ADJOINT;
      N_Param++;
    }

    if (!strcmp(line, "BFGS_VECTORS_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->BFGS_VECTORS_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "PENALTY_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->PENALTY_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "PENALTY_VALUE_AT_ZERO_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->PENALTY_VALUE_AT_ZERO_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "PENALTY_SMALLEST_PARAM_FAC_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->PENALTY_SMALLEST_PARAM_FAC_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "PENALTY_LARGEST_PARAM_FAC_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->PENALTY_LARGEST_PARAM_FAC_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "RELATIVE_DECREASE_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->RELATIVE_DECREASE_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_RESIDUAL_L1_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_RESIDUAL_L2_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_GRADIENT_L1_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_GRADIENT_L1_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_GRADIENT_L2_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_GRADIENT_L2_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_STREAM_DER_L1_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_MAX_EXPONENT_ONE_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_MAX_EXPONENT_TWO_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_MAX_FACTOR_ONE_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_MAX_FACTOR_TWO_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_MAX_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_MAX_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MIN_VAL_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MIN_VAL_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "MAX_VAL_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->MAX_VAL_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_RESIDUAL_LP_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_RESIDUAL_EXP_LP_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "RESIDUAL_LP_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->RESIDUAL_LP_ADJOINT;
      N_Param++;
    }
    if (!strcmp(line, "WEIGHT_RESIDUAL_CW_ADJOINT:"))
    {
      dat >> TDatabase::ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT;
      N_Param++;
    }

    

    if (!strcmp(line, "SC_MG_TYPE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_MG_CYCLE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_CYCLE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTHER_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTHER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRE_SMOOTH_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_POST_SMOOTH_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
      flag[0] =1;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR;
      flag[0] =2;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_SMOOTHER_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_MAXIT_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_RED_FACTOR_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
      flag[1]=1;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR_FINE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSEST_LEVEL_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_FIRST_SOLUTION_LEVEL_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
      N_Param++;
    }

    if (!strcmp(line, "OMPNUMTHREADS:"))
    {
      dat >> TDatabase::ParamDB->OMPNUMTHREADS;
      N_Param++;
    } 

    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_FINE_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_ALL_SCALAR:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR;
      N_Param++;
    }

    // *********** PARAMETERS FOR SADDLE POINT SOLVER ***********

    if (!strcmp(line, "SC_NONLIN_ITE_TYPE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_MAXIT_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_RES_NORM_MIN_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_DAMP_FACTOR_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SOLVER_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_SOLVER_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRECONDITIONER_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_PRECONDITIONER_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_LIN_MAXIT_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_LIN_RED_FACTOR_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
      N_Param++;
    }
     if (!strcmp(line, "SC_LIN_RES_NORM_MIN_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
      N_Param++;
    }

    if (!strcmp(line, "SC_MG_TYPE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_TYPE_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_MG_CYCLE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_CYCLE_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTHER_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTHER_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRE_SMOOTH_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_POST_SMOOTH_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE;
      flag[2] =1;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE;
      flag[2] =2;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_SMOOTHER_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_MAXIT_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_RED_FACTOR_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR_FINE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE;
      flag[3]=1;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSEST_LEVEL_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_FIRST_SOLUTION_LEVEL_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE;
      N_Param++;
    } 
    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_FINE_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE;
      N_Param++;
    }
    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_ALL_SADDLE:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE;
      N_Param++;
    }

    if (!strcmp(line, "SC_LARGEST_DIRECT_SOLVE:"))
    {
      dat >> TDatabase::ParamDB->SC_LARGEST_DIRECT_SOLVE;
      N_Param++;
    }
    if (!strcmp(line, "SC_DOWNWIND_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SC_DOWNWIND_TYPE;
      N_Param++;
    }
    // *********** PARAMETERS FOR BOTH SOLVERS SOLVER ***********

    if (!strcmp(line, "SC_NONLIN_ITE_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SCALAR;
      TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE =
        TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE =
        TDatabase::ParamDB->SC_NONLIN_MAXIT_SCALAR; 
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_RES_NORM_MIN:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
      TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE
        = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_NONLIN_DAMP_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE
        = TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR;
     N_Param++;
    }
    if (!strcmp(line, "SC_LIN_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE
        = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_LIN_RED_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE
        = TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
      N_Param++;
    }
     if (!strcmp(line, "SC_LIN_RES_NORM_MIN:"))
    {
      dat >> TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE 
        = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
      N_Param++;
    }

    if (!strcmp(line, "SC_MG_TYPE:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
      TDatabase::ParamDB->SC_MG_TYPE_SADDLE 
        = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_MG_CYCLE:"))
    {
      dat >> TDatabase::ParamDB->SC_MG_CYCLE_SCALAR;
      TDatabase::ParamDB->SC_MG_CYCLE_SADDLE
        =  TDatabase::ParamDB->SC_MG_CYCLE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTHER:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTHER_SCALAR;
      TDatabase::ParamDB->SC_SMOOTHER_SADDLE
        =  TDatabase::ParamDB->SC_SMOOTHER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_PRE_SMOOTH:"))
    {
      dat >> TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;
      TDatabase::ParamDB->SC_PRE_SMOOTH_SADDLE
        =  TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_POST_SMOOTH:"))
    {
      dat >> TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;
      TDatabase::ParamDB->SC_POST_SMOOTH_SADDLE
        = TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE =
        TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_SMOOTH_DAMP_FACTOR_FINE:"))
    {
      dat >> TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE
        = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
      flag[0] = flag[2] = 1;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_SMOOTHER:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR;
      TDatabase::ParamDB->SC_COARSE_SMOOTHER_SADDLE
        = TDatabase::ParamDB->SC_COARSE_SMOOTHER_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_MAXIT:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
      TDatabase::ParamDB->SC_COARSE_MAXIT_SADDLE
        = TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSE_RED_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SADDLE
        = TDatabase::ParamDB->SC_COARSE_RED_FACTOR_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
      TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE
        = TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
      flag[1]= flag[3] = 1;
      N_Param++;
    }
    if (!strcmp(line, "SC_GMG_DAMP_FACTOR_FINE:"))
    {
      dat >> TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR;
      TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE
        =  TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_COARSEST_LEVEL:"))
    {
      dat >> TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;
      TDatabase::ParamDB->SC_COARSEST_LEVEL_SADDLE =
        TDatabase::ParamDB->SC_COARSEST_LEVEL_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_FIRST_SOLUTION_LEVEL:"))
    {
      dat >> TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
      TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE
        = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
      N_Param++;
    } 
    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_FINE:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR;
      TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE
        = TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR;
      N_Param++;
    }
    if (!strcmp(line, "SC_STEP_LENGTH_CONTROL_ALL:"))
    {
      dat >> TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR;
      TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE
        = TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR;
      N_Param++;
    }


    // read in parameter for time discretization
    if (!strcmp(line, "STEPLENGTH:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH;
      N_Param++;
    }
    if (!strcmp(line, "MIN_TIMESTEPLENGTH:"))
    {
      dat >> TDatabase::TimeDB->MIN_TIMESTEPLENGTH;
      N_Param++;
    }
    if (!strcmp(line, "MAX_TIMESTEPLENGTH:"))
    {
      dat >> TDatabase::TimeDB->MAX_TIMESTEPLENGTH;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_TOL:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_TOL;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_CONTROL:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_CONTROL;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_CONTROLLER:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_CONTROLLER;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_KK_I:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_I;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_KK_P:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_P;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_KK_E:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_E;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_KK_R:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_R;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_KK_D:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_D;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_FAC:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_FAC_MAX:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_FAC_MIN:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_TOL:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_TOL;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_ATOL:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_ATOL;
      N_Param++;
    }
    if (!strcmp(line, "TIMESTEPLENGTH_PARA_RTOL:"))
    {
      dat >> TDatabase::TimeDB->TIMESTEPLENGTH_PARA_RTOL;
      N_Param++;
    }

    if (!strcmp(line, "STARTTIME:"))
    {
      dat >> TDatabase::TimeDB->STARTTIME;
      N_Param++;
    }
    if (!strcmp(line, "ENDTIME:"))
    {
      dat >> TDatabase::TimeDB->ENDTIME;
      N_Param++;
    }
    if (!strcmp(line, "RESET_CURRENTTIME:"))
    {
      dat >> TDatabase::TimeDB->RESET_CURRENTTIME;
      N_Param++;
    }
    if (!strcmp(line, "RESET_CURRENTTIME_STARTTIME:"))
    {
      dat >> TDatabase::TimeDB->RESET_CURRENTTIME_STARTTIME;
      N_Param++;
    }
    if (!strcmp(line, "STEADY_STATE_TOL:"))
    {
      dat >> TDatabase::TimeDB->STEADY_STATE_TOL;
      N_Param++;
    }
    if (!strcmp(line, "SCALE_DIVERGENCE_CONSTRAINT:"))
    {
      dat >> TDatabase::TimeDB->SCALE_DIVERGENCE_CONSTRAINT;
      N_Param++;
    }
    if (!strcmp(line, "EXTRAPOLATE_WEIGHT:"))
    {
      dat >> TDatabase::TimeDB->EXTRAPOLATE_WEIGHT;
      N_Param++;
    }
    if (!strcmp(line, "EXTRAPOLATE_STEPS:"))
    {
      dat >> TDatabase::TimeDB->EXTRAPOLATE_STEPS;
      N_Param++;
    }
    if (!strcmp(line, "EXTRAPOLATE_PRESSURE:"))
    {
      dat >> TDatabase::TimeDB->EXTRAPOLATE_PRESSURE;
      N_Param++;
    }
    if (!strcmp(line, "EXTRAPOLATE_VELOCITY:"))
    {
      dat >> TDatabase::TimeDB->EXTRAPOLATE_VELOCITY;
      N_Param++;
    }
    if (!strcmp(line, "TIME_DISC:"))
    {
      dat >> TDatabase::TimeDB->TIME_DISC;
      N_Param++;
    }
    if (!strcmp(line, "TIME_DISC2:"))
    {
      dat >> TDatabase::TimeDB->TIME_DISC2;
      N_Param++;
    }
    if (!strcmp(line, "FIRST_SSC_STEP:"))
    {
      dat >> TDatabase::TimeDB->FIRST_SSC_STEP;
      N_Param++;
    }
    if (!strcmp(line, "T0:"))
    {
      dat >> TDatabase::TimeDB->T0;
      N_Param++;
    }
    if (!strcmp(line, "T1:"))
    {
      dat >> TDatabase::TimeDB->T1;
      N_Param++;
    }
    if (!strcmp(line, "T2:"))
    {
      dat >> TDatabase::TimeDB->T2;
      N_Param++;
    }
    if (!strcmp(line, "T3:"))
    {
      dat >> TDatabase::TimeDB->T3;
      N_Param++;
    }
    if (!strcmp(line, "T4:"))
    {
      dat >> TDatabase::TimeDB->T4;
      N_Param++;
    }
    if (!strcmp(line, "T5:"))
    {
      dat >> TDatabase::TimeDB->T5;
      N_Param++;
    }
    if (!strcmp(line, "T6:"))
    {
      dat >> TDatabase::TimeDB->T6;
      N_Param++;
    }
    if (!strcmp(line, "T7:"))
    {
      dat >> TDatabase::TimeDB->T7;
      N_Param++;
    }
    if (!strcmp(line, "T8:"))
    {
      dat >> TDatabase::TimeDB->T8;
      N_Param++;
    }
    if (!strcmp(line, "T9:"))
    {
      dat >> TDatabase::TimeDB->T9;
      N_Param++;
    }
    if (!strcmp(line, "STEPS_PER_IMAGE:"))
    {
      dat >> TDatabase::TimeDB->STEPS_PER_IMAGE;
      N_Param++;
    }
    if (!strcmp(line, "STEPS_PER_SNAP:"))
    {
      dat >> TDatabase::TimeDB->STEPS_PER_SNAP;
      N_Param++;
    }
    if (!strcmp(line, "DG_TimeDisc:"))
    {
      dat >> TDatabase::TimeDB->DG_TimeDisc;
      N_Param++;
    }    
    
     if (!strcmp(line, "DG_Order:"))
    {
      dat >> TDatabase::TimeDB->DG_Order;
      N_Param++;
    }   
    if (!strcmp(line, "RB_TYPE:"))
    {
      dat >> TDatabase::TimeDB->RB_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "RB_TYPE2:"))
    {
      dat >> TDatabase::TimeDB->RB_TYPE2;
      N_Param++;
    }
    if (!strcmp(line, "STEPSIZECONTROL:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC;
      N_Param++;
    }
    if (!strcmp(line, "RB_SSC_TOL:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC_TOL;
      N_Param++;
    }
    if (!strcmp(line, "RB_SSC_ALPHA:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC_ALPHA;
      N_Param++;
    }
    if (!strcmp(line, "RB_SSC_ALPHA_MAX:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC_ALPHA_MAX;
      N_Param++;
    }
    if (!strcmp(line, "RB_SSC_ALPHA_MIN:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC_ALPHA_MIN;
      N_Param++;
    }
    if (!strcmp(line, "RB_SSC_MAX_ERROR:"))
    {
      dat >> TDatabase::TimeDB->RB_SSC_MAX_ERROR;
      N_Param++;
    }
    if (!strcmp(line, "RB_APPROX_J:"))
    {
      dat >> TDatabase::TimeDB->RB_APPROX_J;
      N_Param++;
    }
    if (!strcmp(line, "RB_APPROX_C:"))
    {
      dat >> TDatabase::TimeDB->RB_APPROX_C;
      N_Param++;
    }
    if (!strcmp(line, "RB_APPROX_STEPS:"))
    {
      dat >> TDatabase::TimeDB->RB_APPROX_STEPS;
      N_Param++;
    }

    // read in parameter for surface calculations
    if (!strcmp(line, "FS_MAGNETLAW:"))
    {
      dat >> TDatabase::ParamDB->FS_MAGNETLAW;
      N_Param++;
    }
    if (!strcmp(line, "FS_ETA:"))
    {
      dat >> TDatabase::ParamDB->FS_ETA;
      N_Param++;
    }
    if (!strcmp(line, "FS_RHO:"))
    {
      dat >> TDatabase::ParamDB->FS_RHO;
      N_Param++;
    }
    if (!strcmp(line, "FS_ALPHA:"))
    {
      dat >> TDatabase::ParamDB->FS_ALPHA;
      N_Param++;
    }
    if (!strcmp(line, "FS_G:"))
    {
      dat >> TDatabase::ParamDB->FS_G;
      N_Param++;
    }
    if (!strcmp(line, "FS_T:"))
    {
      dat >> TDatabase::ParamDB->FS_T;
      N_Param++;
    }
    if (!strcmp(line, "FS_HM:"))
    {
      dat >> TDatabase::ParamDB->FS_HM;
      N_Param++;
    }
    if (!strcmp(line, "FS_DELTA_H:"))
    {
      dat >> TDatabase::ParamDB->FS_DELTA_H;
      N_Param++;
    }
    if (!strcmp(line, "FS_F:"))
    {
      dat >> TDatabase::ParamDB->FS_F;
      N_Param++;
    }
    if (!strcmp(line, "FS_LH:"))
    {
      dat >> TDatabase::ParamDB->FS_LH;
      N_Param++;
    }
    if (!strcmp(line, "FS_MS:"))
    {
      dat >> TDatabase::ParamDB->FS_MS;
      N_Param++;
    }
    if (!strcmp(line, "FS_CHI0:"))
    {
      dat >> TDatabase::ParamDB->FS_CHI0;
      N_Param++;
    }
    if (!strcmp(line, "FS_WRITE:"))
    {
      dat >> TDatabase::ParamDB->FS_WRITE;
      N_Param++;
    }
    if (!strcmp(line, "FS_READ:"))
    {
      dat >> TDatabase::ParamDB->FS_READ;
      N_Param++;
    }

    if (!strcmp(line, "FS_INNAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->FS_INNAME = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "FS_OUTNAME:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->FS_OUTNAME = aux_char;
      N_Param++;
    }

    if (!strcmp(line, "HEAT_TANGENTIAL_STRESS_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR;
      N_Param++;
    }

    if (!strcmp(line, "HEAT_SOLID_SURFACE_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->HEAT_SOLID_SURFACE_FACTOR;
      N_Param++;
    }

    if (!strcmp(line, "EQ_CONTACT_ANGLE:"))
    {
      dat >> TDatabase::ParamDB->EQ_CONTACT_ANGLE;
      N_Param++;
    }

    if (!strcmp(line, "AD_CONTACT_ANGLE:"))
    {
      dat >> TDatabase::ParamDB->AD_CONTACT_ANGLE;
      N_Param++;
    }
    
    if (!strcmp(line, "RE_CONTACT_ANGLE:"))
    {
      dat >> TDatabase::ParamDB->RE_CONTACT_ANGLE;
      N_Param++;
    }
    
        
    if (!strcmp(line, "DY_CONTACT_ANGLE:"))
    {
      dat >> TDatabase::ParamDB->DY_CONTACT_ANGLE;
      N_Param++;
    }
 
    if (!strcmp(line, "CONTACT_ANGLE_TYPE:"))
    {
      dat >> TDatabase::ParamDB->CONTACT_ANGLE_TYPE;
      N_Param++;
    }   

    if (!strcmp(line, "VMS_SMALL_VELOCITY_SPACE:"))
    {
      dat >> TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE;
#ifdef _MPI
  if(rank==0)
#endif
      OutPut("The name of this parameter is now VMS_LARGE_VELOCITY_SPACE !!!" << endl);
      N_Param++;
    }

    if (!strcmp(line, "VMS_LARGE_VELOCITY_SPACE:"))
    {
      dat >> TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE;
      N_Param++;
    }

    if (!strcmp(line, "VMS_COARSE_MG_SMAGO:"))
    {
	dat >> TDatabase::ParamDB->VMS_COARSE_MG_SMAGO;
      N_Param++;
    }
    if (!strcmp(line, "VMS_ADAPT_LOWER:"))
    {
      dat >> TDatabase::ParamDB->VMS_ADAPT_LOWER;
      N_Param++;
    }
    if (!strcmp(line, "VMS_ADAPT_MIDDLE:"))
    {
      dat >> TDatabase::ParamDB->VMS_ADAPT_MIDDLE;
      N_Param++;
    }
    if (!strcmp(line, "VMS_ADAPT_UPPER:"))
    {
      dat >> TDatabase::ParamDB->VMS_ADAPT_UPPER;
      N_Param++;
    }
    if (!strcmp(line, "VMS_ADAPT_STEPS:"))
    {
      dat >> TDatabase::ParamDB->VMS_ADAPT_STEPS;
      N_Param++;
    }
    if (!strcmp(line, "VMS_ADAPT_COMP:"))
    {
      dat >> TDatabase::ParamDB->VMS_ADAPT_COMP;
      N_Param++;
    }
    if (!strcmp(line, "SUPERCONVERGENCE_ORDER:"))
    {
      dat >> TDatabase::ParamDB->SUPERCONVERGENCE_ORDER;
      N_Param++;
    }
    if (!strcmp(line, "REACTOR_P0:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P0;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P1:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P1;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P2:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P2;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P3:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P3;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P4:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P4;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P5:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P5;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P6:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P6;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P7:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P7;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P8:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P8;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P9:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P9;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P10:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P10;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P11:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P11;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P12:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P12;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P13:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P13;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P14:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P14;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P15:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P15;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P16:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P16;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P17:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P17;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P18:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P18;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P19:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P19;
      N_Param++;
    }
    
    if (!strcmp(line, "REACTOR_P20:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P20;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P21:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P21;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P22:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P22;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P23:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P23;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P24:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P24;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P25:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P25;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P26:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P26;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P27:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P27;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P28:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P28;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P29:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P29;
      N_Param++;
    }

    if (!strcmp(line, "REACTOR_P30:"))
    {
      dat >> TDatabase::ParamDB->REACTOR_P30;
      N_Param++;
    }

    if (!strcmp(line, "FEM_FCT_LINEAR_TYPE:"))
    {
      dat >> TDatabase::ParamDB->FEM_FCT_LINEAR_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "FEM_FCT_PRELIMITING:"))
    {
      dat >> TDatabase::ParamDB->FEM_FCT_PRELIMITING;
      N_Param++;
    }
    if (!strcmp(line, "FEM_FCT_GROUP_FEM:"))
    {
      dat >> TDatabase::ParamDB->FEM_FCT_GROUP_FEM;
      N_Param++;
    }
      
    if (!strcmp(line, "GROUP_FEM:"))
    {
      dat >> TDatabase::ParamDB->GROUP_FEM;
      N_Param++;
    }
    
    if (!strcmp(line, "WENO_TYPE:"))
    {
      dat >> TDatabase::ParamDB->WENO_TYPE;
      N_Param++;
    }


   if (!strcmp(line, "CHANNEL_STATISTICS2_WITH_MODEL:"))
    {
	dat >> TDatabase::ParamDB->CHANNEL_STATISTICS2_WITH_MODEL;
      N_Param++;
    }
   if (!strcmp(line, "CYLINDER_22000_YPLUS_SIDES:"))
    {
	dat >> TDatabase::ParamDB->CYLINDER_22000_YPLUS_SIDES;
      N_Param++;
    }
   if (!strcmp(line, "CYLINDER_22000_YPLUS_FRONT:"))
    {
	dat >> TDatabase::ParamDB->CYLINDER_22000_YPLUS_FRONT;
      N_Param++;
    }
   if (!strcmp(line, "CYLINDER_22000_YPLUS_BACK:"))
    {
	dat >> TDatabase::ParamDB->CYLINDER_22000_YPLUS_BACK;
      N_Param++;
    }

   if (!strcmp(line, "BULK_REACTION_DISC:"))
    {
      dat >> TDatabase::ParamDB->BULK_REACTION_DISC;
      N_Param++;
    }
    if (!strcmp(line, "BULK_SOLD_PARAMETER_TYPE:"))
    {
	dat >> TDatabase::ParamDB->BULK_SOLD_PARAMETER_TYPE;
      N_Param++;
    }
   if (!strcmp(line, "BULK_REACTION_MASS_LUMPING:"))
    {
      dat >> TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING;
      N_Param++;
    }
   if (!strcmp(line, "BULK_METHODS_OF_MOMENTS:"))
    {
	dat >> TDatabase::ParamDB->BULK_METHODS_OF_MOMENTS;
      N_Param++;
    }
   if (!strcmp(line, "BULK_MOM_DISC:"))
    {
	dat >> TDatabase::ParamDB->BULK_MOM_DISC;
      N_Param++;
    }
   if (!strcmp(line, "BULK_REACTION_C_CUT:"))
    {
      dat >> TDatabase::ParamDB->BULK_REACTION_C_CUT;
      N_Param++;
    }
   if (!strcmp(line, "BULK_PB_DISC:"))
    {
      dat >> TDatabase::ParamDB->BULK_PB_DISC;
      N_Param++;
    }
   if (!strcmp(line, "BULK_GROWTH_RATE:"))
    {
      dat >> TDatabase::ParamDB->BULK_GROWTH_RATE;
      N_Param++;
    }
   if (!strcmp(line, "BULK_PB_DISC_STAB:"))
    {
      dat >> TDatabase::ParamDB->BULK_PB_DISC_STAB;
      N_Param++;
    }
    if (!strcmp(line, "BULK_PB_DISC_FCT_GROUP:"))
    {
      dat >> TDatabase::ParamDB->BULK_PB_DISC_FCT_GROUP;
      N_Param++;
    }
    if (!strcmp(line, "N_CELL_LAYERS_PSD:"))
    {
      dat >> TDatabase::ParamDB->N_CELL_LAYERS_PSD;
      N_Param++;
    }
    if (!strcmp(line, "N_CELL_LAYERS_PSD_2:"))
    {
      dat >> TDatabase::ParamDB->N_CELL_LAYERS_PSD_2;
      N_Param++;
    }
    if (!strcmp(line, "OUTPUT_NODE_LAYER_PSD:"))
    {
      dat >> TDatabase::ParamDB->OUTPUT_NODE_LAYER_PSD;
      N_Param++;
    }
    if (!strcmp(line, "BULK_D_P_0:"))
    {
      dat >> TDatabase::ParamDB->BULK_D_P_0;
      N_Param++;
    }
    if (!strcmp(line, "BULK_D_P_MAX:"))
    {
      dat >> TDatabase::ParamDB->BULK_D_P_MAX;
      N_Param++;
    }
    if (!strcmp(line, "BULK_k_g:"))
    {
      dat >> TDatabase::ParamDB->BULK_k_g;
      N_Param++;
    }
    if (!strcmp(line, "BULK_C_g:"))
    {
      dat >> TDatabase::ParamDB->BULK_C_g;
      N_Param++;
    }
    if (!strcmp(line, "BULK_c_C_infty_sat:"))
    {
      dat >> TDatabase::ParamDB->BULK_c_C_infty_sat;
      N_Param++;
    }
    if (!strcmp(line, "BULK_l_infty:"))
    {
      dat >> TDatabase::ParamDB->BULK_l_infty;
      N_Param++;
    }
    if (!strcmp(line, "BULK_u_infty:"))
    {
      dat >> TDatabase::ParamDB->BULK_u_infty;
      N_Param++;
    }
    if (!strcmp(line, "BULK_k_nuc:"))
    {
      dat >> TDatabase::ParamDB->BULK_k_nuc;
      N_Param++;
    }
    if (!strcmp(line, "BULK_C_2:"))
    {
      dat >> TDatabase::ParamDB->BULK_C_2;
      N_Param++;
    }
    if (!strcmp(line, "WINDTUNNEL_dynamic_viscosity:"))
    {
      dat >>  TDatabase::ParamDB->WINDTUNNEL_dynamic_viscosity;
      N_Param++;
    }

    if (!strcmp(line, "WINDTUNNEL_CONFIGURATION:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_CONFIGURATION;
      N_Param++;
    }
    
    if (!strcmp(line, "WINDTUNNEL_STEADY:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_STEADY;
      N_Param++;
    }
    
    if (!strcmp(line, "WINDTUNNEL_INTERPOLATION:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_INTERPOLATION;
      N_Param++;
    }

    if (!strcmp(line, "WINDTUNNEL_BROWNIAN:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_BROWNIAN;
      N_Param++;
    }

    if (!strcmp(line, "WINDTUNNEL_POL_ORDER:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_POL_ORDER;
      N_Param++;
    }

   if (!strcmp(line, "WINDTUNNEL_SHEAR_FACTOR_TYPE:"))
   {
     dat >> TDatabase::ParamDB->WINDTUNNEL_SHEAR_FACTOR_TYPE;
     N_Param++;
   }

    if (!strcmp(line, "WINDTUNNEL_SHEAR_FACTOR:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_SHEAR_FACTOR;
      N_Param++;
    }

    if (!strcmp(line, "WINDTUNNEL_QUAD_METHOD:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_QUAD_METHOD;
      N_Param++;
    }
    if (!strcmp(line, "WINDTUNNEL_MEASURE_MASS:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_MEASURE_MASS;
      N_Param++;
    }
 
    if (!strcmp(line, "WINDTUNNEL_SHIFT:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_SHIFT;
      N_Param++;
    }
    if (!strcmp(line, "WINDTUNNEL_SUPERSAT:"))
    {
      dat >> TDatabase::ParamDB->WINDTUNNEL_SUPERSAT;
      N_Param++;
    }


    if (!strcmp(line, "SSMUM_MP_X:"))
    {
      dat >> TDatabase::ParamDB->SSMUM_MP_X;
      N_Param++;
    }

    if (!strcmp(line, "SSMUM_MP_Y:"))
    {
      dat >> TDatabase::ParamDB->SSMUM_MP_Y;
      N_Param++;
    }

    if (!strcmp(line, "SSMUM_OUTER_RADIUS:"))
    {
	dat >> TDatabase::ParamDB->SSMUM_OUTER_RADIUS;
      N_Param++;
    }

    if (!strcmp(line, "SSMUM_INNER_RADIUS:"))
    {
	dat >> TDatabase::ParamDB->SSMUM_INNER_RADIUS;
      N_Param++;
    }

    if (!strcmp(line, "SSMUM_ROT_PER_SECOND:"))
    {
      dat >> TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;
      N_Param++;
    }

    if (!strcmp(line, "SSMUM_INTERPOLATION:"))
    {
	dat >> TDatabase::ParamDB->SSMUM_INTERPOLATION;
      N_Param++;
    }

    if (!strcmp(line, "UREA_INFLOW_SCALE:"))
    {
	dat >> TDatabase::ParamDB->UREA_INFLOW_SCALE;
      N_Param++;
    }

    if (!strcmp(line, "UREA_PB_DISC:"))
    {
	dat >> TDatabase::ParamDB->UREA_PB_DISC;
      N_Param++;
    }
    
   if (!strcmp(line, "UREA_MODEL:"))
    {
        dat >> TDatabase::ParamDB->UREA_MODEL;
      N_Param++;
    }

   if (!strcmp(line, "UREA_CONC_MAXIT:"))
    {
        dat >> TDatabase::ParamDB->UREA_CONC_MAXIT;
      N_Param++;
    }
    
   if (!strcmp(line, "UREA_inflow_time:"))
    {
        dat >> TDatabase::ParamDB->UREA_inflow_time;
      N_Param++;
    }

       
    if (!strcmp(line, "UREA_u_infty:"))
    {
        dat >> TDatabase::ParamDB->UREA_u_infty;
      N_Param++;
    }    
  
    if (!strcmp(line, "UREA_D_P_MAX:"))
    {
        dat >> TDatabase::ParamDB->UREA_D_P_MAX;
      N_Param++;
    }    
  
  
  

   if (!strcmp(line, "UREA_AGGR_SPATIAL:"))
    {
        dat >> TDatabase::ParamDB->UREA_AGGR_SPATIAL;
      N_Param++;
    }

   if (!strcmp(line, "UREA_AGGR_BROWNIAN:"))
    {
        dat >> TDatabase::ParamDB->UREA_AGGR_BROWNIAN;
      N_Param++;
    }

   if (!strcmp(line, "UREA_AGGR_POL_ORDER:"))
    {
        dat >> TDatabase::ParamDB->UREA_AGGR_POL_ORDER;
      N_Param++;
    }

   if (!strcmp(line, "UREA_AGGR_SHEAR_FACTOR_TYPE:"))
    {
        dat >> TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE;
      N_Param++;
    }

   if (!strcmp(line, "UREA_AGGR_SHEAR_FACTOR:"))
    {
        dat >> TDatabase::ParamDB->UREA_AGGR_SHEAR_FACTOR;
      N_Param++;
    }
    if (!strcmp(line, "UREA_AGGR_BROWNIAN_TEMP:"))
    {
      dat >> TDatabase::ParamDB->UREA_AGGR_BROWNIAN_TEMP;
      N_Param++;
    }
    if (!strcmp(line, "UREA_AGGR_BROWNIAN_SCAL:"))
    {
      dat >> TDatabase::ParamDB->UREA_AGGR_BROWNIAN_SCAL;
      N_Param++;
    }
    if (!strcmp(line, "UREA_SOLD_PARAMETER_TYPE:"))
    {
        dat >> TDatabase::ParamDB->UREA_SOLD_PARAMETER_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "UREA_REACTION_DISC:"))
    {
      dat >> TDatabase::ParamDB->UREA_REACTION_DISC;
      N_Param++;
    }
    if (!strcmp(line, "PB_DISC_TYPE:"))
    {
      dat >> TDatabase::ParamDB->PB_DISC_TYPE;
      N_Param++;
    }
    if (!strcmp(line, "PB_TIME_DISC:"))
    {
      dat >> TDatabase::ParamDB->PB_TIME_DISC;
      N_Param++;
    }





//     if (!strcmp(line, "UREA_D_P_MAX:"))
//     {
//         dat >> TDatabase::ParamDB->UREA_D_P_MAX;
//       N_Param++;
//     }    
//   
//     if (!strcmp(line, "UREA_f_infty:"))
//     {
//         dat >> TDatabase::ParamDB->UREA_f_infty;
//       N_Param++;
//     }    
//   
  
    if (!strcmp(line, "MATLAB_MATRIX:"))
    {
      dat >> line;
      aux_char = new char[strlen(line) + 1];
      strcpy(aux_char, line);
      TDatabase::ParamDB->MATLAB_MATRIX = aux_char;
      N_Param++;
    }
    
    if (!strcmp(line, "WRITE_MATLAB:"))
    {
      dat >> TDatabase::ParamDB->WRITE_MATLAB;
      N_Param++;
    }

    if (!strcmp(line, "WRITE_MATLAB_MATRIX:"))
    {
      dat >> TDatabase::ParamDB->WRITE_MATLAB_MATRIX;
      N_Param++;
    }

    if (!strcmp(line, "NC_TYPE:"))
    {
      dat >> TDatabase::ParamDB->NC_TYPE;
      N_Param++;
    }
    
    if (!strcmp(line, "INPUT_QUAD_RULE:"))
    {
      dat >> TDatabase::ParamDB->INPUT_QUAD_RULE;
      N_Param++;
    }
    
    //======================================================================
    /** parameters for Stokes--Darcy (StoDa) coupling */
    //======================================================================
    if (!strcmp(line, "StoDa_interfaceType:"))
    { //Beavers-Joseph-Saffman or u.t=0
      dat >> TDatabase::ParamDB->StoDa_interfaceType;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_alpha:"))
    { // from Beavers-Joseph-Saffman condition on interface
      dat >> TDatabase::ParamDB->StoDa_alpha;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_problemType:"))
    { // Neumann--Neumann, Robin--Robin, ...
      dat >> TDatabase::ParamDB->StoDa_problemType;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_updatingStrategy:"))
    { // update of the etas
      dat >> TDatabase::ParamDB->StoDa_updatingStrategy;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_theta_f:"))
    { //damping in Stokes (flow) part
      dat >> TDatabase::ParamDB->StoDa_theta_f;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_theta_p:"))
    { //damping in Darcy (porous) part
      dat >> TDatabase::ParamDB->StoDa_theta_p;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_gamma_f:"))
    { // parameter for Robin condition on interface
      dat >> TDatabase::ParamDB->StoDa_gamma_f;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_gamma_p:"))
    { // parameter for Robin condition on interface
      dat >> TDatabase::ParamDB->StoDa_gamma_p;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_weakGamma:"))
    { // parameter for Robin condition on interface
      dat >> TDatabase::ParamDB->StoDa_weakGamma;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_solutionStrategy:"))
    { // only iterative (0) or iterative and one big matrix (1)
      dat >> TDatabase::ParamDB->StoDa_solutionStrategy;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_algorithm:"))
    { // Gauss--Seidel, Jacobi, ...
      dat >> TDatabase::ParamDB->StoDa_algorithm;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_StokesFirst:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_StokesFirst;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_nIterations:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_nIterations;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_relDiff_interfaceError:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_relDiff_interfaceError;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_relDiff_factor1:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_relDiff_factor1;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_relDiff_factor2:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_relDiff_factor2;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_relDiff_factor3:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_relDiff_factor3;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_relDiff_solution:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_relDiff_solution;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_bigResidual:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_bigResidual;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_periodicBoundary:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_periodicBoundary;
      N_Param++;
    }
    if (!strcmp(line, "StoDa_periodicBoundaryPressureDrop:"))
    { // for Gauss--Seidel type method.
      dat >> TDatabase::ParamDB->StoDa_periodicBoundaryPressureDrop;
      N_Param++;
    }

    

    if (!strcmp(line, "Par_P0:"))
    {
      dat >> TDatabase::ParamDB->Par_P0;
      N_Param++;
    }

    if (!strcmp(line, "Par_P1:"))
    {
      dat >> TDatabase::ParamDB->Par_P1;
      N_Param++;
    }

    if (!strcmp(line, "Par_P2:"))
    {
      dat >> TDatabase::ParamDB->Par_P2;
      N_Param++;
    }

    if (!strcmp(line, "Par_P3:"))
    {
      dat >> TDatabase::ParamDB->Par_P3;
      N_Param++;
    }

    if (!strcmp(line, "Par_P4:"))
    {
      dat >> TDatabase::ParamDB->Par_P4;
      N_Param++;
    }

    if (!strcmp(line, "Par_P5:"))
    {
      dat >> TDatabase::ParamDB->Par_P5;
      N_Param++;
    }

    if (!strcmp(line, "Par_P6:"))
    {
      dat >> TDatabase::ParamDB->Par_P6;
      N_Param++;
    }

    if (!strcmp(line, "Par_P7:"))
    {
      dat >> TDatabase::ParamDB->Par_P7;
      N_Param++;
    }

    if (!strcmp(line, "Par_P8:"))
    {
      dat >> TDatabase::ParamDB->Par_P8;
      N_Param++;
    }

    if (!strcmp(line, "Par_P9:"))
    {
      dat >> TDatabase::ParamDB->Par_P9;
      N_Param++;
    }

    if (!strcmp(line, "Par_P10:"))
    {
      dat >> TDatabase::ParamDB->Par_P10;
      N_Param++;
    }

    if (!strcmp(line, "Par_P11:"))
    {
      dat >> TDatabase::ParamDB->Par_P11;
      N_Param++;
    }

    if (!strcmp(line, "Par_P12:"))
    {
      dat >> TDatabase::ParamDB->Par_P12;
      N_Param++;
    }

    if (!strcmp(line, "Par_P13:"))
    {
      dat >> TDatabase::ParamDB->Par_P13;
      N_Param++;
    }

    if (!strcmp(line, "Par_P14:"))
    {
      dat >> TDatabase::ParamDB->Par_P14;
      N_Param++;
    }

    if (!strcmp(line, "Par_P15:"))
    {
      dat >> TDatabase::ParamDB->Par_P15;
      N_Param++;
    }

    if (!strcmp(line, "Par_P16:"))
    {
      dat >> TDatabase::ParamDB->Par_P16;
      N_Param++;
    }

    if (!strcmp(line, "Par_P17:"))
    {
      dat >> TDatabase::ParamDB->Par_P17;
      N_Param++;
    }

    if (!strcmp(line, "Par_P18:"))
    {
      dat >> TDatabase::ParamDB->Par_P18;
      N_Param++;
    }

    if (!strcmp(line, "Par_P19:"))
    {
      dat >> TDatabase::ParamDB->Par_P19;
      N_Param++;
    }

    if (!strcmp(line, "Par_P20:"))
    {
      dat >> TDatabase::ParamDB->Par_P20;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P0:"))
    {
      dat >> TDatabase::ParamDB->PBE_P0;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P1:"))
    {
      dat >> TDatabase::ParamDB->PBE_P1;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P2:"))
    {
      dat >> TDatabase::ParamDB->PBE_P2;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P3:"))
    {
      dat >> TDatabase::ParamDB->PBE_P3;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P4:"))
    {
      dat >> TDatabase::ParamDB->PBE_P4;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P5:"))
    {
      dat >> TDatabase::ParamDB->PBE_P5;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P6:"))
    {
      dat >> TDatabase::ParamDB->PBE_P6;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P7:"))
    {
      dat >> TDatabase::ParamDB->PBE_P7;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P8:"))
    {
      dat >> TDatabase::ParamDB->PBE_P8;
      N_Param++;
    }

    if (!strcmp(line, "PBE_P9:"))
    {
      dat >> TDatabase::ParamDB->PBE_P9;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P0:"))
    {
      dat >> TDatabase::ParamDB->DG_P0;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P1:"))
    {
      dat >> TDatabase::ParamDB->DG_P1;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P2:"))
    {
      dat >> TDatabase::ParamDB->DG_P2;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P3:"))
    {
      dat >> TDatabase::ParamDB->DG_P3;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P4:"))
    {
      dat >> TDatabase::ParamDB->DG_P4;
      N_Param++;
    }
    if (!strcmp(line, "DG_P5:"))
    {
      dat >> TDatabase::ParamDB->DG_P5;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P6:"))
    {
      dat >> TDatabase::ParamDB->DG_P6;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P7:"))
    {
      dat >> TDatabase::ParamDB->DG_P7;
      N_Param++;
    }
    
    if (!strcmp(line, "DG_P8:"))
    {
      dat >> TDatabase::ParamDB->DG_P8;
      N_Param++;
    } 
    
    if (!strcmp(line, "DG_P9:"))
    {
      dat >> TDatabase::ParamDB->DG_P9;
      N_Param++;
    }    
    
    if (!strcmp(line, "MOVING_BOUNDARY:"))
    {
      dat >> TDatabase::ParamDB->MOVING_BOUNDARY;
      N_Param++;
    }
    
    if (!strcmp(line, "DEPENDENT_BASIS:"))
    {
      dat >> TDatabase::ParamDB->DEPENDENT_BASIS;
      N_Param++;
    }
    
    if (!strcmp(line, "DEPENDENT_BASIS_Q1:"))
    {
      dat >> TDatabase::ParamDB->DEPENDENT_BASIS_Q1;
      N_Param++;
    }
    
    if (!strcmp(line, ":"))
    {
      dat >> TDatabase::ParamDB->DEPENDENT_BASIS_Q2;
      N_Param++;
    }
    
    if (!strcmp(line, "timeprofiling:"))
    {
      dat >> TDatabase::ParamDB->timeprofiling;
      N_Param++;
    }
    if (!strcmp(line, "MapperType:"))
    {
      dat >> TDatabase::ParamDB->MapperType;
      N_Param++;
    }
    if (!strcmp(line, "DSType:"))
    {
      dat >> TDatabase::ParamDB->DSType;
      N_Param++;
    }
        
    
 
   // read until end of line
    dat.getline (line, 99);
  }

  if (!flag[0])
  {
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR =
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR =
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
  }
  if (!flag[1])
    TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR =
      TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR;
  if (flag[2]!=1)
  {
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE =
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  }
  if (flag[2]!=2)
  {
    TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE =
      TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE;
  }
  if (!flag[3])
    TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE =
      TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SADDLE;

  if (TDatabase::ParamDB->START_RE_NR < 0)
    TDatabase::ParamDB->START_RE_NR = TDatabase::ParamDB->RE_NR;

  if (TDatabase::TimeDB->MIN_TIMESTEPLENGTH<0)
    TDatabase::TimeDB->MIN_TIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH/10.0;
  if (TDatabase::TimeDB->MAX_TIMESTEPLENGTH<0)
    TDatabase::TimeDB->MAX_TIMESTEPLENGTH = TDatabase::TimeDB->TIMESTEPLENGTH*10.0;

#ifdef _MPI
  if(rank==0)
#endif
  cout << "Parameter file version " << TDatabase::ParamDB->VERSION <<
          " read with " << N_Param << " parameters." << endl;

  dat.close();

  return 0;
}

int TDomain::ReadBdParam(char *ParamFile, int &Flag)
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  char line[100];
  int i, j, CurrentBdPart, N_BdComp, CompType, CompID = 0, N_Spls;
  std::ifstream dat(ParamFile);
#ifdef __2D__
  TBoundComp2D *BdComp;
#else
  TBoundComp3D *BdComp;
  TBoundComp2D *BdComp2D;
#endif

  Flag = 0;

  if (!dat)
  {
#ifdef _MPI
  if(rank==0)
#endif 
    cerr << "cannot open '" << ParamFile << "' for input" << endl;
    exit(-1);
  }

  // determine dimensions for creating arrays
  
  // get the number of boundaries (inner and outer) of the domain
  dat.getline (line, 99);
  dat >> N_BoundParts;
  dat.getline (line, 99);

  BdParts = new TBoundPart*[N_BoundParts];
  Interfaces = new int[N_BoundParts];
  N_BoundComps = 0;
  // StartBdCompID: index in the bdpart list where the new BdPart starts
  StartBdCompID = new int[N_BoundParts + 1];
  StartBdCompID[0] = 0; // set the first to 0

  for (i=0;i<N_BoundParts;i++)
  {
    dat.getline (line, 99); // IBCT
    dat >> CurrentBdPart;
    dat.getline (line, 99);
    Interfaces[i] = CurrentBdPart;
    CurrentBdPart = ABS(CurrentBdPart); // it can be negative (for orientation)
    if (i+1 != CurrentBdPart)
    {
#ifdef _MPI
  if(rank==0)
#endif
      cerr << "Error: different number of boundary part" << endl;
      cerr << "CurrentBdPart " << i << "  " << CurrentBdPart << endl;
      exit(-1);
    }

    // get number of components f the bdpart i
    dat.getline (line, 99);
    dat >> N_BdComp;
    dat.getline (line, 99);
    
    BdParts[i] = new TBoundPart(N_BdComp);
    N_BoundComps += N_BdComp;
    SetStartBdCompID(N_BoundComps, i+1); // set the start ID of the next bdcomp

    dat.getline (line, 99);
    for (j=0;j<N_BdComp;j++)
    {
      dat >> CompType >> N_Spls;  //ITYP NSPLINE NPAR
      dat.getline (line, 99);

#ifdef __2D__
      // 2D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5)
      switch (abs(CompType))
      {
        case 1: BdComp = new TBdLine(CompID++);
                break;
        case 2: BdComp = new TBdCircle(CompID++);
                break;
        case 3: BdComp = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        default:
#ifdef _MPI
  if(rank==0)
#endif
          OutPut("ReadParam.C: Boundary type not implemented" << endl);
          exit(-1);
      }
#else
      // 3D types: Line (1), Circle (2), Spline (3), Poygon (4), NonUnif Spline (5),
      //           Plane (10), Sphere (11)
      switch (abs(CompType))
      {
        case 1: BdComp2D = new TBdLine(CompID++);
                break;
        case 2: BdComp2D = new TBdCircle(CompID++);
                break;
        case 3: BdComp2D = new TBdSpline(CompID++, N_Spls);
                break;
        case 4: BdComp2D = new TBdPolygon(CompID++, N_Spls);
                break;
        case 5: BdComp2D = new TBdNonUniformSpline(CompID++, N_Spls);
                break;
        case 10: BdComp = new TBdPlane(CompID++);
                break;
        case 11: BdComp = new TBdSphere(CompID++);
                break;
        case 4711: BdComp = new TBdNoPRM(CompID++); // create grid without PRM file
	    break;
        default:
#ifdef _MPI
  if(rank==0)
#endif
          OutPut("ReadParam.C: Boundary type (3D) not implemented" << endl);
          exit(-1);
      }

      if(abs(CompType)<10)
      {
        BdComp = new TBdWall(CompID-1, BdComp2D);
        Flag = 1;
      }
#endif // 3D
      BdParts[i]->SetBdComp(j, BdComp);

      if(CompType<0)
      {
        BdComp->SetFreeBoundaryStatus(true);
        cout <<i<< " ReadBdParam : " << j << endl;
      }
    }
  }
    
  dat.getline (line, 99);
  for (i=0;i<N_BoundParts;i++)
  {
    N_BdComp = BdParts[i]->GetN_BdComps();
    for (j=0;j<N_BdComp;j++)
    {
      BdParts[i]->GetBdComp(j)->ReadIn(dat);
    }
  }

  // read HOLES (if any)
  dat.getline (line, 99);
  N_Holes = -12345;
  if (dat.eof())
    N_Holes = 0;
  else
    dat >> N_Holes;

  if(N_Holes == -12345)
    N_Holes = 0;

  dat.getline (line, 99);

  if (N_Holes)
  {
    // coordinates of a point in a hole
    PointInHole = new double[2*N_Holes];

    dat.getline (line, 99);
    for (i=0;i<N_Holes;i++)
    {
      dat >> PointInHole[2*i] >> PointInHole[2*i+1];
      dat.getline (line, 99);
    }
  }
  else
    PointInHole = NULL;

  dat.getline (line, 99);
  N_Regions = -12345;
  if (dat.eof())
    N_Regions = 0;
  else
    dat >> N_Regions;

  if(N_Regions == -12345)
    N_Regions = 0;
    
  dat.getline (line, 99);

  // read REGIONS (if any)
  if (N_Regions)
  {
    PointInRegion = new double[4*N_Regions];
  
    dat.getline (line, 99);
    for (i=0;i<N_Regions;i++)
    {
      dat >> PointInRegion[4*i] >> PointInRegion[4*i+1];
      PointInRegion[4*i+2] = i;
      PointInRegion[4*i+3] = 10000;
      dat.getline (line, 99);
    }
  }
  else
    PointInRegion = NULL;

  dat.close();
  return 0;
}

int TDomain::ReadMapFile(char *MapFile, TDatabase *Database)
{
#ifdef _MPI
  int rank, out_rank=int(TDatabase::ParamDB->Par_P0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  char line[100];
  int i, j, N_, N_Cells, ID, NMortarCell, LocEdge;
  int N_MortarRefs, CellID, N_MortarRefDesc = 0;
  int *MortarRefDesc;
  bool Exist;
  std::ifstream dat(MapFile);
  TBaseCell *CurrCell, *NeighbCell;
  TJoint *CurrJoint;

  if (!dat)
  {
#ifdef _MPI
  if(rank==0)
#endif
    cerr << "cannot open '" << MapFile << "' for input" << endl;
    exit(-1);
  }

  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);
  dat.getline (line, 99);

  // get number of elements
  dat >> N_Cells;
  dat.getline (line, 99);

  if (N_Cells != N_RootCells)
  {
#ifdef _MPI
  if(rank==0)
#endif
    cerr << "Error in ReadMapFile: wrong number of elements!!!" << endl;
    exit(-1);
  }

  dat.getline (line, 99);
  for (i=0;i<N_Cells;i++)
  {
    dat >> ID;
    ((TMacroCell *) CellTree[i])->SetSubGridID(ID);
    dat.getline (line, 99);
  }

#ifdef __MORTAR__

  MortarRefDesc=new int[2 * N_MORTARDESC]; 

  dat.getline (line, 99);
  dat.getline (line, 99);

  // get number of MortarBaseJoints
  dat >> N_MortarFaces;
  dat.getline (line, 99);

  dat.getline (line, 99);
  if (N_MortarFaces)
  {
    MortarFaces = new TMortarFace[N_MortarFaces];

    for (i=0;i<N_MortarFaces;i++)
    {
      dat >> NMortarCell >> LocEdge;
      dat.getline (line, 99);

      CurrCell = CellTree[NMortarCell];
      CurrJoint = CurrCell->GetJoint(LocEdge);
      NeighbCell = CurrJoint->GetNeighbour(CurrCell);

      N_ = NeighbCell->GetRefDesc()->GetN_OrigEdges();
      for (j=0;j<N_;j++)
        if (NeighbCell->GetJoint(j) == CurrJoint) break;

      MortarFaces[i].Cell = CurrCell;
      MortarFaces[i].LocFaceNumber[0] = LocEdge;
      MortarFaces[i].LocFaceNumber[1] = j;
    }
  }

  for (i=0;i<N_MortarFaces;i++)
  {
    j = MortarFaces[i].LocFaceNumber[0];
    CurrCell = MortarFaces[i].Cell;
    CurrJoint = CurrCell->GetJoint(j);
    NeighbCell = CurrJoint->GetNeighbour(CurrCell);

    delete CurrJoint;

    CurrJoint = new TMortarBaseJoint(CurrCell, NeighbCell);
    CurrCell->SetJoint(j, CurrJoint);
    NeighbCell->SetJoint(MortarFaces[i].LocFaceNumber[1], CurrJoint);
  }

  dat.getline (line, 99);
  dat.getline (line, 99);

  // get number of mortar refinements
  dat >> N_MortarRefs;
  dat.getline (line, 99);

  if (N_MortarRefs)
  {
    dat.getline (line, 99);
    for (i=0;i<N_MortarRefs;i++)
    {
      dat >> CellID >> LocEdge >> N_Cells;
      dat.getline (line, 99);

      Exist = false;
      for (j=0;j<N_MortarRefDesc;j++)
        if (MortarRefDesc[2 * j] == LocEdge)
          if (MortarRefDesc[2 * j + 1] == N_Cells)
          {
            Exist = true;
            break;
          }

      if (!Exist)
      {
        if (N_MortarRefDesc == N_MORTARDESC)
        {
#ifdef _MPI
  if(rank==0)
#endif
          cerr << "Error in ReadMapFile: not enough MortarRefDescs" << endl;
          exit(-1);
        }

        if (LocEdge)
          Database->AddMortar1(2 * N_MortarRefDesc, N_Cells);
        else
          Database->AddMortar0(2 * N_MortarRefDesc, N_Cells);

        MortarRefDesc[2 * N_MortarRefDesc] = LocEdge;
        MortarRefDesc[2 * N_MortarRefDesc++ + 1] = N_Cells;
      }

      CellTree[CellID]->SetRefDesc(TDatabase::RefDescDB[N_SHAPES +
                                   Mortar + 2*j]);
    }

    Refine();
  }
#endif

  dat.close();

  return 0;
}
