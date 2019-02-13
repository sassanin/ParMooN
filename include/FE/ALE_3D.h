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
   

#include <stdlib.h>
#include <math.h>
#include <Database.h>


#include <Convolution.h>
#include <MooNMD_Io.h>
#include <ConvDiff.h>

#include <TNSE3D_Routines.h>
#include <MainUtilities.h>

void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF);
// ======================================================================
// declaration for InitializeDiscreteFormGrid()
// ======================================================================

static int GridN_Terms_3D = 3;
static MultiIndex3D GridDerivatives_3D[3] = { D100, D010, D001 };
static int GridSpaceNumbers_3D[3] = { 0, 0 , 0 };
static int GridN_Matrices_3D = 9;
static int GridRowSpace_3D[9] = { 0, 0 , 0,  0, 0 , 0 ,  0, 0 , 0  };
static int GridColumnSpace_3D[9] = { 0, 0 , 0,  0, 0 , 0 ,  0, 0 , 0  };
static int GridN_Rhs_3D = 0;
static int *GridRhsSpace_3D = NULL;

void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs);
