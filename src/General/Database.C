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
// @(#)Database.C        1.37 06/27/00
// 
// Class:       TDatabase
// Purpose:     database of needed refinement, mapping and
//              shape descriptors
//              parameter database
//
// Author:      Volker Behns  29.07.97
//
// =======================================================================
#if defined(_MPI) || defined(_SMPI)
#  include "mpi.h"
#endif

#include <Database.h>
#include <Constants.h>
#include <Mapper.h>
#include <Line.h>
#include <Triangle.h>
#include <Quadrangle.h>
#include <Parallelogram.h>
#include <Rectangle.h>
#include <RefNoRef.h>
#include <RefLineDesc.h>
#include <RefQuadRegDesc.h>
#include <RefQuadBis0Desc.h>
#include <RefQuadBis1Desc.h>
#include <RefQuad1Conf0Desc.h>
#include <RefQuad1Conf1Desc.h>
#include <RefQuad1Conf2Desc.h>
#include <RefQuad1Conf3Desc.h>
#include <RefQuad2Conf0Desc.h>
#include <RefQuad2Conf1Desc.h>
#include <RefQuad2Conf2Desc.h>
#include <RefQuad2Conf3Desc.h>
#include <RefQuadToTri0Desc.h>
#include <RefQuadToTri1Desc.h>
#include <RefTriRegDesc.h>
#include <RefTriBis0Desc.h>
#include <RefTriBis1Desc.h>
#include <RefTriBis2Desc.h>
#include <RefTriBis01Desc.h>
#include <RefTriBis02Desc.h>
#include <RefTriBis10Desc.h>
#include <RefTriBis12Desc.h>
#include <RefTriBis20Desc.h>
#include <RefTriBis21Desc.h>
#include <It_Between.h>
#include <It_EQ.h>
#include <It_Finest.h>
#include <It_EQLevel.h>
#include <It_LELevel.h>
#include <It_OCAF.h>

#include <MooNMD_Io.h>

#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
extern "C" {
#include <amg_solve_main.h>
}
#ifdef __MORTAR__
  #include <It_Mortar.h>
  #include <RefMortar0Desc.h>
  #include <RefMortar1Desc.h>
  #include <RefMortarLineDesc.h>
#endif

#ifdef __3D__
  #include <Tetrahedron.h>
  #include <Hexahedron.h>
  #include <Brick.h>
  #include <RefTetraRegDesc.h>
  #include <RefTetraReg0Desc.h>
  #include <RefTetraReg1Desc.h>
  #include <RefTetraReg2Desc.h>
  #include <RefTetraBis0Desc.h>
  #include <RefTetraBis1Desc.h>
  #include <RefTetraBis2Desc.h>
  #include <RefTetraBis3Desc.h>
  #include <RefTetraBis4Desc.h>
  #include <RefTetraBis5Desc.h>
  #include <RefTetraBis01Desc.h>
  #include <RefTetraBis02Desc.h>
  #include <RefTetraBis03Desc.h>
  #include <RefTetraBis04Desc.h>
  #include <RefTetraBis05Desc.h>
  #include <RefTetraBis10Desc.h>
  #include <RefTetraBis12Desc.h>
  #include <RefTetraBis13Desc.h>
  #include <RefTetraBis14Desc.h>
  #include <RefTetraBis15Desc.h>
  #include <RefTetraBis20Desc.h>
  #include <RefTetraBis21Desc.h>
  #include <RefTetraBis23Desc.h>
  #include <RefTetraBis24Desc.h>
  #include <RefTetraBis25Desc.h>
  #include <RefTetraBis30Desc.h>
  #include <RefTetraBis32Desc.h>
  #include <RefTetraBis34Desc.h>
  #include <RefTetraBis35Desc.h>
  #include <RefTetraBis40Desc.h>
  #include <RefTetraBis41Desc.h>
  #include <RefTetraBis43Desc.h>
  #include <RefTetraBis45Desc.h>
  #include <RefTetraBis51Desc.h>
  #include <RefTetraBis52Desc.h>
  #include <RefTetraBis53Desc.h>
  #include <RefTetraBis54Desc.h>
  #include <RefTetraQuad0Desc.h>
  #include <RefTetraQuad1Desc.h>
  #include <RefTetraQuad2Desc.h>
  #include <RefTetraQuad3Desc.h>
  #include <RefHexaRegDesc.h>
#endif

// Constructors
TDatabase::TDatabase()
{
  // allocate databases
  ShapeDB = new TShapeDesc*[N_SHAPES];
  RefDescDB = new TRefDesc*[N_SHAPES + N_REFDESC + 2*N_MORTARDESC];
  MapperDB = new TMapper*[N_MAPPER];
  IteratorDB = new TIterator*[N_ITERATORS];
  ParamDB = new TParamDB;
  TimeDB = new TTimeDB;

  // initialize shape descriptors
  ShapeDB[S_Line] = new TLine();
  RefDescDB[S_Line] = new TRefNoRef(ShapeDB[S_Line]);

  ShapeDB[Triangle] = new TTriangle();
  RefDescDB[Triangle] = new TRefNoRef(ShapeDB[Triangle]);

  ShapeDB[Quadrangle] = new TQuadrangle();
  RefDescDB[Quadrangle] = new TRefNoRef(ShapeDB[Quadrangle]);

  ShapeDB[Parallelogram] = new TParallelogram();
  RefDescDB[Parallelogram] = new TRefNoRef(ShapeDB[Parallelogram]);

  ShapeDB[Rectangle] = new TRectangle();
  RefDescDB[Rectangle] = new TRefNoRef(ShapeDB[Rectangle]);

  #ifdef __3D__
    ShapeDB[Tetrahedron] = new TTetrahedron();
    RefDescDB[Tetrahedron] = new TRefNoRef(ShapeDB[Tetrahedron]);

    ShapeDB[Hexahedron] = new THexahedron();
    RefDescDB[Hexahedron] = new TRefNoRef(ShapeDB[Hexahedron]);

    ShapeDB[Brick] = new TBrick();
    RefDescDB[Brick] = new TRefNoRef(ShapeDB[Brick]);
  #endif

  // initialize refinement descriptors
  RefDescDB[N_SHAPES + LineReg] = new TRefLineDesc(ShapeDB[S_Line]);
  RefDescDB[N_SHAPES + TriReg]  = new TRefTriRegDesc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis0] = new TRefTriBis0Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis1] = new TRefTriBis1Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis2] = new TRefTriBis2Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis01]= new TRefTriBis01Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis02]= new TRefTriBis02Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis10]= new TRefTriBis10Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis12]= new TRefTriBis12Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis20]= new TRefTriBis20Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis21]= new TRefTriBis21Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + QuadReg] = new TRefQuadRegDesc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + ParallReg] = new TRefQuadRegDesc(ShapeDB[Parallelogram]);
  RefDescDB[N_SHAPES + RectReg] = new TRefQuadRegDesc(ShapeDB[Rectangle]);
  RefDescDB[N_SHAPES + QuadBis0] = new TRefQuadBis0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadBis1] = new TRefQuadBis1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf0] = new TRefQuad1Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf1] = new TRefQuad1Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf2] = new TRefQuad1Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf3] = new TRefQuad1Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf0] = new TRefQuad2Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf1] = new TRefQuad2Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf2] = new TRefQuad2Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf3] = new TRefQuad2Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri0] = new
      TRefQuadToTri0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri1] = new
      TRefQuadToTri1Desc(ShapeDB[Quadrangle]);

  #ifdef __3D__
    RefDescDB[N_SHAPES + TetraReg] =
         new TRefTetraRegDesc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg0] =
         new TRefTetraReg0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg1] =
         new TRefTetraReg1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg2] =
         new TRefTetraReg2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis0] = new TRefTetraBis0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis1] = new TRefTetraBis1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis2] = new TRefTetraBis2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis3] = new TRefTetraBis3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis4] = new TRefTetraBis4Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis5] = new TRefTetraBis5Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis01] = new TRefTetraBis01Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis02] = new TRefTetraBis02Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis03] = new TRefTetraBis03Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis04] = new TRefTetraBis04Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis05] = new TRefTetraBis05Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis10] = new TRefTetraBis10Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis12] = new TRefTetraBis12Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis13] = new TRefTetraBis13Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis14] = new TRefTetraBis14Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis15] = new TRefTetraBis15Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis20] = new TRefTetraBis20Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis21] = new TRefTetraBis21Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis23] = new TRefTetraBis23Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis24] = new TRefTetraBis24Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis25] = new TRefTetraBis25Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis30] = new TRefTetraBis30Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis32] = new TRefTetraBis32Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis34] = new TRefTetraBis34Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis35] = new TRefTetraBis35Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis40] = new TRefTetraBis40Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis41] = new TRefTetraBis41Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis43] = new TRefTetraBis43Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis45] = new TRefTetraBis45Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis51] = new TRefTetraBis51Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis52] = new TRefTetraBis52Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis53] = new TRefTetraBis53Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis54] = new TRefTetraBis54Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad0] = new TRefTetraQuad0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad1] = new TRefTetraQuad1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad2] = new TRefTetraQuad2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad3] = new TRefTetraQuad3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + HexaReg]  =
         new TRefHexaRegDesc(ShapeDB[Hexahedron]);
    RefDescDB[N_SHAPES + BrickReg]  =
         new TRefHexaRegDesc(ShapeDB[Brick]);
  #endif

  #ifdef __3D__
    //initialize mapper
    MapperDB[MapTriReg0] = new TMapper(MapTriReg0);
    MapperDB[MapTriReg1] = new TMapper(MapTriReg1);
    MapperDB[MapTriReg2] = new TMapper(MapTriReg2);

    MapperDB[MapTriBis00] = new TMapper(MapTriBis00);
    MapperDB[MapTriBis01] = new TMapper(MapTriBis01);
    MapperDB[MapTriBis02] = new TMapper(MapTriBis02);
    MapperDB[MapTriBis10] = new TMapper(MapTriBis10);
    MapperDB[MapTriBis11] = new TMapper(MapTriBis11);
    MapperDB[MapTriBis12] = new TMapper(MapTriBis12);
    MapperDB[MapTriBis20] = new TMapper(MapTriBis20);
    MapperDB[MapTriBis21] = new TMapper(MapTriBis21);
    MapperDB[MapTriBis22] = new TMapper(MapTriBis22);
    MapperDB[MapTriBis010] = new TMapper(MapTriBis010);
    MapperDB[MapTriBis011] = new TMapper(MapTriBis011);
    MapperDB[MapTriBis012] = new TMapper(MapTriBis012);
    MapperDB[MapTriBis020] = new TMapper(MapTriBis020);
    MapperDB[MapTriBis021] = new TMapper(MapTriBis021);
    MapperDB[MapTriBis022] = new TMapper(MapTriBis022);
    MapperDB[MapTriBis100] = new TMapper(MapTriBis100);
    MapperDB[MapTriBis101] = new TMapper(MapTriBis101);
    MapperDB[MapTriBis102] = new TMapper(MapTriBis102);
    MapperDB[MapTriBis120] = new TMapper(MapTriBis120);
    MapperDB[MapTriBis121] = new TMapper(MapTriBis121);
    MapperDB[MapTriBis122] = new TMapper(MapTriBis122);
    MapperDB[MapTriBis200] = new TMapper(MapTriBis200);
    MapperDB[MapTriBis201] = new TMapper(MapTriBis201);
    MapperDB[MapTriBis202] = new TMapper(MapTriBis202);
    MapperDB[MapTriBis210] = new TMapper(MapTriBis210);
    MapperDB[MapTriBis211] = new TMapper(MapTriBis211);
    MapperDB[MapTriBis212] = new TMapper(MapTriBis212);

    MapperDB[MapQuadReg0] = new TMapper(MapQuadReg0);
    MapperDB[MapQuadReg1] = new TMapper(MapQuadReg1);
    MapperDB[MapQuadReg2] = new TMapper(MapQuadReg2);
    MapperDB[MapQuadReg3] = new TMapper(MapQuadReg3);
  #endif

  // initialize iterators
  IteratorDB[It_EQ] = new TIt_EQ();
  IteratorDB[It_LE] = new TIt_LE();
  IteratorDB[It_Finest] = new TIt_Finest();
  IteratorDB[It_EQLevel] = new TIt_EQLevel();
  IteratorDB[It_LELevel] = new TIt_LELevel();
  IteratorDB[It_Between] = new TIt_Between();
  IteratorDB[It_OCAF] = new TIt_OCAF();

  #ifdef __MORTAR__
    IteratorDB[It_Mortar1] = new TIt_Mortar();
    IteratorDB[It_Mortar2] = new TIt_Mortar();
  #endif
}

TShapeDesc **TDatabase::ShapeDB = NULL;
TRefDesc   **TDatabase::RefDescDB = NULL;
TMapper    **TDatabase::MapperDB = NULL;
TIterator  **TDatabase::IteratorDB = NULL;
TParamDB   *TDatabase::ParamDB = NULL;
TTimeDB    *TDatabase::TimeDB = NULL;

// Methods

#ifdef __MORTAR__

void TDatabase::AddMortar0(int Mortar_Ni, int N)
{
  RefDescDB[N_SHAPES + Mortar + Mortar_Ni] = new
               TRefMortar0Desc(ShapeDB[Quadrangle], Mortar_Ni, N);
  RefDescDB[N_SHAPES + MortarLine + Mortar_Ni] = new
               TRefMortarLineDesc(ShapeDB[S_Line], N);
}

void TDatabase::AddMortar1(int Mortar_Ni, int N)
{
  RefDescDB[N_SHAPES + Mortar + Mortar_Ni] = new
               TRefMortar1Desc(ShapeDB[Quadrangle], Mortar_Ni, N);
  RefDescDB[N_SHAPES + MortarLine + Mortar_Ni] = new
               TRefMortarLineDesc(ShapeDB[S_Line], N);
}
#endif

void TDatabase::SetDefaultParameters()
{
  char *tmp;
  ParamDB->VERSION = 1;

  tmp = new char[12];
  strcpy(tmp,"NO_GEO_FILE");
  ParamDB->GEOFILE=tmp;
  tmp = new char[17];
  strcpy(tmp,"NO_GEO_FILE_INTL");
  ParamDB->GEOFILE_INTL=tmp;
  
  tmp = new char[12];
  strcpy(tmp,"NO_BND_FILE");
  ParamDB->BNDFILE=tmp;
    
  tmp = new char[17];
  strcpy(tmp,"NO_BND_FILE_INTL");
  ParamDB->BNDFILE_INTL=tmp;
  
  tmp = new char[12];
  strcpy(tmp,"NO_MAP_FILE");
  ParamDB->MAPFILE=tmp;
  tmp = new char[25];
  strcpy(tmp,"MooN_MD_default_outfile");
  ParamDB->OUTFILE=tmp;
  
  ParamDB->PROBLEM_TYPE = 0;
  ParamDB->EXAMPLE = -1; // has to be set to some number >=0
  
  ParamDB->timeprofiling = 0; //time profiling
  ParamDB->MapperType = 1;
  ParamDB->DSType = 1;		//Parallel Direct Solver Type

  ParamDB->MESHGEN_ALLOW_EDGE_REF=0;
  ParamDB->MESHGEN_REF_QUALITY=30;
 
  ParamDB->RE_NR=1.0;
  ParamDB->RA_NR=1.0;
  ParamDB->ROSSBY_NR=0.0;
  ParamDB->START_RE_NR= -4711;
  ParamDB->RE_NR_INCREMENT=1.0;
  ParamDB->FLOW_PROBLEM_TYPE = 0;
  ParamDB->OSEEN_ZERO_ORDER_COEFF = 0.0;

  ParamDB->FR_NR=1.0;
  ParamDB->WB_NR= 1.0;
  ParamDB->PR_NR=1.0;
  ParamDB->PE_NR=1.0;  
  ParamDB->BI_NR = 0;
  ParamDB->WEI_NR=1.0;
  ParamDB->Axial3D = 0;
  ParamDB->Axial3DAxis = 0;  
  
  // ------------------time parameters
  ParamDB->time_system_assemble =0.0;
  ParamDB->time_solve=0.0;
  ParamDB->time_communication=0.0;
  ParamDB->time_GMRES=0.0;
  ParamDB->time_MG=0.0;
  ParamDB->time_projection=0.0;
  ParamDB->time_restriction=0.0;
  ParamDB->time_vanka = 0.0;
  ParamDB->time_vanka_solve = 0.0;
  // ------------------ end of time parameters


  ParamDB->ANSATZ_ORDER = 2;
  ParamDB->TEST_ORDER = 2;

  ParamDB->VELOCITY_SPACE = 22;
  ParamDB->PRESSURE_SPACE = -4711;
  ParamDB->PRESSURE_SEPARATION = 0;
  ParamDB->OMPNUMTHREADS=1;

  ParamDB->LEVELS = 1000;
  ParamDB->UNIFORM_STEPS = 1000;
  ParamDB->DRIFT_X = 0;
  ParamDB->DRIFT_Y = 0;
  ParamDB->DRIFT_Z = 0.41;
  ParamDB->NONLINEARIT_TYPE_NEWTON = 0;

 
  ParamDB->GRID_TYPE = 0;
  ParamDB->GRID_TYPE_1 = 0;
  ParamDB->GRID_TYPE_2 = 0;
  ParamDB->CHANNEL_GRID_STRETCH = 2.75;
  ParamDB->ADAPTIVE_REFINEMENT_CRITERION = 3;
  ParamDB->ERROR_CONTROL = 1;
  ParamDB->REFINE_STRATEGY = 0;
  ParamDB->MAX_CELL_LEVEL = 1000;
  ParamDB->REFTOL = 0.5;
  ParamDB->COARSETOL = 0.0;  
  ParamDB->MIN_FRACTION_TO_CHANGE = 0.1;
  ParamDB->DECREASE_REFTOL_FACTOR = 0.8;
  ParamDB->INCREASE_COARSETOL_FACTOR = 1.1;
  ParamDB->FRACTION_OF_ERROR = 0.25;
  ParamDB->CONVERT_QUAD_TO_TRI = 0;
  ParamDB->N_CELL_LAYERS = 1;
  
  ParamDB->DISCTYPE = 1;
  ParamDB->INTL_DISCTYPE = 1;
  ParamDB->UPWIND_ORDER = 1;
  ParamDB->UPWIND_FLUX_DAMP = 1;
  ParamDB->UPWIND_APPLICATION = 0;
  ParamDB->SHISHKIN_MESH = 0;
  ParamDB->SHISHKIN_DIAM = 1.0;
  ParamDB->NSTYPE = 1;
  ParamDB->DARCYTYPE = 1;
  ParamDB->SIGMA_PERM = 1;
  ParamDB->LAPLACETYPE = 0;
  ParamDB->TENSOR_TYPE = 0;
  ParamDB->TWO_PHASE_FLOW = 0;
  ParamDB->FREE_SURFACE_FLOW = 0;
  ParamDB->PHASE1_TYPE = 1;
  ParamDB->PHASE2_TYPE = 1;
  ParamDB->USE_ISOPARAMETRIC = 1;
  ParamDB->VMM_COARSE_LEVEL = 4711;
  ParamDB->VMM_COARSE_SPACE_ORDER = 1;
  ParamDB->RFB_SUBMESH_LAYERS = 3;
  ParamDB->DEFECT_CORRECTION_TYPE = 0;
  ParamDB->CELL_MEASURE = 0;

  ParamDB->FACE_SIGMA = 1;
  ParamDB->WEAK_BC_SIGMA = 1;
  ParamDB->WEAK_BC = 0;
  ParamDB->TAU=1.0;
  ParamDB->TAU2=1.0;
  ParamDB->TAU3=1.0;
 
  ParamDB->DELTA0=1.0;
  ParamDB->SDFEM_POWER0=1.0;
  ParamDB->DELTA1=1.0;
  ParamDB->SDFEM_TYPE=2;
  ParamDB->SDFEM_NORM_B=0; // l_infty
  ParamDB->CIP_TYPE=0;
  ParamDB->ADJOINT_FACTOR_4_OMEGA_EQ_0 = 10.0;
  ParamDB->DELTA2=1.0;
  
  ParamDB->FILTER_WIDTH_CONSTANT = 2;
  ParamDB->FILTER_WIDTH_POWER = 1;
  ParamDB->GAUSSIAN_GAMMA = 6;
  ParamDB->CONVOLUTE_SOLUTION = 0;
  
  ParamDB->TURBULENT_VISCOSITY_TYPE = 1;
  ParamDB->TURBULENT_VISCOSITY_TENSOR = 0;
  ParamDB->TURBULENT_VISCOSITY_CONSTANT = 0.01;
  ParamDB->TURBULENT_VISCOSITY_POWER = 1;
  ParamDB->TURBULENT_VISCOSITY_SIGMA = 6;
  ParamDB->TURBULENT_MOD_TYPE = 1;
  
  ParamDB->viscosity_max = -1;
  ParamDB->viscosity_min = 100;
  

  ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT = 1; // parameters for VMS
  ParamDB->ARTIFICIAL_VISCOSITY_POWER = 1;

  ParamDB->FRICTION_CONSTANT = 0.0;      // free slip 
  ParamDB->FRICTION_POWER = 0.0;         // free slip
  ParamDB->FRICTION_TYPE = 0;            // friction type
  ParamDB->FRICTION_U0 = 1.0;            // U_0
  ParamDB->PENETRATION_CONSTANT = 1e12;  // no penetration
  ParamDB->PENETRATION_POWER = -2;        // no penetration

  ParamDB->DIV_DIV_STAB_TYPE = 0;        // stabilization for div-div term 
  ParamDB->DIV_DIV_STAB_C1 = 2;
  ParamDB->DIV_DIV_STAB_C2 = 1;

  ParamDB->NSE_NONLINEAR_FORM = 0;       // skew symmetric convective term in NSE

  ParamDB->LP_FULL_GRADIENT = 1;
  ParamDB->LP_STREAMLINE = 0;
  ParamDB->LP_DIVERGENCE = 0;
  ParamDB->LP_PRESSURE = 0;
  ParamDB->LP_COEFF_TYPE = 0;

  ParamDB->LP_FULL_GRADIENT_COEFF = 1.0;
  ParamDB->LP_STREAMLINE_COEFF= 1.0;
  ParamDB->LP_DIVERGENCE_COEFF = 1.0;
  ParamDB->LP_PRESSURE_COEFF = 1.0;

  ParamDB->LP_FULL_GRADIENT_EXPONENT = 1.0;
  ParamDB->LP_STREAMLINE_EXPONENT = 1.0;
  ParamDB->LP_DIVERGENCE_EXPONENT = 1.0;
  ParamDB->LP_PRESSURE_EXPONENT = 1.0;

  ParamDB->LP_ORDER_DIFFERENCE = 1;
  ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = -123;
  ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = -123;
  ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = -123;
  ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = -123;
  
  ParamDB->LP_CROSSWIND = 0;
  ParamDB->LP_CROSSWIND_COEFF_TYPE = 1;
  ParamDB->LP_CROSSWIND_COEFF = 1.0;
  ParamDB->LP_CROSSWIND_EXPONENT = 1.0;

  //======================================================================
  /** parameter for a posteriori parameter computation with adjoint problem */
  //======================================================================
  ParamDB->SOLVE_ADJOINT_PROBLEM = FALSE; 
  ParamDB->SOLD_ADJOINT = 0;
  ParamDB->N_STAGES_ADJOINT = 1;
  ParamDB->SC_NONLIN_ITE_ADJOINT = 1000;
  ParamDB->OPTIMIZATION_ITE_TYPE_ADJOINT = 0;
  ParamDB->BFGS_VECTORS_ADJOINT = 25;
  ParamDB->RELATIVE_DECREASE_ADJOINT = 1e-4;
  ParamDB->PENALTY_ADJOINT = 0;
  ParamDB->PENALTY_VALUE_AT_ZERO_ADJOINT = 1e6;
  ParamDB->PENALTY_SMALLEST_PARAM_FAC_ADJOINT = 1e-2;
  ParamDB->PENALTY_LARGEST_PARAM_FAC_ADJOINT = 1e1;
  ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT = 1;
  ParamDB->WEIGHT_GRADIENT_L1_ADJOINT = 0;
  ParamDB->WEIGHT_GRADIENT_L2_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT = 0;
  ParamDB->REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT = 1.0;
  ParamDB->MIN_VAL_ADJOINT = -1e20;
  ParamDB->MAX_VAL_ADJOINT = 1e20;
  ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT = 1;
  ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT = 1;
  ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT = 1;
  ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT = 1;
  ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT = 0;
  ParamDB->RESIDUAL_LP_ADJOINT = 2.0;

  ParamDB->MIN_MAX_ADJOINT = 0;
  ParamDB->INITIAL_STEEPEST_DESCENT_ADJOINT = 0;

  tmp = new char[30];
  strcpy(tmp,"MooN_MD_default_basefile");
  ParamDB->BASENAME = tmp;
  ParamDB->VTKBASENAME = tmp;
  ParamDB->PSBASENAME = tmp;
  ParamDB->GRAPEBASENAME = tmp;
  ParamDB->GNUBASENAME = tmp;
  ParamDB->READGRAPEBASENAME = tmp;
  ParamDB->GMVBASENAME = tmp;
  ParamDB->MATLABBASENAME = tmp;
 
  
  tmp = new char[2];
  strcpy(tmp,"."); // current directory
  ParamDB->OUTPUTDIR = tmp;
  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_save_data_filename");
  ParamDB->SAVE_DATA_FILENAME=tmp;
  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_read_data_filename");
  ParamDB->READ_DATA_FILENAME=tmp;
  tmp = new char[40];
  strcpy(tmp, "NO_SMESH_FILE");
  ParamDB->SMESHFILE = tmp;
  tmp = new char[40];
  strcpy(tmp,"MooNMD_default_pod_filename");
  ParamDB->POD_FILENAME=tmp;
  //file for storing snapshots (ROM, reduced order modeling)
  tmp = new char[40];
  strcpy(tmp,"MooNMD_default_snap_filename");
  ParamDB->SNAP_FILENAME=tmp;


   /** parameters for SOLD schemes */
  ParamDB->SOLD_TYPE = 0;
  ParamDB->SOLD_PARAMETER_TYPE = 11;
  ParamDB->SOLD_CONST = 1.0;
  ParamDB->SOLD_POWER = 1.0;
  ParamDB->SOLD_S = 1.0;
  ParamDB->SOLD_U0 = 1.0;
  ParamDB->SOLD_PARAMETER_SCALING = 0;
  ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;

  /** parameters for controling the program */
  ParamDB->WRITE_PS = FALSE; 
  ParamDB->WRITE_GRAPE = FALSE; 
  ParamDB->WRITE_GMV = FALSE; 
  ParamDB->WRITE_AMIRA = FALSE; 
  ParamDB->WRITE_VTK = FALSE; 
  ParamDB->WRITE_GNU = FALSE; 
  ParamDB->SAVE_DATA = FALSE; 
  ParamDB->READ_DATA = FALSE; 
  ParamDB->READ_GRAPE_FILE = FALSE; 
  ParamDB->MEASURE_ERRORS = FALSE; 
  ParamDB->ESTIMATE_ERRORS = FALSE; 
  ParamDB->SOLVE_ADJOINT_PROBLEM = FALSE; 
  ParamDB->COMPUTE_VORTICITY_DIVERGENCE = FALSE;
  ParamDB->MESH_TYPE = 0; 
  ParamDB->USE_PRM = 1; 
   
  /** the following parameters are for individual use */
  ParamDB->P2 = 1.0;

  // ******** parameters for scalar system *********//
  ParamDB->SOLVER_TYPE = 1; 

  // parameters for nonlinear iteration
  ParamDB->SC_NONLIN_ITE_TYPE_SCALAR = 0;
  ParamDB->SC_NONLIN_MAXIT_SCALAR = 10000;
  ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR = 1e-10;
  ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR = 1.0;

  // parameters for linear iteration
  ParamDB->SC_SOLVER_SCALAR=AMG_GMRES_FLEX;
  ParamDB->SC_PRECONDITIONER_SCALAR=AMG_MGC;
  ParamDB->SC_LIN_MAXIT_SCALAR = 10000;
  ParamDB->SC_LIN_RED_FACTOR_SCALAR = 0.0;
  ParamDB->SC_LIN_RES_NORM_MIN_SCALAR = 1e-10;
  ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD = ParamDB->SC_LIN_RED_FACTOR_SCALAR;
  ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD = ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
  ParamDB->SC_LIN_MAXIT_SCALAR_SOLD = ParamDB->SC_LIN_MAXIT_SCALAR;
  ParamDB->SC_NONLIN_ITE_ADJOINT = 1000;
  ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER = 1;


  // parameters which are used in scalar multigrid
  ParamDB->SC_MG_TYPE_SCALAR = 0; 
  ParamDB->SC_MG_CYCLE_SCALAR = 1; 
  ParamDB->SC_SMOOTHER_SCALAR = 3;
  ParamDB->SC_PRE_SMOOTH_SCALAR= 2;
  ParamDB->SC_POST_SMOOTH_SCALAR = 2;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR = 1.0;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR = 1.0;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR = 1.0;
  ParamDB->SC_COARSE_SMOOTHER_SCALAR = 3;
  ParamDB->SC_COARSE_MAXIT_SCALAR = 10;
  ParamDB->SC_COARSE_RED_FACTOR_SCALAR =0.1;
  ParamDB->SC_GMG_DAMP_FACTOR_SCALAR = 1.0;
  ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR = 1.0;

  ParamDB->SC_COARSEST_LEVEL_SCALAR = 0;
  ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR = 0;

  ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR = 0;
  ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR = 1;

  // ******** parameters for saddle point system *********//

  // parameters for nonlinear iteration
  ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;
  ParamDB->SC_NONLIN_MAXIT_SADDLE = 1000;
  ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE = 1e-10;
  ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE = 1.0;
  ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE = 0;

  // parameters for linear iteration
  ParamDB->SC_SOLVER_SADDLE=AMG_GMRES_FLEX;
  ParamDB->SC_PRECONDITIONER_SADDLE=AMG_MGC;
  ParamDB->SC_LIN_MAXIT_SADDLE = 10000;
  ParamDB->SC_LIN_RED_FACTOR_SADDLE = 0.0;
  ParamDB->SC_LIN_RES_NORM_MIN_SADDLE = 1e-10;

  // parameters which are used in scalar multigrid
  ParamDB->SC_MG_TYPE_SADDLE = 0; 
  ParamDB->SC_MG_CYCLE_SADDLE = 1; 
  ParamDB->SC_SMOOTHER_SADDLE = 2;
  ParamDB->SC_PRE_SMOOTH_SADDLE= 2;
  ParamDB->SC_POST_SMOOTH_SADDLE = 2;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE = 1.0;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE = 1.0;
  ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE = 1.0;
  ParamDB->SC_COARSE_SMOOTHER_SADDLE = 2;
  ParamDB->SC_COARSE_MAXIT_SADDLE = 10;
  ParamDB->SC_COARSE_RED_FACTOR_SADDLE = 0.1;
  ParamDB->SC_GMG_DAMP_FACTOR_SADDLE = 1.0;
  ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE = 1.0;

  ParamDB->SC_COARSEST_LEVEL_SADDLE = 0;
  ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE = 0;

  ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE = 0;
  ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE = 0;

  ParamDB->SC_LARGEST_DIRECT_SOLVE = 203;
  ParamDB->SC_LARGEST_DIRECT_SOLVE = 100;
  ParamDB->SC_DOWNWIND_TYPE = 0;

  /** AMG solver parameters */
  // coarsen context 
  ParamDB->CC_ALPHA =  0.33333333;
  ParamDB->CC_BETA = 1.0E-5;
  ParamDB->CC_MINCLUSTER=4;
  ParamDB->CC_MAXCLUSTER=6;
  ParamDB->CC_MAXDISTANCE=2;
  ParamDB->CC_MAXCONNECTIVITY=15;
  ParamDB->CC_DEPTHTARGET=20;
  ParamDB->CC_COARSENTARGET=200;
  ParamDB->CC_COARSENRATE=1.2;
  ParamDB->CC_MAJOR=-1;
  ParamDB->CC_DEPENDENCY=AMG_UNSYM;
  ParamDB->CC_RESCALE=1.8;
  ParamDB->CC_VERBOSE=1;
  
  // THESE ARE THE DEFAULTS, DO NOT CHANGE 
  // solver context 
  ParamDB->SC_SYSTEM_TYPE=SCALAR;
  ParamDB->SC_AMG_PREC_IT= 1;
  ParamDB->SC_AMG_PREC_RED_FACTOR= 0.5;
  ParamDB->SC_EX_MAXIT = 0;
  ParamDB->SC_GMRES_RESTART = 10;
  ParamDB->SC_LCD_START_VECTOR = 0;
  ParamDB->SC_ILU_BETA=0.0;
  ParamDB->SC_SOR_OMEGA=1.5;
  ParamDB->SC_SMOOTHER_RED_FACTOR= 0.1;
  ParamDB->SC_OMEGA_COARSE_0=1.0;
  ParamDB->SC_OMEGA_P_0=1.0;
  ParamDB->SC_ILUT_TOL=0.01;
  ParamDB->SC_ILUT_ABSOLUTE_FILLIN=1;
  ParamDB->SC_ILUT_RELATIVE_FILLIN=1.0;
  ParamDB->SC_ILUT_SORT=ILUT_QUICK_SPLIT_0;
  ParamDB->SC_SCHUR_INV_OF_A= AMG_SSOR;
  ParamDB->SC_SCHUR_INV_OF_A_MAXIT= 1;
  ParamDB->SC_SCHUR_ITERATION_DAMP = 0.5;
  ParamDB->SC_SCHUR_ITERATION_MAXIT = 100;
  ParamDB->SC_SCHUR_STEP_LENGTH_CONTROL =0;
  ParamDB->SC_MIXED_BCGS_CGS_SWITCH_TOL=100;
  ParamDB->SC_DIV_FACTOR=1e10;
  ParamDB->SC_NONLIN_DIV_FACTOR=1e10;
  ParamDB->SC_SMOOTHING_STEPS=0;
  ParamDB->SC_N1_PARAM=1;
  ParamDB->SC_N2_PARAM=1;
  ParamDB->SC_MINIT=0;
  ParamDB->SC_VAS_LAZ_DELTA=1.0;
  ParamDB->SC_VERBOSE=1;
  ParamDB->SC_VERBOSE_AMG=1;
  ParamDB->SC_ROW_EQUILIBRATION = 0;

  ParamDB->SC_BRAESS_SARAZIN_MATRIX = 2;
  ParamDB->SC_BRAESS_SARAZIN_ALPHA = 1.5;

  ParamDB->TETGEN_QUALITY = 0.0;
  ParamDB->TETGEN_VOLUMEN = 0.0;
  ParamDB->TETGEN_STEINER = 0;

  ParamDB->CHAR_L0=1.;
  ParamDB->D_VISCOSITY=1.0;
  ParamDB->SURF_TENSION=0.;
  ParamDB->IMPACT_ANGLE=90.;
  ParamDB->Area = 1.;

  // initialize TimeDB
  TimeDB->CURRENTTIME = 0;
  TimeDB->CURRENTTIMESTEPLENGTH = 1;
  TimeDB->TIMESTEPLENGTH = 1;
  TimeDB->MIN_TIMESTEPLENGTH = 1E-4;
  TimeDB->MAX_TIMESTEPLENGTH = 0.5;
  TimeDB->TIMESTEPLENGTH_TOL = 1e-3;
  TimeDB->TIMESTEPLENGTH_CONTROL = 0;
  TimeDB->TIMESTEPLENGTH_CONTROLLER = 0;  // mlh
  TimeDB->TIMESTEPLENGTH_PARA_KK_I = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_P = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_E = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_R = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_D = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_FAC = 0.8;
  TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX = 5.0;
  TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN = 0.2;
  TimeDB->TIMESTEPLENGTH_PARA_TOL = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_ATOL = 0.001;
  TimeDB->TIMESTEPLENGTH_PARA_RTOL = 0.001;
  TimeDB->RESET_CURRENTTIME = 0;
  TimeDB->RESET_CURRENTTIME_STARTTIME = 0.0;
  TimeDB->STEADY_STATE_TOL = 1e-3;
  TimeDB->SCALE_DIVERGENCE_CONSTRAINT = -1.0;

  TimeDB->CONTROL=0;
  TimeDB->CONTROL_ALPHA=1;
  TimeDB->CONTROL_BETA=1;
  TimeDB->CONTROL_GAMMA=1;
  TimeDB->CONTROL_SAFTY=0.9;
  TimeDB->CONTROL_MINSCALE=0.1;
  TimeDB->CONTROL_MAXSCALE=5.0;
  
  // parameters for implicit Euler method
  TimeDB->THETA1 = 1;
  TimeDB->THETA2 = 0;
  TimeDB->THETA3 = 0;
  TimeDB->THETA4 = 1;
  TimeDB->TIME_DISC = 2;
  TimeDB->TIME_DISC2 = -1;

  TimeDB->STARTTIME = 0;
  TimeDB->ENDTIME = 1;
  TimeDB->EXTRAPOLATE_WEIGHT = 1;
  TimeDB->EXTRAPOLATE_STEPS = 0;
  TimeDB->EXTRAPOLATE_PRESSURE = 0;
  TimeDB->EXTRAPOLATE_VELOCITY = 0;

  TimeDB->T0 = 0;
  TimeDB->T1 = 0;
  TimeDB->T2 = 0;
  TimeDB->T3 = 0;
  TimeDB->T4 = 0;
  TimeDB->T5 = 0;
  TimeDB->T6 = 0;
  TimeDB->T7 = 0;
  TimeDB->T8 = 0;
  TimeDB->T9 = 0;

  TimeDB->STEPS_PER_IMAGE = 1;
  TimeDB->STEPS_PER_SNAP = 1;

  TimeDB->RB_TYPE = 3;
  TimeDB->RB_TYPE2 = -1;
  TimeDB->RB_SSC = 0;
  TimeDB->RB_SSC_TOL = 1;
  TimeDB->RB_SSC_ALPHA = 0.9;
  TimeDB->RB_SSC_ALPHA_MIN = 1.0;
  TimeDB->RB_SSC_ALPHA_MAX = 1.0;
  TimeDB->RB_SSC_MAX_ERROR = 1.0;

  TimeDB->RB_APPROX_J = 0;
  TimeDB->RB_APPROX_C = 0;
  
  // parameters for higher order Galerkin-type methods
  TimeDB->INTERNAL_SYSTEMSIZE = 0;
  TimeDB->INTERNAL_ALPHA = NULL;
  TimeDB->INTERNAL_BETA = NULL;
  TimeDB->VALUE_AT_ONE = NULL;
  TimeDB->VAL_AT_QUAD_POINTS = NULL;
  TimeDB->DER_AT_QUAD_POINTS = NULL;
  TimeDB->CORR_AT_QUAD_POINTS = NULL;
  TimeDB->DER_CORR_AT_QUAD_POINTS = NULL;
  TimeDB->DER_AT_START = NULL;
  TimeDB->DER_AT_ONE = NULL;
  TimeDB->DER_COR_AT_ONE = NULL;
  TimeDB->NORMW = 0;
  TimeDB->N_QUADPOINTS = 0;
  TimeDB->N_DEGREES = 0;
  TimeDB->ZETA = NULL;
  TimeDB->WEIGHTS = NULL;
  TimeDB->ALPHA0 = NULL;
  TimeDB->BETA0 = NULL;
  TimeDB->GAMMA0 = NULL;
  TimeDB->CORRECTION = NULL;
  TimeDB->POINTS = NULL;

  TimeDB->DG_TimeDisc = 0; 
  TimeDB->DG_Order = 0; 
    
  ParamDB->INPUT_QUAD_RULE = 0;
  ParamDB->INTERNAL_PROBLEM_LINEAR = 0;
  ParamDB->INTERNAL_PRESSURE_SPACE = 0;
  ParamDB->INTERNAL_PROJECT_PRESSURE = 1;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION = 0;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 0;
  ParamDB->INTERNAL_QUAD_HEXA = 0;
  ParamDB->INTERNAL_QUAD_TETRA = 0;
  ParamDB->INTERNAL_QUAD_QUAD = 0;
  ParamDB->INTERNAL_QUAD_TRIA = 0;
  ParamDB->INTERNAL_QUAD_RULE = 0;
  ParamDB->INTERNAL_LOCAL_DOF = 0;
  ParamDB->INTERNAL_PERIODIC_IDENTITY = 0;
  ParamDB->INTERNAL_PROBLEM_IDENTITY = 0;
  ParamDB->INTERNAL_LEVEL = 0;
  ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 0;
  ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
  ParamDB->INTERNAL_AMG_SOLVES = 0;
  ParamDB->INTERNAL_AMG_PREPARE_TIME = 0.0;
  ParamDB->INTERNAL_GMRES_INFO = 1;
  ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
  ParamDB->INTERNAL_LINEAR_SCHEME = 1;
  ParamDB->INTERNAL_SOLD_ACTIVE = 0;
  ParamDB->INTERNAL_SORT_AMG = 1;
  ParamDB->INTERNAL_COERCIVITY = -4711.0;
  ParamDB->INTERNAL_FACE_INTEGRALS = 0;
  ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS = 0;
  ParamDB->INTERNAL_WRONG_NEUMANN_CHECKED = 0;
  ParamDB->INTERNAL_BFGS_RESTART_ADJOINT = 0;
  ParamDB->INTERNAL_NEW_MATRICES_B = 1;
  ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
  ParamDB->INTERNAL_DISC_FLAG = 0;

  ParamDB->INTERNAL_FESPACE_CONSTRUCT = 0;
  ParamDB->INTERNAL_DO_NOT_RESPECT_DIRICHLET_BC = 0;
  ParamDB->INTERNAL_P1_Array = NULL;
  ParamDB->INTERNAL_START_PARAM = 0;

  ParamDB->MESH_SLIP_WITH_FRICTION = 0;

  /** parameters for free surface calculation */
  ParamDB->INTERFACE_FLOW = FALSE;
  
  ParamDB->FS_MAGNETLAW = 0;

  ParamDB->FS_L = 1;
  ParamDB->FS_U = 1;
  ParamDB->FS_T = 1;

  ParamDB->FS_ETA = 1;
  ParamDB->FS_RHO = 1;
  ParamDB->FS_ALPHA = 1;
  ParamDB->FS_G = 1;

  ParamDB->FS_MS = 1;
  ParamDB->FS_CHI0 = 1;

  ParamDB->FS_HM = 1;
  ParamDB->FS_DELTA_H = 0;
  ParamDB->FS_F = 0;
  
  ParamDB->FS_H0 = 1;
  ParamDB->FS_H1 = 1;
  ParamDB->FS_H2 = 1;
  
  ParamDB->FS_LH = 1;
  ParamDB->FS_GAMMA = 1;
  ParamDB->FS_HT = 1;

  ParamDB->FS_WE = 1;

  ParamDB->FS_WRITE = 0;
  ParamDB->FS_READ = 0;

  tmp = new char[12];
  strcpy(tmp,"FS_INNAME");
  ParamDB->FS_INNAME = tmp;

  tmp = new char[12];
  strcpy(tmp,"FS_OUTNAME");
  ParamDB->FS_OUTNAME = tmp;

  ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR = 0;  
  ParamDB->HEAT_SOLID_SURFACE_FACTOR = 1.;    
  ParamDB->EQ_CONTACT_ANGLE = 0;   
  ParamDB->AD_CONTACT_ANGLE = 0;   
  ParamDB->RE_CONTACT_ANGLE = 0;  
  ParamDB->DY_CONTACT_ANGLE = 0;     
  ParamDB->CONTACT_ANGLE_TYPE = 0;  
   
  // ******** parameters for VMS *********//
  ParamDB->VMS_LARGE_VELOCITY_SPACE = 0;
  ParamDB->VMS_COARSE_MG_SMAGO = 1;
  // constants in AdaptProjectionSpace 
  ParamDB->VMS_ADAPT_LOWER = 0.25;
  ParamDB->VMS_ADAPT_MIDDLE = 1.0;
  ParamDB->VMS_ADAPT_UPPER = 3.0;
  ParamDB->VMS_ADAPT_STEPS = 1;
  ParamDB->VMS_ADAPT_COMP = 1;

  ParamDB->SUPERCONVERGENCE_ORDER = 0;
  ParamDB->FEM_FCT_LINEAR_TYPE = 1;
  ParamDB->FEM_FCT_PRELIMITING = 0;
  ParamDB->FEM_FCT_GROUP_FEM = 0;
  ParamDB->GROUP_FEM = 0;
  ParamDB->WENO_TYPE = 0;

  /* the following parameters are for membrane REACTOR */
  ParamDB->REACTOR_P0 = 0.0;
  ParamDB->REACTOR_P1 = 0.0;
  ParamDB->REACTOR_P2 = 0.0;
  ParamDB->REACTOR_P3 = 0.0;
  ParamDB->REACTOR_P4 = 0.0;
  ParamDB->REACTOR_P5 = 0.0;
  ParamDB->REACTOR_P6 = 0.0;
  ParamDB->REACTOR_P7 = 0.0;
  ParamDB->REACTOR_P8 = 0.0;
  ParamDB->REACTOR_P9 = 0.0;
  ParamDB->REACTOR_P10 = 0.0;
  ParamDB->REACTOR_P11 = 0.0;
  ParamDB->REACTOR_P12 = 0.0;
  ParamDB->REACTOR_P13 = 0.0;
  ParamDB->REACTOR_P14 = 0.0;
  ParamDB->REACTOR_P15 = 0.0;
  ParamDB->REACTOR_P16 = 0.0;
  ParamDB->REACTOR_P17 = 0.0;
  ParamDB->REACTOR_P18 = 0.0;
  ParamDB->REACTOR_P19 = 0.0;
  ParamDB->REACTOR_P20 = 0.0;
  ParamDB->REACTOR_P21 = 0.0;
  ParamDB->REACTOR_P22 = 0.0;
  ParamDB->REACTOR_P23 = 0.0;
  ParamDB->REACTOR_P24 = 0.0;
  ParamDB->REACTOR_P25 = 0.0;
  ParamDB->REACTOR_P26 = 0.0;
  ParamDB->REACTOR_P27 = 0.0;
  ParamDB->REACTOR_P28 = 0.0;
  ParamDB->REACTOR_P29 = 0.0;
  ParamDB->REACTOR_P30 = 0.0;

// parameters for turbulent channel flow
  ParamDB->CHANNEL_STATISTICS2_WITH_MODEL = 0;

// parameters for turbulent flow around a squared cylinder
  ParamDB->CYLINDER_22000_YPLUS_SIDES = 1500;
  ParamDB->CYLINDER_22000_YPLUS_FRONT = 2300;
  ParamDB->CYLINDER_22000_YPLUS_BACK  = 1000;
  
// parameters for BULK computations
  ParamDB->BULK_REACTION_DISC = 1;
  ParamDB->BULK_PB_DISC = 3;
  ParamDB->BULK_PB_DISC_STAB = 1;
  ParamDB->BULK_PB_DISC_FCT_GROUP = 1;
  ParamDB->BULK_COUPLING = 1;
  ParamDB->BULK_GROWTH_RATE = 1;
  ParamDB->BULK_REACTION_MASS_LUMPING = 0;
  ParamDB->BULK_REACTION_C_CUT = 6;
  ParamDB->BULK_METHODS_OF_MOMENTS = 0;
  ParamDB->BULK_MOM_DISC = 2;
  ParamDB->BULK_SOLD_PARAMETER_TYPE = 51;
  ParamDB->N_CELL_LAYERS_PSD = 16;
  ParamDB->N_CELL_LAYERS_PSD_2 = 16;
  ParamDB->OUTPUT_NODE_LAYER_PSD = 1;

// coefficients for the fluid
  ParamDB->BULK_density = 1000;
  ParamDB->BULK_dynamic_viscosity = 1e-3;  

  ParamDB->BULK_l_infty = 1;
  ParamDB->BULK_u_infty = 1e-3;
  ParamDB->BULK_c_infty = 1;
  ParamDB->BULK_c_C_infty_sat = 1.37e-4;

  ParamDB->BULK_C_g = 45.98;
  ParamDB->BULK_C_nuc = 15.33;
  ParamDB->BULK_C_sat = 1.37e-4;
  ParamDB->BULK_C_2 = 7.2e-9;
  ParamDB->BULK_D_A = 1.5e-9;
  ParamDB->BULK_D_P_0 = 1e-9;
  ParamDB->BULK_D_P_MAX = 1e-3;
  ParamDB->BULK_k_g = 1e-7;
  ParamDB->BULK_k_nuc = 1e24;
  ParamDB->BULK_k_r = 1e-2;

// parameters for shear slip mesh update method computations
  ParamDB->SSMUM_MP_X = 0.5;
  ParamDB->SSMUM_MP_Y = 0.5;
  ParamDB->SSMUM_INNER_RADIUS = 0.25;
  ParamDB->SSMUM_OUTER_RADIUS = 0.35;
  ParamDB->SSMUM_ROT_PER_SECOND = 0.5;
  ParamDB->SSMUM_ANGLE = 0.0;
  ParamDB->SSMUM_MAX_CELLS_LAYERS = 1024;
  ParamDB->SSMUM_INTERPOLATION = 0;

  // parameters for WINDTUNNEL computations 
  ParamDB->WINDTUNNEL_CONFIGURATION = 1;
  ParamDB->WINDTUNNEL_INTERPOLATION = 0;
  ParamDB->WINDTUNNEL_STEADY = 0;
  ParamDB->WINDTUNNEL_SPATIAL = 3;
  ParamDB->WINDTUNNEL_BROWNIAN = 0;
  ParamDB->WINDTUNNEL_POL_ORDER = 0;
  ParamDB->WINDTUNNEL_SHEAR_FACTOR_TYPE = 0;
  ParamDB->WINDTUNNEL_SHEAR_FACTOR = 0.05;
  ParamDB->WINDTUNNEL_QUAD_METHOD=0;
  ParamDB->WINDTUNNEL_MEASURE_MASS=0;
  ParamDB->WINDTUNNEL_SHIFT=0.;
  
  ParamDB->WINDTUNNEL_LAYER_NUMBER_X = WINDTUNNEL_LAYER_NUMBER_X_CONST-1;
  ParamDB->WINDTUNNEL_DIM_Y = WINDTUNNEL_DIM_Y_CONST-1;
  ParamDB->WINDTUNNEL_DIM_Z = WINDTUNNEL_DIM_Z_CONST-1;
  ParamDB->WINDTUNNEL_DIM_R = WINDTUNNEL_DIM_R_CONST-1;
  ParamDB->WINDTUNNEL_ENVIR_COND = 5.0613e-10;
  ParamDB->WINDTUNNEL_SUPERSAT = 0.01 ; //test normally supersaturation =0.01
  ParamDB->WINDTUNNEL_U_INFTY = 1; // m/s
  ParamDB->WINDTUNNEL_L_INFTY = 1; // m
  ParamDB->WINDTUNNEL_R_MIN = 0e-6; // = r_min in m
  ParamDB->WINDTUNNEL_R_INFTY = 175e-6;  // log -normal250e-6; // = r_max in m
  ParamDB->WINDTUNNEL_F_INFTY = 1e12; // #/m^4
  //ParamDB->WINDTUNNEL_kinematic_viscosity = 15.68e-6;
  ParamDB->WINDTUNNEL_dynamic_viscosity = 18.15e-6;  // air (Rogers, Yau p. 103)
  ParamDB->WINDTUNNEL_density = 1.2041;   // air
  // ParamDB->WINDTUNNEL_BOUND_KOEFF=0;


  ParamDB->UREA_REACTION_DISC = 1;
  ParamDB->UREA_PB_DISC = 3;
  ParamDB->UREA_MODEL = 0;
  ParamDB->UREA_PB_DISC_STAB = 1;
  ParamDB->UREA_SOLD_PARAMETER_TYPE = 51;
  ParamDB->UREA_PIPE = 0;

  ParamDB->UREA_l_infty = 0.01;
  ParamDB->UREA_u_infty = 0.01;
  ParamDB->UREA_c_infty = 1000;
  ParamDB->UREA_temp_infty = 1;
  ParamDB->UREA_f_infty = 1e13;
  ParamDB->UREA_nu = 1.3612e-6;
  ParamDB->UREA_rho = 789;
  ParamDB->UREA_c_p =  2441.3;
  ParamDB->UREA_lambda = 0.167;
  ParamDB->UREA_D_P_0 = 2.5e-6;
  ParamDB->UREA_D_P_MAX = 5000e-6;
  ParamDB->UREA_k_v = Pi/6.0;
  ParamDB->UREA_m_mol = 60.06e-3;
  ParamDB->UREA_D_J = 1.35e-9;
  ParamDB->UREA_rho_d = 1323;
  ParamDB->UREA_delta_h_cryst = 0.21645e3;
  ParamDB->UREA_k_g = 1e-7;
  //ParamDB->UREA_k_g = 0.;
  ParamDB->UREA_g = 0.5;
  ParamDB->UREA_rho_sat_1 = 35.364;
  ParamDB->UREA_rho_sat_2 = 1.305;
  ParamDB->UREA_beta_nuc = 0.166667e-5;
  //with nucleation
  ParamDB->UREA_alfa_nuc = 1.e8;
  //without nucleation
  //ParamDB->UREA_alfa_nuc = 0.;
  ParamDB->UREA_INFLOW_SCALE = 4.5e-2; // centimeter
  ParamDB->UREA_CONC_TOL = 1e-6;
  ParamDB->UREA_CONC_MAXIT = 1;
  ParamDB->UREA_inflow_time = 5;

  ParamDB->UREA_AGGR_SPATIAL = 3;              //spatial dimension
  ParamDB->UREA_AGGR_BROWNIAN = 0.;            //include brownian kernel
  ParamDB->UREA_AGGR_BROWNIAN_TEMP = 0;        //include temp in brownian kernel
  ParamDB->UREA_AGGR_BROWNIAN_SCAL = 0.;       //scal for brownian kernel
  ParamDB->UREA_AGGR_POL_ORDER = 0.;           //degree of the polynomial basis (at the moment 0 (constant) or 1 (linear))
  ParamDB->UREA_AGGR_SHEAR_FACTOR_TYPE = 0.;   //shear induced kernel factor type (0 - constant factor, 1 - depends on the velocity)
  ParamDB->UREA_AGGR_SHEAR_FACTOR = 0.05;      //the factor itself

  //Param KDP model

  ParamDB->KDP_MODEL = 0;
  ParamDB->KDP_l_infty = 0.01;
  ParamDB->KDP_u_infty = 0.01;
  ParamDB->KDP_c_infty = 1;
  ParamDB->KDP_temp_infty = 1;
  ParamDB->KDP_f_infty = 1e13;
  ParamDB->KDP_nu = 1.2931e-6;
  ParamDB->KDP_rho = 1160;
  ParamDB->KDP_c_p =  4181.3;
  ParamDB->KDP_lambda = 0.602;
  ParamDB->KDP_D_P_0 = 0.0;//1.25e-6;//2.5e-6;
  ParamDB->KDP_D_P_0_2 =0.0;//1.25e-6;//2.5e-6;
  // ParamDB->KDP_D_P_MAX = 2500e-6;
  ParamDB->KDP_D_P_MAX = 1e-3;
  //ParamDB->KDP_D_P_MAX_2 = 5000e-6;
  ParamDB->KDP_D_P_MAX_2 = 1e-3;
  ParamDB->KDP_m_mol = 136.08e-3;
  ParamDB->KDP_D_J = 5.5e-10;
  ParamDB->KDP_rho_d = 2338;
  ParamDB->KDP_delta_h_cryst = 0.119e3;
  // ParamDB->KDP_delta_h_cryst = 119e3;
  ParamDB->KDP_k_g_1 = 1.221e-5;
  ParamDB->KDP_k_g_2 = 10.075e-5;
  //nach christian
  ParamDB->KDP_k_b = 7.875e9;//7.49e10;
  //ParamDB->KDP_k_b = 3.75e13;
  ParamDB->KDP_g_1 = 1.48;
  ParamDB->KDP_g_2 = 1.74;
  ParamDB->KDP_b = 2.04;
  ParamDB->KDP_w_sat_1 = 5.5843e-5;
  ParamDB->KDP_w_sat_1_Ma = 9.3027e-5;
  ParamDB->KDP_w_sat_2 = 2.8159e-2;
  ParamDB->KDP_w_sat_2_Ma = 9.7629e-5;
  ParamDB->KDP_w_sat_3 = 3.6832;
  ParamDB->KDP_w_sat_3_Ma = 0.2087;
  ParamDB->KDP_INTERNAL_NUC_A =0.0;
  ParamDB->KDP_INTERNAL_NUC_B =0.0;
  ParamDB->KDP_INFLOW_SCALE = 4.5e-2; // centimeter
  ParamDB->KDP_CONC_TOL = 1e-6;
  ParamDB->KDP_CONC_MAXIT = 1;
  ParamDB->KDP_inflow_time = 10;
    
    
  //======================================================================
  /** parameters for Stokes--Darcy (StoDa) coupling */
  //======================================================================
  ParamDB->StoDa_interfaceType = 0; //Beavers-Joseph-Saffman or u.t=0
  ParamDB->StoDa_alpha = 1; // from Beavers-Joseph-Saffman condition on interface
  ParamDB->StoDa_problemType = 1;// Neumann--Neumann, Robin--Robin, ...
  ParamDB->StoDa_updatingStrategy = 0; // update of the etas
  ParamDB->StoDa_theta_f = 1; //damping in Stokes (flow) part
  ParamDB->StoDa_theta_p = 1; //damping in Darcy (porous) part
  ParamDB->StoDa_gamma_f = 1; // parameter for Robin condition on interface
  ParamDB->StoDa_gamma_p = 1; // parameter for Robin condition on interface
  ParamDB->StoDa_weakGamma = 1; // parameter for enforcing weak boundary conditions
  ParamDB->StoDa_solutionStrategy = 1; // only iterative (0) or iterative and one big matrix (1)
  ParamDB->StoDa_algorithm = 1; // Gauss--Seidel, Jacobi, ...
  ParamDB->StoDa_StokesFirst = 0; // for Gauss--Seidel type method.
  ParamDB->StoDa_nIterations = 100; // maximum number of iterations
  ParamDB->StoDa_relDiff_interfaceError = 1e-10; // (e_k - e_{k+1})/e_k < this number
  ParamDB->StoDa_relDiff_factor1 = 1; // factor of StokesU in computing E_k
  ParamDB->StoDa_relDiff_factor2 = 1; // factor of StokesP in computing E_k
  ParamDB->StoDa_relDiff_factor3 = 1; // factor of DarcyP  in computing E_k
  ParamDB->StoDa_relDiff_solution = 1e-10; // E_k < this number
  ParamDB->StoDa_bigResidual = 1e-10; // residual of big System < this number
  ParamDB->StoDa_periodicBoundary = 0; // true if there is a periodic boundary
  ParamDB->StoDa_periodicBoundaryPressureDrop = 1.0; // pressure drop at periodic boundary
  
  /** general parameters for population balances */
  ParamDB->PB_DISC_TYPE = 3;
  ParamDB->PB_TIME_DISC = 100;

  
  /** parameters for matlab output */
  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_matlab_matrix");
  ParamDB->MATLAB_MATRIX=tmp;
  
  ParamDB->WRITE_MATLAB = FALSE;
  ParamDB->WRITE_MATLAB_MATRIX = FALSE;

  /** parameter for non-conforming elements
      number corresponds to Apel-Matthies paper 2006 */
  ParamDB->NC_TYPE = 3;
  
  /** parameters for ROM*/
  ParamDB->WRITE_SNAPSHOTS = FALSE;
  ParamDB->DO_ROM = FALSE;
  ParamDB->DO_ROM_P = FALSE;
  ParamDB->RANK_OF_BASIS = 0;
  ParamDB->RANK_OF_BASIS_P = 0;
  ParamDB->POD_INNER_PRODUCT = 1;
  ParamDB->POD_INNER_PRODUCT_P = 1;
  
  ParamDB->BUILD_PODFILE = FALSE;
  ParamDB->POD_FLUCT_FIELD = FALSE;
  ParamDB->POD_FLUCT_FIELD_P = FALSE;
  ParamDB->P_ROM_METHOD = 1;

  tmp = new char[12];
  strcpy(tmp,"NO_PODFILE");
  ParamDB->PODFILE= tmp;
  
  /** parameter for type of projection method (NSE)**/
  ParamDB->PROJECTION_METHOD = 1;

  /** parameters for individual use in parallel computations */
  ParamDB->Par_P0 = 0;
  ParamDB->Par_P1 = 0;
  ParamDB->Par_P2 = 0;
  ParamDB->Par_P3 = 1;
  ParamDB->Par_P4 = 0;
  ParamDB->Par_P5 = 10;
  ParamDB->Par_P6 = 1;
  ParamDB->Par_P7 = 0;
  ParamDB->Par_P8 = 0;
  ParamDB->Par_P9 = 0;
  ParamDB->Par_P10 = -1.0;
  ParamDB->Par_P11 = 0;
  ParamDB->Par_P12 = 0;
  ParamDB->Par_P13 = 0;
  ParamDB->Par_P14 = 0;
  ParamDB->Par_P15 = 0;
  ParamDB->Par_P16 = 0;
  ParamDB->Par_P17 = 0;
  ParamDB->Par_P18 = 0;
  ParamDB->Par_P19 = 0;
  ParamDB->Par_P20 = 0;
  
  ParamDB->MG_DEBUG = 0;
  ParamDB->DOF_Reorder =0;
  ParamDB->DOF_Average =1;
  ParamDB->SC_LOCAL_SMOOTH=1;
  
  
  /** parameters for population balance computations */
  ParamDB->PBE_P0 = 0;
  ParamDB->PBE_P1 = 0;
  ParamDB->PBE_P2 = 0;
  ParamDB->PBE_P3 = 1;
  ParamDB->PBE_P4 = 0;
  ParamDB->PBE_P5 = 0;
  ParamDB->PBE_P6 = 0;
  ParamDB->PBE_P7 = 0;
  ParamDB->PBE_P8 = 0;
  ParamDB->PBE_P9 = 0;
  
   /** parameters for population balance computations */
  ParamDB->DG_P0 = 0.;
  ParamDB->DG_P1 = 0.;
  ParamDB->DG_P2 = 0.;
  ParamDB->DG_P3 = 0.;
  ParamDB->DG_P4 = 0.;
  ParamDB->DG_P5 = 0.;
  ParamDB->DG_P6 = 0.;
  ParamDB->DG_P7 = 0.;
  ParamDB->DG_P8 = 0.;
  ParamDB->DG_P9 = 0.;
  
  /** parameters for moving domains */
  ParamDB->MOVING_BOUNDARY = 0;
  ParamDB->ASSEMBLEMESHMAT = FALSE;
  ParamDB->LameC = 1.0;
  
  /** parameters for moving domains */
  ParamDB->DEPENDENT_BASIS = 0;
  ParamDB->DEPENDENT_BASIS_Q1 = 0;
  ParamDB->DEPENDENT_BASIS_Q2 = 0;
  
  #if defined(_MPI) || defined(_SMPI)
  ParamDB->Comm = MPI_COMM_WORLD;    
  #endif
  return;
}

void TDatabase::WriteParamDB(char *ExecutedFile)
{
  char buf[80];
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
 
  gethostname(buf,80);
  OutFile << "HOSTNAME: " << buf << " started on " << asctime (timeinfo)  << endl;
  OutFile << "EXECUTED FILE: " <<  ExecutedFile << endl;
  OutFile << "VERSION: " << ParamDB->VERSION << endl;
  OutFile << "GEOFILE: " << ParamDB->GEOFILE << endl;
  OutFile << "BNDFILE: " << ParamDB->BNDFILE << endl;
  OutFile << "MAPFILE: " << ParamDB->MAPFILE << endl;
  OutFile << "OUTFILE: " << ParamDB->OUTFILE << endl;
  OutFile << "PROBLEM_TYPE: " << ParamDB->PROBLEM_TYPE << endl;
  OutFile << "EXAMPLE: " << ParamDB->EXAMPLE << endl;
  OutFile << "profiling: " << ParamDB->timeprofiling << endl;
  OutFile << "MapperType: " << ParamDB->MapperType << endl;
  OutFile << "DSType: " << ParamDB->DSType << endl;
  
  OutFile << "MESHGEN_ALLOW_EDGE_REF: " << ParamDB->MESHGEN_ALLOW_EDGE_REF << endl;
  OutFile << "MESHGEN_REF_QUALITY: " << ParamDB->MESHGEN_REF_QUALITY << endl;

  OutFile << "ANSATZ_ORDER: " << ParamDB->ANSATZ_ORDER << endl;
  OutFile << "TEST_ORDER: " << ParamDB->TEST_ORDER << endl;
  
  OutFile << "VELOCITY_SPACE: " << ParamDB->VELOCITY_SPACE << endl;
  OutFile << "PRESSURE_SPACE: " << ParamDB->PRESSURE_SPACE << endl;
  OutFile << "PRESSURE_SEPARATION: " << ParamDB->PRESSURE_SEPARATION << endl;

  OutFile << "OMPNUMTHREADS: " << ParamDB->OMPNUMTHREADS << endl;

  OutFile << "LEVELS: " << ParamDB->LEVELS << endl;
  OutFile << "N_CELL_LAYERS: " << ParamDB->N_CELL_LAYERS << endl;
  OutFile << "UNIFORM_STEPS: " << ParamDB->UNIFORM_STEPS << endl;
  OutFile << "DRIFT_X: " << ParamDB->DRIFT_X << endl;
  OutFile << "DRIFT_Y: " << ParamDB->DRIFT_Y << endl;
  OutFile << "DRIFT_Z: " << ParamDB->DRIFT_Z << endl;
  
  OutFile << "NONLINEARIT_TYPE_NEWTON: " << ParamDB->NONLINEARIT_TYPE_NEWTON << endl;
  
  OutFile << "REFINEMENT: " << ParamDB->REFINEMENT << endl;
  OutFile << "GRID_TYPE: " << ParamDB->GRID_TYPE << endl;
  OutFile << "GRID_TYPE_1: " << ParamDB->GRID_TYPE_1 << endl;
  OutFile << "GRID_TYPE_2: " << ParamDB->GRID_TYPE_2 << endl;
  OutFile << "CHANNEL_GRID_STRETCH: " << ParamDB->CHANNEL_GRID_STRETCH << endl;

  OutFile << "ADAPTIVE_REFINEMENT_CRITERION: " << ParamDB->ADAPTIVE_REFINEMENT_CRITERION << endl;
  OutFile << "ERROR_CONTROL: " << ParamDB->ERROR_CONTROL << endl;
  OutFile << "REFINE_STRATEGY: " << ParamDB->REFINE_STRATEGY << endl;
  OutFile << "MAX_CELL_LEVEL: " << ParamDB->MAX_CELL_LEVEL << endl;
  OutFile << "REFTOL: " << ParamDB->REFTOL << endl;
  OutFile << "COARSETOL: " << ParamDB->COARSETOL << endl;  
  OutFile << "MIN_FRACTION_TO_CHANGE: " << ParamDB->MIN_FRACTION_TO_CHANGE << endl;
  OutFile << "DECREASE_REFTOL_FACTOR: " << ParamDB->DECREASE_REFTOL_FACTOR << endl;
  OutFile << "INCREASE_COARSETOL_FACTOR: " << ParamDB->INCREASE_COARSETOL_FACTOR << endl;
  OutFile << "FRACTION_OF_ERROR: " << ParamDB->FRACTION_OF_ERROR << endl;
  OutFile << "CONVERT_QUAD_TO_TRI: " << ParamDB->CONVERT_QUAD_TO_TRI << endl;
  
  OutFile << "DISCTYPE: " << ParamDB->DISCTYPE << endl;
  OutFile << "INTL_DISCTYPE: " << ParamDB->INTL_DISCTYPE << endl;
  OutFile << "UPWIND_ORDER: " << ParamDB->UPWIND_ORDER << endl;
  OutFile << "UPWIND_FLUX_DAMP: " << ParamDB->UPWIND_FLUX_DAMP << endl;
  OutFile << "UPWIND_APPLICATION: " << ParamDB->UPWIND_APPLICATION<< endl;
  OutFile << "SHISHKIN_MESH: " << ParamDB->SHISHKIN_MESH << endl;
  OutFile << "SHISHKIN_DIAM: " << ParamDB->SHISHKIN_DIAM << endl;
  OutFile << "NSTYPE: " << ParamDB->NSTYPE << endl;
  OutFile << "DARCYTYPE: " << ParamDB->DARCYTYPE << endl;
  OutFile << "SIGMA_PERM: " << ParamDB->SIGMA_PERM << endl;
  OutFile << "LAPLACETYPE: " << ParamDB->LAPLACETYPE << endl;
  OutFile << "TENSOR_TYPE: " << ParamDB->TENSOR_TYPE << endl;
  OutFile << "TWO_PHASE_FLOW: " << ParamDB->TWO_PHASE_FLOW << endl;
  OutFile << "PHASE1_TYPE: " << ParamDB->PHASE1_TYPE << endl;
  OutFile << "PHASE2_TYPE: " << ParamDB->PHASE2_TYPE << endl;
  OutFile << "USE_ISOPARAMETRIC: " << ParamDB->USE_ISOPARAMETRIC << endl;
  OutFile << "VMM_COARSE_LEVEL: " << ParamDB->VMM_COARSE_LEVEL << endl;
  OutFile << "VMM_COARSE_SPACE_ORDER: " << ParamDB->VMM_COARSE_SPACE_ORDER << endl;
  OutFile << "RFB_SUBMESH_LAYERS: " << ParamDB->RFB_SUBMESH_LAYERS << endl;
  OutFile << "DEFECT_CORRECTION_TYPE: " << ParamDB->DEFECT_CORRECTION_TYPE<< endl;
  OutFile << "CELL_MEASURE: " << ParamDB->CELL_MEASURE<< endl;
  OutFile << "FACE_SIGMA: " << ParamDB->FACE_SIGMA << endl;
  OutFile << "WEAK_BC_SIGMA: " << ParamDB->WEAK_BC_SIGMA << endl;
  OutFile << "WEAK_BC: " << ParamDB->WEAK_BC << endl;

  OutFile << "RE_NR: " << ParamDB->RE_NR << endl;
  OutFile << "RA_NR: " << ParamDB->RA_NR << endl;
  OutFile << "ROSSBY_NR: " << ParamDB->ROSSBY_NR << endl;
  OutFile << "START_RE_NR: " << ParamDB->START_RE_NR << endl;
  OutFile << "RE_NR_INCREMENT: " << ParamDB->RE_NR_INCREMENT << endl;
  OutFile << "FLOW_PROBLEM_TYPE: " << ParamDB->FLOW_PROBLEM_TYPE << endl;
  OutFile << "FR_NR: " << ParamDB->FR_NR << endl;
  OutFile << "WB_NR: " << ParamDB->WB_NR << endl;
  OutFile << "PR_NR: " << ParamDB->PR_NR << endl;
  OutFile << "PE_NR: " << ParamDB->PE_NR << endl;  
  OutFile << "BI_NR: " << ParamDB->BI_NR << endl;
  OutFile << "WEI_NR: " << ParamDB->WEI_NR << endl;
  OutFile << "LameC: " << ParamDB->LameC << endl;
  OutFile << "Axial3D: " << ParamDB->Axial3D << endl;
  OutFile << "Axial3DAxis: " << ParamDB->Axial3DAxis << endl;  
  OutFile << "DELTA0: " << ParamDB->DELTA0 << endl;
  OutFile << "SDFEM_POWER0: " << ParamDB->SDFEM_POWER0 << endl;
  OutFile << "SDFEM_TYPE: " << ParamDB->SDFEM_TYPE << endl;
  OutFile << "SDFEM_NORM_B: " << ParamDB->SDFEM_NORM_B << endl;
  OutFile << "DELTA1: " << ParamDB->DELTA1 << endl;
  OutFile << "CIP_TYPE: " << ParamDB->CIP_TYPE << endl;
  OutFile << "ADJOINT_FACTOR_4_OMEGA_EQ_0: " << ParamDB->ADJOINT_FACTOR_4_OMEGA_EQ_0 << endl;
  OutFile << "DELTA2: " << ParamDB->DELTA2 << endl;

  OutFile << "SOLD_TYPE: " << ParamDB->SOLD_TYPE << endl;
  OutFile << "SOLD_PARAMETER_TYPE: " << ParamDB->SOLD_PARAMETER_TYPE << endl;
  OutFile << "SOLD_CONST: " << ParamDB->SOLD_CONST << endl;
  OutFile << "SOLD_POWER: " << ParamDB->SOLD_POWER << endl;
  OutFile << "SOLD_S: " << ParamDB->SOLD_S << endl;
  OutFile << "SOLD_U0: " << ParamDB->SOLD_U0 << endl;
  OutFile << "SOLD_PARAMETER_SCALING: " << ParamDB->SOLD_PARAMETER_SCALING << endl;
  OutFile << "SOLD_PARAMETER_SCALING_FACTOR: " << ParamDB->SOLD_PARAMETER_SCALING_FACTOR << endl;
  

  OutFile << "FILTER_WIDTH_CONSTANT: " << ParamDB->FILTER_WIDTH_CONSTANT << endl;
  OutFile << "FILTER_WIDTH_POWER: " << ParamDB->FILTER_WIDTH_POWER << endl;
  OutFile << "GAUSSIAN_GAMMA: " << ParamDB->GAUSSIAN_GAMMA  << endl;
  OutFile << "CONVOLUTE_SOLUTION: " << ParamDB->CONVOLUTE_SOLUTION << endl;

  OutFile << "TURBULENT_VISCOSITY_TYPE: " << ParamDB->TURBULENT_VISCOSITY_TYPE << endl;
  OutFile << "TURBULENT_VISCOSITY_TENSOR: " << ParamDB->TURBULENT_VISCOSITY_TENSOR << endl;
  OutFile << "TURBULENT_VISCOSITY_CONSTANT: " << ParamDB->TURBULENT_VISCOSITY_CONSTANT << endl;
  OutFile << "TURBULENT_VISCOSITY_POWER: " << ParamDB->TURBULENT_VISCOSITY_POWER << endl;
  OutFile << "TURBULENT_VISCOSITY_SIGMA: " << ParamDB->TURBULENT_VISCOSITY_SIGMA << endl;
  OutFile << "TURBULENT_MOD_TYPE: " << ParamDB->TURBULENT_MOD_TYPE << endl;

  OutFile << "ARTIFICIAL_VISCOSITY_CONSTANT: " << ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT << endl;
  OutFile << "ARTIFICIAL_VISCOSITY_POWER: " << ParamDB->ARTIFICIAL_VISCOSITY_POWER << endl;

  OutFile << "FRICTION_CONSTANT: " << ParamDB->FRICTION_CONSTANT << endl;
  OutFile << "FRICTION_POWER: " << ParamDB->FRICTION_POWER << endl;
  OutFile << "FRICTION_TYPE: " << ParamDB->FRICTION_TYPE << endl;
  OutFile << "FRICTION_U0: " << ParamDB->FRICTION_U0 << endl;
  OutFile << "PENETRATION_CONSTANT: " << ParamDB->PENETRATION_CONSTANT << endl;
  OutFile << "PENETRATION_POWER: " << ParamDB->PENETRATION_POWER << endl;
  OutFile << "DIV_DIV_STAB_TYPE: " << ParamDB->DIV_DIV_STAB_TYPE << endl; 
  OutFile << "DIV_DIV_STAB_C1: " << ParamDB->DIV_DIV_STAB_C1 << endl; 
  OutFile << "DIV_DIV_STAB_C2: " << ParamDB->DIV_DIV_STAB_C2 << endl; 
  OutFile << "NSE_NONLINEAR_FORM: " << ParamDB->NSE_NONLINEAR_FORM << endl;
  OutFile << "OSEEN_ZERO_ORDER_COEFF: " << ParamDB->OSEEN_ZERO_ORDER_COEFF << endl;

  OutFile << "LP_FULL_GRADIENT: " << ParamDB->LP_FULL_GRADIENT << endl;
  OutFile << "LP_STREAMLINE: " << ParamDB->LP_STREAMLINE << endl;
  OutFile << "LP_DIVERGENCE: " << ParamDB->LP_DIVERGENCE << endl;
  OutFile << "LP_PRESSURE: " << ParamDB->LP_PRESSURE << endl;
  OutFile << "LP_COEFF_TYPE: " << ParamDB->LP_COEFF_TYPE << endl;
  OutFile << "LP_FULL_GRADIENT_COEFF: " << ParamDB->LP_FULL_GRADIENT_COEFF << endl;
  OutFile << "LP_STREAMLINE_COEFF: " << ParamDB->LP_STREAMLINE_COEFF << endl;
  OutFile << "LP_DIVERGENCE_COEFF: " << ParamDB->LP_DIVERGENCE_COEFF << endl;
  OutFile << "LP_PRESSURE_COEFF: " << ParamDB->LP_PRESSURE_COEFF << endl;
  OutFile << "LP_FULL_GRADIENT_EXPONENT: " << ParamDB->LP_FULL_GRADIENT_EXPONENT << endl;
  OutFile << "LP_STREAMLINE_EXPONENT: " << ParamDB->LP_STREAMLINE_EXPONENT << endl;
  OutFile << "LP_DIVERGENCE_EXPONENT: " << ParamDB->LP_DIVERGENCE_EXPONENT << endl;
  OutFile << "LP_PRESSURE_EXPONENT: " << ParamDB->LP_PRESSURE_EXPONENT << endl;
  OutFile << "LP_ORDER_DIFFERENCE: " << ParamDB->LP_ORDER_DIFFERENCE << endl;
  OutFile << "LP_FULL_GRADIENT_ORDER_DIFFERENCE: " << ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE << endl;
  OutFile << "LP_STREAMLINE_ORDER_DIFFERENCE: " << ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE << endl;
  OutFile << "LP_DIVERGENCE_ORDER_DIFFERENCE: " << ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE << endl;
  OutFile << "LP_PRESSURE_ORDER_DIFFERENCE: " << ParamDB->LP_PRESSURE_ORDER_DIFFERENCE << endl;

  OutFile << "SOLVE_ADJOINT_PROBLEM: " << ParamDB->SOLVE_ADJOINT_PROBLEM << endl;
  OutFile << "SOLD_ADJOINT: " << ParamDB->SOLD_ADJOINT << endl;
  OutFile << "N_STAGES_ADJOINT: " << ParamDB->N_STAGES_ADJOINT << endl;
  OutFile << "SC_NONLIN_ITE_ADJOINT: " << ParamDB->SC_NONLIN_ITE_ADJOINT << endl;
  OutFile << "OPTIMIZATION_ITE_TYPE_ADJOINT: " << ParamDB->OPTIMIZATION_ITE_TYPE_ADJOINT << endl;
  OutFile << "BFGS_VECTORS_ADJOINT: " << ParamDB->BFGS_VECTORS_ADJOINT << endl;
  OutFile << "RELATIVE_DECREASE_ADJOINT: " << ParamDB->RELATIVE_DECREASE_ADJOINT << endl;
  OutFile << "PENALTY_ADJOINT: " << ParamDB->PENALTY_ADJOINT << endl;
  OutFile << "PENALTY_VALUE_AT_ZERO_ADJOINT: " << ParamDB->PENALTY_VALUE_AT_ZERO_ADJOINT << endl;
  OutFile << "PENALTY_SMALLEST_PARAM_FAC_ADJOINT: " << ParamDB->PENALTY_SMALLEST_PARAM_FAC_ADJOINT << endl;
  OutFile << "PENALTY_LARGEST_PARAM_FAC_ADJOINT: " << ParamDB->PENALTY_LARGEST_PARAM_FAC_ADJOINT << endl;
  OutFile << "WEIGHT_RESIDUAL_L1_ADJOINT: " << ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT << endl;
  OutFile << "WEIGHT_RESIDUAL_L2_ADJOINT: " << ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT << endl;
  OutFile << "WEIGHT_GRADIENT_L1_ADJOINT: " << ParamDB->WEIGHT_GRADIENT_L1_ADJOINT << endl;
  OutFile << "WEIGHT_GRADIENT_L2_ADJOINT: " << ParamDB->WEIGHT_GRADIENT_L2_ADJOINT << endl;
  OutFile << "WEIGHT_STREAM_DER_L1_ADJOINT: " << ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT << endl;
  OutFile << "WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT: " << ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT << endl;
  OutFile << "WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT: " << ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT << endl;
  OutFile << "REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT: " << ParamDB->REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT << endl;
  OutFile << "WEIGHT_RESIDUAL_LP_ADJONT: " << ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT << endl;
  OutFile << "WEIGHT_RESIDUAL_EXP_LP_ADJONT: " << ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT << endl;
  OutFile << "WEIGHT_RESIDUAL_CW_ADJOINT: " << ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT << endl;
  OutFile << "RESIDUAL_LP_ADJONT: " << ParamDB->RESIDUAL_LP_ADJOINT << endl;
  OutFile << "MIN_VAL_ADJOINT: " <<ParamDB->MIN_VAL_ADJOINT <<endl;
  OutFile << "MAX_VAL_ADJOINT: " <<ParamDB->MAX_VAL_ADJOINT << endl;
  OutFile << "MIN_MAX_EXPONENT_ONE_ADJOINT: " <<ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT<< endl;
  OutFile << "MIN_MAX_EXPONENT_TWO_ADJOINT: " <<ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT<< endl;
  OutFile << "MIN_MAX_FACTOR_ONE_ADJOINT: " << ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT<< endl;
  OutFile << "MIN_MAX_FACTOR_TWO_ADJOINT: " << ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT<< endl;
  OutFile << "MIN_MAX_ADJOINT: " << ParamDB->MIN_MAX_ADJOINT<< endl;
  
  OutFile << "BASENAME: " << ParamDB->BASENAME << endl;
  OutFile << "VTKBASENAME: " << ParamDB->VTKBASENAME << endl;
  OutFile << "PSBASENAME: " << ParamDB->PSBASENAME << endl;
  OutFile << "GRAPEBASENAME: " << ParamDB->GRAPEBASENAME << endl;
  OutFile << "GNUBASENAME: " << ParamDB->GNUBASENAME << endl;
  OutFile << "READGRAPEBASENAME: " << ParamDB->READGRAPEBASENAME << endl;
  OutFile << "GMVBASENAME: " << ParamDB->GMVBASENAME << endl;
  OutFile << "MATLABBASENAME: " << ParamDB->MATLABBASENAME << endl;
  OutFile << "SAVE_DATA_FILENAME: " << ParamDB->SAVE_DATA_FILENAME << endl;
  OutFile << "READ_DATA_FILENAME: " << ParamDB->READ_DATA_FILENAME << endl;
  OutFile << "POD_FILENAME: " << ParamDB->POD_FILENAME << endl;
  OutFile << "SNAP_FILENAME: " << ParamDB->SNAP_FILENAME << endl;

  OutFile << "SOLVER_TYPE: " << ParamDB->SOLVER_TYPE << endl;
  OutFile << "WRITE_PS: " << ParamDB->WRITE_PS << endl;
  OutFile << "WRITE_GRAPE: " << ParamDB->WRITE_GRAPE << endl;
  OutFile << "WRITE_GMV: " << ParamDB->WRITE_GMV << endl;
  OutFile << "WRITE_AMIRA: " << ParamDB->WRITE_AMIRA << endl;
  OutFile << "WRITE_GNU: " << ParamDB->WRITE_GNU << endl;
  OutFile << "WRITE_AMIRA: " << ParamDB->WRITE_AMIRA << endl;
  OutFile << "WRITE_VTK: " << ParamDB->WRITE_VTK << endl;
  OutFile << "MESH_TYPE: " << ParamDB->MESH_TYPE << endl;  
  OutFile << "WRITE_SNAPSHOTS: " << ParamDB->WRITE_SNAPSHOTS << endl;
  OutFile << "WRITE_MATLAB_MATRIX: " << ParamDB->WRITE_MATLAB_MATRIX << endl; 
  OutFile << "WRITE_MATLAB: " << ParamDB->WRITE_MATLAB << endl;
  OutFile << "SAVE_DATA: " << ParamDB->SAVE_DATA << endl;
  OutFile << "READ_DATA: " << ParamDB->READ_DATA << endl;
  OutFile << "MEASURE_ERRORS: " << ParamDB->MEASURE_ERRORS << endl;
  OutFile << "ESTIMATE_ERRORS: " << ParamDB->ESTIMATE_ERRORS << endl;
  OutFile << "SOLVE_ADJOINT_PROBLEM: " << ParamDB->SOLVE_ADJOINT_PROBLEM << endl;
  OutFile << "COMPUTE_VORTICITY_DIVERGENCE: " << ParamDB->COMPUTE_VORTICITY_DIVERGENCE << endl;
  OutFile << "READ_GRAPE_FILE: " << ParamDB->READ_GRAPE_FILE  << endl;
  OutFile << "BUILD_PODFILE:" << ParamDB->BUILD_PODFILE << endl;
  OutFile << "POD_FLUCT_FIELD:" << ParamDB->POD_FLUCT_FIELD << endl;
  OutFile << "POD_FLUCT_FIELD_P:" << ParamDB->POD_FLUCT_FIELD_P << endl;
  OutFile << "PROJECTION_METHOD:" << ParamDB->PROJECTION_METHOD << endl;

  OutFile << "P0: " << ParamDB->P0 << endl;
  OutFile << "P1: " << ParamDB->P1 << endl;
  OutFile << "P2: " << ParamDB->P2 << endl;
  OutFile << "P3: " << ParamDB->P3 << endl;
  OutFile << "P4: " << ParamDB->P4 << endl;
  OutFile << "P5: " << ParamDB->P5 << endl;
  OutFile << "P6: " << ParamDB->P6 << endl;
  OutFile << "P7: " << ParamDB->P7 << endl;
  OutFile << "P8: " << ParamDB->P8 << endl;
  OutFile << "P9: " << ParamDB->P9 << endl;
  OutFile << "P10: " << ParamDB->P10 << endl;
  OutFile << "P11: " << ParamDB->P11 << endl;
  OutFile << "P12: " << ParamDB->P12 << endl;
  OutFile << "P13: " << ParamDB->P13 << endl;
  OutFile << "P14: " << ParamDB->P14 << endl;
  OutFile << "P15: " << ParamDB->P15 << endl;

  OutFile << "PBE_P0: " << ParamDB->PBE_P0 << endl;
  OutFile << "PBE_P1: " << ParamDB->PBE_P1 << endl;
  OutFile << "PBE_P2: " << ParamDB->PBE_P2 << endl;
  OutFile << "PBE_P3: " << ParamDB->PBE_P3 << endl;
  OutFile << "PBE_P4: " << ParamDB->PBE_P4 << endl;
  OutFile << "PBE_P5: " << ParamDB->PBE_P5 << endl;
  OutFile << "PBE_P6: " << ParamDB->PBE_P6 << endl;
  OutFile << "PBE_P7: " << ParamDB->PBE_P7 << endl;
  OutFile << "PBE_P8: " << ParamDB->PBE_P8 << endl;
  OutFile << "PBE_P9: " << ParamDB->PBE_P9 << endl;

  OutFile << "DG_P0: " << ParamDB->DG_P0 << endl;
  OutFile << "DG_P1: " << ParamDB->DG_P1 << endl;
  OutFile << "DG_P2: " << ParamDB->DG_P2 << endl;
  OutFile << "DG_P3: " << ParamDB->DG_P3 << endl;
  OutFile << "DG_P4: " << ParamDB->DG_P4 << endl;
  OutFile << "DG_P5: " << ParamDB->DG_P5 << endl;
  OutFile << "DG_P6: " << ParamDB->DG_P6 << endl;
  OutFile << "DG_P7: " << ParamDB->DG_P7 << endl;
  OutFile << "DG_P8: " << ParamDB->DG_P8 << endl;
  OutFile << "DG_P9: " << ParamDB->DG_P9 << endl;

  
  OutFile << "VORTICITY THICKNESS FOR MIXING LAYER (P8): " << ParamDB->P8 << endl;

  OutFile << "*********** PARAMETERS FOR SCALAR SOLVER ***********" << endl;
  OutFile << "SC_NONLIN_ITE_TYPE_SCALAR: " << ParamDB->SC_NONLIN_ITE_TYPE_SCALAR << endl;
  OutFile << "SC_NONLIN_MAXIT_SCALAR: " << ParamDB->SC_NONLIN_MAXIT_SCALAR << endl;
  OutFile << "SC_NONLIN_RES_NORM_MIN_SCALAR: " << ParamDB->SC_NONLIN_RES_NORM_MIN_SCALAR << endl;
  OutFile << "SC_NONLIN_DAMP_FACTOR_SCALAR: " << ParamDB->SC_NONLIN_DAMP_FACTOR_SCALAR << endl;

  OutFile << "SC_SOLVER_SCALAR: " << ParamDB->SC_SOLVER_SCALAR << endl;
  OutFile << "SC_PRECONDITIONER_SCALAR: " << ParamDB->SC_PRECONDITIONER_SCALAR << endl;
  OutFile << "SC_LIN_RED_FACTOR_SCALAR: " << ParamDB->SC_LIN_RED_FACTOR_SCALAR << endl;
  OutFile << "SC_LIN_RES_NORM_MIN_SCALAR: " << ParamDB->SC_LIN_RES_NORM_MIN_SCALAR << endl;
  OutFile << "SC_LIN_MAXIT_SCALAR: " << ParamDB->SC_LIN_MAXIT_SCALAR << endl;
  OutFile << "SC_LIN_RED_FACTOR_SCALAR_SOLD: " << ParamDB->SC_LIN_RED_FACTOR_SCALAR_SOLD << endl;
  OutFile << "SC_LIN_RES_NORM_MIN_SCALAR_SOLD: " << ParamDB->SC_LIN_RES_NORM_MIN_SCALAR_SOLD << endl;
  OutFile << "SC_LIN_MAXIT_SCALAR_SOLD: " << ParamDB->SC_LIN_MAXIT_SCALAR_SOLD << endl;
  OutFile << "SC_FLEXIBLE_KRYLOV_SPACE_SOLVER: " << ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER << endl;
  OutFile << "SC_NONLIN_ITE_ADJOINT: " << ParamDB->SC_NONLIN_ITE_ADJOINT << endl;

  OutFile << "SC_MG_TYPE_SCALAR: " << ParamDB->SC_MG_TYPE_SCALAR << endl; 
  OutFile << "SC_MG_CYCLE_SCALAR: " << ParamDB->SC_MG_CYCLE_SCALAR << endl; 
  OutFile << "SC_SMOOTHER_SCALAR: " << ParamDB->SC_SMOOTHER_SCALAR << endl;
  OutFile << "SC_PRE_SMOOTH_SCALAR: " << ParamDB->SC_PRE_SMOOTH_SCALAR << endl;
  OutFile << "SC_POST_SMOOTH_SCALAR: " << ParamDB->SC_POST_SMOOTH_SCALAR << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_SCALAR: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SCALAR << endl;
  OutFile << "SC_COARSE_SMOOTHER_SCALAR: " << ParamDB->SC_COARSE_SMOOTHER_SCALAR << endl;
  OutFile << "SC_COARSE_MAXIT_SCALAR: " << ParamDB->SC_COARSE_MAXIT_SCALAR << endl;
  OutFile << "SC_COARSE_RED_FACTOR_SCALAR: " << ParamDB->SC_COARSE_RED_FACTOR_SCALAR << endl;

  OutFile << "SC_GMG_DAMP_FACTOR_SCALAR: " << ParamDB->SC_GMG_DAMP_FACTOR_SCALAR << endl;
  OutFile << "SC_GMG_DAMP_FACTOR_FINE_SCALAR: " << ParamDB->SC_GMG_DAMP_FACTOR_FINE_SCALAR << endl;

  OutFile << "SC_COARSEST_LEVEL_SCALAR: " << ParamDB->SC_COARSEST_LEVEL_SCALAR << endl;
  OutFile << "SC_FIRST_SOLUTION_LEVEL_SCALAR: " << ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR << endl;
  
  OutFile << "SC_STEP_LENGTH_CONTROL_FINE_SCALAR: " << ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR << endl;
  OutFile << "SC_STEP_LENGTH_CONTROL_ALL_SCALAR: " << ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR << endl;

  OutFile << "*********** PARAMETERS FOR SADDLE POINT SOLVER ***********" << endl;
  OutFile << "SC_NONLIN_ITE_TYPE_SADDLE: " << ParamDB->SC_NONLIN_ITE_TYPE_SADDLE << endl;
  OutFile << "SC_NONLIN_MAXIT_SADDLE: " << ParamDB->SC_NONLIN_MAXIT_SADDLE << endl;
  OutFile << "SC_NONLIN_RES_NORM_MIN_SADDLE: " << ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE << endl;
  OutFile << "SC_NONLIN_DAMP_FACTOR_SADDLE: " << ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE << endl;
  OutFile << "SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE: " << ParamDB->SC_NONLIN_RES_NORM_MIN_SCALE_SADDLE << endl;

  OutFile << "SC_SOLVER_SADDLE: " << ParamDB->SC_SOLVER_SADDLE << endl;
  OutFile << "SC_PRECONDITIONER_SADDLE: " << ParamDB->SC_PRECONDITIONER_SADDLE << endl;
  OutFile << "SC_LIN_RED_FACTOR_SADDLE: " << ParamDB->SC_LIN_RED_FACTOR_SADDLE << endl;
  OutFile << "SC_LIN_RES_NORM_MIN_SADDLE: " << ParamDB->SC_LIN_RES_NORM_MIN_SADDLE << endl;
  OutFile << "SC_LIN_MAXIT_SADDLE: " << ParamDB->SC_LIN_MAXIT_SADDLE << endl;

  OutFile << "SC_MG_TYPE_SADDLE: " << ParamDB->SC_MG_TYPE_SADDLE << endl; 
  OutFile << "SC_MG_CYCLE_SADDLE: " << ParamDB->SC_MG_CYCLE_SADDLE << endl; 
  OutFile << "SC_SMOOTHER_SADDLE: " << ParamDB->SC_SMOOTHER_SADDLE << endl;
  OutFile << "SC_PRE_SMOOTH_SADDLE: " << ParamDB->SC_PRE_SMOOTH_SADDLE << endl;
  OutFile << "SC_POST_SMOOTH_SADDLE: " << ParamDB->SC_POST_SMOOTH_SADDLE << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_SADDLE: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_SADDLE << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SADDLE << endl;
  OutFile << "SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE: " << ParamDB->SC_SMOOTH_DAMP_FACTOR_COARSE_SADDLE << endl;
  OutFile << "SC_COARSE_SMOOTHER_SADDLE: " << ParamDB->SC_COARSE_SMOOTHER_SADDLE << endl;
  OutFile << "SC_COARSE_MAXIT_SADDLE: " << ParamDB->SC_COARSE_MAXIT_SADDLE << endl;
  OutFile << "SC_COARSE_RED_FACTOR_SADDLE: " << ParamDB->SC_COARSE_RED_FACTOR_SADDLE << endl;
  OutFile << "SC_GMG_DAMP_FACTOR_SADDLE: " << ParamDB->SC_GMG_DAMP_FACTOR_SADDLE << endl;
  OutFile << "SC_GMG_DAMP_FACTOR_FINE_SADDLE: " << ParamDB->SC_GMG_DAMP_FACTOR_FINE_SADDLE << endl;
  
  OutFile << "SC_COARSEST_LEVEL_SADDLE: " << ParamDB->SC_COARSEST_LEVEL_SADDLE << endl;
  OutFile << "SC_FIRST_SOLUTION_LEVEL_SADDLE: " << ParamDB->SC_FIRST_SOLUTION_LEVEL_SADDLE << endl;

  OutFile << "SC_STEP_LENGTH_CONTROL_FINE_SADDLE: " << ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SADDLE << endl;
  OutFile << "SC_STEP_LENGTH_CONTROL_ALL_SADDLE: " << ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SADDLE << endl;
  OutFile << "SC_LARGEST_DIRECT_SOLVE: " << ParamDB->SC_LARGEST_DIRECT_SOLVE << endl;
  OutFile << "SC_DOWNWIND_TYPE: " << ParamDB->SC_DOWNWIND_TYPE << endl;

  OutFile << "CC_ALPHA: " << ParamDB->CC_ALPHA << endl;
  OutFile << "CC_BETA: " << ParamDB->CC_BETA << endl;
  OutFile << "CC_MINCLUSTER: " << ParamDB->CC_MINCLUSTER << endl;
  OutFile << "CC_MAXCLUSTER: " << ParamDB->CC_MAXCLUSTER << endl;
  OutFile << "CC_MAXDISTANCE: " << ParamDB->CC_MAXDISTANCE << endl;
  OutFile << "CC_MAXCONNECTIVITY: " << ParamDB->CC_MAXCONNECTIVITY << endl;
  OutFile << "CC_DEPTHTARGET: " << ParamDB->CC_DEPTHTARGET << endl;
  OutFile << "CC_COARSENTARGET: " << ParamDB->CC_COARSENTARGET << endl;
  OutFile << "CC_COARSENRATE: " << ParamDB->CC_COARSENRATE << endl;
  OutFile << "CC_MAJOR: " << ParamDB->CC_MAJOR << endl;
  OutFile << "CC_DEPENDENCY: " << ParamDB->CC_DEPENDENCY << endl;
  OutFile << "CC_RESCALE: " << ParamDB->CC_RESCALE << endl;
  OutFile << "CC_VERBOSE: " << ParamDB->CC_VERBOSE << endl;
  
  OutFile << "SC_SYSTEM_TYPE: " << ParamDB->SC_SYSTEM_TYPE << endl;
  OutFile << "SC_AMG_PREC_IT: " << ParamDB->SC_AMG_PREC_IT << endl;
  OutFile << "SC_AMG_PREC_RED_FACTOR: " << ParamDB->SC_AMG_PREC_RED_FACTOR << endl;
  OutFile << "SC_SMOOTHER_RED_FACTOR: " << ParamDB->SC_SMOOTHER_RED_FACTOR << endl;
  OutFile << "SC_EX_MAXIT: " << ParamDB->SC_EX_MAXIT << endl;
  OutFile << "SC_GMRES_RESTART: " << ParamDB->SC_GMRES_RESTART << endl;
  OutFile << "SC_LCD_START_VECTOR: " << ParamDB->SC_LCD_START_VECTOR << endl;
  OutFile << "SC_ILU_BETA: " << ParamDB->SC_ILU_BETA << endl;
  OutFile << "SC_SOR_OMEGA: " << ParamDB->SC_SOR_OMEGA << endl;
  OutFile << "SC_OMEGA_COARSE_0: " << ParamDB->SC_OMEGA_COARSE_0 << endl;
  OutFile << "SC_OMEGA_P_0: " << ParamDB->SC_OMEGA_P_0 << endl;
  OutFile << "SC_ILUT_TOL: " << ParamDB->SC_ILUT_TOL << endl;
  OutFile << "SC_ILUT_ABSOLUTE_FILLIN: " << ParamDB->SC_ILUT_ABSOLUTE_FILLIN << endl;
  OutFile << "SC_ILUT_RELATIVE_FILLIN: " << ParamDB->SC_ILUT_RELATIVE_FILLIN << endl;
  OutFile << "SC_ILUT_SORT: " << ParamDB->SC_ILUT_SORT << endl;
  OutFile << "SC_SCHUR_INV_OF_A: " << ParamDB->SC_SCHUR_INV_OF_A << endl;
  OutFile << "SC_SCHUR_INV_OF_A_MAXIT: " << ParamDB->SC_SCHUR_INV_OF_A_MAXIT << endl;
  OutFile << "SC_SCHUR_ITERATION_DAMP: " << ParamDB->SC_SCHUR_ITERATION_DAMP << endl;
  OutFile << "SC_SCHUR_ITERATION_MAXIT: " << ParamDB->SC_SCHUR_ITERATION_MAXIT << endl;
  OutFile << "SC_SCHUR_STEP_LENGTH_CONTROL: " << ParamDB->SC_SCHUR_STEP_LENGTH_CONTROL << endl;
  OutFile << "SC_MIXED_BCGS_CGS_SWITCH_TOL: " << ParamDB->SC_MIXED_BCGS_CGS_SWITCH_TOL << endl;
  OutFile << "SC_DIV_FACTOR: " << ParamDB->SC_DIV_FACTOR << endl;
  OutFile << "SC_NONLIN_DIV_FACTOR: " << ParamDB->SC_NONLIN_DIV_FACTOR << endl;
  OutFile << "SC_SMOOTHING_STEPS: " << ParamDB->SC_SMOOTHING_STEPS << endl;
  OutFile << "SC_N1_PARAM: " << ParamDB->SC_N1_PARAM << endl;
  OutFile << "SC_N2_PARAM: " << ParamDB->SC_N2_PARAM << endl;
  OutFile << "SC_MINIT: " << ParamDB->SC_MINIT << endl;
  OutFile << "SC_VAS_LAZ_DELTA: " << ParamDB->SC_VAS_LAZ_DELTA << endl;
  OutFile << "SC_ROW_EQUILIBRATION: " << ParamDB->SC_ROW_EQUILIBRATION << endl;
  OutFile << "SC_BRAESS_SARAZIN_MATRIX: " << ParamDB->SC_BRAESS_SARAZIN_MATRIX << endl;
  OutFile << "SC_BRAESS_SARAZIN_ALPHA: " << ParamDB->SC_BRAESS_SARAZIN_ALPHA << endl;
  OutFile << "SC_VERBOSE: " << ParamDB->SC_VERBOSE << endl;
  OutFile << "SC_VERBOSE_AMG: " << ParamDB->SC_VERBOSE_AMG << endl;

  OutFile << "CHAR_L0: " << ParamDB->CHAR_L0 << endl;
  OutFile << "D_VISCOSITY: " << ParamDB->D_VISCOSITY << endl;
  OutFile << "SURF_TENSION: " << ParamDB->SURF_TENSION << endl;
  OutFile << "IMPACT_ANGLE: " << ParamDB->IMPACT_ANGLE << endl;
  OutFile << "Area: " << ParamDB->Area << endl;

  // ******** parameters for VMS *********//
  OutFile << "VMS_LARGE_VELOCITY_SPACE: " <<  ParamDB->VMS_LARGE_VELOCITY_SPACE<< endl;
  OutFile << "VMS_COARSE_MG_SMAGO: " <<  ParamDB->VMS_COARSE_MG_SMAGO<< endl;
  // constants in AdaptProjectionSpace 
  OutFile << "VMS_ADAPT_LOWER: " << ParamDB->VMS_ADAPT_LOWER << endl; 
  OutFile << "VMS_ADAPT_MIDDLE: " << ParamDB->VMS_ADAPT_MIDDLE << endl; 
  OutFile << "VMS_ADAPT_UPPER: " << ParamDB->VMS_ADAPT_UPPER << endl; 
  OutFile << "VMS_ADAPT_STEPS: " << ParamDB->VMS_ADAPT_STEPS << endl; 
  OutFile << "VMS_ADAPT_COMP: " << ParamDB->VMS_ADAPT_COMP << endl; 

  OutFile << "SUPERCONVERGENCE_ORDER: " << ParamDB->SUPERCONVERGENCE_ORDER << endl;
  OutFile << "FEM_FCT_LINEAR_TYPE: " << ParamDB->FEM_FCT_LINEAR_TYPE << endl;
  OutFile << "FEM_FCT_PRELIMITING: " << ParamDB->FEM_FCT_PRELIMITING << endl;
  OutFile << "FEM_FCT_GROUP_FEM: " << ParamDB->FEM_FCT_GROUP_FEM << endl;
  OutFile << "GROUP_FEM: " << ParamDB->GROUP_FEM << endl;
  OutFile << "WENO_TYPE: " << ParamDB->WENO_TYPE << endl;
  
  OutFile << "WRITE_SNAPSHOTS: " << ParamDB->WRITE_SNAPSHOTS << endl;
  OutFile << "DO_ROM: " << ParamDB->DO_ROM << endl;
  OutFile << "DO_ROM_P: " << ParamDB->DO_ROM_P << endl;
  OutFile << "RANK_OF_BASIS: " << ParamDB->RANK_OF_BASIS << endl;
  OutFile << "RANK_OF_BASIS_P: " << ParamDB->RANK_OF_BASIS_P << endl;
  OutFile << "POD_INNER_PRODUCT: " << ParamDB->POD_INNER_PRODUCT << endl;
  OutFile << "POD_INNER_PRODUCT_P: " << ParamDB->POD_INNER_PRODUCT_P << endl;
  OutFile << "BUILD_PODFILE: " << ParamDB->BUILD_PODFILE << endl;
  OutFile << "POD_FILENAME: " << ParamDB->POD_FILENAME << endl;
  OutFile << "SNAP_FILENAME: " << ParamDB->SNAP_FILENAME << endl;
  OutFile << "POD_FLUCT_FIELD: " << ParamDB->POD_FLUCT_FIELD << endl;
  OutFile << "POD_FLUCT_FIELD_P: " << ParamDB->POD_FLUCT_FIELD_P << endl;
  OutFile << "P_ROM_METHOD: " << ParamDB->P_ROM_METHOD << endl;
  OutFile << "PROJECTION_METHOD: " << ParamDB->PROJECTION_METHOD << endl;

  /** write parameter for non-conforming elements */
  OutFile << "NC_TYPE: " << ParamDB->NC_TYPE << endl;

  OutFile << "CHANNEL_STATISTICS2_WITH_MODEL: " << ParamDB->CHANNEL_STATISTICS2_WITH_MODEL <<endl;
  OutFile << "BULK_REACTION_DISC: " << ParamDB->BULK_REACTION_DISC <<endl;
  OutFile << "BULK_PB_DISC: " << ParamDB->BULK_PB_DISC <<endl;
  OutFile << "BULK_PB_DISC_STAB: " << ParamDB->BULK_PB_DISC_STAB <<endl;
  OutFile << "BULK_COUPLING: " << ParamDB->BULK_COUPLING <<endl;
  OutFile << "BULK_GROWTH_RATE: " << ParamDB->BULK_GROWTH_RATE <<endl;
  OutFile << "BULK_REACTION_MASS_LUMPING: " << ParamDB->BULK_REACTION_MASS_LUMPING <<endl;
  OutFile << "BULK_REACTION_C_CUT: " << ParamDB->BULK_REACTION_C_CUT <<endl;
  OutFile << "BULK_METHODS_OF_MOMENTS: " << ParamDB->BULK_METHODS_OF_MOMENTS <<endl;
  OutFile << "BULK_MOM_DISC: " << ParamDB->BULK_MOM_DISC <<endl;
  OutFile << "BULK_SOLD_PARAMETER_TYPE: " << ParamDB->BULK_SOLD_PARAMETER_TYPE <<endl;
  OutFile << "N_CELL_LAYERS_PSD: " << ParamDB->N_CELL_LAYERS_PSD <<endl;
  OutFile << "N_CELL_LAYERS_PSD_2: " << ParamDB->N_CELL_LAYERS_PSD <<endl;
  OutFile << "OUTPUT_NODE_LAYER_PSD: " << ParamDB->OUTPUT_NODE_LAYER_PSD <<endl;
  OutFile << "BULK_l_infty: " << ParamDB->BULK_l_infty << endl;
  OutFile << "BULK_u_infty: " << ParamDB->BULK_u_infty << endl;
  OutFile << "BULK_c_infty: " << ParamDB->BULK_c_infty << endl;
  OutFile << "BULK_c_C_infty_sat: " << ParamDB->BULK_c_C_infty_sat << endl;
  OutFile << "BULK_C_g: " << ParamDB->BULK_C_g << endl;
  OutFile << "BULK_C_nuc: " << ParamDB->BULK_C_nuc << endl;
  OutFile << "BULK_C_sat: " << ParamDB->BULK_C_sat << endl;
  OutFile << "BULK_C_2: " << ParamDB->BULK_C_2 << endl;
  OutFile << "BULK_D_A: " << ParamDB->BULK_D_A << endl;
  OutFile << "BULK_D_P_0: " << ParamDB->BULK_D_P_0 << endl;
  OutFile << "BULK_D_P_MAX: " << ParamDB->BULK_D_P_MAX << endl;
  OutFile << "BULK_k_g: " << ParamDB->BULK_k_g << endl;
  OutFile << "BULK_k_r: " << ParamDB->BULK_k_r << endl;
  OutFile << "BULK_k_nuc: " << ParamDB->BULK_k_nuc << endl;
  OutFile << "SSMUM_MP_X: " << ParamDB->SSMUM_MP_X << endl;
  OutFile << "SSMUM_MP_Y: " << ParamDB->SSMUM_MP_Y << endl;
  OutFile << "SSMUM_INNER_RADIUS: " << ParamDB->SSMUM_INNER_RADIUS << endl;
  OutFile << "SSMUM_OUTER_RADIUS: " << ParamDB->SSMUM_OUTER_RADIUS << endl;
  OutFile << "SSMUM_ROT_PER_SECOND: " << ParamDB->SSMUM_ROT_PER_SECOND << endl;
  OutFile << "SSMUM_MAX_CELLS_LAYERS: " << ParamDB->SSMUM_MAX_CELLS_LAYERS << endl;
  OutFile << "SSMUM_INTERPOLATION: " << ParamDB->SSMUM_INTERPOLATION << endl;
  
  OutFile << "INPUT_QUAD_RULE: " << ParamDB->INPUT_QUAD_RULE << endl;
  
  /** Parameters for Stokes--Darcy */
  OutFile << "StoDa_interfaceType: " << ParamDB->StoDa_interfaceType << endl;
  OutFile << "StoDa_alpha: " << ParamDB->StoDa_alpha << endl;
  OutFile << "StoDa_problemType: " << ParamDB->StoDa_problemType << endl;
  OutFile << "StoDa_updatingStrategy: " << ParamDB->StoDa_updatingStrategy << endl;
  OutFile << "StoDa_theta_f: " << ParamDB->StoDa_theta_f << endl;
  OutFile << "StoDa_theta_p: " << ParamDB->StoDa_theta_p << endl;
  OutFile << "StoDa_gamma_f: " << ParamDB->StoDa_gamma_f << endl;
  OutFile << "StoDa_gamma_p: " << ParamDB->StoDa_gamma_p << endl;
  OutFile << "StoDa_weakGamma: " << ParamDB->StoDa_weakGamma << endl;
  OutFile << "StoDa_solutionStrategy: " << ParamDB->StoDa_solutionStrategy << endl;
  OutFile << "StoDa_StokesFirst: " << ParamDB->StoDa_StokesFirst << endl;
  OutFile << "StoDa_algorithm: " << ParamDB->StoDa_algorithm << endl;
  OutFile << "StoDa_relDiff_interfaceError: " << ParamDB->StoDa_relDiff_interfaceError << endl;
  OutFile << "StoDa_relDiff_factor1: " << ParamDB->StoDa_relDiff_factor1 << endl;
  OutFile << "StoDa_relDiff_factor2: " << ParamDB->StoDa_relDiff_factor2 << endl;
  OutFile << "StoDa_relDiff_factor3: " << ParamDB->StoDa_relDiff_factor3 << endl;
  OutFile << "StoDa_relDiff_solution: " << ParamDB->StoDa_relDiff_solution << endl;
  OutFile << "StoDa_bigResidual: " << ParamDB->StoDa_bigResidual << endl;
  OutFile << "StoDa_periodicBoundary: " << ParamDB->StoDa_periodicBoundary << endl;
  OutFile << "StoDa_periodicBoundaryPressureDrop: " << ParamDB->StoDa_periodicBoundaryPressureDrop << endl;
  OutFile << "StoDa_nIterations: " << ParamDB->StoDa_nIterations << endl;

 
  OutFile << "HEAT_TANGENTIAL_STRESS_FACTOR: " << ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR << endl;
  OutFile << "HEAT_SOLID_SURFACE_FACTOR: " << ParamDB->HEAT_SOLID_SURFACE_FACTOR << endl;
  OutFile << "EQ_CONTACT_ANGLE: " << ParamDB->EQ_CONTACT_ANGLE << endl;
  OutFile << "AD_CONTACT_ANGLE: " << ParamDB->AD_CONTACT_ANGLE << endl;
  OutFile << "RE_CONTACT_ANGLE: " << ParamDB->RE_CONTACT_ANGLE << endl;
  OutFile << "DY_CONTACT_ANGLE: " << ParamDB->DY_CONTACT_ANGLE << endl;  
  OutFile << "CONTACT_ANGLE_TYPE: " << ParamDB->CONTACT_ANGLE_TYPE << endl; 
  
}

void TDatabase::WriteTimeDB()
{
  OutFile << "CURRENTTIME: " << TimeDB->CURRENTTIME << endl;
  OutFile << "CURRENTTIMESTEPLENGTH: " << TimeDB->CURRENTTIMESTEPLENGTH << endl;
  OutFile << "TIMESTEPLENGTH: " << TimeDB->TIMESTEPLENGTH << endl;
  OutFile << "MIN_TIMESTEPLENGTH: " << TimeDB->MIN_TIMESTEPLENGTH << endl;
  OutFile << "MAX_TIMESTEPLENGTH: " << TimeDB->MAX_TIMESTEPLENGTH << endl;
  OutFile << "TIMESTEPLENGTH_TOL: " << TimeDB->TIMESTEPLENGTH_TOL << endl;
  OutFile << "TIMESTEPLENGTH_CONTROL: " << TimeDB->TIMESTEPLENGTH_CONTROL << endl;
  OutFile << "TIMESTEPLENGTH_CONTROLLER: " <<  TimeDB->TIMESTEPLENGTH_CONTROLLER <<endl;    
  OutFile << "TIMESTEPLENGTH_PARA_KK_I: " <<  TimeDB->TIMESTEPLENGTH_PARA_KK_I <<endl;
  OutFile << "TIMESTEPLENGTH_PARA_KK_P: " <<  TimeDB->TIMESTEPLENGTH_PARA_KK_P << endl;
  OutFile << "TIMESTEPLENGTH_PARA_KK_E: " <<  TimeDB->TIMESTEPLENGTH_PARA_KK_E << endl;
  OutFile << "TIMESTEPLENGTH_PARA_KK_R: " <<  TimeDB->TIMESTEPLENGTH_PARA_KK_R << endl;
  OutFile << "TIMESTEPLENGTH_PARA_KK_D: " << TimeDB->TIMESTEPLENGTH_PARA_KK_D << endl;
  OutFile << "TIMESTEPLENGTH_PARA_FAC: " <<  TimeDB->TIMESTEPLENGTH_PARA_FAC << endl;
  OutFile << "TIMESTEPLENGTH_PARA_FAC_MAX: " <<  TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX<< endl;
  OutFile << "TIMESTEPLENGTH_PARA_FAC_MIN: " <<  TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN << endl;
  OutFile << "TIMESTEPLENGTH_PARA_TOL: " <<  TimeDB->TIMESTEPLENGTH_PARA_TOL << endl;
  OutFile << "TIMESTEPLENGTH_PARA_ATOL: " <<  TimeDB->TIMESTEPLENGTH_PARA_ATOL << endl;
  OutFile << "TIMESTEPLENGTH_PARA_RTOL: " <<  TimeDB->TIMESTEPLENGTH_PARA_RTOL << endl;
  OutFile << "RESET_CURRENTTIME: " << TimeDB->RESET_CURRENTTIME << endl;
  OutFile << "RESET_CURRENTTIME_STARTTIME: " << TimeDB->RESET_CURRENTTIME_STARTTIME << endl;
  OutFile << "STEADY_STATE_TOL: "<< TimeDB->STEADY_STATE_TOL << endl;
  OutFile << "SCALE_DIVERGENCE_CONSTRAINT: "<< TimeDB->SCALE_DIVERGENCE_CONSTRAINT << endl;

  OutFile << "CONTROL: "<< TimeDB->CONTROL << endl;
  OutFile << "CONTROL_ALPHA: "<< TimeDB->CONTROL_ALPHA << endl;
  OutFile << "CONTROL_BETA: "<< TimeDB->CONTROL_BETA << endl;
  OutFile << "CONTROL_GAMMA: "<< TimeDB->CONTROL_GAMMA << endl;
  OutFile << "CONTROL_SAFTY: "<< TimeDB->CONTROL_SAFTY << endl;
  OutFile << "CONTROL_MAXSCALE: "<< TimeDB->CONTROL_MAXSCALE << endl;
  OutFile << "CONTROL_MINSCALE: "<< TimeDB->CONTROL_MINSCALE << endl;
  
  OutFile << "THETA1: " << TimeDB->THETA1 << endl;
  OutFile << "THETA2: " << TimeDB->THETA2 << endl;
  OutFile << "THETA3: " << TimeDB->THETA3 << endl;
  OutFile << "THETA4: " << TimeDB->THETA4 << endl;

  OutFile << "TIME_DISC: " << TimeDB->TIME_DISC << endl;
  OutFile << "TIME_DISC2: " << TimeDB->TIME_DISC2 << endl;

  OutFile << "STARTTIME: " << TimeDB->STARTTIME << endl;
  OutFile << "ENDTIME: " << TimeDB->ENDTIME << endl;

  OutFile << "T0: " << TimeDB->T0 << endl;
  OutFile << "T1: " << TimeDB->T1 << endl;
  OutFile << "T2: " << TimeDB->T2 << endl;
  OutFile << "T3: " << TimeDB->T3 << endl;
  OutFile << "T4: " << TimeDB->T4 << endl;
  OutFile << "T5: " << TimeDB->T5 << endl;
  OutFile << "T6: " << TimeDB->T6 << endl;
  OutFile << "T7: " << TimeDB->T7 << endl;
  OutFile << "T8: " << TimeDB->T8 << endl;
  OutFile << "T9: " << TimeDB->T9 << endl;

  OutFile << "STEPS_PER_IMAGE: " << TimeDB->STEPS_PER_IMAGE << endl;
  OutFile << "STEPS_PER_SNAP: " << TimeDB->STEPS_PER_SNAP << endl;

  OutFile << "RB_TYPE: " << TimeDB->RB_TYPE << endl;
  OutFile << "RB_TYPE2: " << TimeDB->RB_TYPE2 << endl;

  OutFile << "EXTRAPOLATE_VELOCITY: " << TimeDB->EXTRAPOLATE_VELOCITY << endl;
  OutFile << "EXTRAPOLATE_PRESSURE: " << TimeDB->EXTRAPOLATE_PRESSURE << endl;
  OutFile << "EXTRAPOLATE_STEPS: " << TimeDB->EXTRAPOLATE_STEPS << endl;
  OutFile << "EXTRAPOLATE_WEIGHT: " << TimeDB->EXTRAPOLATE_WEIGHT << endl;
} 

void TDatabase::CheckParameterConsistencyNSE()
{
  // Newton method
  if ((ParamDB->SC_NONLIN_ITE_TYPE_SADDLE)&&(ParamDB->NSTYPE<=2))
  {
    ParamDB->NSTYPE+=2;
    OutPut("NSTYPE changed to " << ParamDB->NSTYPE);
    OutPut(" because of SC_NONLIN_ITE_TYPE_SADDLE  = " <<ParamDB->SC_NONLIN_ITE_TYPE_SADDLE << endl);
  }

  if(ParamDB->PRESSURE_SPACE == -4711 || ParamDB->PRESSURE_SPACE>0)
  {
    // continuous pressure and cell Vanka do not work
    if (((ParamDB->VELOCITY_SPACE==2) || (ParamDB->VELOCITY_SPACE==3)
	 ||(ParamDB->VELOCITY_SPACE==101))
        &&(ParamDB->SC_SMOOTHER_SADDLE<3))
    {
      ParamDB->SC_SMOOTHER_SADDLE+=2;
      OutPut("SC_SMOOTHER_SADDLE changed to " << ParamDB->SC_SMOOTHER_SADDLE);
      OutPut(" because of continuous pressure"<< endl);
    }

    if (((ParamDB->VELOCITY_SPACE==2) || (ParamDB->VELOCITY_SPACE==3)||
	 (ParamDB->VELOCITY_SPACE==101))  
        &&(ParamDB->SC_COARSE_SMOOTHER_SADDLE<3))
    {
      ParamDB->SC_COARSE_SMOOTHER_SADDLE+=2;
      OutPut("SC_COARSE_SMOOTHER_SADDLE changed to " << ParamDB->SC_COARSE_SMOOTHER_SADDLE);
      OutPut(" because of continuous pressure"<< endl);
    }
  }
  if (ParamDB->GROUP_FEM)
  {
    if (ParamDB->DISCTYPE != GALERKIN)
    {
      ParamDB->DISCTYPE = GALERKIN;
      OutPut("GROUP_FEM: changed DISCTYPE to " << ParamDB->DISCTYPE << endl);
    }
    if (ParamDB->NSTYPE != 1)
    {
      //ParamDB->NSTYPE = 1;
      //OutPut("GROUP_FEM: changed NSTYPE to " << ParamDB->NSTYPE << endl);
      OutPut("WARNING: GROUP_FEM works properly only with NSTYPE = 1" << endl);
    }
    if (ParamDB->SC_MG_TYPE_SADDLE != 0)
    {
      ParamDB->SC_MG_TYPE_SADDLE = 0;
      OutPut("GROUP_FEM: changed SC_MG_TYPE_SADDLE to " << ParamDB->SC_MG_TYPE_SADDLE << endl);
    }
  }


  if ((ParamDB->DISCTYPE == SDFEM) && (ParamDB->NSTYPE==1))
  {
      //ParamDB->NSTYPE = 2;
      //OutPut("NSTYPE changed from 1 to 2 because of SDFEM discretization "<< endl);
      OutPut("NSTYPE 1: only reduced SDFEM, only for 2D, fixed point, not skew !!!" << endl);
  }

  if ((ParamDB->DISCTYPE == SDFEM) && (ParamDB->NSTYPE==3))
  {
    ParamDB->NSTYPE = 4;
    OutPut("NSTYPE changed from 3 to 4 because of SDFEM discretization "<< endl);
  }

  if ((ParamDB->LAPLACETYPE == 1) && (ParamDB->NSTYPE ==1))
  {
    ParamDB->NSTYPE = 3 ;
    OutPut("NSTYPE changed from 1 to 3 because of LAPLACETYPE "<< endl);
  }

  if ((ParamDB->LAPLACETYPE == 1) && (ParamDB->NSTYPE ==2))
  {
    ParamDB->NSTYPE = 4 ;
    OutPut("NSTYPE changed from 2 to 4 because of LAPLACETYPE "<< endl);
  }

  // equal order
  if (ParamDB->NSTYPE == 14)
  {
      if (!(ParamDB->DISCTYPE == SDFEM))
      {
	  ParamDB->DISCTYPE = SDFEM;
	  OutPut("DISCTYPE changed to SDFEM !!!"<<endl);
      }
/*
      if (ParamDB->SC_SMOOTHER_SADDLE<3)
      {
	  ParamDB->SC_SMOOTHER_SADDLE+=2;
	  OutPut("SC_SMOOTHER_SADDLE changed to " << ParamDB->SC_SMOOTHER_SADDLE);
	  OutPut(" because of continuous pressure"<< endl);
      }
      if (ParamDB->SC_COARSE_SMOOTHER_SADDLE<3)
      {
	  ParamDB->SC_COARSE_SMOOTHER_SADDLE+=2;
	  OutPut("SC_COARSE_SMOOTHER_SADDLE changed to " << ParamDB->SC_COARSE_SMOOTHER_SADDLE);
	  OutPut(" because of continuous pressure"<< endl);
      }
*/
  }

  if (ParamDB->SOLVER_TYPE == 0) // AMG
  {
    ParamDB->SC_MG_TYPE_SADDLE=0;
  }
  
  if (ParamDB->PROBLEM_TYPE == 3)
  {
     if (ParamDB->PRESSURE_SEPARATION==1)
     {
        TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 1;
     }
     else
     {
        TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE = 0;
     }
     TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE = 1.0;    
  }

  // rotational form
  if (ParamDB->NSE_NONLINEAR_FORM==2||(ParamDB->NSE_NONLINEAR_FORM==4))
  {
      if (ParamDB->NSTYPE<=2)
      {
	  ParamDB->NSTYPE+=2;
	  OutPut("NSTYPE changed to " << ParamDB->NSTYPE);
	  OutPut(" because of NSE_NONLINEAR_FORM = " <<ParamDB->NSE_NONLINEAR_FORM << endl);
      }
      // change DISCTYPE for internal reasons
      if (ParamDB->DISCTYPE == 1)
      {
	  ParamDB->DISCTYPE = 4;
	  TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE = 0;
	  OutPut("DISCTYPE changed to 4 for internal reasons, turbulent viscosity is switched off."<< endl);
      }	  
  }

  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==5)
  {
      if (!((TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION)||
	    (TDatabase::ParamDB->DISCTYPE == VMS_PROJECTION_EXPL)))
      {
	  OutPut("TURBULENT_VISCOSITY_TYPE = 5 only defined for projection-based VMS methods"<<endl);
	  OutPut("Set different TURBULENT_VISCOSITY_TYPE !!!"<<endl);
	  exit(4711);
      }	  
  }

  // LOCAL_PROJECTION
  if (ParamDB->DISCTYPE == LOCAL_PROJECTION)
  {
    if (ParamDB->TENSOR_TYPE == 0)    
    {
    if (ParamDB->LP_FULL_GRADIENT)
    {
      if (ParamDB->LP_STREAMLINE)
      {
        ParamDB->LP_STREAMLINE = 0;
        OutPut("LP_STREAMLINE changed to " << ParamDB->LP_STREAMLINE);
        OutPut(" due to LP_FULL_GRADIENT = " << ParamDB->LP_FULL_GRADIENT << endl);
      }

      if (ParamDB->LP_DIVERGENCE)
      {
        ParamDB->LP_DIVERGENCE = 0;
        OutPut("LP_DIVERGENCE changed to " << ParamDB->LP_DIVERGENCE);
        OutPut(" due to LP_FULL_GRADIENT = " << ParamDB->LP_FULL_GRADIENT << endl);
      }
    } // end LP_FULL_GRADIENT

    if (ParamDB->LP_DIVERGENCE)
    {
      if (ParamDB->NSTYPE<=2)
      {
        ParamDB->NSTYPE+=2;
        OutPut("NSTYPE changed to " << ParamDB->NSTYPE);
        OutPut("LP_DIVERGENCE = " << ParamDB->LP_DIVERGENCE << endl);
      }
    }

    if(ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
      ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = ParamDB->LP_ORDER_DIFFERENCE;

    if(ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
      ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = ParamDB->LP_ORDER_DIFFERENCE;

    if(ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE == -123)
      ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = ParamDB->LP_ORDER_DIFFERENCE;

    if(ParamDB->LP_PRESSURE_ORDER_DIFFERENCE == -123)
      ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = ParamDB->LP_ORDER_DIFFERENCE;
    }

  } // end DISCTYPE == LOCAL_PROJECTION
  else
  {
    // switch off all local projection terms
    ParamDB->LP_FULL_GRADIENT = 0;
    ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    ParamDB->LP_STREAMLINE = 0;
    ParamDB->LP_STREAMLINE_COEFF = 0;
    ParamDB->LP_STREAMLINE_EXPONENT = 1;

    ParamDB->LP_DIVERGENCE = 0;
    ParamDB->LP_DIVERGENCE_COEFF = 0;
    ParamDB->LP_DIVERGENCE_EXPONENT = 1;

    ParamDB->LP_PRESSURE = 0;
    ParamDB->LP_PRESSURE_COEFF = 0;
    ParamDB->LP_PRESSURE_EXPONENT = 1;

    ParamDB->LP_ORDER_DIFFERENCE = 1;
    ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = 1;
    ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = 1;
    ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = 1;
    ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = 1;
  }
  
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)  
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
    switch (TDatabase::ParamDB->NSTYPE)
    {
      case 1: OutPut("Galerkin discretization for Oseen because of NSTYPE "
  << TDatabase::ParamDB->NSTYPE << endl);
  TDatabase::ParamDB->DISCTYPE =  1;
  break;
      case 14:  OutPut("SUPG/PSPG/grad-div discretization for Oseen because of NSTYPE "
  << TDatabase::ParamDB->NSTYPE << endl);
  TDatabase::ParamDB->DISCTYPE =  2;
  break;
      default:
  OutPut("No method for Oseen implemented for NSTYPE " << TDatabase::ParamDB->NSTYPE << endl);
  exit(4711);
    }
    if (ParamDB->SC_NONLIN_MAXIT_SADDLE > 1)
    {
      ParamDB->SC_NONLIN_MAXIT_SADDLE = 1;
      OutPut("Set SC_NONLIN_MAXIT_SADDLE " << ParamDB->SC_NONLIN_MAXIT_SADDLE <<
      " for Oseen, further assembling not implemented"<<endl);
      TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE = 1.0;
      OutPut("Set TDatabase::ParamDB->SC_NONLIN_DAMP_FACTOR_SADDLE to 1.0" << endl);
    }
  }

} // end CheckParameterConsistencyNSE
