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
   
/** ************************************************************************ 
* @brief     source file for TSystemDarcy2D
* @author    Ulrich Wilbrandt,
* @date      15.03.15
 ************************************************************************  */
#include <Database.h>
#include <SystemDarcy2D.h>
#include <Darcy2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>

TSystemDarcy2D::TSystemDarcy2D(TFESpace2D **fespaces)
{
  //store the FEspace
  fe_spaces[0] = fespaces[0]; // velocity
  fe_spaces[1] = fespaces[1]; // pressure
  
  // build matrices
  // first build matrix structures
  // velocity-velocty coupling
  TSquareStructure2D *sqstructureA = new TSquareStructure2D(fe_spaces[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  // pressure-pressure coupling
  TSquareStructure2D *sqstructureC = new TSquareStructure2D(fe_spaces[1]);
  sqstructureC->Sort();  // sort column numbers: numbers are in increasing order
  // velocity-pressure and pressure-velocity coupling
  TStructure2D *structureB = new TStructure2D(fe_spaces[1], fe_spaces[0]);
  TStructure2D *structureBT = new TStructure2D(fe_spaces[0], fe_spaces[1]);
  
  /** create the velocity-velocity coupling matrix */
  sq_matrices[0] = new TSquareMatrix2D(sqstructureA);
  /** create the pressure-pressure coupling matrix */
  sq_matrices[1] = new TSquareMatrix2D(sqstructureC);
  /** create velocity-pressure and pressure-velocity coupling matrices */
  rect_matrices[0] = new TMatrix2D(structureBT);
  rect_matrices[1] = new TMatrix2D(structureB);
  
  N_Matrices = 4;
}

TSystemDarcy2D::~TSystemDarcy2D()
{
  delete sq_matrices[0]->GetStructure();
  delete sq_matrices[1]->GetStructure();
  delete rect_matrices[0]->GetStructure();
  delete rect_matrices[1]->GetStructure();
  delete sq_matrices[0];
  delete sq_matrices[1];
  delete rect_matrices[0];
  delete rect_matrices[1];
}
  
  
void TSystemDarcy2D::Init(BoundCondFunct2D **BoundCond, 
                             BoundValueFunct2D **BoundValue)
{
  BoundaryConditions[0] =  BoundCond[0];
  BoundaryConditions[1] =  BoundCond[1];
  BoundaryValues[0] = BoundValue[0];
  BoundaryValues[1] = BoundValue[1];
} // TSystemDarcy2D::Init


void TSystemDarcy2D::Assemble(LocalAssembling2D& la, double *sol,
                                 double *rhs)
{
  int N_U = fe_spaces[0]->GetN_DegreesOfFreedom();
  int N_U_Active = fe_spaces[0]->GetActiveBound();
  int N_P = fe_spaces[1]->GetN_DegreesOfFreedom();
  int N_DOF = N_U + N_P;
  
  memset(rhs, 0, N_DOF*SizeOfDouble);
  double *rhs_blocks[2] = { rhs, rhs+N_U };
 
  // initialize matrices
  sq_matrices[0]->Reset();
  sq_matrices[1]->Reset();
  rect_matrices[0]->Reset();
  rect_matrices[1]->Reset();
  
  MultiIndex2D derivatives[6] = { D00, D00, D10, D01, D10, D01 };
  int spacesNumbers[6] = { 0, 1, 0, 0, 1, 1};
  int rowSpace[4] = {0, 1, 0, 1};
  int columnSpace[4] = { 0, 1, 1, 0};
  int rhsSpace[2] = { 0, 1 };
  
  // assemble
  TDiscreteForm2D discreteForm((char*)"dummy",(char*)"dummy", 6, derivatives,
                               spacesNumbers, 4, 2, rowSpace, columnSpace,
                               rhsSpace, BilinearAssembleDarcyGalerkin,
                               la.GetCoeffFct(), NULL);
  Assemble2D_VectFE(2, fe_spaces, 2, sq_matrices, 2, rect_matrices, 2, 
                    rhs_blocks, fe_spaces, &discreteForm, BoundaryConditions,
                    BoundaryValues);
  
    
  // copy Dirichlet values from rhs to solution vector (this is not really 
  // necessary in case of a direct solver)
  memcpy(sol+N_U_Active, rhs+N_U_Active, (N_U-N_U_Active)*SizeOfDouble);
} // void TSystemDarcy2D::Assemble


void TSystemDarcy2D::Solve(double *sol, double *rhs)
{ 
  switch(TDatabase::ParamDB->SOLVER_TYPE)
  {
   case AMG_SOLVE:
     cout << "AMG_SOLVE not yet implemented " <<endl;
   break;

   case GMG:
     cout << "GMG solver not yet implemented " <<endl;
   break;

   case DIRECT:
     DirectSolver(sq_matrices[0], sq_matrices[1], rect_matrices[0], 
                  rect_matrices[1], rhs, sol);
   break;      
 
   default:
     OutPut("Unknown Solver" << endl);
     exit(4711);;
  }
}
