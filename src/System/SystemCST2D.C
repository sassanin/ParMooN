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
* @brief     source file for TSystemCST2D
* @author    Jagannath Venkatesan, 
* @date      18.08.15
* @History   Moved aux to main program, Removed FE_Vel in constructor
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <SystemCST2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <stdlib.h>
#include <string.h>
#include <MainUtilities.h>


TSystemCST2D::TSystemCST2D(TFESpace2D *stress_fespace, TFEVectFunct2D *Stress, int tensortype, int disctype, int solver, 
			   TFESpace2D *Velocity_FeSpace, TFESpace2D* Pressure_FeSpace, TFESpace2D* Deformation_FeSpace, TFEVectFunct2D *Velocity)
{
  //store the FEspace
  FeSpace[0] = Velocity_FeSpace;
  FeSpace[1] = Pressure_FeSpace;
  
    //set the tensor type
  Tensortype = tensortype;
  
  if (Tensortype == 1)
  {
  FeSpace[2] = stress_fespace;
  if(TDatabase::ParamDB->TENSOR_TYPE ==2)
  FeSpace[3] = Deformation_FeSpace;
  }
  else if(Tensortype == 2 )
  {
  FeSpace[2] = Deformation_FeSpace;
  FeSpace[3] = stress_fespace;
  }
  FeFct[0] = Stress->GetComponent(0);
  FeFct[1] = Stress->GetComponent(1); 
  FeFct[2] = Stress->GetComponent(2);
  
  // no. of DOFs
  N_S = stress_fespace->GetN_DegreesOfFreedom();
  N_Active =  stress_fespace->GetActiveBound();
  N_DirichletDof = N_S - N_Active;  
  

  
  //set the discretization type
  Disctype = disctype;
  
  //set the solver type
  SOLVER = solver;
  
  // build matrices
  // first build matrix structure
  
  if (Disctype == LOCAL_PROJECTION)
  {
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
  }

  if (Tensortype == 1)
  {
  sqstructure = new TSquareStructure2D(FeSpace[2]);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order
  /** S is the stiffness/system mat for stationary problem   */
  SqmatrixS11 = new TSquareMatrix2D(sqstructure);  
  SqmatrixS12 = new TSquareMatrix2D(sqstructure); 
  SqmatrixS21 = new TSquareMatrix2D(sqstructure); 
  SqmatrixS22 = new TSquareMatrix2D(sqstructure); 
  SqmatrixS23 = new TSquareMatrix2D(sqstructure);
  SqmatrixS32 = new TSquareMatrix2D(sqstructure); 
  SqmatrixS33 = new TSquareMatrix2D(sqstructure);
  
  N_Matrices = 7;
  }
  else if (Tensortype == 2)
  {
      sqstructure = new TSquareStructure2D(FeSpace[3]);
  sqstructure->Sort();
      /** S is the stiffness/system mat for stationary problem   */
  SqmatrixS11 = new TSquareMatrix2D(sqstructure);  
  SqmatrixS22 = new TSquareMatrix2D(sqstructure); 
  SqmatrixS33 = new TSquareMatrix2D(sqstructure);
  
  N_Matrices = 3;
  }
  else
  {
    cout<<"Invalid tensor type !!!!!"<<endl;
    exit(0);
  }
    
     // matrices for methods
   sqmatrices = (TSquareMatrix **)SQMATRICES;

   NSEaux_error = NULL;
   NSEaux = NULL;
}

TSystemCST2D::~TSystemCST2D()
{
   if (Tensortype == 1)
   {
   delete sqstructure;
   delete SqmatrixS11;
   delete SqmatrixS12;
   delete SqmatrixS21;
   delete SqmatrixS22;
   delete SqmatrixS23;
   delete SqmatrixS32;
   delete SqmatrixS33;
  }
  else if(Tensortype == 2)
  {
    delete sqstructure;
   delete SqmatrixS11;
   delete SqmatrixS22;
   delete SqmatrixS33;
  }
   
}
  
  
void TSystemCST2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *S1BoundValue, BoundValueFunct2D *S2BoundValue, BoundValueFunct2D *S3BoundValue, TAuxParam2D *aux, TAuxParam2D *auxerror)
{
  TDiscreteForm2D *DiscreteFormGalerkin, *DiscreteFormSDFEM;
  
    // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  
  BoundaryConditions[2] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = S1BoundValue;
  BoundaryValues[1] = S2BoundValue;
  BoundaryValues[2] = S3BoundValue;
 
  // save the CST bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  

  NSEaux = aux;
 
  //  aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {    
     NSEaux_error =  auxerror;             
    }
  
    if (Tensortype == 1)
   {
  // set the Discreteforms
  InitializeDiscreteForms_CST(DiscreteFormGalerkin, DiscreteFormSDFEM, LinCoeffs[0]);
   }
   else if (Tensortype ==2)
   {
     // set the Discreteforms
  InitializeDiscreteForms_DFT(DiscreteFormGalerkin, LinCoeffs[0]);
   }
  
    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormGalerkin;
           break;
      case SDFEM:
           DiscreteFormARhs = DiscreteFormSDFEM;
           break;
      
      default:
            OutPut("Unknown DISCTYPE " << endl);
            exit(4711);;
     }  
     
     
} 


void TSystemCST2D::Assemble(double *sol, double *rhs)
{

  double *RHSs[3];
  TFESpace2D *ferhs[3], *fesp[3];
  int N_FESpace;
     memset(rhs, 0, (3*N_S)*SizeOfDouble);
    RHSs[0] = rhs;
    RHSs[1] = rhs + N_S;
    RHSs[2] = rhs + 2*N_S;
 
//     fesp[0] = FeSpace;


   
        if (Tensortype == 1)
   {
     ferhs[0] = FeSpace[2];
    ferhs[1] = FeSpace[2];
    ferhs[2] = FeSpace[2];
    // initialize matrices
    SQMATRICES[0] = SqmatrixS11;
    SQMATRICES[1] = SqmatrixS12;
    SQMATRICES[2] = SqmatrixS21;
    SQMATRICES[3] = SqmatrixS22;
    SQMATRICES[4] = SqmatrixS23;
    SQMATRICES[5] = SqmatrixS32;
    SQMATRICES[6] = SqmatrixS33;

    SQMATRICES[0]->Reset(); 
    SQMATRICES[1]->Reset(); 
    SQMATRICES[2]->Reset(); 
    SQMATRICES[3]->Reset(); 
    SQMATRICES[4]->Reset(); 
    SQMATRICES[5]->Reset(); 
    SQMATRICES[6]->Reset(); 
   }
      else if (Tensortype == 2)
   {
     ferhs[0] = FeSpace[3];
    ferhs[1] = FeSpace[3];
    ferhs[2] = FeSpace[3];
         // initialize matrices
    SQMATRICES[0] = SqmatrixS11;
    SQMATRICES[1] = SqmatrixS22;
    SQMATRICES[2] = SqmatrixS33;

    SQMATRICES[0]->Reset(); 
    SQMATRICES[1]->Reset(); 
    SQMATRICES[2]->Reset(); 
   }
   if(TDatabase::ParamDB->TENSOR_TYPE ==1)
     N_FESpace = 3;
   else if (TDatabase::ParamDB->TENSOR_TYPE ==2)
     N_FESpace = 4;

    // assemble
    Assemble2D(N_FESpace, FeSpace,
               N_Matrices, SQMATRICES,
               0, NULL,
               3, RHSs, ferhs,
               DiscreteFormARhs,
               BoundaryConditions,
               BoundaryValues,
               NSEaux);
 
     // apply local projection stabilization method
     if(Disctype==LOCAL_PROJECTION)
      {
	// cout<<"LPS!!!!!!\n";
	double delta0 = TDatabase::ParamDB->DELTA0;
         AddStreamlineTerm(SQMATRICES[0], SQMATRICES[3], 
			      SQMATRICES[6],
			      u1, u2, delta0, 1.0, 1);
	 
	 AddStretchingTerm(SQMATRICES, u1, u2, delta0, 1.0, 1);
      }
    
    
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);     
    memcpy(sol+N_S+N_Active, rhs+N_S+N_Active, N_DirichletDof*SizeOfDouble);  
    memcpy(sol+(2*N_S)+N_Active, rhs+(2*N_S)+N_Active, N_DirichletDof*SizeOfDouble);  
     
} // void TSystemCD2D::Assemble(T


void TSystemCST2D::Solve(double *sol, double *rhs)
{
  
    switch(SOLVER)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        
	  if (Tensortype == 1)
	DirectSolver(SqmatrixS11, SqmatrixS12, SqmatrixS21, SqmatrixS22, SqmatrixS23, SqmatrixS32, SqmatrixS33, rhs, sol);
         else if (Tensortype == 2)
        DirectSolver(SqmatrixS11, SqmatrixS22, SqmatrixS33, rhs, sol);
       
	break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}

void TSystemCST2D::MeasureErrors(DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *s_error)
{
  double errors[4];
  TFESpace2D *fesp[1];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    if (Tensortype == 1)
  fesp[0] = FeSpace[2];
else if (Tensortype == 2)
  fesp[0] = FeSpace[3];
     // errors in first stress component
     FeFct[0]->GetErrors(ExactS1, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, fesp, errors);
      s_error[0] = errors[0];
      s_error[1] = errors[1];
      
     // errors in second stress component
     FeFct[1]->GetErrors(ExactS2, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, fesp, errors);
     s_error[2] = errors[0];
     s_error[3] = errors[1];      
      
      // errors in third stress component
     FeFct[2]->GetErrors(ExactS3, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, fesp, errors);     
     s_error[4] = errors[0];
     s_error[5] = errors[1];        
}


#endif // #ifdef __2D__
