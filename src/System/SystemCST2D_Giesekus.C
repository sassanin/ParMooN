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
* @brief     source file for TSystemCST2D_Giesekus
* @author    Jagannath Venkatesan, 
* @date      03.09.15
* @History 
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <SystemCST2D_Giesekus.h>
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


TSystemCST2D_Giesekus::TSystemCST2D_Giesekus(TFESpace2D *stress_fespace, TFEVectFunct2D *Stress,  int disctype, int solver)
{
  //store the FEspace
  FeSpace = stress_fespace;
  
  // store the FeFunctions
  FeFct[0] = Stress->GetComponent(0);
  FeFct[1] = Stress->GetComponent(1); 
  FeFct[2] = Stress->GetComponent(2);
  
 
  N_S = stress_fespace->GetN_DegreesOfFreedom();
  N_Active =  stress_fespace->GetActiveBound();
  N_DirichletDof = N_S - N_Active;  
  
  //set the discretization type
  Disctype = disctype;
  
  //set the solver type
  SOLVER = solver;
  
  // build matrices
  // first build matrix structure
  sqstructure = new TSquareStructure2D(FeSpace);
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
  
     // matrices for methods
   sqmatrices = (TSquareMatrix **)SQMATRICES;

   
   CSTaux_error = NULL;
   CSTNSEaux = NULL;
}

TSystemCST2D_Giesekus::~TSystemCST2D_Giesekus()
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
  
  
void TSystemCST2D_Giesekus::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *S1BoundValue, BoundValueFunct2D *S2BoundValue, BoundValueFunct2D *S3BoundValue, TAuxParam2D *aux, TAuxParam2D *auxerror)
{
  TDiscreteForm2D *DiscreteFormGalerkin;
  
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
  
   CSTNSEaux = aux;
 
    // aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
        CSTaux_error =  auxerror;
    }
  
  
  // set the Discreteforms
  InitializeDiscreteForms_CST_Giesekus(DiscreteFormGalerkin, LinCoeffs[0]);
  
  
    switch(Disctype)
     {
      case GALERKIN:
           DiscreteFormARhs = DiscreteFormGalerkin;
           break;
      
      default:
            OutPut("Unknown DISCTYPE " << endl);
            exit(4711);;
     }  
     
     
} 


void TSystemCST2D_Giesekus::Assemble(double *sol, double *rhs)
{

  double *RHSs[3];
  TFESpace2D *ferhs[3], *fesp[3];
  
     memset(rhs, 0, (3*N_S)*SizeOfDouble);
    RHSs[0] = rhs;
    RHSs[1] = rhs + N_S;
    RHSs[2] = rhs + 2*N_S;

    fesp[0] = FeSpace;

    ferhs[0] = FeSpace;
    ferhs[1] = FeSpace;
    ferhs[2] = FeSpace;
    
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
 
    // assemble
    Assemble2D(1, fesp,
               N_Matrices, SQMATRICES,
               0, NULL,
               3, RHSs, ferhs,
               DiscreteFormARhs,
               BoundaryConditions,
               BoundaryValues,
               CSTNSEaux);
 
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);     
    memcpy(sol+N_S+N_Active, rhs+N_S+N_Active, N_DirichletDof*SizeOfDouble);  
    memcpy(sol+(2*N_S)+N_Active, rhs+(2*N_S)+N_Active, N_DirichletDof*SizeOfDouble);  
     
} // void TSystemCD2D::Assemble(T


void TSystemCST2D_Giesekus::GetResidual(double *sol, double *rhs, double *res)
{
  int *KColS, *RowPtrS, N_S, i, j, begin, end, posi, posj;
  double *EntriesS11, *EntriesS12, *EntriesS21, *EntriesS22, *EntriesS23, *EntriesS32, *EntriesS33;
  
  KColS = SqmatrixS11->GetKCol();
  RowPtrS = SqmatrixS11->GetRowPtr();
  N_S = SqmatrixS11->GetN_Rows();
  EntriesS11 = SqmatrixS11->GetEntries();
  EntriesS12 = SqmatrixS12->GetEntries();
  EntriesS21 = SqmatrixS21->GetEntries();
  EntriesS22 = SqmatrixS22->GetEntries();
  EntriesS23 = SqmatrixS23->GetEntries();
  EntriesS32 = SqmatrixS32->GetEntries();
  EntriesS33 = SqmatrixS33->GetEntries();

  for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
       posi = i ;
       posj = KColS[j];
       res[posi] += EntriesS11[j] * sol[posj];
       if(i<N_Active)
       {
       posj = KColS[j] + N_S;
       res[posi] += EntriesS12[j] * sol[posj];
       }
       
       posi = i + N_S;
       if(i<N_Active)
       {
       posj = KColS[j];
       res[posi] += EntriesS21[j] * sol[posj];
       }
       posj = KColS[j] + N_S;
       res[posi] += EntriesS22[j] * sol[posj];
       if(i<N_Active)
       {
       posj = KColS[j] + (2.0*N_S);
       res[posi] += EntriesS23[j] * sol[posj];
       }
       
       posi = i + (2*N_S);
        if(i<N_Active)
       {
       posj = KColS[j] + N_S;
       res[posi] += EntriesS32[j] * sol[posj];
       }
       posj = KColS[j] + (2*N_S);
       res[posi] += EntriesS33[j] * sol[posj];     
    }
  }
  
  for(i=0;i<(3*N_S);i++)
    res[i] = rhs[i] - res[i];

}


void TSystemCST2D_Giesekus::Solve(double *sol, double *rhs)
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
        DirectSolver(SqmatrixS11, SqmatrixS12, SqmatrixS21, SqmatrixS22, SqmatrixS23, SqmatrixS32, SqmatrixS33, rhs, sol);
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
  
}

void TSystemCST2D_Giesekus::MeasureErrors(DoubleFunct2D *ExactS1, DoubleFunct2D *ExactS2, DoubleFunct2D *ExactS3, double *s_error)
{
  double errors[4];
  TFESpace2D *fesp[1];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
    fesp[0] = FeSpace;


     // errors in first stress component
     FeFct[0]->GetErrors(ExactS1, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, CSTaux_error, 1, fesp, errors);

      s_error[0] = errors[0];
      s_error[1] = errors[1];


     // errors in second stress component
     FeFct[1]->GetErrors(ExactS2, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, CSTaux_error, 1, fesp, errors);
     s_error[2] = errors[0];
     s_error[3] = errors[1];      
      
      // errors in third stress component
     FeFct[2]->GetErrors(ExactS3, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, CSTaux_error, 1, fesp, errors);     
     s_error[4] = errors[0];
     s_error[5] = errors[1];        
}



#endif // #ifdef __2D__
