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
* @brief     source file for TSystemNSE2D
* @author    Sashikumaar Ganesan, 
* @date      23.08.14
* @History   06.02.16 : Updated vms and stress fespaces
 ************************************************************************  */

#ifdef __2D__

#include <Database.h>
#include <SystemNSE2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <FEVectFunct2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
// #include <NSE2D_ParamRout.h>
#include <MainUtilities.h>
#include <Upwind.h>

#include <stdlib.h>
#include <string.h>

#ifdef __PRIVATE__  
 #include <VMS.h>
#endif

TSystemNSE2D::TSystemNSE2D(TFESpace2D *velocity_fespace, TFESpace2D *presssure_fespace, TFEVectFunct2D *Velocity, 
                           TFEFunction2D *p, double *sol, double *rhs, int disctype, int nsetype, int solver
#ifdef __PRIVATE__  
                           ,TFESpace2D *Projection_space, TFESpace2D *Stress_FeSpace,   TFESpace2D *Deformation_FeSpace

 #endif    
)
{
  //store the FEspaces and fefunct
  FeSpaces[0] = velocity_fespace;
  FeSpaces[1] = presssure_fespace;

#ifdef __PRIVATE__   
  if (TDatabase::ParamDB->TENSOR_TYPE ==1 || TDatabase::ParamDB->TENSOR_TYPE ==2)
   { 
    FeSpaces[2] = Stress_FeSpace;
    
    if (TDatabase::ParamDB->TENSOR_TYPE ==2) 
     {
      FeSpaces[3] = Deformation_FeSpace;
     }
  }
#endif

  N_U = velocity_fespace->GetN_DegreesOfFreedom();
  N_P = presssure_fespace->GetN_DegreesOfFreedom();

  N_Active =  velocity_fespace->GetActiveBound();

  N_DirichletDof = N_U - N_Active;  
  
  VelocityFct = Velocity;
  
  FeFct[0] = Velocity->GetComponent(0);
  FeFct[1] = Velocity->GetComponent(1); 
  FeFct[2] = p;
   
  Sol = sol;
  Rhs = rhs;
  
  RHSs[0] = Rhs;
  RHSs[1] = Rhs + N_U;
  RHSs[2] = Rhs + 2*N_U;
 

  //set the discretization type
  Disctype = disctype;
  
  // NSE type
  NSEType = nsetype;
  
  //set the solver type
  Solver = solver;
 
  // build matrices
  // first build matrix structure
  sqstructureA = new TSquareStructure2D(FeSpaces[0]);
  sqstructureA->Sort();  // sort column numbers: numbers are in increasing order
  
  structureB = new TStructure2D(FeSpaces[1], FeSpaces[0]);
  structureBT = new TStructure2D(FeSpaces[0], FeSpaces[1]);
    
    switch(NSEType)
     {
      case 1:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE1;
      break;

      case 2:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE2;
      break;

      case 3:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE3;
      break;

      case 4:
        MatrixB1 = new TMatrix2D(structureB);
        MatrixB2 = new TMatrix2D(structureB);
        MatrixB1T = new TMatrix2D(structureBT);
        MatrixB2T = new TMatrix2D(structureBT);

        SqmatrixA11 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA12 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA21 = new TSquareMatrix2D(sqstructureA);
        SqmatrixA22 = new TSquareMatrix2D(sqstructureA);
        Defect = Defect_NSE4;
      break;
      
      default:
            OutPut("Unknown NSETYPE, it must be 1 to 4" << endl);
            exit(4711);;      
      
     }  
 
 #ifdef __PRIVATE__ 
   if (Disctype == VMS_PROJECTION)
   { 
     if(NSEType==1 || NSEType==2)
     {
      Error("NSETYPE should be 3 or 4 for VMS_PROJECTION !!!" << endl);
      exit(-1);          
     }
    
    if (TDatabase::ParamDB->TENSOR_TYPE ==1)
    {
    FeSpaces[3] = velocity_fespace; //  to be included the convolution space if needed
    FeSpaces[4] = Projection_space; 
    }
    else 
    {
    FeSpaces[2] = velocity_fespace; //  to be included the convolution space if needed
    FeSpaces[3] = Projection_space;   
    }
 
    sqstructureL = new TSquareStructure2D(Projection_space);
    sqstructureL->Sort();
    structure_tilde_G = new TStructure2D(velocity_fespace, Projection_space);
    structure_G = new TStructure2D(Projection_space, velocity_fespace);    

    MatricesL = new TSquareMatrix2D(sqstructureL);
    Matrices_tilde_G11 = new TMatrix2D(structure_tilde_G);
    Matrices_tilde_G22 = new TMatrix2D(structure_tilde_G);    
    Matrices_G11 = new TMatrix2D(structure_G);
    Matrices_G22 = new TMatrix2D(structure_G);
   }
#endif

   // matrices for methods
   sqmatrices = (TSquareMatrix **)SQMATRICES;
   matrices = (TMatrix **)MATRICES;
   
   NSEaux_error = NULL;
   NSEaux = NULL;
}

TSystemNSE2D::~TSystemNSE2D()
{
  delete sqstructureA; delete structureB; delete structureBT;
  delete MatrixB1; delete MatrixB2; delete MatrixB1T; delete MatrixB2T;
  delete SqmatrixA11; delete SqmatrixA12; delete SqmatrixA21; delete SqmatrixA22;
}


void TSystemNSE2D::Init(CoeffFct2D *lincoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *U1BoundValue, BoundValueFunct2D *U2BoundValue, TAuxParam2D *aux, TAuxParam2D *auxerror)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces;

 TFESpace2D *fesprhs[3];
  
  TDiscreteForm2D *DiscreteFormGalerkin, *DiscreteFormSDFEM, *DiscreteFormUpwind, *DiscreteFormSmagorinsky;
  TDiscreteForm2D *DiscreteFormVMSProjection, *DiscreteFormNLGalerkin, *DiscreteFormNLSDFEM, *DiscreteFormNLUpwind;
  TDiscreteForm2D *DiscreteFormNLSmagorinsky, *DiscreteFormNLVMSProjection, *DiscreteFormPressSep, *DiscreteFormAuxProbPressSep;
  TDiscreteForm2D *DiscreteFormNSRFBRhs, *DiscreteFormNSCSTRhs;
    
  // save the boundary condition
  BoundaryConditions[0] = BoundCond;
  BoundaryConditions[1] = BoundCond;  

  // save the boundary values  
  BoundaryValues[0] = U1BoundValue;
  BoundaryValues[1] = U2BoundValue;
 
  // save the nse bilinear coefficient   
  LinCoeffs[0] = lincoeffs;
  
  
   NSEaux = aux;
 
  //  aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {    
     NSEaux_error =  auxerror;             
    }
    
  // set the Discreteforms
  InitializeDiscreteForms(
    DiscreteFormGalerkin, DiscreteFormSDFEM,
    DiscreteFormUpwind, DiscreteFormSmagorinsky,
    DiscreteFormVMSProjection,
    DiscreteFormNLGalerkin, DiscreteFormNLSDFEM,
    DiscreteFormNLUpwind, DiscreteFormNLSmagorinsky,
    DiscreteFormNLVMSProjection,
    DiscreteFormPressSep,
    DiscreteFormAuxProbPressSep,
    DiscreteFormNSRFBRhs,
    DiscreteFormNSCSTRhs,
    LinCoeffs[0], NSEType);

    // find discrete form
    switch(Disctype)
       {
	  case LOCAL_PROJECTION:
          case GALERKIN:
          DiscreteFormARhs = DiscreteFormGalerkin;
            DiscreteFormNL = DiscreteFormNLGalerkin;
	    DiscreteFormRhs = DiscreteFormNSCSTRhs;
          break;

          case SDFEM:
            DiscreteFormARhs = DiscreteFormSDFEM;
            DiscreteFormNL = DiscreteFormNLSDFEM; 
          break;

          case UPWIND:
            DiscreteFormARhs = DiscreteFormUpwind;
            DiscreteFormNL = DiscreteFormNLUpwind;    
            break;

          case SMAGORINSKY:
            DiscreteFormARhs = DiscreteFormSmagorinsky;
            DiscreteFormNL = DiscreteFormNLSmagorinsky;              
            break;

          case VMS_PROJECTION:
	      DiscreteFormARhs = DiscreteFormVMSProjection;
	      DiscreteFormNL = DiscreteFormNLVMSProjection;
	      DiscreteFormRhs = DiscreteFormNSCSTRhs;
	      if (TDatabase::ParamDB->NSTYPE != 1 && TDatabase::ParamDB->NSTYPE != 4)
	       {
                OutPut("VMS only for NSTYPE 1 and 4 implemented !!!"<<endl);
		exit(4711);
	       }
            break;

          default:
            Error("Unknown DISCTYPE" << endl);
            exit(-1);
        } 
     
     // set the discrete form for the Stokes equation
      if (TDatabase::ParamDB->PROBLEM_TYPE==STOKES)
       {
        DiscreteFormARhs = DiscreteFormUpwind;     
        DiscreteFormNL = NULL;
       }
       
       
      // initilize the assemble          
      switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 2;
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 1;
          N_RectMatrices = 4;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 2;
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;

          SQMATRICES[0]->Reset();
          SQMATRICES[1]->Reset();
          SQMATRICES[2]->Reset();
          SQMATRICES[3]->Reset();
          MATRICES[0]->Reset();
          MATRICES[1]->Reset();
          MATRICES[2]->Reset();
          MATRICES[3]->Reset();

          N_SquareMatrices = 4;
          N_RectMatrices = 4;

#ifdef __PRIVATE__  
        if(Disctype == VMS_PROJECTION)
        {
          N_SquareMatrices = 5;
          SQMATRICES[4] =  MatricesL;
          SQMATRICES[4]->Reset();
 
          N_RectMatrices = 8;
          MATRICES[4] = Matrices_tilde_G11;
          MATRICES[5] = Matrices_tilde_G22;
          MATRICES[6] = Matrices_G11;
          MATRICES[7] = Matrices_G22;
          MATRICES[4]->Reset();
          MATRICES[5]->Reset();
          MATRICES[6]->Reset();
          MATRICES[7]->Reset();

          N_FESpaces = 4;
        }
#endif 
  
          break;
      } //  switch(NSEType)      

    N_Rhs = 2;
    N_FESpaces = 2;   
#ifdef __PRIVATE__       
  if (TDatabase::ParamDB->TENSOR_TYPE ==2)
   { N_FESpaces+=2; }
  else if (TDatabase::ParamDB->TENSOR_TYPE ==1)
   { N_FESpaces+=1; }   
#endif     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];      
   
     //assemble object
     AMatRhsAssemble = new TAssembleMat2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES, N_RectMatrices, MATRICES,
                              N_Rhs, RHSs, fesprhs, DiscreteFormARhs, BoundaryConditions, BoundaryValues, NSEaux);   
    
     AMatRhsAssemble->Init();     
     
  // ===============================================================================================================
  // initilize the nonliner assemble    
     switch(TDatabase::ParamDB->NSTYPE)
       {
        case 1:
        case 2:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[0]->Reset();

          N_SquareMatrices = 1;
        break;

        case 3:
        case 4:
          if (TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE==0)
           {
            SQMATRICES[0] = SqmatrixA11;
            SQMATRICES[1] = SqmatrixA22;
            SQMATRICES[0]->Reset();
            SQMATRICES[1]->Reset();

            N_SquareMatrices = 2;
           }
          else
           {
            // Newton method
            cout<< "Newton method not tested " <<endl;
            exit(0);
           }

         break;
        } // switch(TDatabase::ParamDB->NSTYPE)
            
      N_RectMatrices = 0;          
      N_Rhs = 0;
      N_FESpaces = 1;     
     
     AMatAssembleNonLinear = new TAssembleMat2D(N_FESpaces, FeSpaces, N_SquareMatrices, SQMATRICES, N_RectMatrices, NULL,
                              N_Rhs, NULL, NULL, DiscreteFormNL, BoundaryConditions, BoundaryValues, NSEaux);   
    
     AMatAssembleNonLinear->Init();    
     
/*    cout << "TAssembleMat2D  () done ! " << endl;
  exit(0);   */   
} // TSystemNSE2D::Init

 
void TSystemNSE2D::UpdateUpwind()
{
        switch(NSEType)
         {
          case 1:
          case 2:
            // do upwinding with one matrix
            UpwindForNavierStokes(LinCoeffs[0], SqmatrixA11, FeFct[0], FeFct[1]);
            cout << "UPWINDING DONE : level " << endl;
            break;

          case 3:
          case 4:
            // do upwinding with two matrices
            cout << "UPWINDING DONE : level " << endl;
            UpwindForNavierStokes(LinCoeffs[0], SqmatrixA11, FeFct[0], FeFct[1]);
            UpwindForNavierStokes(LinCoeffs[0], SqmatrixA22, FeFct[0], FeFct[1]);
	    break;
         }                        // endswitch 
} // TSystemNSE2D::UpdateUpwind()

  
void TSystemNSE2D::UpdateLPS()
{
   switch(NSEType)
      {
	   case 1:
	   case 3:
	     OutPut("LPS only for NSTYPE 2 and 4 implemented !!!"<<endl);
		exit(4711);
	     break;
	   case 2:
	     UltraLocalProjection(SqmatrixA11, 0);
	     break;
	   case 4:
	     
	     AddDeformationTensorTerm(SqmatrixA11,SqmatrixA12,
				      SqmatrixA21,SqmatrixA22,
				      TDatabase::ParamDB->DELTA1, 2.0, 1);
	     break;
	     default:
            Error("Unknown NSTYPE" << endl);
            exit(-1);
        }  
  
}


void TSystemNSE2D::Assemble(double *sol, double *rhs)
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces; 
 double *RHSs[3];  
 TFESpace2D *fesprhs[3];

  /** initialize matrices */
  AMatRhsAssemble->Reset();
  
  /** assemble */
  AMatRhsAssemble->Assemble2D();  
  
  /** apply local projection stabilization method */
  if(Disctype==LOCAL_PROJECTION)
   {
    cout<<"LPS!!!!!!\n";
    this->UpdateLPS();
   }
 
  /** upwind */ 
  if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
   {
    this->UpdateUpwind(); 
   }
            
  // slip with friction boundary condition
  if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
   {
    if(NSEType <4)
     {
      OutPut("For slip with friction bc NSTYPE 4 is ");
      OutPut("necessary !!!!! " << endl);
      exit(4711);
     }
 
     N_FESpaces = 1;
     N_SquareMatrices = 4;
     N_RectMatrices = 2;
     N_Rhs = 2;
     
     RHSs[0] = rhs;
     RHSs[1] = rhs + N_U;
     RHSs[2] = rhs + 2*N_U;

        fesprhs[0] = FeSpaces[0];
        fesprhs[1] = FeSpaces[0];
        fesprhs[2] = FeSpaces[1];
      
        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
        SQMATRICES[2] = SqmatrixA12;
        SQMATRICES[3] = SqmatrixA21;

        MATRICES[0] = MatrixB1T;
        MATRICES[1] = MatrixB2T;

        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, MATRICES,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=      
    
     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
      
     
#ifdef __PRIVATE__   
      // update matrices
      if (Disctype == VMS_PROJECTION)
        {
         SQMATRICES[0] = SqmatrixA11;
         SQMATRICES[1] = SqmatrixA12;
         SQMATRICES[2] = SqmatrixA21;
         SQMATRICES[3] = SqmatrixA22;
         SQMATRICES[6] =  MatricesL;
         MATRICES[2] = Matrices_tilde_G11;
         MATRICES[3] = Matrices_tilde_G22;
         MATRICES[4] = Matrices_G11;
         MATRICES[5] = Matrices_G22;
         if (TDatabase::ParamDB->TENSOR_TYPE ==2)
         VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[5]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
	 else if(TDatabase::ParamDB->TENSOR_TYPE ==1)
	 VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[4]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
	 else
	 {
	 VMSProjectionUpdateMatrices(N_U, FeSpaces[0]->GetActiveBound(), FeSpaces[3]->GetN_DegreesOfFreedom(),
                                     SQMATRICES, MATRICES);
        }
	}
#endif  
     
//     cout << "Test Assemble " << endl; 
} // TSystemNSE2D::Assemble(T

void TSystemNSE2D::AssembleNonLinear(double *sol, double *rhs)   // yet to be updated with stress and deformation spaces in case of DEVSS
{
 int N_SquareMatrices, N_RectMatrices, N_Rhs, N_FESpaces, last_sq;

 TAuxParam2D *NLaux; 
  double *RHSs[3];

  TFESpace2D *fesprhs[3];

      // reset the nonliner matrices
      AMatAssembleNonLinear->Reset();
            
      // assemble the nonlinear matrix */      
      AMatAssembleNonLinear->Assemble2D(); 
      
     /** apply local projection stabilization method */
     if(Disctype==LOCAL_PROJECTION)
      {
       cout<<"LPS!!!!!!\n";
       this->UpdateLPS();
      }
 
     /** upwind */ 
     if( (Disctype==UPWIND) && !(TDatabase::ParamDB->PROBLEM_TYPE == 3) )
      {
       this->UpdateUpwind(); 
      }
       
      // slip with boundary condition
      if (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >= 1)
      { 
        N_FESpaces = 1;
        N_SquareMatrices = 2;
        N_RectMatrices = 0;
        N_Rhs = 2;

        SQMATRICES[0] = SqmatrixA11;
        SQMATRICES[1] = SqmatrixA22;
	
        RHSs[0] = rhs;
        RHSs[1] = rhs + N_U;
	
        fesprhs[0] = FeSpaces[0];
        fesprhs[1] = FeSpaces[0];
      
        NLaux = new TAuxParam2D(1, 0, 0, 0, FeSpaces, NULL, NULL, NULL, NULL, 0, NULL);

        Assemble2DSlipBC(N_FESpaces, FeSpaces,
                         N_SquareMatrices, SQMATRICES,
                         N_RectMatrices, NULL,
                         N_Rhs, RHSs, fesprhs,
                         NULL,
                         BoundaryConditions,
                         BoundaryValues,
                         NSEaux,
                         FeFct[0], FeFct[1]);

	
	delete NLaux;
      }// (TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION >=         
     
      
     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble);       
      
      
      
} //TSystemNSE2D::AssembleNonLinear(

// Used for the case of conformation stress in rhs
  void TSystemNSE2D::AssembleRhsOnly(double *sol, double *rhs)
{
  int N_Rhs, N_FESpaces;
  
  double *RHSs[3];

  TFESpace2D *fesprhs[3];

      N_Rhs = 2;
      N_FESpaces = 2;   
     
      RHSs[0] = rhs;
      RHSs[1] = rhs + N_U;
      RHSs[2] = rhs + 2*N_U;
      memset(rhs, 0, (2*N_U+N_P)*SizeOfDouble);
     
      fesprhs[0] = FeSpaces[0];
      fesprhs[1] = FeSpaces[0];
      fesprhs[2] = FeSpaces[1];
      
       if (TDatabase::ParamDB->TENSOR_TYPE == 1 || TDatabase::ParamDB->TENSOR_TYPE ==2)
      { 
	++N_FESpaces;
	
	if(TDatabase::ParamDB->TENSOR_TYPE ==2)
	  ++N_FESpaces;
    
  }
#ifdef __PRIVATE__   
    if (Disctype == VMS_PROJECTION)
     {
      N_FESpaces += 2;
     }
#endif  

      // assemble
      Assemble2D(N_FESpaces, FeSpaces,
        0, NULL,
        0, NULL,
        N_Rhs, RHSs, fesprhs,
        DiscreteFormRhs,
        BoundaryConditions,
        BoundaryValues,
        NSEaux);
     
     // set rhs for Dirichlet nodes
     memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);
     memcpy(sol+N_U+N_Active, rhs+N_U+N_Active, N_DirichletDof*SizeOfDouble); 
      
}
  
void TSystemNSE2D::GetResidual(double *sol, double *rhs, double *res)
{
  
     switch(NSEType)
      {
        case 1:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
        break;

        case 2:
          SQMATRICES[0] = SqmatrixA11;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
        break;

        case 3:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
        break;

        case 4:
          SQMATRICES[0] = SqmatrixA11;
          SQMATRICES[1] = SqmatrixA12;
          SQMATRICES[2] = SqmatrixA21;
          SQMATRICES[3] = SqmatrixA22;
          MATRICES[0] = MatrixB1;
          MATRICES[1] = MatrixB2;
          MATRICES[2] = MatrixB1T;
          MATRICES[3] = MatrixB2T;
          break;
      } //  switch(NSEType)  
  
  
   Defect(sqmatrices, matrices, sol, rhs, res); 

   
   
//    double *defect = new double[2*N_U+N_P];
//    
//    memset(defect, 0, 2*N_U+N_P*SizeOfDouble);
//    memset(res, 0, 2*N_U+N_P*SizeOfDouble);
// 
//       // Block A
//    MatVect(SqmatrixA11, sol, defect);
//    Daxpy(N_U, 1, defect, res);  
//    MatVectActive(SqmatrixA12, sol+N_U, defect);
//    Daxpy(N_Active, 1, defect, res);
//    MatVectActive(SqmatrixA21, sol, defect+N_U);
//    Daxpy(N_Active, 1, defect+N_U, res+N_U); 
//    MatVect(SqmatrixA22, sol+N_U, defect+N_U);
//    Daxpy(N_U, 1, defect+N_U, res+N_U);
//    
//       // Block BT
//       MatVect1(MatrixB1T, sol+2*N_U, defect);
//       Daxpy(N_Active, 1, defect, res);    
//       MatVect1(MatrixB2T, sol+2*N_U, defect+N_U);
//       Daxpy(N_Active, 1, defect+N_U, res+N_U); 
//      
//       // Block B
//       MatVect1(MatrixB1, sol, defect+2*N_U);
//       Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 
//       MatVect1(MatrixB2, sol+N_U, defect+2*N_U);
//       Daxpy(N_P, 1, defect+2*N_U, res+2*N_U); 
// 
//    
//     Daxpy(2*N_U+N_P, -1., rhs, res);


} // TSystemNSE2D::GetResidual

void TSystemNSE2D::Solve(double *sol, double *rhs)
{
  
    switch(Solver)
     {
      case AMG_SOLVE:
        cout << "AMG_SOLVE not yet implemented " <<endl;
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        switch(NSEType)
         {
          case 1:
            DirectSolver(SqmatrixA11, MatrixB1,  MatrixB2, rhs, sol);
          break;

          case 2:
             DirectSolver(SqmatrixA11, MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol);
          break;

          case 3:
           cout << "Solver not included for NSTYPE 3 in this version" <<endl;
            cout << "try NSTYPE 4 " <<endl;   
	    exit(0);
          break;

          case 4:    
             DirectSolver(SqmatrixA11, SqmatrixA12, SqmatrixA21, SqmatrixA22, 
                          MatrixB1T, MatrixB2T, MatrixB1,  MatrixB2, rhs, sol); 
	     
          break;
      } //  switch(NSEType) 

      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  
}

void TSystemNSE2D::MeasureErrors(DoubleFunct2D *ExactU1, DoubleFunct2D *ExactU2, DoubleFunct2D *ExactP,
                                    double *u_error, double *p_error)
{
  double errors[4];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
     // errors in first velocity component
     FeFct[0]->GetErrors(ExactU1, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
      u_error[0] = errors[0];
      u_error[1] = errors[1];
      
     // errors in second velocity component
     FeFct[1]->GetErrors(ExactU2, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces, errors);
     u_error[2] = errors[0];
     u_error[3] = errors[1];      
      
      // errors in pressure
     FeFct[2]->GetErrors(ExactP, 3, AllDerivatives, 2,
                         L2H1Errors,
                         NULL, NSEaux_error, 1, FeSpaces+1, errors);     
     p_error[0] = errors[0];
     p_error[1] = errors[1];        
}
    
    
#endif // #ifdef __3D__
