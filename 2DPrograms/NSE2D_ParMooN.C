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
//
// Purpose:     main program for solving a stationary NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.08.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
#include <NSE2D_ParamRout.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
// =======================================================================
// include current example
// =======================================================================
//#include "../Examples/NSE_2D/DrivenCavity.h" // Lid Driven Cavity 
#include "../Examples/NSE_2D/flow_past_cylinder.h" // Flow past cylinder using Gmsh
// #include "../Examples/NSE_2D/Poiseuille.h" 
// #include "../Examples/NSE_2D/Benchmark.h"
// #include "../Examples/NSE_2D/Hemker_Wide.h"

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, N_Cells, ORDER, N_U, N_P, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code;
  
  double *sol, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
  double limit, u_error[4], p_error[2];
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *mortarcoll = NULL;
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
  TFEVectFunct2D *Velocity;
  TFEFunction2D *u1, *u2, *Pressure, *fefct[3];
  TOutput2D *Output;
  TSystemNSE2D *SystemMatrix;
  TAuxParam2D *aux, *auxerror;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
   
  const char vtkdir[] = "VTK";
  char *PsBaseName, *VtkBaseName, *GEO;
  char UString[] = "u";
  char PString[] = "p";
  
  std::ostringstream os;
  os << " ";
      
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);
  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  ExampleFile();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh
  cout << " MESH TYPE :  " << TDatabase::ParamDB->MESH_TYPE << endl;
   if(TDatabase::ParamDB->MESH_TYPE==0)
   {
    GEO = TDatabase::ParamDB->GEOFILE;
    Domain->Init(NULL, GEO);
   }
     else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
         cout << " MESH - GMSH GEN " << endl;
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE); 
       OutPut("GMSH used for meshing !!!" << endl);
    }
       else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     } 
     
    //  TriaReMeshGen(Domain);
   
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();  
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
   
  if(TDatabase::ParamDB->WRITE_VTK)
   { mkdir(vtkdir, 0777); }   
  
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);

  // fespaces for velocity and pressure  
  GetVelocityAndPressureSpace(coll, BoundCondition, mortarcoll, Velocity_FeSpace,
                              Pressure_FeSpace, &pressure_space_code,
                              TDatabase::ParamDB->VELOCITY_SPACE,
                              TDatabase::ParamDB->PRESSURE_SPACE);
  
  // defaulty inf-sup pressure space will be selected based on the velocity space, so update it in database
  TDatabase::ParamDB->INTERNAL_PRESSURE_SPACE = pressure_space_code;
  velocity_space_code =   TDatabase::ParamDB->VELOCITY_SPACE;
    
  N_U = Velocity_FeSpace->GetN_DegreesOfFreedom();
  N_P = Pressure_FeSpace->GetN_DegreesOfFreedom();    
  N_TotalDOF = 2*N_U + N_P;

  OutPut("Dof Velocity : "<< setw(10) << 2* N_U << endl);
  OutPut("Dof Pressure : "<< setw(10) << N_P << endl);
  OutPut("Total Dof all: "<< setw(10) << N_TotalDOF  << endl);
//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_TotalDOF];
    rhs = new double[N_TotalDOF];

    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
    
    fesp[0] = Velocity_FeSpace;
    fefct[0] = u1;
    fefct[1] = u2;
    aux = new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
			   NSN_ParamFctVelo,
                           NSN_FEValuesVelo,
                           fesp, fefct,
                           NSFctVelo,
                           NSFEFctIndexVelo, NSFEMultiIndexVelo,
                           NSN_ParamsVelo, NSBeginParamVelo);
    
      //  aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {    
     auxerror =  new TAuxParam2D(NSN_FESpacesVelo, NSN_FctVelo, 
			   NSN_ParamFctVelo,
                           NSN_FEValuesVelo,
                           fesp, fefct,
                           NSFctVelo,
                           NSFEFctIndexVelo, NSFEMultiIndexVelo,
                           NSN_ParamsVelo, NSBeginParamVelo);             
    }
    else
    { auxerror = NULL;
    }
    
    
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    NSEType = TDatabase::ParamDB->NSTYPE;
    
    SystemMatrix = new TSystemNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, GALERKIN, NSEType, DIRECT);
 
    // initilize the system matrix with the functions defined in Example file
    //last two arguments is aux that is used to pass additional fe functions (eg. mesh velocity)
    // otherwise, just pass with NULL 
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, auxerror);

    // assemble the system matrix with given sol and rhs 
    SystemMatrix->Assemble(sol, rhs);
 
    // calculate the residual
    defect = new double[N_TotalDOF];
    memset(defect,0,N_TotalDOF*SizeOfDouble);
    
    SystemMatrix->GetResidual(sol, rhs, defect);
    
    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
      residual =  Ddot(N_TotalDOF, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect);  

      OutPut("Nonlinear iteration step   0");
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);     

//====================================================================== 
//Solve the system
// the nonlinear iteration
//======================================================================
    limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
    Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;
    
    for(j=1;j<=Max_It;j++)
     {   
      // Solve the NSE system
      SystemMatrix->Solve(sol, rhs);
      
      //no nonlinear iteration for Stokes problem  
      if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
       break;
      
      // assemble the system matrix with given aux, sol and rhs 
      SystemMatrix->AssembleNonLinear(sol, rhs);     
   
      // get the residual
      memset(defect,0,N_TotalDOF*SizeOfDouble);
      SystemMatrix->GetResidual(sol, rhs, defect);

    //correction due to L^2_O Pressure space 
     if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
      residual =  Ddot(N_TotalDOF, defect, defect);
      impuls_residual = Ddot(2*N_U, defect, defect); 

      OutPut("nonlinear iteration step " << setw(3) << j);
      OutPut(setw(14) << impuls_residual);
      OutPut(setw(14) << residual-impuls_residual);
      OutPut(setw(14) << sqrt(residual) << endl);
 
      if ((sqrt(residual)<=limit)||(j==Max_It))
       {
        break;
       }//if ((sqrt(residual)<=limit)||(j==Max_It))
       
     } //for(j=1;j<=Max_It;j
    
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
       IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);
    
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput2D(2, 1, 1, 2, Domain);

   Output->AddFEVectFunct(Velocity);
   Output->AddFEFunction(Pressure);
   
//    u1->Interpolate(ExactU1);
//    u2->Interpolate(ExactU2);
//    Pressure->Interpolate(ExactP);
   
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      img++;
     }   
     
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
      SystemMatrix->MeasureErrors(ExactU1, ExactU2, ExactP, u_error, p_error);

       OutPut("L2(u): " <<  sqrt(u_error[0]*u_error[0]+u_error[2]*u_error[2]) << endl);
       OutPut("H1-semi(u): " <<  sqrt(u_error[1]*u_error[1]+u_error[3]*u_error[3]) << endl);
       OutPut("L2(p): " <<  p_error[0] << endl);
       OutPut("H1-semi(p): " <<  p_error[1] << endl); 
     } //if(TDatabase::ParamDB->MEASURE_ERRORS)
     
 // Calculation of drag and lift coefficients in the case of flow past a cylinder
// double cd=0.0,cl=0.0;
//  GetCdCl(u1,u2,Pressure,u1,u2,cd,cl);
// cout<<"cd : "<<cd<<"\t"<<"cl : "<<cl<<"\n";
// 
  CloseFiles();
  
  return 0;
} // end main
