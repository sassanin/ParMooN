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
// Purpose:     main program for solving a time-dependent NSE equation in ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 03.09.2014

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SystemTNSE2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
// #include <TimeUtilities.h>
#include <TNSE2D_ParamRout.h>
#include <TimeDiscRout.h>

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
// #include "../Examples/TNSE_2D/DrivenCavity.h" //   in unit square
// #include "../Examples/TNSE_2D/Bsp3.h" // smooth sol in unit square
// #include "../Examples_All/TNSE_2D/Benchmark2.h"  
// #include "../Examples/TNSE_2D/SinCos.h" // smooth sol in unit square
#include "../Examples/TNSE_2D/SinCos2.h" // smooth sol in unit square
// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, j, l, m, N_Cells, ORDER, N_U, N_P, N_L, N_TotalDOF, img=1, pressure_space_code;
  int Max_It, NSEType, velocity_space_code, N_SubSteps, Disctype;
  
  double *sol, *rhs, *oldrhs, *defect, t1, t2, residual, impuls_residual;
  double limit, AllErrors[7], end_time, oldtau, tau;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll, *mortarcoll = NULL;
  TFESpace2D *Velocity_FeSpace, *Pressure_FeSpace, *fesp[2];
  TFEVectFunct2D *Velocity;
  TFEFunction2D *u1, *u2, *Pressure, *fefct[2];
  TOutput2D *Output;
  TSystemTNSE2D *SystemMatrix;
  TAuxParam2D *aux, *NSEaux_error;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
   
#ifdef __PRIVATE__ 
  TFESpace2D *Projection_space;
#endif
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO;
  char UString[] = "u";
  char PString[] = "p";
  char NameString[] = "VMS";
  
  std::ostringstream os;
  os << " ";   
  
  mkdir(vtkdir, 0777);
    
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->CheckParameterConsistencyNSE();
  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh
   GEO = TDatabase::ParamDB->GEOFILE;
   Domain->Init(NULL, GEO);
   
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
   
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  NSEType = TDatabase::ParamDB->NSTYPE;
  Disctype = TDatabase::ParamDB->DISCTYPE;
    
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
    
#ifdef __PRIVATE__
   if (Disctype == VMS_PROJECTION)
    {
      if (TDatabase::ParamDB->VMS_LARGE_VELOCITY_SPACE==0)
        Projection_space = new TFESpace2D(coll,NameString, UString, BoundCondition, DiscP_PSpace, 0, mortarcoll);
      else
        Projection_space = new TFESpace2D(coll,NameString, UString, BoundCondition, DiscP_PSpace, 1, mortarcoll);
      
      N_L = Projection_space->GetN_DegreesOfFreedom();
      OutPut("Dof Projection : " << setw(10) << N_L << endl);      
    } 
#endif     

//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_TotalDOF]; 
    rhs = new double[N_TotalDOF];
    oldrhs = new double[N_TotalDOF];
    
    memset(sol, 0, N_TotalDOF*SizeOfDouble);
    memset(rhs, 0, N_TotalDOF*SizeOfDouble);

    Velocity = new TFEVectFunct2D(Velocity_FeSpace, UString,  UString,  sol, N_U, 2);
    u1 = Velocity->GetComponent(0);
    u2 = Velocity->GetComponent(1);
    Pressure = new TFEFunction2D(Pressure_FeSpace, PString,  PString,  sol+2*N_U, N_P);
    
//  interpolate the initial solution
    u1->Interpolate(InitialU1);
    u2->Interpolate(InitialU2);
    Pressure->Interpolate(InitialP);    

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN 
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT
    SystemMatrix = new TSystemTNSE2D(Velocity_FeSpace, Pressure_FeSpace, Velocity, Pressure, sol, rhs, Disctype, NSEType, DIRECT
#ifdef __PRIVATE__  
                   , Projection_space, NULL, NULL
#endif     
    );

    //define the aux
    fesp[0] = Velocity_FeSpace;

    fefct[0] = u1;
    fefct[1] = u2;  
    
    switch(Disctype)
     {
      // turbulent viscosity must be computed
      case SMAGORINSKY:
      case VMS_PROJECTION:
      case CLASSICAL_LES:
      case GL00_CONVOLUTION:
      case GL00_AUX_PROBLEM:

       aux =  new TAuxParam2D(TimeNSN_FESpacesVelo_GradVelo, TimeNSN_FctVelo_GradVelo,
                        TimeNSN_ParamFctVelo_GradVelo,
                        TimeNSN_FEValuesVelo_GradVelo,
                        fesp, fefct,
                        TimeNSFctVelo_GradVelo,
                        TimeNSFEFctIndexVelo_GradVelo,
                        TimeNSFEMultiIndexVelo_GradVelo,
                        TimeNSN_ParamsVelo_GradVelo,
                        TimeNSBeginParamVelo_GradVelo);

        break;

      default:
        // 2 parameters are needed for assembling (u1_old, u2_old)
        aux =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2, TimeNSN_ParamFct2,
                           TimeNSN_FEValues2,
                           fesp, fefct,
                           TimeNSFct2,
                           TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                           TimeNSN_Params2, TimeNSBeginParam2);  
    }
    
   // aux for calculating the error
   if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
       NSEaux_error =  new TAuxParam2D(TimeNSN_FESpaces2, TimeNSN_Fct2,
                             TimeNSN_ParamFct2,
                             TimeNSN_FEValues2,
                             fesp, fefct,
                             TimeNSFct2,
                             TimeNSFEFctIndex2, TimeNSFEMultiIndex2,
                             TimeNSN_Params2, TimeNSBeginParam2);     
    }
      double hmin, hmax;
         coll->GetHminHmax(&hmin,&hmax);
   OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
//      TDatabase::TimeDB->TIMESTEPLENGTH = hmin;
//       cout<<TDatabase::TimeDB->TIMESTEPLENGTH<<"\n"; 
   
  //======================================================================  

    // initilize the system matrix with the functions defined in Example file
    // last argument is aux that is used to pass additional fe functions (eg. mesh velocity)    
    SystemMatrix->Init(LinCoeffs, BoundCondition, U1BoundValue, U2BoundValue, aux, NSEaux_error);
    
    // assemble M, A matrices and rhs 
    SystemMatrix->Assemble(sol, rhs);
  
//======================================================================
// produce outout
//======================================================================
   VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
   Output = new TOutput2D(2, 2, 1, 1, Domain);

   Output->AddFEVectFunct(Velocity);
   Output->AddFEFunction(Pressure);
   
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
// time disc loop
//======================================================================    
   // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   oldtau = 1.;
   end_time = TDatabase::TimeDB->ENDTIME; 
   limit = TDatabase::ParamDB->SC_NONLIN_RES_NORM_MIN_SADDLE;
   Max_It = TDatabase::ParamDB->SC_NONLIN_MAXIT_SADDLE;        
   memset(AllErrors, 0, 7.*SizeOfDouble);
   
   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {                                               // time cycle
     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;
     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters(1);

      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
   
      OutPut(endl << "CURRENT TIME: ");
      OutPut(TDatabase::TimeDB->CURRENTTIME << endl);          
       
      //copy sol, rhs to olssol, oldrhs
      memcpy(oldrhs, rhs, N_TotalDOF*SizeOfDouble);        
    
      // assemble only rhs, nonlinear matrix for NSE will be assemble in fixed point iteration
      // not needed if rhs is not time-dependent
      if(m!=1)
       { SystemMatrix->AssembleRhs(sol, rhs); }
      else
       { SystemMatrix->Assemble(sol, rhs);  }
       
      //scale B matices and assemble NSE-rhs based on the \theta time stepping scheme 
      SystemMatrix->AssembleSystMat(tau/oldtau, oldrhs, rhs, sol); 
      oldtau = tau;
         
      // calculate the residual
      defect = new double[N_TotalDOF];
      memset(defect,0,N_TotalDOF*SizeOfDouble);

      SystemMatrix->GetTNSEResidual(sol, defect);

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
//Nonlinear iteration of fixed point type
//======================================================================
     for(j=1;j<=Max_It;j++)
      {   
       // Solve the NSE system
       SystemMatrix->Solve(sol);
      
       if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
        IntoL20FEFunction(sol+2*N_U, N_P, Pressure_FeSpace, velocity_space_code, pressure_space_code);   

       //no nonlinear iteration for Stokes problem  
       if(TDatabase::ParamDB->FLOW_PROBLEM_TYPE==STOKES)
        break;
     
        // restore the mass matrix for the next nonlinear iteration      
        SystemMatrix->RestoreMassMat();     
    
        // assemble the system matrix with given aux, sol and rhs 
        SystemMatrix->AssembleANonLinear(sol, rhs);    
          
        // assemble system mat, S = M + dt\theta_1*A
        SystemMatrix->AssembleSystMatNonLinear();        

        // get the residual
        memset(defect, 0, N_TotalDOF*SizeOfDouble);
        SystemMatrix->GetTNSEResidual(sol, defect);       
             
        //correction due to L^2_O Pressure space 
        if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
         IntoL20Vector2D(defect+2*N_U, N_P, pressure_space_code);
    
        residual =  Ddot(N_TotalDOF, defect, defect);
        impuls_residual = Ddot(2*N_U, defect, defect);
        OutPut("nonlinear iteration step " << setw(3) << j);
        OutPut(setw(14) << impuls_residual);
        OutPut(setw(14) << residual-impuls_residual);
        OutPut(setw(14) << sqrt(residual) << endl); 
	
        if(sqrt(residual)<=limit)
         break;

       } // for(j=1;j<=Max_It;j++)
/*           cout << " test VHM main " << endl;
exit(0);      */  
       // restore the mass matrix for the next time step          
       SystemMatrix->RestoreMassMat();     
       
      } // for(l=0;l<N_SubSteps;
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {   
//        u1->Interpolate(ExactU1);
//        u2->Interpolate(ExactU2);
//        Pressure->Interpolate(ExactP); 
       
      SystemMatrix->MeasureTNSEErrors(ExactU1, ExactU2, ExactP, AllErrors);

       OutPut("L2(u): " <<   AllErrors[0] << endl);
       OutPut("H1-semi(u): " <<  AllErrors[1] << endl);
       OutPut("L2(p): " <<  AllErrors[2] << endl);
       OutPut("H1-semi(p): " <<  AllErrors[3]   << endl); 
       OutPut(AllErrors[4] <<  " l_infty(L2(u)) " <<AllErrors[5] << endl);
       OutPut(TDatabase::TimeDB->CURRENTTIME << " L2(0,t,L2)(u) : " <<   sqrt(AllErrors[6]) << endl); 
      
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

//======================================================================
// produce outout
//======================================================================     
     if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)      
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
                
    } // while(TDatabase::TimeDB->CURRENTTIME< e

//======================================================================
// produce final outout
//======================================================================
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

  CloseFiles();
  
  return 0;
} // end main
