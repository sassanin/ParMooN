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
// Purpose:     main program for scalar equations with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 08.08.2014

// =======================================================================

#include <Domain.h>
#include <Database.h>
#include <SystemTCD2D.h>
#include <SystemTCD2D_ALE.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>
#include <MainUtilities.h>
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
#include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
// #include "../Main_Users/Sashi/TCD_2D/Hemker.h"

int main(int argc, char* argv[])
{
  int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img=1, N_G;
  int N_Active;

  double *sol, *rhs, *oldrhs, t1, t2, errors[5], Linfty;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  char *VtkBaseName;
  const char vtkdir[] = "VTK"; 
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TFESpace2D *Scalar_FeSpace, *fesp[1];
  TFEFunction2D *Scalar_FeFunction;
  TOutput2D *Output;
  TSystemTCD2D *SystemMatrix;  
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = {D00, D10, D01}; 
  
  std::ostringstream os;
  os << " ";   
  
  // ======================================================================
  // set the database values and generate mesh
  // ======================================================================
  // set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based 
  Domain = new TDomain(argv[1]);  
  
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 2;
   OpenFiles();

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
  
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */ 
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
      Domain->ReadGeo(TDatabase::ParamDB->GEOFILE); 
      OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE);
       OutPut("GMSH used for meshing !!!" << endl);
    }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)    //triangle mesh
     {  
      OutPut("Triangle.h used for meshing !!!" << endl);
      TriaReMeshGen(Domain);
     } 
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }
     
#if defined(__HEMKER__) || defined(__BEAM__) 
  TriaReMeshGen(Domain);  
  TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif   
  
  // refine grid up to the coarsest level
  for(i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain->RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain->PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(vtkdir, 0777);

  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
  
  // fespaces for scalar equation 
  Scalar_FeSpace = new TFESpace2D(coll, (char*)"fe space", (char*)"solution space", 
                                          BoundCondition, ORDER, NULL);
   
  N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  N_Active =  Scalar_FeSpace->GetActiveBound();
  OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
  OutPut("dof active   : "<< setw(10) << N_Active << endl);
   
//======================================================================
// construct all finite element functions
//======================================================================
  sol = new double[N_DOF];
  rhs = new double[N_DOF];
  oldrhs = new double[N_DOF];
    
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char*)"sol", (char*)"sol", sol, N_DOF); 
  
  //interpolate the initial value
 Scalar_FeFunction->Interpolate(InitialCondition);

  //======================================================================
  // SystemMatrix construction and solution
  //======================================================================    
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    SystemMatrix = new TSystemTCD2D(Scalar_FeSpace, GALERKIN, DIRECT);
    
    // initilize the system matrix with the functions defined in Example file
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);
       
    // assemble the system matrix with given aux, sol and rhs 
    // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass with NULL 
    SystemMatrix->AssembleMRhs(NULL, sol, rhs);   

  
  //======================================================================
  // produce outout at t=0
  //======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);

//     Scalar_FeFunction->Interpolate(Exact);   
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

     
    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpace;       
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
     for(j=0;j<5;j++)
       errors[j] = 0;
     
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);
     
      olderror = errors[0];
      olderror1 = errors[1]; 
     
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);     
     Linfty=errors[0];
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
       
       
   coll->GetHminHmax(&hmin,&hmax);
   OutPut("h_min : " << hmin << " h_max : " << hmax << endl);

//    TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;

 
  // TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
  //======================================================================
  // time disc loop
  //======================================================================    
  // parameters for time stepping scheme
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;

   // time loop starts
   while(TDatabase::TimeDB->CURRENTTIME < end_time)
   {
    m++;
    TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

    for(l=0; l<N_SubSteps; l++) // sub steps of fractional step theta
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

      
      //copy rhs to oldrhs
      memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 

      // unless the stiffness matrix or rhs change in time, it is enough to 
      // assemble only once at the begning
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
      {
        SystemMatrix->AssembleARhs(NULL, sol, rhs);

        // M:= M + (tau*THETA1)*A
        // rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs +[M-(tau*THETA2)A]*oldsol
        // note! sol contains only the previous time step value, so just pass 
        // sol for oldsol
        SystemMatrix->AssembleSystMat(oldrhs, sol, rhs, sol);
        ConvectionFirstTime = FALSE;
      }
     
      // solve the system matrix 
      SystemMatrix->Solve(sol, rhs);
    
      // restore the mass matrix for the next time step    
      // unless the stiffness matrix or rhs change in time, it is not necessary to assemble the system matrix in every time step
      if(UpdateStiffnessMat || UpdateRhs)
      {
        SystemMatrix->RestoreMassMat();
      }   
    } // for(l=0;l<N_SubSteps;l++) 
    
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

    //======================================================================
    // measure errors to known solution
    //======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
    {
      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors, BilinearCoeffs, aux, 1, fesp, errors);


      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
      OutPut(" L2: " << errors[0]);
      OutPut(" H1-semi: " << errors[1] << endl);

      errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      olderror = errors[0];
      OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

      errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
      OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
      olderror1 = errors[1];   
      
      
      if(Linfty<errors[0])
      Linfty=errors[0];

      OutPut( "Linfty " << Linfty << endl);      
    } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  
  
   } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

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
