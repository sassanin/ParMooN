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
// Purpose:     main program for time-dependent scalar equation with new kernels of ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.01.15

// =======================================================================
 
#include <Domain.h>
#include <Database.h>
#include <SystemTCD3D.h>
#include <FEDatabase3D.h>
#include <FESpace3D.h>
#include <SquareStructure3D.h>
#include <Structure3D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble3D.h>
#include <Output3D.h>
#include <LinAlg.h>
#include <MainUtilities.h>
#include <TimeDiscRout.h>
// #include <TimeUtilities.h>
#include <Solver.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <tetgen.h>
#include <GridCell.h>
#include <MacroCell.h>
#include <BdPlane.h>
#include <BdSphere.h>
#include <IsoBoundFace.h>
#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
// #include "../Examples/TCD_3D/Sin4.h"
// #include "../Examples/TCD_3D/ConstTSmooth.h"
//  #include "../Examples/TCD_3D/ConstT.h"
#include "../Examples/TCD_3D/amc.h"
// #include "../Examples/TCD_3D/HeatChanel.h"
// #include "../Main_Users/Sashi/TCD_3D/TeraHertzBreast.h"


int main(int argc, char* argv[])
{
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, LEVELS, mg_level, N_Cells, N_DOF, img=1;
  int N_Active, mg_type;
  
  double *sol, *oldsol, *rhs, *oldrhs, t1, t2, errors[5], **Sol_array, **Rhs_array;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax;
  double start_time, stop_time, 
	 start_vtk=0, end_vtk=0, total_vtk=0,
	 start_assembling=0, end_assembling=0, total_assembling=0,
	 start_solve=0, end_solve=0, total_solve=0,
	 start_int=0, end_int=0, total_int=0;

  TDomain *Domain;
  TDatabase *Database = new TDatabase();

  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName, *GEO, *PRM;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  double Linfty=0;
  char SubID[] = "";
  
  int profiling;
#ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2, temp, reduced_errors[4];

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;

  int Refine;
  int metisType[2] = {0,0};
#endif
  
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions;
  TOutput3D *Output;
  TSystemTCD3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = {D000, D100, D010, D001};

  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);
  
// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  
  profiling = TDatabase::ParamDB->timeprofiling;
  if(profiling)
  {
#ifdef _MPI
    start_time = MPI_Wtime();
#else
    start_time = GetTime();
#endif
  }  
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  #ifdef _MPI
  out_rank=TDatabase::ParamDB->Par_P0;
  //out_rank = 0;
  if(rank == out_rank)
  #endif
   {
    Database->WriteParamDB(argv[0]);
    Database->WriteTimeDB();
    ExampleFile();
   }
  
  /* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
  // standard mesh  
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
  /* meshgenerator */
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)   
    {Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE); } //tetgen mesh
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }

  LEVELS = TDatabase::ParamDB->LEVELS;
  if(TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
  {
    TDatabase::ParamDB->UNIFORM_STEPS += (LEVELS-1);
    LEVELS = 1;
  }
  // refine grid up to the coarsest level
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
    Domain->RegRefineAll();
  
#ifdef _MPI
  Domain->GenerateEdgeInfo();
  
  if(profiling)  t1 = MPI_Wtime();
  
  if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("metis type set to %d\n",TDatabase::ParamDB->Par_P2);
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
  //this loop checks if number of cells are sufficient in the coarsest level, such that each 
  //rank get some own cells to work on
  //it does so by changing the metis type first, if not possible then refine and start again
  do
  {
    metisType[TDatabase::ParamDB->Par_P2] = 1;
    Refine = Partition_Mesh3D(Comm, Domain, MaxCpV);	//MaxCpV=maximum cell per vertex
    
    if(metisType[0]*metisType[1] == 1 && Refine)
    {
      metisType[0] = 0;      metisType[1] = 0;
      TDatabase::ParamDB->Par_P2 = 0;
      if(rank == 0)
       {
	  printf("\n----------------------------------------------------------------------------------------\n");
	  printf("Warning :: both metisType used. Now refining the mesh by one step \n");
	  printf("metis type set to 0\n");
	  printf("----------------------------------------------------------------------------------------\n\n");
       }
      Domain->RegRefineAll();
      Domain->GenerateEdgeInfo();
      TDatabase::ParamDB->UNIFORM_STEPS +=1;
    }
  }while(Refine);
  
  if(profiling)  t2 = MPI_Wtime(); 
  
  if(profiling){
    time2 = t2-t1;
    MPI_Reduce(&time2, &time1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
    if(rank == out_rank)
      printf("Time taken for Domain Decomposition is %e\n", time1);
  }

  Domain->GenerateEdgeInfo();
  MaxSubDomainPerDof = MIN(MaxCpV, size);
  TDatabase::ParamDB->WRITE_PS = 0;
#endif   
  
  if(TDatabase::ParamDB->WRITE_PS)
   {
    // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
   }
  

//=========================================================================
// set data for multigrid
//=========================================================================  

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG_SOLVE|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR = mg_type;
    }
  
  if(mg_type)
   {
    mg_level =  LEVELS + 1;
    ORDER = -1;
   }
  else
   {
    mg_level = LEVELS;
    ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
   }
   
  if(TDatabase::ParamDB->SOLVER_TYPE==GMG)
   {
#ifdef _MPI  
    if(rank == out_rank)
    {
#endif 
    OutPut("=======================================================" << endl);
    OutPut("======           GEOMETRY  LEVEL ");
    OutPut(LEVELS << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level << "              ======" << endl);
    OutPut("=======================================================" << endl);   
#ifdef _MPI 
    }
#endif 
   }
    
  Scalar_FeSpaces = new TFESpace3D*[mg_level];  
  Scalar_FeFunctions = new TFEFunction3D*[mg_level]; 
  Sol_array = new double*[mg_level];
  Rhs_array = new double*[mg_level];
 
//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {   
    if(i)
     { Domain->RegRefineAll(); }
     
#ifdef _MPI
     if(rank == out_rank)
       printf("Level :: %d\n\n",i);
     if(i)
     {
       Domain->GenerateEdgeInfo();
       Domain_Crop(Comm, Domain);       // removing unwanted cells in the hallo after refinement 
     }
#endif

     coll=Domain->GetCollection(It_Finest, 0);
     
     // fespaces for scalar equation
     
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     

#ifdef _MPI
     Scalar_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
     
     //multilevel multigrid disc
     if(i==LEVELS-1 && mg_type==1) 
      {
       ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
       Scalar_FeSpaces[mg_level-1] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);
#ifdef _MPI
       Scalar_FeSpaces[mg_level-1]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
      } //  if(i==LEVELS-1 && i!=mg_level-1) 
     
//======================================================================
// construct all finite element functions
//======================================================================
    N_DOF = Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    Sol_array[i] = sol;
    Rhs_array[i] = rhs;   
    
    Scalar_FeFunction  = new TFEFunction3D(Scalar_FeSpaces[i], CString, CString, sol, N_DOF);  
    Scalar_FeFunctions[i] = Scalar_FeFunction;
     
    if(i==LEVELS-1 && mg_type==1) 
     {  
      N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      sol = new double[N_DOF];
      rhs = new double[N_DOF];
      Sol_array[mg_level-1] = sol;
      Rhs_array[mg_level-1] = rhs;

      Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
      Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
     }//   if(i==LEVELS-1 && mg_type==1) 
    
#ifdef _MPI
     N_Cells = coll->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\n",rank,N_Cells,N_DOF);
#endif
   }// for(i=0;i<LEVELS;i++)


   oldrhs = new double[N_DOF];
   oldsol = new double[N_DOF];
   
//    int t;
//    for(t=0;t<N_DOF;t++){
//      oldsol[t]=0.;
//      sol[t]=0.;
//      oldrhs[t]=0.;
//    }
#ifndef _MPI   
   N_Cells = coll->GetN_Cells();
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
   OutPut(endl);
#endif

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    /** interpolate the initial value */
    Scalar_FeFunction->Interpolate(InitialCondition);   
    
    /** Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) SUPG (or) LOCAL_PROJECTION
     *  Solver: AMG_SOLVE (or) GMG  (or) DIRECT */
    if(profiling){
#ifdef _MPI
      start_int = MPI_Wtime();
#else
      start_int = GetTime();
#endif
    }
    
    SystemMatrix = new TSystemTCD3D(mg_level, Scalar_FeSpaces, Sol_array, Rhs_array,
                                              TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE);
    
#ifdef _MPI
    if(rank==0)
#endif
    printf("SystemMatrix constructed\n");
    /** initilize the system matrix with the functions defined in Example file */
    // last argument aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass it with NULL 
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue, NULL);
    if(profiling){
#ifdef _MPI
      end_int = MPI_Wtime();

#else
      end_int = GetTime();
#endif
      total_int = end_int-start_int;
    }
//exit(0);
    /** assemble the system matrix with given aux, sol and rhs 
     *  aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
     *  otherwise, just pass it with NULL  */
    if(profiling){
#ifdef _MPI
      start_assembling = MPI_Wtime();
#else
      start_assembling = GetTime();
#endif
    }
    SystemMatrix->AssembleMRhs(); 
    
    if(profiling){
#ifdef _MPI
      end_assembling = MPI_Wtime();
#else
      end_assembling = GetTime();
#endif
      total_assembling += (end_assembling-start_assembling);
    }
  //exit(0);  
   /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
    memcpy(oldrhs, rhs, N_DOF*SizeOfDouble);     
//======================================================================
// produce outout at t=0
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);

//     Scalar_FeFunction->Interpolate(Exact);
#ifdef _MPI
    if(profiling)	start_vtk = MPI_Wtime();

    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling)	end_vtk = MPI_Wtime();
#else
    if(profiling)	start_vtk = GetTime();
   
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
    if(profiling)	end_vtk = GetTime();
#endif

    
 //Scalar_FeFunction->Interpolate(InitialCondition);   
    
    /** measure errors to known solution */
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpaces[mg_level-1];       
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
      
      for(j=0;j<5;j++)
       errors[j] = 0;

      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);

#ifdef _MPI
      MPI_Allreduce(errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);
      for(i=0;i<2;i++)
	errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2: " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi: " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
#endif           
//      Linfty=errors[0];
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)  

//   exit(0);
  
// #ifdef _MPI
//      if(rank == out_rank)
// #endif     
//      cout << "time " << TDatabase::TimeDB->CURRENTTIME << endl;
     
     coll->GetHminHmax(&hmin,&hmax);
 
#ifdef _MPI
     MPI_Allreduce(&hmax, &temp, 1, MPI_DOUBLE, MPI_MAX, Comm);
     hmax = temp;
     temp = 0;
#endif

#ifdef _MPI
     if(rank == out_rank)
#endif
      OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
      
//      TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
  
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
   
//    cout << "init " << Ddot(N_DOF, sol, sol)<< endl;
   /** time loop starts */
   while(TDatabase::TimeDB->CURRENTTIME< end_time)
    {

     m++;
     TDatabase::TimeDB->INTERNAL_STARTTIME = TDatabase::TimeDB->CURRENTTIME;

     for(l=0;l<N_SubSteps;l++) // sub steps of fractional step theta
      {
       SetTimeDiscParameters(1);

#ifdef _MPI
      if(rank == out_rank)
#endif
      if(m==1)
       {
        OutPut("Theta1: " << TDatabase::TimeDB->THETA1<< endl);
        OutPut("Theta2: " << TDatabase::TimeDB->THETA2<< endl);
        OutPut("Theta3: " << TDatabase::TimeDB->THETA3<< endl);
        OutPut("Theta4: " << TDatabase::TimeDB->THETA4<< endl);
       }

      /** copy the sol to old sol */  
      memcpy(oldsol, sol, N_DOF*SizeOfDouble);     
       
      tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
      TDatabase::TimeDB->CURRENTTIME += tau;
      
#ifdef _MPI
      if(rank == out_rank)
#endif
      OutPut(endl << "CURRENT TIME: "<<TDatabase::TimeDB->CURRENTTIME << endl);
      
      /** unless the stiffness matrix or rhs change in time, it is enough to assemble only once at the begning   */
      if(profiling){ 
#ifdef _MPI
	start_assembling = MPI_Wtime();
#else
	start_assembling = GetTime();
#endif
      }
      
      if(UpdateStiffnessMat || UpdateRhs ||  ConvectionFirstTime)
       {  
        if(UpdateRhs)
        { SystemMatrix->AssembleARhs(); }
        else
         { SystemMatrix->AssembleARhs(); }
         
//          for(i=0;i<N_DOF;i++)
// 	   rhs[i]=0;
         
        /**  M:= M + (tau*TDatabase::TimeDB->THETA1)*A
         *   rhs: =(tau*THETA4)*rhs +(tau*THETA3)*oldrhs + [ M - (tau*THETA2)A]*oldsol **/
	SystemMatrix->AssembleSystMat(oldrhs, oldsol, rhs, sol
#ifdef _MPI
                                      , Rhs_array
#endif
                                     );
        /** copy rhs to oldrhs before calling the solver, as rhs will change in multigrid solver */
        memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
        ConvectionFirstTime = FALSE;
       }
      
      if(profiling){
#ifdef _MPI
	end_assembling = MPI_Wtime();
#else
	end_assembling = GetTime();
#endif
	total_assembling += (end_assembling-start_assembling);
      }
      // solve the system matrix 
      if(profiling){
#ifdef _MPI
	start_solve = MPI_Wtime();
#else
	start_solve = GetTime();
#endif
      }

      SystemMatrix->Solve(sol);
      
      if(profiling){
#ifdef _MPI
	end_solve = MPI_Wtime();
#else
	end_solve = GetTime();
#endif
	total_solve += (end_solve-start_solve);
      }
      
      /** restore the mass matrix for the next time step unless the stiffness matrix 
       * or rhs change in time, it is not necessary to assemble the system matrix in every time step */
      if(UpdateStiffnessMat || UpdateRhs)
       {         
        SystemMatrix->RestoreMassMat();
       }
      
     } // for(l=0;l<N_SubSteps;l++) 
//======================================================================
// produce output
//======================================================================
#ifdef _MPI
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = MPI_Wtime();
    }
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling)	end_vtk = MPI_Wtime();
#else
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = GetTime();
    }
   
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
    if(profiling)	end_vtk = GetTime();
#endif

//======================================================================
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {      
      olderror  = errors[0];
      olderror1 = errors[1];
 
 
      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);
 

#ifdef _MPI
      MPI_Allreduce(errors, reduced_errors, 2, MPI_DOUBLE, MPI_SUM, Comm);      
      for(i=0;i<2;i++)
	errors[i] = sqrt(reduced_errors[i]);
      if(rank == out_rank){
	OutPut(endl);
	OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
	OutPut( "L2: " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi: " << sqrt(reduced_errors[1]) << endl);
      }
#else
      OutPut(endl);
      OutPut("time: " << TDatabase::TimeDB->CURRENTTIME<<endl);
      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
#endif 
       
      if(m>1)
      {      
        errors[3] += (errors[0]*errors[0] + olderror * olderror)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	errors[4] += (errors[1]*errors[1] +olderror1 * olderror1)*TDatabase::TimeDB->TIMESTEPLENGTH/2.0;
	
#ifdef _MPI
	if(rank == out_rank)
#endif
	{
	  OutPut(TDatabase::TimeDB->CURRENTTIME <<  " L2(0,T;L2) " << sqrt(errors[3]) << " ");

          OutPut( "L2(0,T;H1) " << sqrt(errors[4]) << endl);
        }
        

        if(Linfty<errors[0])
          Linfty=errors[0];

#ifdef _MPI
        if(rank == out_rank)
#endif
	OutPut( "Linfty " << Linfty << endl);      
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS) 

  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

//======================================================================
// produce final output
//======================================================================
  #ifdef _MPI
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = MPI_Wtime();
    }
  if(m==1 || m % TDatabase::TimeDB->STEPS_PER_IMAGE == 0)
    if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;

    if(profiling){
      end_vtk = MPI_Wtime();
      total_vtk += (end_vtk-start_vtk);
    }
#else
    if(profiling){
      total_vtk += (end_vtk-start_vtk);
      start_vtk = GetTime();
    }
   
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
    if(profiling){
      end_vtk = GetTime();
      total_vtk += (end_vtk-start_vtk);
    }
#endif
    }
  if(profiling){
#ifdef _MPI
    stop_time = MPI_Wtime();
#else
    stop_time = GetTime();
#endif
  }
//======================================================================
// Time profiling Output
//======================================================================  
  if(profiling){
#ifdef _MPI
    
//     int min_Ncells;
//     MPI_Allreduce(&N_Cells, &min_Ncells, 1, MPI_INT, MPI_MIN, Comm);
//     if(min_Ncells == N_Cells)
//     {
//       OutPut( "min NCells : " << N_Cells << endl);
//       OutPut( "corresponding min Ndof : " << N_DOF << endl);
//     }
//     int max_Ncells;
//     MPI_Allreduce(&N_Cells, &max_Ncells, 1, MPI_INT, MPI_MAX, Comm);
//     if(max_Ncells == N_Cells)
//     {
//       OutPut( "max NCells : " << N_Cells << endl);
//       OutPut( "corresponding max Ndof : " << N_DOF << endl);
//     }

    int Total_cells, Total_dof;
    MPI_Reduce(&N_Cells, &Total_cells, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    MPI_Reduce(&N_DOF, &Total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    N_Cells = Total_cells;
    N_DOF   = Total_dof;
    if(rank == out_rank){
#endif
    OutPut(endl<<"#Levels :: "<<LEVELS<<"  #Uniform refinement :: "<<TDatabase::ParamDB->UNIFORM_STEPS <<"  Order :: "<<TDatabase::ParamDB->ANSATZ_ORDER<<endl);
    OutPut("Total Cells :: "<<N_Cells<<"     Total_dof :: "<<N_DOF<<endl<<endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
    OutPut( "Total time taken for initializing System Matrix : " << (total_int) << "("<<100*(total_int)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for vtk writing : " << (total_vtk) << "("<<100*(total_vtk)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for assembling : " << (total_assembling) << "("<<100*(total_assembling)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for solving : " << (total_solve) << "("<<100*(total_solve)/(stop_time-start_time)<<"%)"<<endl);
    OutPut( "Total time taken for communication : " << timeC << "(" <<100*timeC/(stop_time-start_time) <<"%)"<< endl);
    OutPut( "Total time taken throughout : " << (stop_time-start_time) << endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl);
#ifdef _MPI
    }
#endif
  }

  CloseFiles();
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
} // end main
