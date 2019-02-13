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
// Purpose:     main program for solving a 3D stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 23.01.2015
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <FESpace3D.h>
#include <SystemCD3D.h>
#include <Output3D.h>
#include <MainUtilities.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _MPI
#include "mpi.h"
#include <MeshPartition.h>
#endif

#define AMG 0
#define GMG 1
#define DIRECT 2

double bound = 0;
double timeC = 0;
// =======================================================================
// include current example
// =======================================================================
#include "../Examples/CD_3D/Laplace.h"
// #include "../Examples/CD_3D/HeatChanel.h"
// =======================================================================

// main program
// =======================================================================
int main(int argc, char* argv[])
{
  // ======================================================================
  //  declaration of variables
  // ======================================================================
  int i, N_Cells, ORDER, N_DOF,  img=1;
  int mg_type, mg_level, LEVELS;
  double total_dof = 0;
  
  double *sol, *rhs, **Sol_array, **Rhs_array, t1, t2, errors[4];
  double start_time, end_time;
  double construction, assembling, solving, vtk, total;                         //time calc
  
  TDatabase *Database = new TDatabase();
  const char vtkdir[] = "VTK"; 
  char *PsBaseName, *VtkBaseName;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  char SubID[] = "";
  
  int profiling;
  #ifdef _MPI
  int rank, size, out_rank;
  int MaxCpV, MaxSubDomainPerDof;
  
  MPI_Comm Comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  double time1, time2;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  TDatabase::ParamDB->Comm = Comm;
  
  bool Refine;
  bool metisType[2] = {false,false};
  #endif 
  
  TDomain *Domain;
  TFEDatabase3D *FEDatabase = new TFEDatabase3D(); 
  TCollection *coll;
  TFESpace3D **Scalar_FeSpaces, *fesp[1];
  TFEFunction3D *Scalar_FeFunction, **Scalar_FeFunctions;
  TOutput3D *Output;
  TSystemCD3D *SystemMatrix;
  TAuxParam3D *aux;
  MultiIndex3D AllDerivatives[4] = { D000, D100, D010, D001 };

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
  
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
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
  PsBaseName = TDatabase::ParamDB->BASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    
  /* meshgenerator */
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
      Domain->Init(TDatabase::ParamDB->BNDFILE, TDatabase::ParamDB->GEOFILE); 
      OutPut("PRM-GEO used for meshing !!!" << endl);
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {
       Domain->GmshGen(TDatabase::ParamDB->GEOFILE); 
       OutPut("GMSH used for meshing !!!" << endl);
    }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)   
    {
      Domain->TetrameshGen(TDatabase::ParamDB->GEOFILE); 
      OutPut("Tetgen used for meshing !!!" << endl);
    } //tetgen mesh
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
 
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG|| TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   { 
     mg_type=0; 
     TDatabase::ParamDB->SC_MG_TYPE_SCALAR =  mg_type;
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
    OutPut(LEVELS-1 << "              ======" << endl);
    OutPut("======           MULTIGRID LEVEL ");
    OutPut(mg_level-1 << "              ======" << endl);
    OutPut("=======================================================" << endl);   
#ifdef _MPI 
    }
#endif 
   }
    
  Scalar_FeSpaces = new TFESpace3D*[LEVELS+1];  
  Scalar_FeFunctions = new TFEFunction3D*[LEVELS+1]; 
  Sol_array = new double*[LEVELS+1];
  Rhs_array = new double*[LEVELS+1];

//=========================================================================
// construct all finite element spaces
// loop over all levels (not a multigrid level but for convergence study)  
//=========================================================================
  for(i=0;i<LEVELS;i++)
   {  
     //if(i=LEVELS-1) exit(0);
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
     
     coll = Domain->GetCollection(It_Finest, 0);
     OutPut("ORDER   : " << ORDER <<endl);
     // fespaces for scalar equation 
     Scalar_FeSpaces[i] =  new TFESpace3D(coll, Name, Description, BoundCondition, ORDER);     
      
#ifdef _MPI
     Scalar_FeSpaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
#endif
     
     //multilevel multigrid disc
     if(i==LEVELS-1 && i!=mg_level-1) 
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
     
    if(i==LEVELS-1 && i!=mg_level-1) 
     {  
      N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
      sol = new double[N_DOF];
      rhs = new double[N_DOF];
      Sol_array[mg_level-1] = sol;
      Rhs_array[mg_level-1] = rhs;

      Scalar_FeFunction = new TFEFunction3D(Scalar_FeSpaces[mg_level-1], CString, CString, sol, N_DOF);   
      Scalar_FeFunctions[mg_level-1] = Scalar_FeFunction;
     }//  if(i==LEVELS-1 && i!=mg_level-1) 
     
#ifdef _MPI
     N_Cells = coll->GetN_Cells();
     printf("rank=%d\t N_Cells   : %d\t Dof all   :%d\n",rank,N_Cells,N_DOF);
#endif

   }// for(i=0;i<LEVELS;i++)

#ifndef _MPI   
   N_Cells = coll->GetN_Cells();
   OutPut("N_Cells   : " << N_Cells <<endl);
   OutPut("Dof all   : " << N_DOF  << endl);  
   OutPut(endl);
#endif 
//    OutPut(endl);  
//    N_DOF = Scalar_FeSpaces[mg_level-1]->GetN_DegreesOfFreedom();
//    OutPut("Dof all   : " << N_DOF  << endl);  
// 
//    int *BeginIndex = Scalar_FeSpaces[mg_level-1]->GetBeginIndex();
//    int *GlobalNumbers = Scalar_FeSpaces[mg_level-1]->GetGlobalNumbers();
//    int j, *DOF;
//    TBaseCell *cell;
//    double x, y, z;
//    
//    TVertex *vert;
//    for(i=0;i<N_Cells;i++)
//    {
//     DOF =GlobalNumbers+BeginIndex[i];
//     cell = coll->GetCell(i);
//     cout << endl;
//     for(j=0;j<4;j++)
//      { 
//       cell->GetVertex(j)->GetCoords(x, y, z);
//       cout << i << " x " << x << " y " << y << " z " << z << " L_DOF: " <<  j<< " G_DOF: " <<  DOF[j]<<endl;
//      }
//            exit(0);
//    }
   

//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
    // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
    if(profiling){
#ifdef _MPI
      t1 = MPI_Wtime();
#else
      t1 = GetTime();
#endif
    }
    
    SystemMatrix = new TSystemCD3D(mg_level, Scalar_FeSpaces, Sol_array, Rhs_array,
                                          TDatabase::ParamDB->DISCTYPE, TDatabase::ParamDB->SOLVER_TYPE);

    // initilize the system matrix with the functions defined in Example file
    // last argument aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
    // otherwise, just pass it with NULL 
    SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue, NULL);
    
    if(profiling){
#ifdef _MPI
      t2 = MPI_Wtime();
#else
      t2 = GetTime();
#endif
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
#endif
      construction = t2;
    } 

    if(profiling){
#ifdef _MPI
      t1 = MPI_Wtime();
#else
      t1 = GetTime();
#endif
    }
  
    SystemMatrix->Assemble();
    
    if(profiling){
#ifdef _MPI
      t2 = MPI_Wtime();
#else
      t2 = GetTime();
#endif
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
#endif
      assembling = t2;
    }
    
    //Solve the system
    if(profiling){
#ifdef _MPI
      t1 = MPI_Wtime();
#else
      t1 = GetTime();
#endif
    }
    
    SystemMatrix->Solve(sol, rhs);

    if(profiling){
#ifdef _MPI
      t2 = MPI_Wtime();
#else
      t2 = GetTime();
#endif
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
#endif
      solving = t2;
    }
   
//======================================================================
// produce outout
//======================================================================
    if(profiling){
#ifdef _MPI
      t1 = MPI_Wtime();
#else
      t1 = GetTime();
#endif
    }
    
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput3D(2, 2, 1, 1, Domain);

    Output->AddFEFunction(Scalar_FeFunction);

     //Scalar_FeFunction->Interpolate(Exact);   
#ifdef _MPI
     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;       
#else
    if(TDatabase::ParamDB->WRITE_VTK)
     {
      os.seekp(std::ios::beg);
       if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".vtk" << ends;
         else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".vtk" << ends;
          else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".vtk" << ends;
           else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".vtk" << ends;
            else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".vtk" << ends;
      Output->WriteVtk(os.str().c_str());
      
     }
       img++;    
#endif 

    if(profiling){
#ifdef _MPI
      t2 = MPI_Wtime();
#else
      t2 = GetTime();
#endif
      t2 = t2-t1;
#ifdef _MPI
      MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
#endif
      vtk = t2;
    }       
//====================================================================== 
// measure errors to known solution
//======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      Scalar_FeFunction->Interpolate(Exact); 
      
      fesp[0] = Scalar_FeSpaces[mg_level-1];
      aux =  new TAuxParam3D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
     
      Scalar_FeFunction->GetErrors(Exact, 4, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);
      
      delete aux;
      
#ifdef _MPI
      double reduced_errors[4];
      MPI_Reduce(errors, reduced_errors, 4, MPI_DOUBLE, MPI_SUM, out_rank, Comm);
      if(rank == out_rank){
	OutPut(endl);
	OutPut( "L2: " << sqrt(reduced_errors[0]) << endl);
	OutPut( "H1-semi: " << sqrt(reduced_errors[1]) << endl);
	//OutPut( "SD: " << sqrt(reduced_errors[2]) << endl);
      }
#else
      OutPut(endl);
      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
      //OutPut( "SD: " << errors[2] << endl);
#endif
     } // if(TDatabase::ParamDB->MEASURE_ERRORS)

  
    if(profiling){
#ifdef _MPI
      end_time = MPI_Wtime();
      MPI_Reduce(&start_time, &t1, 1, MPI_DOUBLE, MPI_MIN, out_rank, Comm);
      MPI_Reduce(&end_time,   &t2, 1, MPI_DOUBLE, MPI_MAX, out_rank, Comm);
      total = t2 - t1;
#else
      end_time = GetTime();
      total = end_time - start_time;
#endif
    } 

//======================================================================
// Time profiling Output
//======================================================================  
  if(profiling){
#ifdef _MPI
    int Total_cells, Total_dof;
    MPI_Reduce(&N_Cells, &Total_cells, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    MPI_Reduce(&N_DOF, &Total_dof, 1, MPI_INT, MPI_SUM, out_rank, Comm);
    N_Cells = Total_cells;
    N_DOF   = Total_dof;
    if(rank == out_rank){
#endif
    OutPut("#Levels :: "<<LEVELS<<"  #Uniform refinement :: "<<TDatabase::ParamDB->UNIFORM_STEPS <<"  Order :: "<<TDatabase::ParamDB->ANSATZ_ORDER<<endl);  
    OutPut("Total Cells :: "<<N_Cells<<"     Total_dof :: "<<N_DOF<<endl<<endl);  
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl); 
    OutPut( "Total time taken for initializing System Matrix : " << (construction) << "("<<100*(construction)/(total)<<"%)"<<endl);
    OutPut( "Total time taken for vtk writing : " << (vtk) << "("<<100*(vtk)/(total)<<"%)"<<endl);
    OutPut( "Total time taken for assembling : " << (assembling) << "("<<100*(assembling)/(total)<<"%)"<<endl);
    OutPut( "Total time taken for solving : " << (solving) << "("<<100*(solving)/(total)<<"%)"<<endl);
    OutPut( "Total time taken for communication : " << timeC << "(" <<100*timeC/(total) <<"%)"<< endl);
    OutPut( "Total time taken throughout : " << (total) << endl);
    OutPut("----------------------------------------------------------------------------------------------------------------------"<<endl);
    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl); 
#ifdef _MPI
    }
#endif
  }
  
//   for(i=0;i<mg_level;i++)
//    {
//      total_dof += Scalar_FeSpaces[i]->GetN_DegreesOfFreedom();
//    }
//    printf("total dof over all levels :: %lf\n",total_dof);
  
  CloseFiles();
#ifdef _MPI
  MPI_Finalize(); 
#endif
  
  return 0;
} // end main
