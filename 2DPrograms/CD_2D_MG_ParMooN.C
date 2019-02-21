// =======================================================================
//
// Purpose:     main program with parallel solver  
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 24.07.2009
//              MultiGrid implementation - 03.02.2013
// =======================================================================
#ifdef _MPI
# include "mpi.h"
#endif

#ifdef _OMPONLY
#include <omp.h>
#endif

#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <AuxParam2D.h>
#include <QuadAffin.h>
#include <DirectSolver.h>
#include <Assemble2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <CD2DErrorEstimator.h>

double bound = 0;

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#include <MainUtilities.h>
#include <Upwind.h>
#include <FixedPointIte.h>

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>
#include <LocalProjection.h>
#include <gridgen.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>


#include <ItMethod.h>
#include <FixedPointIte.h>
#include <MultiGridScaIte.h>
#include <FgmresIte.h>
#include <MultiGrid2D.h>
#include <MGLevel2D.h>
#include <Solver.h>

#include <sys/stat.h>
#include <sys/types.h>


#ifdef _MPI
 #include <MeshPartition.h>
 #include <MeshPartition2D.h>
 #include <ParFECommunicator2D.h>
//  #include <MumpsSolver.h>
// #include <ParVectorNSE.h>
#endif


#define AMG 0
#define GMG 1
#define DIRECT 2

// // for external mesh generator
// extern "C"
// {
//   void triangulate(char*, struct triangulateio*,
//                    struct triangulateio*, struct triangulateio*);
// }

// ======================================================================
// include the required example file
// ======================================================================
// #include "../Examples/CD_2D/ReactionDominate.h" 
// #include "../Examples/CD_2D/Plane.h"  // unit square
#include "../Examples/CD_2D/Hemker.h" // circle in a channel
// #include "../Examples/CD_2D/TwoBoundaryLayers.h"  // unit square
// #include "../Examples/CD_2D/TwoBoundaryLayersLaplace.h"  // unit square
// #include "../Examples/CD_2D/Smooth.h"  // unit square
// #include "../Examples/CD_2D/TwoInteriorLayers.h"  // unit square
// #include "../Examples/CD_2D/SineLaplaceDiriHom.h" 
// #include "../Examples/CD_2D/SineLaplace.h" 
// #include "../Examples/CD_2D/Hemker1996.h" 
// #include "../Examples/CD_2D/Smooth.h"
// #include "../Examples/CD_2D/Circle02.h"
// #include "../Examples/CD_2D/Circle03.h"
//    #include "../Examples/CD_2D/Terahertz.h"
//    #include "../Examples/CD_2D/Sine.h"
//    #include "../Examples/CD_2D/DiffOptTomography.h"
// ======================================================================
// utilities for main program
// ======================================================================

/** printing scalar at y=0 **/ 
void PrintScalar(TFEFunction2D *Fefunct, int &N_BData)
{
 int i, N;
 
 double dx, x, y, Val[3];
 
 char *VtkBaseName;
 
 VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

 std::ostringstream os;
 os << " ";

 os.seekp(std::ios::beg);
//       os << "surfact"<< i << ".dat" << ends;
  if(N_BData<10) os << "BDData/"<<VtkBaseName<<"Gamma_0000"<<N_BData<<".data" << ends;
  else if(N_BData<100) os <<"BDData/"<<VtkBaseName<<"Gamma_000"<<N_BData<<".data" << ends;
  else if(N_BData<1000) os <<"BDData/"<<VtkBaseName<<"Gamma_00"<<N_BData<<".data" << ends;
  else if(N_BData<10000) os <<"BDData/"<<VtkBaseName<<"Gamma_0"<<N_BData<<".data" << ends;
  else  os <<"BDData/"<<VtkBaseName<<"Gamma_"<<N_BData<<".data" << ends;
  
 std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }
  dat << "%% Scalar data created by ParMooN" << endl;
//   dat << "%% Current Reference Time :" << TDatabase::TimeDB->CURRENTTIME << endl;
  dat << "%% x, y, value" << endl; 
  
  N = 50;
  
  y=0;  
  x=0.1;
  dx = -60./(double)N;
  
  for(i=0;i<N;i++)
   {
    Fefunct->FindGradient(x,  y, Val);
    dat << x<< " " << y << "  " << Val[0] <<endl;

//     x += dx;
    y += dx;   
   }  
   
  N_BData++;
}




int main(int argc, char* argv[])
{
  TDatabase *Database = new TDatabase();
  
#ifdef _MPI
  const int root = 0;
  int rank, size, len, MaxCpV, MaxSubDomainPerDof, out_rank;
  double t_par1, t_par2;
  char  name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm Comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(name, &len);
  
  TCollection *own_coll;  
  TFESpace2D **OwnScalar_Spaces, *ownscalar_space; 
  TParFECommunicator2D **ParComm;
#endif 
   
  TDomain *Domain = new TDomain();   
  TFEDatabase2D *FEDatabase = new TFEDatabase2D();
  TCollection *coll;
  TBaseCell *cell;
  TFESpace2D *scalar_space, *scalar_space_low, *RefinedOutput_space[2], **Scalar_Spaces;
  TOutput2D *Output;
  TAuxParam2D *aux;
  TSquareStructure2D *sqstructureA, *sqstructureA_low;
  TSquareMatrix2D *sqmatrixA, *sqmatrixA_low, *SQMATRICES[1];
  TSquareMatrix2D **MatricesA, *Tmp_MatricesA;
  TSquareMatrix **sqmatrices = (TSquareMatrix **)SQMATRICES;
  TFEFunction2D *u1, *u2;
  TFESpace2D  *velocity_space;
  TFEVectFunct2D *u;  
  
  TDiscreteForm2D *DiscreteForm;
  TDiscreteForm2D *DiscreteFormUpwind; 
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormGLS;
  TFEFunction2D *C, *C_low, *Refinmed_C, **FEFunctArray;
  TFESpace2D *fesp[2], *ferhs[3];
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  TMGLevel2D *MGLevel, *MGLevel_low;
  TMultiGrid2D *MG;
  TItMethod *itmethod, *prec;  
  
  int i,j,k,l,ret, ORDER, N_DOF, N_DOF_low, N_Active, N_NonActive, *N_DOFs;
  int N_Unknowns, img=1, N_RDOF[0], zerostart;
  int N_LinIter, LEVELS, N_U_DOF, N_BData=0;
  int mg_level, mg_type, n_aux, N_Paramters=1, FirstSolve, low;
  
  double solver_time, hmin, hmax;
  double *sol, *oldsol, *rhs, *defect, *PostSol, *refined_sol[2], *auxsol;
  double *rhs_low, *sol_low,  *RHSs[1], t1, t2, errors[4], **RhsArray;
  double *l2, *h1, *sd;
  double *u_sol, *itmethod_sol, *itmethod_rhs;
  double Parameters[2];

  char *PRM, *GEO;
  char *PsBaseName, *VtkBaseName;
  char CString[] = "c";
  char UString[] = "u";
  char RUString[] = "ref_u";
  char PString[] = "p";
  char ReadinDat [] = "readin.dat";
  char MassMatrix[] = "Mass matrix";
  char Mass[] = "Mass";
  char Name[] = "Temp";
  char UName[] = "u";
  char Description[] = "description";
  char CdString[] = "Conv-Diff";
  char GalString[] = "Galerkin";
  char SDFEMString[] = "SDFEM";
  char UpwString[] = "Upwind";
  const char BDdir[] = "BDData";
  const char vtkdir[] = "VTK";
  char SubID[] = "";
  
  boolean ForPressure=FALSE, ReadVelo=FALSE;

  std::ostringstream os, opts;
  os << " ";
  opts << " ";
  
  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);  
  

//======================================================================
// read parameter file Readin.dat
//======================================================================
  if(argc>=2)
    ret=Domain->ReadParam(argv[1]);
  else
    ret=Domain->ReadParam(ReadinDat);
 
  if(ret==-1)
    exit(-1);

  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  ExampleFile();
  
//======================================================================
// copy read parameters into local variables
//======================================================================
  PRM = TDatabase::ParamDB->BNDFILE;
  GEO = TDatabase::ParamDB->GEOFILE;
  PsBaseName = TDatabase::ParamDB->PSBASENAME;
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;

  // set type of multilevel
  mg_type = TDatabase::ParamDB->SC_MG_TYPE_SCALAR;
  if(TDatabase::ParamDB->SOLVER_TYPE==AMG||
      TDatabase::ParamDB->SOLVER_TYPE==DIRECT)
   mg_type = 0;  
  
  
  if(mg_type)
    mg_level = 0;
  else
    mg_level = -1;
    
  LEVELS = TDatabase::ParamDB->LEVELS;
//   mg_level += LEVELS;

  // initialize multigrid
  if(TDatabase::ParamDB->SOLVER_TYPE == GMG)
  {
    Parameters[0] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_SCALAR;
    Parameters[1] = TDatabase::ParamDB->SC_SMOOTH_DAMP_FACTOR_FINE_SCALAR;
    i=1;
    MG = new TMultiGrid2D(i, N_Paramters, Parameters);
    
    switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
      {
          case 11:
           zerostart = 1;
          break;
          case 16:
           zerostart = 0;
          break;
      }    
    }  
  
  FirstSolve = TDatabase::ParamDB->SC_FIRST_SOLUTION_LEVEL_SCALAR;
  
  
  // finite element spaces on the different levels of the multigrid

  Scalar_Spaces = new TFESpace2D*[LEVELS+1];  
  FEFunctArray = new TFEFunction2D*[LEVELS+1];
  MatricesA = new TSquareMatrix2D*[LEVELS+1];
  
  N_DOFs = new int[LEVELS+1];  
  l2 = new double[LEVELS+1];
  h1 = new double[LEVELS+1];
  sd = new double[LEVELS+1];
  RhsArray = new double*[LEVELS+1];  
  
#ifdef _MPI    
  OwnScalar_Spaces = new TFESpace2D*[LEVELS+1];   
  ParComm = new TParFECommunicator2D*[LEVELS+1];  
  
  out_rank=TDatabase::ParamDB->Par_P0;  
#endif    
//======================================================================
// define discrete form
//======================================================================


//   DiscreteFormGalerkin = new TDiscreteForm2D(CdString, GalString,
//                                   N_Terms, Derivatives, SpacesNumbers,
//                                   N_Matrices, N_Rhs, RowSpace,
//                                   ColumnSpace, RhsSpace, BilinearAssemble,
//                                   BilinearCoeffs, NULL);

  InitializeDiscreteForms_Stationary(DiscreteFormUpwind, DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormGLS,
                                     BilinearCoeffs);
  
//======================================================================
// generate mesh and refun (if needed)
//======================================================================
    Domain->Init(PRM, GEO);
    
#ifdef _MPI  
   if(rank==0)
   {
#endif      
//  write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain_Coarse" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
#ifdef _MPI 
   }
#endif   
   
   
/* generate special mesh for Hemker example */
#ifdef __HMM_1986__
    if(!strcmp(GEO, "InitGrid"))
     if(TDatabase::ParamDB->REACTOR_P25)
         MeshReGen_HemkerResolved(Domain);
     else
         MeshReGen_Hemker(Domain);
#endif

  // refine grid
  for(i=0;i<TDatabase::ParamDB->UNIFORM_STEPS;i++)
	      Domain->RegRefineAll();
//        Domain->AdaptRefineAll();  
 
#ifdef _MPI  
   if(rank==0)
   {
#endif  
      // write grid into an Postscript file
    os.seekp(std::ios::beg);
    os << "Domain" << ".ps" << ends;
    Domain->PS(os.str().c_str(),It_Finest,0);
#ifdef _MPI 
   }
#endif  
    
//=================================================
// STEP >>>>>>> PARTITION GRID USING METIS
//=================================================
#ifdef _MPI
    t1 = MPI_Wtime();
    Partition_Mesh2D(Comm, Domain, MaxCpV); 
    t2 = MPI_Wtime();
                
    if(rank==0)
     printf("Domain Decomposition Complete, time taken for it  is %e\n", (t2-t1));

  os.seekp(std::ios::beg);
  os << rank << ".ps" << ends;
  Domain->PS(os.str().c_str(),It_Finest,0);
  MPI_Barrier(Comm);
 
  MaxSubDomainPerDof = MIN(MaxCpV, size);
#endif        
 
//======================================================================
// include the boundary condition and boundary values from the example file
//======================================================================
  BoundCondFunct2D *BoundaryConditions[1] = { BoundCondition };
  BoundValueFunct2D *BoundaryValues[1] = { BoundValue };

  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
//======================================================================
// loop over all levels (not a multigrid level but for convergence study)
//======================================================================
//  OutPut("Fixed Point Comm " << Comm << endl);
  for(i=0;i<LEVELS;i++)
   {
    mg_level++;
    
#ifdef _MPI  
   if(rank==0)
   {
#endif      
    OutPut("*******************************************************" << endl);
    OutPut("******           GEOMETRY  LEVEL ");
    OutPut(i << "              *******" << endl);
    OutPut("******           MULTIGRID LEVEL ");
    OutPut(mg_level << "              *******" << endl);
    OutPut("*******************************************************" << endl);
//     OutPut("memory before: " << setw(10) << GetMemory() << endl);
#ifdef _MPI 
   }
#endif   
    solver_time = 0.0;
    N_LinIter = 0;  
    
    
    if(i)
     Domain->RegRefineAll();    
#ifdef _MPI
     // removing unwanted cells in the hallo after refinement
     if(i)
      Domain_Crop(Comm, Domain);
#endif
 
    coll=Domain->GetCollection(It_Finest, 0);
//     OutPut( "number of cells: " << coll->GetN_Cells() << endl);
//     coll->GetHminHmax(&hmin,&hmax);
//     OutPut("h_min : " << hmin << " h_max : " << hmax << endl);
//     cout << endl << endl; 

#ifdef _MPI  
    own_coll = Domain->GetOwnCollection(It_Finest, 0, rank);	 
#endif
   
    
//======================================================================
// construct all finite element spaces
//======================================================================
    
    // if multiple discretization multilevel method is used
    // get space for low order disc on finest geo grid
    if(mg_type==1)
    {
     // nonconforming space, including hallo in parallel
     scalar_space_low = new TFESpace2D(coll, Name, Description,BoundCondition, -1, NULL);
     
     Scalar_Spaces[i] = scalar_space_low;
     N_DOF_low = scalar_space_low->GetN_DegreesOfFreedom();
     N_DOFs[i] = N_DOF_low;         
        
#ifdef _MPI
    t1 = MPI_Wtime();

    Scalar_Spaces[i]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[i] = new TParFECommunicator2D(Comm, Scalar_Spaces[i]);
  
    t2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace_low SubDomain dof mapping %e\n", (t2-t1));
      printf("DOF of FeSpace  space_low : %d \n", ParComm[i]->GetN_GlobalDegreesOfFreedom());      
     }
 
     // fespace on own cels
     ownscalar_space = new TFESpace2D(own_coll,  Name, Description, BoundCondition, -1, NULL); 
     OwnScalar_Spaces[i] = ownscalar_space;
#else
     OutPut("N_Dof_low   : "<<  setw(10) << N_DOF_low  << endl);  
    
#endif    
    } // if(mg_type==1)


   
    // get fe space of high order disc on finest geo grid
    if((i>=FirstSolve)||(mg_type==0))
    { 
     scalar_space =  new TFESpace2D(coll, Name, Description, BoundCondition, ORDER, NULL);
      
     Scalar_Spaces[mg_level] = scalar_space;
     N_DOF = scalar_space->GetN_DegreesOfFreedom();
     N_DOFs[mg_level] = N_DOF;  

#ifdef _MPI
    t1 = MPI_Wtime();

    Scalar_Spaces[mg_level]->SetMaxSubDomainPerDof(MaxSubDomainPerDof);
    ParComm[mg_level] = new TParFECommunicator2D(Comm, Scalar_Spaces[mg_level]);
  
    t2 = MPI_Wtime();
    if(rank==out_rank)
     {
      printf("Time taken for FeSpace SubDomain dof mapping %e\n", (t2-t1));
      printf("DOF of FeSpace  space  : %d \n", ParComm[mg_level]->GetN_GlobalDegreesOfFreedom());      
     }
 
     // fespace on own cels
     ownscalar_space = new TFESpace2D(own_coll,  Name, Description, BoundCondition, ORDER, NULL); 
     OwnScalar_Spaces[mg_level] = ownscalar_space;
#else
      OutPut("N_Dof       : "<<  setw(10) << N_DOF  << endl);     
#endif      
     
    }  

    if(ReadVelo)
     {
      velocity_space =  new TFESpace2D(coll, UName, UString, BoundCondition, 2, NULL);
     }
     
//     	double x,y;
// 	int GCell, *BeginIndex, *GlobalNumbers, *DOF;
//       if(rank==1)
//       {
// 	
//        BeginIndex = ownscalar_space->GetBeginIndex();
//         GlobalNumbers = ownscalar_space->GetGlobalNumbers();
//  
//   
// 	 cell = own_coll->GetCell(0);
// 	 DOF = GlobalNumbers+BeginIndex[0];
// 	  GCell = cell->GetGlobalCellNo();
// 	 
// 	  for(j=0;j<4;j++)
// 	  {
// 	    cell->GetVertex(j)->GetCoords(x, y);
// 	    
// 	    
// 	    printf("GCellNo %d DOF: %d X %e and Y %e\n", GCell,  DOF[j],  x, y);  
// 	    
// 	  }
//  
// //        ownscalar_space->GetDOFPosition(71, x, y);
// //        printf("DOF: %d X %e and Y %e\n", 71, x, y);  
//     
//       }
/*       MPI_Finalize();  
      exit(0); */  
          
//======================================================================
// construct all finite element functions
//======================================================================
    if (mg_type==1)
     {
       
//       OutPut("Dof_low       : "<<  setw(10) << N_DOF_low  << endl);

      rhs_low = new double[N_DOF_low];
      memset(rhs_low, 0, N_DOF_low*SizeOfDouble);
      RhsArray[i] = rhs_low;
      
      sol_low = new double[N_DOF_low];
      memset(sol_low, 0, N_DOF_low*SizeOfDouble);
      C_low = new TFEFunction2D(scalar_space_low, UString, UString, sol_low, N_DOF_low);   
      FEFunctArray[i] = C_low;      
     }
       
     
    // high order disc
    N_Unknowns = N_DOF;  
//     OutPut("dof all      : "<< setw(10) << N_Unknowns  << endl);
     
    if ((i>=FirstSolve)||(mg_type==0))
    {     
     sol = new double[N_Unknowns];
     memset(sol, 0, N_Unknowns*SizeOfDouble);      
     rhs = new double[N_Unknowns];
     memset(rhs, 0, N_Unknowns*SizeOfDouble);   
     RhsArray[mg_level] =  rhs;        

     C = new TFEFunction2D(scalar_space, UString, UString, sol, N_DOF);
     FEFunctArray[mg_level] = C;
 
//======================================================================
// produce outout
//======================================================================
   // prepare output, only the concentration will be saved

    Output = new TOutput2D(2, 2, 1, 1,Domain);
    
    Output->AddFEFunction(C);
    
//     C->Interpolate(Exact);       
#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(Comm, img, SubID);
      img++;      
     t_par2 = MPI_Wtime();

     if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
     
//      exit(0);
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
    } // if ((i>=FirstSolve)||(mg_type==0))   
 
 
 
//======================================================================
// allocate memory for all matrices
//======================================================================
    // build matrices for low order disc
    if(mg_type==1)
    {
      // matrix structures
      sqstructureA_low = new TSquareStructure2D(scalar_space_low);
      sqstructureA_low->Sort();
      
      sqmatrixA_low = new TSquareMatrix2D(sqstructureA_low);
      MatricesA[i] = sqmatrixA_low;
    } 

    // build matrices for high order disc
    if ((i>=FirstSolve)||(mg_type==0))
    {
     sqstructureA = new TSquareStructure2D(scalar_space);
     sqstructureA->Sort();

     sqmatrixA = new TSquareMatrix2D(sqstructureA);
     MatricesA[mg_level] = sqmatrixA;
    }


    
//======================================================================    
// prepare multigrid method
//======================================================================   
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {
      case AMG:
      case DIRECT:
        low = mg_level;
        break;

      case GMG:
        // coarsest grid number
        low = 0;
        // determine number of auxiliary arrays
         if ((TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_ALL_SCALAR)
              || (TDatabase::ParamDB->SC_STEP_LENGTH_CONTROL_FINE_SCALAR))
          n_aux=4;
         else
          n_aux=2;
          
        if (mg_type==1)
        {    
#ifdef _MPI  
         MGLevel_low = new TMGLevel2D(i, sqmatrixA_low, rhs_low, sol_low, FEFunctArray[i], ParComm[i], OwnScalar_Spaces[i], n_aux, NULL);
#else
         MGLevel_low = new TMGLevel2D(i, sqmatrixA_low, rhs_low, sol_low, n_aux, NULL);
#endif 
         if (i==0)
          MG->AddLevel(MGLevel_low); 
         else
          MG->ReplaceLevel(i, MGLevel_low);  
        }
        
        // high order disc
        if ((i>=FirstSolve)||(mg_type==0))
        {  
#ifdef _MPI  
         MGLevel = new TMGLevel2D(mg_level, sqmatrixA, rhs, sol, FEFunctArray[mg_level], ParComm[mg_level], OwnScalar_Spaces[mg_level], n_aux, NULL);
#else
         MGLevel = new TMGLevel2D(mg_level, sqmatrixA, rhs, sol, n_aux, NULL);
#endif
         MG->AddLevel(MGLevel);
        }
        break;
    }  // switch(TDatabase::ParamDB->SOLVER_TYPE)
    

    // restrict solution to all grids
    if (TDatabase::ParamDB->SOLVER_TYPE==GMG)
      MG->RestrictToAllGrids();


       
    // if no solution on this grid, continue
    if(FirstSolve>i)
      continue;    
//======================================================================
// assemble all matrices
//======================================================================
      
     // build the discretizations
    for(k=low;k<=mg_level;k++)
    {   
     N_DOF = N_DOFs[k];      
     N_Active = Scalar_Spaces[k]->GetActiveBound();
     N_NonActive = N_DOF - N_Active;
 
     rhs = RhsArray[k];  
     RHSs[0] = rhs;   
    
     memset(rhs, 0, N_DOF*SizeOfDouble);
 
     fesp[0] = Scalar_Spaces[k];;
     ferhs[0] = Scalar_Spaces[k];;
 
     // find discrete form
//      if ((mg_type==1) && (k<i+1))
//       {
//        DiscreteForm = DiscreteFormUpwind;
//        OutPut("UPWIND"<<endl);
//       }
//      else
      switch(TDatabase::ParamDB->DISCTYPE)
      {
       case GALERKIN:
           DiscreteForm = DiscreteFormGalerkin;
       break;

       case SDFEM:
           DiscreteForm = DiscreteFormSDFEM;
       break;

       case GLS:
           DiscreteForm = DiscreteFormGLS;
       break;

       default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
      }

      // initialize matrices
     SQMATRICES[0] = MatricesA[k];
     SQMATRICES[0]->Reset();
     aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);
    
      // assemble
      Assemble2D(1, fesp,
               1, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteForm,
               BoundaryConditions,
               BoundaryValues,
               aux);

#ifdef _MPI  
       ParComm[k]->SetSlaveDofRows(SQMATRICES[0]->GetRowPtr(), SQMATRICES[0]->GetKCol(), SQMATRICES[0]->GetEntries(), rhs);     
#endif
      
      
#ifdef _OPTTOMOGRAPHY_          
   RobinInt(SQMATRICES[0], rhs, BoundCondition);  
#endif      
  

     delete aux;
    } // for(k=low;k<=mg_level;k++)   // endfor, assembling done
 
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_NonActive*SizeOfDouble);
 
//======================================================================
// solve the system
//======================================================================
    N_Unknowns = N_DOF;  
    PostSol = new double[N_Unknowns];
    memset(PostSol, 0, N_Unknowns*SizeOfDouble);     
    oldsol = new double[N_Unknowns]; 
    memset(oldsol, 0, N_Unknowns*SizeOfDouble);     
    defect = new double[N_Unknowns]; 
    memset(defect,0,N_Unknowns*SizeOfDouble);

//     OutPut("norm of solution " <<  sqrt(Ddot(N_Unknowns,sol,sol))  << endl);
//     OutPut("norm of rhs " <<  sqrt(Ddot(N_Unknowns,rhs,rhs))  << endl);
    
//     exit(0);    
    // solve system
    switch(TDatabase::ParamDB->SOLVER_TYPE)
    {  

     case AMG:
        t1 = GetTime();
        Solver(MatricesA[mg_level], RhsArray[mg_level], sol);
        t2 = GetTime();
        OutPut( "time for AMG solving: " << t2-t1 << endl);
        OutPut("solution " << sqrt(Ddot(N_Unknowns,sol,sol)) << endl);
      break;

      case GMG:
        t1 = GetTime();
        // build preconditioner
        switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
        {
         case 5:
           prec = new TMultiGridScaIte(MatVect_Scalar, Defect_Scalar, NULL, 0, N_Unknowns, MG, zerostart);
         break;
         default:
          OutPut("Unknown preconditioner !!!" << endl);
          exit(4711);
        } //switch (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR)
        
        
       switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
       {
        // fixed point iteration
        case 11:
         itmethod = new TFixedPointIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_Unknowns, 1
#ifdef _MPI   
                               , ParComm[mg_level]
#endif                        
                                                );
         
         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
          memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
          memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble); 
         }
         else
         {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
         }
        break;

        case 16:
         // FGMRES
         itmethod = new TFgmresIte(MatVect_Scalar, Defect_Scalar, prec, 0, N_Unknowns, 1);
 
         if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
         {
          itmethod_sol = new double[N_Unknowns];
          itmethod_rhs = new double[N_Unknowns];
          memcpy(itmethod_sol, sol, N_Unknowns*SizeOfDouble);
          memcpy(itmethod_rhs, rhs, N_Unknowns*SizeOfDouble);
         }
         else
         {
          itmethod_sol = sol;
          itmethod_rhs = rhs;
         }
         
        break;

        default:
        OutPut("Unknown solver !!!" << endl);
        exit(4711);
       }  //  switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)    
        
       // solve linear system
       itmethod->Iterate(sqmatrices, NULL, itmethod_sol, itmethod_rhs);  
       
       switch (TDatabase::ParamDB->SC_SOLVER_SCALAR)
        {
          case 11:
          case 16:
            if (TDatabase::ParamDB->SC_PRECONDITIONER_SCALAR == 5)
            {
              memcpy(sol, itmethod_sol, N_Unknowns*SizeOfDouble);
              memcpy(rhs, itmethod_rhs, N_Unknowns*SizeOfDouble);
            }
            break;
        }        
        
        
     case DIRECT:     
      t1 = GetTime();
      DirectSolver(MatricesA[mg_level], rhs, sol);
      t2 = GetTime();
      OutPut( "time for solving: " << t2-t1 << endl);
//     OutPut("solution " << sqrt(Ddot(N_Unknowns,rhs,sol)) << endl);
     break;        
        
    } //  switch(TDatabase::ParamDB->SOLVER_TYPE)
    
    OutPut("MEMORY: " << setw(10) << GetMemory()/(1048576.0));
    OutPut(" MB" << endl); 
        
     
//        MPI_Finalize();  
//       exit(0);   
#ifdef _MPI
     t_par1 = MPI_Wtime();

     if(TDatabase::ParamDB->WRITE_VTK)
      Output->Write_ParVTK(MPI_COMM_WORLD, img, SubID);
        img++;       
     t_par2 = MPI_Wtime();

     if(rank==out_rank)
     printf("Time taken for writing the parvtk file %e\n", (t_par2-t_par1));
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
      
//       os.seekp(std::ios::beg);
//        if(img<10) os <<  "VTK/"<<VtkBaseName<<".0000"<<img<<".plt" << ends;
//          else if(img<100) os <<  "VTK/"<<VtkBaseName<<".000"<<img<<".plt" << ends;
//           else if(img<1000) os <<  "VTK/"<<VtkBaseName<<".00"<<img<<".plt" << ends;
//            else if(img<10000) os <<  "VTK/"<<VtkBaseName<<".0"<<img<<".plt" << ends;
//             else  os <<  "VTK/"<<VtkBaseName<<"."<<img<<".plt" << ends;
//       
//       Output->WriteBinaryPlt(os.str().c_str());

     }
       img++;    
#endif

    // measure errors to known solution
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      C->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                   BilinearCoeffs, aux, 1, fesp, errors);

      delete aux;

      l2[i] = errors[0];
      h1[i] = errors[1];
      sd[i] = errors[2];

      if (i)
       {
        OutPut( "L2: " << errors[0]<< " order " <<  log(l2[i-1]/l2[i])/ln2 << endl);
        OutPut( "H1-semi: " << errors[1] << " order " << log(h1[i-1]/h1[i])/ln2 << endl);
        OutPut( "SD: " << errors[2] << endl);
       }
      else
       {
        OutPut( "L2: " << errors[0] << endl);
        OutPut( "H1-semi: " << errors[1] << endl);
        OutPut( "SD: " << errors[2] << endl);
       }


     } // if(TDatabase::ParamDB->MEASURE_ERRORS)


    OutPut("memory after: " << setw(10) << GetMemory() << endl);
//     exit(0);
    
   } // for(i=0;i<LEVELS;i++)


#ifdef _MPI    
   MPI_Finalize();  
#endif  
  CloseFiles();
  OutPut("used time: " << GetTime() << endl);
  OutPut("used bytes: " << GetMemory() << endl);
  return 0;
}
