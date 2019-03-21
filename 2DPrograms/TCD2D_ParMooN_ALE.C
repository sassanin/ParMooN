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

#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
// #include <gridgen.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>
#include <TimeConvDiff2D.h>

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
// #include "../Examples/TCD_2D/exp.h"
// #include "../Examples/TCD_2D/SinCos1.h"
// #include "../Examples_All/TCD_2D/Time3.h"
// #include "../Examples/TCD_2D/exp_0.h"
//    #include "../Examples/TCD_2D/exp_2.h"
// #include "../Examples_All/TCD_2D/exp_1.h"
// #include  "../Examples/TCD_2D/TimeDomianNoSource.h"
#include  "../Examples/TCD_2D/TimeTwoInteriorLayers.h"
// #include  "../TCD_2D/TimeDomianNoSource_error.h"
//  #include  "../TCD_2D/Bonito_errorchebyshev.h"
//#include  "../TCD_2D/Hemker.h"
// #include  "../TCD_2D/TimeDomian_beam.h"
// #include  "../TCD_2D/ALEOsciDrop.h"
// #include  "../TCD_2D/2D_Smmoth_dGmovin.h" 

int main(int argc, char* argv[])
{
 
// ======================================================================
// variable declaration
// ======================================================================
  int i, j, l, m, N_SubSteps, ORDER, N_Cells, N_DOF, img=1, N_G;
  int N_Active, N_MovVert[5];
  int **IsoCellEdgeNos, N_Data=1;
  
  double *sol, *rhs, *oldrhs, *MeshSol, t1, t2, errors[5], *Iso_refX, Params[10], InitialMass;
  double tau, end_time, *defect, olderror, olderror1, hmin, hmax, CurrTime, tau1, oldt;
  
  bool UpdateStiffnessMat, UpdateRhs,  ConvectionFirstTime;
  
  TDomain *Domain;
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TFESpace2D *Scalar_FeSpace, *fesp[1], *GridFESpace;
  TFEFunction2D *Scalar_FeFunction;
  TOutput2D *Output;
  TSystemTCD2D_ALE *SystemMatrix_ALE; 
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  TFEVectFunct2D *MeshVelocity; 
  TIsoBoundEdge **Free_Joint, *IsoJoint;
  TVertex **MovBoundVert;
  TBaseCell **Free_Cells;
  TAuxParam2D *Aux_ALE;
  TFEFunction2D  *fefct[4];
  
  const char vtkdir[] = "VTK"; 
  const char BDdir[] = "BDData";
  
  char *PsBaseName, *VtkBaseName, *GEO;
  char Name[] = "name";
  char Description[] = "description";
  char CString[] = "C";
  char WString[] = "w";  
  double Linfty;

  std::ostringstream os;
  os << " ";   
    
  mkdir(vtkdir, 0777);
  mkdir(BDdir, 0777);
  
// ======================================================================
// set the database values and generate mesh
// ======================================================================    
  /** set variables' value in TDatabase using argv[1] (*.dat file), and generate the MESH based */
  Domain = new TDomain(argv[1]);  
  TDatabase::ParamDB->Domain = Domain;
  
  OpenFiles();
  OutFile.setf(std::ios::scientific);

  Database->WriteParamDB(argv[0]);
  Database->WriteTimeDB();
  ExampleFile();
 
  /* include the mesh from a mesh generator, for a standard mesh use the 
   * build-in function. The GEOFILE describes the boundary of the domain. */ 
  if(TDatabase::ParamDB->MESH_TYPE==0)
    {
     Domain->ReadGeo(TDatabase::ParamDB->GEOFILE); 
     Domain->ConvertQuadToTri(TDatabase::ParamDB->CONVERT_QUAD_TO_TRI);  
    } // ParMooN  build-in Geo mesh
  else if(TDatabase::ParamDB->MESH_TYPE==1)  
     {Domain->GmshGen(TDatabase::ParamDB->GEOFILE); }//gmsh mesh
  else if(TDatabase::ParamDB->MESH_TYPE==2)    //triangle mesh
     {
#if defined(__ALEDROP__) || defined(__HEMKER__) ||  defined(__BEAM__)         
       TriaReMeshGen(Domain); 
#endif         
    } 
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }
        
  
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
  
//   exit(0);
  
//=========================================================================
// construct all finite element spaces
//=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  coll=Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells (space) : " << N_Cells <<endl);
//   exit(0);
  
  // fespaces for scalar equation 
   Scalar_FeSpace  =  new TFESpace2D(coll, Name, Description, BoundCondition, ORDER, NULL);
   
   N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
   N_Active =  Scalar_FeSpace->GetActiveBound();
   OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
        
  // mesh velocity space 
   GridFESpace = new TFESpace2D(coll, Name, WString, GridBoundCondition, 1, NULL);
 
   N_G = GridFESpace->GetN_DegreesOfFreedom();
   OutPut("N_G    : "<< setw(10) << N_G  << endl);     
//======================================================================
// construct all finite element functions
//======================================================================
    sol = new double[N_DOF];
    rhs = new double[N_DOF];
    oldrhs = new double[N_DOF];
    defect = new double[N_DOF];
    
    memset(sol, 0, N_DOF*SizeOfDouble);
    memset(rhs, 0, N_DOF*SizeOfDouble);

    Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, CString, CString, sol, N_DOF); 
    
    //interpolate the initial value
    Scalar_FeFunction->Interpolate(InitialCondition);
         
    MeshSol = new double[2*N_G];
    memset(MeshSol, 0, 2*N_G*SizeOfDouble);
    MeshVelocity = new TFEVectFunct2D(GridFESpace, WString, WString, MeshSol, N_G, 2);         

//======================================================================
// information for moving mesh based on boundary displacement
//====================================================================== 
#if defined(__ALEDROP__) || defined(__HEMKER__) || defined(__BEAM__) 
   IsoCellEdgeNos = new int *[2];
   GetMovingBoundData(coll, N_MovVert, MovBoundVert, Free_Joint, Free_Cells, IsoCellEdgeNos, Iso_refX); 
   OutPut("N_MovVert  : " << N_MovVert[1] <<endl);  
#endif      
//    exit(0);
//======================================================================
// SystemMatrix construction and solution
//======================================================================  
    SystemMatrix_ALE = new TSystemTCD2D_ALE(Scalar_FeSpace, TDatabase::ParamDB->DISCTYPE, DIRECT, GridFESpace, MeshVelocity,
#ifdef __CONSERVATIVEALE__       
                                                      TRUE
#else         
                                                      FALSE
#endif                                                      
    );
    
   fesp[0] = GridFESpace;
   fefct[0] = MeshVelocity->GetComponent(0);
   fefct[1] = MeshVelocity->GetComponent(1);
   fefct[2] = MeshVelocity->GetComponent(0);
   fefct[3] = MeshVelocity->GetComponent(1);
    
   Aux_ALE =  new TAuxParam2D(TimeCDParamsVeloFieldN_FESpaces,
                            TimeCDParamsVeloFieldN_Fct,
                            TimeCDParamsVeloFieldN_ParamFct,
                            TimeCDParamsVeloFieldN_FEValues_ALE,
                            fesp, fefct,
                            TimeCDParamsVeloFieldFct_ALE,
                            TimeCDParamsVeloFieldFEFctIndex_ALE,
                            TimeCDParamsVeloFieldFEMultiIndex_ALE,
                            TimeCDParamsVeloFieldN_Params_ALE,
                            TimeCDParamsVeloFieldBeginParam);  
    
    
    SystemMatrix_ALE->Init(BilinearCoeffs, BoundCondition, BoundValue, GridCoeffs, GridBoundCondition, GridBoundValue, Aux_ALE);
  
    
#if defined(__ALEDROP__) || defined(__HEMKER__) || defined(__BEAM__)   
    SystemMatrix_ALE->AddBoundModifyFunction(ModifyBdCoords, N_MovVert, MovBoundVert, Free_Joint, Iso_refX);
#else
    SystemMatrix_ALE->AddMeshModifyFunction(ModifyCoords);
#endif    
    
    SystemMatrix_ALE->AssembleMRhs(sol, rhs); 
    
// #ifdef __CONSERVATIVEALE__      
    SystemMatrix_ALE->StoreMmat();
/*#else        
    if(TDatabase::TimeDB->TIME_DISC==1)
     {
      ReAssembleM = FALSE;
     } */  
// #endif  

// exit(0);

//======================================================================
// produce outout at t=0
//======================================================================
    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
    Output = new TOutput2D(2, 2, 1, 1, Domain);
   Output->AddFEVectFunct(MeshVelocity);
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
  
    
   TDatabase::TimeDB->TIMESTEPLENGTH *=  hmin;
  
  
  SystemMatrix_ALE->GetMassAndArea(Scalar_FeFunction, Params);   
  InitialMass = Params[0];
  OutPut("C-Mass: " << Params[0] << " C-Mass Diff: " << InitialMass  - Params[0]<< " Relative C-Mass diff: " << ( InitialMass  - Params[0])/InitialMass << " C " << Params[2] << endl);  
       
  // TDatabase::TimeDB->TIMESTEPLENGTH =  hmax;
//======================================================================
// time disc loop
//======================================================================    
   // parameters for time stepping scheme
   oldt = 0;
   m = 0;
   N_SubSteps = GetN_SubSteps();
   end_time = TDatabase::TimeDB->ENDTIME; 

   UpdateStiffnessMat = FALSE; //check BilinearCoeffs in example file
   UpdateRhs = TRUE; //check BilinearCoeffs in example file
   ConvectionFirstTime=TRUE;

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
//         cout<<Ddot(N_DOF, sol, sol)<<endl;   
      //copy rhs to oldrhs
      memcpy(oldrhs, rhs, N_DOF*SizeOfDouble); 
           
      // move the mesh and assemble all matrices dep. on time discretization
      SystemMatrix_ALE->AssembleALEMat(sol, rhs, tau);
      
      //assemble the system Matrix
      SystemMatrix_ALE->AssembleSystMat(oldrhs, sol, rhs, sol);   
//       cout<<Ddot(N_DOF, sol, sol)<<endl;
      // solve the system matrix 
      SystemMatrix_ALE->Solve(sol, rhs);
    
//       cout<<Ddot(N_DOF, sol, sol)<<endl;
//       cout << " tau " << tau << endl;
//       exit(0);
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

       OutPut("time: " << TDatabase::TimeDB->CURRENTTIME);
       OutPut(" min Div W: " << TDatabase::ParamDB->REACTOR_P29);
       OutPut(" max Div W: " << TDatabase::ParamDB->REACTOR_P30);
      OutPut( ", " << errors[0] << ": Linfty " << Linfty << endl);  
      
      SystemMatrix_ALE->GetMassAndArea(Scalar_FeFunction, Params);   
       OutPut("C-Mass: " << Params[0] << " C-Mass Diff: " << InitialMass  - Params[0]<< " Relative C-Mass diff: " << ( InitialMass  - Params[0])/InitialMass << " C " << Params[2] << endl);     
      
     } //  if(TDatabase::ParamDB->MEASURE_ERRORS)
//      while(TDatabase::TimeDB->CURRENTTIME< end_time){
//      if(Linfty<errors[0])
//        Linfty=errors[0];}
//  OutPut( ", " << errors[0] << ": Linfty " << Linfty << endl);  
//       exit(0);
  } // while(TDatabase::TimeDB->CURRENTTIME< end_time)

#ifdef  __TWO_INTERIOR_LAYERS__
   double h, x=0,values[3],y;
   int  bound_points = 201;
   h = 1.0/(bound_points-1);
   
   os.seekp(std::ios::beg);  
   os << "BDData/Boundary."<<N_Data<<".data" << ends;   
   N_Data++;
   std::ofstream dat(os.str().c_str());
    if (!dat)
     {
      cerr << "cannot open file for output" << endl;
      return -1;
     }
    dat << "%% Boundary data created by ParMooN" << endl;
    dat << "%% Current Time :" << TDatabase::TimeDB->CURRENTTIME << endl;   
   
  
    for (i=0;i<bound_points; i++)
     {
      y = i*h;
      Scalar_FeFunction->FindGradient(x,y,values);
//       OutPut("outflow " << level << " " <<  y <<  " " <<  values[0] << endl);
      dat << y << " " <<  values[0]<< endl;
     }
     
     dat.close();
     cout << endl;
     cout << "Boundary wrote output into file " << endl;    
#endif
    
    
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









