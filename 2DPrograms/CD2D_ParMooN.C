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
// Purpose:     main program for solving a stationary scalar equation using ParMooN
//
// Author:      Sashikumaar Ganesan
//
// History:     Implementation started on 22.08.2014
// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <SystemCD2D.h>
#include <Output2D.h>
#include <MainUtilities.h>

#include <stdlib.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>

// =======================================================================
// include current example
// =======================================================================
#include "../Examples/CD_2D/Hemker.h" // circle in a channel
// #include "../Examples/CD_2D/SineLaplace.h" // smooth sol in unitsquares
// #include "../Examples/CD_2D/TwoInteriorLayers.h" // smooth sol in unitsquares
// #include "../Examples/CD_2D/furnace.h"  
//#include "../Main_Users/Sashi/TCD_2D/Hemker.h"

// =======================================================================
// main program
class TFEFunction2D;
class TFEFunction2D;
// =======================================================================
int main(int argc, char* argv[])
{
  //  declaration of database, you need this in every program
  int i, ORDER, N_Cells, N_DOF, img=0, Disctype;
  
  double *sol, *rhs, t1, t2, errors[4];
     
  char *VtkBaseName;
     
  TDatabase *Database = new TDatabase();
  TFEDatabase2D *FEDatabase = new TFEDatabase2D(); 
  TCollection *coll;
  TDomain *Domain;
  TFESpace2D *Scalar_FeSpace, *fesp[1];
  TFEFunction2D *Scalar_FeFunction;
  TSystemCD2D *SystemMatrix;
  TOutput2D *Output;
  TAuxParam2D *aux;
  MultiIndex2D AllDerivatives[3] = { D00, D10, D01 };
  
  std::ostringstream os;
  os << " ";     
  
  // set variables' value in TDatabase using argv[1] (*.dat file) 
  Domain = new TDomain(argv[1]);
  
  //set PROBLEM_TYPE to CD if not yet set
  if(TDatabase::ParamDB->PROBLEM_TYPE == 0)
    TDatabase::ParamDB->PROBLEM_TYPE = 1;
  //open OUTFILE, this is where all output is written to (addionally to console)
  OpenFiles();
 
  // write all Parameters to the OUTFILE (not to console) for later reference
  Database->WriteParamDB(argv[0]);
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
    else
     {  
      OutPut("Mesh Type not known, set MESH_TYPE correctly!!!" << endl);
      exit(0);
     }    
     
#if defined(__HEMKER__) || defined(__BEAM__) 
//   TriaReMeshGen(Domain);  
//   TDatabase::ParamDB->UNIFORM_STEPS = 0;
#endif   
          
  // refine grid up to the coarsest level
  for(i=0; i<TDatabase::ParamDB->UNIFORM_STEPS; i++)
    Domain->RegRefineAll();  
  
  // write grid into an Postscript file
  if(TDatabase::ParamDB->WRITE_PS)
    Domain->PS("Domain.ps", It_Finest, 0);
  
  // create output directory, if not already existing
  if(TDatabase::ParamDB->WRITE_VTK)
    mkdir(TDatabase::ParamDB->OUTPUTDIR, 0777);

  //=========================================================================
  // construct all finite element spaces
  //=========================================================================
  ORDER = TDatabase::ParamDB->ANSATZ_ORDER;
  
  // a collection is basically only an array of cells, which is needed to create
  // a finite element space
  coll = Domain->GetCollection(It_Finest, 0);
  // print out some information about the mesh
  N_Cells = coll->GetN_Cells();
  OutPut("N_Cells : " << N_Cells <<endl);
  
  // create fespace for scalar equation
  Scalar_FeSpace = new TFESpace2D(coll, (char*)"name", (char*)"description", BoundCondition, ORDER, NULL);
  // print out some information on the finite element space
  N_DOF = Scalar_FeSpace->GetN_DegreesOfFreedom();
  OutPut("dof all      : "<< setw(10) << N_DOF  << endl);
 
  //======================================================================
  // construct all finite element functions
  //======================================================================
  sol = new double[N_DOF];
  rhs = new double[N_DOF];
  // set solution and right hand side vectors to zero
  memset(sol, 0, N_DOF*SizeOfDouble);
  memset(rhs, 0, N_DOF*SizeOfDouble);

  // create a finite element function
  Scalar_FeFunction = new TFEFunction2D(Scalar_FeSpace, (char*)"C", (char*)"C", sol, N_DOF);
    
  //======================================================================
  // SystemMatrix construction and solution
  //====================================================================== 
  // Disc type: GALERKIN (or) SDFEM  (or) UPWIND (or) GLS (or) SUPG (or) LOCAL_PROJECTION
  // Solver: AMG_SOLVE (or) GMG  (or) DIRECT 
  Disctype = TDatabase::ParamDB->DISCTYPE;
  SystemMatrix = new TSystemCD2D(Scalar_FeSpace, Disctype, TDatabase::ParamDB->SOLVER_TYPE);
  
  // initilize the system matrix with the functions defined in the example
  SystemMatrix->Init(BilinearCoeffs, BoundCondition, BoundValue);

  // assemble the system matrix with given aux, sol and rhs 
  // aux is used to pass  addition fe functions (eg. mesh velocity) that is nedded for assembling,
  // otherwise, just pass with NULL 
  SystemMatrix->Assemble(NULL, sol, rhs);

  
  //Solve the system
    t1 = GetTime();
    SystemMatrix->Solve(sol, rhs);
    t2 = GetTime();
    OutPut( "time for solving: " << t2-t1 << endl);
    
  //======================================================================
  // produce outout
  //======================================================================
  VtkBaseName = TDatabase::ParamDB->VTKBASENAME;    
  Output = new TOutput2D  (1, 1, 0, 0, Domain);
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
   
  //====================================================================== 
  // measure errors to known solution
  // If an exact solution is not known, it is usually set to be zero, so that
  // in such a case here only integrals of the solution are computed.
  //======================================================================    
    if(TDatabase::ParamDB->MEASURE_ERRORS)
     {
      fesp[0] = Scalar_FeSpace;
      aux =  new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL);

      Scalar_FeFunction->GetErrors(Exact, 3, AllDerivatives, 2, L2H1Errors,
                                   BilinearCoeffs, aux, 1, fesp, errors);

      delete aux;

      OutPut( "L2: " << errors[0] << endl);
      OutPut( "H1-semi: " << errors[1] << endl);
      OutPut( "SD: " << errors[2] << endl);

     } // if(TDatabase::ParamDB->MEASURE_ERRORS)
 
  delete [] sol;
  delete [] rhs;
  delete coll;
  CloseFiles();
  return 0;
} // end main
