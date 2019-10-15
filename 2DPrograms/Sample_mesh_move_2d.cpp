#include <Domain.h>
#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <FESpace2D.h>
#include <SquareStructure2D.h>
#include <Structure2D.h>
#include <Output2D.h>
#include <LinAlg.h>
#include <Assemble2D.h>
#include <DirectSolver.h>
#include <MainUtilities.h>
#include <DeformMesh2D.h>
#include <LinAlg.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <numeric>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>


void BoundCondition_Y(int i, double t, BoundCond &cond)
{

      cond = DIRICHLET;

}

void BoundCondition_X(int i, double t, BoundCond &cond)
{
  if(i == 1 || i ==3 )
    cond = DIRICHLET;
  else
  {
      cond = NEUMANN;
  } 
}

void BilinearCoeffs(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  double angle = 0, v1, v2;
  int i;
  double *coeff;

  //v1 = cos(angle);
  //v2 = sin(angle);
}
double boundval = 0;
void BoundValue_X(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0:
      value = 0.;
      break;
    case 2:
      value = 0.;
      break;
    case 1:
      value = 4. * (1-Param)*Param; // boundval *Param;
      break;
    case 3:
      value = 2. * (1-Param)*Param;; // boundval * (1 - Param);
      break;
    }
}



void BoundValue_Y(int BdComp, double Param, double &value)
{
  //cout << "Bd Comp = " << BdComp << ", Param = " << Param << ", value = " << value << endl;
	value = 0;
}

int main (int argc, char** argv)
{

int N_Cells, ORDER, N_U, N_P, N_DOF, N_Active;

	double *sol, *sol_buffer, *rhs, *defect, t1, t2, errors[4], residual, impuls_residual;
	double limit, u_error[4], p_error[2];
	srand(time(NULL));
	TDomain *Domain;
	TDatabase *Database = new TDatabase();
	TFEDatabase2D *FEDatabase = new TFEDatabase2D();


	/* ---------------  Collection of all the cells --------------/
	* -- Connectivity Information , Global to local Node NUmbering -- */
	TCollection *coll;


	/* ----------- Finite element Space Used ----------- /
	*  Order of the Element ,Shape functions , Derivative matices  *
	*  For navier Stokes ,we need two Fespaces , for Velocity and pressure */
	TFESpace2D *fespace, *fesp[1], *ferhs[2];
	/* ------------- Vector Unknown -----------*/
	/* will npot be used for scalar variables ( i.e ) scalar Equation  */
	TFEVectFunct2D *displacements;
	/* -------------- Object to Output Class ----------------*/
	TOutput2D *Output;
	TAuxParam2D *aux;
	/* ------------------ Object for Galerkin Discrete Form --------------------- */
	TDiscreteForm2D *discreteform;
	/* --------------------- List of derivatives used  ------------------*/
	/* 00 - the actual unknown variable ( without derivative )
	* 10 - Derivative with respect to the first variable
	* 10 - Derivative with respect to second variable        */
	MultiIndex2D AllDerivatives[3] = {D00, D10, D01};

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
	// Enigma ----- //
	// Database->CheckParameterConsistencyNSE();
	Database->WriteParamDB(argv[0]);

	/* include the mesh from a meshgenerator, for a standard mesh use the build-in function */
	// standard mesh
	GEO = TDatabase::ParamDB->GEOFILE;
	Domain->Init(NULL, GEO);

	  //  TriaReMeshGen(Domain);
  // refine grid up to the Finest  level
	for (int i = 0; i < TDatabase::ParamDB->UNIFORM_STEPS; i++)
	Domain->RegRefineAll();

	ORDER = TDatabase::ParamDB->ANSATZ_ORDER;

	coll = Domain->GetCollection(It_Finest, 0);
	N_Cells = coll->GetN_Cells();
	OutPut("N_Cells : " << N_Cells << endl);

	cout << " Before FE SPACE "<< endl;
	fespace = new TFESpace2D(coll, (char *)"name", (char *)"description", BoundCondition_X, ORDER, NULL);
	cout << " After FE SPAVCE "<< endl;



        
    N_U = fespace->GetN_DegreesOfFreedom();
    N_Active = fespace->GetActiveBound();
    N_DOF = 2 * N_U;
    cout<<" N- DOF : "<< N_DOF<< endl;
	cout<<" N- Active : "<< N_Active << endl;



    VtkBaseName = TDatabase::ParamDB->VTKBASENAME;
    Output = new TOutput2D(1, 0, 1, 1, Domain);
	sol = new double[N_DOF]();
	
	displacements = new TFEVectFunct2D(fespace, (char *)"disp", (char *)"description", sol, N_U, 2);
	Output->AddFEVectFunct(displacements);
    mkdir(vtkdir, 0777);
    os.seekp(std::ios::beg);
    int output_write_counter = 0;
    os <<  "VTK/"<<VtkBaseName<<".0000" <<output_write_counter <<".vtk" << ends;
    output_write_counter =  output_write_counter + 1;
    Output->WriteVtk(os.str().c_str());



	 double* gridpos_1 = new double[N_DOF]();
	 TFEVectFunct2D* gridCoordinates1 = new TFEVectFunct2D(fespace, (char*)"C", (char*)"C", gridpos_1, N_U, 2);
	 gridCoordinates1 -> GridToData();




	cout << " Before Calling class " << endl;
    // Call the class for mesh movement
    deformMesh2D* mesh_deform = new deformMesh2D(coll, BoundCondition_X, 
	BoundCondition_Y, BoundValue_X,  BoundValue_Y);

	cout << " After  Calling class --------------- Success" << endl;
	

    os.seekp(std::ios::beg);
    os <<  "VTK/"<<VtkBaseName<<".123"  <<".vtk" << ends;
    output_write_counter =  output_write_counter + 1;
    Output->WriteVtk(os.str().c_str());


    return 0;
}




