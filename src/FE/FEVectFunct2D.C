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
// @(#)FEVectFunct2D.C        1.3 07/22/99
// 
// Class:       TFEVectFunct2D
// Purpose:     a function from a finite element space in 2D
//
// Author:      Gunar Matthies (17.01.98)
//
// History:     start of implementation 17.01.98 (Gunar Matthies)
//
//              start of reimplementation 06.08.1998 (GM)
//
//              WriteSol/ReadSol    08.09.09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef __2D__

#ifdef _MPI
# include "mpi.h"
#endif

#include <FEVectFunct2D.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <Database.h>


#include <Joint.h>
#include <BoundEdge.h>
#include <InterfaceJoint.h>

#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <MooNMD_Io.h>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/** constructor with vector initialization */
TFEVectFunct2D::TFEVectFunct2D(TFESpace2D *fespace2D, char *name, 
                             char *description, double *values, 
                             int length, int n_components)
  : TFEFunction2D(fespace2D, name, description, values, length)
{
  N_Components = n_components;
}

/** convert current grid to vector-values FE function */
void TFEVectFunct2D::GridToData()
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  int N_Cells;
  int N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, *xi, *eta;
  double Val[MaxN_BaseFunctions2D];
  double OutVal[MaxN_BaseFunctions2D];
  int *DOF, Index;
  RefTrans2D F_K;
  TRefTrans2D *rt;
  double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_PointsForNodal2D];
  double FunctionalValuesX[MaxN_PointsForNodal2D];
  double FunctionalValuesY[MaxN_PointsForNodal2D];
  double FctVal[4];

  // begin code
  
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  for(i=0;i<N_Cells;i++)
  {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    F_K = Element->GetRefTransID();

    switch(F_K)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        break;
    }
    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta,
                                X, Y, AbsDetjk);

    // for(j=0;j<4;j++)
    //   cout << X[j] << " " << Y[j] << endl;

    nf->GetAllFunctionals(Coll, cell, X, FunctionalValuesX);
    nf->GetAllFunctionals(Coll, cell, Y, FunctionalValuesY);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
    {
      k = DOF[j];
      Values[k] = FunctionalValuesX[j];
      Values[k+Length] = FunctionalValuesY[j];
      // cout << k << " " << Values[k] << " " << Values[k+Length] << endl;
    }
  } // endfor i
}

/** use current data for grid replacement */
void TFEVectFunct2D::DataToGrid()
{
  int i,j,k,l;
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FEId;
  TFE2D *Element;
  TBaseFunct2D *bf;
  int N_Cells, N_Vertices;
  int N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_, N_Points;
  double s, t;
  double Val[MaxN_BaseFunctions2D];
  double OutVal[MaxN_BaseFunctions2D];
  int *DOF, Index;
  RefTrans2D F_K;
  double xi[4], eta[4];
  double X[4], Y[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];

  // begin code
  
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    bf = Element->GetBaseFunct2D();
    N_LocalDOFs = Element->GetN_DOF();

    F_K = Element->GetRefTransID();

    switch(F_K)
    {
      case QuadAffin:
      case QuadBilinear:
      case QuadIsoparametric:
        N_Vertices = 4;
        xi[0]  = -1; xi[1]  =  1; xi[2]  =  1; xi[3]  = -1;
        eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
        X[0] = X[1] = X[2] = X[3] = 0;
        Y[0] = Y[1] = Y[2] = Y[3] = 0;
      break;

      case TriaAffin:
      case TriaIsoparametric:
        N_Vertices = 3;
        xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
        eta[0] = 0; eta[1] = 0; eta[2] = 1;
        X[0] = X[1] = X[2] = 0;
        Y[0] = Y[1] = Y[2] = 0;
      break;
    }

    for(j=0;j<N_Vertices;j++)
    {
      bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
      bf->ChangeBF(Coll, cell, FunctValues[j]);
    }

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
    {
      k = DOF[j];
      s = Values[k];
      t = Values[k+Length];
      for(k=0;k<N_Vertices;k++)
      {
        X[k] += FunctValues[k][j]*s;
        Y[k] += FunctValues[k][j]*t;
      } // endfor k
    } // endfor j

    for(j=0;j<N_Vertices;j++)
    {
      // cout << "j: " << j << " x: " << X[j] << " y: " << Y[j] << endl;
      cell->GetVertex(j)->SetCoords(X[j], Y[j]);
    }
  } // endfor i
}

/** calculate errors to given vector function */
void TFEVectFunct2D::GetDeformationTensorErrors( 
  DoubleFunct2D *Exact, DoubleFunct2D *Exact1,
  int N_Derivatives,
  MultiIndex2D *NeededDerivatives,
  int N_Errors, ErrorMethod2D *ErrorMeth, 
  CoeffFct2D *Coeff, 
  TAuxParam2D *Aux,
  int n_fespaces, TFESpace2D **fespaces,
  double *errors)
{
  int i,j,k,l,n,m, N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_, N_U;
  int Used[N_FEs2D], *N_BaseFunct;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  TFE2D *ele;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux, *aux1, *aux2, *aux3;
  double *Derivatives[2*MaxN_QuadPoints_2D];
  double *ExactVal[2*MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  double **OrigFEValues, *Orig, value, value1;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex;
  double LocError[4], *Values0,*Values1;
  double hK;
  bool *SecondDer;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  SecondDer = new bool[n_fespaces];
  for(i=0;i<n_fespaces;i++)
    SecondDer[i] = FALSE;

  N_Parameters = Aux->GetN_Parameters();
  aux1 = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux1 + j*N_Parameters;

  aux2 = new double [2*MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux2 + j*N_Derivatives;
  
  aux3 = new double [2*MaxN_QuadPoints_2D * 4];
  for(j=0;j<2*MaxN_QuadPoints_2D;j++)
    ExactVal[j] = aux3 + j*4;

  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20]; 
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  fespace = fespaces[0];
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();
  N_U = Length;
  Values0 = Values;
  Values1 = Values+Length;

  for(i=0;i<N_Errors;i++)
    errors[i] = 0.0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = fespaces[0]->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    memset(Used, 0, N_FEs2D*SizeOfInt);
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      Used[CurrentElement] = 1;
    }

    N_LocalUsedElements = 0;
    memset(LocalUsedElements, 0, SizeOfInt*N_FEs2D);
    j = 0;
    for(k=0;k<N_FEs2D;k++)
      if(Used[k])
      {
        LocalUsedElements[j] = (FE2D)k;
        j++;
      }
    N_LocalUsedElements = j;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetjk);

    if(N_Parameters>0)
      Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param); 

    // calculate all needed derivatives of this FE function
    CurrentElement = fespace->GetFE2D(i, cell);
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = GlobalNumbers + BeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues[l] = Values0[DOF[l]];
      FEFunctValues1[l] = Values1[DOF[l]];
    }

    // for all needed derivatives
    for(k=0;k<N_Derivatives;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      // for all quadrature points
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value = 0;
        value1 = 0;
        for(l=0;l<N_;l++)
        {
          value += FEFunctValues[l] * Orig[l];
          value1 += FEFunctValues1[l] * Orig[l];
        } // endfor l
        Derivatives[j][k] = value;
        Derivatives[j+N_Points][k] = value1;
      } // endfor j
    } // endfor k

    // exact value for first component
    for(j=0;j<N_Points;j++)
      Exact(X[j], Y[j], ExactVal[j]);

    // exact value for second component
    for(j=0;j<N_Points;j++)
      Exact1(X[j], Y[j], ExactVal[j+N_Points]);

    if(Coeff)
      Coeff(N_Points, X, Y, Param, AuxArray);      

    ErrorMeth(N_Points, X, Y, AbsDetjk, weights, hK, Derivatives, 
              ExactVal, AuxArray, LocError);

    for(j=0;j<N_Errors;j++)
      errors[j] += LocError[j];

  } // endfor i

  for(j=0;j<N_Errors;j++)
  {
    if (errors[j] > 0)
    errors[j] = sqrt(errors[j]);
  }

  delete aux;
  delete aux1;
  delete aux2;
  delete aux3;
  delete SecondDer;

} // TFEVectFunct2D::GetDeformationTensorErrors

/** calculate L2-nrom of divergence */
double TFEVectFunct2D::GetL2NormDivergence()
{
  int i,j,k,l;
  int N_Cells;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D UsedElements[1], FEid;
  int N_UsedElements = 1;
  TNodalFunctional2D *nf;
  BaseFunct2D BaseFunct, *BaseFuncts;
  int *N_BaseFunct, N_Bf;
  bool SecondDer[1] = { FALSE };
  double diverror, locdiv;
  int N_Points;
  double *xi, *eta, *weights;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetJK[MaxN_QuadPoints_2D];
  RefTrans2D F_K;
  TRefTrans2D *rt;
  double *Values0, *Values1;
  double FEFunctValues0[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *DOF;
  double **OrigFEValuesX, *OrigX, value;
  double **OrigFEValuesY, *OrigY;
  
  diverror = 0.0;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  Values0 = Values;
  Values1 = Values+Length;

  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  BeginIndex = FESpace2D->GetBeginIndex();

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEid = FESpace2D->GetFE2D(i, cell);
    UsedElements[0] = FEid;

    TFEDatabase2D::GetOrig(N_UsedElements, UsedElements, 
                         Coll, cell, SecondDer,
                         N_Points, xi, eta, weights, X, Y, AbsDetJK);

    // calculate all needed derivatives of this FE function
    BaseFunct = BaseFuncts[FEid];
    N_Bf = N_BaseFunct[FEid];

    DOF = GlobalNumbers + BeginIndex[i];
    for(j=0;j<N_Bf;j++)
    {
      k = DOF[j];
      FEFunctValues0[j] = Values0[k];
      FEFunctValues1[j] = Values1[k];
    }

    OrigFEValuesX = TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    OrigFEValuesY = TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    locdiv = 0;
    // for all quadrature points
    for(j=0;j<N_Points;j++)
    {
      OrigX = OrigFEValuesX[j];
      OrigY = OrigFEValuesY[j];
      value = 0;
      for(l=0;l<N_Bf;l++)
      {
        value += FEFunctValues0[l] * OrigX[l] + FEFunctValues1[l] * OrigY[l];
      } // endfor l
      locdiv += AbsDetJK[j]*weights[j]*(value*value);
    } // endfor j

    diverror += locdiv;

  } // endfor i

  diverror = sqrt(diverror);

  return diverror;
} // TFEVectFunct2D::GetL2NormDivergence


/** write the solution into a data file - written by Sashi **/
void TFEVectFunct2D::WriteSol(double t)
{
  int i, N_Joints, N_Cells;

  char *BaseName, Dquot;

  #ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  #endif

  TCollection *Coll;
  TBaseCell *cell;

  Dquot = 34; //  see ASCII Chart
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();
  BaseName = TDatabase::ParamDB->BASENAME;
  char *output_directory = TDatabase::ParamDB->OUTPUTDIR;

  std::ostringstream os;
  os << " ";

  #ifdef _MPI
  OutPut("Writing solution into "<< output_directory << "/" << BaseName << rank
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << output_directory << "/" << BaseName<<rank<<".Sol" << ends;
  #else
  OutPut("Writing solution into "<< output_directory << "/" << BaseName << t
         << ".Sol MooNMD file"<< endl);
  os.seekp(std::ios::beg);
  os << output_directory << "/" << BaseName << t<<".Sol" << ends;
  #endif

  std::ofstream dat(os.str().c_str());

  if (!dat)
   {
    cerr << "cannot open file for output" << endl;
    exit(0);
   }

    dat << "# Solution of the vector "<<Dquot<<Name<<Dquot<<", written by MooNMD"<< endl;
    dat << "# N_Cells, Cell_Type, N_Dim, N_Dof" <<  endl;
    dat <<N_Cells << " " << N_Joints << " " << N_Components << " " << Length<< endl;
    dat <<  endl;

    dat << "# Dof "<< " Nodal Values"<< endl;

    for(i=0;i<Length;i++)
     dat << i << " " << Values[i] << " " << Values[Length + i] << endl;
  dat.close();
}


/** Read the solution from a given data file - written by Sashi **/
void TFEVectFunct2D::ReadSol(char *BaseName)
{
 int i, j, N_Joints, N_Cells, N_cells, N_joints, N_components, length;
 char line[100];

 TCollection *Coll;
 TBaseCell *cell;

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();

  i=0;
  cell =  Coll->GetCell(i);
  N_Joints = cell->GetN_Joints();

  std::ifstream dat(BaseName);
  if (!dat)
   {
    cerr << "cannot open '" <<  BaseName << "' for input" << endl;
    exit(-1);
   }

  dat.getline (line, 99);
  dat.getline (line, 99);
  dat >> N_cells >> N_joints >> N_components >> length;
  dat.getline (line, 99);

  if(N_cells!=N_Cells || N_joints!=N_Joints || N_components!=N_Components || length!=Length )
   {
    OutPut("Given data file does not match with this FE Vector function !"<<endl);
    OutPut(N_cells <<", "<< N_joints<< ", " << N_components<< " , "<< length<< " , " <<endl);
    OutPut(N_Cells <<", "<< N_Joints<< ", " << N_Components<< " , "<< Length<< " , " <<endl);
    exit(0);
   }

  dat.getline (line, 99);
  OutPut("Reading nodal values of the FE Vector function !"<<endl);

  for(i=0;i<Length;i++)
   {
    dat.getline (line, 99);
    dat >> j >> Values[i] >> Values[Length + i];
   }

  dat.close();
}

/** interpolate the old vect value to the new function */
void TFEVectFunct2D::Interpolate(TFEVectFunct2D *OldVectFunct)
{
 int i, j, N_Cells, N_DOFs, N_LocalDOFs, N_Points, N_Edges;
 int *BeginIndex, *GlobalNumbers, *DOF;
 int PolynomialDegree, ApproxOrder;

 double s, *xi, *eta;
 double X[MaxN_PointsForNodal2D], Y[MaxN_PointsForNodal2D];
 double AbsDetjk[MaxN_PointsForNodal2D];
 double PointValues1[MaxN_PointsForNodal2D];
 double PointValues2[MaxN_PointsForNodal2D];
 double FunctionalValues[MaxN_PointsForNodal2D];
 double FctVal[4];
 double val1[3], val2[3];
  
 bool IsIsoparametric;
   
 TCollection *Coll;
 TBaseCell *cell;
 RefTrans2D RefTrans, *RefTransArray;
 FE2D FEId;
 TFE2D *Element;
 TNodalFunctional2D *nf;
 QuadFormula2D QuadFormula; 
 BF2DRefElements RefElement; 
 TJoint *joint;
 JointType jointtype;
 BoundTypes bdtype;
 RefTrans2D F_K;
 TRefTrans2D *rt;

  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
  
  memset(Values, 0, SizeOfDouble*N_Components*N_DOFs);
  RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();
    
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace2D->GetFE2D(i, cell);
    Element = TFEDatabase2D::GetFE2D(FEId);
    nf = Element->GetNodalFunctional2D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
    ApproxOrder = TFEDatabase2D::GetAccuracyFromFE2D(FEId);

    RefElement = Element->GetBaseFunct2D()->GetRefElement();

    switch(RefElement)
    {
      case BFUnitSquare:
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree
                         (3*PolynomialDegree);
        N_Edges = 4;
      break;

      case BFUnitTriangle:
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree
                         (3*PolynomialDegree-1);
        N_Edges = 3;
      break;
    }
    
    RefTrans = RefTransArray[FEId];

    IsIsoparametric = FALSE;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = TRUE;
      }
    } // endif    
 
    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric
 

    switch(RefTrans)
    {
      case QuadAffin:
        rt = TFEDatabase2D::GetRefTrans2D(QuadAffin);
        ((TQuadAffin *)rt)->SetCell(cell);
        F_K = QuadAffin;
        break;
      case QuadBilinear:
        rt = TFEDatabase2D::GetRefTrans2D(QuadBilinear);
        ((TQuadBilinear *)rt)->SetCell(cell);
        F_K = QuadBilinear;
        break;
      case QuadIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(QuadIsoparametric);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        F_K = QuadIsoparametric;
        break;
      case TriaAffin:
        rt = TFEDatabase2D::GetRefTrans2D(TriaAffin);
        ((TTriaAffin *)rt)->SetCell(cell);
        F_K = TriaAffin;
        break;
      case TriaIsoparametric:
        rt = TFEDatabase2D::GetRefTrans2D(TriaIsoparametric);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(ApproxOrder);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        F_K = TriaIsoparametric;
        break;
      default:
        cout << "unknown reftrans id: " << RefTrans << endl;
    }
    TFEDatabase2D::GetOrigFromRef(F_K, N_Points, xi, eta, X, Y, AbsDetjk);  

    for(j=0;j<N_Points;j++)
    { 
     OldVectFunct->FindVectGradient(X[j], Y[j], val1, val2);
     PointValues1[j] = val1[0];
     PointValues2[j] = val2[0];
    }// for(j=0;j<N_Points;j++)


    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues1, FunctionalValues);    

    DOF = GlobalNumbers+BeginIndex[i];
    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];

    nf->GetAllFunctionals(Coll, (TGridCell *)cell, PointValues2, FunctionalValues);     

    for(j=0;j<N_LocalDOFs;j++)
      Values[N_DOFs + DOF[j]] = FunctionalValues[j];
  } // for(i=0;i<N_Cells;i++)  
} // Interpolate
   
   
/** determine the value of a vect function and its first derivatives at
 the given point */
void TFEVectFunct2D::FindVectGradient(double x, double y, double *val1, double *val2)
{
  int i, j, k, N_Cells, N_Found, N_BaseFunct, N_DOFs;
  int *BeginIndex, *GlobalNumbers, *Numbers, N_Edges;
  
  double xv, yv, xi, eta, U1, U2;  
  double *uorig, *uxorig, *uyorig, *uref, *uxiref, *uetaref;
  double u1, u1x, u1y, u2, u2x, u2y;
 
  TBaseCell *cell;
  TCollection *Coll;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;
  JointType jointtype;
  TJoint *joint;
  BoundTypes bdtype;
  BF2DRefElements RefElement; 
 
  bool IsIsoparametric;
  
  N_Found = 0;
  
  val1[0] = 0;
  val1[1] = 0;
  val1[2] = 0;

  val2[0] = 0;
  val2[1] = 0;
  val2[2] = 0;  


  BeginIndex = FESpace2D->GetBeginIndex();
  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  N_DOFs = FESpace2D->GetN_DegreesOfFreedom();
  
  Coll = FESpace2D->GetCollection();
  N_Cells = Coll->GetN_Cells();  

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    
    if(cell->PointInCell(x,y))
    {
      N_Found++;          
      // cout << "point found" << endl;
      FE_ID = FESpace2D->GetFE2D(i, cell);
      FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
      RefTrans = FE_Obj->GetRefTransID();

      // set cell for reference transformation
      TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

      // find local coordinates of the given point
      // mainly used for remeshing, so Affine trans os enough
      TFEDatabase2D::GetRefFromOrig(RefTrans, x, y, xi, eta);
      
      // values are calculated in oldfunction, so isoparam        
      IsIsoparametric = FALSE;           
      if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
      {
       RefElement = FE_Obj->GetBaseFunct2D()->GetRefElement();     

       switch(RefTrans)
       {
         case BFUnitSquare:    
          N_Edges = 4;
         break;

         case BFUnitTriangle:
          N_Edges = 3;
         break;
        }      

       for(j=0;j<N_Edges;j++)
       {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)
          IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
       {
        bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
        if(bdtype != Line)
          IsIsoparametric = TRUE;
       }
       
       if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
        IsIsoparametric = TRUE;
       }
      }// endif 

     if(IsIsoparametric)
      {
       switch(RefElement)
       {
         case BFUnitSquare:
           RefTrans = QuadIsoparametric;
         break;

         case BFUnitTriangle:
           RefTrans = TriaIsoparametric;
         break;
       }
     } // endif IsIsoparametric           
           
      // set cell for reference transformation
      TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
      
      // cout << " xi: " << xi << endl;
      // cout << "eta: " << eta << endl;      
 
      // get base function object
      bf = FE_Obj->GetBaseFunct2D();      
      N_BaseFunct = bf->GetDimension();      

      uorig = new double[N_BaseFunct];
      uxorig = new double[N_BaseFunct];
      uyorig = new double[N_BaseFunct];

      uref = new double[N_BaseFunct];
      uxiref = new double[N_BaseFunct];
      uetaref = new double[N_BaseFunct];
      
      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);          
      
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
                                   uref, uxiref, uetaref, uorig, uxorig, uyorig);     

      u1 = 0;
      u1x = 0;
      u1y = 0;
      
      u2 = 0;
      u2x = 0;
      u2y = 0;

      Numbers = GlobalNumbers + BeginIndex[i];
      for(j=0;j<N_BaseFunct;j++)
      {
        U1 = Values[Numbers[j]];
        u1 += uorig[j]*U1;
        u1x += uxorig[j]*U1;
        u1y += uyorig[j]*U1;

        U2 = Values[N_DOFs + Numbers[j]];        
        u2 += uorig[j]*U2;
        u2x += uxorig[j]*U2;
        u2y += uyorig[j]*U2;        
      }
  
      val1[0] += u1;
      val1[1] += u1x;
      val1[2] += u1y;
      
      val2[0] += u2;
      val2[1] += u2x;
      val2[2] += u2y;     

      delete [] uorig;
      delete [] uxorig;
      delete [] uyorig;
      delete [] uref;
      delete [] uxiref;
      delete [] uetaref;      
      
    } // if(cell->PointInCell(x,y))
  } // for(i=0;i<N_Cells;i+



  if(N_Found>0)
   {
    val1[0] /= (double)N_Found;
    val1[1] /= (double)N_Found;
    val1[2] /= (double)N_Found;
    
    val2[0] /= (double)N_Found;
    val2[1] /= (double)N_Found;
    val2[2] /= (double)N_Found;
   }
  else 
   {
    cout<<"("<<x<<" , " <<y<<" ) Point not found !!!!!"<<endl;
    exit(0);
   }
}


TFEVectFunct2D& TFEVectFunct2D::operator*=(double alpha)
{
  int N_Active = FESpace2D->GetActiveBound();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] *= alpha;
    }
  }
  return *this;
}

TFEVectFunct2D & TFEVectFunct2D::operator+=(const TFEVectFunct2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different fe spaces. Exiting" << endl);
    exit(1);
  }
  if(Length != rhs.Length)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different lengths. Exiting" << endl);
    exit(1);
  }
  if(Values == rhs.Values)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have the same solution vector. This operation would be "
           << "equivalent to a multiplication by 2! Exiting" << endl);
    exit(1);
  }
  if(N_Components != rhs.N_Components)
  {
    OutPut("ERROR: TFEVectFunct2D::operator+=() The two arguments "
           << "have different number of components. You are trying to add a "
           << rhs.N_Components << "-dimensional vector to a " << N_Components
           << "-dimensional vector! Exiting" << endl);
    exit(1);
  }
  int N_Active = FESpace2D->GetActiveBound();
  for(int n=0; n<N_Components; n++)
  {
    for (int i=0; i<N_Active; i++)
    {
      Values[i+n*Length] += rhs.Values[i+n*Length];
    }
  }
  return *this;
}

TFEVectFunct2D & TFEVectFunct2D::operator=(const TFEVectFunct2D & rhs)
{
  if(FESpace2D != rhs.FESpace2D)
  {
    OutPut("ERROR: TFEVectFunct2D::operator=() The two arguments "
           << "have different fe spaces. Exiting" << endl);
    exit(1);
  }
  if(Length != rhs.Length)
  {
    OutPut("ERROR: TFEVectFunct2D::operator=() The two arguments "
           << "have different lengths. Exiting" << endl);
    exit(1);
  }
  if(N_Components != rhs.N_Components)
  {
    OutPut("ERROR: TFEVectFunct2D::operator=() The two arguments "
           << "have different number of components. You are trying to set a "
           << rhs.N_Components << "-dimensional vector to a " << N_Components
           << "-dimensional vector! Exiting" << endl);
    exit(1);
  }
  if(Values == rhs.Values)
  {
    OutPut("WARNING: TFEVectFunct2D::operator=() The two arguments "
           << "have the same solution vector already. No action taken.");
    return *this;
  }
  // copy the values from rhs to *this
  for (int i=0;i<Length*N_Components; i++)
  {
    Values[i] = rhs.Values[i];
  }
  return *this;
}

#endif // #ifdef __2D__
