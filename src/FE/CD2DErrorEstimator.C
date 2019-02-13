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
// @(#)CD2DErrorEstimator.C        1.6 08/17/99
//
// Class:       TCD2DErrorEstimator
//
// Purpose:     compute a posteriori error estimate
//
// Author:      Volker John
//
// History:     04.02.1998 start implementation
//
// =======================================================================

#ifdef __2D__

#include <CD2DErrorEstimator.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <Enumerations.h>
#include <BaseFunct2D.h>
#include <Database.h>
#include <ConvDiff.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifndef __MAC64__
#include <malloc.h>
#endif

/**********************************************************************/
// 0 - gradient indicator
// 1 - H^1 estimator
// 2 - L^2 estimator
// 3 - energy norm + dual norm, Verf"urth 2005
// 4 - energy norm estimator without jumps
// 5 - supg estimator John/Novo, upper estimate
// 6 - supg estimator John/Novo, lower estimate
/**********************************************************************/
#define N_estimators 7

TCD2DErrorEstimator::TCD2DErrorEstimator(int fe_local_estimator,
TFEFunction2D *fe_function2D,
int error_control)
{
  FELocalEstimator = fe_local_estimator;
  FESpace2D =  (TFESpace2D *) fe_function2D->GetFESpace2D();
  Collection = FESpace2D->GetCollection();
  FEFunction2D = fe_function2D;
  ErrorControl = error_control;
}

void TCD2DErrorEstimator::GetErrorEstimate(int N_Derivatives,
MultiIndex2D *NeededDerivatives,
CoeffFct2D *Coeff,
BoundCondFunct2D **BoundaryConds,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Aux,
int n_fespaces,
TFESpace2D **fespaces,
double *eta_K,
double *maximal_local_error,
double *estimated_global_error)
{
    //const int MaxN_BaseFunctions2D_loc = 16;
  const int N_BaseFuncts2D_loc = 5;
  const int MaxN_QuadPoints_1D_loc = 5;
  int i,j,k,l,n,N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_;
  int Used[N_FEs2D],current_estimator, MaxN_BaseFunctions2D_loc;
  int *N_BaseFunct, BaseFunct_loc;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunct;
  TBaseFunct2D *bf;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta,*weights1D, *zeta;
  double xi1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double eta1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *xietaval_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *xideriv_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
  double *etaderiv_ref1D[N_BaseFuncts2D_loc][4][MaxN_QuadPoints_1D_loc];
//  double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
//  double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
//  double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D_loc];
  double *xyval_1D[4];
  double *xderiv_1D[4];
  double *yderiv_1D[4];
  double *X1D[4], *Y1D[4], val[3];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux,*aux1,*aux2;
  double *Derivatives[MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, N_DOF;
  double **OrigFEValues, *Orig, value;
  double *FEFunctValues;
  double *Values,max_loc_err;
  int *GlobalNumbers, *BeginIndex;
  double xc, yc;
  double estimated_global_errors[N_estimators], estimated_local_errors[N_estimators];
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  bool *SecondDer;

  int ee_verbose=1;                               // verbosity

  int memory[3],data_base_memory;
#ifndef __MAC64__  
  struct mallinfo MALLINFO;

  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
  data_base_memory=0;
#endif 
  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  SecondDer = new bool[n_fespaces];

  memset(Used, 0, N_FEs2D*SizeOfInt);
  MaxN_BaseFunctions2D_loc = 0;
  for(i=0;i<n_fespaces;i++)
  {
    fespace = fespaces[i];                        // fe space
    n = fespace->GetN_UsedElements();             // # used finite elements
    UsedElements = fespace->GetUsedElements();    // used finite elements
    for(j=0;j<n;j++)                              // for all finite elements
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
      k = TFEDatabase2D::GetN_BaseFunctFromFE2D(CurrentElement);
      if (k>MaxN_BaseFunctions2D_loc)
      {
	  MaxN_BaseFunctions2D_loc = k;
      }
    }                                             // enfor j
  }                                               // endfor i

  FEFunctValues = new double[MaxN_BaseFunctions2D_loc];
  for (i=0;i<N_BaseFuncts2D_loc;i++)
  {
      for (j=0;j<4;j++)
      {
	  for (k=0;k<MaxN_QuadPoints_1D_loc;k++)
	  {
	      xietaval_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
	      xideriv_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
	      etaderiv_ref1D[i][j][k] = new double[MaxN_BaseFunctions2D_loc];
	  }
      }
  }

  N_UsedElements = 0;                             // compute number of used elements
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  if (N_UsedElements>N_BaseFuncts2D_loc)
  {
      OutPut("CD2DErrorEstimator: too many finite elements " << N_UsedElements << endl);
      OutPut("Increase N_BaseFuncts2D_loc !!!"<<endl);
      exit(4711);
  }

  UsedElements = new FE2D[N_UsedElements];        // store used finite elements 
  j=0;                                            // in array 
  for(i=0;i<N_FEs2D;i++)
  {
    if(Used[i])
    {
	UsedElements[j] = (FE2D)i;
	j++;
    }                                               // endif
  }


  if (ee_verbose>1)
  {
    cout << "estimator number of used elements: " << N_UsedElements << endl;
    for(i=0;i<N_UsedElements;i++)
      cout << "UsedElements[" << i << "]: " << UsedElements[i] << endl;
  }

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################

 for(i=0;i<N_UsedElements;i++)                   // for used finite elements
  {
    CurrentElement = UsedElements[i];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    if (N_Points1D > MaxN_QuadPoints_1D_loc)
    {
	OutPut("CD2DErrorEstimator: too many 1D quadrature points " << N_Points1D << endl);
	OutPut("Increase  MaxN_QuadPoints_1D_loc !!!"<<endl);
	exit(4711);
    }
    BaseFunct = BaseFuncts[CurrentElement];
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunct);// get base functions
    bf2Drefelements = bf->GetRefElement();
    for (j=0;j<N_UsedElements;j++)
    {
	if ((int) BaseFunct == (int) UsedElements[j])
	    break;
    }
    BaseFunct_loc = j;

    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
	// quadrilateral cell
	case BFUnitSquare :                         // edge 0
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][0][j] = zeta[j];
		eta1D[BaseFunct_loc][0][j] = -1;
		bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunct_loc][0][j]);
		bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunct_loc][0][j]);
		bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunct_loc][0][j]);
	    }                                         // edge 1
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][1][j] = 1;
		eta1D[BaseFunct_loc][1][j] = zeta[j];
		bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunct_loc][1][j]);
		bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunct_loc][1][j]);
		bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunct_loc][1][j]);
	    }                                         // edge 2
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][2][j] = -zeta[j];
		eta1D[BaseFunct_loc][2][j] = 1;
		bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunct_loc][2][j]);
		bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunct_loc][2][j]);
		bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunct_loc][2][j]);
	    }                                         // edge 3
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][3][j] = -1;
		eta1D[BaseFunct_loc][3][j] = -zeta[j];
		bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunct_loc][3][j]);
		bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunct_loc][3][j]);
		bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunct_loc][3][j]);
	    }
	    break;
	    
	case BFUnitTriangle :                       // triangular cell
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][0][j] = (zeta[j]+1)/2;
		eta1D[BaseFunct_loc][0][j] = 0;
		bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunct_loc][0][j]);
		bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunct_loc][0][j]);
		bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunct_loc][0][j]);
	    }                                         // edge 1
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][1][j] = (-zeta[j]+1)/2;
		eta1D[BaseFunct_loc][1][j] = (zeta[j]+1)/2;
		bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunct_loc][1][j]);
		bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunct_loc][1][j]);
		bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunct_loc][1][j]);
	    }                                         // edge 2
	    for (j=0;j<N_Points1D;j++)                // for all quadrature points
	    {
		xi1D[BaseFunct_loc][2][j] = 0;
		eta1D[BaseFunct_loc][2][j] = (-zeta[j] +1)/2;
		bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunct_loc][2][j]);
		bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunct_loc][2][j]);
		bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunct_loc][2][j]);
	    }
	    break;
    }
  }                                               // endfor i
 
  for (i=0;i<4;i++)                               // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[i] = new double[N_Points1D];              // coordinates of edge i
    Y1D[i] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[i] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
      xderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
      yderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D_loc];
    }
    xyval_1D[i] = new double[N_Points1D];         // arrays for values in original cell
    xderiv_1D[i] = new double[N_Points1D];
    yderiv_1D[i] = new double[N_Points1D];
  }

  N_Parameters = Aux->GetN_Parameters();          // get number of parameters of equation
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  aux1 = new double [MaxN_QuadPoints_2D*N_Derivatives];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux1 + j*N_Derivatives;

  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux2 + j*20;

  GlobalNumbers = FESpace2D->GetGlobalNumbers();
  BeginIndex = FESpace2D->GetBeginIndex();

  // ########################################################################
  // prepare error estimates
  // ########################################################################

  // all spaces use same Coll
  // Coll = FESpace2D->GetCollection();            // collection of mesh cells
  Coll = Collection;                              // collection of mesh cells
  N_Cells = Coll->GetN_Cells();                   // number of mesh cells
  N_DOF = FEFunction2D->GetLength();              // number of global dof
  Values = FEFunction2D->GetValues();             // values of fe function

  for(i=0;i<N_Cells;i++)                          // do for all mesh cells
  {                                               // on the finest level
    cell=Coll->GetCell(i);
    k=cell->GetN_Edges();                         // # edges
    for(j=0;j<k;j++)                              // for all edges
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);// neighbour cell
      if(neigh) neigh->SetClipBoard(-1);          // set clipboard to -1
    }
    cell->SetClipBoard(-1);
  }                                               // endfor i
  // non finest neighbours of finest cells have clipboard -1

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  for (i=0;i<N_estimators;i++)                               // initialize some quantities
    estimated_global_errors[i]=0.0;
  max_loc_err = 0;
  current_estimator = GetFELocalEstimator();

  // ########################################################################
  // loop over all cells
  // ########################################################################
#ifndef __MAC64__
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    cell = Coll->GetCell(i);                      // next cell
    eta_K[i] = 0.0;                               // initialize local estimate

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i,cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];   // local basis functions
      LocBF[j] = BaseFuncts[CurrentElement];
      SecondDer[j] = TRUE;                        // with 2nd derivative
    }
    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################

    // get reference transformation
    RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);
    if(N_Parameters>0)                            // get parameters of equ.
      Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

    // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D->GetFE2D(i,cell);  // finite element on cell
    BaseFunct = BaseFuncts[CurrentElement];       // basis functions
    N_ = N_BaseFunct[CurrentElement];             // # basis functions
    DOF = GlobalNumbers + BeginIndex[i];          // dof of current mesh cell

    for(l=0;l<N_;l++)
      FEFunctValues[l] = Values[DOF[l]];          // fe values of dofs

    // compute values for all derivatives
    // in all quadrature points
    // in original mesh cell
    for(k=0;k<N_Derivatives;k++)                  // for all derivatives
    {                                             // get values in original cell
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)                     // for all quadrature points
      {
        Orig = OrigFEValues[j];                   // value in original cell
        value = 0;
        for(l=0;l<N_;l++)                         // for all basis functions
          value += FEFunctValues[l] * Orig[l];    // accumulate value of derivative in point j
        Derivatives[j][k] = value;                // for k-th derivative
      }                                           // endfor j
    }                                             // endfor k

    if(Coeff)                                     // get coefficients of pde
      Coeff(N_Points, X, Y, Param, AuxArray);
    // prepare 1D quadrature formula
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

                                                  // update data base
    TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);

    for (j=0;j<N_UsedElements;j++)
    {
	if ((int) BaseFunct == (int) UsedElements[j])
	    break;
    }
    BaseFunct_loc = j;

    N_Edges=cell->GetN_Edges();                   // # edges
    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {                                             // get original coordinates of edge quad. points
      TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunct_loc][j],
        eta1D[BaseFunct_loc][j],
        X1D[j], Y1D[j], AbsDetjk1D[j]);
      for(k=0;k<N_Points1D;k++)                   // get values and derivatives in original cell
      {
        TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunct_loc][j][k],
          eta1D[BaseFunct_loc][j][k],
          TFEDatabase2D::GetBaseFunct2D(BaseFunct),
          Coll, (TGridCell *)cell,
          xietaval_ref1D[BaseFunct_loc][j][k],
          xideriv_ref1D[BaseFunct_loc][j][k],
          etaderiv_ref1D[BaseFunct_loc][j][k],
          xyval_ref1D[j][k],
          xderiv_ref1D[j][k],
          yderiv_ref1D[j][k]);
      }

      for(k=0;k<N_Points1D;k++)                   // for all quadrature points
      {
        val[0]=val[1]=val[2] = 0;
        for(l=0;l<N_;l++)                         // for all basis functions
        {
                                                  // accumulate value of derivative
          val[0] += FEFunctValues[l] * xyval_ref1D[j][k][l];
                                                  // accumulate value of derivative
          val[1] += FEFunctValues[l] * xderiv_ref1D[j][k][l];
                                                  // accumulate value of derivative
          val[2] += FEFunctValues[l] * yderiv_ref1D[j][k][l];
        }                                         // endfor l
        xyval_1D[j][k]= val[0];                   // for k-th
        xderiv_1D[j][k]= val[1];                  // for k-th
        yderiv_1D[j][k]= val[2];                  // for k-th
      }                                           // endfor k
    }                                             // endfor j
    
    if((TDatabase::ParamDB->DISCTYPE == SDFEM)
        || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM))
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
      N_Edges = cell->GetN_Edges();
      for (int ij=0;ij<N_Edges;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if (N_Edges==3)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }

    // estimate local errors
    EstimateCellError_new(fespace,cell,N_Points, X, Y, AbsDetjk, weights, Derivatives,
		      AuxArray,BoundaryConds,BoundaryValues,
		      N_Points1D, zeta, X1D, Y1D, weights1D,
		      xyval_1D, xderiv_1D, yderiv_1D, MaxN_BaseFunctions2D_loc,
		      GlobalNumbers,BeginIndex,DOF,Values,
		      estimated_local_errors);

    for(k=0;k<N_estimators;k++) // update global error estimates
      estimated_global_errors[k]+=estimated_local_errors[k];
    // update maximal local error estimate
    if (estimated_local_errors[current_estimator]>max_loc_err)
      max_loc_err = estimated_local_errors[current_estimator];

    eta_K[i] = estimated_local_errors[current_estimator];
  }                                               // endfor i (loop over cells)
#ifndef __MAC64__
  MALLINFO = mallinfo();
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
  data_base_memory+= memory[2]-memory[1];
#endif
  for(i=1;i<N_estimators;i++)                                // compute global error estimates
    estimated_global_error[i]=sqrt(estimated_global_errors[i]);
  estimated_global_error[0]=estimated_global_errors[0];
  // compute maximal local error estimate
  *maximal_local_error = sqrt(max_loc_err);
  // set memory free
  delete Param[0];
  delete AuxArray[0];
  delete SecondDer;
  delete UsedElements;
  for (i=0;i<4;i++)
  {
    delete X1D[i];
    delete Y1D[i];
    delete AbsDetjk1D[i];
    for (j=0;j<N_Points1D;j++)
    {
      delete xyval_ref1D[i][j];
      delete xderiv_ref1D[i][j];
      delete yderiv_ref1D[i][j];
    }
    delete xyval_1D[i];
    delete xderiv_1D[i];
    delete yderiv_1D[i];
  }
  delete aux1;

  for (i=0;i<N_BaseFuncts2D_loc;i++)
  {
      for (j=0;j<4;j++)
      {
	  for (k=0;k<MaxN_QuadPoints_1D_loc;k++)
	  {
	      delete xietaval_ref1D[i][j][k];
	      delete xideriv_ref1D[i][j][k];
	      delete etaderiv_ref1D[i][j][k];
	  }
      }
  }

  delete FEFunctValues;
#ifndef __MAC64__  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;

  if ((memory[1]- memory[0])!=data_base_memory)
    cout << "WARNING : Error Estimator did not set all memory free !!!" <<  memory[1]- memory[0] << endl;
#endif
}                                                 // TCD2DErrorEstimator::GetErrorEstimate

/*******************************************************************************/
// EstimateCellError_new
// compute the local contribution to the error estimator
// * gradient indicator - estimated_error[0]
// * residual based estimators - estimated_error[i], i>0
// alpha[i] - weight for cell residual
// beta[i]  - weight for edge residuals
// 0 - H^1 estimator
// 1 - L^2 estimator
// 2 - estimator for dual norm of convection, Verf"urth 2005
// 3 - estimator for dual norm of convection, Verf"urth 2005
//     without jumps
// 4 - estimator for SUPG norm, version 1, Novo
// 5 - estimator for SUPG norm, version 2, Novo
/*******************************************************************************/

void  TCD2DErrorEstimator::EstimateCellError_new(TFESpace2D *fespace,
					     TBaseCell *cell,
					     int N_Points,
					     double *X,
					     double *Y,
					     double *AbsDetjk,
					     double *weights,
					     double **Derivatives,
					     double **coeffs,
					     BoundCondFunct2D **BoundaryConds,
					     BoundValueFunct2D **BoundaryValues,
					     int N_Points1D,
					     double *zeta,
					     double *X1D[4],
					     double *Y1D[4],
					     double *weights1D,
					     double *xyval_1D[4],
					     double *xderiv_1D[4],
					     double *yderiv_1D[4],
					     int MaxN_BaseFunctions2D_loc,
					     int *GlobalNumbers,
					     int *BeginIndex,
					     int *DOF,
					     double *Values,
					     double *estimated_local_error)
{
  int i,j,k,l,n,N_Edges,comp,parent_edge,MaxLen1,MaxLen2,MaxLen3,N_child,neigh_edge;
  int chnum1,l_child,child_N_,edge1,neigh_N_,N_Neigh;
  double *deriv, w,e1,e2,*coeff,strong_residual,alpha[N_estimators-1],beta[N_estimators-1],hK,hE,val[3],hE2;
  double estimated_error[N_estimators],t0,t1,nx,ny,x0,x1,y0,y1,neumann_data,jump,absdetjk1D, meas;
  double absdet1D[MaxN_QuadPoints_2D];
  double cell_x0,cell_y0,cell_x1,cell_y1,*xi1DNeigh,*eta1DNeigh, *FEFunctValuesNeigh;
  double *X1DNeigh,*Y1DNeigh,*X1DCell,*Y1DCell;
  double delta_K;
  const int *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpLen3;
  const int  *TmpoEnlE, *TmpEdVer, *TmpECI, *TmpCE, *TmpEdVerParent, *TmpEdVerNeigh;
  TJoint *joint,*parent_joint;
  TBoundEdge *boundedge;
  TBoundComp2D *BoundComp;
  BoundCond Cond0;
  TRefDesc *refdesc,*refdesc_child;
  TBaseCell *neigh, *child, *parent;
  TVertex *ver0,*ver1, *ver2,*ver3;
  FE2D CurrEleNeigh;
  BaseFunct2D BaseFunctNeigh;
  QuadFormula2D QuadFormulaNeigh;
  TQuadFormula2D *qfNeigh;
  QuadFormula1D LineQuadFormulaNeigh;
  TQuadFormula1D *qf1DNeigh;
  TBaseFunct2D *bfNeigh;
  int N_Points1DNeigh,N_PointsNeigh;
  double *weights1DNeigh,*zetaNeigh,*weightsNeigh,*xiNeigh,*etaNeigh;
  TFE2D *eleNeigh;
  RefTrans2D RefTransNeigh;
  BF2DRefElements bf2DrefelementsNeigh;
  double *xietaval_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D];
  double *xideriv_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D];
  double *etaderiv_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D];
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;
  double *xderiv_Cell1D, *yderiv_Cell1D, *xyval_Cell1D;
  int part,edge2neigh;
  int ee_verbose = 1, check_cont, conform_grid=TDatabase::ParamDB->GRID_TYPE;
  TCollection *Coll;

  Coll = fespace->GetCollection();
  // this is for the parameter optimization
  if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
  {
    N_Edges=cell->GetN_Edges();
    for(j=0;j<N_Edges;j++) // loop over all edges of cell
    {
      ver0 = cell->GetVertex(j);
      // clip board contains information on the position of the vertex
      w =  ver0->GetClipBoard();
      // vertex not on the boundary
      if (w<-1e-8)
        continue;
      // component of boundary
      comp = floor(w+1e-8);
      // parameter
      w -= comp;
      // get boundary condition
      //OutPut(comp << " " << w << endl);
      BoundaryConds[0](comp, w, Cond0);
      // Dirichlet
      if (Cond0== DIRICHLET)
      {
        for (i=0;i<N_estimators;i++)
          estimated_error[i] = 0;
	    return;
      }
    }
  }
  k = MaxN_BaseFunctions2D*MaxN_QuadPoints_1D*MaxN_BaseFunctions2D_loc;
  xietaval_refNeigh1D[0][0] = new double[3*k];
  xideriv_refNeigh1D[0][0] = xietaval_refNeigh1D[0][0] + k;
  etaderiv_refNeigh1D[0][0] = xideriv_refNeigh1D[0][0] + k;
  
  for (i=0;i<MaxN_BaseFunctions2D;i++)
  {
      n = i*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D_loc;
      for (j=0;j<MaxN_QuadPoints_1D;j++)
      {
	  l = n + j*MaxN_BaseFunctions2D_loc;
	  xietaval_refNeigh1D[i][j] = xietaval_refNeigh1D[0][0] + l;
	  xideriv_refNeigh1D[i][j] = xideriv_refNeigh1D[0][0] + l;
	  etaderiv_refNeigh1D[i][j] = etaderiv_refNeigh1D[0][0] + l;
      }
  }

  if (ee_verbose>1)
  {
    for (i=0;i<4;i++)
    {
      cout << "1DDD " ;
      for (k=0;k<N_Points1D;k++)
        cout << " x " << X1D[i][k] << " y " << Y1D[i][k];
      cout << endl;
    }

    for (j=0;j<4;j++)
      for (k=0;k<N_Points1D;k++)
        cout << "xx " << xyval_1D[j][k] <<  "dx " << xderiv_1D[j][k] << "dy " <<
          yderiv_1D[j][k] << endl;
  }

  for (i=0;i<N_estimators;i++)                               // initalize local error estimates
    estimated_error[i]=0.0;

  /*************************************************************************/
  /*                                                                       */
  /*   GRADIENT INDICATOR                                                  */
  /*                                                                       */
  /*************************************************************************/
  if (FELocalEstimator==cd_gradient_indicator)    // gradient indicator
  {
    for(i=0;i<N_Points;i++)                       // for all quadrature points
    {
      deriv = Derivatives[i];                     // all derivatives in quadrature points
      w = weights[i]*AbsDetjk[i];

      e1 = deriv[0];                              // x derivative
      e2 = deriv[1];                              // y derivative
      estimated_error[0] += w*(e1*e1+e2*e2);
    }                                             // endfor i
  }

  /*************************************************************************/
  /*                                                                       */
  /*   RESIDUAL BASED EXPLICIT ERROR ESTIMATORS                            */
  /*                                                                       */
  /*************************************************************************/
  if ((FELocalEstimator>0)||(ErrorControl))
  {
    xi1DNeigh = new double[12*N_Points1D];
    eta1DNeigh = xi1DNeigh + N_Points1D;
    X1DNeigh = eta1DNeigh + N_Points1D;
    Y1DNeigh = X1DNeigh + N_Points1D;
    X1DCell = Y1DNeigh + N_Points1D;
    Y1DCell = X1DCell + N_Points1D;
    xderiv_Neigh1D = Y1DCell + N_Points1D;
    yderiv_Neigh1D = xderiv_Neigh1D + N_Points1D;
    xyval_Neigh1D = yderiv_Neigh1D + N_Points1D;
    xderiv_Cell1D = xyval_Neigh1D + N_Points1D;
    yderiv_Cell1D = xderiv_Cell1D + N_Points1D;
    xyval_Cell1D = yderiv_Cell1D + N_Points1D;
    FEFunctValuesNeigh = new double[(3*N_Points1D+1)*MaxN_BaseFunctions2D];
    k = MaxN_BaseFunctions2D;
    for (i=0;i<N_Points1D;i++)
    {
	xyval_refNeigh1D[i]= FEFunctValuesNeigh+k;
	k += MaxN_BaseFunctions2D;
	xderiv_refNeigh1D[i]= FEFunctValuesNeigh+k;
	k += MaxN_BaseFunctions2D;
	yderiv_refNeigh1D[i]= FEFunctValuesNeigh+k;
	k += MaxN_BaseFunctions2D;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER>0)
      check_cont=1;
    else
      check_cont=0;

    /*************************************************************************/
    /*  strong residual                                                      */
    /*************************************************************************/
    hK = cell->GetDiameter();
    N_Edges = cell->GetN_Edges();
     for (i=0;i<N_Edges;i++)
    {
      TDatabase::ParamDB->INTERNAL_VERTEX_X[i] = cell->GetVertex(i)->GetX();
      TDatabase::ParamDB->INTERNAL_VERTEX_Y[i] = cell->GetVertex(i)->GetY();
    }
    if (N_Edges==3)
      TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;

    meas = cell->GetMeasure();
    strong_residual = 0;
    for(i=0;i<N_Points;i++)                       // for all quadrature points
    {
      coeff = coeffs[i];
      //cout << " eps " << coeff[0] << " b1 " << coeff[1] << " b2 " << coeff[2]
      //   << " c " << coeff[3] << " f " << coeff[4] << endl;

      deriv = Derivatives[i];                     // all derivatives in quadrature points
      w = weights[i]*AbsDetjk[i];

      //      cout << " xx " << deriv[3] << " yy " << deriv[4] << " x " << deriv[0]
      //   << " y " << deriv[1] << " c " << deriv[2] << endl;

      e1 = -coeff[0]*(deriv[3]+deriv[4])+coeff[1]*deriv[0]+coeff[2]*deriv[1]
        +coeff[3]*deriv[2]
        -coeff[4];                                // strong residual
      if(  TDatabase::ParamDB->P15 == -4711.0)    // time-dependent problem
        e1 += coeff[6];

      strong_residual += w*e1*e1;                 // L^2 norm
    }                                             // endfor i
    double hK_tilde = Mesh_size_in_convection_direction(hK, coeff[1], coeff[2]);

    alpha[0] = hK*hK;                             // weight for H^1 estimator
    alpha[1] = hK*hK*hK*hK;                       // weight for L^2 estimator
    alpha[2] = hK*hK/coeff[0];                    // weight for energy norm estimator
    if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
    {
	if (1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY<alpha[2])
	    alpha[2] = 1.0/TDatabase::ParamDB->INTERNAL_COERCIVITY; // update weight for energy norm estimator
    }
    alpha[3] = alpha[2];
    alpha[4] = alpha[3];
    double linfb = fabs(coeff[1]);
    if (fabs(coeff[2]) > linfb)
      linfb = fabs(coeff[2]);
    if (TDatabase::ParamDB->DISCTYPE == SUPG)
    {
      // compute stabilization parameter
      delta_K = Compute_SDFEM_delta(hK, coeff[0], coeff[1], coeff[2], coeff[3], linfb);
      if (alpha[4] > 24 * delta_K)
         alpha[4] =  24 * delta_K;
      // second contribution
      alpha[4] +=  24 * delta_K;
    }
    alpha[5] = hK;
    if(TDatabase::ParamDB->DISCTYPE == SUPG)
    {
      alpha[5] = hK*hK/(3*sqrt(10.0)*coeff[0]);
      if (coeff[3]>0)
      {
        if (1.0/coeff[3]<alpha[5])
          alpha[5] = 1.0/coeff[3];
      }
      // compute stabilization parameter
      if (hK_tilde/(sqrt(2.0)*linfb)<alpha[5])
        alpha[5] = hK/(sqrt(2.0)*linfb);
      alpha[5] *= 2*alpha[5];
    }
    for (i=1;i<N_estimators;i++)
      estimated_error[i] = alpha[i-1]*strong_residual;
    /*************************************************************************/
    /*  compute jumps across the edges                                       */
    /*************************************************************************/
    // no jumps in this estimator
    beta[3] = 0;
    N_Edges=cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)                        // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
        (joint->GetType() == IsoBoundEdge))       // boundary edge
      {
        boundedge=(TBoundEdge *)joint;
        BoundComp=boundedge->GetBoundComp();      // get boundary component
        boundedge->GetParameters(t0, t1);         // parameter interval
        comp=BoundComp->GetID();                  // boundary id
                                                  // type of boundary condition
        BoundaryConds[0](comp, (t0+t1)/2.0, Cond0);
                                                  // at midpoint of boundary
        switch(Cond0)
        {
          case DIRICHLET:                         // boundary is Dirichlet
            if(TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
            {
              for (i=0;i<N_estimators;i++)
                estimated_error[i] = 0;
              delete xietaval_refNeigh1D[0][0];
              delete xi1DNeigh;
              delete FEFunctValuesNeigh;
              return;
            }
            break;                                // no error
          case NEUMANN:                           // boundary is Neumann
            boundedge->GetXYofT(t0,x0,y0);        // coordinates at begin of parameter interval
            boundedge->GetXYofT(t1,x1,y1);        // coordinates at end of parameter interval
            nx = y1-y0;                           // outer normal vector
            ny = x0-x1;
            hE = sqrt(nx*nx+ny*ny);               // length of edge
            nx /= hE;                             // normalized normal vector
            ny /= hE;
            jump=0;
            for (i=0;i<N_Points1D;i++)            // compute difference to Neumann condition
            {
              x0 =  X1D[j][i];
              y0 =  Y1D[j][i];
              boundedge->GetXYofT(t0,x0,y0);      // coordinates at quadrature points
                                                  // Neumann data
              BoundaryValues[0](comp,t0,neumann_data);
              e1 = coeff[0] *(xderiv_1D[j][i]*nx+yderiv_1D[j][i]*ny) - neumann_data;
              w = weights1D[i]*hE/2.0;
              jump+= w*e1*e1;                     // integral on the edge
            }
            beta[0]=hE;                           // weight for H^1 estimator
            beta[1]=hE*hE*hE;                     // weight for L^2 estimator
            beta[2]=hE/coeff[0];                  // weight for energy norm estimator
	    if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
	    {
		w = 1/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
		if (w < beta[2])
		    beta[2]= w;
	    }
            if (24.0/linfb< beta[2])
            {
              beta[4] = 24.0/linfb;
            }
            else
            {
              beta[4] = beta[2];
            }
            /*
            if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
            {
              linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
              if (1.0/linfb < beta[4])
                beta[4] = 1.0/linfb;
            }
            linfb = hE / coeff[0];
            if (linfb < beta[4])
              beta[4] = linfb;
            */
            beta[5] = 1.0;
            beta[5] = alpha[5]*hE/(4.0*meas);
            for (i=1;i<N_estimators;i++)
              estimated_error[i] += beta[i-1]*jump;
            break;
        default:
         cerr << "Only few BC implementation done "  << endl;
          exit (-1);
         break;    
        }                                         // endswitch
      }                                           // end boundary edge
      else                                        // begin inner edge
      {
        refdesc=cell->GetRefDesc();               // get refinement descriptor
        refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
        ver0=cell->GetVertex(TmpEdVer[2*j]);      // get vertices of face j
        ver1=cell->GetVertex(TmpEdVer[2*j+1]);
                                                  // coordinates of face j
        cell_x0 = cell->GetVertex(TmpEdVer[2*j])->GetX();
        cell_y0 = cell->GetVertex(TmpEdVer[2*j])->GetY();
        cell_x1 = cell->GetVertex(TmpEdVer[2*j+1])->GetX();
        cell_y1 = cell->GetVertex(TmpEdVer[2*j+1])->GetY();
        jump=0;
        neigh=joint->GetNeighbour(cell);
        nx = cell_y1 - cell_y0;                   // compute normal
        ny = cell_x0 - cell_x1;
        hE = sqrt(nx*nx+ny*ny);                   // length of edge
        nx /= hE;                                 // normalized normal vector
        ny /= hE;
        if (ee_verbose>1)
        {
          cout << " A " << cell_x0 << " " << cell_y0;
          cout << " B " << cell_x1 << " " << cell_y1;
          cout << " n " << nx << " " << ny;
        }

        /*************************************************************************/
        /*  no neighbour, find neighbour of parent                               */
        /*************************************************************************/
        if(!neigh)
        {
          // there is no neighbour on the same level
          //  => finer cell in 1 regularity
          parent = cell->GetParent();             // parent cell
          refdesc= parent->GetRefDesc();
          refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerParent);
          refdesc->GetChildEdge(TmpCE,MaxLen1);
          refdesc->GetNewEdgeOldEdge(TmpoEnlE);
          l=0;
          while(parent->GetChild(l)!=cell) l++;   // local child number
          parent_edge = TmpCE[l*MaxLen1+j];       // number of father edge
          parent_edge = TmpoEnlE[parent_edge];    // number of father edge

          parent_joint= parent->GetJoint(parent_edge);
                                                  // neighbour to parent
          neigh=parent_joint->GetNeighbour(parent);
          if (!neigh)
          {
            cout << "Hier sollte man aber nicht hinkommen 2 !"<< endl;
          }

          neigh_edge=0;
          while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=parent) neigh_edge ++;
                                                  // vertices of edge
          ver2 = neigh->GetVertex(TmpEdVerParent[2*neigh_edge]);
          ver3 = neigh->GetVertex(TmpEdVerParent[2*neigh_edge+1]);
          if (ver1==ver2)                       // first part of long edge
          {
            part = -1;
          }
          else if (ver0==ver3)                 // second part of long edge
          {
            part = 1;
          }
          else
          {
            cout << "Hier sollte man aber nicht hinkommen 4 !"<< endl;
          }

          neigh_N_ = neigh->GetClipBoard();       // number of neighbour in iterator
          if(neigh_N_==-1)
          {
            cout << "Hier sollte man aber nicht hinkommen 3 !"<< endl;
          }
                                                  // finite element on neighbour
          CurrEleNeigh = fespace->GetFE2D(neigh_N_,neigh);
          eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);

                                                  // basis functions on neighbout
          BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
          N_Neigh = eleNeigh->GetN_DOF();         // number of basis functions

          bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
                                                  // referenz cell of neighbour
          bf2DrefelementsNeigh = bfNeigh->GetRefElement();

          if (conform_grid)
          {
            switch(bf2DrefelementsNeigh)          // compute coordinates of line quadrature
            {                                     // points in reference cell
              case BFUnitSquare :                 // edge 0
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = -zeta[i];
                    eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge==1)
                {                                 // edge 1
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = 1;
                    eta1DNeigh[i] = -zeta[i];
                  }
                }
                if (neigh_edge==2)
                {                                 // edge 2
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = zeta[i];
                    eta1DNeigh[i]= 1;
                  }
                }

                if (neigh_edge==3)
                {                                 // edge 3
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = -1;
                    eta1DNeigh[i]= zeta[i];
                  }
                }
                break;

              case BFUnitTriangle :
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (-zeta[i]+1)/2;
                    eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge==1)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (zeta[i]+1)/2;
                    eta1DNeigh[i] = (-zeta[i]+1)/2;
                  }
                }
                if (neigh_edge==2)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = 0;
                    eta1DNeigh[i] = (zeta[i] +1)/2;
                  }
                }
                break;
            }
          }
          else
          {
            switch(bf2DrefelementsNeigh)          // compute coordinates of line quadrature
            // this is only for 1-regular triangulations
            {                                     // points in reference cell
              case BFUnitSquare :                 // edge 0
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (-zeta[i]+part)/2;
                    eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge==1)
                {                                 // edge 1
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = 1;
                    eta1DNeigh[i] = (-zeta[i]+part)/2;
                  }
                }
                if (neigh_edge==2)
                {                                 // edge 2
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (zeta[i]-part)/2;
                    eta1DNeigh[i]= 1;
                  }
                }

                if (neigh_edge==3)
                {                                 // edge 3
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = -1;
                    eta1DNeigh[i]= (zeta[i]-part)/2;
                  }
                }
                break;

              case BFUnitTriangle :
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    if (part==-1)
                      part=0;
                    xi1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                    eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge==1)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    if (part==1)
                      part=0;
                    xi1DNeigh[i] = ((zeta[i]+1)/2-part)/2;
                    if (part==0)
                      part=1;
                    if (part==-1)
                      part=0;
                    eta1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                    if (part==0)
                      part=-1;
                    //cout << "part " << part << endl;
                  }
                }
                if (neigh_edge==2)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    if (part==1)
                      part=0;
                    xi1DNeigh[i] = 0;
                    eta1DNeigh[i] = ((zeta[i] +1)/2-part)/2;
                  }
                }
                break;
            }
          }
          if (ee_verbose>1)
            for (i=0;i<N_Points1D;i++)            // for all quadrature points
              cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

                                                  // compute gradients in reference cell of the neighbour
          for (i=0;i<N_Points1D;i++)              // for all quadrature points
          {
            bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
            bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
            bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
          }
                                                  // reftrafo of neighbour
          RefTransNeigh= eleNeigh->GetRefTransID();
          TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);

          DOF = GlobalNumbers + BeginIndex[neigh_N_];
          for(i=0;i<N_Neigh;i++)
          {
            FEFunctValuesNeigh[i] = Values[DOF[i]];
            if (ee_verbose>1)
              cout << " value " <<  FEFunctValuesNeigh[i] << endl;
          }
          for (i=0;i<N_Points1D;i++)              // get values and derivatives in original cell
          {
            TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
              eta1DNeigh[i],
              TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
              Coll, (TGridCell *)neigh,
              xietaval_refNeigh1D[BaseFunctNeigh][i],
              xideriv_refNeigh1D[BaseFunctNeigh][i],
              etaderiv_refNeigh1D[BaseFunctNeigh][i],
              xyval_refNeigh1D[i],
              xderiv_refNeigh1D[i],
              yderiv_refNeigh1D[i]);
          }

          for(i=0;i<N_Points1D;i++)               // for all quadrature points
          {
            val[0]=val[1]= val[2]=0;
            for(l=0;l<N_Neigh;l++)                // for all basis functions
            {
                                                  // accumulate value of derivative
              val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
              val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
              val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l];
              if (ee_verbose>1)
                cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " <<
                  yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
            }                                     // endfor l
            xderiv_Neigh1D[i]= val[0];            // for k-th
            yderiv_Neigh1D[i]= val[1];            // for k-th
            xyval_Neigh1D[i]= val[2];             // for k-th
          }                                       // endfor i

          TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh,
            eta1DNeigh,
            X1DNeigh, Y1DNeigh, absdet1D);
          jump=0.0;
          absdetjk1D = hE/2.0;
          for (i=0;i<N_Points1D;i++)              // compute jump
          {
            if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
              cout << " wrong quad points 1 " << X1D[j][i] << " , " << Y1D[j][i]
                << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
            if (check_cont)
              if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
            {
              cout << "quad points a " << X1D[j][i] << " , " << Y1D[j][i] << endl;
              cout << " i " << i << " vala " << xyval_1D[j][i]<< " neigha " << xyval_Neigh1D[i]<< " " << fabs(xyval_1D[j][i]-xyval_Neigh1D[i]) << endl;
            }
            e1 =coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
              +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny);
            if (ee_verbose>1)
              cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
            w = weights1D[i]*absdetjk1D;
            jump+= w*e1*e1;                       // integral on the edge
          }
          if (ee_verbose>1)
            cout << "jump " << jump << endl;

          beta[0]=hE;                             // weight for H^1 estimator
          beta[1]=hE*hE*hE;                       // weight for L^2 estimator
          beta[2]=hE/coeff[0];                    // weight for energy norm estimator
	  if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
	  {
	      w = 1/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
	      if (w < beta[2])
		  beta[2]= w;
	  }
	  //beta[2] *= 2.0;
          if (24.0/linfb < beta[2])
          {
            beta[4] = 24.0/linfb;
          }
          else
          {
            beta[4] = beta[2];
          }
          /*
          beta[4] = 24;
          if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
          {
            linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
            if (1.0/linfb < beta[4])
              beta[4] = 1.0/linfb;
          }
          linfb = hE / coeff[0];
          if (linfb < beta[4])
            beta[4] = linfb;*/
          beta[5] = 1.0;
          beta[5] = alpha[5]*hE/(4.0*meas);
          for (i=1;i<N_estimators;i++)
            estimated_error[i] += beta[i-1]*jump/2.0;
        }                                         // end no neighbour
        /*************************************************************************/
        /*  neighbour is not on the finest level, find children of neighbour     */
        /*************************************************************************/
        else                                      // there is a neighbour on the same level
        {
          n=neigh->GetClipBoard();
          if(n==-1)
          {
            // the neighbour is no member of the collection
            // check whether the children of neigh are in collection
            // find the local edge of neigh on which cell is -> l

            edge2neigh=0;
            while(neigh->GetJoint(edge2neigh)->GetNeighbour(neigh)!=cell)
              edge2neigh++;                       // find connections between cells
            refdesc=neigh->GetRefDesc();          // ref desc of neighbour
                                                  // get edges
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);
                                                  // get connection to child edges
            refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);
                                                  // get cell belonging to child edge (TmpEC)
            refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);
                                                  // get local no.s of child edge
            refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);
            if (conform_grid)
              N_child = 1;
            else
              N_child = 2;                        // not general  !!!

            for(k=0;k < N_child ; k++)            // find children of neigh on face l -> child
            {
              edge1=TmpoEnE[edge2neigh*MaxLen1+k];// edge child, not general !!!
              chnum1=TmpEC[edge1*MaxLen2];        // local number of child cell
              child =neigh->GetChild(chnum1);     // child cell
              child_N_=child->GetClipBoard();     // id of child cell
                                                  // get local indices of child edge
              refdesc->GetEdgeChildIndex(TmpECI,TmpLen3, MaxLen3);
              l_child = TmpECI[edge1*MaxLen3];    // local index of child edge

              refdesc_child=child->GetRefDesc();  // ref desc of child
                                                  // conn. edge -> vertices
              refdesc_child->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
                                                  // vertices of edge
              ver2 = child->GetVertex(TmpEdVer[2*l_child]);
              ver3 = child->GetVertex(TmpEdVer[2*l_child+1]);

              if (ee_verbose>1)
              {
                cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
              }

              if (ver1==ver2)
              {
                part=1;
              }
              else if (ver0==ver3)
              {
                part=-1;
              }
              else
              {
                cout << " something wrong 5 " <<  endl;
                cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
              }
              // now from point of view of child cell -> cell becomes the neighbour
              // prepare intergration for the half part of edge j

              neigh_N_ = cell->GetClipBoard();    // number of original cell  in iterator
              if(neigh_N_==-1)
              {
                cout << "Hier sollte man aber nicht hinkommen 33 !"<< endl;
              }
                                                  // finite element on neighbour
              CurrEleNeigh = fespace->GetFE2D(neigh_N_,cell);
              eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);

                                                  // basis functions on neighbout
              BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
              N_Neigh = eleNeigh->GetN_DOF();     // number of basis functions

              bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
                                                  // referenz cell of neighbour
              bf2DrefelementsNeigh = bfNeigh->GetRefElement();

              neigh_edge = j;
              if (conform_grid)
              {
                switch(bf2DrefelementsNeigh)      // compute coordinates of line quadrature
                {                                 // points in reference cell
                  case BFUnitSquare :             // edge 0
                    if (neigh_edge==0)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = -zeta[i];
                        eta1DNeigh[i] = -1;
                      }
                    }
                    if (neigh_edge==1)
                    {                             // edge 1
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = 1;
                        eta1DNeigh[i] = -zeta[i];
                      }
                    }
                    if (neigh_edge==2)
                    {                             // edge 2
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = zeta[i];
                        eta1DNeigh[i]= 1;
                      }
                    }

                    if (neigh_edge==3)
                    {                             // edge 3
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = -1;
                        eta1DNeigh[i]= zeta[i];
                      }
                    }
                    break;

                  case BFUnitTriangle :
                    if (neigh_edge==0)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = (-zeta[i]+1)/2;
                        eta1DNeigh[i] = 0;
                      }
                    }
                    if (neigh_edge==1)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = (zeta[i]+1)/2;
                        eta1DNeigh[i] = (-zeta[i]+1)/2;
                      }
                    }
                    if (neigh_edge==2)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = 0;
                        eta1DNeigh[i] = (zeta[i] +1)/2;
                      }
                    }
                    break;
                }
              }
              else
              {
                switch(bf2DrefelementsNeigh)      // compute coordinates of line quadrature
                // this is only for 1-regular triangulations
                {                                 // points in reference cell
                  case BFUnitSquare :             // edge 0
                    if (neigh_edge==0)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = (-zeta[i]+part)/2;
                        eta1DNeigh[i] = -1;
                      }
                    }
                    if (neigh_edge==1)
                    {                             // edge 1
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = 1;
                        eta1DNeigh[i] = (-zeta[i]+part)/2;
                      }
                    }
                    if (neigh_edge==2)
                    {                             // edge 2
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = (zeta[i]-part)/2;
                        eta1DNeigh[i]= 1;
                      }
                    }

                    if (neigh_edge==3)
                    {                             // edge 3
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        xi1DNeigh[i] = -1;
                        eta1DNeigh[i]= (zeta[i]-part)/2;
                      }
                    }
                    break;

                  case BFUnitTriangle :
                    if (neigh_edge==0)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        if (part==-1)
                          part=0;
                        xi1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                        eta1DNeigh[i] = 0;
                      }
                    }
                    if (neigh_edge==1)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        if (part==1)
                          part=0;
                        xi1DNeigh[i] = ((zeta[i]+1)/2-part)/2;
                        if (part==0)
                          part=1;
                        if (part==-1)
                          part=0;
                        eta1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                        if (part==0)
                          part=-1;
                      }
                    }
                    if (neigh_edge==2)
                    {
                      for (i=0;i<N_Points1D;i++)  // for all quadrature points
                      {
                        if (part==1)
                          part=0;
                        xi1DNeigh[i] = 0;
                        eta1DNeigh[i] = ((zeta[i] +1)/2-part)/2;
                      }
                    }
                    break;
                }
              }
              if (ee_verbose>1)
                for (i=0;i<N_Points1D;i++)        // for all quadrature points
                  cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

              // compute gradients in reference cell of the neighbour
              for (i=0;i<N_Points1D;i++)          // for all quadrature points
              {
                bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
                                                  // reftrafo of neighbour
              RefTransNeigh= eleNeigh->GetRefTransID();
              TFEDatabase2D::SetCellForRefTrans(cell,RefTransNeigh);

              DOF = GlobalNumbers + BeginIndex[neigh_N_];
              for(i=0;i<N_Neigh;i++)
              {
                FEFunctValuesNeigh[i] = Values[DOF[i]];
                if (ee_verbose>1)
                  cout << " value " <<  FEFunctValuesNeigh[i] << endl;
              }

              for (i=0;i<N_Points1D;i++)          // get values and derivatives in original cell
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                  eta1DNeigh[i],
                  TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                  Coll, (TGridCell *)neigh,
                  xietaval_refNeigh1D[BaseFunctNeigh][i],
                  xideriv_refNeigh1D[BaseFunctNeigh][i],
                  etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  xyval_refNeigh1D[i],
                  xderiv_refNeigh1D[i],
                  yderiv_refNeigh1D[i]);
              }

              for(i=0;i<N_Points1D;i++)           // for all quadrature points
              {
                val[0]=val[1]= val[2]=0;
                for(l=0;l<N_Neigh;l++)            // for all basis functions
                {
                                                  // accumulate value of derivative
                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l];
                  if (ee_verbose>1)
                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " <<
                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                }                                 // endfor l
                xderiv_Cell1D[i]= val[0];         // for k-th
                yderiv_Cell1D[i]= val[1];         // for k-th
                xyval_Cell1D[i]= val[2];          // for k-th
              }                                   // endfor i

              TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh,
                eta1DNeigh,
                X1DCell, Y1DCell, absdet1D);
              // prepare integration for the child of the neighbour belong to the half part
              // of edge j

              neigh_N_ = child->GetClipBoard();   // number of neighbour in iterator
                                                  // finite element on neighbour
              CurrEleNeigh = fespace->GetFE2D(neigh_N_,child);
              eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);

                                                  // basis functions on neighbout
              BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
              N_Neigh = eleNeigh->GetN_DOF();     // number of basis functions

              bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
                                                  // referenz cell of neighbour
              bf2DrefelementsNeigh = bfNeigh->GetRefElement();

              neigh_edge = l_child;
              switch(bf2DrefelementsNeigh)        // compute coordinates of line quadrature
              {                                   // points in reference cell
                case BFUnitSquare :               // edge 0
                  if (neigh_edge==0)
                  {
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = zeta[i];
                      eta1DNeigh[i] = -1;
                    }
                  }
                  if (neigh_edge==1)
                  {                               // edge 1
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = 1;
                      eta1DNeigh[i] = zeta[i];
                    }
                  }
                  if (neigh_edge==2)
                  {                               // edge 2
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = -zeta[i];
                      eta1DNeigh[i]= 1;
                    }
                  }

                  if (neigh_edge==3)
                  {                               // edge 3
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = -1;
                      eta1DNeigh[i]= -zeta[i];
                    }
                  }
                  break;

                case BFUnitTriangle :
                  if (neigh_edge==0)
                  {
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = (zeta[i]+1)/2;
                      eta1DNeigh[i] = 0;
                    }
                  }
                  if (neigh_edge==1)
                  {
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = (-zeta[i]+1)/2;
                      eta1DNeigh[i] = (zeta[i]+1)/2;
                    }
                  }
                  if (neigh_edge==2)
                  {
                    for (i=0;i<N_Points1D;i++)    // for all quadrature points
                    {
                      xi1DNeigh[i] = 0;
                      eta1DNeigh[i] = (-zeta[i] +1)/2;
                    }
                  }
                  break;
              }

              if (ee_verbose>1)
                for (i=0;i<N_Points1D;i++)        // for all quadrature points
                  cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

              // compute gradients in reference cell of the neighbour
              for (i=0;i<N_Points1D;i++)          // for all quadrature points
              {
                bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
              }
                                                  // reftrafo of neighbour
              RefTransNeigh= eleNeigh->GetRefTransID();
              TFEDatabase2D::SetCellForRefTrans(child,RefTransNeigh);

              DOF = GlobalNumbers + BeginIndex[neigh_N_];
              for(i=0;i<N_Neigh;i++)
              {
                FEFunctValuesNeigh[i] = Values[DOF[i]];
                if (ee_verbose>1)
                  cout << " value " <<  FEFunctValuesNeigh[i] << endl;
              }

              for (i=0;i<N_Points1D;i++)          // get values and derivatives in original cell
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                  eta1DNeigh[i],
                  TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                  Coll, (TGridCell *)neigh,
                  xietaval_refNeigh1D[BaseFunctNeigh][i],
                  xideriv_refNeigh1D[BaseFunctNeigh][i],
                  etaderiv_refNeigh1D[BaseFunctNeigh][i],
                  xyval_refNeigh1D[i],
                  xderiv_refNeigh1D[i],
                  yderiv_refNeigh1D[i]);
              }

              for(i=0;i<N_Points1D;i++)           // for all quadrature points
              {
                val[0]=val[1]= val[2]=0;
                for(l=0;l<N_Neigh;l++)            // for all basis functions
                {
                                                  // accumulate value of derivative
                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l];
                  if (ee_verbose>1)
                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " <<
                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                }                                 // endfor l
                xderiv_Neigh1D[i]= val[0];        // for k-th
                yderiv_Neigh1D[i]= val[1];        // for k-th
                xyval_Neigh1D[i]= val[2];         // for k-th
              }                                   // endfor i

              TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh,
                eta1DNeigh,
                X1DNeigh, Y1DNeigh, absdet1D);
              jump=0.0;
              absdetjk1D = hE/(2.0*N_child);      // only half edge is considered
              for (i=0;i<N_Points1D;i++)          // compute jump
              {
                if ((fabs(X1DCell[i]-X1DNeigh[i])+fabs(Y1DCell[i]-Y1DNeigh[i]))>1e-8)
                  cout << " wrong quad points 2 " << X1DCell[i] << " , " << Y1DCell[i]
                    << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                if (check_cont)
                  if (fabs(xyval_Neigh1D[i]-xyval_Cell1D[i])>1e-8)
                {
                  cout << "quad points b " << X1DCell[i] << " , " << Y1DCell[i] << endl;
                  cout << " i " << i << " valb " << xyval_Cell1D[i]<< " neighb " << xyval_Neigh1D[i] << " " << fabs(xyval_Cell1D[i]-xyval_Neigh1D[i]) << endl;
                }
                e1 =coeff[0] *((xderiv_Cell1D[i]-xderiv_Neigh1D[i])*nx
                  +(yderiv_Cell1D[i]-yderiv_Neigh1D[i])*ny);
                if (ee_verbose>1)
                  cout << i<< " jumpx " << xderiv_Cell1D[i] << " " << xderiv_Neigh1D[i] << endl;
                w = weights1D[i]*absdetjk1D;
                jump+= w*e1*e1;                   // integral on the edge
              }
              if (ee_verbose>1)
                cout << "jump " << jump << endl;
              hE2=hE/N_child;
              beta[0]=hE2;                        // weight for H^1 estimator
              beta[1]=hE2*hE2*hE2;                // weight for L^2 estimator
              beta[2]=hE2/coeff[0];               // weight for energy norm estimator
	      if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
	      {
		  w = 1/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
		  if (w < beta[2])
		      beta[2]= w;
	      }
              //beta[2] *= 2.0;
              if (24.0/linfb< beta[2])
              {
                beta[4] = 24.0/linfb;
              }
              else
              {
                beta[4] = beta[2];
              }
              /*
              beta[4] = 24;
              if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
              {
                linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
                if (1.0/linfb < beta[4])
                  beta[4] = 1.0/linfb;
              }
              linfb = hE / coeff[0];
              if (linfb < beta[4])
                beta[4] = linfb;
              */
              beta[5] = 1.0;
              beta[5] = alpha[5]*hE/(4.0*meas);
              for (i=1;i<N_estimators;i++)
                estimated_error[i] += beta[i-1]*jump/2.0;
            }
          }                                       // end clipboard==-1
          else
            /*************************************************************************/
            /*  neighbour is on the finest level                                     */
            /*************************************************************************/

          {                                       // the neighbour is a member of the collection
            // find the finite element on the other side
            //                      cout << " neighbour found " << endl;
            // find the local edge of neigh on which cell is -> l
            neigh_edge=0;
            while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;
            refdesc=neigh->GetRefDesc();
            refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);
            ver0=  cell->GetVertex(TmpEdVer[2*j]);
            ver1=  cell->GetVertex(TmpEdVer[2*j+1]);
                                                  // vertices of edge
            ver2 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge]);
            ver3 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge+1]);
            if (((ver0==ver2)&&(ver1==ver3))||((ver0==ver3)&&(ver1==ver2)))
              ;
            else
            {
              cout << "wrong edge " << endl;
            }

            // compute gradient at the quadrature points on the edge of
            // the neighbour element

            neigh_N_ = neigh->GetClipBoard();     // number of neighbour in iterator
                                                  // finite element on neighbour
            CurrEleNeigh = fespace->GetFE2D(neigh_N_,neigh);
            eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);

                                                  // basis functions on neighbout
            BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();
            N_Neigh = eleNeigh->GetN_DOF();       // number of basis functions

            bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);
                                                  // referenz cell of neighbour
            bf2DrefelementsNeigh = bfNeigh->GetRefElement();

            switch(bf2DrefelementsNeigh)          // compute coordinates of line quadrature
            {                                     // points in reference cell
              case BFUnitSquare :                 // edge 0
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = -zeta[i];
                    eta1DNeigh[i] = -1;
                  }
                }
                if (neigh_edge==1)
                {                                 // edge 1
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = 1;
                    eta1DNeigh[i] = -zeta[i];
                  }
                }
                if (neigh_edge==2)
                {                                 // edge 2
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = zeta[i];
                    eta1DNeigh[i]= 1;
                  }
                }

                if (neigh_edge==3)
                {                                 // edge 3
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = -1;
                    eta1DNeigh[i]= zeta[i];
                  }
                }
                break;

              case BFUnitTriangle :
                if (neigh_edge==0)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (-zeta[i]+1)/2;
                    eta1DNeigh[i] = 0;
                  }
                }
                if (neigh_edge==1)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = (zeta[i]+1)/2;
                    eta1DNeigh[i] = (-zeta[i]+1)/2;
                  }
                }
                if (neigh_edge==2)
                {
                  for (i=0;i<N_Points1D;i++)      // for all quadrature points
                  {
                    xi1DNeigh[i] = 0;
                    eta1DNeigh[i] = (zeta[i] +1)/2;
                  }
                }
                break;
            }

            if (ee_verbose>1)
              for (i=0;i<N_Points1D;i++)          // for all quadrature points
                cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

            // compute gradients in reference cell of the neighbour
            for (i=0;i<N_Points1D;i++)            // for all quadrature points
            {
              bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
              bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
            }
                                                  // reftrafo of neighbour
            RefTransNeigh= eleNeigh->GetRefTransID();
            TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);

            DOF = GlobalNumbers + BeginIndex[neigh_N_];
            for(i=0;i<N_Neigh;i++)
            {
              FEFunctValuesNeigh[i] = Values[DOF[i]];
              if (ee_verbose>1)
                cout << " value " <<  FEFunctValuesNeigh[i] << endl;
            }
            for (i=0;i<N_Points1D;i++)            // get values and derivatives in original cell
            {
              TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                eta1DNeigh[i],
                TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                Coll, (TGridCell *)neigh,
                xietaval_refNeigh1D[BaseFunctNeigh][i],
                xideriv_refNeigh1D[BaseFunctNeigh][i],
                etaderiv_refNeigh1D[BaseFunctNeigh][i],
                xyval_refNeigh1D[i],
                xderiv_refNeigh1D[i],
                yderiv_refNeigh1D[i]);
            }

            for(i=0;i<N_Points1D;i++)             // for all quadrature points
            {
              val[0]=val[1]= val[2]=0;
              for(l=0;l<N_Neigh;l++)              // for all basis functions
              {
                                                  // accumulate value of derivative
                val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l];
                                                  // accumulate value of derivative
                val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l];
                if (ee_verbose>1)
                  cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " <<
                    yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
              }                                   // endfor l
              xderiv_Neigh1D[i]= val[0];          // for k-th
              yderiv_Neigh1D[i]= val[1];          // for k-th
              xyval_Neigh1D[i]= val[2];           // for k-th
            }                                     // endfor i

            TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh,
              eta1DNeigh,
              X1DNeigh, Y1DNeigh, absdet1D);

            jump=0.0;
            absdetjk1D = hE/2.0;
            for (i=0;i<N_Points1D;i++)            // compute jump
            {
              if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
                cout << " wrong quad points 0 " << X1D[j][i] << " , " << Y1D[j][i]
                  << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
              if (check_cont)
                if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
                  cout << " i " << i << " valc " << xyval_1D[j][i]<< " neighc " << xyval_Neigh1D[i]<< endl;
              e1 =coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
                +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny);
              if (ee_verbose>1)
                cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
              w = weights1D[i]*absdetjk1D;
              jump+= w*e1*e1;                     // integral on the edge
            }
            if (ee_verbose>1)
              cout << "jump " << jump << endl;
            beta[0]=hE;                           // weight for H^1 estimator
            beta[1]=hE*hE*hE;                     // weight for L^2 estimator
            beta[2]=hE/coeff[0];                  // weight for energy norm estimator
	    if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
	    {
		w = 1/sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY * coeff[0]);
		if (w < beta[2])
		    beta[2]= w;
	    }
	    //beta[2] *= 2.0;
            if (24.0/linfb< beta[2])
            {
              beta[4] = 24.0/linfb;
            }
            else
            {
              beta[4] = beta[2];
            }
            /*beta[4] = 24;
            if (TDatabase::ParamDB->INTERNAL_COERCIVITY>0)
            {
              linfb = sqrt(TDatabase::ParamDB->INTERNAL_COERCIVITY) * sqrt(coeff[0]);
              if (1.0/linfb < beta[4])
                beta[4] = 1.0/linfb;
            }
            linfb = hE / coeff[0];
            if (linfb < beta[4])
              beta[4] = linfb;
            */
            beta[5] = 1.0;
            beta[5] = alpha[5]*hE/(4.0*meas);
            for (i=1;i<N_estimators;i++)
              estimated_error[i] += beta[i-1]*jump/2.0;
          }                                       // end neighbour is member of the collection
        }                                         // end neighbour on the finer level
      }                                           // end inner edge
    }                                             // end for j
    delete xi1DNeigh;
    delete FEFunctValuesNeigh;
  }                                               // end residual based estimators

  delete xietaval_refNeigh1D[0][0];

  //for (i=0;i<4;i++)
  //  cout << i << " estimated error " << estimated_error[i] << endl;
  if (delta_K!=4711)
  {
      for (i=0;i<N_estimators;i++)
    estimated_local_error[i]=estimated_error[i];
  }
  else
  {
    for (i=0;i<5;i++)
      estimated_local_error[i]=estimated_error[i];
    // no SUPG, i.e. no delta_K
    estimated_local_error[5] = estimated_local_error[6] = 0;
  }
}

void  TCD2DErrorEstimator::EstimateCellError(TFESpace2D *fespace,
                                           TBaseCell *cell,
                                           int N_Points, 
                                           double *X, 
                                           double *Y, 
                                           double *AbsDetjk, 
                                           double *weights,
                                           double **Derivatives, 
                                           double **coeffs, 
                                           BoundCondFunct2D **BoundaryConds,
                                           BoundValueFunct2D **BoundaryValues,
                                           int N_Points1D, 
                                           double *zeta,
                                           double *X1D[4], 
                                           double *Y1D[4], 
                                           double *weights1D,
                                           double *xyval_1D[4],
                                           double *xderiv_1D[4],
                                           double *yderiv_1D[4],
                                           int *GlobalNumbers, 
                                           int *BeginIndex,
                                           int *DOF,
                                           double *Values,                                           
                                           double *estimated_local_error)
{
  int i,j,k,l,n,N_Edges,comp,parent_edge,MaxLen1,MaxLen2,MaxLen3,N_child,neigh_edge;
  int chnum1,l_child,child_N_,edge1,neigh_N_,N_Neigh;
  double *deriv, w,e1,e2,*coeff,strong_residual,alpha[N_estimators-1],beta[N_estimators-1],hK,hE,val[3],hE2;
  double estimated_error[4],t0,t1,nx,ny,x0,x1,y0,y1,neumann_data,jump,absdetjk1D;
  double absdet1D[MaxN_QuadPoints_2D];
  double cell_x0,cell_y0,cell_x1,cell_y1,*xi1DNeigh,*eta1DNeigh, *FEFunctValuesNeigh;
  double *X1DNeigh,*Y1DNeigh,*X1DCell,*Y1DCell;
  const int *TmpoEnE, *TmpLen1, *TmpEC, *TmpLen2, *TmpLen3; 
  const int  *TmpoEnlE, *TmpEdVer, *TmpECI, *TmpCE, *TmpEdVerParent, *TmpEdVerNeigh;
  TJoint *joint,*parent_joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0;
  TRefDesc *refdesc,*refdesc_child;
  TBaseCell *neigh, *child, *parent;
  TVertex *ver0,*ver1, *ver2,*ver3;
  FE2D CurrEleNeigh;
  BaseFunct2D BaseFunctNeigh;
  QuadFormula2D QuadFormulaNeigh;
  TQuadFormula2D *qfNeigh;
  QuadFormula1D LineQuadFormulaNeigh;
  TQuadFormula1D *qf1DNeigh;
  TBaseFunct2D *bfNeigh;
  int N_Points1DNeigh,N_PointsNeigh;
  double *weights1DNeigh,*zetaNeigh,*weightsNeigh,*xiNeigh,*etaNeigh;
  TFE2D *eleNeigh;
  RefTrans2D RefTransNeigh;
  BF2DRefElements bf2DrefelementsNeigh;
  double xietaval_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double xideriv_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double etaderiv_refNeigh1D[MaxN_BaseFunctions2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;
  double *xderiv_Cell1D, *yderiv_Cell1D, *xyval_Cell1D;
  int part,edge2neigh;
  int ee_verbose = 1, check_cont, conform_grid=TDatabase::ParamDB->GRID_TYPE;
  TCollection *Coll;

  Coll = fespace->GetCollection();
  
  if (ee_verbose>1)
    {
      for (i=0;i<4;i++)
        {
          cout << "1DDD " ;
          for (k=0;k<N_Points1D;k++)
            cout << " x " << X1D[i][k] << " y " << Y1D[i][k];
          cout << endl;
        }
 
      for (j=0;j<4;j++)
        for (k=0;k<N_Points1D;k++)
          cout << "xx " << xyval_1D[j][k] <<  "dx " << xderiv_1D[j][k] << "dy " <<
            yderiv_1D[j][k] << endl;
    }

  for (i=0;i<4;i++)             // initalize local error estimates
    estimated_error[i]=0.0;

/*************************************************************************/
/*                                                                       */
/*   GRADIENT INDICATOR                                                  */
/*                                                                       */
/*************************************************************************/
  if (FELocalEstimator==cd_gradient_indicator)      // gradient indicator
    {
      for(i=0;i<N_Points;i++)                // for all quadrature points
        {
          deriv = Derivatives[i];            // all derivatives in quadrature points 
          w = weights[i]*AbsDetjk[i];

          e1 = deriv[0];                     // x derivative
          e2 = deriv[1];                     // y derivative
          estimated_error[0] += w*(e1*e1+e2*e2);
        } // endfor i
    }
  
/*************************************************************************/
/*                                                                       */
/*   RESIDUAL BASED EXPLICIT ERROR ESTIMATORS                            */
/*                                                                       */
/*************************************************************************/
  if ((FELocalEstimator>0)||(ErrorControl))
  {
    xi1DNeigh = new double[N_Points1D];
    eta1DNeigh = new double[N_Points1D];
    X1DNeigh = new double[N_Points1D];
    Y1DNeigh = new double[N_Points1D];
    X1DCell = new double[N_Points1D];
    Y1DCell = new double[N_Points1D];
    FEFunctValuesNeigh = new double[MaxN_BaseFunctions2D];
    for (i=0;i<N_Points1D;i++)
    {
      xyval_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
      xderiv_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
      yderiv_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
    }
    xderiv_Neigh1D = new double[N_Points1D];
    yderiv_Neigh1D = new double[N_Points1D];
    xyval_Neigh1D = new double[N_Points1D];
    xderiv_Cell1D = new double[N_Points1D];
    yderiv_Cell1D = new double[N_Points1D];
    xyval_Cell1D = new double[N_Points1D];
    if (TDatabase::ParamDB->ANSATZ_ORDER>0)
      check_cont=1;
    else
      check_cont=0;

/*************************************************************************/
/*  strong residual                                                      */
/*************************************************************************/
    hK = cell->GetDiameter();
    strong_residual = 0;
    for(i=0;i<N_Points;i++)                // for all quadrature points
    {
      coeff = coeffs[i];
      //cout << " eps " << coeff[0] << " b1 " << coeff[1] << " b2 " << coeff[2]
      //   << " c " << coeff[3] << " f " << coeff[4] << endl;
      
      deriv = Derivatives[i];            // all derivatives in quadrature points 
      w = weights[i]*AbsDetjk[i];
      
      //      cout << " xx " << deriv[3] << " yy " << deriv[4] << " x " << deriv[0]
      //   << " y " << deriv[1] << " c " << deriv[2] << endl;

      e1 = -coeff[0]*(deriv[3]+deriv[4])+coeff[1]*deriv[0]+coeff[2]*deriv[1] 
        +coeff[3]*deriv[2]
        -coeff[4];              // strong residual
      strong_residual += w*e1*e1;        // L^2 norm
    }                                     // endfor i
    alpha[0] = hK*hK;                      // weight for H^1 estimator
    alpha[1] = hK*hK*hK*hK;                // weight for L^2 estimator
    alpha[2] = 1;                          // weight for energy norm estimator
    if (hK*hK/coeff[0]<1)
      alpha[2] =hK*hK/coeff[0];           // update weight for energy norm estimator  
    alpha[3] = alpha[2];
    for (i=1;i<5;i++)
      estimated_error[i] = alpha[i-1]*strong_residual;
/*************************************************************************/
/*  compute jumps across the edges                                       */
/*************************************************************************/
       N_Edges=cell->GetN_Edges();
       for(j=0;j<N_Edges;j++)              // loop over all edges of cell
        {
          joint=cell->GetJoint(j);
          if ((joint->GetType() == BoundaryEdge)||
                (joint->GetType() == IsoBoundEdge)) // boundary edge 
            {
              boundedge=(TBoundEdge *)joint;  
              BoundComp=boundedge->GetBoundComp();  // get boundary component
              boundedge->GetParameters(t0, t1);     // parameter interval
              comp=BoundComp->GetID();              // boundary id 
              BoundaryConds[0](comp, (t0+t1)/2.0, Cond0); // type of boundary condition
                                                    // at midpoint of boundary 
              switch(Cond0)
                {
                case DIRICHLET:                         // boundary is Dirichlet
                  break;                                // no error
                case NEUMANN:                           // boundary is Neumann 
                  boundedge->GetXYofT(t0,x0,y0);       // coordinates at begin of parameter interval 
                  boundedge->GetXYofT(t1,x1,y1);        // coordinates at end of parameter interval
                  nx = y1-y0;                           // outer normal vector
                  ny = x0-x1;           
                  hE = sqrt(nx*nx+ny*ny);               // length of edge
                  nx /= hE;                             // normalized normal vector
                  ny /= hE;
                  jump=0;
                  for (i=0;i<N_Points1D;i++)           // compute difference to Neumann condition
                    {
                      x0 =  X1D[j][i];
                      y0 =  Y1D[j][i];                      
                      boundedge->GetXYofT(t0,x0,y0);   // coordinates at quadrature points 
                      BoundaryValues[0](comp,t0,neumann_data);  // Neumann data
                      e1 = coeff[0] *(xderiv_1D[j][i]*nx+yderiv_1D[j][i]*ny) - neumann_data;
                      w = weights1D[i]*hE/2.0;
                      jump+= w*e1*e1;                       // integral on the edge
                    }
                  beta[0]=hE;                          // weight for H^1 estimator
                  beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                  beta[2]=hE/coeff[0];                 // weight for energy norm estimator
                  if (1.0/sqrt(coeff[0])< beta[2])
                    beta[2]=1.0/sqrt(coeff[0]);
                  for (i=1;i<4;i++)
                    estimated_error[i] += beta[i-1]*jump;                  
                  break;
        default:
         cerr << "Only few BC implementation done "  << endl;
          exit (-1);
         break;    
                }                                   // endswitch 
             }                                       // end boundary edge
          else                                      // begin inner edge 
            {
              refdesc=cell->GetRefDesc();           // get refinement descriptor
              refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer); 
              ver0=cell->GetVertex(TmpEdVer[2*j]);  // get vertices of face j
              ver1=cell->GetVertex(TmpEdVer[2*j+1]);
              cell_x0 = cell->GetVertex(TmpEdVer[2*j])->GetX();  // coordinates of face j
              cell_y0 = cell->GetVertex(TmpEdVer[2*j])->GetY();
              cell_x1 = cell->GetVertex(TmpEdVer[2*j+1])->GetX();
              cell_y1 = cell->GetVertex(TmpEdVer[2*j+1])->GetY();
              jump=0;
              neigh=joint->GetNeighbour(cell);
              nx = cell_y1 - cell_y0;              // compute normal
              ny = cell_x0 - cell_x1;
              hE = sqrt(nx*nx+ny*ny);               // length of edge
              nx /= hE;                             // normalized normal vector
              ny /= hE;
              if (ee_verbose>1)
                {
                  cout << " A " << cell_x0 << " " << cell_y0;
                  cout << " B " << cell_x1 << " " << cell_y1;
                  cout << " n " << nx << " " << ny;
                }

/*************************************************************************/
/*  no neighbour, find neighbour of parent                               */
/*************************************************************************/
              if(!neigh)
                {
                  // there is no neighbour on the same level
                  //  => finer cell in 1 regularity
                  parent = cell->GetParent();                      // parent cell
                  refdesc= parent->GetRefDesc();
                  refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerParent); 
                  refdesc->GetChildEdge(TmpCE,MaxLen1);
                  refdesc->GetNewEdgeOldEdge(TmpoEnlE);
                  l=0;
                  while(parent->GetChild(l)!=cell) l++;           // local child number
                  parent_edge = TmpCE[l*MaxLen1+j];               // number of father edge
                  parent_edge = TmpoEnlE[parent_edge];            // number of father edge
                  
                  parent_joint= parent->GetJoint(parent_edge);
                  neigh=parent_joint->GetNeighbour(parent);        // neighbour to parent
                  if (!neigh)
                    {
                      cout << "Hier sollte man aber nicht hinkommen 2 !"<< endl;
                    }                                             
 
                  neigh_edge=0;
                  while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=parent) neigh_edge ++;
                  ver2 = neigh->GetVertex(TmpEdVerParent[2*neigh_edge]);          // vertices of edge
                  ver3 = neigh->GetVertex(TmpEdVerParent[2*neigh_edge+1]);
                  if (ver1==ver2)                              // first part of long edge
                    {
                      part = -1;
                    }
                  else if (ver0==ver3)                          // second part of long edge
                    {
                      part = 1;
                    }
                  else 
                    {
                      cout << "Hier sollte man aber nicht hinkommen 4 !"<< endl;
                    }                                             

                  neigh_N_ = neigh->GetClipBoard();                  // number of neighbour in iterator
                  if(neigh_N_==-1)
                    {
                      cout << "Hier sollte man aber nicht hinkommen 3 !"<< endl;
                    }    
                  CurrEleNeigh = fespace->GetFE2D(neigh_N_,neigh);   // finite element on neighbour
                  eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                  
                  BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                  N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions
                  
                  bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                  bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                  if (conform_grid)
                  {
                      switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                        {                                    // points in reference cell
                        case BFUnitSquare :                  // edge 0
                          if (neigh_edge==0)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = -zeta[i];
                                  eta1DNeigh[i] = -1;
                                }
                            }
                          if (neigh_edge==1)
                            {                               // edge 1
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = 1;
                                  eta1DNeigh[i] = -zeta[i];
                                }
                            }
                          if (neigh_edge==2)
                            {                               // edge 2
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = zeta[i];
                                  eta1DNeigh[i]= 1;
                                }
                            }

                          if (neigh_edge==3)
                            {                               // edge 3
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = -1;
                                  eta1DNeigh[i]= zeta[i];
                                }
                            }
                          break;

                        case BFUnitTriangle :
                           if (neigh_edge==0)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = (-zeta[i]+1)/2;
                                  eta1DNeigh[i] = 0;
                                }
                            }
                          if (neigh_edge==1)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = (zeta[i]+1)/2;
                                  eta1DNeigh[i] = (-zeta[i]+1)/2;
                                }
                            }
                          if (neigh_edge==2)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = 0;
                                  eta1DNeigh[i] = (zeta[i] +1)/2;
                                }
                            }
                          break;
                        }
                  }
                  else
                  {
                  switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                    // this is only for 1-regular triangulations
                    {                                    // points in reference cell
                    case BFUnitSquare :                  // edge 0
                      if (neigh_edge==0)
                        {
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              xi1DNeigh[i] = (-zeta[i]+part)/2;
                              eta1DNeigh[i] = -1;
                            }
                        }
                      if (neigh_edge==1)
                        {                               // edge 1
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              xi1DNeigh[i] = 1;
                              eta1DNeigh[i] = (-zeta[i]+part)/2;
                            }
                        }
                      if (neigh_edge==2)
                        {                               // edge 2
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              xi1DNeigh[i] = (zeta[i]-part)/2;
                              eta1DNeigh[i]= 1;
                            }
                        }
                      
                      if (neigh_edge==3)
                        {                               // edge 3
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              xi1DNeigh[i] = -1;
                              eta1DNeigh[i]= (zeta[i]-part)/2;
                            }
                        }
                      break;
                      
                    case BFUnitTriangle :
                      if (neigh_edge==0)
                        {
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              if (part==-1)
                                part=0;
                              xi1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                              eta1DNeigh[i] = 0;
                            }
                        }
                      if (neigh_edge==1)
                        {
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              if (part==1)
                                part=0;                              
                              xi1DNeigh[i] = ((zeta[i]+1)/2-part)/2;
                              if (part==0)
                                part=1;                              
                              if (part==-1)
                                part=0;
                              eta1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                              if (part==0)
                                part=-1;                              
                              //cout << "part " << part << endl;
                            }
                        }
                      if (neigh_edge==2)
                        {
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {
                              if (part==1)
                                part=0;
                              xi1DNeigh[i] = 0;
                              eta1DNeigh[i] = ((zeta[i] +1)/2-part)/2;
                            }
                        }
                      break;
                    }
                  }
                  if (ee_verbose>1)
                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                      cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

                                                          // compute gradients in reference cell of the neighbour
                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                    {                      
                      bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                      bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                      bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                    }
                  RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                  TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);
                  
                  DOF = GlobalNumbers + BeginIndex[neigh_N_];
                  for(i=0;i<N_Neigh;i++)
                    {
                      FEFunctValuesNeigh[i] = Values[DOF[i]];
                      if (ee_verbose>1)
                        cout << " value " <<  FEFunctValuesNeigh[i] << endl;
                    }
                  for (i=0;i<N_Points1D;i++)  // get values and derivatives in original cell
                    {
                      TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                                                 eta1DNeigh[i],
                                                 TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                                 Coll, (TGridCell *)neigh,
                                                 xietaval_refNeigh1D[BaseFunctNeigh][i],
                                                 xideriv_refNeigh1D[BaseFunctNeigh][i],
                                                 etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                                 xyval_refNeigh1D[i],
                                                 xderiv_refNeigh1D[i],
                                                 yderiv_refNeigh1D[i]);
                    }
                  
                  for(i=0;i<N_Points1D;i++)     // for all quadrature points
                    {
                      val[0]=val[1]= val[2]=0;
                      for(l=0;l<N_Neigh;l++)       // for all basis functions 
                        {
                          val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                          if (ee_verbose>1)
                            cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                              yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                        } // endfor l
                      xderiv_Neigh1D[i]= val[0]; // for k-th 
                      yderiv_Neigh1D[i]= val[1]; // for k-th 
                      xyval_Neigh1D[i]= val[2]; // for k-th 
                    } // endfor i                                        
                  
                  
                  TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                              eta1DNeigh, 
                                              X1DNeigh, Y1DNeigh, absdet1D);   
                  jump=0.0;
                  absdetjk1D = hE/2.0;
                  for (i=0;i<N_Points1D;i++)           // compute jump
                    {
                      if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
                        cout << " wrong quad points 1 " << X1D[j][i] << " , " << Y1D[j][i] 
                             << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                      if (check_cont)
                        if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
                          {
                            cout << "quad points a " << X1D[j][i] << " , " << Y1D[j][i] << endl;
                            cout << " i " << i << " vala " << xyval_1D[j][i]<< " neigha " << xyval_Neigh1D[i]<< " " << fabs(xyval_1D[j][i]-xyval_Neigh1D[i]) << endl;
                          }
                      e1 =coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
                                     +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny);
                      if (ee_verbose>1)
                        cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
                      w = weights1D[i]*absdetjk1D;
                      jump+= w*e1*e1;                       // integral on the edge
                    }
                  if (ee_verbose>1)
                    cout << "jump " << jump << endl;
                  
                  beta[0]=hE;                          // weight for H^1 estimator
                  beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                  beta[2]=hE/coeff[0];                 // weight for energy norm estimator
                  if (1.0/sqrt(coeff[0])< beta[2])
                    beta[2]=1.0/sqrt(coeff[0]);
                  for (i=1;i<4;i++)
                    estimated_error[i] += beta[i-1]*jump/2.0;
                                    
                } // end no neighbour                                        
/*************************************************************************/
/*  neighbour is not on the finest level, find children of neighbour     */
/*************************************************************************/
              else                                // there is a neighbour on the same level
                {
                  n=neigh->GetClipBoard();
                  if(n==-1)
                    {
                      // the neighbour is no member of the collection
                      // check whether the children of neigh are in collection
                      // find the local edge of neigh on which cell is -> l

                      edge2neigh=0;
                      while(neigh->GetJoint(edge2neigh)->GetNeighbour(neigh)!=cell) 
                        edge2neigh++; // find connections between cells
                      refdesc=neigh->GetRefDesc();                          // ref desc of neighbour
                      refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh);// get edges   
                      refdesc->GetOldEdgeNewEdge(TmpoEnE, TmpLen1, MaxLen1);// get connection to child edges
                      refdesc->GetEdgeChild(TmpEC, TmpLen2, MaxLen2);       // get cell belonging to child edge (TmpEC)
                      refdesc->GetOldEdgeNewLocEdge(TmpoEnlE);              // get local no.s of child edge
                      if (conform_grid)
                        N_child = 1;
                      else
                        N_child = 2;                       // not general  !!! 

                      for(k=0;k < N_child ; k++) // find children of neigh on face l -> child 
                        {
                          edge1=TmpoEnE[edge2neigh*MaxLen1+k];                 // edge child, not general !!!               
                          chnum1=TmpEC[edge1*MaxLen2];                // local number of child cell 
                          child =neigh->GetChild(chnum1);             // child cell   
                          child_N_=child->GetClipBoard();             // id of child cell
                          refdesc->GetEdgeChildIndex(TmpECI,TmpLen3, MaxLen3); // get local indices of child edge
                          l_child = TmpECI[edge1*MaxLen3];                        // local index of child edge

                          refdesc_child=child->GetRefDesc();          // ref desc of child 
                          refdesc_child->GetShapeDesc()->GetEdgeVertex(TmpEdVer);// conn. edge -> vertices
                          ver2 = child->GetVertex(TmpEdVer[2*l_child]);          // vertices of edge
                          ver3 = child->GetVertex(TmpEdVer[2*l_child+1]);
                          
                          if (ee_verbose>1)
                            {
                              cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                              cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                              cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                              cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
                            }
                          
                          if (ver1==ver2)
                            {
                              part=1;
                            }
                          else if (ver0==ver3)
                            {
                              part=-1;
                            }
                          else 
                            {
                              cout << " something wrong 5 " <<  endl;
                              cout << "ver 0 " << ver0->GetX() << "  " << ver0->GetY() << endl;
                              cout << "ver 1 " << ver1->GetX() << "  " << ver1->GetY() << endl;
                              cout << "ver 2 " << ver2->GetX() << "  " << ver2->GetY() << endl;
                              cout << "ver 3 " << ver3->GetX() << "  " << ver3->GetY() << endl;
                            }
                          // now from point of view of child cell -> cell becomes the neighbour
                          // prepare intergration for the half part of edge j


                          neigh_N_ = cell->GetClipBoard();           // number of original cell  in iterator
                          if(neigh_N_==-1)
                            {
                              cout << "Hier sollte man aber nicht hinkommen 33 !"<< endl;
                            }    
                          CurrEleNeigh = fespace->GetFE2D(neigh_N_,cell);   // finite element on neighbour
                          eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                  
                          BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                          N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions
                  
                          bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                          bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                          neigh_edge = j;
                          if (conform_grid)
                          {
                            switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                              {                                    // points in reference cell
                              case BFUnitSquare :                  // edge 0
                                if (neigh_edge==0)
                                  {
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = -zeta[i];
                                        eta1DNeigh[i] = -1;
                                      }
                                  }
                                if (neigh_edge==1)
                                  {                               // edge 1
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = 1;
                                        eta1DNeigh[i] = -zeta[i];
                                      }
                                  }
                                if (neigh_edge==2)
                                  {                               // edge 2
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = zeta[i];
                                        eta1DNeigh[i]= 1;
                                      }
                                  }
                                
                                if (neigh_edge==3)
                                  {                               // edge 3
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = -1;
                                        eta1DNeigh[i]= zeta[i];
                                      }
                                  }
                                break;
                                
                              case BFUnitTriangle :
                                if (neigh_edge==0)
                                  {
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = (-zeta[i]+1)/2;
                                        eta1DNeigh[i] = 0;
                                      }
                                  }
                                if (neigh_edge==1)
                                  {
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = (zeta[i]+1)/2;
                                        eta1DNeigh[i] = (-zeta[i]+1)/2;
                                      }
                                  }
                                if (neigh_edge==2)
                                  {
                                    for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                      {
                                        xi1DNeigh[i] = 0;
                                        eta1DNeigh[i] = (zeta[i] +1)/2;
                                      }
                                  }
                                break;
                              }
                          }
                          else
                          {
                            switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                              // this is only for 1-regular triangulations
                              {                                    // points in reference cell
                              case BFUnitSquare :                  // edge 0
                              if (neigh_edge==0)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = (-zeta[i]+part)/2;
                                      eta1DNeigh[i] = -1;
                                    }
                                }
                              if (neigh_edge==1)
                                {                               // edge 1
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = 1;
                                      eta1DNeigh[i] = (-zeta[i]+part)/2;
                                    }
                                }
                              if (neigh_edge==2)
                                {                               // edge 2
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = (zeta[i]-part)/2;
                                      eta1DNeigh[i]= 1;
                                    }
                                }
                              
                              if (neigh_edge==3)
                                {                               // edge 3
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = -1;
                                      eta1DNeigh[i]= (zeta[i]-part)/2;
                                    }
                                }
                              break;
                      
                            case BFUnitTriangle :
                              if (neigh_edge==0)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      if (part==-1)
                                        part=0;
                                      xi1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                                      eta1DNeigh[i] = 0;
                                    }
                                }
                              if (neigh_edge==1)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      if (part==1)
                                        part=0;                              
                                      xi1DNeigh[i] = ((zeta[i]+1)/2-part)/2;
                                      if (part==0)
                                        part=1;                              
                                      if (part==-1)
                                        part=0;
                                      eta1DNeigh[i] = ((-zeta[i]+1)/2+part)/2;
                                      if (part==0)
                                        part=-1;                              
                                    }
                                }
                              if (neigh_edge==2)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      if (part==1)
                                        part=0;
                                      xi1DNeigh[i] = 0;
                                      eta1DNeigh[i] = ((zeta[i] +1)/2-part)/2;
                                    }
                                }
                              break;
                            }
                          }
                          if (ee_verbose>1)
                            for (i=0;i<N_Points1D;i++)         // for all quadrature points
                              cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;


                          // compute gradients in reference cell of the neighbour
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {                      
                              bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                              bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                              bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                            }
                          RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                          TFEDatabase2D::SetCellForRefTrans(cell,RefTransNeigh);

                          DOF = GlobalNumbers + BeginIndex[neigh_N_];
                          for(i=0;i<N_Neigh;i++)
                            {
                              FEFunctValuesNeigh[i] = Values[DOF[i]];
                              if (ee_verbose>1)
                                cout << " value " <<  FEFunctValuesNeigh[i] << endl;
                            }

                          for (i=0;i<N_Points1D;i++)  // get values and derivatives in original cell
                            {
                              TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                                                         eta1DNeigh[i],
                                                         TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                                         Coll, (TGridCell *)neigh,
                                                         xietaval_refNeigh1D[BaseFunctNeigh][i],
                                                         xideriv_refNeigh1D[BaseFunctNeigh][i],
                                                         etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                                         xyval_refNeigh1D[i],
                                                         xderiv_refNeigh1D[i],
                                                         yderiv_refNeigh1D[i]);
                            }
                          
                          for(i=0;i<N_Points1D;i++)     // for all quadrature points
                            {
                              val[0]=val[1]= val[2]=0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  if (ee_verbose>1)
                                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                                } // endfor l
                              xderiv_Cell1D[i]= val[0]; // for k-th 
                              yderiv_Cell1D[i]= val[1]; // for k-th 
                              xyval_Cell1D[i]= val[2]; // for k-th 
                            } // endfor i                                        

                          TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                      eta1DNeigh, 
                                                      X1DCell, Y1DCell, absdet1D);   
                          // prepare integration for the child of the neighbour belong to the half part
                          // of edge j

                           neigh_N_ = child->GetClipBoard();                  // number of neighbour in iterator
                          CurrEleNeigh = fespace->GetFE2D(neigh_N_,child);   // finite element on neighbour
                          eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                          
                          BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                          N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions

                          bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                          bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour


                          neigh_edge = l_child;
                          switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                            {                                    // points in reference cell
                            case BFUnitSquare :                  // edge 0
                              if (neigh_edge==0)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = zeta[i];
                                      eta1DNeigh[i] = -1;
                                    }
                                }
                              if (neigh_edge==1)
                                {                               // edge 1
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = 1;
                                      eta1DNeigh[i] = zeta[i];
                                    }
                                }
                              if (neigh_edge==2)
                                {                               // edge 2
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = -zeta[i];
                                      eta1DNeigh[i]= 1;
                                    }
                                }
                              
                              if (neigh_edge==3)
                                {                               // edge 3
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = -1;
                                      eta1DNeigh[i]= -zeta[i];
                                    }
                                }
                              break;
                              
                            case BFUnitTriangle :
                              if (neigh_edge==0)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = (zeta[i]+1)/2;
                                      eta1DNeigh[i] = 0;
                                    }
                                }
                              if (neigh_edge==1)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = (-zeta[i]+1)/2;
                                      eta1DNeigh[i] = (zeta[i]+1)/2;
                                    }
                                }
                              if (neigh_edge==2)
                                {
                                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                    {
                                      xi1DNeigh[i] = 0;
                                      eta1DNeigh[i] = (-zeta[i] +1)/2;
                                    }
                                }
                              break;
                            }
                          
                          if (ee_verbose>1)
                            for (i=0;i<N_Points1D;i++)         // for all quadrature points
                              cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;


                          // compute gradients in reference cell of the neighbour
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {                      
                              bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                              bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                              bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                            }
                          RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                          TFEDatabase2D::SetCellForRefTrans(child,RefTransNeigh);
                          
                          DOF = GlobalNumbers + BeginIndex[neigh_N_];
                          for(i=0;i<N_Neigh;i++)
                            {
                              FEFunctValuesNeigh[i] = Values[DOF[i]];
                              if (ee_verbose>1)
                                cout << " value " <<  FEFunctValuesNeigh[i] << endl;
                            }

                          for (i=0;i<N_Points1D;i++)  // get values and derivatives in original cell
                            {
                              TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                                                         eta1DNeigh[i],
                                                         TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                                         Coll, (TGridCell *)neigh,
                                                         xietaval_refNeigh1D[BaseFunctNeigh][i],
                                                         xideriv_refNeigh1D[BaseFunctNeigh][i],
                                                         etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                                         xyval_refNeigh1D[i],
                                                         xderiv_refNeigh1D[i],
                                                         yderiv_refNeigh1D[i]);
                            }
                          

                          for(i=0;i<N_Points1D;i++)     // for all quadrature points
                            {
                              val[0]=val[1]= val[2]=0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  if (ee_verbose>1)
                                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                                } // endfor l
                              xderiv_Neigh1D[i]= val[0]; // for k-th 
                              yderiv_Neigh1D[i]= val[1]; // for k-th 
                              xyval_Neigh1D[i]= val[2]; // for k-th 
                            } // endfor i                                        
                          
                          TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                      eta1DNeigh, 
                                                      X1DNeigh, Y1DNeigh, absdet1D);   
                          jump=0.0;
                          absdetjk1D = hE/(2.0*N_child);       // only half edge is considered
                           for (i=0;i<N_Points1D;i++)           // compute jump
                            {                              
                              if ((fabs(X1DCell[i]-X1DNeigh[i])+fabs(Y1DCell[i]-Y1DNeigh[i]))>1e-8)
                                cout << " wrong quad points 2 " << X1DCell[i] << " , " << Y1DCell[i] 
                                     << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                              if (check_cont)
                                if (fabs(xyval_Neigh1D[i]-xyval_Cell1D[i])>1e-8)
                               {
                                 cout << "quad points b " << X1DCell[i] << " , " << Y1DCell[i] << endl;
                                 cout << " i " << i << " valb " << xyval_Cell1D[i]<< " neighb " << xyval_Neigh1D[i] << " " << fabs(xyval_Cell1D[i]-xyval_Neigh1D[i]) << endl;
                                }
                              e1 =coeff[0] *((xderiv_Cell1D[i]-xderiv_Neigh1D[i])*nx
                                             +(yderiv_Cell1D[i]-yderiv_Neigh1D[i])*ny);
                              if (ee_verbose>1)
                                cout << i<< " jumpx " << xderiv_Cell1D[i] << " " << xderiv_Neigh1D[i] << endl;
                              w = weights1D[i]*absdetjk1D;
                              jump+= w*e1*e1;                       // integral on the edge
                            }
                          if (ee_verbose>1)
                            cout << "jump " << jump << endl;
                          hE2=hE/N_child;
                          beta[0]=hE2;                          // weight for H^1 estimator
                          beta[1]=hE2*hE2*hE2;                    // weight for L^2 estimator
                          beta[2]=hE2/coeff[0];                 // weight for energy norm estimator
                          if (1.0/sqrt(coeff[0])< beta[2])
                            beta[2]=1.0/sqrt(coeff[0]);
                          for (i=1;i<4;i++)
                            estimated_error[i] += beta[i-1]*jump/2.0;

                        }
                    } // end clipboard==-1
                  else
/*************************************************************************/
/*  neighbour is on the finest level                                     */
/*************************************************************************/
                    
                    {              // the neighbour is a member of the collection
                                   // find the finite element on the other side
                      //                      cout << " neighbour found " << endl; 
                      // find the local edge of neigh on which cell is -> l
                      neigh_edge=0;
                      while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;
                      refdesc=neigh->GetRefDesc();
                      refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh); 
                      ver0=  cell->GetVertex(TmpEdVer[2*j]);
                      ver1=  cell->GetVertex(TmpEdVer[2*j+1]);
                      ver2 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge]);          // vertices of edge
                      ver3 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge+1]);
                      if(!(((ver0==ver2)&&(ver1==ver3))||((ver0==ver3)&&(ver1==ver2))))
                      {
                        OutPut("wrong edge 2" << endl);
                        exit(4711);
                      }
                      
                      // compute gradient at the quadrature points on the edge of 
                      // the neighbour element
                      
                      neigh_N_ = neigh->GetClipBoard();                  // number of neighbour in iterator
                      CurrEleNeigh = fespace->GetFE2D(neigh_N_,neigh);   // finite element on neighbour
                      eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                      
                      BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                      N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions

                      bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                      bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                      switch(bf2DrefelementsNeigh)                // compute coordinates of line quadrature
                        {                                    // points in reference cell
                        case BFUnitSquare :                  // edge 0
                          if (neigh_edge==0)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = -zeta[i];
                                  eta1DNeigh[i] = -1;
                                }
                            }
                          if (neigh_edge==1)
                            {                               // edge 1
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = 1;
                                  eta1DNeigh[i] = -zeta[i];
                                }
                            }
                          if (neigh_edge==2)
                            {                               // edge 2
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = zeta[i];
                                  eta1DNeigh[i]= 1;
                                }
                            }

                          if (neigh_edge==3)
                            {                               // edge 3
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = -1;
                                  eta1DNeigh[i]= zeta[i];
                                }
                            }
                          break;

                        case BFUnitTriangle :
                           if (neigh_edge==0)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = (-zeta[i]+1)/2;
                                  eta1DNeigh[i] = 0;
                                }
                            }
                          if (neigh_edge==1)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = (zeta[i]+1)/2;
                                  eta1DNeigh[i] = (-zeta[i]+1)/2;
                                }
                            }
                          if (neigh_edge==2)
                            {
                              for (i=0;i<N_Points1D;i++)         // for all quadrature points
                                {
                                  xi1DNeigh[i] = 0;
                                  eta1DNeigh[i] = (zeta[i] +1)/2;
                                }
                            }
                          break;
                        }
                      
                      if (ee_verbose>1)
                        for (i=0;i<N_Points1D;i++)         // for all quadrature points
                          cout << "xiN " << xi1DNeigh[i] << " etaN " << eta1DNeigh[i] << endl;

                      // compute gradients in reference cell of the neighbour
                      for (i=0;i<N_Points1D;i++)         // for all quadrature points
                        {                      
                          bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                          bfNeigh->GetDerivatives(D10,xi1DNeigh[i],eta1DNeigh[i],xideriv_refNeigh1D[BaseFunctNeigh][i]);
                          bfNeigh->GetDerivatives(D01,xi1DNeigh[i],eta1DNeigh[i],etaderiv_refNeigh1D[BaseFunctNeigh][i]);
                        }
                      RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                      TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);

                      DOF = GlobalNumbers + BeginIndex[neigh_N_];
                      for(i=0;i<N_Neigh;i++)
                        {
                          FEFunctValuesNeigh[i] = Values[DOF[i]];
                          if (ee_verbose>1)
                            cout << " value " <<  FEFunctValuesNeigh[i] << endl;
                        }
                      for (i=0;i<N_Points1D;i++)  // get values and derivatives in original cell
                        {
                          TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1DNeigh[i],
                                                     eta1DNeigh[i],
                                                     TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                                                     Coll, (TGridCell *)neigh,
                                                     xietaval_refNeigh1D[BaseFunctNeigh][i],
                                                     xideriv_refNeigh1D[BaseFunctNeigh][i],
                                                     etaderiv_refNeigh1D[BaseFunctNeigh][i],
                                                     xyval_refNeigh1D[i],
                                                     xderiv_refNeigh1D[i],
                                                     yderiv_refNeigh1D[i]);
                        }

                      for(i=0;i<N_Points1D;i++)     // for all quadrature points
                        {
                          val[0]=val[1]= val[2]=0;
                          for(l=0;l<N_Neigh;l++)       // for all basis functions 
                            {
                              val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                              if (ee_verbose>1)
                                cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                  yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                            } // endfor l
                          xderiv_Neigh1D[i]= val[0]; // for k-th 
                          yderiv_Neigh1D[i]= val[1]; // for k-th 
                          xyval_Neigh1D[i]= val[2]; // for k-th 
                        } // endfor i                                        
                      
                      TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                  eta1DNeigh, 
                                                  X1DNeigh, Y1DNeigh, absdet1D);   
                  
                      jump=0.0;
                      absdetjk1D = hE/2.0;
                      for (i=0;i<N_Points1D;i++)           // compute jump
                        {
                          if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
                            cout << " wrong quad points 0 " << X1D[j][i] << " , " << Y1D[j][i] 
                                 << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                          if (check_cont)
                            if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
                              cout << " i " << i << " valc " << xyval_1D[j][i]<< " neighc " << xyval_Neigh1D[i]<< endl;
                          e1 =coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
                                          +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny);
                          if (ee_verbose>1)
                            cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
                          w = weights1D[i]*absdetjk1D;
                         jump+= w*e1*e1;                       // integral on the edge
                        }
                      if (ee_verbose>1)
                        cout << "jump " << jump << endl;
                      beta[0]=hE;                          // weight for H^1 estimator
                      beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                      beta[2]=hE/coeff[0];                 // weight for energy norm estimator
                      if (1.0/sqrt(coeff[0])< beta[2])
                        beta[2]=1.0/sqrt(coeff[0]);
                      for (i=1;i<4;i++)
                        estimated_error[i] += beta[i-1]*jump/2.0;
                                          
                      
                    }  // end neighbour is member of the collection
                } // end neighbour on the finer level 
            }     // end inner edge
        }         // end for j
       delete xi1DNeigh; 
       delete eta1DNeigh;
       delete X1DNeigh; 
       delete Y1DNeigh;
       delete X1DCell; 
       delete Y1DCell;
       delete FEFunctValuesNeigh;
       for (i=0;i<N_Points1D;i++)
         {
           delete xyval_refNeigh1D[i];
           delete  xderiv_refNeigh1D[i];
           delete   yderiv_refNeigh1D[i];
         }
       delete xderiv_Neigh1D;
       delete  yderiv_Neigh1D;
       delete  xyval_Neigh1D;
       delete xderiv_Cell1D;
       delete yderiv_Cell1D;
       delete xyval_Cell1D;
   }             // end residual based estimators
  

  //for (i=0;i<4;i++)
  //  cout << i << " estimated error " << estimated_error[i] << endl;
  for (i=0;i<5;i++)
    estimated_local_error[i]=estimated_error[i];
}

#endif
