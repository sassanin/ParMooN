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
// @(#)NS2DErrorEstimator.C        1.4 10/26/99
// 
// Class:       TNS2DErrorEstimator
//
// Purpose:     compute a posteriori error estimate
//
// Author:      Volker John
//
// History:     04.02.1998 start implementation
//
// =======================================================================
#ifdef __2D__

#include <NSE2DErrorEstimator.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <Enumerations.h>
#include <BaseFunct2D.h>
#include <Database.h>
#include <MooNMD_Io.h>

#include <math.h>
#include <string.h>
#include <stdlib.h>
#ifdef __MAC64__
#include <malloc/malloc.h>
#else
#include <malloc.h>
#endif

TNS2DErrorEstimator::TNS2DErrorEstimator( int fe_local_estimator,
                                      TFEVectFunct2D *u,
                                      TFEFunction2D *p,
                                      int error_control,
                                      int navierstokes)
{
  FELocalEstimator = fe_local_estimator;
  FESpace2D_U =  (TFESpace2D *) u->GetFESpace2D();
  Collection_U = FESpace2D_U->GetCollection();
  FESpace2D_P =  (TFESpace2D *) p->GetFESpace2D();
  Collection_P = FESpace2D_P->GetCollection();
  U = u;
  P = p;
  ErrorControl = error_control;
  NavierStokes = navierstokes;
}

void TNS2DErrorEstimator::GetErrorEstimate(int N_Derivatives,
                                         MultiIndex2D *NeededDerivatives,
                                         int N_DerivativesP,
                                         MultiIndex2D *NeededDerivativesP,
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
  int i,j,k,l,m,n,N_UsedElements, N_LocalUsedElements;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_U, N_P;
  int Used[N_FEs2D],current_estimator;
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement,CurrentElementP ;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunct,BaseFunctP;
  TBaseFunct2D *bf;
  TCollection *Coll;
  TBaseCell *cell, *neigh;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta,*weights1D, *zeta;
  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *xyval_1D[4];
  double *xderiv_1D[4];
  double *yderiv_1D[4];
  double *X1D[4], *Y1D[4], val[6];
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  RefTrans2D RefTrans;
  double *Param[MaxN_QuadPoints_2D], *aux;
  double *Derivatives[3*MaxN_QuadPoints_2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, N_DOF, *DOFP;
  double **OrigFEValues, *Orig, value[2];
  double FEFunctValues[3*MaxN_BaseFunctions2D];
  double *Values,*ValuesP,max_loc_err;
  int *GlobalNumbers, *GlobalNumbersP, *BeginIndex, *BeginIndexP;
  double LocError[4];
  double estimated_global_errors[4], estimated_local_errors[4];
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  bool *SecondDer;

  int ee_verbose=2;                             // verbosity

  int memory[3],data_base_memory;
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;
 info = mstats();
 
 memory[0]=memory[1]=memory[2]=0.;
#else    
  struct mallinfo MALLINFO;
  MALLINFO = mallinfo();
  memory[0]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif
// ########################################################################
// store information in local arrays
// ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  SecondDer = new bool[n_fespaces];

  memset(Used, 0, N_FEs2D*SizeOfInt);

  for(i=0;i<n_fespaces;i++)
  {
    fespace = fespaces[i];                      /* fe space */
    n = fespace->GetN_UsedElements();           /* # used finite elements */
    UsedElements = fespace->GetUsedElements();  /* used finite elements */
    for(j=0;j<n;j++)                            /* for all finite elements */
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
    }                                           // enfor j
  }                                             // endfor i

  N_UsedElements = 0;                           /* compute number of used elements */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];      /* store used finite elements */
  j=0;                                          /* in array */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
    {
      UsedElements[j] = (FE2D)i;
      j++;
    }                                           // endif

  if (ee_verbose>1)
    {
      cout << "estimator number of used elements: " << N_UsedElements << endl;
      for(i=0;i<N_UsedElements;i++)
        cout << "UsedElements[" << i << "]: " << UsedElements[i] << endl;
    }

// ########################################################################
// calculate values of base functions and derivatives on ref element
// ########################################################################
  CurrentElement = (fespaces[0]->GetUsedElements())[0]; // fe for velocity
  l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement); // get 1d quad fromula
  LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
  qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
  qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

  for(i=0;i<N_UsedElements;i++)                 // for used finite elements
  {
    CurrentElement = UsedElements[i];
    /*l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);*/
    BaseFunct = BaseFuncts[CurrentElement];     

    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunct); // get base functions
    bf2Drefelements = bf->GetRefElement();
    
    switch(bf2Drefelements)                      // compute coordinates of line quadrature
    {                                            // points in reference cell
                                                 // quadrilateral cell 
    case BFUnitSquare :                          // edge 0
      for (j=0;j<N_Points1D;j++)                 // for all quadrature points
        {
          xi1D[BaseFunct][0][j] = zeta[j];
          eta1D[BaseFunct][0][j] = -1;
          bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunct][0][j]);
        }                                        // edge 1
      for (j=0;j<N_Points1D;j++)               // for all quadrature points
        {       
          xi1D[BaseFunct][1][j] = 1;
          eta1D[BaseFunct][1][j] = zeta[j];
          bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunct][1][j]);
        }                                        // edge 2
      for (j=0;j<N_Points1D;j++)                 // for all quadrature points
        {
          xi1D[BaseFunct][2][j] = -zeta[j];
          eta1D[BaseFunct][2][j] = 1;
          bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunct][2][j]);
        }                                         // edge 3
      for (j=0;j<N_Points1D;j++)                  // for all quadrature points
        {
          xi1D[BaseFunct][3][j] = -1;
          eta1D[BaseFunct][3][j] = -zeta[j];
          bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunct][3][j]);
          bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunct][3][j]);
          bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunct][3][j]);
        }
      break;
      
    case BFUnitTriangle :                        // triangular cell 
      for (j=0;j<N_Points1D;j++)                 // for all quadrature points
        {
          xi1D[BaseFunct][0][j] = (zeta[j]+1)/2;
          eta1D[BaseFunct][0][j] = 0;
          bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunct][0][j]);
          bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunct][0][j]);
        }                                       // edge 1
      for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct][1][j] = (-zeta[j]+1)/2;
          eta1D[BaseFunct][1][j] = (zeta[j]+1)/2;
          bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunct][1][j]);
          bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunct][1][j]);
        }                                       // edge 2
      for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunct][2][j] = 0;
          eta1D[BaseFunct][2][j] = (-zeta[j] +1)/2;
          bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunct][2][j]);
          bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunct][2][j]);
        }
      break;
    }
  }                                               // endfor i

  for (i=0;i<4;i++)                               // arrays for coordinates, values and  
    {                                             // determinant for 1D quadrature
      X1D[i] = new double[N_Points1D];            // coordinates of edge i  
      Y1D[i] = new double[N_Points1D];
      AbsDetjk1D[i] = new double[MaxN_QuadPoints_2D];     // determinant of affine mapping 
      for (j=0;j<N_Points1D;j++)                  // arrays for values in reference cell
        {
          xyval_ref1D[i][j] = new double[MaxN_BaseFunctions2D];
          xderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D];
          yderiv_ref1D[i][j] = new double[MaxN_BaseFunctions2D];
        }
      xyval_1D[i] = new double[3*N_Points1D];       // arrays for values in original cell 
      xderiv_1D[i] = new double[3*N_Points1D];
      yderiv_1D[i] = new double[3*N_Points1D];
    }

  N_Parameters = Aux->GetN_Parameters();          // get number of parameters of equation
  aux = new double [MaxN_QuadPoints_2D*N_Parameters]; // allocate memory for Param arrays
  for(j=0;j<MaxN_QuadPoints_2D;j++)               // set pointers
    Param[j] = aux + j*N_Parameters;

  aux = new double [3*MaxN_QuadPoints_2D*N_Derivatives]; // allocate memory for Derivatives arrays, hier kann man was sparen
  for(j=0;j<3*MaxN_QuadPoints_2D;j++)               // set pointers
    Derivatives[j] = aux + j*N_Derivatives;
 
  // 20 <= number of term
  aux = new double [MaxN_QuadPoints_2D*20]; // allocate memory for AuxArray
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;
 
  GlobalNumbers = FESpace2D_U->GetGlobalNumbers();
  BeginIndex = FESpace2D_U->GetBeginIndex();
  GlobalNumbersP = FESpace2D_P->GetGlobalNumbers();
  BeginIndexP = FESpace2D_P->GetBeginIndex();

// ########################################################################
// prepare error estimates
// ########################################################################

  // all spaces use same Coll
  Coll = FESpace2D_U->GetCollection();          // collection of mesh cells
  N_Cells = Coll->GetN_Cells();                 // number of mesh cells
  N_U = U->GetLength();
  N_DOF = 2*U->GetLength()+P->GetLength();      // number of global dof 
  Values = U->GetValues();                      // values of fe function
  ValuesP = P->GetValues();

  for(i=0;i<N_Cells;i++)                        // do for all mesh cells
  {                                           // on the finest level 
    cell = Coll->GetCell(i);
    k=cell->GetN_Edges();                       // # edges
    for(j=0;j<k;j++)                            // for all edges
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell); // neighbour cell
      if(neigh) neigh->SetClipBoard(-1);        // set clipboard to -1
    }
    cell->SetClipBoard(-1);          
  }                                             // endfor i
  // non finest neighbours of finest cells have clipboard -1 

  for(i=0;i<N_Cells;i++)                        // set clipboard of cells on finest  
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  for (i=0;i<4;i++)                           // initialize some quantities
    estimated_global_errors[i]=0.0;
  max_loc_err = 0;
  data_base_memory =0;
  current_estimator = GetFELocalEstimator();


// ########################################################################
// loop over all cells
// ########################################################################
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
 
#else  
  MALLINFO = mallinfo();
  memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  for(i=0;i<N_Cells;i++)                      // for all cells on the finest level
  {
     cell = Coll->GetCell(i);                        // next cell 
    eta_K[i] = 0.0;                           // initialize local estimate

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
   for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i,cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement]; // local basis functions
      LocBF[j] = BaseFuncts[CurrentElement];
      if (j==0)
        SecondDer[j] = TRUE;                      // with 2nd derivative
      else                                        // for velo
        SecondDer[j] = FALSE;
    }
    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    
                                      // get reference transformation
    RefTrans = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, weights, X, Y, AbsDetjk);
    if(N_Parameters>0)                // get parameters of equ.
        Aux->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param); 
    

    // velocity first
                                       // calculate all needed derivatives of this FE function
    CurrentElement = FESpace2D_U->GetFE2D(i,cell);  // finite element on cell
    BaseFunct = BaseFuncts[CurrentElement];       // basis functions
    N_ = N_BaseFunct[CurrentElement];             // # basis functions    
    DOF = GlobalNumbers + BeginIndex[i];    // dof of current mesh cell
    // cout << "cell no " << i << " begin " << *DOF << endl; 
    for(l=0;l<N_;l++)
      {
        FEFunctValues[l] = Values[DOF[l]];    // u fe values of dofs
        FEFunctValues[l+MaxN_BaseFunctions2D] = Values[DOF[l]+N_U];    // v fe values of dofs
        //cout << FEFunctValues[l] << " " << FEFunctValues[l+MaxN_BaseFunctions2D] << endl;
      }

                                         // compute values for all derivatives 
                                            // in all quadrature points
                                            // in original mesh cell 
    for(k=0;k<N_Derivatives;k++)            // for all derivatives 
    {                                       // get values in original cell   
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)               // for all quadrature points
      {
        Orig = OrigFEValues[j];             // value in original cell 
        value[0] = value[1] = 0;
        for(l=0;l<N_;l++)                   // for all basis functions 
          {
            value[0] += FEFunctValues[l] * Orig[l]; // accumulate u value of derivative in point j
            value[1] += FEFunctValues[l+MaxN_BaseFunctions2D] * Orig[l]; // accumulate v value of derivative in point j
          }
        Derivatives[j][k] = value[0];       // for k-th u derivative
        Derivatives[j+MaxN_QuadPoints_2D][k] = value[1];// for k-th v derivative
      }                                     // endfor j
    }                                       // endfor k
          
    if(Coeff)                               // get coefficients of pde 
      Coeff(N_Points, X, Y, Param, AuxArray); // 0 - eps, 1,2 - rhs
                                            // prepare 1D quadrature formula
    // pressure first
    CurrentElementP = FESpace2D_P->GetFE2D(i,cell);  // finite element on cell
    BaseFunctP = BaseFuncts[CurrentElementP];       // basis functions
    N_P = N_BaseFunct[CurrentElementP];             // # basis functions    
    DOFP = GlobalNumbersP + BeginIndexP[i];    // dof of current mesh cell
    //cout << "cell no " << i << " begin " << *DOFP << 
    // " " << *GlobalNumbersP << " " << BeginIndexP[i] << endl; 
    for(l=0;l<N_P;l++)
      {
        FEFunctValues[l+2*MaxN_BaseFunctions2D] = ValuesP[DOFP[l]];    // p fe values of dofs
        //cout <<  FEFunctValues[l+2*MaxN_BaseFunctions2D] << endl;
      }
   
    for(k=0;k<N_DerivativesP;k++)            // for all derivatives 
    {                                       // get values in original cell   
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunctP,NeededDerivativesP[k]);
      for(j=0;j<N_Points;j++)               // for all quadrature points
      {
        Orig = OrigFEValues[j];             // value in original cell 
        value[0] = 0;
        for(l=0;l<N_P;l++)                   // for all basis functions 
          {
            value[0] += FEFunctValues[l+2*MaxN_BaseFunctions2D] * Orig[l]; // accumulate p value of derivative in point j
          }
        Derivatives[j+2*MaxN_QuadPoints_2D][k] = value[0];// for k-th v derivative
      }                                     // endfor j
    }                                       // endfor k

    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);

    N_Edges=cell->GetN_Edges();             // # edges
     // pressure second 
    
    for(j=0;j<N_Edges;j++)                  // loop over all edges of cell
      {                                     // get original coordinates of edge quad. points
        TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctP][j], 
                                    eta1D[BaseFunctP][j], 
                                    X1D[j], Y1D[j], AbsDetjk1D[j]);   
        for(k=0;k<N_Points1D;k++)  // get values and derivatives in original cell
          {
            TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctP][j][k],
                                       eta1D[BaseFunctP][j][k],
                                       TFEDatabase2D::GetBaseFunct2D(BaseFunctP),
                                       Coll, (TGridCell *)cell,
                                       xietaval_ref1D[BaseFunctP][j][k],
                                       xideriv_ref1D[BaseFunctP][j][k],
                                       etaderiv_ref1D[BaseFunctP][j][k],
                                       xyval_ref1D[j][k],
                                       xderiv_ref1D[j][k],
                                       yderiv_ref1D[j][k]);
          }
        
        for(k=0;k<N_Points1D;k++)     // for all quadrature points
          {
            val[0]= 0;
            for(l=0;l<N_P;l++)       // for all basis functions 
              {
                m = l+2*MaxN_BaseFunctions2D;
                val[0] += FEFunctValues[m] * xyval_ref1D[j][k][l]; // accumulate value of derivative 
              } // endfor l
            m = k + 2*N_Points1D;
            xyval_1D[j][m]= val[0]; // for k-th
          } // endfor k                                        
      }     // endfor j
   // velocity second
    for(j=0;j<N_Edges;j++)                  // loop over all edges of cell
      {                                     // get original coordinates of edge quad. points
        TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunct][j], 
                                    eta1D[BaseFunct][j], 
                                    X1D[j], Y1D[j], AbsDetjk1D[j]);   
        for(k=0;k<N_Points1D;k++)  // get values and derivatives in original cell
          {
            TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunct][j][k],
                                       eta1D[BaseFunct][j][k],
                                       TFEDatabase2D::GetBaseFunct2D(BaseFunct),
                                       Coll, (TGridCell *)cell,
                                       xietaval_ref1D[BaseFunct][j][k],
                                       xideriv_ref1D[BaseFunct][j][k],
                                       etaderiv_ref1D[BaseFunct][j][k],
                                       xyval_ref1D[j][k],
                                       xderiv_ref1D[j][k],
                                       yderiv_ref1D[j][k]);
          }
        
        for(k=0;k<N_Points1D;k++)     // for all quadrature points
          {
            val[0]=val[1]=val[2] = val[3]=val[4]=val[5]= 0;
            for(l=0;l<N_;l++)       // for all basis functions 
              {
                m = l+MaxN_BaseFunctions2D;
                val[0] += FEFunctValues[l] * xyval_ref1D[j][k][l]; // accumulate value of derivative 
                val[1] += FEFunctValues[l] * xderiv_ref1D[j][k][l]; // accumulate value of derivative
                val[2] += FEFunctValues[l] * yderiv_ref1D[j][k][l]; // accumulate value of derivative
                val[3] += FEFunctValues[m] * xyval_ref1D[j][k][l]; // accumulate value of derivative 
                val[4] += FEFunctValues[m] * xderiv_ref1D[j][k][l]; // accumulate value of derivative
                val[5] += FEFunctValues[m] * yderiv_ref1D[j][k][l]; // accumulate value of derivative
              } // endfor l
            m = k + N_Points1D;
            xyval_1D[j][k]= val[0]; // for k-th
            xderiv_1D[j][k]= val[1]; // for k-th 
            yderiv_1D[j][k]= val[2]; // for k-th 
            xyval_1D[j][m]= val[3]; // for k-th
            xderiv_1D[j][m]= val[4]; // for k-th 
            yderiv_1D[j][m]= val[5]; // for k-th 
          } // endfor k                                        
      }     // endfor j


    // estimate local errors 
    EstimateCellError(fespaces,cell,N_Points, X, Y, AbsDetjk, weights, Derivatives, 
                      AuxArray,BoundaryConds,BoundaryValues,
                      N_Points1D, zeta, X1D, Y1D, weights1D,
                      xyval_1D, xderiv_1D, yderiv_1D,
                      GlobalNumbers,BeginIndex,DOF,Values,
                      GlobalNumbersP,BeginIndexP,DOF,ValuesP,
                      estimated_local_errors);
    
    for(k=0;k<4;k++) // update global error estimates
      estimated_global_errors[k]+=estimated_local_errors[k];
                     // update maximal local error estimate
    if (estimated_local_errors[current_estimator]>max_loc_err)
      max_loc_err = estimated_local_errors[current_estimator];
    
    eta_K[i] = estimated_local_errors[current_estimator];
  }                 // endfor i (loop over cells)
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
 
#else    
  MALLINFO = mallinfo();   
  memory[2]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif  
  data_base_memory+= memory[2]-memory[1];
  
  for(i=1;i<4;i++)  // compute global error estimates
    estimated_global_error[i]=sqrt(estimated_global_errors[i]);
  estimated_global_error[0]=estimated_global_errors[0];
                    // compute maximal local error estimate
  *maximal_local_error = sqrt(max_loc_err);
                    // set memory free  

  delete UsedElements;
  delete SecondDer;
 
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
  delete Param[0];
  delete Derivatives[0];
   delete AuxArray[0];
  
#ifdef _MALLOC_MALLOC_H_
 info = mstats();
 
#else  
   MALLINFO = mallinfo();
   memory[1]=MALLINFO.usmblks+MALLINFO.uordblks;
#endif   
   if (memory[1]- memory[0]!=data_base_memory)
     cout << "WARNING : Error Estimator did not set all memory free !!!" <<  memory[1]- memory[0] << 
       "  " << data_base_memory << endl;
  
} // TCDErrorEstimator::GetErrorEstimate


void  TNS2DErrorEstimator::EstimateCellError(TFESpace2D **fespaces,
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
                                           int *GlobalNumbersP, 
                                           int *BeginIndexP,
                                           int *DOFP,
                                           double *ValuesP,                                           
                                           double *estimated_local_error)
{
  int i,j,k,l,n,N_Edges,comp,parent_edge,MaxLen1,MaxLen2,MaxLen3,N_child,neigh_edge;
  int chnum1,l_child,child_N_,edge1,neigh_N_,N_Neigh,N_,N_U,m;
  double *deriv, w,e1,e2,e3,e4,*coeff,strong_residual,alpha[3],beta[3],hK,hE,val[6],hE2;
  double *deriv1,*deriv2;
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
  TFESpace2D *fespace, *fespaceP;
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
  double *xyval_refNeigh1D[3*MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[3*MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[3*MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;
  double *xderiv_Cell1D, *yderiv_Cell1D, *xyval_Cell1D;
  int part,edge2neigh;
  int ee_verbose = 1,check_cont_u,check_cont_p,conform_grid=TDatabase::ParamDB->GRID_TYPE;
  double delta= TDatabase::ParamDB->FILTER_WIDTH_CONSTANT;
  TCollection *Coll;

  Coll = fespaces[0]->GetCollection();
  
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
  if (FELocalEstimator==ns_gradient_indicator)      // gradient indicator
    {
      for(i=0;i<N_Points;i++)                // for all quadrature points
        {
          w = weights[i]*AbsDetjk[i];

          deriv = Derivatives[i];            // all u derivatives in quadrature points 
          e1 = deriv[0];                     // x derivative
          e2 = deriv[1];                     // y derivative
          deriv = Derivatives[i+MaxN_QuadPoints_2D];// all v derivatives in quadrature points 
          e3 = deriv[0];
          e4 = deriv[1];
          estimated_error[0] += w*(e1*e1+e2*e2+e3*e3+e4*e4);
        } // endfor i
    }
  
/*************************************************************************/
/*                                                                       */
/*   RESIDUAL BASED EXPLICIT ERROR ESTIMATORS                            */
/*                                                                       */
/*************************************************************************/
  if ((FELocalEstimator>0)||(ErrorControl))
  {
    fespace = fespaces[0];
    fespaceP = fespaces[1];
    N_U =  U->GetLength();
    xi1DNeigh = new double[N_Points1D];
    eta1DNeigh = new double[N_Points1D];
    X1DNeigh = new double[N_Points1D];
    Y1DNeigh = new double[N_Points1D];
    X1DCell = new double[N_Points1D];
    Y1DCell = new double[N_Points1D];
    FEFunctValuesNeigh = new double[3*MaxN_BaseFunctions2D];
    for (i=0;i<N_Points1D;i++)
    {
      xyval_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
      xderiv_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
      yderiv_refNeigh1D[i]= new double[MaxN_BaseFunctions2D];
    }
    xderiv_Neigh1D = new double[3*N_Points1D];
    yderiv_Neigh1D = new double[3*N_Points1D];
    xyval_Neigh1D = new double[3*N_Points1D];
    xderiv_Cell1D = new double[3*N_Points1D];
    yderiv_Cell1D = new double[3*N_Points1D];
    xyval_Cell1D = new double[3*N_Points1D];

    check_cont_u=check_cont_p=0;
    if (TDatabase::ParamDB->ANSATZ_ORDER>0)
      check_cont_u=check_cont_p=1;
    if (TDatabase::ParamDB->ANSATZ_ORDER<-1)
      check_cont_u=1;
    // this has to be set for all continuous pressure spaces
    if ((TDatabase::ParamDB->VELOCITY_SPACE>0)&&(TDatabase::ParamDB->VELOCITY_SPACE<10))
      check_cont_p = 1;
    else
      check_cont_p = 0;
	
      

/*************************************************************************/
/*  strong residual                                                      */
/*************************************************************************/
    hK = cell->GetDiameter();
    strong_residual = 0;
    for(i=0;i<N_Points;i++)                // for all quadrature points
    {
      coeff = coeffs[i];
      w = weights[i]*AbsDetjk[i];
      //cout << " eps " << coeff[0] << " f0 " << coeff[1] << " f1 " << coeff[2] << endl;
      
      deriv = Derivatives[i];            // all u derivatives in quadrature points 
      deriv1 = Derivatives[i+MaxN_QuadPoints_2D];
      deriv2 = Derivatives[i+2*MaxN_QuadPoints_2D];
      w = weights[i]*AbsDetjk[i];
      if (NavierStokes)
        {
          e1 = -coeff[0]*(deriv[3]+deriv[4])+deriv[2]*deriv[0]+deriv1[2]*deriv[1]+deriv2[0]
            -coeff[1];                   // strong residual, u - comp
          e2 = -coeff[0]*(deriv1[3]+deriv1[4])+deriv[2]*deriv1[0]+deriv1[2]*deriv1[1]+deriv2[1]
            -coeff[2];                   // strong residual, v - comp
        }
      else
        {
          e1 = -coeff[0]*(deriv[3]+deriv[4])+deriv2[0]
            -coeff[1];                   // strong residual, u - comp
          e2 = -coeff[0]*(deriv1[3]+deriv1[4])+deriv2[1]
            -coeff[2];                   // strong residual, v - comp
        }
      strong_residual += w*(e1*e1+e2*e2);        // L^2 norm
    }                                     // endfor i
    alpha[0] = hK*hK;                      // weight for H^1 estimator
    alpha[1] = hK*hK*hK*hK;                // weight for L^2 estimator
    alpha[2] = 1;                          // weight for energy norm estimator
    if (TDatabase::ParamDB->P4==123456789)
      alpha[1]*=hK*hK/(delta*delta);
    if (hK*hK/coeff[0]<1)
      alpha[2] =hK*hK/coeff[0];           // update weight for energy norm estimator  
    for (i=1;i<4;i++)
       estimated_error[i] = alpha[i-1]*strong_residual;
    //estimated_error[i] = 0;
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
                case ROBIN:
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
                      m = i+N_Points1D;  
                      x0 =  X1D[j][i];
                      y0 =  Y1D[j][i];                      
                      boundedge->GetXYofT(t0,x0,y0);   // coordinates at quadrature points 
                      BoundaryValues[0](comp,t0,neumann_data);  // Neumann data
                      e1 = coeff[0] *(xderiv_1D[j][i]*nx+yderiv_1D[j][i]*ny) - neumann_data;
                      BoundaryValues[1](comp,t0,neumann_data);  // Neumann data
                      e2 = coeff[0] *(xderiv_1D[j][m]*nx+yderiv_1D[j][m]*ny) - neumann_data;                       
                      w = weights1D[i]*hE/2.0;
                      jump+= w*(e1*e1+e2*e2);          // integral on the edge
                    }
                  beta[0]=hE;                          // weight for H^1 estimator
                  beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                  if (TDatabase::ParamDB->P4==123456789)
                    beta[1]*=hE*hE/(delta*delta);
                  beta[2]=hE/coeff[0];                 // weight for energy norm estimator
                  if (1.0/sqrt(coeff[0])< beta[2])
                    beta[2]=1.0/sqrt(coeff[0]);
                  for (i=1;i<4;i++)
                    estimated_error[i] += beta[i-1]*jump;                  
                  break;
                 default:
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
                  if(ver1==ver2)                               // first part of long edge
                    {
                      part = -1;
                    }
                  else if(ver0==ver3)                          // second part of long edge
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
                  
                  BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbour    
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
                      FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] = Values[DOF[i]+N_U]; // v values                          
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
                      val[0]=val[1]= val[2]= val[3]=val[4]= val[5] = 0;
                       for(l=0;l<N_Neigh;l++)       // for all basis functions 
                        {
                          m = l+MaxN_BaseFunctions2D;
                          val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                          val[3] += FEFunctValuesNeigh[m] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[4] += FEFunctValuesNeigh[m] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                          val[5] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                          if (ee_verbose>1)
                            cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                              yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                        } // endfor l
                      m = i+N_Points1D; 
                      xderiv_Neigh1D[i]= val[0]; // for k-th 
                      yderiv_Neigh1D[i]= val[1]; // for k-th 
                      xyval_Neigh1D[i]= val[2]; // for k-th 
                      xderiv_Neigh1D[m]= val[3]; // for k-th 
                      yderiv_Neigh1D[m]= val[4]; // for k-th 
                      xyval_Neigh1D[m]= val[5]; // for k-th 
                    } // endfor i                                        
                  
                  
                  TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                              eta1DNeigh, 
                                              X1DNeigh, Y1DNeigh, absdet1D);   
                  // pressure second                   
                  // compute values at the quadrature points on the edge of 
                  // the neighbour element
                  
                  CurrEleNeigh = fespaceP->GetFE2D(neigh_N_,neigh);   // finite element on neighbour
                  eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                  
                  BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                  N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions
                  
                  bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                  bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour
                  
                  // compute gradients in reference cell of the neighbour
                  for (i=0;i<N_Points1D;i++)         // for all quadrature points
                    {                      
                      bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                    }
                  
                  RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                  TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);
                  
                  DOFP = GlobalNumbersP + BeginIndexP[neigh_N_];
                  
                  for(i=0;i<N_Neigh;i++)
                    {
                      FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] = ValuesP[DOFP[i]]; // p values
                      if (ee_verbose>1)
                        cout << " p value "  <<  FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] << endl;
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
                  
                  for(i=0;i<N_Points1D;i++)     // for all quadrature points on edge
                    {
                      val[0]= 0;
                      for(l=0;l<N_Neigh;l++)       // for all basis functions 
                        {
                          m = l+2*MaxN_BaseFunctions2D;
                          val[0] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                        } // endfor l
                      m = i+2*N_Points1D; 
                      xyval_Neigh1D[m]= val[0]; // for k-th 
                      //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
                    } // endfor i                                        
                      
                  jump=0.0;
                  absdetjk1D = hE/2;
                  for (i=0;i<N_Points1D;i++)           // compute jump
                    {
                      m = i+N_Points1D; 
                      l = m+N_Points1D; 
                      if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
                        cout << " wrong quad points_a " << X1D[j][i] << " , " << Y1D[j][i] 
                             << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                      if (check_cont_u)
                        {
                          if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
                            cout << " i " << i << " uval_a " << xyval_1D[j][i]<< " uneigh_a " << xyval_Neigh1D[i]<< endl;
                          if (fabs(xyval_Neigh1D[m]-xyval_1D[j][m])>1e-8)
                            cout << " i " << i << " vval_a " << xyval_1D[j][m]<< " vneigh_a " << xyval_Neigh1D[m]<< endl;
                        }
                      if (check_cont_p)
                        {
                          if (fabs(xyval_Neigh1D[l]-xyval_1D[j][l])>1e-8)
                            cout << " i " << i << " pval_a " << xyval_1D[j][l]<< " pneigh_a " << xyval_Neigh1D[l]<< endl;
                        }
                      e1 = coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
                                      +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny)
                        -(xyval_1D[j][l]-xyval_Neigh1D[l])*nx;
                      e2 = coeff[0] *((xderiv_1D[j][m]-xderiv_Neigh1D[m])*nx
                                      +(yderiv_1D[j][m]-yderiv_Neigh1D[m])*ny)
                        -(xyval_1D[j][l]-xyval_Neigh1D[l])*ny;
                      if (ee_verbose>1)
                        cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
                      w = weights1D[i]*absdetjk1D;
                      jump+= w*(e1*e1+e2*e2);                       // integral on the edge
                    }
                  if (ee_verbose>1)
                    cout << "jump " << jump << endl;
                  
                  beta[0]=hE;                          // weight for H^1 estimator
                  beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                  if (TDatabase::ParamDB->P4==123456789)
                    beta[1]*=hE*hE/(delta*delta);
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
                            cout << " something wrong 5 " <<  endl;
                          
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
                              FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] = Values[DOF[i]+N_U]; // v values                          
                              if (ee_verbose>1)
                                cout << " value " <<  FEFunctValuesNeigh[i] << 
                                  " " <<  FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] << endl;
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
                              val[0]=val[1]= val[2]= val[3] =val[4] = val[5] = 0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  m = l+MaxN_BaseFunctions2D;
                                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[3] += FEFunctValuesNeigh[m] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[4] += FEFunctValuesNeigh[m] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[5] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  if (ee_verbose>1)
                                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                                } // endfor l
                              m = i+N_Points1D; 
                              xderiv_Cell1D[i]= val[0]; // for k-th 
                              yderiv_Cell1D[i]= val[1]; // for k-th 
                              xyval_Cell1D[i]= val[2]; // for k-th 
                              xderiv_Cell1D[m]= val[3]; // for k-th 
                              yderiv_Cell1D[m]= val[4]; // for k-th 
                              xyval_Cell1D[m]= val[5]; // for k-th 
                            } // endfor i                                        

                          TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                      eta1DNeigh, 
                                                      X1DCell, Y1DCell, absdet1D);   

                          // pressure second                   
                          // compute values at the quadrature points on the edge of 
                          // the cell 
                      
                          CurrEleNeigh = fespaceP->GetFE2D(neigh_N_,cell);   // pressure space
                          eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                      
                          BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                          N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions

                          bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                          bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                          // compute gradients in reference cell of the neighbour
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {                      
                              bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                            }

                          RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                          TFEDatabase2D::SetCellForRefTrans(cell,RefTransNeigh);

                          DOFP = GlobalNumbersP + BeginIndexP[neigh_N_];

                          for(i=0;i<N_Neigh;i++)
                            {
                              FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] = ValuesP[DOFP[i]]; // p values
                              if (ee_verbose>1)
                                cout << " p value "  <<  FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] << endl;
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

                          for(i=0;i<N_Points1D;i++)     // for all quadrature points on edge
                            {
                              val[0]= 0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  m = l+2*MaxN_BaseFunctions2D;
                                  val[0] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                } // endfor l
                              m = i+2*N_Points1D; 
                              xyval_Cell1D[m]= val[0]; // for k-th 
                              //cout << "p value " << xyval_Cell1D[m] << endl;
                            } // endfor i                                        
                          
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
                              FEFunctValuesNeigh[i] = Values[DOF[i]];       // u values
                              FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] = Values[DOF[i]+N_U]; // v values                          
                              if (ee_verbose>1)
                                cout << " value " <<  FEFunctValuesNeigh[i] << 
                                  " " <<  FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] << endl;
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
                              val[0]=val[1]= val[2]=val[3]=val[4]=val[5]=0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  m = l+MaxN_BaseFunctions2D;
                                  val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[3] += FEFunctValuesNeigh[m] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[4] += FEFunctValuesNeigh[m] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                                  val[5] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                  if (ee_verbose>1)
                                    cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                      yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                                } // endfor l
                              m = i+N_Points1D; 
                              xderiv_Neigh1D[i]= val[0]; // for k-th 
                              yderiv_Neigh1D[i]= val[1]; // for k-th 
                              xyval_Neigh1D[i]= val[2]; // for k-th 
                              xderiv_Neigh1D[m]= val[3]; // for k-th 
                              yderiv_Neigh1D[m]= val[4]; // for k-th 
                              xyval_Neigh1D[m]= val[5]; // for k-th 
                            } // endfor i                                        
                          
                          TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                      eta1DNeigh, 
                                                      X1DNeigh, Y1DNeigh, absdet1D);   
                          
                          // pressure second                   
                          // compute values at the quadrature points on the edge of 
                          // the neighbour element
                          
                          CurrEleNeigh = fespaceP->GetFE2D(neigh_N_,child);   // finite element on neighbour
                          eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                      
                          BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                          N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions

                          bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                          bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                          // compute gradients in reference cell of the neighbour
                          for (i=0;i<N_Points1D;i++)         // for all quadrature points
                            {                      
                              bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                            }

                          RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                          TFEDatabase2D::SetCellForRefTrans(child,RefTransNeigh);

                          DOFP = GlobalNumbersP + BeginIndexP[neigh_N_];

                          for(i=0;i<N_Neigh;i++)
                            {
                              FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] = ValuesP[DOFP[i]]; // p values
                              if (ee_verbose>1)
                                cout << " p value "  <<  FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] << endl;
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
                          
                          for(i=0;i<N_Points1D;i++)     // for all quadrature points on edge
                            {
                              val[0]= 0;
                              for(l=0;l<N_Neigh;l++)       // for all basis functions 
                                {
                                  m = l+2*MaxN_BaseFunctions2D;
                                  val[0] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                                } // endfor l
                              m = i+2*N_Points1D; 
                              xyval_Neigh1D[m]= val[0]; // for k-th 
                              //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
                            } // endfor i
                                        
                          jump=0.0;
                          absdetjk1D = hE/(2.0*N_child);       // only half edge is considered
                          for (i=0;i<N_Points1D;i++)           // compute jump
                            {                              
                              m = i+N_Points1D; 
                              l = m+N_Points1D; 
                              if ((fabs(X1DCell[i]-X1DNeigh[i])+fabs(Y1DCell[i]-Y1DNeigh[i]))>1e-8)
                                 cout << " wrong quad points_c " << X1D[j][i] << " , " << Y1D[j][i] 
                                      << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                              if (check_cont_u)
                                {
                                  if (fabs(xyval_Neigh1D[i]-xyval_Cell1D[i])>1e-8)
                                    cout << " i " << i << " uval_c " << xyval_Cell1D[i]<< " uneigh_c " << xyval_Neigh1D[i]<< endl;
                                  if (fabs(xyval_Neigh1D[m]-xyval_Cell1D[m])>1e-8)
                                    cout << " i " << i << " vval_c " << xyval_Cell1D[m]<< " vneigh_c " << xyval_Neigh1D[m]<< endl;
                                }
                              if (check_cont_p)
                                {
                                  if (fabs(xyval_Neigh1D[l]-xyval_Cell1D[l])>1e-8)
                                    cout << " i " << i << " pval_c " << xyval_Cell1D[l]<< " pneigh_c " << xyval_Neigh1D[l]<< endl;
                                }
                              e1 = coeff[0] *((xderiv_Cell1D[i]-xderiv_Neigh1D[i])*nx
                                               +(yderiv_Cell1D[i]-yderiv_Neigh1D[i])*ny)
                                 -(xyval_Cell1D[l]-xyval_Neigh1D[l])*nx;
                               e2 = coeff[0] *((xderiv_Cell1D[m]-xderiv_Neigh1D[m])*nx
                                               +(yderiv_Cell1D[m]-yderiv_Neigh1D[m])*ny)
                                 -(xyval_Cell1D[l]-xyval_Neigh1D[l])*ny;
                            }
                          if (ee_verbose>1)
                            cout << "jump_c " << jump << endl;
                          hE2=hE/N_child;
                          beta[0]=hE2;                          // weight for H^1 estimator
                          beta[1]=hE2*hE2*hE2;                    // weight for L^2 estimator
                          if (TDatabase::ParamDB->P4==123456789)
                            beta[1]*=hE2*hE2/(delta*delta);
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
                      // find the local edge of neigh on which cell is -> l
                      neigh_edge=0;
                      while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;
                      refdesc=neigh->GetRefDesc();
                      refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVerNeigh); 
                      ver0=  cell->GetVertex(TmpEdVer[2*j]);
                      ver1=  cell->GetVertex(TmpEdVer[2*j+1]);
                      ver2 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge]);          // vertices of edge
                      ver3 = neigh->GetVertex(TmpEdVerNeigh[2*neigh_edge+1]);
                      if (((ver0==ver2)&&(ver1==ver3))||((ver0==ver3)&&(ver1==ver2)))
                        ;
                      else
                        {
                          cout << "wrong edge " << endl;                          
                        }
                      // velocity first                      
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
                          FEFunctValuesNeigh[i] = Values[DOF[i]];       // u values
                          FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] = Values[DOF[i]+N_U]; // v values                          
                          if (ee_verbose>1)
                            cout << " value " <<  FEFunctValuesNeigh[i] << 
                              " " <<  FEFunctValuesNeigh[i+MaxN_BaseFunctions2D] << endl;
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

                      for(i=0;i<N_Points1D;i++)     // for all quadrature points on edge
                        {
                          val[0]=val[1]= val[2]= val[3]=val[4]= val[5] = 0;
                          for(l=0;l<N_Neigh;l++)       // for all basis functions 
                            {
                              m = l+MaxN_BaseFunctions2D;
                              val[0] += FEFunctValuesNeigh[l] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[1] += FEFunctValuesNeigh[l] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[2] += FEFunctValuesNeigh[l] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                              val[3] += FEFunctValuesNeigh[m] * xderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[4] += FEFunctValuesNeigh[m] * yderiv_refNeigh1D[i][l]; // accumulate value of derivative
                              val[5] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                              if (ee_verbose>1)
                                cout << l << "  " << xderiv_refNeigh1D[i][l] << "  " << 
                                  yderiv_refNeigh1D[i][l] <<  "  " << FEFunctValuesNeigh[l] << endl;
                            } // endfor l
                          m = i+N_Points1D; 
                          xderiv_Neigh1D[i]= val[0]; // for k-th 
                          yderiv_Neigh1D[i]= val[1]; // for k-th 
                          xyval_Neigh1D[i]= val[2]; // for k-th 
                          xderiv_Neigh1D[m]= val[3]; // for k-th 
                          yderiv_Neigh1D[m]= val[4]; // for k-th 
                          xyval_Neigh1D[m]= val[5]; // for k-th 
                        } // endfor i                                        
                      
                      // just for testing quad points, may be deleted later
                      TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1DNeigh, 
                                                  eta1DNeigh, 
                                                  X1DNeigh, Y1DNeigh, absdet1D);
                      // pressure second                   
                      // compute values at the quadrature points on the edge of 
                      // the neighbour element
                      
                      CurrEleNeigh = fespaceP->GetFE2D(neigh_N_,neigh);   // finite element on neighbour
                      eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh); 
                      
                      BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbout    
                      N_Neigh = eleNeigh->GetN_DOF();                    // number of basis functions

                      bfNeigh = TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh);        
                      bf2DrefelementsNeigh = bfNeigh->GetRefElement();   // referenz cell of neighbour

                      // compute gradients in reference cell of the neighbour
                      for (i=0;i<N_Points1D;i++)         // for all quadrature points
                        {                      
                          bfNeigh->GetDerivatives(D00,xi1DNeigh[i],eta1DNeigh[i],xietaval_refNeigh1D[BaseFunctNeigh][i]);
                        }

                      RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
                      TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);

                      DOFP = GlobalNumbersP + BeginIndexP[neigh_N_];

                      for(i=0;i<N_Neigh;i++)
                        {
                          FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] = ValuesP[DOFP[i]]; // p values
                          if (ee_verbose>1)
                            cout << " p value "  <<  FEFunctValuesNeigh[i+2*MaxN_BaseFunctions2D] << endl;
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

                      for(i=0;i<N_Points1D;i++)     // for all quadrature points on edge
                        {
                          val[0]= 0;
                          for(l=0;l<N_Neigh;l++)       // for all basis functions 
                            {
                              m = l+2*MaxN_BaseFunctions2D;
                              val[0] += FEFunctValuesNeigh[m] * xyval_refNeigh1D[i][l]; // accumulate value of derivative
                            } // endfor l
                          m = i+2*N_Points1D; 
                          xyval_Neigh1D[m]= val[0]; // for k-th 
                          //cout << "p value " << xyval_Neigh1D[m] << " " << xyval_1D[j][m] << endl;
                        } // endfor i                                        
                      
                      jump=0.0;
                      absdetjk1D = hE/2;
                      for (i=0;i<N_Points1D;i++)           // compute jump
                        {
                          m = i+N_Points1D; 
                          l = m+N_Points1D;
                          
                          if ((fabs(X1D[j][i]-X1DNeigh[i])+fabs(Y1D[j][i]-Y1DNeigh[i]))>1e-8)
                            cout << " wrong quad points_b " << X1D[j][i] << " , " << Y1D[j][i] 
                                 << "   " << X1DNeigh[i] << " , " << Y1DNeigh[i]  << endl;
                          if (check_cont_u)
                            {
                              if (fabs(xyval_Neigh1D[i]-xyval_1D[j][i])>1e-8)
                                cout << " i " << i << " uval_b " << xyval_1D[j][i]<< " uneigh_b " << xyval_Neigh1D[i]<< endl;
                              if (fabs(xyval_Neigh1D[m]-xyval_1D[j][m])>1e-8)
                                cout << " i " << i << " vval_b " << xyval_1D[j][m]<< " vneigh_b " << xyval_Neigh1D[m]<< endl;
                            }
                          if (check_cont_p)
                            {
                              if (fabs(xyval_Neigh1D[l]-xyval_1D[j][l])>1e-8)
                                cout << " i " << i << " pval_b " << xyval_1D[j][l]<< " pneigh_b " << xyval_Neigh1D[l]<< endl;
                            }
                          e1 = coeff[0] *((xderiv_1D[j][i]-xderiv_Neigh1D[i])*nx
                                          +(yderiv_1D[j][i]-yderiv_Neigh1D[i])*ny)
                            -(xyval_1D[j][l]-xyval_Neigh1D[l])*nx;
                          e2 = coeff[0] *((xderiv_1D[j][m]-xderiv_Neigh1D[m])*nx
                                          +(yderiv_1D[j][m]-yderiv_Neigh1D[m])*ny)
                            -(xyval_1D[j][l]-xyval_Neigh1D[l])*ny;
                          if (ee_verbose>1)
                            cout << i<< " jumpx " << xderiv_1D[j][i] << " " << xderiv_Neigh1D[i] << endl;
                          w = weights1D[i]*absdetjk1D;
                          jump+= w*(e1*e1+e2*e2);                       // integral on the edge
                        }
                      if (ee_verbose>1)
                        cout << "jump " << jump << endl;
                      beta[0]=hE;                          // weight for H^1 estimator
                      beta[1]=hE*hE*hE;                    // weight for L^2 estimator
                      if (TDatabase::ParamDB->P4==123456789)
                        beta[1]*=hE*hE/(delta*delta);
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
  for (i=0;i<4;i++)
    estimated_local_error[i]=estimated_error[i];
}

#endif // #ifdef __2D__
