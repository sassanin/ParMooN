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
// @(#)Assemble2D.C        06.07 04/13/00
//
// Purpose:     bilinear form (discretized and edge-stabilized assemble)
//
// Author:      Rudolf Umla (07.07.2009)
//
// History:     
//
// =======================================================================

#include <DefineParams.h>

//#include <Assemble2D_edge_convdiv.h>
#include <Enumerations.h>
#include <Matrix2D.h>
#include <AuxParam2D.h>
#include <DiscreteForm2D.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>
#include <FEDatabase2D.h>
#include <NodalFunctional2D.h>
#include <SquareMatrix2D.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <Convolution.h>

#include <string.h>
#include <stdlib.h>


/******************************************************************************/
//
// compute edge terms for Oseen equations
//
/******************************************************************************/
#ifdef __2D__
void Assemble2D_edge_Oseen(CoeffFct2D *Coeff,int n_fespaces, TFESpace2D **fespaces,
			   int n_sqmatrices, TSquareMatrix2D **sqmatrices,
			   int n_matrices, TMatrix2D **matrices,
			   int n_rhs, double **rhs, TFESpace2D **ferhs,
			   BoundCondFunct2D **BoundaryConditions,
			   BoundValueFunct2D **BoundaryValues,
			   TAuxParam2D *Parameters)
{
  const int MaxN_BaseFunctions2D_Ersatz =100;

  double hk,nu,w,integrant,tau_par,tau_par2,tau_par3;
  int N_AllMatrices = n_sqmatrices+n_matrices,verbose;
  int i,j,k,l,l1,l2,l3,n,n_neigh,m,r,q,dummy,N_UsedElements,N_LocalUsedElements,ii,jj,ll,l_test,l_ansatz;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_Hanging,N_Rows;
  int N_Test, N_Ansatz, N_Joints, N_TestNeighbour, N_AnsatzNeighbour;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TBaseFunct2D *bf;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement, TestElementNeighbour, AnsatzElementNeighbour;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunctCell, TestFunctCell, AnsatzFunctCell, TestFunctNeigh, AnsatzFunctNeighbour;
  TCollection *Coll;
  TBaseCell *cell;
  TJoint *joint;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  int **GlobalNumbers, **BeginIndex;
  int **RhsGlobalNumbers, **RhsBeginIndex;
  int **TestGlobalNumbers, **TestBeginIndex;
  int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
  TFE2D *ele;
  TFEDesc2D *FEDesc_Obj;
  BF2DRefElements bf2Drefelements;
  double *weights, *xi, *eta, *weights1D, *weights_neigh, *xi_neigh, *eta_neigh, *weights1D_neigh;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], X_neigh[MaxN_QuadPoints_2D], Y_neigh[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D], AbsDetjk_neigh[MaxN_QuadPoints_2D],*AbsDetjk1D[4];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux, *aux2, *aux3, *aux4 ;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *Coeffs[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF, *TestDOF_neigh, *AnsatzDOF_neigh;
  double **SqEntries,**Entries;
  int **SqColInd, **ColInd, **SqRowPtr, **RowPtr;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s,integral;
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D_Ersatz];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double delta;
  bool *SecondDer;
  double max_b,dummy3,Re,Re_K;
  double max_nb;
  double *Coefficients1D[MaxN_QuadPoints_2D];
  double *Parameters1D[MaxN_QuadPoints_2D];

  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  //double xietaval_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double xideriv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  //double etaderiv_ref1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];

  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *xyval_ref1D_test[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D_test[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D_test[4][MaxN_QuadPoints_1D];
  double *xyval_ref1D_ansatz[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D_ansatz[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D_ansatz[4][MaxN_QuadPoints_1D];

  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  RefTrans2D RefTrans;
  int N_DOF;
  double *Values;
  //double value_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],value_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double xderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],xderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  //double yderiv_basefunct_ref1D[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz],yderiv_basefunct_ori[N_BaseFuncts2D][6][MaxN_BaseFunctions2D_Ersatz];
  double*** value_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** xderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** yderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double *value_basefunct_ori[6];
  double *xderiv_basefunct_ori[6];
  double *yderiv_basefunct_ori[6];
  double *x_pos_ref= new double[6];
  double *y_pos_ref=new double[6];
  double *x_pos=new double[6];
  double *y_pos=new double[6];
  double *value_basefunct_ori_neigh[6];
  double *xderiv_basefunct_ori_neigh[6];
  double *yderiv_basefunct_ori_neigh[6];
  double *x_pos_neigh=new double[6];
  double *y_pos_neigh=new double[6];
  double *dummy2=new double[6];

  int ref_n;
  int neigh_edge;
  int neigh_N_,N_Neigh;
  double absdet1D_neigh[MaxN_QuadPoints_2D];
  double xi1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D], eta1DNeigh[N_BaseFuncts2D][MaxN_QuadPoints_1D];
  double *X1DNeigh,*Y1DNeigh;
  TBaseCell *neigh;
  FE2D LocalUsedElements_neigh[N_FEs2D], CurrEleNeigh;
  BaseFunct2D BaseFunctNeigh;
  QuadFormula2D QuadFormulaNeigh;
  TQuadFormula2D *qfNeigh;
  QuadFormula1D LineQuadFormulaNeigh;
  TQuadFormula1D *qf1DNeigh;
  int LocN_BF_neigh[N_BaseFuncts2D];
  BaseFunct2D LocBF_neigh[N_BaseFuncts2D];
  int N_Points1DNeigh,N_PointsNeigh;
  double *weights1DNeigh,*zetaNeigh,*weightsNeigh,*xiNeigh,*etaNeigh;
  TFE2D *eleNeigh;
  RefTrans2D RefTransNeigh;
  BF2DRefElements bf2DrefelementsNeigh;
  int *DOF_neigh;
  double xietaval_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double xideriv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];
  double etaderiv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D_Ersatz];

  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xyval_refNeigh1D_test[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D_test[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D_test[MaxN_QuadPoints_1D];
  double *xyval_refNeigh1D_ansatz[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D_ansatz[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D_ansatz[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xyval_test[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv_test[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv_test[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xyval_ansatz[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv_ansatz[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv_ansatz[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

  OutPut("egde integrals "<<endl);

  //Ouput instruction
  verbose=TDatabase::ParamDB->SC_VERBOSE_AMG;

#ifdef __3D__
  double z0, z1;
#endif

  OutPut("egde integrals "<<endl);

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // prepare objects which are defined by the square matrices
  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];

    SqEntries = new double* [n_sqmatrices];
    SqColInd = new int* [n_sqmatrices];
    SqRowPtr =  new int* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
      j = sqmatrices[i]->GetHangingN_Entries();
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    }                                             // endfor
  }                                               // endif n_sqmatrices

  // print square matrices 
  if(verbose==2)
  {
    OutPut("Squarematrices before Edge assembling" << endl);
    for(k=0;k<n_sqmatrices;k++)
    {
      OutPut(endl);
      OutPut("sqmatrix: " << k << endl);
      SqRowPtr[k] = sqmatrices[k]->GetRowPtr();
      SqEntries[k] = sqmatrices[k]->GetEntries();
      SqColInd[k] = sqmatrices[k]->GetKCol();
      N_Rows = sqmatrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=SqRowPtr[k][i+1];
        for(j=SqRowPtr[k][i];j<end;j++)
        {
          OutPut("Matrix: " << setw(5) << i << setw(5) << SqColInd[k][j] << "   ");
          OutPut(setw(10) << SqEntries[k][j] << endl);
        }
      }
      OutPut(endl);
    }                         
  }

  // prepare objects which are defined by the rectangular matrices  
  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
    for(i=0;i<n_matrices;i++)
    {
      fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetTestSpace();
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();

      fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_matrices

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D_Ersatz*MaxN_BaseFunctions2D_Ersatz];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D_Ersatz];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D_Ersatz;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D_Ersatz;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D_Ersatz;
  }                                               // endif N_AllMatrices

  SecondDer = new bool[n_fespaces];
  SecondDer[0] = FALSE;

  // allocate arrays for storing information for the integration on the edges
  // for all basis functions 
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    value_basefunct_ref1D[i] = new double* [6];
    xderiv_basefunct_ref1D[i] = new double* [6];
    yderiv_basefunct_ref1D[i] = new double* [6];
    for (j=0;j<6;j++)
    {
      value_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      xderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];
      yderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D_Ersatz];

      memset( value_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( xderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      memset( yderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
    }
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
    xietaval_ref1D[i] = new double** [4];
    xideriv_ref1D[i] = new double** [4];
    etaderiv_ref1D[i] = new double** [4];
    for (j=0;j<4;j++)
    {
      xietaval_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      xideriv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      etaderiv_ref1D[i][j] = new double* [MaxN_QuadPoints_1D];
      for (n=0;n<MaxN_QuadPoints_1D;n++)
      {
        xietaval_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        xideriv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];
        etaderiv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D_Ersatz];

        memset( xietaval_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( xideriv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
        memset( etaderiv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D_Ersatz );
      }
    }
  }

  memset(Used, 0, N_FEs2D*SizeOfInt);

  for(i=0;i<n_fespaces;i++)
  {
    fespace = fespaces[i];                        /* fe space */
    n = fespace->GetN_UsedElements();             /* # used finite elements */
    UsedElements = fespace->GetUsedElements();    /* used finite elements */
    for(j=0;j<n;j++)                              /* for all finite elements */
    {
      CurrentElement = UsedElements[j];
      Used[CurrentElement] = 1;
    }                                             // enfor j
  }                                               // endfor i

  N_UsedElements = 0;                             /* compute number of used elements */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i]) N_UsedElements++;

  UsedElements = new FE2D[N_UsedElements];        /* store used finite elements */
  j=0;                                            /* in array */
  for(i=0;i<N_FEs2D;i++)
    if(Used[i])
  {
    UsedElements[j] = (FE2D)i;
    j++;
  }                                               // endif

  // ########################################################################
  // calculate values of base functions and derivatives on ref element
  // ########################################################################
  if(verbose == 1)
  {
    OutPut("# square matrices: " << n_sqmatrices << endl);
    OutPut(" # fespaces: " << n_fespaces<<endl);
    OutPut("N_UsedElements:" << N_UsedElements << " " << endl); fflush(0);
    OutPut("N_BaseFuncts2D:" << N_BaseFuncts2D << " " << endl);
    OutPut("MaxN_QuadPoints_1D:" << MaxN_QuadPoints_1D << " " << endl);
    OutPut("MaxN_BaseFunctions2D_Ersatz:" << MaxN_BaseFunctions2D_Ersatz << " " << endl);
  }
  // prepare integration on the edges of the reference cell
  for(n=0;n<N_UsedElements;n++)                   // for used finite elements
  {
    CurrentElement = UsedElements[n];
    l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
    LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1D->GetFormulaData(N_Points1D, weights1D, zeta);
    BaseFunctCell = BaseFuncts[CurrentElement];
                                                  // get base functions
    bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctCell);
    bf2Drefelements = bf->GetRefElement();
    switch(bf2Drefelements)                       // compute coordinates of line quadrature
    {                                             // points in reference cell
      // quadrilateral cell
      case BFUnitSquare :                         // edge 0

        bf->GetDerivatives(D00, -1, 1, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D10, -1, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, -1, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = -1;
        y_pos_ref[0] = 1;
        bf->GetDerivatives(D00, 1, -1, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(D10, 1, -1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 1, -1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] =-1;
        bf->GetDerivatives(D00, 1, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(D10, 1, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 1, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 1;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(D00, -1, -1, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(D10, -1, -1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, -1, -1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = -1;
        y_pos_ref[3] = -1;

        ref_n=4;

        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][0][j] = zeta[j];
          eta1D[BaseFunctCell][0][j] = -1;
          bf->GetDerivatives(D00, zeta[j], -1, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, zeta[j], -1, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, zeta[j], -1, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = 1;
          eta1D[BaseFunctCell][1][j] = zeta[j];
          bf->GetDerivatives(D00, 1, zeta[j], xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, 1, zeta[j], xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, 1, zeta[j], etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = -zeta[j];
          eta1D[BaseFunctCell][2][j] = 1;
          bf->GetDerivatives(D00, -zeta[j], 1, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, -zeta[j], 1, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, -zeta[j], 1, etaderiv_ref1D[BaseFunctCell][2][j]);
        }                                         // edge 3
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][3][j] = -1;
          eta1D[BaseFunctCell][3][j] = -zeta[j];
          bf->GetDerivatives(D00, -1, -zeta[j], xietaval_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D10, -1, -zeta[j], xideriv_ref1D[BaseFunctCell][3][j]);
          bf->GetDerivatives(D01, -1, -zeta[j], etaderiv_ref1D[BaseFunctCell][3][j]);
        }
        break;

      case BFUnitTriangle :                       // triangular cell

        bf->GetDerivatives(D00, 0, 0, value_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D10, 0, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[0] = 0;
        y_pos_ref[0] = 0;
        bf->GetDerivatives(D00, 1, 0, value_basefunct_ref1D[BaseFunctCell][1]);
        bf->GetDerivatives(D10, 1, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 1, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[1] = 1;
        y_pos_ref[1] = 0;
        bf->GetDerivatives(D00, 0, 1, value_basefunct_ref1D[BaseFunctCell][2]);
        bf->GetDerivatives(D10, 0, 1, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 1, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[2] = 0;
        y_pos_ref[2] = 1;

        bf->GetDerivatives(D00, 0.5, 0, value_basefunct_ref1D[BaseFunctCell][3]);
        bf->GetDerivatives(D10, 0.5, 0, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[3] = 0.5;
        y_pos_ref[3] = 0;
        bf->GetDerivatives(D00, 0.5, 0.5, value_basefunct_ref1D[BaseFunctCell][4]);
        bf->GetDerivatives(D10, 0.5, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0.5, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[4] = 0.5;
        y_pos_ref[4] = 0.5;
        bf->GetDerivatives(D00, 0, 0.5, value_basefunct_ref1D[BaseFunctCell][5]);
        bf->GetDerivatives(D10, 0, 0.5, xderiv_basefunct_ref1D[BaseFunctCell][0]);
        bf->GetDerivatives(D01, 0, 0.5, yderiv_basefunct_ref1D[BaseFunctCell][0]);
        x_pos_ref[5] = 0;
        y_pos_ref[5] = 0.5;

        ref_n=6;

        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][0][j] = (zeta[j]+1)/2;
          eta1D[BaseFunctCell][0][j] = 0;
          bf->GetDerivatives(D00, (zeta[j]+1)/2, 0, xietaval_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D10, (zeta[j]+1)/2, 0, xideriv_ref1D[BaseFunctCell][0][j]);
          bf->GetDerivatives(D01, (zeta[j]+1)/2, 0, etaderiv_ref1D[BaseFunctCell][0][j]);
        }                                         // edge 1
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][1][j] = (-zeta[j]+1)/2;
          eta1D[BaseFunctCell][1][j] = (zeta[j]+1)/2;
          bf->GetDerivatives(D00, (-zeta[j]+1)/2, (zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D10, (-zeta[j]+1)/2, (zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][1][j]);
          bf->GetDerivatives(D01, (-zeta[j]+1)/2, (zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][1][j]);
        }                                         // edge 2
        for (j=0;j<N_Points1D;j++)                // for all quadrature points
        {
          xi1D[BaseFunctCell][2][j] = 0;
          eta1D[BaseFunctCell][2][j] = (-zeta[j] +1)/2;
          bf->GetDerivatives(D00, 0, (-zeta[j]+1)/2, xietaval_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D10, 0, (-zeta[j]+1)/2, xideriv_ref1D[BaseFunctCell][2][j]);
          bf->GetDerivatives(D01, 0, (-zeta[j]+1)/2, etaderiv_ref1D[BaseFunctCell][2][j]);
        }
        break;
    }
  }                                               // endfor n

  for(l=0;l<ref_n;l++)
  {
    value_basefunct_ori[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    value_basefunct_ori_neigh[l] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  for(m=0;m<4;m++)                                // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[m] = new double[N_Points1D];              // coordinates of edge i
    Y1D[m] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[m] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xyval_ref1D_test[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D_test[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D_test[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xyval_ref1D_ansatz[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      xderiv_ref1D_ansatz[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];
      yderiv_ref1D_ansatz[m][j] = new double[MaxN_BaseFunctions2D_Ersatz];

    }
  }                                               // endfor m

  for (j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xyval_refNeigh1D_test[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D_test[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D_test[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xyval_refNeigh1D_ansatz[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D_ansatz[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D_ansatz[j] = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################
  // get number of parameters of equation
  N_Parameters = Parameters->GetN_Parameters();   
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coeffs[j] = aux2 + j*20;

  aux3 = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Parameters1D[j] = aux3 + j*N_Parameters;

  aux4 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coefficients1D[j] = aux4 + j*20;

  // ########################################################################
  // prepare loop over cells
  // ########################################################################

  // all spaces use same Coll
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)                          // set clipboard of cells on finest
  {
    cell=Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  // ########################################################################
  // loop over all cells
  // ########################################################################
  for(i=0;i<N_Cells;i++)                          // for all cells on the finest level
  {
    cell = Coll->GetCell(i);                      // next cell

    for(n=0;n<n_sqmatrices;n++)
    {
      if ( verbose == 1)
      {
        OutPut("-------------------------------------------" << endl);
        OutPut(" Enter squarematrix " << n << endl);
        OutPut("-------------------------------------------" << endl);
      }

      if(n>4)
      {
	  OutPut("Error in Assemble2D_edge_Oseen: only sqmatrices 0,...,4 are implemented so far" 
		 << endl);  
	  exit(4711);
      }
      // calculate all needed derivatives of this FE function
      fespace = sqmatrices[n]->GetFESpace();
      CurrentElement = fespace->GetFE2D(i,cell);  // finite element on cell

      BaseFunctCell = BaseFuncts[CurrentElement]; // basis functions
      N_ = N_BaseFunct[CurrentElement];           // # basis functions
      DOF = GlobalNumbers[n] + BeginIndex[n][i];  // dof of current mesh cell

      LocalUsedElements[0] = CurrentElement;
      LocN_BF[0] = N_BaseFunct[CurrentElement];   // local basis functions
      LocBF[0] = BaseFuncts[CurrentElement];
      SecondDer[0] = FALSE;
      RefTrans = TFEDatabase2D::GetOrig(1, LocalUsedElements,
					Coll, cell, SecondDer,
					N_Points, xi, eta, weights, X, Y, AbsDetjk);
      if(N_Parameters>0)                          // get parameters of equ.
        Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

                                                  // get coefficients of pde
      if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);
      // start to compute L-infty norm of convection in given cell
      // Oseen: convection is given (= solution)
      max_b=0;
      for(j=0;j<N_Points;j++)
      {
        dummy3 = sqrt(Coeffs[j][3]*Coeffs[j][3]+Coeffs[j][4]*Coeffs[j][4]);
        if (verbose == 1) 
	    OutPut("X[" << j <<"], " <<  "Y[" << j <<"]: " << X[j] << " " 
		   << Y[j]<< " |b|= " << dummy3 << endl);
        if(dummy3 > max_b) max_b= dummy3;
      }

      // prepare 1D quadrature formula
      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
      if(verbose ==1){ OutPut("Polynomial degree on cell: " << l << endl);}
      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

      if(verbose ==1)
      {
        for(j=0;j<N_Points1D; j++)
        {
          OutPut("weights1D["<<j<<"]:" <<  weights1D[j] << endl);
        }
        OutPut(endl);
      }

                                                  // update data base
      TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
        ->MakeRefElementData(LineQuadFormula);
      N_Edges=cell->GetN_Edges();                 // # edges

      if(verbose ==1)
      {
        for(r=0;r<N_Edges;r++)
        {
          cell->GetVertex(r)->GetCoords(x0, y0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
          OutPut("Local edge r: " << r << " Ecke A " << x0 << " " << y0 << 
		 " Ecke B " << x1 << " " << y1 << endl);
        }
      }

      for(r=0;r<N_Edges;r++)                      // loop over all edges of cell
      {                                           // get original coordinates of edge quad. points
        TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][r],
				      eta1D[BaseFunctCell][r],
				      X1D[r], Y1D[r], AbsDetjk1D[r]);

        for(j=0;j<N_Points1D;j++)                 // get values and derivatives in original cell
        {
          TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctCell][r][j],
				       eta1D[BaseFunctCell][r][j],
				       TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
				       Coll, (TGridCell *)cell,
				       xietaval_ref1D[BaseFunctCell][r][j],
				       xideriv_ref1D[BaseFunctCell][r][j],
				       etaderiv_ref1D[BaseFunctCell][r][j],
				       xyval_ref1D[r][j],
				       xderiv_ref1D[r][j],
				       yderiv_ref1D[r][j]);
        }
      }                                           // endfor r

      TFEDatabase2D::GetOrigFromRef(RefTrans,ref_n,x_pos_ref,y_pos_ref,x_pos,y_pos,dummy2);
      for(l=0;l<ref_n;l++)
      {  
	  //OutPut("x_pos_ref[l]: " << x_pos_ref[l] << "value_basefunct_ref1D[BaseFunctCell][l]: " << value_basefunct_ref1D[BaseFunctCell][l]);
	TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
				     y_pos_ref[l],
				     TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
				     Coll, (TGridCell *)cell,
				     value_basefunct_ref1D[BaseFunctCell][l],
				     xderiv_basefunct_ref1D[BaseFunctCell][l],
				     yderiv_basefunct_ref1D[BaseFunctCell][l],
				     value_basefunct_ori[l],
				     xderiv_basefunct_ori[l],
				     yderiv_basefunct_ori[l]);
      }
      // start to compute L-infty norm of convection times normal
      // Oseen: convection is given (= solution)
      max_nb=0;
      for(r=0;r<N_Edges;r++)
      {
        if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
        cell->GetVertex(r)->GetCoords(x0, y0);
        cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
        // compute length of the boundary edge
        hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
        // compute normal vector to this boundary (normalized)
        nx = (y1-y0)/hE;
        ny = (x0-x1)/hE;

        for(j=0;j<N_Points1D; j++)
        {
          dummy3= fabs(Coefficients1D[j][3]*nx+Coefficients1D[j][4]*ny);
          if (verbose == 1) 
	      OutPut("X1D, Y1D on Edge "<< r << " in quadpoint " << j << 
		     ": " << X1D[r][j] << " " << Y1D[r][j]<< " |b*n|= " << dummy3 << endl);
          if(dummy3>max_nb) max_nb=dummy3;
        }

      }

      if (verbose == 1)OutPut(endl);
      // continue to compute L-infty norm of convection
      for(r=0;r<N_Edges;r++)
      {
        if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
        for(j=0;j<N_Points1D; j++)
        {
          dummy3= sqrt(Coeffs[j][3]*Coeffs[j][3]+Coeffs[j][4]*Coeffs[j][4]);
          if (verbose == 1) OutPut("X1D, Y1D on Edge "<< r << " in quadpoint " 
				   << j << ": " << X1D[r][j] << " " << Y1D[r][j]<< 
				   " |b|= " << dummy3 << endl);
          if(dummy3 > max_b) max_b= dummy3;
        }
      }

      for(r=0;r<N_Edges;r++)
      {

        if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
        neigh=cell->GetJoint(r)->GetNeighbour(cell);
        // If there is a neighbour to the edge, do...
        if(neigh)
        {                                         
          // Get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
          if(i<q||i>=q)
          {
            // calculate all needed derivatives of this FE function
                                                  // finite element on cell
            CurrentElement = fespace->GetFE2D(i,cell);
                                                  // finite element on neighbour
            CurrEleNeigh = fespace ->GetFE2D(q,neigh);
            BaseFunctNeigh = BaseFuncts[CurrEleNeigh];
            eleNeigh =  TFEDatabase2D::GetFE2D(CurrEleNeigh);
            //BaseFunctNeigh = eleNeigh->GetBaseFunct2D_ID();    // basis functions on neighbour
            N_Neigh = eleNeigh->GetN_DOF();       // number of basis functions on neighbour
                                                  // dof of current mesh cell on neighbour cell
            DOF_neigh = GlobalNumbers[n] + BeginIndex[n][q];

            LocalUsedElements_neigh[0] = CurrEleNeigh;
                                                  // local basis functions
            LocN_BF_neigh[0] = N_BaseFunct[CurrEleNeigh];
            LocBF_neigh[0] = BaseFuncts[CurrEleNeigh];

            RefTransNeigh = TFEDatabase2D::GetOrig(1, LocalUsedElements_neigh,
              Coll, neigh, SecondDer,
              N_Points, xi_neigh, eta_neigh, weights_neigh, X_neigh, Y_neigh, AbsDetjk_neigh);

            /* To do: Include this for FE-Spaces of which the polynomial degree
            of the basis fuctions depend on the cell
            RefTrans_neigh = TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements_neigh,
                                                    Coll, neigh, SecondDer,
                                                    N_Points_neigh, xi_neigh, eta_neigh, weights_neigh,
                                                X_neigh, Y_neigh,AbsDetjk_neigh);
            if(N_Parameters>0)                // get parameters of equ.
            Parameters->GetParameters(N_Points_neigh, Coll, neigh, q, xi_neigh, eta_neigh, X_neigh, Y_neigh, Param_neigh);

            if(Coeff)                               // get coefficients of pde in the neighbour cell
            Coeff(N_Points_neigh, X_neigh, Y_neigh, Param_neigh, );

            // prepare 1D quadrature formula in the neighbour cell
            l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrEleNeigh);
            LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
            qf1D_neigh = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
            qf1D_neigh->GetFormulaData(N_Points1D_neigh, weights1D_neigh, zeta_neigh)
            */
            // Get edge of the neighbour cell which is the edge r
            neigh_edge=0;
            while(neigh->GetJoint(neigh_edge)->GetNeighbour(neigh)!=cell) neigh_edge ++;

            //RefTransNeigh= eleNeigh->GetRefTransID();          // reftrafo of neighbour
            //TFEDatabase2D::SetCellForRefTrans(neigh,RefTransNeigh);
            //TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrEleNeigh)  // update data base
            //->MakeRefElementData(LineQuadFormula);

            for(m=0;m<4;m++)                      // arrays for coordinates on neighbour cell
            {
              X1D_neigh[m] = new double[N_Points1D];
              Y1D_neigh[m] = new double[N_Points1D];
            }

            // get original coordinates of edge quad. points of neighbour cell
            TFEDatabase2D::GetOrigFromRef(RefTransNeigh,N_Points1D, xi1D[BaseFunctNeigh][neigh_edge],
              eta1D[BaseFunctNeigh][neigh_edge],X1D_neigh[neigh_edge], Y1D_neigh[neigh_edge], AbsDetjk1D[neigh_edge]);

            // get values and derivatives on original neighbour cell on edge neigh_edge
            for (j=0;j<N_Points1D;j++)
            {
              // OutPut(" xi1D " << xi1D[BaseFunctNeigh][neigh_edge][j] << " eta1D " << eta1D[BaseFunctCell][r][j] << " BaseFunctNeigh: " << BaseFunctNeigh << "\n");
              if(verbose ==1){  OutPut("X1D[r][j]: " << X1D[r][j] << " Y1D[r][j]: " << Y1D[r][j] << " X1D[neigh_edge][j] " <<  X1D[neigh_edge][j] << "  Y1D[neigh_edge][j]: " <<  Y1D[neigh_edge][j] <<  endl);}
            }

            if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0] )
            {
              if(verbose ==1)
              {
                OutPut("Quadrature points on neighbour edge in the correct order." << endl);
              }
              for (j=0;j<N_Points1D;j++)
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
                  eta1D[BaseFunctNeigh][neigh_edge][j],
                  TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                  Coll, (TGridCell *)neigh,
                  xietaval_ref1D[BaseFunctNeigh][neigh_edge][j],
                  xideriv_ref1D[BaseFunctNeigh][neigh_edge][j],
                  etaderiv_ref1D[BaseFunctNeigh][neigh_edge][j],
                  xyval_refNeigh1D[j],
                  xderiv_refNeigh1D[j],
                  yderiv_refNeigh1D[j]);
              }                                   //endfor j
            }                                     //endif
            else
            {
              if(verbose ==1)
              {
                OutPut("Inverse the order of the quadrature oints on neighbour edge !" << endl);
              }
              for (j=0;j<N_Points1D;j++)
              {
                TFEDatabase2D::GetOrigValues(RefTransNeigh, xi1D[BaseFunctNeigh][neigh_edge][j],
                  eta1D[BaseFunctNeigh][neigh_edge][j],
                  TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                  Coll, (TGridCell *)neigh,
                  xietaval_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  xideriv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  etaderiv_ref1D[BaseFunctNeigh][neigh_edge][N_Points1D-j-1],
                  xyval_refNeigh1D[j],
                  xderiv_refNeigh1D[j],
                  yderiv_refNeigh1D[j]);

              }                                   //endfor j
            }                                     //endelse

            TFEDatabase2D::GetOrigFromRef(RefTransNeigh,ref_n,x_pos_ref,y_pos_ref,x_pos_neigh,y_pos_neigh,dummy2);
            for(l=0;l<ref_n;l++)
            {
              TFEDatabase2D::GetOrigValues(RefTrans, x_pos_ref[l],
                y_pos_ref[l],
                TFEDatabase2D::GetBaseFunct2D(BaseFunctNeigh),
                Coll, (TGridCell *)neigh,
                value_basefunct_ref1D[BaseFunctNeigh][l],
                xderiv_basefunct_ref1D[BaseFunctNeigh][l],
                yderiv_basefunct_ref1D[BaseFunctNeigh][l],
                value_basefunct_ori_neigh[l],
                xderiv_basefunct_ori_neigh[l],
                yderiv_basefunct_ori_neigh[l]);
            }

            for(k=0;k<N_; k++)
            {
              dummy = 0;
              l=0;
              // Check if basis function k of cell i is in the FE-Space of neighbour cell q
              while(l<N_Neigh && dummy == 0)
              {
                if(DOF[k] == DOF_neigh[l])dummy=1;
                l++;
              }
              l = l-1;
              // if basis function k of cell i is in the local FE-Space of neighbour cell q do
              if(dummy ==1 )
              {                                   // Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
                for(j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][k] = xyval_ref1D[r][j][k]  -  xyval_refNeigh1D[j][l];
                  jump_xderiv[j][k] = xderiv_ref1D[r][j][k] - xderiv_refNeigh1D[j][l];
                  jump_yderiv[j][k] = yderiv_ref1D[r][j][k] - yderiv_refNeigh1D[j][l];

                  // OutPut(" k: " << k << " l: " << l << " DOF[k]: " << DOF[k] << " DOF_neigh[l]: " << DOF_neigh[l] << endl);
                  if(verbose ==1)
                  {
                    OutPut("xyval_Zelle: "<< xyval_ref1D[r][j][k] <<" xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("xderiv_Zelle: "<< xderiv_ref1D[r][j][k] <<  " xderiv_Neigh: " << xderiv_refNeigh1D[j][l] << " jump= " << jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("yderiv_Zelle: "<< yderiv_ref1D[r][j][k] <<  " yderiv_Neigh: " << yderiv_refNeigh1D[j][l] << " jump= " << jump_yderiv[j][k] << " of basefunction: "<< DOF[k]  << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut(endl);
                  }
                  if(fabs(jump_xyval[j][k])>=0.0000000000001){OutPut(" Error in Assemble Oseen: jump over edge although function is continuous "<< endl); exit(4711); }
                }
              }                                   //endif
              // if basis function k of cell i is NOT in the local FE-Space of neighbour cell q do
              if (dummy == 0)
              {
                for(j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][k]  = xyval_ref1D[r][j][k] ;
                  jump_xderiv[j][k] = xderiv_ref1D[r][j][k];
                  jump_yderiv[j][k] = yderiv_ref1D[r][j][k];

                  if(verbose ==1)
                  {
                    OutPut(" No Neighbour: xyval_Zelle: "<< xyval_ref1D[r][j][k] << " jump= " <<  jump_xyval[j][k] <<  " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("No Neighbour: x-deriv-jump= " <<  jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("No Neighbour: y-deriv-jump= " <<  jump_yderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
                    OutPut(endl);
                  }
                  if(fabs(jump_xyval[j][k])>=0.0000000000001){OutPut(" Error in Assemble Oseen: jump over edge although function is continuous "<< endl); exit(4711); }
                }                                 //endfor j
              }                                   //endif
            }                                     //endfor k

            // Then for the basis functions of neighbour cell q
            //which are not in the local FE-Space of cell i
            for(l=0;l<N_Neigh; l++)
            {
              dummy = 0;
              k=0;
              while(k<N_ && dummy == 0 )
              {
                if(DOF_neigh[l] == DOF[k]) dummy=1 ;
                k++;
              }
              k=k-1;
              // If basis function l of neighbour cell q is NOT  in the local FE-Space of cell i do
              if( dummy == 0)
              {

                for(j=0;j<N_Points1D;j++)
                {
                  jump_xyval[j][l+N_] = -xyval_refNeigh1D[j][l] ;
                  jump_xderiv[j][l+N_]= -xderiv_refNeigh1D[j][l];
                  jump_yderiv[j][l+N_]= -yderiv_refNeigh1D[j][l];

                  if(verbose ==1)
                  {
                    OutPut("Neighbour!!" << "xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][l+N_]<<  " of basefunction: "<< DOF_neigh[l] << " in cell: " << q << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("Neighbour!! " << "x-deriv-jump: " << jump_xderiv[j][l+N_] << " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("Neighbour!! y-deriv-jump= " <<  jump_yderiv[j][l+N_]<< " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
                    OutPut(endl);
                  }
                  if(fabs(jump_xyval[j][k])>=0.0000000000001){OutPut(" Error in Assemble Oseen: jump over edge although function is continuous "<< endl); exit(4711); }
                }                                 //endfor j
              }                                   //endif
            }                                     //endfor l

            // #################################################################################
            // Compute the edge integrals with the jumps of the basis functions and their derivatives
            // #################################################################################

            // get vertices of boundary edge
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            hk = cell->GetDiameter();
            if(verbose ==1){  OutPut("Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
            // compute length of the boundary edge
            hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            if(verbose ==1) OutPut("nx: " << nx << " ny: "<< ny<<endl);
            // tangential normal vector to this boundary (normalized)
            tx = (x1-x0)/hE;
            ty = (y1-y0)/hE;
            tau_par = TDatabase::ParamDB->TAU;
            tau_par2 = TDatabase::ParamDB->TAU2;
            tau_par3 = TDatabase::ParamDB->TAU3;
            Re =  TDatabase::ParamDB->RE_NR;	    
            Re *= max_b * hk;
            if (Re > 1) Re=1;

            /*  max_nb=0;
              for(j=0;j<N_Points1D; j++)
              {dummy3= fabs(Coefficients1D[j][3]*nx+Coefficients1D[j][4]*ny);
                if(dummy3>max_nb) max_nb=dummy3;
              }*/
            // max_b=1;
            //  max_nb=1;
            if (verbose == 1)
            {
              OutPut("tau edge stabilization: " << tau_par << endl);
              OutPut("tau2 edge stabilization: " << tau_par2 << endl);
              OutPut("tau3 edge stabilization: " << tau_par3 << endl);
              OutPut("Supremumsnorm of b in cell i: " << max_b << endl);
              OutPut("Supremumsnorm of b*n on edge r : " << max_nb << endl);
              OutPut("Reynoldsnumber in cell i : " << Re << endl);
              OutPut("Diameter of cell i:" << hk << endl);
            }

            fespace = sqmatrices[n]->GetFESpace();
            ActiveBound = fespace->GetActiveBound();
            SqEntries[n] = sqmatrices[n]->GetEntries();
            SqRowPtr[n] = sqmatrices[n]->GetRowPtr();
            SqColInd[n] = sqmatrices[n]->GetKCol();

            if(verbose ==1)
            {
              OutPut("ActiveBound von sqmatrix " << n << "  :" << ActiveBound << endl);
              for(j=0;j<N_Points1D; j++)
              {
                OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<< Coefficients1D[j][3] <<" b2: " << Coefficients1D[j][4] << endl);
              }
            }
            //edge integrals: test function from cell <-> ansatz function from cell
            if(verbose ==1){  OutPut("testfct. cell<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              //OutPut("dof_ii" << dof_ii << endl);
              // Dirichlet node
              if (dof_ii>=ActiveBound) continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch(n)
                  {
                    case 0:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*nx*
			  (jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj]*jump_xderiv[j][ii]);
                      break;

                    case 1:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*
			  (jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj]*jump_xderiv[j][ii]);
                      break;

                    case 2:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj]*jump_yderiv[j][ii]);
                      break;

                    case 3:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*ny*ny*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj]*jump_yderiv[j][ii]);
                      break;

                    case 4:
                      //pressure stabilization term
                      //  integrant = tau_par3*Re*hk*hk*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      if(max_b>1e-10) 
		      { if(1/max_b < hk*TDatabase::ParamDB->RE_NR)
			  integrant=tau_par3/max_b*hk*hk*(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			      *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);

                        else  integrant=tau_par3*TDatabase::ParamDB->RE_NR*hk*hk*hk*(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			      *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);


                      }
                      else
                        integrant = tau_par3*hk*hk*hk*TDatabase::ParamDB->RE_NR
			    *(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			    *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);
                      break;
                  }

                  w = weights1D[j]*hE/2;//AbsDetjk1D[r][j];
                  //w = weights1D[j]*AbsDetjk1D[r][j];
                  integral += w*integrant;        // integral on the edge
                }

                // update matrix
                found = 0;
                for (ll=SqRowPtr[n][dof_ii];ll < SqRowPtr[n][dof_ii+1]; ll++)
                {
                  if (SqColInd[n][ll] == dof_jj)
                  {
                    if(verbose ==1)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    SqEntries[n][ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  OutPut("ERROR in Assemble_edge_integrals (test function cell - ansatz function cell) " << endl);
                  exit(4711);
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)
            if (verbose ==1){  OutPut(endl);}

            //edge integrals: test function from cell <-> ansatz function from neighbour
            if(verbose ==1){  OutPut("testfct. cell<-> ansatzfct. neighbour Integral:" << endl);}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound)
                continue;
              for (jj=0;jj<N_Neigh;jj++)
              {
                dof_jj = DOF_neigh[jj];

                //################################################################################################
                //#### Check if ansatz  function jj  of the neighbour cell is in the FE-Space of the cell  ####
                dummy = 0;
                l=0;
                while(l<N_ && dummy == 0)
                {
                  if(dof_jj == DOF[l])dummy=1;
                  l++;
                }
                if (dummy == 1) continue;
                //################################################################################################
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch(n)
                  {
                    case 0:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*nx*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii]);
                      break;

                    case 1:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj+N_]*jump_xderiv[j][ii]);
                      break;

                    case 2:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj+N_]*jump_yderiv[j][ii]);
                      break;

                    case 3:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*ny*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj+N_]*jump_yderiv[j][ii]);
                      break;

                   case 4:
                      //pressure stabilization term
                      //  integrant = tau_par3*Re*hk*hk*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      if(max_b>1e-10) 
		      { if(1/max_b < hk*TDatabase::ParamDB->RE_NR)
			  integrant=tau_par3/max_b*hk*hk*(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			      *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);

                        else  integrant=tau_par3*hk*TDatabase::ParamDB->RE_NR*hk*hk*(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			      *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);


                      }
                      else
                        integrant = tau_par3*hk*hk*hk*TDatabase::ParamDB->RE_NR
			    *(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			    *(nx*jump_xderiv[j][ii]+ny*jump_yderiv[j][ii]);
                      break;
                      // integral on the edge
                  }
                  w = weights1D[j]* hE/2;//AbsDetjk1D[r][j];
                  //w = weights1D[j]* AbsDetjk1D[r][j];
                  integral+= w*integrant;
                }
                // update first matrix
                found = 0;
                for (ll=SqRowPtr[n][dof_ii];ll < SqRowPtr[n][dof_ii+1]; ll++)
                {
                  if (SqColInd[n][ll] == dof_jj)
                  {
                    if(verbose ==1)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    SqEntries[n][ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  //OutPut("ERROR in Assemble_edge_integrals (test function cell - ansatz function neighbour)" << endl);
                  exit(4711);
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)
            if (verbose ==1){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from cell
            if(verbose ==1){ OutPut("testfct. neighbour<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound)
                continue;

              //################################################################################################
              //#### Check if test function ii  of the neighbour cell is in the FE-Space of the cell    ####
              dummy = 0;
              l=0;

              while(l<N_ && dummy == 0)
              {
                if(dof_ii == DOF[l])dummy=1;
                l++;
              }
              if (dummy == 1) continue;
              //################################################################################################

              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];

                if(dummy == 0)                    //
                  // OutPut("jj " << dof_jj << endl);
                  // initialize the boundary integrals
                  integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch(n)
                  {
                    case 0:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*nx*(jump_xderiv[j][jj]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj]*jump_xderiv[j][ii+N_]);
                      break;

                    case 1:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj]*jump_xderiv[j][ii+N_]);
                      break;

                    case 2:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj]*jump_yderiv[j][ii+N_]);
                      break;

                    case 3:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*ny*ny*(jump_xderiv[j][jj]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj]*jump_yderiv[j][ii+N_]);
                      break;

                   case 4:
                      //pressure stabilization term
                      //  integrant = tau_par3*Re*hk*hk*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      if(max_b>1e-10) 
		      { if(1/max_b < hk*TDatabase::ParamDB->RE_NR)
			  integrant=tau_par3/max_b*hk*hk*(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			      *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);

                        else  integrant=tau_par3*hk*TDatabase::ParamDB->RE_NR*hk*hk*(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			      *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);


                      }
                      else
                        integrant = tau_par3*hk*hk*hk*TDatabase::ParamDB->RE_NR
			    *(nx*jump_xderiv[j][jj]+ny*jump_yderiv[j][jj])
			    *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);
                      break;
                  }

                  w = weights1D[j]* hE/2;//AbsDetjk1D[r][j];
                  //w = weights1D[j]* AbsDetjk1D[r][j];
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=SqRowPtr[n][dof_ii];ll < SqRowPtr[n][dof_ii+1]; ll++)
                {
                  if (SqColInd[n][ll] == dof_jj)
                  {
                    if(verbose ==1)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    SqEntries[n][ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  OutPut("ERROR in Assemble_edge_integrals (test function neighbour - ansatz function cell)" << endl);
                  //exit(4711);
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)
            if (verbose ==1){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from neighbour
            if(verbose ==1){     OutPut("testfct. neighbour <-> ansatzfct. neighbour Integral:" << endl);}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound) continue;

              //################################################################################################
              //#### Check if test function ii  of the neighbour cell is in the FE-Space of the cell     ####
              dummy = 0;
              l=0;
              while(l<N_ && dummy == 0)
              {
                if(dof_ii == DOF[l])dummy=1;
                l++;
              }
              if (dummy == 1) continue;
              //################################################################################################

              for (jj=0;jj<N_Neigh;jj++)
              {
                dof_jj = DOF_neigh[jj];

                //################################################################################################
                //####  Check if ansatz  function jj  of the neighbour cell is in the FE-Space of the cell  ####
                dummy = 0;
                l=0;
                while(l<N_ && dummy == 0)
                {
                  if(dof_jj == DOF[l])dummy=1;
                  l++;
                }
                if (dummy == 1) continue;
                //################################################################################################

                if(dummy == 0)                    //
                  // OutPut("jj " << dof_jj << endl);
                  // initialize the boundary integrals
                  integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch(n)
                  {
                    case 0:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*nx*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii+N_]);
                      break;

                    case 1:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj+N_]*jump_xderiv[j][ii+N_]);
                      break;

                    case 2:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*nx*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_xderiv[j][jj+N_]*jump_yderiv[j][ii+N_]);
                      break;

                    case 3:
                      //stabilization term with n*grad(u)
                      integrant = tau_par*Re*hk*hk*max_nb*ny*ny*(jump_xderiv[j][jj+N_]*jump_xderiv[j][ii+N_] + jump_yderiv[j][ii+N_] * jump_yderiv[j][jj+N_]);
                      //stabilization term with
                      integrant+=tau_par2*Re*hk*hk*max_b*(jump_yderiv[j][jj+N_]*jump_yderiv[j][ii+N_]);
                      break;

                   case 4:
                      //pressure stabilization term
                      //  integrant = tau_par3*Re*hk*hk*(jump_xderiv[j][jj]*jump_xderiv[j][ii] + jump_yderiv[j][ii] * jump_yderiv[j][jj]);
                      if(max_b>1e-10) 
		      { if(1/max_b < hk*TDatabase::ParamDB->RE_NR)
			  integrant=tau_par3/max_b*hk*hk*(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			      *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);

                        else  integrant=tau_par3*hk*TDatabase::ParamDB->RE_NR*hk*hk*(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			      *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);


                      }
                      else
                        integrant = tau_par3*hk*hk*hk*TDatabase::ParamDB->RE_NR
			    *(nx*jump_xderiv[j][jj+N_]+ny*jump_yderiv[j][jj+N_])
			    *(nx*jump_xderiv[j][ii+N_]+ny*jump_yderiv[j][ii+N_]);
                      break;
                  }

                  w = weights1D[j]* hE/2;//AbsDetjk1D[r][j];
                  //w = weights1D[j]* AbsDetjk1D[r][j];
		  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=SqRowPtr[n][dof_ii];ll < SqRowPtr[n][dof_ii+1]; ll++)
                {
                  if (SqColInd[n][ll] == dof_jj)
                  {
                    if(verbose ==1)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    SqEntries[n][ll] += integral;
                    found = 1;
                    break;
                  }
                }
                if (!found)
                {
                  OutPut("ERROR in Assemble_edge_integrals (test function neighbour - ansatz function neighbour)" << endl);
                  //exit(4711);
                }
                // update second matrix
              }                                   // end inner loop over dof (jj)
            }                                     // end outer loop over dof (ii)

            for (m=0;m<4;m++)
            {
              delete X1D_neigh[m];
              delete Y1D_neigh[m];
            }                                     //endfor m

          }                                       //endif i<q
        }                                         //endif neigh
      }                                           //endfor r (loop over edges)
    }                                             // endfor n (loop over sqmatrices)
  }                                               // endfor i (loop over cells)

  if(n_sqmatrices)
  {
    delete GlobalNumbers;
    delete BeginIndex;

    for(i=0;i<n_sqmatrices;i++)
      delete HangingEntries[i];
    delete HangingEntries;
  }

  if(n_matrices)
  {
    delete AnsatzGlobalNumbers;
    delete AnsatzBeginIndex;
    delete TestGlobalNumbers;
    delete TestBeginIndex;
  }

  /*if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete HangingRhs[i];
    delete HangingRhs;

    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }*/

  if(N_Parameters)
  {
    delete Param[0];
  }

  if(N_AllMatrices)
  {
    delete LocMatrices;
    delete Matrices[0];
    delete Matrices;
  }

  delete Coeffs[0];
  delete Coefficients1D[0];
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
      delete xyval_ref1D_test[i][j];
      delete xderiv_ref1D_test[i][j];
      delete yderiv_ref1D_test[i][j];
      delete xyval_ref1D_ansatz[i][j];
      delete xderiv_ref1D_ansatz[i][j];
      delete yderiv_ref1D_ansatz[i][j];
    }
  }

  /*for (n=0;n<N_BaseFuncts2D;n++)
    { delete xi1DNeigh[n];
      delete eta1DNeigh [n];
    }*/

  //delete X1DNeigh;
  //delete Y1DNeigh;
  for (i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete  xderiv_refNeigh1D[i];
    delete   yderiv_refNeigh1D[i];
  }
  /*        delete xderiv_Neigh1D;
          delete  yderiv_Neigh1D;
          delete  xyval_Neigh1D;*/
  // delete aux;
  //delete aux3;

  /*for(i=0; i < N_BaseFuncts2D; i++)
  {
      for(j=0; j < 4; j++)
      {
          for(m=0; m < MaxN_QuadPoints_1D; m++)
          {
              //for(n=0; n < MaxN_BaseFunctions2D_Ersatz; n++)
              //{
                 delete xietaval_ref1D[i][j][m];
                 delete xideriv_ref1D[i][j][m];
                 delete etaderiv_ref1D[i][j][m];
  // }
  }
  }
  }*/

  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  if(verbose==2)
  {
    OutPut("Squarematrices after Edge assembling" << endl);
    for(k=0;k<n_sqmatrices;k++)
    {
      cout << endl;
      cout << "sqmatrix: " << k << endl;
      SqRowPtr[k] = sqmatrices[k]->GetRowPtr();
      SqEntries[k] = sqmatrices[k]->GetEntries();
      SqColInd[k] = sqmatrices[k]->GetKCol();
      N_Rows = sqmatrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=SqRowPtr[k][i+1];
        for(j=SqRowPtr[k][i];j<end;j++)
        {
          // cout << j << endl;
          //cout << "Matrix: " << setw(5) << i << setw(5) << SqColInd[k][j] << "   ";
          //cout << setw(10) << SqEntries[k][j] << endl;
          cout << "b(" << i+1 << "," << SqColInd[k][j]+1 << ") = ";
          cout << SqEntries[k][j] << ";" << endl;
        }
      }
      cout << endl;
    }                                             // endfor k
  }
}                                                 // end of Assemble

#endif
