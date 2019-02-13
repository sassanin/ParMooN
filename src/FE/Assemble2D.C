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
// @(#)Assemble2D.C        1.16 04/13/00
//
// Purpose:     bilinear form (discretized and stabilized assemble)
//
// Author:      Gunar Matthies (10.08.98)
//
// History:     start of implementation 10.08.98 (Gunar Matthies)
//
// =======================================================================

#include <DefineParams.h>

#include <Assemble2D.h>
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
#include <vector>


void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
                int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                int n_matrices, TMatrix2D **matrices,
                int n_rhs, double **rhs, TFESpace2D **ferhs,
                TDiscreteForm2D *DiscreteForm,
                BoundCondFunct2D **BoundaryConditions,
                BoundValueFunct2D **BoundaryValues,
                TAuxParam2D *Parameters
#ifdef __3D__
                , TAux2D3D *Aux2D3D
#endif
                , int AssemblePhaseID
)
{
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ii,jj,ll,ij;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, N_Edges;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
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
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s,integral[2];
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-4;
#ifdef __3D__
  double z0, z1;
#endif
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue;
  bool *SecondDer;
  int *Bounds, NeumannBound, RobinBound;
  int lr, mr;
  double RobinScale = 16.0;
  TVertex *ver0;
  TRefTrans2D *rt;

  int N_Rows;
  int AnsatzActiveBound, AnsatzHangingBound;
 
  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
 
  if(n_sqmatrices+n_matrices)
  {
    HangingEntries = new double* [n_sqmatrices+n_matrices];

    for(i=0;i<n_sqmatrices;i++)
    {
      j = sqmatrices[i]->GetHangingN_Entries();
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    }

    for(i=0;i<n_matrices;i++)
    {
      j = matrices[i]->GetHangingN_Entries();
      HangingEntries[i+n_sqmatrices] = new double [j];
      memset(HangingEntries[i+n_sqmatrices], 0, SizeOfDouble*j);
    }
  }
 
  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                              // endif n_sqmatrices

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

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  N_Parameters = Parameters->GetN_Parameters();
  
#ifdef __3D__
  N_Parameters += 7;                              // (u, ux, uy, uz, nx, ny, nz)
#endif

  if(N_Parameters)
  {
    aux = new double [MaxN_QuadPoints_2D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_2D;j++)
      Param[j] = aux + j*N_Parameters;
  }

  // 40 <= number of terms in bilinear form
  // DUE NOTE CHANGE BELOW 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
  aux = new double [MaxN_QuadPoints_2D*40];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*40;

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices
 
  SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++) // set the cell indices
  {
    cell = Coll->GetCell(i);
    cell->SetCellIndex(i);
  }

  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    TDatabase::ParamDB->INTERNAL_CELL = i;

    // only for multiphase flows
    if(AssemblePhaseID >= 0)
     {
      if(AssemblePhaseID != cell->GetRegionID() )
       continue;
     }

    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
      case 0:                                     // diameter
        hK = cell->GetDiameter();
        break;
        //case 1: // with reference map
        //OutPut("cell measure " << endl);
        //hK = cell->GetLengthWithReferenceMap();
        //break;
      case 2:                                     // shortest edge
        hK = cell->GetShortestEdge();
        break;
      case 1:                                     // with reference map
      case 3:                                     // measure
        hK = cell->GetMeasure();
        hK = sqrt(hK);
        break;
      case 4:                                     // mesh size in convection direction, this is just a dummy
        hK = cell->GetDiameter();
        break;
      case 5:                                     // take value from an array
        // this is in general not the diameter but a pw constant value
        // which is needed for some reasons
        hK = cell->GetDiameter();
        break;
      default:                                    // diameter
        hK = cell->GetDiameter();
        OutPut("CELL_MEASURE " << TDatabase::ParamDB->CELL_MEASURE <<
             " not available, set CELL_MEASURE: 0 !!!" << endl);
        TDatabase::ParamDB->CELL_MEASURE = 0;
        break;
    }   
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
 
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
    
    Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);
 
    if((TDatabase::ParamDB->DISCTYPE == SDFEM)
      || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
      || (TDatabase::ParamDB->CELL_MEASURE == 4))
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
      N_Edges = cell->GetN_Edges();
      for (ij=0;ij<N_Edges;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if (N_Edges==3)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }

#ifdef __3D__
    if(Aux2D3D)
      Aux2D3D->GetGradient(i, N_Points, xi, eta, Param);
#endif
    // use DiscreteForm to assemble a few matrices and
    // right-hand sides at once

    if(DiscreteForm)
    {
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk,
        hK, X, Y,
        LocN_BF, LocBF,
        Param, AuxArray,
        cell,
        N_AllMatrices, n_rhs,
        LocMatrices, LocRhs);
    }

    N_Joints = cell->GetN_Joints();
    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################
    for(j=0;j<n_sqmatrices;j++) 
    {
      // find space for this bilinear form
      fespace = sqmatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE2D(i, cell);
      N_ = N_BaseFunct[CurrentElement];

      Matrix = Matrices+j*MaxN_BaseFunctions2D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr = sqmatrices[j]->GetRowPtr();
      ColInd = sqmatrices[j]->GetKCol();

      CurrentHangingEntries = HangingEntries[j];
      HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
      HangingColInd = sqmatrices[j]->GetHangingKCol();

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = GlobalNumbers[j] + BeginIndex[j][i];
      Bounds = fespace->GetBoundaryNodesBound();
      NeumannBound = Bounds[0];
      RobinBound = Bounds[1]; 

      /*
      BoundaryCondition = BoundaryConditions[j];
      for(m=0;m<N_Joints;m++)
      {
      joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge ||
           joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge)
          {
            boundedge = (TBoundEdge *)joint;
      BoundComp = boundedge->GetBoundComp();
      boundedge->GetParameters(t0, t1);
      }
      else
      {
      isoboundedge = (TIsoBoundEdge *)joint;
      BoundComp = isoboundedge->GetBoundComp();
      isoboundedge->GetParameters(t0, t1);
      }
      // get id of the boundary component
      comp = BoundComp->GetID();
      // get type of the boundary condition at the beginning
      // and at the end of the current edge
      BoundaryCondition(comp, t0, Cond0);
      //      cout << "bound1" << endl;
      if(Cond0 == ROBIN)
      {
      #ifdef __2D__
      //cout << "robin" << endl;
      // Robin
      lr = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);

      // get a suitable line quadrature formula
      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*lr);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

      TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);

      JointValues=TFEDatabase2D::GetJointValues2D(
      BaseFuncts[CurrentElement], LineQuadFormula, m);
      // get vertices of boundary edge
      cell->GetVertex(m)->GetCoords(x0, y0);
      cell->GetVertex((m+1) % 4)->GetCoords(x1, y1);
      // compute (half of the) length of the boundary edge
      hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
      // cout << "x0: " << x0 << " y0: " << y0 << endl;
      // cout << "x1: " << x1 << " y1: " << y1 << endl;
      // compute boundary integral
      for(n=0;n<N_LinePoints;n++)
      {
      // values of test functions in this quadrature point
      JointValue = JointValues[n];
      // cout << "Zeta  :" << zeta[n] << endl;

      // get quadrature point on the boundary
      for(l=0;l<N_;l++)
      {
      MatrixRow = Matrix[l];
      s = JointValue[l];
      // multiply value with weights from quadrature formula
      // and determinant from integral transformation to the
      // unit edge (-1,1)
      s *= hE * LineWeights[n];
      // !! hold the robin boundary values of
      // !! the function alpha from parameters
      // s *= alpha;
      // update rhs for all test functions
      for(k=0;k<N_;k++)
      MatrixRow[k] += s*JointValue[k]*RobinScale;
      } // endfor l
      } // endfor n
      #endif
      } // end Robin
      } // endif BoundEdge
      } // endfor m
      */

      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        if(l<ActiveBound)
        {
          for(k=0;k<N_;k++)
          {
            // DOF[k] is the global index of the k-th local degree of freedom
            // MatrixRow[k] is the assembled value corresponding to the m-th
            // local test function and k-th local ansatz function. That means it
            // corresponds to the l=DOF[m]-th global test function and the 
            // DOF[k]-th global ansatz function
	
             sqmatrices[j]->add(l, DOF[k], MatrixRow[k]); 
	 
          }
          /*
          // node l is inner or Neumann node
          end=RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            for(k=0;k<N_;k++)
            {
              if(DOF[k] == ColInd[n])
              {
                //cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
          */
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            end = HangingRowPtr[l+1];
            for(n=HangingRowPtr[l];n<end;n++)
            {
              for(k=0;k<N_;k++)
              {
                if(DOF[k] == HangingColInd[n])
                {
                  CurrentHangingEntries[n] += MatrixRow[k];
                  break;
                }                                 // endif
              }                                   // endfor k
            }                                     // endfor n
          }
          else
          {
            // Dirichlet node
            sqmatrices[j]->set(l,l,1.0);
            /*n=RowPtr[l];
            if(ColInd[n]==l)
            {
              Entries[n]=1.0;
            }*/
          }
        }
      }                                           // endfor m
    }                                             // endfor j

    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(j=0;j<n_matrices;j++)   
    {
      TestElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetTestSpace())->GetFE2D(i, cell);
      AnsatzElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetAnsatzSpace())->GetFE2D(i, cell);
      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->GetRowPtr();
      ColInd = matrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      fespace = (TFESpace2D *)(matrices[j]->GetStructure()->GetTestSpace());
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
      HangingRowPtr = matrices[j]->GetHangingRowPtr();
      HangingColInd = matrices[j]->GetHangingKCol();

      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound || l>=DirichletBound)
        {
          end=RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            for(k=0;k<N_Ansatz;k++)
            {
              if(AnsatzDOF[k] == ColInd[n])
              {
                // cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }
        else
        {
          // hanging node
          l -= ActiveBound;
          end = HangingRowPtr[l+1];
          for(n=HangingRowPtr[l];n<end;n++)
          {
            // cout << l << " HangingColInd: " << HangingColInd[n] << endl;
            for(k=0;k<N_Ansatz;k++)
            {
              // cout << "AnsatzDOF: " << AnsatzDOF[k] << endl;
              if(AnsatzDOF[k] == HangingColInd[n])
              {
                CurrentHangingEntries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)

    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      fespace = ferhs[j];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE2D(i, cell);

      N_ = N_BaseFunct[CurrentElement];

      local_rhs = righthand+j*MaxN_BaseFunctions2D;
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];
      // find space for this linear form

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
	if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
	  {
	    l = 0;
	    N_Edges=cell->GetN_Edges();
	    for(jj=0;jj<N_Edges;jj++)                        // loop over all edges of cell
	      {
		ver0 = cell->GetVertex(jj);
		t0 =  ver0->GetClipBoard();
		// vertex not on the boundary
		if (t0<-1e-8)
		  continue;
		// component of boundary
		comp = floor(t0+1e-8);
		// parameter
		t0 -= comp;
		// get boundary condition
		BoundaryConditions[0](comp, t0, Cond0);
		// Dirichlet
		if (Cond0== DIRICHLET)
		  {
		    l = -4711;
		  }
		/*joint=cell->GetJoint(jj);
		if ((joint->GetType() == BoundaryEdge)||
		    (joint->GetType() == IsoBoundEdge))       // boundary edge
		  {
		    boundedge=(TBoundEdge *)joint;
		    BoundComp=boundedge->GetBoundComp();      // get boundary component
		    boundedge->GetParameters(t0, t1);         // parameter interval
		    comp=BoundComp->GetID();                  // boundary id
		    BoundaryConditions[0](comp, (t0+t1)/2.0, Cond0);
		    // at midpoint of boundary
		    if (Cond0==DIRICHLET)
		      l = -4711;
		      }*/
	      }
	    if (l==-4711) 
	      break; // do nothing for this mesh cell
	  }
		
        l=DOF[m];
        //cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
          // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      }                                           // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase2D::GetFE2D(CurrentElement);
      //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
      //  continue;

      nf = ele->GetNodalFunctional2D();

      if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
      {
        /* Superconvergence boundary interpolation */
        if(nf->GetID() == NF_C_Q_Q2_2D)
          nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
      }

      nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      FEDesc_Obj = ele->GetFEDesc2D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Edges();

      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);

        if(joint->GetType() == BoundaryEdge ||
           joint->GetType() == IsoBoundEdge ||
           joint->GetType() == InterfaceJoint)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            switch(Cond0)
            {
              case DIRICHLET:
                // if DG
                if (N_EdgeDOF==0)
                  break;
                // read boundary values for each quadrature point
                for(l=0;l<N_EdgePoints;l++)
                {
                  s = EdgePoints[l];
                  t = 0.5*(t0*(1-s) + t1*(1+s));
                  BoundaryValue(comp, t, PointValues[l]);
                }                                 // endfor l
                // compute boundary values for each dof on the
                // boundary edge with the nodal functionals

                nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                  FunctionalValues);
                EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                // save boundary values of each dof on the boundary
                // edge in the rhs
                for(l=0;l<N_EdgeDOF;l++)
                {
                  RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                }
                break;

              case NEUMANN:
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[l];
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                break;
              case ROBIN:
#ifdef __2D__
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                // get vertices of boundary edge
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[l];
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
#endif
                break;
              case SLIP:
                OutPut("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !"<< endl);
                exit(4711);
                break;
              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble2DSlipBC, see below
                break;

            case FREESURF:
            case INTERFACE:
              // do nothing here
              // everything is done in Freesurfint, see Freesurface2D.C
              break;
              default :
                OutPut("Unknown boundary condition !"<< endl);
                exit(4711);

            }                                     // endswitch Cond0
          }                                       // endif (Cond0==Cond1)
          else
          {
            OutPut("different boundary condition on one edge ");
            OutPut("are not allowed!" << endl);
            exit(4711);
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)


  
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j];
    HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
    HangingColInd = sqmatrices[j]->GetHangingKCol();

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      end = HangingRowPtr[i+1];
      for(n=HangingRowPtr[i];n<end;n++)
      {
        v = CurrentHangingEntries[n];
        m = HangingColInd[n];
        for(l=0;l<k;l++)
        {
          l1 = DOF[l];
          if(l1<ActiveBound)
          {
            last=RowPtr[l1+1];
            for(l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i
  }                                               // endfor j

  for(j=0;j<n_matrices;j++)
  {
    // hanging nodes in test space
    fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetTestSpace());
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = matrices[j]->GetEntries();
    RowPtr = matrices[j]->GetRowPtr();
    ColInd = matrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
    HangingRowPtr = matrices[j]->GetHangingRowPtr();
    HangingColInd = matrices[j]->GetHangingKCol();

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      end = HangingRowPtr[i+1];
      for(n=HangingRowPtr[i];n<end;n++)
      {
        v = CurrentHangingEntries[n];
        m = HangingColInd[n];
        for(l=0;l<k;l++)
        {
          l1 = DOF[l];
          if(l1<ActiveBound)
          {
            last=RowPtr[l1+1];
            for(l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i

    // hanging nodes in ansatz space
    N_Rows =  matrices[j]->GetN_Rows();
    fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetAnsatzSpace());

    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();
    AnsatzActiveBound = fespace->GetActiveBound();
    AnsatzHangingBound = fespace->GetHangingBound();
    for(i=0;i<N_Rows;i++)
    {
      end = RowPtr[i+1];
      for(k=RowPtr[i];k<end;k++)
      {
        l = ColInd[k];
        if(l>=AnsatzActiveBound && l<AnsatzHangingBound)
        {
          // l is hanging node in ansatz space
          hn = HangingNodes[l-AnsatzActiveBound];
          HNDescr = hn->GetType();
          HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
          m = HNDescr_Obj->GetN_Nodes();
          Coupling = HNDescr_Obj->GetCoeff();
          DOF = hn->GetDOF();
          v = Entries[k];
          for(n=0;n<m;n++)
          {
            for(l1=RowPtr[i];l1<end;l1++)
            {
              if(ColInd[l1] == DOF[n])
                Entries[l1] += v*Coupling[n];
            }
          }
          Entries[k] = 0;
        }                                         // endif l
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  for(j=0;j<n_rhs;j++)
  {
    fespace = ferhs[j];
    N_Hanging = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    RHS = rhs[j];
    CurrentHangingRhs = HangingRhs[j];

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      N_ = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      for(k=0;k<N_;k++)
      {
        l = DOF[k];
        if(l<ActiveBound)
        {
          RHS[l] += Coupling[k] * CurrentHangingRhs[i];
        }
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    ActiveBound = fespace->GetActiveBound();

    n = RowPtr[ActiveBound];

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      Entries[n] = 1.0;
      n++;

      for(l=0;l<k;l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      }                                           // endfor l

    }                                             // endfor i
  }

  if(n_sqmatrices)
  {
    delete [] GlobalNumbers;
    delete [] BeginIndex;
  }

  if(n_matrices)
  {
    delete [] AnsatzGlobalNumbers;
    delete [] AnsatzBeginIndex;
    delete [] TestGlobalNumbers;
    delete [] TestBeginIndex;
  }

  if(n_sqmatrices+n_matrices)
  {
    for(i=0;i<n_sqmatrices+n_matrices;i++)
      delete [] HangingEntries[i];

    delete [] HangingEntries;
  }

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete [] HangingRhs[i];
    delete [] HangingRhs;

    delete [] righthand;
    delete [] LocRhs;
    delete [] RhsBeginIndex;
    delete [] RhsGlobalNumbers;
  }

  if(N_Parameters)
  {
    delete [] Param[0];
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];

  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  /*
    for(k=0;k<n_sqmatrices;k++)
    {
      cout << endl;
      cout << "sqmatrix: " << k << endl;
      RowPtr = sqmatrices[k]->GetRowPtr();
      Entries = sqmatrices[k]->GetEntries();
      ColInd = sqmatrices[k]->GetKCol();
      N_Rows = sqmatrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  OutPut("a"<<k<<"("<<i+1<<","<< ColInd[j]+1 << ")=   ");
  if (fabs(Entries[j])<1e-12)
  Entries[j] = 0;
  OutPut(Entries[j] << ";" << endl);
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_matrices;k++)
  {
  cout << endl;
  cout << "matrix: " << k << endl;
  RowPtr = matrices[k]->GetRowPtr();
  Entries = matrices[k]->GetEntries();
  ColInd = matrices[k]->GetKCol();
  N_Rows = matrices[k]->GetN_Rows();
  for(i=0;i<N_Rows;i++)
  {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  OutPut("b"<<k<<"("<<i+1<<","<< ColInd[j]+1 << ")=   ");
  if (fabs(Entries[j])<1e-12)
  Entries[j] = 0;
  OutPut(Entries[j] << ";" << endl);
  //cout << setw(10) << i << setw(10) << ColInd[j] << "   ";
  //cout << setw(20) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

//   for(k=0;k<n_rhs;k++)
//   {
//   cout << "rhs: " << k << endl;
//   N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
//   RHS=rhs[k];
//   for(i=0;i<N_Rows;i++)
//   //cout << setw(10) << i << setw(20) << RHS[i] << endl;
//   OutPut("f"<<k<<"("<< i+1  << ")=   " << RHS[i] << endl);
//   }
  */
}                                                 // end of Assemble


////////////////////////////////////////////////////////////////////////
//
// Assembling of matrices multiplied by a factor
//
////////////////////////////////////////////////////////////////////////

void Assemble2D_FCT(int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
TDiscreteForm2D *DiscreteForm,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters,
double factor
#ifdef __3D__
, TAux2D3D *Aux2D3D
#endif
)
{
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ij;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, N_Edges;
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
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
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries;
  int *ColInd, *RowPtr;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s;
  int comp;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, eps=1e-4;
  #ifdef __3D__
  double z0, z1;
  #endif
  double **JointValues, *JointValue;
  bool *SecondDer;
  //int *Bounds, NeumannBound, RobinBound;
  //double RobinScale = 16.0;
  
  int N_Rows;
  int AnsatzActiveBound, AnsatzHangingBound;

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices+n_matrices)
  {
    HangingEntries = new double* [n_sqmatrices+n_matrices];

    for(i=0;i<n_sqmatrices;i++)
    {
      j = sqmatrices[i]->GetHangingN_Entries();
      HangingEntries[i] = new double [j];
      memset(HangingEntries[i], 0, SizeOfDouble*j);
    }

    for(i=0;i<n_matrices;i++)
    {
      j = matrices[i]->GetHangingN_Entries();
      HangingEntries[i+n_sqmatrices] = new double [j];
      memset(HangingEntries[i+n_sqmatrices], 0, SizeOfDouble*j);
    }
  }

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_sqmatrices

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

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  N_Parameters = Parameters->GetN_Parameters();

  #ifdef __3D__
  N_Parameters += 7;                              // (u, ux, uy, uz, nx, ny, nz)
  #endif

  if(N_Parameters)
  {
    aux = new double [MaxN_QuadPoints_2D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_2D;j++)
      Param[j] = aux + j*N_Parameters;
  }

  // 40 <= number of terms in bilinear form
  // DUE NOTE CHANGE BELOW 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
  aux = new double [MaxN_QuadPoints_2D*40];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*40;

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    cell->SetCellIndex(i);
  }
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    TDatabase::ParamDB->INTERNAL_CELL = i;
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
      case 0:                                     // diameter
        hK = cell->GetDiameter();
        break;
        //case 1: // with reference map
        //OutPut("cell measure " << endl);
        //hK = cell->GetLengthWithReferenceMap();
        //break;
      case 2:                                     // shortest edge
        hK = cell->GetShortestEdge();
        break;
      case 1:                                     // with reference map
      case 3:                                     // measure
        hK = cell->GetMeasure();
        hK = sqrt(hK);
        break;
      case 4:                                     // mesh size in convection direction, this is just a dummy
        hK = cell->GetDiameter();
        break;
      case 5:                                     // take value from an array
        // this is in general not the diameter but a pw constant value
        // which is needed for some reasons
        hK = cell->GetDiameter();
        break;
      default:                                    // diameter
        hK = cell->GetDiameter();
        break;
    }
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################

    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

    Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

    if ((TDatabase::ParamDB->DISCTYPE == SDFEM)||
      (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM))
    {
      TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
      N_Edges = cell->GetN_Edges();
      for (ij=0;ij<N_Edges;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
      }
      if (N_Edges==3)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }

    #ifdef __3D__
    if(Aux2D3D)
      Aux2D3D->GetGradient(i, N_Points, xi, eta, Param);
    #endif

    // use DiscreteForm to assemble a few matrices and
    // right-hand sides at once

    if(DiscreteForm)
    {
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk,
        hK, X, Y,
        LocN_BF, LocBF,
        Param, AuxArray,
        cell,
        N_AllMatrices, n_rhs,
        LocMatrices, LocRhs, factor);
    }
    N_Joints = cell->GetN_Joints();

    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################
    for(j=0;j<n_sqmatrices;j++)
    {
      // find space for this bilinear form
      fespace = sqmatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE2D(i, cell);
      N_ = N_BaseFunct[CurrentElement];

      Matrix = Matrices+j*MaxN_BaseFunctions2D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr = sqmatrices[j]->GetRowPtr();
      ColInd = sqmatrices[j]->GetKCol();

      CurrentHangingEntries = HangingEntries[j];
      HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
      HangingColInd = sqmatrices[j]->GetHangingKCol();

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = GlobalNumbers[j] + BeginIndex[j][i];
      //Bounds = fespace->GetBoundaryNodesBound();
      //NeumannBound = Bounds[0];
      //RobinBound = Bounds[1];

      /*
        BoundaryCondition = BoundaryConditions[j];
        for(m=0;m<N_Joints;m++)
        {
        joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge ||
        joint->GetType() == IsoBoundEdge)
        {
        if(joint->GetType() == BoundaryEdge)
        {
              boundedge = (TBoundEdge *)joint;
      BoundComp = boundedge->GetBoundComp();
      boundedge->GetParameters(t0, t1);
      }
      else
      {
      isoboundedge = (TIsoBoundEdge *)joint;
      BoundComp = isoboundedge->GetBoundComp();
      isoboundedge->GetParameters(t0, t1);
      }
      // get id of the boundary component
      comp = BoundComp->GetID();
      // get type of the boundary condition at the beginning
      // and at the end of the current edge
      BoundaryCondition(comp, t0, Cond0);
      //      cout << "bound1" << endl;
      if(Cond0 == ROBIN)
      {
      #ifdef __2D__
      //cout << "robin" << endl;
      // Robin
      lr = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);

      // get a suitable line quadrature formula
      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*lr);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

      TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
      ->MakeRefElementData(LineQuadFormula);

      JointValues=TFEDatabase2D::GetJointValues2D(
      BaseFuncts[CurrentElement], LineQuadFormula, m);
      // get vertices of boundary edge
      cell->GetVertex(m)->GetCoords(x0, y0);
      cell->GetVertex((m+1) % 4)->GetCoords(x1, y1);
      // compute (half of the) length of the boundary edge
      hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
      // cout << "x0: " << x0 << " y0: " << y0 << endl;
      // cout << "x1: " << x1 << " y1: " << y1 << endl;
      // compute boundary integral
      for(n=0;n<N_LinePoints;n++)
      {
      // values of test functions in this quadrature point
      JointValue = JointValues[n];
      // cout << "Zeta  :" << zeta[n] << endl;

      // get quadrature point on the boundary
      for(l=0;l<N_;l++)
      {
      MatrixRow = Matrix[l];
      s = JointValue[l];
      // multiply value with weights from quadrature formula
      // and determinant from integral transformation to the
      // unit edge (-1,1)
      s *= hE * LineWeights[n];
      // !! hold the robin boundary values of
      // !! the function alpha from parameters
      // s *= alpha;
      // update rhs for all test functions
      for(k=0;k<N_;k++)
      MatrixRow[k] += s*JointValue[k]*RobinScale;
      } // endfor l
      } // endfor n
      #endif
      } // end Robin
      } // endif BoundEdge
      } // endfor m
      */

      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        //cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          end=RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            for(k=0;k<N_;k++)
            {
              if(DOF[k] == ColInd[n])
              {
                //cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            end = HangingRowPtr[l+1];
            for(n=HangingRowPtr[l];n<end;n++)
            {
              for(k=0;k<N_;k++)
              {
                if(DOF[k] == HangingColInd[n])
                {
                  CurrentHangingEntries[n] += MatrixRow[k];
                  break;
                }                                 // endif
              }                                   // endfor k
            }                                     // endfor n
          }
          else
          {
            // Dirichlet node
            n=RowPtr[l];
            if(ColInd[n]==l)
            {
              Entries[n]=1.0;
            }
          }
        }
      }                                           // endfor m
    }                                             // endfor j

    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(j=0;j<n_matrices;j++)
    {
      TestElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetTestSpace())->GetFE2D(i, cell);
      AnsatzElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetAnsatzSpace())->GetFE2D(i, cell);
      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->GetRowPtr();
      ColInd = matrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      fespace = (TFESpace2D *)(matrices[j]->GetStructure()->GetTestSpace());
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
      HangingRowPtr = matrices[j]->GetHangingRowPtr();
      HangingColInd = matrices[j]->GetHangingKCol();

      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound || l>=DirichletBound)
        {
          end=RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            for(k=0;k<N_Ansatz;k++)
            {
              if(AnsatzDOF[k] == ColInd[n])
              {
                // cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }
        else
        {
          // hanging node
          l -= ActiveBound;
          end = HangingRowPtr[l+1];
          for(n=HangingRowPtr[l];n<end;n++)
          {
            // cout << l << " HangingColInd: " << HangingColInd[n] << endl;
            for(k=0;k<N_Ansatz;k++)
            {
              // cout << "AnsatzDOF: " << AnsatzDOF[k] << endl;
              if(AnsatzDOF[k] == HangingColInd[n])
              {
                CurrentHangingEntries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)
    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      fespace = ferhs[j];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE2D(i, cell);

      N_ = N_BaseFunct[CurrentElement];

      local_rhs = righthand+j*MaxN_BaseFunctions2D;
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];
      // find space for this linear form

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
        if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
        {
          // do nothing for this cell
          if (TDatabase::ParamDB->INTERNAL_INDICATOR_Array[i]!=1)
          {
            break;
          }
        }

        l=DOF[m];
        //cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
          // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      }                                           // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase2D::GetFE2D(CurrentElement);
      //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
      //  continue;

      nf = ele->GetNodalFunctional2D();

      if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
      {
        /* Superconvergence boundary interpolation */
        if(nf->GetID() == NF_C_Q_Q2_2D)
          nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
      }

      nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      FEDesc_Obj = ele->GetFEDesc2D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Edges();

      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);

        if(joint->GetType() == BoundaryEdge ||
          joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            switch(Cond0)
            {
              case DIRICHLET:
                // if DG
                if (N_EdgeDOF==0)
                  break;
                // read boundary values for each quadrature point
                for(l=0;l<N_EdgePoints;l++)
                {
                  s = EdgePoints[l];
                  t = 0.5*(t0*(1-s) + t1*(1+s));
                  BoundaryValue(comp, t, PointValues[l]);
                }                                 // endfor l
                // compute boundary values for each dof on the
                // boundary edge with the nodal functionals

                nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                  FunctionalValues);
                EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                // save boundary values of each dof on the boundary
                // edge in the rhs
                for(l=0;l<N_EdgeDOF;l++)
                {
                  RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                }
                break;

              case NEUMANN:
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                // get vertices of boundary edge
                #ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
                #else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                #endif
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[l];
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                break;
              case ROBIN:
                #ifdef __2D__
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                // get vertices of boundary edge
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(l=0;l<N_LinePoints;l++)
                {
                  // values of test functions in this quadrature point
                  JointValue = JointValues[l];
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[l];
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    if((l3 = DOF[k])<ActiveBound)
                      RHS[l3] += s*JointValue[k];
                }
                #endif
                break;
              case SLIP:
                OutPut("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !"<< endl);
                exit(4711);
                break;
              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble2DSlipBC, see below
                break;
              default :
                OutPut("Unknown boundary condition !"<< endl);
                exit(4711);

            }                                     // endswitch Cond0
          }                                       // endif (Cond0==Cond1)
          else
          {
            OutPut("different boundary condition on one edge ");
            OutPut("are not allowed!" << endl);
            exit(4711);
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)
  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j];
    HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
    HangingColInd = sqmatrices[j]->GetHangingKCol();

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      end = HangingRowPtr[i+1];
      for(n=HangingRowPtr[i];n<end;n++)
      {
        v = CurrentHangingEntries[n];
        m = HangingColInd[n];
        for(l=0;l<k;l++)
        {
          l1 = DOF[l];
          if(l1<ActiveBound)
          {
            last=RowPtr[l1+1];
            for(l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i
  }                                               // endfor j

  for(j=0;j<n_matrices;j++)
  {
    // hanging nodes in test space
    fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetTestSpace());
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = matrices[j]->GetEntries();
    RowPtr = matrices[j]->GetRowPtr();
    ColInd = matrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
    HangingRowPtr = matrices[j]->GetHangingRowPtr();
    HangingColInd = matrices[j]->GetHangingKCol();

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      end = HangingRowPtr[i+1];
      for(n=HangingRowPtr[i];n<end;n++)
      {
        v = CurrentHangingEntries[n];
        m = HangingColInd[n];
        for(l=0;l<k;l++)
        {
          l1 = DOF[l];
          if(l1<ActiveBound)
          {
            last=RowPtr[l1+1];
            for(l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i

    // hanging nodes in ansatz space
    N_Rows =  matrices[j]->GetN_Rows();
    fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetAnsatzSpace());

    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();
    AnsatzActiveBound = fespace->GetActiveBound();
    AnsatzHangingBound = fespace->GetHangingBound();
    for(i=0;i<N_Rows;i++)
    {
      end = RowPtr[i+1];
      for(k=RowPtr[i];k<end;k++)
      {
        l = ColInd[k];
        if(l>=AnsatzActiveBound && l<AnsatzHangingBound)
        {
          // l is hanging node in ansatz space
          hn = HangingNodes[l-AnsatzActiveBound];
          HNDescr = hn->GetType();
          HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
          m = HNDescr_Obj->GetN_Nodes();
          Coupling = HNDescr_Obj->GetCoeff();
          DOF = hn->GetDOF();
          v = Entries[k];
          for(n=0;n<m;n++)
          {
            for(l1=RowPtr[i];l1<end;l1++)
            {
              if(ColInd[l1] == DOF[n])
                Entries[l1] += v*Coupling[n];
            }
          }
          Entries[k] = 0;
        }                                         // endif l
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  for(j=0;j<n_rhs;j++)
  {
    fespace = ferhs[j];
    N_Hanging = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    RHS = rhs[j];
    CurrentHangingRhs = HangingRhs[j];

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      N_ = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      for(k=0;k<N_;k++)
      {
        l = DOF[k];
        if(l<ActiveBound)
        {
          RHS[l] += Coupling[k] * CurrentHangingRhs[i];
        }
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    ActiveBound = fespace->GetActiveBound();

    n = RowPtr[ActiveBound];

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      Entries[n] = 1.0;
      n++;

      for(l=0;l<k;l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      }                                           // endfor l

    }                                             // endfor i
  }

  if(n_sqmatrices)
  {
    delete GlobalNumbers;
    delete BeginIndex;
  }

  if(n_matrices)
  {
    delete AnsatzGlobalNumbers;
    delete AnsatzBeginIndex;
    delete TestGlobalNumbers;
    delete TestBeginIndex;
  }

  if(n_sqmatrices+n_matrices)
  {
    for(i=0;i<n_sqmatrices+n_matrices;i++)
      delete HangingEntries[i];

    delete HangingEntries;
  }

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete HangingRhs[i];
    delete HangingRhs;

    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }

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

  delete AuxArray[0];

} // end of Assemble2D_FCT



// =======================================================================
//
// Assemble2DSlipBC
//
// some manipulations in matrices and the rhs are necessary
//
// =======================================================================

void Assemble2DSlipBC(int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
TDiscreteForm2D *DiscreteForm,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters,
TFEFunction2D *u1, TFEFunction2D *u2)
{
  
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ii,jj,ll;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
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
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  double *Entries6, *Entries7, *Entries8, *Entries9;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5, *ColInd6, *RowPtr6, *ColInd7, *RowPtr7;
  int *ColInd8, *RowPtr8, *ColInd9, *RowPtr9;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s,integral[2];
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double friction_constant= TDatabase::ParamDB->FRICTION_CONSTANT;
  double friction_power = TDatabase::ParamDB->FRICTION_POWER;
  double penetration_constant = TDatabase::ParamDB->PENETRATION_CONSTANT;
  double penetration_power = TDatabase::ParamDB->PENETRATION_POWER;
  int friction_type = TDatabase::ParamDB->FRICTION_TYPE;
  double RE_NR, tangential_velo, U0, denominator;
  double delta, r_axial, beta_h;
//   double Frict_Force1=0., Frict_Force2=0.;
  bool *SecondDer;
#ifdef __3D__
  double z0, z1;
#endif

  int axial3D = TDatabase::ParamDB->Axial3D;

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];
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

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;

  }                                               // endif n_rhs

  N_Parameters = Parameters->GetN_Parameters();
  if(N_Parameters)
  {
    aux = new double [MaxN_QuadPoints_2D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_2D;j++)
      Param[j] = aux + j*N_Parameters;
  }

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_2D*40];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*40;

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  SecondDer = new bool[n_fespaces];
  SecondDer[0] = FALSE;

  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    hK = cell->GetDiameter();
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // ####################################################################
    // manipulate global matrices
    // manipulate global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      fespace = ferhs[j];
      CurrentElement = fespace->GetFE2D(i, cell);

      N_ = N_BaseFunct[CurrentElement];
      

      local_rhs = righthand+j*MaxN_BaseFunctions2D;
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];

      // find bounds in fe space
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
      
      // only for edges on the boundary
      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase2D::GetFE2D(CurrentElement);
      // nf = ele->GetNodalFunctional2D();
      // nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      FEDesc_Obj = ele->GetFEDesc2D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

      // find slip type bc
      N_Joints = cell->GetN_Edges();
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
           joint->GetType() == IsoBoundEdge)
         {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();                  
          
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }

          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            switch(Cond0)
            {
              case DIRICHLET:
                break;

              case NEUMANN:
                break;

            case FREESURF:
            case INTERFACE:
              // do nothing here
              // everything is done in Freesurfint, see Freesurface2D.C
              break;

              case SLIP:
                exit(4711);
                break;

              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // edge is assumed to be straight line
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                // get values of test functions in all quadrature points
                // on joint m
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                // compute normal vector to this boundary (normalized)
                nx = (y1-y0)/hE;
                ny = (x0-x1)/hE;
                // tangential normal vector to this boundary (normalized)
                tx = (x1-x0)/hE;
                ty = (y1-y0)/hE;  

//                   OutPut("A " << x0 << " " << y0 << " B " << x1 << " " << y1 << endl);
//                   OutPut("t " << tx << " " << ty << " n " << nx << " " << ny << endl);

//                 delta = CharacteristicFilterWidth(hK);
#ifdef  __CHANNELSTEPSLIP__
                // upper boundary - free slip
                if (comp==6)
                  friction_constant = 0;
#endif
                // penalty value for weak imposition of no penetration bc
                penetration_penalty = penetration_constant*pow(hE,penetration_power);;

                // parameter for friction
                //BoundaryValue(comp, (t0+t1)/2.0,friction_parameter);
//                 switch(friction_type)
//                 {
//                   case 1:
//                     // linear friction
//                     friction_parameter = friction_constant * pow(hE,friction_power);
//                     //OutPut("fric " <<  friction_parameter << endl);
//                     break;
//                   case 2:
//                     // nonlinear type 1
//                     // beta = sqrt(Re)(w\cdot tau)/(U_0/2-(w\cdot tau)
//                     // centre of the boundary edge
//                     x = (x1+x0)/2.0;
//                     y = (y1+y0)/2.0;
//                     // compute velocity in (x,y);
//                     u1->FindGradientLocal(cell,i,x,y,u1_values);
//                     u2->FindGradientLocal(cell,i,x,y,u2_values);
//                     // compute tangential velocity
//                     tangential_velo = u1_values[0]*tx + u2_values[0]*ty;
//                     // get Reynolds number
//                     RE_NR = TDatabase::ParamDB->RE_NR;
//                     U0 = TDatabase::ParamDB->FRICTION_U0;
//                     denominator = U0/2-tangential_velo;
//                     if (fabs(denominator)<1e-8)
//                     {
//                       OutPut("nonlinear slip bc type 1, denominator zero !!!" << endl);
//                       exit(4711);
//                     }
//                     friction_parameter = sqrt(RE_NR) * tangential_velo/denominator;
//                     friction_parameter = fabs(friction_parameter);
//                     OutPut("x " << x << " y " << y);
//                     OutPut(" friction paramater " << friction_parameter << endl);
//                     break;
//                 }

                hE = hE/2;

                EdgeDOF = FEDesc_Obj->GetJointDOF(m);

                // compute additional matrix entries
                // for all velo dof in the mesh cell
                // ii - test function
				
                for (ii=0;ii<N_;ii++)
                {
                  // look for 'ii'-th row in all matrices
                  dof_ii = DOF[ii];
                  // Dirichlet node
                  if (dof_ii>=ActiveBound)
                    continue;

                  // !!!!!! assumed that A_11 - A_22 are in sqmatrices[0] - [2]
                  // first velocity component -> matrices A_11 and A_12 (and M_11)
                  if (j==0)
                  {
                    Entries1 = sqmatrices[0]->GetEntries();
                    RowPtr1 = sqmatrices[0]->GetRowPtr();
                    ColInd1 = sqmatrices[0]->GetKCol();

                    if (n_sqmatrices>2)
                    {
                      Entries2 = sqmatrices[2]->GetEntries();
                      RowPtr2 = sqmatrices[2]->GetRowPtr();
                      ColInd2 = sqmatrices[2]->GetKCol();
                    }

                    // time dependent problem and NSTYPE 4
                    if (n_sqmatrices==8)
                    {
                      Entries3 = sqmatrices[4]->GetEntries();
                      RowPtr3 = sqmatrices[4]->GetRowPtr();
                      ColInd3 = sqmatrices[4]->GetKCol();
                      Entries4 = sqmatrices[6]->GetEntries();
                      RowPtr4 = sqmatrices[6]->GetRowPtr();
                      ColInd4 = sqmatrices[6]->GetKCol();
                    }

                    if (n_matrices>=2)
                    {
                      Entries5 = matrices[0]->GetEntries();
                      RowPtr5 = matrices[0]->GetRowPtr();
                      ColInd5 = matrices[0]->GetKCol();
                    }
                      if (n_matrices>=6)
                    {
                      Entries6 = matrices[2]->GetEntries();
                      RowPtr6 = matrices[2]->GetRowPtr();
                      ColInd6 = matrices[2]->GetKCol();
		      
		      Entries7 = matrices[3]->GetEntries();
                      RowPtr7 = matrices[3]->GetRowPtr();
                      ColInd7 = matrices[3]->GetKCol();
                    }
                    
                      if (n_matrices==10)
                    {
                      Entries8 = matrices[6]->GetEntries();
                      RowPtr8 = matrices[6]->GetRowPtr();
                      ColInd8 = matrices[6]->GetKCol();
		      
		      Entries9 = matrices[7]->GetEntries();
                      RowPtr9 = matrices[7]->GetRowPtr();
                      ColInd9 = matrices[7]->GetKCol();
                    }
                    
                    
                  }
                  // second velocity component -> matrices A_21 and A_22
                  if (j==1)
                  {
                    if (n_sqmatrices>2)
                    {
                      Entries1 = sqmatrices[3]->GetEntries();
                      RowPtr1 = sqmatrices[3]->GetRowPtr();
                      ColInd1 = sqmatrices[3]->GetKCol();
                    }

                    Entries2 = sqmatrices[1]->GetEntries();
                    RowPtr2 = sqmatrices[1]->GetRowPtr();
                    ColInd2 = sqmatrices[1]->GetKCol();

                    // time dependent problem and NSTYPE 4
                    if (n_sqmatrices==8)
                    {
                      Entries3 = sqmatrices[5]->GetEntries();
                      RowPtr3 = sqmatrices[5]->GetRowPtr();
                      ColInd3 = sqmatrices[5]->GetKCol();
                      Entries4 = sqmatrices[7]->GetEntries();
                      RowPtr4 = sqmatrices[7]->GetRowPtr();
                      ColInd4 = sqmatrices[7]->GetKCol();
                    }
                    if (n_matrices>=2)
                    {
                      Entries5 = matrices[1]->GetEntries();
                      RowPtr5 = matrices[1]->GetRowPtr();
                      ColInd5 = matrices[1]->GetKCol();
                    }
                    
                     if (n_matrices>=6)
                    {
                      Entries6 = matrices[4]->GetEntries();
                      RowPtr6 = matrices[4]->GetRowPtr();
                      ColInd6 = matrices[4]->GetKCol();
		      
		      Entries7 = matrices[5]->GetEntries();
                      RowPtr7 = matrices[5]->GetRowPtr();
                      ColInd7 = matrices[5]->GetKCol();
                    }
                    
                     if (n_matrices==10)
                    {
                      Entries8 = matrices[8]->GetEntries();
                      RowPtr8 = matrices[8]->GetRowPtr();
                      ColInd8 = matrices[8]->GetKCol();
		      
		      Entries9 = matrices[9]->GetEntries();
                      RowPtr9 = matrices[9]->GetRowPtr();
                      ColInd9 = matrices[9]->GetKCol();
                    }
                    
                    
                  }

                  //OutPut("ii " << dof_ii << endl);
                  // for all dof in the mesh cell
                  // jj - ansatz function
                  for (jj=0;jj<N_;jj++)
                  {
                    dof_jj = DOF[jj];
                    // OutPut("jj " << dof_jj << endl);
                    // initialize the boundary integrals
                    for (l=0;l<2;l++)
                      integral[l] = 0;

                    // compute boundary integrals
                    // first component of the velocity
                    if (j==0)
                    {
                      for(l=0;l<N_LinePoints;l++)
                      {
                        // values of test functions in this quadrature point
                        JointValue = JointValues[l];
                        // get quadrature point on the boundary
                        x = x0 + 0.5*(x1-x0)*(zeta[l]+1);
                        y = y0 + 0.5*(y1-y0)*(zeta[l]+1);

                        if(!axial3D) 
                         { r_axial = 1.;}
                        else if(axial3D && fabs(x)<1e-12)
                         { r_axial = 0.0; }
                        else
                         { r_axial = fabs(x); }


               // added varient for impinging droplets - Sashi
                
                if(fabs(x)>1e-12)
                 { beta_h = fabs(x); }
                else
                 { beta_h = 1.;}
                 
                switch(friction_type)
                {
                 case 1:
                  friction_parameter = friction_constant * pow(hE,friction_power);
                  break;

                 case 2:  
                  friction_parameter = friction_constant*pow(hE,friction_power)/(beta_h);
                  break;

                 case 3:
                  friction_parameter = friction_constant * pow(hE,friction_power)/(beta_h*beta_h);
                  break;

                 case 4:
                  friction_parameter = friction_constant *pow(hE,friction_power)/(beta_h*beta_h*beta_h);
                  break;

                 case 12: // Volker part
                    // nonlinear type 1
                    // beta = sqrt(Re)(w\cdot tau)/(U_0/2-(w\cdot tau)
                    // centre of the boundary edge
                    x = (x1+x0)/2.0;
                    y = (y1+y0)/2.0;
                    // compute velocity in (x,y);
                    u1->FindGradientLocal(cell,i,x,y,u1_values);
                    u2->FindGradientLocal(cell,i,x,y,u2_values);
                    // compute tangential velocity
                    tangential_velo = u1_values[0]*tx + u2_values[0]*ty;
                    // get Reynolds number
                    RE_NR = TDatabase::ParamDB->RE_NR;
                    U0 = TDatabase::ParamDB->FRICTION_U0;
                    denominator = U0/2-tangential_velo;
                    if (fabs(denominator)<1e-8)
                    {
                      OutPut("nonlinear slip bc type 1, denominator zero !!!" << endl);
                      exit(4711);
                    }
                    friction_parameter = sqrt(RE_NR) * tangential_velo/denominator;
                    friction_parameter = fabs(friction_parameter);
                    OutPut("x " << x << " y " << y);
                    OutPut(" friction paramater " << friction_parameter << endl);
                    break;
    
                default:
                 Error("Slip with friction Type not maching !!!!!!!!!!" << endl);
                 Error("file: " << __FILE__ << " line " << __LINE__ << endl);
                 exit(-1);  

                }
                
                // added by Jagannath - linearly varying slip with friction
//                 if(TDatabase::ParamDB->TENSOR_TYPE==1)
// 		{
// 		  if(TDatabase::ParamDB->FREE_SURFACE_FLOW==1)
// 		  {
// 		    if(friction_type==1)
// 		    {
// 		      if(x>0.9*TDatabase::ParamDB->P7)
// 		      friction_parameter = TDatabase::ParamDB->FRICTION_CONSTANT + (x - 0.9*TDatabase::ParamDB->P7) * (TDatabase::ParamDB->P8 - TDatabase::ParamDB->FRICTION_CONSTANT)/(0.1*TDatabase::ParamDB->P7);
// 		      else
// 		      friction_parameter = TDatabase::ParamDB->FRICTION_CONSTANT;
// 		    }
// 		    else
// 		    {
// 		     cout<<"Not yet implemented for other friction types \n";
// 		     exit(4711);
// 		    }
// 		  }
// 		}
		  
                
                
                       //t = t0 + 0.5*(t1-t0)*(zeta[l]+1);

                        // weight times determinant of reference trafo
                        // r_axial = 0, on axial BD so integrals are always zero on Axial BD
                        s = hE * LineWeights[l]*r_axial;
                        // (A_11)_{ii,jj}
//                         OutPut ("before " << integral[0] << " " << integral[1] << endl);
                        val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nx;
                        val += friction_parameter*JointValue[jj]*tx*JointValue[ii]*tx;
                        integral[0] += val*s;
// 			Frict_Force1 += val*s;

                        // (A_12)_{ii,jj}
                        val =  penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nx;
                        val+= friction_parameter*JointValue[jj]*ty*JointValue[ii]*tx;
                        integral[1] += val*s;
// 			Frict_Force2 += val*s;
                        //OutPut("    s   " <<s << endl);
                        //OutPut("       " << penetration_penalty*JointValue[jj]*nx*JointValue[ii]*nx<< endl);
//                         OutPut("       " << friction_parameter*JointValue[jj]*tx*JointValue[ii]*tx << endl);
                        //OutPut("       " << penetration_penalty*JointValue[jj]*ny*JointValue[ii]*nx << endl);
//                         OutPut("       " << friction_parameter*JointValue[jj]*ty*JointValue[ii]*tx << endl);
                      }
                      //if ((integral[0]==0)&&( integral[1]==0)) continue;
                      //OutPut ("1 " << integral[0] << " " << integral[1] << endl);

                      // edge not parallel to y axis or penetration
//                       if ((fabs(ny)>eps)||(penetration_penalty<1e3))
                      if ((fabs(tx)>eps)||(penetration_penalty>0))
                       {
//                         if(comp==0 || comp==1 || comp==5 )
//                          cout << j <<  " X " << x << " Y " << y <<endl;
                        // update first matrix
                        found = 0;
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                        {
                          if (ColInd1[ll] == dof_jj)
                          {
                            //OutPut("noy0 " << integral[0] << " ");
                            Entries1[ll] += integral[0];
                            found = 1;
                            break;
                          }
                        }
                        if (!found)
                        {
                          OutPut("ERROR A_11 " << endl);
                          exit(4711);
                        }
                        // update second matrix
                        if (n_sqmatrices>2)
                        {
                          found = 0;
                          for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                          {
                            if (ColInd2[ll] == dof_jj)
                            {
                              found = 1;
                              //OutPut("noy1 " << integral[1] << " ");
                              Entries2[ll] += integral[1];
                              break;
                            }
                          }
                          if (!found)
                          {
                            OutPut("ERROR A_12 "<< endl);
                            exit(4711);
                          }
                        }
                      }
                      else   // edge parallel to y-axis and no penetration
                      {
                        found = 0;
                        for (ll=0;ll<N_EdgeDOF; ll++)
                        {
                          if (dof_ii==DOF[EdgeDOF[ll]])
                          {
                            found =1;
                            break;
                          }
                        }
                        if (!found)
                          continue;
                        // OutPut("y" << endl);
                        // update first matrix, set diagonal entry to 1
                        // all other entries to zero
                        for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                        {
                          if (ColInd1[ll] == dof_ii)
                            Entries1[ll] = 1;
                          else
                            Entries1[ll] = 0;
                        }
                        // update second matrix, set all entries to zero
                        if (n_sqmatrices>2)
                        {
                          for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                            Entries2[ll] = 0;
                        }

                        if (n_sqmatrices==8)
                        {                         // M_11, set off diagonal to zero
                          for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          {
                            if (ColInd3[ll] != dof_ii)
                              Entries3[ll] = 0;
                          }
                          // M_12, set row to zero
                          for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                            Entries4[ll] = 0;
                        }

                        if (n_matrices>=2)
                          for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                            Entries5[ll] = 0;
			  
			if (n_matrices>=6)
			{  
			  for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                            Entries6[ll] = 0; 
			  
			   for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                            Entries7[ll] = 0; 
			}
			
	                  if (n_matrices==10)
			{   
			  for (ll=RowPtr8[dof_ii];ll < RowPtr8[dof_ii+1]; ll++)
                            Entries8[ll] = 0; 
			  
			   for (ll=RowPtr9[dof_ii];ll < RowPtr9[dof_ii+1]; ll++)
                            Entries9[ll] = 0; 
			}

                        // set rhs to zero
                        RHS[dof_ii] = 0;
                      }
                    }                             // end first component (j==0)

                    // second component
                    if (j==1)
                    {
                      for(l=0;l<N_LinePoints;l++)
                      {
                        // values of test functions in this quadrature point
                        JointValue = JointValues[l];
                        // get quadrature point on the boundary
                        x = x0 + 0.5*(x1-x0)*(zeta[l]+1);
                        y = y0 + 0.5*(y1-y0)*(zeta[l]+1);
                        
                        if(!axial3D) 
                         { r_axial = 1.;}
                        else if(axial3D && fabs(x)<1e-12)
                         { r_axial = 0.0; }
                        else
                         { r_axial = fabs(x); }
                         
                        //t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                        // get velocity in this quadrature point

                        // weight times determinant of reference trafo
                        s = hE * LineWeights[l]*r_axial;
                        // (A_21)_{ii,jj}
                        val = penetration_penalty*JointValue[jj]*nx*JointValue[ii]*ny;
                        val += friction_parameter*JointValue[jj]*ty*JointValue[ii]*tx;
                        integral[0] += s*val;
	
                        // (A_22)_{ii,jj}
                        val = penetration_penalty*JointValue[jj]*ny*JointValue[ii]*ny;
                        val += friction_parameter*JointValue[jj]*ty*JointValue[ii]*ty;
                        integral[1] += s*val;
                      }
                      //if ((integral[0]==0)&&( integral[1]==0)) continue;
                      //OutPut ("2 " << integral[0] << " " << integral[1] << endl);

                      // edge not parallel to x-axis or pentration

//                      if ((fabs(nx) > eps)|| (penetration_penalty < 1e3))
                     if ((fabs(nx) > eps) || (fabs(penetration_penalty) > 0.0))
                      { 
                        // OutPut("nox" << endl);
                        if (n_sqmatrices>2)
                        {
                          // update first matrix
                          found = 0;
                          for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                          {
                            if (ColInd1[ll] == dof_jj)
                            {
                              Entries1[ll] += integral[0];
                              found =1 ;
                              break;
                            }
                          }
                          if (!found)
                          {
                            OutPut("ERROR A_21 "<< endl);
                            exit(4711);
                          }
                        }

                        // update second matrix
                        found = 0;
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                        {
                          if (ColInd2[ll] == dof_jj)
                          {
                            Entries2[ll] += integral[1];
                            found =1 ;
                            break;
                          }
                        }
                        if (!found)
                        {
                          OutPut("ERROR A_22 "<< endl);
                          exit(4711);
                        }
                      }
                      else                        // edge parallel to x-axis and no penetration
                      {
                        //OutPut("x ");
                        found = 0;
                        for (ll=0;ll<N_EdgeDOF; ll++)
                        {
                          if (dof_ii==DOF[EdgeDOF[ll]])
                          {
                            found =1;
                            break;
                          }
                        }
                        if (!found)
                          continue;

                        // update first matrix, set all entries to zero
                        if (n_sqmatrices>2)
                        {                            
                          for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                            Entries1[ll] = 0;
                        }
                        // update second matrix, set diagonal entry to 1,
                        // all other entries to 0
                        for (ll=RowPtr2[dof_ii];ll < RowPtr2[dof_ii+1]; ll++)
                        {
                          if (ColInd2[ll] == dof_ii)
                          { Entries2[ll] = 1; }
                          else
                           {  Entries2[ll] = 0; }
                        }

                        // set rhs to zero
                        RHS[dof_ii] = 0;

                        // update mass matrix
                        if (n_sqmatrices==8)
                        {
                          for (ll=RowPtr3[dof_ii];ll < RowPtr3[dof_ii+1]; ll++)
                          {                       //M_22
                            if (ColInd3[ll] != dof_ii)
                              Entries3[ll] = 0;
                          }
                          // M_21
                          for (ll=RowPtr4[dof_ii];ll < RowPtr4[dof_ii+1]; ll++)
                            Entries4[ll] = 0;
                        }

                        if (n_matrices>=2)
                          for (ll=RowPtr5[dof_ii];ll < RowPtr5[dof_ii+1]; ll++)
                            Entries5[ll] = 0;
			  
			  if (n_matrices>=6)
			{   
			  for (ll=RowPtr6[dof_ii];ll < RowPtr6[dof_ii+1]; ll++)
                            Entries6[ll] = 0; 
			  
			   for (ll=RowPtr7[dof_ii];ll < RowPtr7[dof_ii+1]; ll++)
                            Entries7[ll] = 0; 
			}
			
		          if (n_matrices==10)
			{   
			  for (ll=RowPtr8[dof_ii];ll < RowPtr8[dof_ii+1]; ll++)
                            Entries8[ll] = 0; 
			  
			   for (ll=RowPtr9[dof_ii];ll < RowPtr9[dof_ii+1]; ll++)
                            Entries9[ll] = 0; 
			}

                      }
                    }                             // end first component (j==1)
                  }                               // end inner loop over dof (jj)
                }                                 // end outer loop over dof (ii)

                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                break;                            // end slip with friction and penetration with resistance bc
            default:
 
            break;   
            }                                     // endswitch Cond0     
          }                                       // endif (Cond0==Cond1)
          else
          {
            OutPut("different boundary condition on one edge ");
            OutPut("are not allowed!" << endl);
            exit(4711);
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)

  
//   OutPut(" Frict_Force1 " << Frict_Force1  << " " << Frict_Force2  <<endl);
  
  if(n_sqmatrices)
  {
    delete [] GlobalNumbers;
    delete [] BeginIndex;

    for(i=0;i<n_sqmatrices;i++)
      delete [] HangingEntries[i];
    delete [] HangingEntries;
  }

  if(n_matrices)
  {
    delete [] AnsatzGlobalNumbers;
    delete [] AnsatzBeginIndex;
    delete [] TestGlobalNumbers;
    delete [] TestBeginIndex;
  }

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete [] HangingRhs[i];
    delete [] HangingRhs;

    delete [] righthand;
    delete [] LocRhs;
    delete [] RhsBeginIndex;
    delete [] RhsGlobalNumbers;
  }

  if(N_Parameters)
  {
    delete [] Param[0];
  }

  if(N_AllMatrices)
  {
    delete [] LocMatrices;
    delete [] Matrices[0];
    delete [] Matrices;
  }

  delete [] AuxArray[0];
  delete [] SecondDer;

  /*
    int N_Rows;
    // ####################################################################
    // print the whole matrix -- SECOND
    // ####################################################################
    for(k=0;k<n_sqmatrices;k++)
    {
      cout << endl;
      cout << "sqmatrix: " << k << endl;
      RowPtr = sqmatrices[k]->GetRowPtr();
      Entries = sqmatrices[k]->GetEntries();
  ColInd = sqmatrices[k]->GetKCol();
  N_Rows = sqmatrices[k]->GetN_Rows();
  for(i=0;i<N_Rows;i++)
  {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_matrices;k++)
  {
  cout << endl;
  cout << "matrix: " << k << endl;
  RowPtr = matrices[k]->GetRowPtr();
  Entries = matrices[k]->GetEntries();
  ColInd = matrices[k]->GetKCol();
  N_Rows = matrices[k]->GetN_Rows();
  for(i=0;i<N_Rows;i++)
  {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_rhs;k++)
  {
  cout << "rhs: " << k << endl;
  N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
  RHS=rhs[k];
  for(i=0;i<N_Rows;i++)
  cout << setw(5) << i << setw(20) << RHS[i] << endl;
  }
  */

}  // end of Assemble


/*
  Assemble2D_VectFE:
    Assemble for vector finite elements (Raviart-Thomas, Brezzi-Douglas-Marini)
    Need the global orientation of normal at each inner edge/face

    implementation: Alfonso (07.09.2010)
*/

void Assemble2D_VectFE(int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
TDiscreteForm2D *DiscreteForm,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters
)
{
#ifdef __2D__
  if(Parameters)
  {
    Error("Assemble2D_VectFE: input 'Parameters' of type 'TAuxParam2D*' is "<<
          "not set to NULL. This is usually done if you want values of FE "<<
          "functions during local assembling, for example in nonlinear "
          "problems. This is not yet supported. Exiting.\n");
    exit(1);
  }
  
  int N_AllMatrices = n_sqmatrices+n_matrices;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *righthand;
  
  
  // check if hanging nodes exist, I don't know how this works, we quit the
  // program in that case.
  for (int iSpace=0; iSpace<n_fespaces; iSpace++)
  {
    if(fespaces[iSpace]->GetN_Hanging()!=0)
    {
      Error("Assemble2D_VectFE: hanging entries not supported. Exiting");
      exit(1);
    }
  }
  bool *SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
  
  // ########################################################################
  // loop over all cells
  // ########################################################################
  TCollection *Coll = fespaces[0]->GetCollection();// all spaces use same Coll
  int N_Cells = Coll->GetN_Cells(); // number of cells in this collection
  // set cell indices for all cells
  for(int icell=0; icell<N_Cells; icell++)
    Coll->GetCell(icell)->SetCellIndex(icell);
  for(int icell=0; icell<N_Cells; icell++)
  {
    double **LocRhs;
    if(n_rhs)
    {
      LocRhs = new double* [n_rhs];
      righthand = new double [n_rhs*MaxN_BaseFunctions2D];
      memset(righthand,0,SizeOfDouble*n_rhs*MaxN_BaseFunctions2D);
      for(int i=0;i<n_rhs;i++)
        LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
    }
    TBaseCell *cell = Coll->GetCell(icell); // current cell
    int N_Edges  = cell->GetN_Edges(); // number of edges of this cell
    
    // set orientation of normal at internal edges
    // this function is not accessable during the local assembling routine
    // so we use a global parameter as well
    cell->SetNormalOrientation();
    for (int ijo = 0; ijo<N_Edges; ijo++)
    {
      if (N_Edges==4)
      {                                       // Raviart-Thomas, quadrilaterals
        TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[ijo] = 
            cell->GetNormalOrientation(ijo);
        
      } 
      else if (N_Edges==3)
      {                                       // Raviart-Thomas, triangles
        TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[ijo] = 
            cell->GetNormalOrientation(ijo);
      }
    }                                             // for ijo
    double hK = 0; // size of cell
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
      case 0:                           // diameter
      case 4:                           // mesh size in convection direction
      case 5:                           // take value from an array
        hK = cell->GetDiameter();
        break;
      case 2:                           // shortest edge
        hK = cell->GetShortestEdge();
        break;
      case 1:                           // with reference map
      case 3:                           // measure
        hK = cell->GetMeasure();
        hK = sqrt(hK);
        break;
      default:                          // diameter
        hK = cell->GetDiameter();
        break;
    }
    // ########################################################################
    // find local used elements on this cell
    // ########################################################################
    std::vector<FE2D> LocalUsedElements(n_fespaces);
    std::vector<int> LocN_BF(n_fespaces);
    std::vector<BaseFunct2D> LocBF(n_fespaces);
    for(int iSpace=0; iSpace<n_fespaces; iSpace++)
    {
      FE2D CurrentElement = fespaces[iSpace]->GetFE2D(icell, cell);
      TFE2D *fe = TFEDatabase2D::GetFE2D(CurrentElement);
      LocalUsedElements[iSpace] = CurrentElement;
      LocN_BF[iSpace] = fe->GetSize();
      LocBF[iSpace] = fe->GetBaseFunct2D_ID();
    }
    
    int N_LocalUsedElements = n_fespaces;
    
    // ########################################################################
    // calculate values on original element
    // ########################################################################
    int N_Points; // number of quadrature points
    TFEDatabase2D::GetOrig(N_LocalUsedElements, &LocalUsedElements[0],
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);
    
    // this could provide values of FE functions during the local assemble
    // routine, not yet supported.
    //Parameters->GetParameters(N_Points,Coll, cell,icell, xi,eta, X,Y, Param);

    // ########################################################################
    // assemble local matrices and right hand sides
    // ########################################################################
    // maximum number of used basis functions
    //int max_n_BF = *max_element(LocN_BF.begin(),LocN_BF.end()); 
    // for every matrix we allocate a local matrix with corresponding number of 
    // rows and columns
    double ***LocMatrices = NULL;
    if(N_AllMatrices)
    {
      LocMatrices = new double**[N_AllMatrices];
      for(int i=0;i<N_AllMatrices;i++)
      {
        int n_rows =LocN_BF[DiscreteForm->rowSpaceOfMat(i)];//number of rows
        int n_cols =LocN_BF[DiscreteForm->colSpaceOfMat(i)];//number of columns
        LocMatrices[i] = new double*[n_rows];
        for(int j=0; j<n_rows; j++)
        {
          LocMatrices[i][j] = new double[n_cols]; 
          memset(LocMatrices[i][j], 0, SizeOfDouble*n_cols);
        }
      }
    }                                               // endif N_AllMatrices
    
    if(DiscreteForm)
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk, hK, X, Y, 
        &LocN_BF[0], &LocBF[0], cell, LocMatrices, LocRhs);
    
    
    // ########################################################################
    // add local matrices to global matrices (ansatz == test)
    // ########################################################################
    for(int iSqMat=0;iSqMat<n_sqmatrices;iSqMat++)
    {
      // fe space for this square matrix
      TFESpace2D *fespace = fespaces[DiscreteForm->rowSpaceOfMat(iSqMat)];
      // the number of local basis functions (= size of local matrix)
      int N_BaseFunctions = LocN_BF[DiscreteForm->rowSpaceOfMat(iSqMat)];
      
      double **Matrix = LocMatrices[iSqMat];
      int ActiveBound = fespace->GetActiveBound();
      int *DOF = fespace->GetGlobalDOF(icell);
      
      // add local matrix to global
      for(int irow = 0; irow < N_BaseFunctions; irow++)
      {
        int RowDOF = DOF[irow];
        if(RowDOF<ActiveBound)
        { // active degree of freedom
          for(int icolumn=0;icolumn<N_BaseFunctions;icolumn++)
          {
            int columnDOF=DOF[icolumn];
            sqmatrices[iSqMat]->add(RowDOF,columnDOF,Matrix[irow][icolumn]);
          }
        }
        else
        { // nonactive degree of freedom (Dirichlet)
          sqmatrices[iSqMat]->set(RowDOF,RowDOF,1.0); // 1 on diagonal
        }
      }                                           // endfor m
    }                                             // endfor j
    
    // ########################################################################
    // add local matrices to global matrices (ansatz != test)
    // ########################################################################
    for(int iMat=0;iMat<n_matrices;iMat++)
    {
      TFESpace2D* testSpace = matrices[iMat]->GetStructure()->GetTestSpace2D();
      TFESpace2D*ansatzSpace=matrices[iMat]->GetStructure()->GetAnsatzSpace2D();
      TFE2D*test_fe = TFEDatabase2D::GetFE2D(testSpace->GetFE2D(icell, cell));
      TFE2D*ansatz_fe=TFEDatabase2D::GetFE2D(ansatzSpace->GetFE2D(icell,cell));
      
      // number of test and ansatz functions
      int N_Test = test_fe->GetSize();
      int N_Ansatz = ansatz_fe->GetSize();
      
      double **Matrix = LocMatrices[iMat+n_sqmatrices];
      
      int *TestDOF = testSpace->GetGlobalDOF(icell);
      int *AnsatzDOF = ansatzSpace->GetGlobalDOF(icell);
      
      int ActiveBound = testSpace->GetActiveBound();
      
      // add local matrix to global
      for(int irow = 0; irow < N_Test; irow++)
      {
        int rowDOF = TestDOF[irow];
        if(rowDOF<ActiveBound)
        {
          for(int icolumn=0; icolumn<N_Ansatz; icolumn++)
          {
            int columnDOF = AnsatzDOF[icolumn];
            matrices[iMat]->add(rowDOF,columnDOF,Matrix[irow][icolumn]);
          }
        }
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)
    
    if(N_AllMatrices)
    {
      for(int i=0;i<N_AllMatrices;i++)
      {
        int n_rows =LocN_BF[DiscreteForm->rowSpaceOfMat(i)];//number of rows
        for(int j=0; j<n_rows; j++)
        {
          delete [] LocMatrices[i][j];
        }
        delete [] LocMatrices[i];
      }
      delete [] LocMatrices;
    }
    // ########################################################################
    // add local right-hand sides to global right-hand side
    // ########################################################################
    for(int irhs=0;irhs<n_rhs;irhs++)
    {
      TFESpace2D *fespace = ferhs[irhs];
      FE2D CurrentElement = fespace->GetFE2D(icell, cell);
      TFE2D *fe = TFEDatabase2D::GetFE2D(CurrentElement);

      int N_BaseFunctions = fe->GetSize();

      double *local_rhs = LocRhs[irhs];
      double *RHS = rhs[irhs];
      int ActiveBound = fespace->GetActiveBound();
      
      // dof of the rhs nodes connected to this cell
      int *DOF = fespace->GetGlobalDOF(icell);

      // add local right-hand side to the global one
      for(int irow = 0; irow < N_BaseFunctions; irow++)
      { 
        int rowDOF = DOF[irow];
        if(rowDOF<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[rowDOF] += local_rhs[irow];
        }                                         // endif l
      }                                           // endfor m
      //////////////////////////////////////////////////////////////////////
      // take care of boundary conditions:      
      BoundCondFunct2D *BoundaryCondition = BoundaryConditions[irhs];
      BoundValueFunct2D *BoundaryValue = BoundaryValues[irhs];
      
      for(int ijoint=0; ijoint<N_Edges; ijoint++)
      {
        TJoint *joint = cell->GetJoint(ijoint);
        if(joint->GetType() == BoundaryEdge)
        {
          TBoundEdge *boundedge = (TBoundEdge *)joint;
          double t0,t1;
          boundedge->GetParameters(t0, t1);
          // get id of the boundary component
          int comp = boundedge->GetBoundComp()->GetID();
          // get type of the boundary condition in the middle of the edge
          BoundCond Cond;
          BoundaryCondition(comp, (t0+t1)/2, Cond);
          switch(Cond)
          {
            case DIRICHLET:
            {
              TFEDesc2D *FEDesc_Obj = fe->GetFEDesc2D();
              int N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              // if DG
              if (N_EdgeDOF==0)
                break;
              TNodalFunctional2D *nf = fe->GetNodalFunctional2D();
              // number of points used for computation of nodal functionals
              int N_EdgePoints; 
              // points used for computation of nodal functionals
              double *EdgePoints; 
              nf->GetPointsForEdge(N_EdgePoints, EdgePoints);
              
              // read boundary values for each point
              double PointValues[N_EdgePoints];
              for(int iedgePoint=0;iedgePoint<N_EdgePoints;iedgePoint++)
              {
                double s = EdgePoints[iedgePoint];
                s = 0.5*(t0*(1-s) + t1*(1+s)); // map s from [-1,1] to [t0,t1]
                BoundaryValue(comp, s, PointValues[iedgePoint]);
              }
              // compute values for each dof on the boundary edge with the 
              // nodal functionals
              double FunctionalValues[N_EdgeDOF];
              nf->GetEdgeFunctionals(Coll, cell, ijoint, PointValues,
                FunctionalValues);
              int *EdgeDOF = FEDesc_Obj->GetJointDOF(ijoint);
              // save boundary values of each dof on the boundary
              // edge in the rhs
              for(int l=0;l<N_EdgeDOF;l++)
              {
                RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
              }
              break;
            }
            case NEUMANN:
            {
              // Basis functions
              TBaseFunct2D* BaseFunct = fe->GetBaseFunct2D();
              // get polynomial degree of fe
              int polynomialDegree = BaseFunct->GetPolynomialDegree();
              // get a suitable line quadrature formula
              QuadFormula1D LineQuadFormula = 
                     TFEDatabase2D::GetQFLineFromDegree(2*polynomialDegree);
              TQuadFormula1D *qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
              int N_LinePoints; // number of quadrature points on this edge
              double *LineWeights, *zeta;// quadrature points and weights
              qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
              BaseFunct->MakeRefElementData(LineQuadFormula);
              double **JointValues=TFEDatabase2D::GetJointValues2D(
                BaseFunct->GetID(), LineQuadFormula, ijoint);
              BaseFunct->ChangeBF(Coll, cell, N_LinePoints, JointValues);
              // get vertices of boundary edge
              double x0, x1, y0, y1;
              cell->GetVertex(ijoint)->GetCoords(x0, y0);
              cell->GetVertex((ijoint+1) % N_Edges)->GetCoords(x1, y1);
              // compute (half of the) length of the boundary edge
              double hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
              double nx = (y1-y0)/(2*hE);
              double ny = (x0-x1)/(2*hE);
              // for Piola transform
              double *JointValuesTransformed = 
                    new double [N_BaseFunctions * BaseFunct->GetBaseVectDim()];
              // compute boundary integral
              for(int ilinePoint=0; ilinePoint<N_LinePoints; ilinePoint++)
              {
                // values of test functions in this quadrature point
                double *JointValue = JointValues[ilinePoint];
                TFEDatabase2D::GetOrigValues(
                    fe->GetRefTransID(), zeta[ilinePoint], BaseFunct, ijoint,
                    JointValue, NULL,NULL, JointValuesTransformed, NULL,NULL); 
                // get quadrature point on the boundary
                double t = t0 + 0.5*(t1-t0)*(zeta[ilinePoint]+1);
                double s;
                // get value in this quadrature point (in s)
                BoundaryValue(comp, t, s);
                // multiply value with weights from quadrature formula
                // and determinant from integral transformation to the
                // unit edge (-1,1)
                s *= hE * LineWeights[ilinePoint];
                // in case of the pressure right hand side, the values in 
                // JointValuesTransformed[k+N_] are not valid. Therefore we 
                // need HOMOGENEOUS Neumann boundary conditions for the 
                // pressure function
                if(s==0.0)
                  continue;
                
                // update rhs for all test functions
                for(int k=0;k<N_BaseFunctions;k++)
                {
                  int rowDOF = DOF[k];
                  if(rowDOF < ActiveBound)
                  {
                    RHS[rowDOF] -= s*(JointValuesTransformed[k]*nx 
                                 +JointValuesTransformed[k+N_BaseFunctions]*ny);
                  }
                }
              }
              
              delete [] JointValuesTransformed;
              break;
            }
            case ROBIN:
              /*// get polynomial degree of fe
              l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                (CurrentElement);
              // get a suitable line quadrature formula
              LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
              qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
              qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
              TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                ->MakeRefElementData(LineQuadFormula);
              JointValues=TFEDatabase2D::GetJointValues2D(
                BaseFuncts[CurrentElement], LineQuadFormula, m);
              // get vertices of boundary edge
              cell->GetVertex(m)->GetCoords(x0, y0);
              cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
              // compute (half of the) length of the boundary edge
              hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
              // compute boundary integral
              for(l=0;l<N_LinePoints;l++)
              {
                // values of test functions in this quadrature point
                JointValue = JointValues[l];
                // get quadrature point on the boundary
                t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                // get value in this quadrature point (in s)
                BoundaryValue(comp, t, s);
                // multiply value with weights from quadrature formula
                // and determinant from integral transformation to the
                // unit edge (-1,1)
                s *= hE * LineWeights[l];
                // update rhs for all test functions
                for(k=0;k<N_BaseFunctions;k++)
                  if((l3 = DOF[k])<ActiveBound)
                    RHS[l3] += s*JointValue[k];
              }*/
              Error("Robin boundary conditions not yet supported. Exiting\n");
              exit(1);
              break;
            default:
              OutPut("Unknown boundary condition !"<< endl);
              exit(4711);
              break;
          }                                     // endswitch Cond0
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
    if(n_rhs)
    {
      delete [] LocRhs[0]; 
      delete [] LocRhs;
    }
  } // endfor i (N_Cells)
#endif // __2D__
} // end of Assemble_VectFE

#ifdef __MORTAR__
// compute zeta \in (-1, 1) as: v2 = v0 + 0.5(zeta + 1) * (v1 - v0)
double GetZeta(TVertex *v0, TVertex *v1, TVertex *v2)
{
  double T, delX, delY;

  delX = v1->GetX() - v0->GetX();
  if (ABS(delX) > 1e-8)
    T = (v2->GetX() - v0->GetX()) / delX;
  else
  {
    delY = v1->GetY() - v0->GetY();
    if (ABS(delY) > 1e-8)
      T = (v2->GetY() - v0->GetY()) / delY;
  }

  return 2*T - 1;
}


double GetLambda(double sX, double sY, TVertex *vert, double dX, double dY)
{
  if (ABS(dX) > 1e-8)
    return (vert->GetX() - sX) / dX;
  else
    return (vert->GetY() - sY) / dY;
}


#include <MortarBaseJoint.h>
#include <It_Mortar.h>

// assemble additional link term
#ifdef __ADD_LINK_SDFEM__
void AddLinkSDFEM(TFESpace2D *Space2D, TSquareMatrix2D *Matrix,
TDiscreteForm2D *DiscreteForm)
{
  // additional SDFEM link term
  int i, j, k, l, N_;
  TIt_Mortar *It1 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar1];
  TIt_Mortar *It2 = (TIt_Mortar *) TDatabase::IteratorDB[It_Mortar2];
  TBaseCell *CellM, *CellNM;
  FE2D CurrElementID;
  TFE2D *CurrElement;
  TFEDesc2D *CurrDesc;
  int lev, infoNM, infoM;
  double lamM0 = 0, lamM1 = 0, lamNM0, lamNM1, startX, startY, endX, endY;
  double lam, lam0, lam1;
  double delX, delY;
  bool NewMortarSideEle;
  int **J_DOF_NM, N_J_DOF_NM, **J_DOF_M, N_J_DOF_M;
  int indexM = 0, indexNM;
  double *coeff = new double[6], *param = new double[1];
  double nx, ny, nlen, val;
  int ll, kk, *NumbersM, *NumbersNM, index;
  int *GlobalNumbers = Space2D->GetGlobalNumbers();
  int *BeginIndex = Space2D->GetBeginIndex();
  QuadFormula1D QuadFormula;
  TQuadFormula1D *Line_Formula;
  double *Entries = Matrix->GetEntries();
  int *RowPtr = Matrix->GetRowPtr(), *KCol = Matrix->GetKCol();
  int NeumannBound = Space2D->GetN_Inner() + Space2D->GetN_BoundaryNodes()[0];
  int PolynomialDegree, N_QuadPoints_Line;
  double *Line_Weights, *Line_Xi, **Ele_NM_Values, **Ele_M_Values, *Flux;
  double zeta0NM, zeta1NM, zeta0M, zeta1M, X, Y;
  double zeta, xi, eta;
  TVertex *AuxVert0 = new TVertex(0, 0), *AuxVert1 = new TVertex(0, 0);
  const int *TmpEV_NM, *TmpEV_M;
  TBaseFunct2D *Ele_NM_BF, *Ele_M_BF;
  int LocDofNM, LocDofM, LocDofNM1, LocDofM1;
  double delX_int, delY_int, len_int;

  CurrElementID = Space2D->GetUsedElements()[0];
  PolynomialDegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrElementID);
  QuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*PolynomialDegree);

  Line_Formula = TFEDatabase2D::GetQuadFormula1D(QuadFormula);
  Line_Formula->GetFormulaData(N_QuadPoints_Line, Line_Weights, Line_Xi);

  Ele_NM_Values = new double*[N_QuadPoints_Line];
  Ele_M_Values = new double*[N_QuadPoints_Line];
  Flux = new double[N_QuadPoints_Line];

  *Ele_NM_Values = new double[MaxN_BaseFunctions2D * N_QuadPoints_Line];
  for (i=1;i<N_QuadPoints_Line;i++)
    Ele_NM_Values[i] = *Ele_NM_Values + i*MaxN_BaseFunctions2D;
  *Ele_M_Values = new double[MaxN_BaseFunctions2D * N_QuadPoints_Line];
  for (i=1;i<N_QuadPoints_Line;i++)
    Ele_M_Values[i] = *Ele_M_Values + i*MaxN_BaseFunctions2D;

  // assemble parts which connect both siedes of a mortar
  // loop over all mortar edges
  N_ = It1->GetN_MortarFace();
  for (i=0;i<N_;i++)
  {
    NewMortarSideEle = TRUE;
    lev = MAX_ItLevel + (i << 8);
    It1->Init(lev);
    It2->Init(-lev);

    It1->GetPoint(startX, startY);
    It2->GetPoint(endX, endY);
    delX = endX - startX;
    delY = endY - startY;

    // Assumption: n = const. (NM -> M) on each mortar edge !!!!!!
    nlen = sqrt(delX*delX + delY*delY);
    nx = delY / nlen;
    ny = -delX / nlen;

    while (CellNM = It1->Next(infoNM))
    {
      CellNM->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV_NM);
      // get a new element on non-mortar side
      lamNM0 = GetLambda(startX, startY, CellNM->GetVertex(infoNM),
        delX, delY);
      lamNM1 = GetLambda(startX, startY, CellNM->GetVertex((infoNM+1)
        % CellNM->GetN_Vertices()), delX, delY);
      if (CellNM->GetClipBoard() == -1)
      {
        cerr << "Error in MatrixStructure: cell out of collection!!!"
          << endl;
        exit(-3);
      }

      indexNM = CellNM->GetClipBoard();
      CurrElementID = Space2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      Ele_NM_BF = CurrElement->GetBaseFunct2D();
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      if (!NewMortarSideEle)
      {
        // assemble link term if the next while clause wont be entered
        NumbersM = GlobalNumbers + BeginIndex[indexM];
        NumbersNM = GlobalNumbers + BeginIndex[indexNM];

        lam0 = MAX(lamM0, lamNM0);
        lam1 = MIN(lamM1, lamNM1);

        AuxVert0->SetCoords(startX + lam0*delX, startY + lam0*delY);
        AuxVert1->SetCoords(startX + lam1*delX, startY + lam1*delY);

        delX_int = (lam1 - lam0)*delX;
        delY_int = (lam1 - lam0)*delY;
        len_int = sqrt(delX_int*delX_int + delY_int*delY_int);

        zeta0NM = GetZeta(CellNM->GetVertex(TmpEV_NM[2*infoNM]),
          CellNM->GetVertex(TmpEV_NM[2*infoNM + 1]),
          AuxVert0);
        zeta1NM = GetZeta(CellNM->GetVertex(TmpEV_NM[2*infoNM]),
          CellNM->GetVertex(TmpEV_NM[2*infoNM + 1]),
          AuxVert1);
        zeta0M = GetZeta(CellM->GetVertex(TmpEV_M[2*infoM]),
          CellM->GetVertex(TmpEV_M[2*infoM + 1]),
          AuxVert0);
        zeta1M = GetZeta(CellM->GetVertex(TmpEV_M[2*infoM]),
          CellM->GetVertex(TmpEV_M[2*infoM + 1]),
          AuxVert1);

        for (k=0;k<N_QuadPoints_Line;k++)
        {
          lam = lam0 + 0.5*(Line_Xi[k] + 1)*(lam1 - lam0);
          X = startX + lam*delX;
          Y = startY + lam*delY;
          DiscreteForm->GetCoeffFct()(1, &X, &Y, &param, &coeff);
          Flux[k] = coeff[1]*nx + coeff[2]*ny;
        }

        switch (CellNM->GetType())
        {
          case Triangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0NM + 0.5*(Line_Xi[k] + 1)*(zeta1NM - zeta0NM);
              zeta = 0.5*(zeta + 1);
              switch(infoNM)
              {
                case 0: xi=zeta;
                eta=0;
                break;
                case 1: xi=1-zeta;
                eta=zeta;
                break;
                case 2: xi=0;
                eta=1-zeta;
                break;
              }
              Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
            }
            break;

          case Rectangle:
          case Parallelogram:
          case Quadrangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0NM + 0.5*(Line_Xi[k] + 1)*(zeta1NM - zeta0NM);
              switch(infoNM)
              {
                case 0: xi=zeta;
                eta=-1;
                break;
                case 1: xi=1;
                eta=zeta;
                break;
                case 2: xi=-zeta;
                eta=1;
                break;
                case 3: xi=-1;
                eta=-zeta;
                break;
              }
              Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
            }
            break;
        }

        switch (CellM->GetType())
        {
          case Triangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0M + 0.5*(Line_Xi[k] + 1)*(zeta1M - zeta0M);
              zeta = 0.5*(zeta + 1);
              switch(infoM)
              {
                case 0: xi=zeta;
                eta=0;
                break;
                case 1: xi=1-zeta;
                eta=zeta;
                break;
                case 2: xi=0;
                eta=1-zeta;
                break;
              }
              Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
            }
            break;

          case Rectangle:
          case Parallelogram:
          case Quadrangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0M + 0.5*(Line_Xi[k] + 1)*(zeta1M - zeta0M);
              switch(infoM)
              {
                case 0: xi=zeta;
                eta=-1;
                break;
                case 1: xi=1;
                eta=zeta;
                break;
                case 2: xi=-zeta;
                eta=1;
                break;
                case 3: xi=-1;
                eta=-zeta;
                break;
              }
              Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
            }
            break;
        }

        for (j=0;j<N_J_DOF_NM;j++)
          for (k=0;k<N_J_DOF_M;k++)
        {
          LocDofNM = J_DOF_NM[infoNM][j];
          kk = NumbersNM[LocDofNM];
          LocDofM = J_DOF_M[infoM][k];
          ll = NumbersM[LocDofM];

          if (kk < NeumannBound || ll < NeumannBound)
          {
            for (val=l=0;l<N_QuadPoints_Line;l++)
              val += Line_Weights[l]*Flux[l] *
                Ele_NM_Values[l][LocDofNM] * Ele_M_Values[l][LocDofM];

            val *= len_int;

            if (kk < NeumannBound)
            {
              index = RowPtr[kk];
              l = KCol[index];
              while (l != ll)
                l = KCol[++index];
              Entries[index] -= 0.5*val;
            }
            if (ll < NeumannBound)
            {
              index = RowPtr[ll];
              l = KCol[index];
              while (l != kk)
                l = KCol[++index];
              Entries[index] += 0.5*val;
            }
          }
        }
      }

      // check which side is next
      if (lamM1 < lamNM1) NewMortarSideEle = TRUE;

      while (NewMortarSideEle)
      {
        // get a new element on mortar side
        if (!(CellM = It2->Next(infoM))) break;
        CellM->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV_M);
        lamM0 = GetLambda(startX, startY, CellM->GetVertex((infoM+1)
          % CellM->GetN_Vertices()), delX, delY);
        lamM1 = GetLambda(startX, startY, CellM->GetVertex(infoM),
          delX, delY);
        if (CellM->GetClipBoard() == -1)
        {
          cerr << "Error in MatrixStructure: cell out of collection!!!"
            << endl;
          exit(-3);
        }

        indexM = CellM->GetClipBoard();
        CurrElementID = Space2D->GetFE2D(indexM, CellM);
        CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
        Ele_M_BF = CurrElement->GetBaseFunct2D();
        CurrDesc = CurrElement->GetFEDesc2D();

        J_DOF_M = CurrDesc->GetJointDOF();
        N_J_DOF_M = CurrDesc->GetN_JointDOF();

        // assemble link term
        NumbersM = GlobalNumbers + BeginIndex[indexM];
        NumbersNM = GlobalNumbers + BeginIndex[indexNM];

        lam0 = MAX(lamM0, lamNM0);
        lam1 = MIN(lamM1, lamNM1);

        AuxVert0->SetCoords(startX + lam0*delX, startY + lam0*delY);
        AuxVert1->SetCoords(startX + lam1*delX, startY + lam1*delY);

        delX_int = (lam1 - lam0)*delX;
        delY_int = (lam1 - lam0)*delY;
        len_int = sqrt(delX_int*delX_int + delY_int*delY_int);

        zeta0NM = GetZeta(CellNM->GetVertex(TmpEV_NM[2*infoNM]),
          CellNM->GetVertex(TmpEV_NM[2*infoNM + 1]),
          AuxVert0);
        zeta1NM = GetZeta(CellNM->GetVertex(TmpEV_NM[2*infoNM]),
          CellNM->GetVertex(TmpEV_NM[2*infoNM + 1]),
          AuxVert1);
        zeta0M = GetZeta(CellM->GetVertex(TmpEV_M[2*infoM]),
          CellM->GetVertex(TmpEV_M[2*infoM + 1]),
          AuxVert0);
        zeta1M = GetZeta(CellM->GetVertex(TmpEV_M[2*infoM]),
          CellM->GetVertex(TmpEV_M[2*infoM + 1]),
          AuxVert1);

        for (k=0;k<N_QuadPoints_Line;k++)
        {
          lam = lam0 + 0.5*(Line_Xi[k] + 1)*(lam1 - lam0);
          X = startX + lam*delX;
          Y = startY + lam*delY;
          DiscreteForm->GetCoeffFct()(1, &X, &Y, &param, &coeff);
          Flux[k] = coeff[1]*nx + coeff[2]*ny;
        }

        switch (CellNM->GetType())
        {
          case Triangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0NM + 0.5*(Line_Xi[k] + 1)*(zeta1NM - zeta0NM);
              zeta = 0.5*(zeta + 1);
              switch(infoNM)
              {
                case 0: xi=zeta;
                eta=0;
                break;
                case 1: xi=1-zeta;
                eta=zeta;
                break;
                case 2: xi=0;
                eta=1-zeta;
                break;
              }
              Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
            }
            break;

          case Rectangle:
          case Parallelogram:
          case Quadrangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0NM + 0.5*(Line_Xi[k] + 1)*(zeta1NM - zeta0NM);
              switch(infoNM)
              {
                case 0: xi=zeta;
                eta=-1;
                break;
                case 1: xi=1;
                eta=zeta;
                break;
                case 2: xi=-zeta;
                eta=1;
                break;
                case 3: xi=-1;
                eta=-zeta;
                break;
              }
              Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
            }
            break;
        }

        switch (CellM->GetType())
        {
          case Triangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0M + 0.5*(Line_Xi[k] + 1)*(zeta1M - zeta0M);
              zeta = 0.5*(zeta + 1);
              switch(infoM)
              {
                case 0: xi=zeta;
                eta=0;
                break;
                case 1: xi=1-zeta;
                eta=zeta;
                break;
                case 2: xi=0;
                eta=1-zeta;
                break;
              }
              Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
            }
            break;

          case Rectangle:
          case Parallelogram:
          case Quadrangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = zeta0M + 0.5*(Line_Xi[k] + 1)*(zeta1M - zeta0M);
              switch(infoM)
              {
                case 0: xi=zeta;
                eta=-1;
                break;
                case 1: xi=1;
                eta=zeta;
                break;
                case 2: xi=-zeta;
                eta=1;
                break;
                case 3: xi=-1;
                eta=-zeta;
                break;
              }
              Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
            }
            break;
        }

        for (j=0;j<N_J_DOF_NM;j++)
          for (k=0;k<N_J_DOF_M;k++)
        {
          LocDofNM = J_DOF_NM[infoNM][j];
          kk = NumbersNM[LocDofNM];
          LocDofM = J_DOF_M[infoM][k];
          ll = NumbersM[LocDofM];

          if (kk < NeumannBound || ll < NeumannBound)
          {
            for (val=l=0;l<N_QuadPoints_Line;l++)
              val += Line_Weights[l]*Flux[l] *
                Ele_NM_Values[l][LocDofNM] * Ele_M_Values[l][LocDofM];

            val *= len_int;

            if (kk < NeumannBound)
            {
              index = RowPtr[kk];
              l = KCol[index];
              while (l != ll)
                l = KCol[++index];
              Entries[index] -= 0.5*val;
            }
            if (ll < NeumannBound)
            {
              index = RowPtr[ll];
              l = KCol[index];
              while (l != kk)
                l = KCol[++index];
              Entries[index] += 0.5*val;
            }
          }
        }

        // check which side is next
        if (lamM0 >= lamNM1 || lamM1 >= lamNM1) NewMortarSideEle = FALSE;
      }
    }
  }

  delete AuxVert0;
  delete AuxVert1;

  // assemble parts on one side of a mortar edge
  // loop over all mortar edges
  N_ = It1->GetN_MortarFace();
  for (i=0;i<N_;i++)
  {
    lev = MAX_ItLevel + (i << 8);
    It1->Init(lev);
    It2->Init(-lev);

    It1->GetPoint(startX, startY);
    It2->GetPoint(endX, endY);
    delX = endX - startX;
    delY = endY - startY;

    // Assumption: n = const. (NM -> M) on each mortar edge !!!!!!
    nlen = sqrt(delX*delX + delY*delY);
    nx = delY / nlen;
    ny = -delX / nlen;

    while (CellNM = It1->Next(infoNM))
    {
      CellNM->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV_NM);

      indexNM = CellNM->GetClipBoard();
      CurrElementID = Space2D->GetFE2D(indexNM, CellNM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      Ele_NM_BF = CurrElement->GetBaseFunct2D();
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_NM = CurrDesc->GetJointDOF();
      N_J_DOF_NM = CurrDesc->GetN_JointDOF();

      NumbersNM = GlobalNumbers + BeginIndex[indexNM];

      AuxVert0 = CellNM->GetVertex(TmpEV_NM[2*infoNM]);
      startX = AuxVert0->GetX();
      startY = AuxVert0->GetY();
      AuxVert0 = CellNM->GetVertex(TmpEV_NM[2*infoNM + 1]);
      endX = AuxVert0->GetX();
      endY = AuxVert0->GetY();

      delX_int = endX - startX;
      delY_int = endY - startY;
      len_int = sqrt(delX_int*delX_int + delY_int*delY_int);

      for (k=0;k<N_QuadPoints_Line;k++)
      {
        lam = 0.5*(Line_Xi[k] + 1);
        X = startX + lam*delX;
        Y = startY + lam*delY;
        DiscreteForm->GetCoeffFct()(1, &X, &Y, &param, &coeff);
        Flux[k] = coeff[1]*nx + coeff[2]*ny;
      }

      switch (CellNM->GetType())
      {
        case Triangle:
          for (k=0;k<N_QuadPoints_Line;k++)
          {
            zeta = 0.5*(Line_Xi[k] + 1);
            switch(infoNM)
            {
              case 0: xi=zeta;
              eta=0;
              break;
              case 1: xi=1-zeta;
              eta=zeta;
              break;
              case 2: xi=0;
              eta=1-zeta;
              break;
            }
            Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
          }
          break;

        case Rectangle:
        case Parallelogram:
        case Quadrangle:
          for (k=0;k<N_QuadPoints_Line;k++)
          {
            zeta = Line_Xi[k];
            switch(infoNM)
            {
              case 0: xi=zeta;
              eta=-1;
              break;
              case 1: xi=1;
              eta=zeta;
              break;
              case 2: xi=-zeta;
              eta=1;
              break;
              case 3: xi=-1;
              eta=-zeta;
              break;
            }
            Ele_NM_BF->GetDerivatives(D00, xi, eta, Ele_NM_Values[k]);
          }
          break;
      }

      for (j=0;j<N_J_DOF_NM;j++)
        for (k=j;k<N_J_DOF_NM;k++)
      {
        LocDofNM = J_DOF_NM[infoNM][j];
        kk = NumbersNM[LocDofNM];
        LocDofNM1 = J_DOF_NM[infoNM][k];
        ll = NumbersNM[LocDofNM1];

        if (kk < NeumannBound || ll < NeumannBound)
        {
          for (val=l=0;l<N_QuadPoints_Line;l++)
            val += Line_Weights[l]*Flux[l] *
              Ele_NM_Values[l][LocDofNM] * Ele_NM_Values[l][LocDofNM1];

          val *= len_int;

          if (kk < NeumannBound)
          {
            index = RowPtr[kk];
            l = KCol[index];
            while (l != ll)
              l = KCol[++index];
            Entries[index] += 0.5*val;
          }
          if (ll < NeumannBound && ll != kk)
          {
            index = RowPtr[ll];
            l = KCol[index];
            while (l != kk)
              l = KCol[++index];
            Entries[index] += 0.5*val;
          }
        }
      }
    }

    while (CellM = It2->Next(infoM))
    {
      CellM->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV_M);

      indexM = CellM->GetClipBoard();
      CurrElementID = Space2D->GetFE2D(indexM, CellM);
      CurrElement = TFEDatabase2D::GetFE2D(CurrElementID);
      Ele_M_BF = CurrElement->GetBaseFunct2D();
      CurrDesc = CurrElement->GetFEDesc2D();

      J_DOF_M = CurrDesc->GetJointDOF();
      N_J_DOF_M = CurrDesc->GetN_JointDOF();

      NumbersM = GlobalNumbers + BeginIndex[indexM];

      AuxVert0 = CellM->GetVertex(TmpEV_M[2*infoM]);
      startX = AuxVert0->GetX();
      startY = AuxVert0->GetY();
      AuxVert0 = CellM->GetVertex(TmpEV_M[2*infoM + 1]);
      endX = AuxVert0->GetX();
      endY = AuxVert0->GetY();

      delX_int = endX - startX;
      delY_int = endY - startY;
      len_int = sqrt(delX_int*delX_int + delY_int*delY_int);

      for (k=0;k<N_QuadPoints_Line;k++)
      {
        lam = 0.5*(Line_Xi[k] + 1);
        X = startX + lam*delX;
        Y = startY + lam*delY;
        DiscreteForm->GetCoeffFct()(1, &X, &Y, &param, &coeff);
        Flux[k] = coeff[1]*nx + coeff[2]*ny;
      }

      switch (CellM->GetType())
      {
        case Triangle:
          for (k=0;k<N_QuadPoints_Line;k++)
          {
            zeta = 0.5*(Line_Xi[k] + 1);
            switch(infoM)
            {
              case 0: xi=zeta;
              eta=0;
              break;
              case 1: xi=1-zeta;
              eta=zeta;
              break;
              case 2: xi=0;
              eta=1-zeta;
              break;
            }
            Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
          }
          break;

        case Rectangle:
        case Parallelogram:
        case Quadrangle:
          for (k=0;k<N_QuadPoints_Line;k++)
          {
            zeta = Line_Xi[k];
            switch(infoM)
            {
              case 0: xi=zeta;
              eta=-1;
              break;
              case 1: xi=1;
              eta=zeta;
              break;
              case 2: xi=-zeta;
              eta=1;
              break;
              case 3: xi=-1;
              eta=-zeta;
              break;
            }
            Ele_M_BF->GetDerivatives(D00, xi, eta, Ele_M_Values[k]);
          }
          break;
      }

      for (j=0;j<N_J_DOF_M;j++)
        for (k=j;k<N_J_DOF_M;k++)
      {
        LocDofM = J_DOF_M[infoM][j];
        kk = NumbersM[LocDofM];
        LocDofM1 = J_DOF_M[infoM][k];
        ll = NumbersM[LocDofM1];

        if (kk < NeumannBound || ll < NeumannBound)
        {
          for (val=l=0;l<N_QuadPoints_Line;l++)
            val += Line_Weights[l]*Flux[l] *
              Ele_M_Values[l][LocDofM] * Ele_M_Values[l][LocDofM1];

          val *= len_int;

          if (kk < NeumannBound)
          {
            index = RowPtr[kk];
            l = KCol[index];
            while (l != ll)
              l = KCol[++index];
            Entries[index] -= 0.5*val;
          }
          if (ll < NeumannBound && ll != kk)
          {
            index = RowPtr[ll];
            l = KCol[index];
            while (l != kk)
              l = KCol[++index];
            Entries[index] -= 0.5*val;
          }
        }
      }
    }
  }

  delete *Ele_NM_Values;
  delete Ele_NM_Values;
  delete *Ele_M_Values;
  delete Ele_M_Values;
  delete Flux;

  delete coeff;
  delete param;

  /*
    int N_Rows, end;
    // print the whole matrix -- MORTAR SDFEM modification
    cout << endl;
    cout << "matrix (MORTAR SDFEM modification):" << endl;
    N_Rows = Matrix->GetN_Rows();
    for(i=0;i<N_Rows;i++)
    {
      end = RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << KCol[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  */
}
#endif                                            // __ADD_LINK_SDFEM__

/** assemble the matrix entries */
void Assemble(TMatrix2D *matrix)
{
  TStructure *structure = matrix->GetStructure();
  TFESpace1D *TestSpace = (TFESpace1D *) structure->GetTestSpace();
  TFESpace2D *AnsatzSpace = (TFESpace2D *) structure->GetAnsatzSpace();
  TCollection *coll2D = AnsatzSpace->GetCollection();
  TCollection *coll1D = TestSpace->GetCollection();
  int N_Cells, N_Edges, i, j, k, l, m, n, o;
  int **JointDOFs, N_JointDOFs;
  const int *TmpEV;
  double Z00, Z01, Z10, Z11;
  double X, Y, X0, Y0, X1, Y1;
  bool MortarSide;
  TBaseCell *CurrCell1D, *CurrCell2D;
  TVertex *CurrVertex;
  TQuadFormula1D *line_formula;
  int N_QuadPoints_Line;
  double *Line_Weights, *Line_Xi;
  TBaseFunct2D *Ele2DBF;
  TBaseFunct1D *LineBF;
  double **LineValues, **Ele2DValues;
  double zeta, xi, eta, sum;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *TestNumbers, *AnsatzNumbers;
  int TestNode, AnsatzNode, AnsatzNode_loc;
  int *RowPtr = matrix->GetRowPtr();
  int *KCol = matrix->GetKCol();
  double *Entries = matrix->GetEntries();
  TFE1D *CurrEle_1D;
  TFE2D *CurrEle_2D;
  TFEDesc2D *CurrDesc;
  FE1D CurrEleID_1D;
  FE2D CurrEleID_2D;
  JointType type;

  line_formula = TFEDatabase2D::GetQuadFormula1D(Gauss3Line);
  line_formula->GetFormulaData(N_QuadPoints_Line, Line_Weights, Line_Xi);

  Ele2DValues = new double*[N_QuadPoints_Line];
  LineValues = new double*[N_QuadPoints_Line];

  *Ele2DValues = new double[MaxN_BaseFunctions2D * N_QuadPoints_Line];
  for (i=1;i<N_QuadPoints_Line;i++)
    Ele2DValues[i] = *Ele2DValues + i*MaxN_BaseFunctions2D;

  *LineValues = new double[MaxN_BaseFunctions1D * N_QuadPoints_Line];
  for (i=1;i<N_QuadPoints_Line;i++)
    LineValues[i] = *LineValues + i*MaxN_BaseFunctions1D;

  TestGlobalNumbers = TestSpace->GetGlobalNumbers();
  AnsatzGlobalNumbers = AnsatzSpace->GetGlobalNumbers();

  TestBeginIndex = TestSpace->GetBeginIndex();
  AnsatzBeginIndex = AnsatzSpace->GetBeginIndex();

  N_Cells = coll2D->GetN_Cells();
  for (i=0;i<N_Cells;i++)
  {
    CurrCell2D = coll2D->GetCell(i);
    N_Edges = CurrCell2D->GetN_Edges();
    for (j=0;j<N_Edges;j++)
      if ((type = CurrCell2D->GetJoint(j)->GetType()) == MortarJoint ||
      type == MortarBaseJoint)
    {
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

      CurrCell2D->GetRefDesc()->GetShapeDesc()->GetEdgeVertex(TmpEV);

      if (type == MortarJoint)
        l = ((TMortarJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();
      else
        l = ((TMortarBaseJoint *) CurrCell2D->GetJoint(j))->GetMEdgeInColl();

      if (l < 0)
      {
        MortarSide = TRUE;
        l = -(l + 1);
      }
      else
        MortarSide = FALSE;

      // !!! important
      // at the moment, this is only made for conform elements

      CurrEleID_2D = AnsatzSpace->GetFE2D(i, CurrCell2D);
      CurrEle_2D = TFEDatabase2D::GetFE2D(CurrEleID_2D);
      Ele2DBF = CurrEle_2D->GetBaseFunct2D();
      CurrDesc = CurrEle_2D->GetFEDesc2D();

      JointDOFs = CurrDesc->GetJointDOF();
      N_JointDOFs = CurrDesc->GetN_JointDOF();

      if (!MortarSide)
      {
        Z00 = -1;
        Z01 = 1;
        Z10 = -1;
        Z11 = 1;

        CurrCell1D = coll1D->GetCell(l);
        CurrEleID_1D = TestSpace->GetFE1D(l, CurrCell1D);
        CurrEle_1D = TFEDatabase2D::GetFE1D(CurrEleID_1D);
        LineBF = CurrEle_1D->GetBaseFunct1D();
      }

      CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j]);
      X0 = CurrVertex->GetX();
      Y0 = CurrVertex->GetY();
      CurrVertex = CurrCell2D->GetVertex(TmpEV[2*j + 1]);
      X1 = CurrVertex->GetX();
      Y1 = CurrVertex->GetY();

      do
      {
        if (MortarSide)
        {
          CurrCell1D = coll1D->GetCell(l);
          CurrEleID_1D = TestSpace->GetFE1D(l, CurrCell1D);
          CurrEle_1D = TFEDatabase2D::GetFE1D(CurrEleID_1D);
          LineBF = CurrEle_1D->GetBaseFunct1D();

          CurrVertex = CurrCell1D->GetVertex(0);
          X = CurrVertex->GetX();
          Y = CurrVertex->GetY();
          if ((X0-X)*(X0-X1) + (Y0-Y)*(Y0-Y1) < 1e-16) break;

          Z01 = GetZeta(CurrCell2D->GetVertex(TmpEV[2*j]),
            CurrCell2D->GetVertex(TmpEV[2*j + 1]),
            CurrCell1D->GetVertex(1));
          Z00 = GetZeta(CurrCell2D->GetVertex(TmpEV[2*j]),
            CurrCell2D->GetVertex(TmpEV[2*j + 1]),
            CurrCell1D->GetVertex(0));
          Z10 = GetZeta(CurrCell1D->GetVertex(0), CurrCell1D->GetVertex(1),
            CurrCell2D->GetVertex(TmpEV[2*j + 1]));
          Z11 = GetZeta(CurrCell1D->GetVertex(0), CurrCell1D->GetVertex(1),
            CurrCell2D->GetVertex(TmpEV[2*j]));

          if (Z01 < -1) Z01 = -1;
          if (Z10 < -1) Z10 = -1;
          if (Z00 > 1) Z00 = 1;
          if (Z11 > 1) Z11 = 1;
        }

        for (k=0;k<N_QuadPoints_Line;k++)
        {
          zeta = Z10 + 0.5*(Line_Xi[k] + 1)*(Z11 - Z10);
          LineBF->GetDerivatives(D0, zeta, LineValues[k]);
        }

        switch (CurrCell2D->GetType())
        {
          case Triangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = Z00 + 0.5*(Line_Xi[k] + 1)*(Z01 - Z00);
              zeta = 0.5*(zeta + 1);
              switch(j)
              {
                case 0: xi=zeta;
                eta=0;
                break;
                case 1: xi=1-zeta;
                eta=zeta;
                break;
                case 2: xi=0;
                eta=1-zeta;
                break;
              }
              Ele2DBF->GetDerivatives(D00, xi, eta, Ele2DValues[k]);
            }
            break;

          case Rectangle:
          case Parallelogram:
          case Quadrangle:
            for (k=0;k<N_QuadPoints_Line;k++)
            {
              zeta = Z00 + 0.5*(Line_Xi[k] + 1)*(Z01 - Z00);
              switch(j)
              {
                case 0: xi=zeta;
                eta=-1;
                break;
                case 1: xi=1;
                eta=zeta;
                break;
                case 2: xi=-zeta;
                eta=1;
                break;
                case 3: xi=-1;
                eta=-zeta;
                break;
              }
              Ele2DBF->GetDerivatives(D00, xi, eta, Ele2DValues[k]);
            }
            break;
        }

        TestNumbers = TestGlobalNumbers + TestBeginIndex[l];
        n = TestBeginIndex[l+1] - TestBeginIndex[l];

        for (k=0;k<n;k++)
          for (m=0;m<N_JointDOFs;m++)
        {
          AnsatzNode_loc = JointDOFs[j][m];
          AnsatzNode = AnsatzNumbers[AnsatzNode_loc];

          for (sum=o=0;o<N_QuadPoints_Line;o++)
            sum += Line_Weights[o] * Ele2DValues[o][AnsatzNode_loc] *
              LineValues[o][k];

          sum *= 0.25*(Z01 - Z00) * sqrt((X1-X0)*(X1-X0) + (Y1-Y0)*(Y1-Y0));

          TestNode = TestNumbers[k];
          o = RowPtr[TestNode];
          while (KCol[o] != AnsatzNode && o < RowPtr[TestNode+1]) o++;

          if (o == RowPtr[TestNode+1])
          {
            cerr << "Error in Assemble(Mortar): " << i << " : "
              << o << endl;
            exit(-1);
          }

          Entries[o] += sum;
        }

        l++;

      } while ((CurrCell1D = CurrCell1D->GetJoint(1)->
        GetNeighbour(CurrCell1D)) && MortarSide);
    }
  }

  /*
    int N_Rows = matrix->GetN_Rows();
    // print the whole matrix
    cout << endl;
    cout << "mortar matrix" << endl;

    for(i=0;i<N_Rows;i++)
    {
      k = RowPtr[i + 1];
      for(j=RowPtr[i];j<k;j++)
      {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << KCol[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  */
}
#endif                                            // MORTAR

void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
TDiscreteForm2D *DiscreteForm,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters,
TAuxParam2D *ParametersBound,
TypeBoundSwitchFunct2D *TypeBoundSwitcher,
int *CounterBoundaryParam
#ifdef __3D__
, TAux2D3D *Aux2D3D
#endif
)
{
  double hK;
  int N_AllMatrices = n_sqmatrices+n_matrices;
  int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ii,jj,ll;
  int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TFESpace2D *fespace;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
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
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double Xi[MaxN_PointsForNodal2D], Eta[MaxN_PointsForNodal2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *Param[MaxN_QuadPoints_2D];
  double *local_rhs;
  double *righthand;
  double **Matrices, *aux;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *AuxArray[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  double *RHS, *MatrixRow;
  double **HangingEntries, **HangingRhs;
  double *CurrentHangingEntries, *CurrentHangingRhs;
  int *HangingRowPtr, *HangingColInd;
  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  double *Coupling, v;
  TBoundComp *BoundComp;
  double t0, t1, t, s,integral[2];
  int comp, dof_ii,dof_jj, found;
  BoundCond Cond0, Cond1;
  BoundCondFunct2D *BoundaryCondition;
  BoundValueFunct2D *BoundaryValue;
  TNodalFunctional2D *nf;
  int N_EdgePoints;
  double *EdgePoints;
  double PointValues[MaxN_PointsForNodal2D];
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-4;
#ifdef __3D__
  double z0, z1;
#endif
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue;
  bool *SecondDer;
  int *Bounds, NeumannBound, RobinBound;
  int lr, mr;

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];
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

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  if(ParametersBound)
    N_Parameters = MAX(Parameters->GetN_Parameters(),
      ParametersBound->GetN_Parameters());
  else
    N_Parameters = Parameters->GetN_Parameters();

#ifdef __3D__
  N_Parameters += 7;                              // (u, ux, uy, uz, nx, ny, nz)
#endif

  if(N_Parameters)
  {
    aux = new double [MaxN_QuadPoints_2D*N_Parameters];
    for(j=0;j<MaxN_QuadPoints_2D;j++)
      Param[j] = aux + j*N_Parameters;
  }

  // 20 <= number of term in bilinear form
  aux = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    AuxArray[j] = aux + j*20;

  if(N_AllMatrices)
  {
    aux = new double
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  SecondDer = DiscreteForm->GetNeeds2ndDerivatives();

  // ########################################################################
  // loop over all cells
  // ########################################################################
  Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    hK = cell->GetDiameter();

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    for(j=0;j<n_fespaces;j++)
    {
      CurrentElement = fespaces[j]->GetFE2D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = n_fespaces;

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
      Coll, cell, SecondDer,
      N_Points, xi, eta, weights, X, Y, AbsDetjk);

    Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);

#ifdef __3D__
    if(Aux2D3D)
      Aux2D3D->GetGradient(i, N_Points, xi, eta, Param);
#endif

    // use DiscreteForm to assemble a few matrices and
    // right-hand sides at once
    if(DiscreteForm)
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk,
        hK, X, Y,
        LocN_BF, LocBF,
        Param, AuxArray,
        cell,
        N_AllMatrices, n_rhs,
        LocMatrices, LocRhs);

    N_Joints=cell->GetN_Joints();

    // ####################################################################
    // add local matrices to global matrices (ansatz == test)
    // ####################################################################
    for(j=0;j<n_sqmatrices;j++)
    {
      // find space for this bilinear form
      fespace = sqmatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE2D(i, cell);
      N_ = N_BaseFunct[CurrentElement];

      Matrix = Matrices+j*MaxN_BaseFunctions2D;
      Entries = sqmatrices[j]->GetEntries();
      RowPtr = sqmatrices[j]->GetRowPtr();
      ColInd = sqmatrices[j]->GetKCol();

      CurrentHangingEntries = HangingEntries[j];
      HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
      HangingColInd = sqmatrices[j]->GetHangingKCol();

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = GlobalNumbers[j] + BeginIndex[j][i];

      Bounds = fespace->GetBoundaryNodesBound();
      NeumannBound = Bounds[0];
      RobinBound = Bounds[1];

      BoundaryCondition = BoundaryConditions[j];
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
           joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }

          // get id of the boundary component
          comp = BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          BoundaryCondition(comp, t0, Cond0);
          if(Cond0 == ROBIN)
          {
            // Robin
            lr = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);

            // get a suitable line quadrature formula
            LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*lr);
            qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
            qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

            TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
              ->MakeRefElementData(LineQuadFormula);

            JointValues=TFEDatabase2D::GetJointValues2D(
              BaseFuncts[CurrentElement], LineQuadFormula, m);
            // get vertices of boundary edge
#ifdef __2D__
            cell->GetVertex(m)->GetCoords(x0, y0);
            cell->GetVertex((m+1) % 4)->GetCoords(x1, y1);
#endif
            // compute (half of the) length of the boundary edge
            hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
            // cout << "x0: " << x0 << " y0: " << y0 << endl;
            // cout << "x1: " << x1 << " y1: " << y1 << endl;
            // compute boundary integral
            if(TypeBoundSwitcher(comp,0.5*(t1+t0))==0)
            {
              // ParametersBound->GetParameters(N_LinePoints, cell, i,
              //				     zeta, m, Param);
              for(n=0;n<N_LinePoints;n++)
              {
                // values of test functions in this quadrature point
                JointValue = JointValues[n];
                // cout << "Zeta  :" << zeta[n] << endl;

                // get quadrature point on the boundary
                for(l=0;l<N_;l++)
                {
                  MatrixRow = Matrix[l];
                  s = JointValue[l];
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[n];
                  // !! hold the robin boundary values of
                  // !! the function alpha from parameters
                  // s *= alpha;
                  // update rhs for all test functions
                  // s *= Param[n][CounterBoundaryParam[j]];
                  for(k=0;k<N_;k++)
                    MatrixRow[k] += s*JointValue[k];
                }                                 // endfor l
              }                                   // endfor n
            }                                     // if switcher
            else
            {
              if(!ParametersBound)
              {
                Error("ParametersBound not set" << endl);
                exit(-1);
              }
              ParametersBound->GetParameters(N_LinePoints, Coll, cell, i,
                zeta, m, Param);
              for(n=0;n<N_LinePoints;n++)
              {
                // values of test functions in this quadrature point
                JointValue = JointValues[n];
                // cout << "Zeta  :" << zeta[n] << endl;

                // get quadrature point on the boundary
                for(l=0;l<N_;l++)
                {
                  MatrixRow = Matrix[l];
                  s = JointValue[l];
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  s *= hE * LineWeights[n];
                  // !! hold the robin boundary values of
                  // !! the function alpha from parameters
                  s *= Param[n][CounterBoundaryParam[j]];
                  // update rhs for all test functions
                  for(k=0;k<N_;k++)
                    MatrixRow[k] += s*JointValue[k];
                }                                 // endfor l
              }                                   // endfor n

            }
          }                                       // end Robin
        }                                         // endif BoundEdge
      }                                           // endfor m

      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          end=RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            for(k=0;k<N_;k++)
            {
              if(DOF[k] == ColInd[n])
              {
                // cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                break;
              }                                   // endif
            }                                     // endfor k
          }                                       // endfor n
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            end = HangingRowPtr[l+1];
            for(n=HangingRowPtr[l];n<end;n++)
            {
              for(k=0;k<N_;k++)
              {
                if(DOF[k] == HangingColInd[n])
                {
                  CurrentHangingEntries[n] += MatrixRow[k];
                  break;
                }                                 // endif
              }                                   // endfor k
            }                                     // endfor n
          }
          else
          {
            // Dirichlet node
            n=RowPtr[l];
            if(ColInd[n]==l)
            {
              Entries[n]=1.0;
            }
          }
        }
      }                                           // endfor m
    }                                             // endfor j
    // ####################################################################
    // add local matrices to global matrices (ansatz != test)
    // ####################################################################
    for(j=0;j<n_matrices;j++)
    {
      TestElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetTestSpace())->GetFE2D(i, cell);
      AnsatzElement = ((TFESpace2D *) matrices[j]->GetStructure()->
        GetAnsatzSpace())->GetFE2D(i, cell);

      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;

      Entries = matrices[j]->GetEntries();
      RowPtr = matrices[j]->GetRowPtr();
      ColInd = matrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        end=RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Ansatz;k++)
          {
            if(AnsatzDOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += MatrixRow[k];
              break;
            }                                     // endif
          }                                       // endfor k
        }                                         // endfor n
      }                                           // endfor m
    }                                             // endfor j  (n_matrices)

    // ####################################################################
    // add local right-hand sides to global right-hand side
    // ####################################################################
    for(j=0;j<n_rhs;j++)
    {
      fespace = ferhs[j];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE2D(i, cell);

      N_ = N_BaseFunct[CurrentElement];

      local_rhs = righthand+j*MaxN_BaseFunctions2D;
      RHS = rhs[j];
      CurrentHangingRhs = HangingRhs[j];
      // find space for this linear form

      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m];
          // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
        }                                         // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            l -= ActiveBound;
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      }                                           // endfor m

      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase2D::GetFE2D(CurrentElement);
      nf = ele->GetNodalFunctional2D();
      nf->GetPointsForEdge(N_EdgePoints, EdgePoints);

      FEDesc_Obj = ele->GetFEDesc2D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();

      // setting Dirichlet boundary condition
      N_Joints = cell->GetN_Edges();
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
           joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }
          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            switch(Cond0)
            {
              case DIRICHLET:
                // analytic boundary condition

                if(TypeBoundSwitcher(comp,0.5*(t1+t0))==0)
                {
                  // read boundary values for each quadrature point
                  for(l=0;l<N_EdgePoints;l++)
                  {
                    s = EdgePoints[l];
                    t = 0.5*(t0*(1-s) + t1*(1+s));
                    BoundaryValue(comp, t, PointValues[l]);
                  }                               // endfor l
                  // compute boundary values for each dof on the
                  // boundary edge with the nodal functionals
                  nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                    FunctionalValues);
                  EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                  // save boundary values of each dof on the boundary
                  // edge in the rhs
                  for(l=0;l<N_EdgeDOF;l++)
                  {
                    RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                  }
                }
                else
                {
                  if(!ParametersBound)
                  {
                    Error("ParametersBound not set" << endl);
                    exit(-1);
                  }
                  ParametersBound->GetParameters(N_EdgePoints, Coll, cell, i,
                    EdgePoints, m, Param);
                  for(l=0;l<N_EdgePoints;l++)
                  {
                    PointValues[l]=Param[l][CounterBoundaryParam[j]];
                    // cout << "test Point Values " << PointValues[l] << endl;
                  }

                  nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
                    FunctionalValues);
                  EdgeDOF = FEDesc_Obj->GetJointDOF(m);
                  // save boundary values of each dof on the boundary
                  // edge in the rhs
                  for(l=0;l<N_EdgeDOF;l++)
                  {
                    RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
                    //  cout << "hallo Functionals: " << FunctionalValues[l] << endl;
                  }
                }                                 // end else switcher
                break;

              case NEUMANN:
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                if(TypeBoundSwitcher(comp,0.5*(t1+t0))==0)
                {
                  for(l=0;l<N_LinePoints;l++)
                  {
                    // values of test functions in this quadrature point
                    JointValue = JointValues[l];
                    // get quadrature point on the boundary
                    t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                    // get value in this quadrature point (in t)
                    // analytic bc
                    BoundaryValue(comp, t, s);
                    // multiply value with weights from quadrature formula
                    // and determinant from integral transformation to the
                    // unit edge (-1,1)
                    s *= hE * LineWeights[l];
                    // update rhs for all test functions
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += s*JointValue[k];
                  }                               // end for l
                }
                else
                {
                  if(!ParametersBound)
                  {
                    Error("ParametersBound not set" << endl);
                    exit(-1);
                  }
                  ParametersBound->GetParameters(N_EdgePoints, Coll, cell, i,
                    EdgePoints, m, Param);
                  for(l=0;l<N_LinePoints;l++)
                  {
                    JointValue = JointValues[l];
                    s = Param[l][CounterBoundaryParam[j]]*hE*LineWeights[l];
                    // cout << "Param Neumann" << Param[l][CounterBoundaryParam] << endl;
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += s*JointValue[k];
                  }
                }
                break;
              case ROBIN:
                // get polynomial degree of fe
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(
                  BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
                  ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
                if(TypeBoundSwitcher(comp,0.5*(t1+t0))==0)
                {
                  for(l=0;l<N_LinePoints;l++)
                  {
                    // values of test functions in this quadrature point
                    JointValue = JointValues[l];
                    // get quadrature point on the boundary
                    t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                    // get value in this quadrature point (in s)
                    // analytic bc
                    BoundaryValue(comp, t, s);
                    // multiply value with weights from quadrature formula
                    // and determinant from integral transformation to the
                    // unit edge (-1,1)
                    s *= hE * LineWeights[l];
                    // update rhs for all test functions
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += s*JointValue[k];
                  }                               // end for l
                }
                else
                {
                  if(!ParametersBound)
                  {
                    Error("ParametersBound not set" << endl);
                    exit(-1);
                  }
                  ParametersBound->GetParameters(N_EdgePoints, Coll, cell, i,
                    EdgePoints, m, Param);
                  for(l=0;l<N_LinePoints;l++)
                  {
                    JointValue = JointValues[l];
                    s = Param[l][CounterBoundaryParam[j]+1]*hE*LineWeights[l];
                    // cout << "Param Robin" << Param[l][CounterBoundaryParam+1] << endl;
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += s*JointValue[k];
                  }
                }
                break;

                /*
                      // get polynomial degree of fe
                      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                        (CurrentElement);
                      // get a suitable line quadrature formula
                      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                      TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                        ->MakeRefElementData(LineQuadFormula);
                      JointValues=TFEDatabase2D::GetJointValues2D(
                BaseFuncts[CurrentElement], LineQuadFormula, m);
                // get vertices of boundary edge
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
                // compute (half of the) length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
                // compute boundary integral
                for(l=0;l<N_LinePoints;l++)
                {
                // values of test functions in this quadrature point
                JointValue = JointValues[l];
                // get quadrature point on the boundary
                t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                // get value in this quadrature point (in s)
                BoundaryValue(comp, t, s);
                // multiply value with weights from quadrature formula
                // and determinant from integral transformation to the
                // unit edge (-1,1)
                s *= hE * LineWeights[l];
                // update rhs for all test functions
                for(k=0;k<N_;k++)
                if((l3 = DOF[k])<ActiveBound)
                RHS[l3] += s*JointValue[k];
                }
                break;
                */
              case SLIP:
                exit(4711);
                break;
              case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble2DSlipBC, see below
                break;
         default:

         break;
            }                                     // endswitch Cond0

          }                                       // endif (Cond0==Cond1)
          else
          {
            OutPut("different boundary condition on one edge ");
            OutPut("are not allowed!" << endl);
            exit(4711);
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor j (n_rhs)
  }                                               // endfor i (N_Cells)

  // ####################################################################
  // modify matrix according to coupling
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j];
    HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
    HangingColInd = sqmatrices[j]->GetHangingKCol();

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      end = HangingRowPtr[i+1];
      for(n=HangingRowPtr[i];n<end;n++)
      {
        v = CurrentHangingEntries[n];
        m = HangingColInd[n];
        for(l=0;l<k;l++)
        {
          l1 = DOF[l];
          if(l1<ActiveBound)
          {
            last=RowPtr[l1+1];
            for(l2=RowPtr[l1];l2<last;l2++)
            {
              if(ColInd[l2] == m)
              {
                Entries[l2] += Coupling[l] * v;
              }
            }                                     // endfor l2
          }                                       // endif
        }                                         // endfor l
      }                                           // endfor n
    }                                             // endfor i
  }                                               // endfor j

  for(j=0;j<n_rhs;j++)
  {
    fespace = ferhs[j];
    N_Hanging = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    RHS = rhs[j];
    CurrentHangingRhs = HangingRhs[j];

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      N_ = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      for(k=0;k<N_;k++)
      {
        l = DOF[k];
        if(l<ActiveBound)
        {
          RHS[l] += Coupling[k] * CurrentHangingRhs[i];
        }
      }                                           // endfor k
    }                                             // endfor i
  }                                               // endfor j

  // ####################################################################
  // write coupling into matrix
  // ####################################################################
  for(j=0;j<n_sqmatrices;j++)
  {
    fespace = sqmatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    HangingNodes = fespace->GetHangingNodes();

    Entries = sqmatrices[j]->GetEntries();
    RowPtr = sqmatrices[j]->GetRowPtr();
    ColInd = sqmatrices[j]->GetKCol();

    ActiveBound = fespace->GetActiveBound();

    n = RowPtr[ActiveBound];

    for(i=0;i<N_;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      Entries[n] = 1.0;
      n++;

      for(l=0;l<k;l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      }                                           // endfor l

    }                                             // endfor i
  }

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

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete HangingRhs[i];
    delete HangingRhs;

    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }

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

  delete AuxArray[0];

  /*
  int N_Rows;
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  for(k=0;k<n_sqmatrices;k++)
  {
    cout << endl;
    cout << "sqmatrix: " << k << endl;
    RowPtr = sqmatrices[k]->GetRowPtr();
    Entries = sqmatrices[k]->GetEntries();
  ColInd = sqmatrices[k]->GetKCol();
  N_Rows = sqmatrices[k]->GetN_Rows();
  for(i=0;i<N_Rows;i++)
  {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(10) << i << setw(10) << ColInd[j] << "   ";
  cout << setw(20) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_matrices;k++)
  {
  cout << endl;
  cout << "matrix: " << k << endl;
  RowPtr = matrices[k]->GetRowPtr();
  Entries = matrices[k]->GetEntries();
  ColInd = matrices[k]->GetKCol();
  N_Rows = matrices[k]->GetN_Rows();
  for(i=0;i<N_Rows;i++)
  {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(10) << i << setw(10) << ColInd[j] << "   ";
  cout << setw(20) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_rhs;k++)
  {
  cout << "rhs: " << k << endl;
  N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
  RHS=rhs[k];
  for(i=0;i<N_Rows;i++)
  cout << setw(10) << i << setw(20) << RHS[i] << endl;
  }

  */

}                                                 // end of Assemble


void ModifyDivergenceMatrices(int N_Active, int N_Dirichlet,
TMatrix2D *B1, TMatrix2D *B2)
{
  int i, entries, *columns;
  double *B1Entries, *B2Entries;

  entries = B1->GetN_Entries();
  columns =   B1->GetKCol();
  B1Entries = B1->GetEntries();
  B2Entries = B2->GetEntries();

  for (i=0;i<entries;i++)
  {
    if ((columns[i]>= N_Active)&&(columns[i] < N_Dirichlet))
    {
      B1Entries[i] = 0;
      B2Entries[i] = 0;
    }
  }

  return;
}


#ifdef __2D__

void ComputeInterfaceConditionDirichlet(TSquareMatrix2D  **sqmatrices,
TMatrix2D **matrices,
TFEVectFunct2D *u_other,
double *rhs,
int *interface_dof_Array,
int *interface_dof_other_Array,
int interface_dof,
double density,
double density_other)
{
  int i,j, N_DOF;
  int start,end,indexColumn,indexRow;
  int *RowPtrA11, *KColA11, *RowPtrA12, *KColA12, *RowPtrA21, *KColA21, *RowPtrA22;
  int *KColA22, *RowPtrB1T, *KColB1T, *RowPtrB2T, *KColB2T;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22, *EntriesB1T, *EntriesB2T;
  TSquareMatrix *A11, *A12, *A21, *A22;
  TMatrix2D *B1T, *B2T;
  double t1=1, t2=7.113e-4,n,n1,n2, *u1, *u2, u, v;
  TFEFunction2D *ufct;

  //OutPut("Dirichlet" << endl);
  n=sqrt(t1*t1+t2*t2);
  t1/=n;
  t2/=n;
  n1=-t2;
  n2=t1;

  ufct = u_other->GetComponent(0);
  u1 = ufct->GetValues();
  delete ufct;
  ufct = u_other->GetComponent(1);
  u2 = ufct->GetValues();
  delete ufct;

  /*u1 = u_other->GetComponent(0)->GetValues();
  u2 = u_other->GetComponent(1)->GetValues();
  */

  A11 = sqmatrices[0];
  A12 = sqmatrices[1];
  A21 = sqmatrices[2];
  A22 = sqmatrices[3];
  B1T = matrices[0];
  B2T = matrices[1];

  RowPtrA11 = A11->GetRowPtr();
  KColA11 = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  N_DOF = A11->GetN_Rows();

  RowPtrA12 = A12->GetRowPtr();
  KColA12 = A12->GetKCol();
  EntriesA12 = A12->GetEntries();

  RowPtrA21 = A21->GetRowPtr();
  KColA21 = A21->GetKCol();
  EntriesA21 = A21->GetEntries();

  RowPtrA22 = A22->GetRowPtr();
  KColA22 = A22->GetKCol();
  EntriesA22 = A22->GetEntries();

  RowPtrB1T = B1T->GetRowPtr();
  KColB1T = B1T->GetKCol();
  EntriesB1T = B1T->GetEntries();

  RowPtrB2T = B2T->GetRowPtr();
  KColB2T = B2T->GetKCol();
  EntriesB2T = B2T->GetEntries();

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA11[indexRow];
    end = RowPtrA11[indexRow+1];
    for (j=start;j<end;j++)
    {
      indexColumn = KColA11[j];
      if(indexColumn == indexRow)
      {
        EntriesA11[j] = 1;
      }
      else
      {
        EntriesA11[j] = 0;
      }
    }
    u = u1[interface_dof_other_Array[i]];
    v = u2[interface_dof_other_Array[i]];
    rhs[indexRow] = ((density_other/density)*n1*n1+n2*n2)*u + (density_other/density-1)*v*n1*n2;
    // OutPut("diri " << i << " " <<indexRow << " " << interface_dof_other_Array[i] << " " << u << " " << v <<
    //   " " <<  rhs[indexRow] << " " << density_other << " " << density << endl);
  }

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA12[indexRow];
    end = RowPtrA12[indexRow+1];
    for (j=start;j<end;j++)
    {
      EntriesA12[j] = 0;
    }
  }

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA22[indexRow];
    end = RowPtrA22[indexRow+1];
    for (j=start;j<end;j++)
    {
      indexColumn = KColA22[j];
      if(indexColumn == indexRow)
      {
        EntriesA22[j] = 1;
      }
      else
      {
        EntriesA22[j] = 0;
      }
    }
    u = u1[interface_dof_other_Array[i]];
    v = u2[interface_dof_other_Array[i]];
    rhs[indexRow+N_DOF] = -(1-density_other/density)*u*n1*n2 + v*(n1*n1+(density_other/density)*n2*n2);
  }

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA21[indexRow];
    end = RowPtrA21[indexRow+1];
    for (j=start;j<end;j++)
    {
      EntriesA21[j] = 0;
    }
  }

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrB1T[indexRow];
    end = RowPtrB1T[indexRow+1];
    for (j=start;j<end;j++)
    {
      EntriesB1T[j] = 0;
    }
  }

  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrB2T[indexRow];
    end = RowPtrB2T[indexRow+1];
    for (j=start;j<end;j++)
    {
      EntriesB2T[j] = 0;
    }
  }
  return;
}


void ComputeInterfaceConditionStress(TSquareMatrix2D  **sqmatrices,
TMatrix2D **matrices,
TFEVectFunct2D *u_other,
double *rhs,
int *interface_dof_Array,
int *interface_dof_other_Array,
int interface_dof, double density,
double density_other,
double *rhs_for_stress_other)
{
  int i,j, N_DOF, N_DOF_other;
  int start,end,indexColumn,indexRow;
  int *RowPtrA11, *KColA11, *RowPtrA12, *KColA12, *RowPtrA21, *KColA21;
  int *RowPtrA22, *KColA22, *RowPtrB1T, *KColB1T, *RowPtrB2T, *KColB2T;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22, *EntriesB1T, *EntriesB2T;
  double entrieA11,entrieA21,entrieA22,entrieA12,entrieB1T,entrieB2T;
  TSquareMatrix *A11, *A12, *A21, *A22;
  TMatrix2D *B1T, *B2T;
  double t1=1, t2=7.113e-4,n,n1,n2, *u1, *u2, u, v, scale = 1.0;
  TFEFunction2D *ufct;

  //OutPut("Stress" << endl);
  // compute normal and tangential
  n=sqrt(t1*t1+t2*t2);
  t1/=n;
  t2/=n;
  n1=-t2;
  n2=t1;
  // get velocity of the other subdomain

  ufct = u_other->GetComponent(0);
  u1 = ufct->GetValues();
  N_DOF_other = ufct->GetLength();
  delete ufct;
  ufct = u_other->GetComponent(1);
  u2 = ufct->GetValues();
  delete ufct;

  /*u1 = u_other->GetComponent(0)->GetValues();
  u2 = u_other->GetComponent(1)->GetValues();
  N_DOF_other = u_other->GetComponent(0)->GetLength();
  */

  // get matrices
  A11 = sqmatrices[0];
  A12 = sqmatrices[1];
  A21 = sqmatrices[2];
  A22 = sqmatrices[3];
  B1T = matrices[0];
  B2T = matrices[1];

  RowPtrA11 = A11->GetRowPtr();
  KColA11 = A11->GetKCol();
  EntriesA11 = A11->GetEntries();
  N_DOF = A11->GetN_Rows();

  RowPtrA12 = A12->GetRowPtr();
  KColA12 = A12->GetKCol();
  EntriesA12 = A12->GetEntries();

  RowPtrA21 = A21->GetRowPtr();
  KColA21 = A21->GetKCol();
  EntriesA21 = A21->GetEntries();

  RowPtrA22 = A22->GetRowPtr();
  KColA22 = A22->GetKCol();
  EntriesA22 = A22->GetEntries();

  RowPtrB1T = B1T->GetRowPtr();
  KColB1T = B1T->GetKCol();
  EntriesB1T = B1T->GetEntries();

  RowPtrB2T = B2T->GetRowPtr();
  KColB2T = B2T->GetKCol();
  EntriesB2T = B2T->GetEntries();

  //  OutPut("dofs " << N_DOF << " " << N_DOF_other << endl);
  // add rows of matrices A11 and A21 to A11
  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA11[indexRow];
    end = RowPtrA11[indexRow+1];
    for (j=start;j<end;j++)
    {
      entrieA11 = EntriesA11[j];
      entrieA21 = EntriesA21[j];
      EntriesA11[j] = n2*entrieA11 - n1*entrieA21;
      EntriesA11[j] *= scale;
    }
    rhs[indexRow] = -n2*rhs_for_stress_other[interface_dof_other_Array[i]]
      + n1*rhs_for_stress_other[interface_dof_other_Array[i]+N_DOF_other];
    rhs[indexRow] = scale*rhs[indexRow];
    //     OutPut("rhs " << indexRow << " " <<  rhs[indexRow] << " ");
    //OutPut(i << " " << interface_dof_other_Array[i] << " " << rhs_for_stress_other[interface_dof_other_Array[i]]
    //	     << " " << rhs_for_stress_other[interface_dof_other_Array[i]+N_DOF_other] << endl);
    // set rhs for additional equation in A21,A22,B2T
    u = u1[interface_dof_other_Array[i]];
    v = u2[interface_dof_other_Array[i]];
    rhs[indexRow+N_DOF] = density_other*(u*n1+v*n2);
    //rhs[indexRow+N_DOF] = -rhs[indexRow+N_DOF];
    //OutPut(indexRow+N_DOF << " " <<  rhs[indexRow+N_DOF] << endl);
    //OutPut(u << " " << n1 << " " << v << " "<< n2 << endl);
  }

  // add rows of matrices B1T and B2T to B1T
  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrB1T[indexRow];
    end = RowPtrB1T[indexRow+1];
    for (j=start;j<end;j++)
    {
      entrieB1T = EntriesB1T[j];
      entrieB2T = EntriesB2T[j];
      EntriesB1T[j] = n2*entrieB1T - n1*entrieB2T;
      EntriesB1T[j] *= scale;
      EntriesB2T[j] = 0;
    }
  }

  //  add rows of matrices A21 and A22 to A21
  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA21[indexRow];
    end = RowPtrA21[indexRow+1];
    for (j=start;j<end;j++)
    {
      entrieA12 = EntriesA12[j];
      entrieA22 = EntriesA22[j];
      EntriesA12[j] = n2*entrieA12 - n1*entrieA22;
      EntriesA12[j] *= scale;
    }
  }

  // modify matrix A21 for additional equations
  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA21[indexRow];
    end = RowPtrA21[indexRow+1];
    for (j=start;j<end;j++)
    {
      indexColumn = KColA21[j];
      if(indexColumn == indexRow)
      {
        EntriesA21[j] = density*n1;
      }
      else
      {
        EntriesA21[j] = 0;
      }
    }
  }

  // modify matrix A22 for additional equations
  for(i=0;i<interface_dof;i++)
  {
    indexRow = interface_dof_Array[i];
    start = RowPtrA22[indexRow];
    end = RowPtrA22[indexRow+1];
    for (j=start;j<end;j++)
    {
      indexColumn = KColA22[j];
      if(indexColumn == indexRow)
      {
        EntriesA22[j] = density*n2;
      }
      else
      {
        EntriesA22[j] = 0;
      }
    }
  }

}
#endif

#ifdef __2D__
// =======================================================================
//
// Assemble2D_DG
//
// assembling for discontinuous Galerkin discretization
//
// =======================================================================

void Assemble2D_DG(CoeffFct2D *Coeff,int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters)
{
  double hK,w,integrant,tau_par;
  int N_AllMatrices = n_sqmatrices+n_matrices,out;
  int i,j,k,l,l1,l2,l3,n,n_neigh,m,r,q,dummy,N_UsedElements,N_LocalUsedElements,ii,jj,ll;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, ref_n;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TBaseFunct2D *bf;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunctCell;
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
  double **Matrices, *aux, *aux2 ;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *Coeffs[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5;
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
  double FunctionalValues[MaxN_BaseFunctions2D];
  int *EdgeDOF, N_EdgeDOF;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-12;
  double penetration_penalty, friction_parameter;
  double **JointValues, *JointValue, u1_values[3], u2_values[3];
  double delta;
  bool *SecondDer;

  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  RefTrans2D RefTrans;
  int N_DOF;
  double *Values;
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
  double xietaval_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double xideriv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double etaderiv_refNeigh1D[N_BaseFuncts2D][MaxN_QuadPoints_1D][MaxN_BaseFunctions2D];
  double *xyval_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *yderiv_refNeigh1D[MaxN_QuadPoints_1D];
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

#ifdef __3D__
  double z0, z1;
#endif
  out =  TDatabase::ParamDB->SC_VERBOSE;
  if(out>2)OutPut("DG "<<endl);

  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];
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
      [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
    Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
    for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
      Matrices[j] = aux+j*MaxN_BaseFunctions2D;

    LocMatrices = new double** [N_AllMatrices];
    for(i=0;i<N_AllMatrices;i++)
      LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
  }                                               // endif N_AllMatrices

  SecondDer = new bool[n_fespaces];
  SecondDer[0] = FALSE;

  if(n_rhs)
  {
    HangingRhs = new double* [n_rhs];
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();

      j = fespace->GetN_Hanging();
      HangingRhs[i] = new double [j];
      memset(HangingRhs[i], 0, SizeOfDouble*j);
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

  // The following arrays are used to store the values of the base functions and their derivatives for each vertex which is used to define the FE-Space.
  // It is only used to test the routine for errors. The index i has to be chosen according to the number of vertices of the used FE Space.
  for (i=0;i<N_BaseFuncts2D;i++)
  {
    value_basefunct_ref1D[i] = new double* [6];
    xderiv_basefunct_ref1D[i] = new double* [6];
    yderiv_basefunct_ref1D[i] = new double* [6];
    for (j=0;j<6;j++)
    {

      value_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D];
      xderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D];
      yderiv_basefunct_ref1D[i][j] = new double [MaxN_BaseFunctions2D];

      memset( value_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D );
      memset( xderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D );
      memset( yderiv_basefunct_ref1D[i][j] , 0 , sizeof(double)* MaxN_BaseFunctions2D );

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
        xietaval_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D];
        xideriv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D];
        etaderiv_ref1D[i][j][n] = new double [MaxN_BaseFunctions2D];

        memset( xietaval_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D );
        memset( xideriv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D );
        memset( etaderiv_ref1D[i][j][n] , 0 , sizeof(double)* MaxN_BaseFunctions2D );
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
  if(out>2)
  {
    OutPut("N_UsedElements:" << N_UsedElements << " " << endl); fflush(0);
    OutPut("N_BaseFuncts2D:" << N_BaseFuncts2D << " " << endl);
    OutPut("MaxN_QuadPoints_1D:" << MaxN_QuadPoints_1D << " " << endl);
    OutPut("MaxN_BaseFunctions2D:" << MaxN_BaseFunctions2D << " " << endl);
  };

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

        // Store values and derivatives of the base functions in vertices
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
    value_basefunct_ori[l] = new double[MaxN_BaseFunctions2D];
    xderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D];
    yderiv_basefunct_ori[l]  = new double[MaxN_BaseFunctions2D];
    value_basefunct_ori_neigh[l] = new double[MaxN_BaseFunctions2D];
    xderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D];
    yderiv_basefunct_ori_neigh[l]  = new double[MaxN_BaseFunctions2D];

  }

  for(m=0;m<4;m++)                                // arrays for coordinates, values and
  {                                               // determinant for 1D quadrature
    X1D[m] = new double[N_Points1D];              // coordinates of edge i
    Y1D[m] = new double[N_Points1D];
                                                  // determinant of affine mapping
    AbsDetjk1D[m] = new double[MaxN_QuadPoints_2D];
    for (j=0;j<N_Points1D;j++)                    // arrays for values in reference cell
    {
      xyval_ref1D[m][j] = new double[MaxN_BaseFunctions2D];
      xderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D];
      yderiv_ref1D[m][j] = new double[MaxN_BaseFunctions2D];
    }

  }                                               // endfor m

  for (j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################

  N_Parameters = Parameters->GetN_Parameters();   // get number of parameters of equation
  aux = new double [MaxN_QuadPoints_2D*N_Parameters];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Param[j] = aux + j*N_Parameters;

  // 20 <= number of term
  aux2 = new double [MaxN_QuadPoints_2D*20];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Coeffs[j] = aux2 + j*20;

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
      // prepare 1D quadrature formula
      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
      if(out>2){ OutPut("Polynomial degree on cell: " << l << endl);}
      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

      if(out>2)
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

      if(out>2)
      {
        for(r=0;r<N_Edges;r++)
        {
          cell->GetVertex(r)->GetCoords(x0, y0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
          if(out>2){OutPut("Local edge r: " << r << " Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
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
        if(out>2)
        {
          for (ii=0;ii<N_;ii++)                   //ii - test function
          {
            for(j=0;j<N_Points1D;j++)             // quadrature points
            {
              OutPut("basefunction: " << ii << " edge "<< r << " quadrature point " << j << " value " << xyval_ref1D[r][j][ii] << " xderiv " << xderiv_ref1D[r][j][ii] << " yderiv " << yderiv_ref1D[r][j][ii] << endl);
            }
          }
        }

      }                                           // endfor r

      // Values of the base functions in vertices. Used only for test purposes.
      TFEDatabase2D::GetOrigFromRef(RefTrans,ref_n,x_pos_ref,y_pos_ref,x_pos,y_pos,dummy2);
      for(l=0;l<ref_n;l++)
      {
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
        if(out>2){OutPut("x_pos_ref[l]: " << x_pos_ref[l] << " value_basefunct_ref1D[BaseFunctCell][l]: " << value_basefunct_ref1D[BaseFunctCell][l] <<  endl);};
      }

      for(r=0;r<N_Edges;r++)
      {                                           // For each edge, get the corresponding neighbour cell.
        neigh=cell->GetJoint(r)->GetNeighbour(cell);
        // If there is a neighbour to the edge, do...
        if(neigh)
        {                                         // Get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
          if(i<q)
          {

            // calculate all needed derivatives of this FE function
                                                  // finite element on neighbour
            CurrEleNeigh = fespaces[n]->GetFE2D(q,neigh);
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
            {                                     //for (k=0; k< MaxN_BaseFunctions2D; k++)
              // OutPut(" xi1D " << xi1D[BaseFunctCell][r][j] << " eta1D " << eta1D[BaseFunctCell][r][j] << " BaseFunct: " << BaseFunctCell << "\n");
              // OutPut(" xi1D " << xi1D[BaseFunctNeigh][neigh_edge][j] << " eta1D " << eta1D[BaseFunctCell][r][j] << " BaseFunctNeigh: " << BaseFunctNeigh << "\n");
              if(out>2){  OutPut("X1D[r][j]: " << X1D[r][j] << " Y1D[r][j]: " << Y1D[r][j] << " X1D_neigh[neigh_edge][j] " <<  X1D_neigh[neigh_edge][j] << "  Y1D_neigh[neigh_edge][j]: " <<  Y1D_neigh[neigh_edge][j] <<  endl);}
            }

            if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0])
            {
              if(out>2)
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

              if(out>2)
              {
                for (jj=0;jj<N_Neigh;jj++)
                {
                  for(j=0;j<N_Points1D;j++)       // quadrature points
                  {
                    OutPut("basefunction: " << jj << " edge " << r << " quadrature point " << j << " value " << xyval_ref1D[j][jj] << " xderiv " << xderiv_ref1D[j][jj] << " yderiv " << yderiv_ref1D[j][jj] << endl);
                  }
                }
              }

            }                                     //endif
            else
            {
              if(out>2)
              {
                OutPut("Inverse the order of the quadrature points on neighbour edge !" << endl);
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
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    OutPut("Basisfkt: "<< DOF[k] <<"  (x,y)-coordinate: " << x_pos[l] << " " << y_pos[l] << " value: " <<  value_basefunct_ori[l][k] << endl);
                  }
                }
              }
            }                                     //endfor k

            for(k=0;k<N_Neigh; k++)
            {
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    OutPut("Basisfkt neigh: "<< DOF_neigh[k] <<"  (x,y)-coordinate: " << x_pos_neigh[l] << " " << y_pos_neigh[l] << " value: " <<  value_basefunct_ori_neigh[l][k] << endl);
                  }
                }
              }
            }                                     //endfor k

            // #################################################################################
            // Compute the edge integrals for discontinuous Galerkin
            // #################################################################################

            //[0][4]: " << Coeffs[0][3]);
            // get vertices of boundary edge
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            if(out>2){  OutPut("Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
            // compute length of the boundary edge
            hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            // tangential normal vector to this boundary (normalized)
            tx = (x1-x0)/hE;
            ty = (y1-y0)/hE;
            tau_par = TDatabase::ParamDB->FACE_SIGMA;
            if(out>2)OutPut(" tau edge stabilization: " << tau_par << endl);

            Entries1 = sqmatrices[n]->GetEntries();
            RowPtr1 = sqmatrices[n]->GetRowPtr();
            ColInd1 = sqmatrices[n]->GetKCol();
            fespace = sqmatrices[n]->GetFESpace();
            ActiveBound = fespace->GetActiveBound();
            if(out>2)
            {
              OutPut("ActiveBound von sqmatrix " << n << "  :" << ActiveBound  << endl);
              for(j=0;j<N_Points1D; j++)
              {
                OutPut("Det["<<j<<"]:" <<  hE/2 <<" b1: "<< Coeffs[j][1] <<" b2: " << Coeffs[j][2] << endl);
              }
            }

            //edge integrals: test function from cell <-> ansatz function from cell
            if(out>2){  OutPut("testfct. cell<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              //OutPut("dof_ii" << dof_ii << endl);
              // Dirichlet node

              if (dof_ii>=ActiveBound)continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral; Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
                {
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += 0.5*Coeffs[j][0] * (-xyval_ref1D[r][j][ii]*(xderiv_ref1D[r][j][jj]*nx + yderiv_ref1D[r][j][jj]*ny) - xyval_ref1D[r][j][jj]*(xderiv_ref1D[r][j][ii]*nx + yderiv_ref1D[r][j][ii]*ny));

                                                  // penalty term to achieve coercitivity
                  integrant += tau_par*Coeffs[j][0]*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) < 0) integrant -= (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                  w = weights1D[j]*hE/2;

                  integral += w*integrant;        // integral on the edge
                }
                // update matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from cell <-> ansatz function from neighbour
            if(out>2){  OutPut("testfct. cell<-> ansatzfct. neighbour Integral:" << endl);}
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
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += 0.5*Coeffs[j][0] * (-xyval_ref1D[r][j][ii]*(xderiv_refNeigh1D[j][jj]*nx + yderiv_refNeigh1D[j][jj]*ny) + xyval_refNeigh1D[j][jj]*(xderiv_ref1D[r][j][ii]*nx + yderiv_ref1D[r][j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += -tau_par*Coeffs[j][0]*xyval_ref1D[r][j][ii]*xyval_refNeigh1D[j][jj];
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) < 0) integrant += (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_refNeigh1D[j][jj]*xyval_ref1D[r][j][ii];
                  w = weights1D[j]*hE/2;
                  integral += w*integrant;        // integral on the edge                  // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from cell
            if(out>2){ OutPut("testfct. neighbour<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound)
                continue;
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];

                if(dummy == 0)                    //
                  // OutPut("jj " << dof_jj << endl);
                  // initialize the boundary integrals
                  integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += 0.5*Coeffs[j][0] * (xyval_refNeigh1D[j][ii]*(xderiv_ref1D[r][j][jj]*nx + yderiv_ref1D[r][j][jj]*ny) - xyval_ref1D[r][j][jj]*(xderiv_refNeigh1D[j][ii]*nx + yderiv_refNeigh1D[j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += -tau_par*Coeffs[j][0]*xyval_refNeigh1D[j][ii]*xyval_ref1D[r][j][jj];
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) > 0) integrant -= (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_ref1D[r][j][jj]*xyval_refNeigh1D[j][ii];
                  w = weights1D[j]*hE/2;
                  integral += w*integrant;        // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from neighbour
            if(out>2){     OutPut("testfct. neighbour <-> ansatzfct. neighbour Integral:" << endl);}
            for (ii=0;ii<N_Neigh;ii++)            //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF_neigh[ii];

              // Dirichlet node
              if (dof_ii>=ActiveBound) continue;
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
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += 0.5*Coeffs[j][0] * (xyval_refNeigh1D[j][ii]*(xderiv_refNeigh1D[j][jj]*nx + yderiv_refNeigh1D[j][jj]*ny) + xyval_refNeigh1D[j][jj]*(xderiv_refNeigh1D[j][ii]*nx + yderiv_refNeigh1D[j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += tau_par*Coeffs[j][0]*xyval_refNeigh1D[j][ii]*xyval_refNeigh1D[j][jj];
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) > 0) integrant += (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_refNeigh1D[j][ii]*xyval_refNeigh1D[j][jj];
                  w = weights1D[j]*hE/2;
                  integral += w*integrant;        // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
        else
        {
#ifdef __3D__
          cell->GetVertex(r)->GetCoords(x0, y0, z0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
          cell->GetVertex(r)->GetCoords(x0, y0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
          if(out>2){  OutPut("Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
          // compute length of the boundary edge
          hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
          // compute normal vector to this boundary (normalized)
          nx = (y1-y0)/hE;
          ny = (x0-x1)/hE;
          // tangential normal vector to this boundary (normalized)
          tx = (x1-x0)/hE;
          ty = (y1-y0)/hE;
          tau_par = TDatabase::ParamDB->FACE_SIGMA;
          if(out>2)OutPut(" tau edge stabilization: " << tau_par << endl);

          Entries1 = sqmatrices[n]->GetEntries();
          RowPtr1 = sqmatrices[n]->GetRowPtr();
          ColInd1 = sqmatrices[n]->GetKCol();
          fespace = sqmatrices[n]->GetFESpace();
          ActiveBound = fespace->GetActiveBound();
          if(out>2)
          {
            OutPut("ActiveBound von sqmatrix " << n << "  :" << ActiveBound  << endl);
            for(j=0;j<N_Points1D; j++)
            {
              OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<< Coeffs[j][1] <<" b2: " << Coeffs[j][2] << endl);
            }
          }

          //edge integrals: test function from cell <-> ansatz function from cell
          if(out>2){  OutPut("testfct. cell<-> ansatzfct. cell Integral:" << endl);}
          for (ii=0;ii<N_;ii++)                   //ii - test function
          {
            // look for 'ii'-th row in all matrices
            dof_ii = DOF[ii];
            //OutPut("dof_ii" << dof_ii << endl);
            // Dirichlet node

            if (dof_ii>=ActiveBound)continue;

            // for all dof in the mesh cell
            // jj - ansatz function
            for (jj=0;jj<N_;jj++)
            {
              dof_jj = DOF[jj];
              // initialize the boundary integrals
              integral=0;
              for (j=0;j<N_Points1D;j++)          // compute edge integral; Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
              {
                integrant = 0;
                                                  // contribution of diffusive term
                integrant += Coeffs[j][0] * (-xyval_ref1D[r][j][ii]*(xderiv_ref1D[r][j][jj]*nx + yderiv_ref1D[r][j][jj]*ny) - xyval_ref1D[r][j][jj]*(xderiv_ref1D[r][j][ii]*nx + yderiv_ref1D[r][j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                integrant += 2*tau_par/hE*Coeffs[j][0]*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                                                  // contribution of convective term
                if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) < 0) integrant -= (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                w = weights1D[j]*hE/2;
                integral += w*integrant;          // integral on the edge
              }

              // update matrix
              found = 0;
              for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
              {
                if (ColInd1[ll] == dof_jj)
                {
                  if(out>2)
                  {
                    OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                      <<  " DOF Ansatzfkt: " << dof_jj << endl);
                  }

                  Entries1[ll] += integral;
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
            }                                     // end inner loop over dof (jj)
          }                                       // end outer loop over dof (ii)
          if (out>2){  OutPut(endl);}
        }                                         //end else for if(neigh)
      }                                           //endfor r (loop over edges)
    }                                             // endfor n (loop over sqmatrices

    if(out>2)
    {
      OutPut("weak boundary conditions" << endl);
      OutPut("cell: " << i << endl);
    }

    for(n=0;n<n_rhs;n++)
    {
      fespace = ferhs[n];
      ActiveBound = fespace->GetActiveBound();
      CurrentElement = fespace->GetFE2D(i,cell);  // finite element on cell

      LocalUsedElements[0] = CurrentElement;
      
      N_ = N_BaseFunct[CurrentElement];
      if(out>2){OutPut("N_BaseFunct " << N_ << endl);}

      local_rhs = righthand+n*MaxN_BaseFunctions2D;
      RHS = rhs[n];
      CurrentHangingRhs = HangingRhs[n];
      // find space for this linear form

      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[n] + RhsBeginIndex[n][i];
      if(out>2)
      {
        for(k=0; k<N_; k++)
        {
          OutPut("DOF[k]: k: " << k << " DOF: " << DOF[k]  << endl);
        }
      }

      BoundaryCondition = BoundaryConditions[n];
      BoundaryValue = BoundaryValues[n];

      N_Joints = cell->GetN_Edges();

      if(out>2)OutPut("N_Joints: " << N_Joints << endl);

      RefTrans = TFEDatabase2D::GetOrig(1, LocalUsedElements,
        Coll, cell, SecondDer,
        N_Points, xi, eta, weights, X, Y, AbsDetjk);

      if(N_Parameters>0)                          // get parameters of equ.
        Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);
                                                  // get coefficients of pde
      if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);

      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
          joint->GetType() == IsoBoundEdge)
        {
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint)
          {
            boundedge = (TBoundEdge *)joint;
            BoundComp = boundedge->GetBoundComp();
            boundedge->GetParameters(t0, t1);
          }
          else
          {
            isoboundedge = (TIsoBoundEdge *)joint;
            BoundComp = isoboundedge->GetBoundComp();
            isoboundedge->GetParameters(t0, t1);
          }
          // get id of the boundary component
          comp=BoundComp->GetID();
          if(out>2)OutPut("boundary comp: " << comp << endl);
          // get type of the boundary condition at the beginning
          // and at the end of the current edge
          if (t0 < t1)
          {
            BoundaryCondition(comp, t0+eps, Cond0);
            BoundaryCondition(comp, t1-eps, Cond1);
          }
          else
          {
            BoundaryCondition(comp, t0-eps, Cond0);
            BoundaryCondition(comp, t1+eps, Cond1);
          }
          // only one boundary condition per edge allowed
          if(Cond0 == Cond1)
          {
            if(Cond0 == DIRICHLET)
            {
              // get polynomial degree of fe
              if(out>2){OutPut("Edge " << m << endl);}
              l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                (CurrentElement);
              // get a suitable line quadrature formula
              LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
              qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
              if(out>2){OutPut("QuadFormula " << qf1D << endl);}
              qf1D->GetFormulaData(N_LinePoints, LineWeights, zeta);
              TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                ->MakeRefElementData(LineQuadFormula);

              TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][m],
                eta1D[BaseFunctCell][m],
                X1D[m], Y1D[m], AbsDetjk1D[m]);

              for(j=0;j<N_Points1D;j++)           // get values and derivatives in original cell
              {
                TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctCell][m][j],
                  eta1D[BaseFunctCell][m][j],
                  TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
                  Coll, (TGridCell *)cell,
                  xietaval_ref1D[BaseFunctCell][m][j],
                  xideriv_ref1D[BaseFunctCell][m][j],
                  etaderiv_ref1D[BaseFunctCell][m][j],
                  xyval_ref1D[m][j],
                  xderiv_ref1D[m][j],
                  yderiv_ref1D[m][j]);
              }

              // get vertices of boundary edge
#ifdef __3D__
              cell->GetVertex(m)->GetCoords(x0, y0, z0);
              cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
              cell->GetVertex(m)->GetCoords(x0, y0);
              cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
              // compute the length of the boundary edge
              hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
              nx = (y1-y0)/hE;
              ny = (x0-x1)/hE;
              // tangential normal vector to this boundary (normalized)
              tx = (x1-x0)/hE;
              ty = (y1-y0)/hE;
              tau_par = TDatabase::ParamDB->FACE_SIGMA;

              // compute boundary integral: ansatz function from cell with Dirichlet boundary value
              for(k=0;k<N_;k++)
              {
                l3 = DOF[k];
                integral=0;
                for(l=0;l<N_LinePoints;l++)
                {
                  if(out>2){OutPut("l:" << l << " normal vector: " << nx <<  " " << ny << endl);};
                  // get quadrature point on the boundary
                  t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                  if(out>2){OutPut("Quad formula data: Quadrature point[l]:" << zeta[l] <<  " weight: " << LineWeights[l] << endl);};
                  // get value in this quadrature point (in s)
                  BoundaryValue(comp, t, s);
                  if(out>2)
                  {
                    OutPut("Parameter Quadrature Point: " << t << " Boundary Value: " << s << endl);
                    OutPut("Function values Quadrature Point: " << k << " xyval " << xyval_ref1D[m][l][k] << " xderiv " << xderiv_ref1D[m][l][k] << " yderiv " << yderiv_ref1D[m][l][k] << endl);

                  };
                  // multiply value with weights from quadrature formula
                  // and determinant from integral transformation to the
                  // unit edge (-1,1)
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += Coeffs[j][0]*(-s*(xderiv_ref1D[m][l][k]*nx + yderiv_ref1D[m][l][k]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += 2*tau_par/hE*Coeffs[j][0]*s*xyval_ref1D[m][l][k];
                  if((Coeffs[l][1]*nx + Coeffs[l][2]*ny) < 0) integrant -= (Coeffs[l][1]*nx + Coeffs[l][2]*ny)*s*xyval_ref1D[m][l][k];
                  integral += hE/2 * LineWeights[l]*integrant;
                }
                if(l3 < ActiveBound) RHS[l3] += integral;
                else
                {
                  OutPut("Index l3 " << l3 << " bigger than ActiveBound " << ActiveBound << " !");
                  exit(4711);
                }
                if(out>2) OutPut("DOF: " << l3 << " integral: " << integral << " RHS: " << RHS[l3] << endl);
              }
            }                                     // endif
          }                                       // endif (Cond0==Cond1)
          else
          {
            OutPut("different boundary condition on one edge ");
            OutPut("are not allowed!" << endl);
            exit(4711);
          }
        }                                         // endif (boundary joint)
      }                                           // endfor m (N_Joints)
    }                                             // endfor n (rhs)
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

  if(n_rhs)
  {
    for(i=0;i<n_rhs;i++)
      delete HangingRhs[i];
    delete HangingRhs;

    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }

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
  delete SecondDer;

  delete UsedElements;
  for(i=0; i < N_BaseFuncts2D; i++)
  {
    for(j=0; j < 4; j++)
    {
      for(m=0; m < MaxN_QuadPoints_1D; m++)
      {
        delete xietaval_ref1D[i][j][m];
        delete xideriv_ref1D[i][j][m];
        delete etaderiv_ref1D[i][j][m];
      }
      delete xietaval_ref1D[i][j];
      delete xideriv_ref1D[i][j];
      delete etaderiv_ref1D[i][j];
    }
    delete xietaval_ref1D[i];
    delete xideriv_ref1D[i];
    delete  etaderiv_ref1D[i];
  }

  delete xietaval_ref1D;
  delete xideriv_ref1D;
  delete etaderiv_ref1D;

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
  }
  
  for(l=0;l<ref_n;l++)
  {
    delete value_basefunct_ori[l];
    delete xderiv_basefunct_ori[l];
    delete yderiv_basefunct_ori[l];
    delete value_basefunct_ori_neigh[l];
    delete xderiv_basefunct_ori_neigh[l];
    delete yderiv_basefunct_ori_neigh[l];
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
    for (j=0;j<ref_n;j++)
    {
      delete value_basefunct_ref1D[i][j];
      delete xderiv_basefunct_ref1D[i][j];
      delete yderiv_basefunct_ref1D[i][j];
    }
    delete value_basefunct_ref1D[i];
    delete xderiv_basefunct_ref1D[i];
    delete yderiv_basefunct_ref1D[i];
  }

  delete value_basefunct_ref1D;
  delete xderiv_basefunct_ref1D;
  delete yderiv_basefunct_ref1D;


  /*for (n=0;n<N_BaseFuncts2D;n++)
    { delete xi1DNeigh[n];
      delete eta1DNeigh [n];
    }*/

  for (i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete  xderiv_refNeigh1D[i];
    delete   yderiv_refNeigh1D[i];
  }
  /*        delete xderiv_Neigh1D;
          delete  yderiv_Neigh1D;
          delete  xyval_Neigh1D;*/
  delete aux;

  int N_Rows;
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  if(out>2)
  {
    for(k=0;k<n_sqmatrices;k++)
    {
      cout << endl;
      cout << "sqmatrix: " << k << endl;
      RowPtr = sqmatrices[k]->GetRowPtr();
      Entries = sqmatrices[k]->GetEntries();
      ColInd = sqmatrices[k]->GetKCol();
      N_Rows = sqmatrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=RowPtr[i+1];
        for(j=RowPtr[i];j<end;j++)
        {
          // cout << j << endl;
          cout << "Matrix: " << setw(5) << i << setw(5) << ColInd[j] << "   ";
          cout << setw(10) << Entries[j] << endl;
        }
      }
      cout << endl;
    }                                             // endfor k

    for(k=0;k<n_rhs;k++)
    {
      cout << "rhs: " << k << endl;
      N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
      RHS=rhs[k];
      for(i=0;i<N_Rows;i++)
        cout << setw(5) << i << setw(20) << RHS[i] << endl;
    }
  }                                               //endif

  /*  for(k=0;k<n_matrices;k++)
    {
      cout << endl;
      cout << "matrix: " << k << endl;
      RowPtr = matrices[k]->GetRowPtr();
      Entries = matrices[k]->GetEntries();
      ColInd = matrices[k]->GetKCol();
      N_Rows = matrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k */

} // end of Assemble2D_DG


// =======================================================================
//
// Assemble2D_CIP
//
// assembling for continuous interior penalty discretization
//
// =======================================================================

void Assemble2D_CIP(CoeffFct2D *Coeff,int n_fespaces, TFESpace2D **fespaces,
int n_sqmatrices, TSquareMatrix2D **sqmatrices,
int n_matrices, TMatrix2D **matrices,
int n_rhs, double **rhs, TFESpace2D **ferhs,
BoundCondFunct2D **BoundaryConditions,
BoundValueFunct2D **BoundaryValues,
TAuxParam2D *Parameters)
{
  const int MaxN_BaseFunctions2D_Ersatz =100;

  double hK,w,integrant,tau_par,sigma_par;
  int N_AllMatrices = n_sqmatrices+n_matrices,out;
  int i,j,k,l,l1,l2,l3,n,n_neigh,m,r,q,dummy,N_UsedElements,N_LocalUsedElements,ii,jj,ll,weak;
  int N_Cells, N_Points, N_Parameters, N_Points1D, N_Edges, N_, N_Hanging;
  int N_Test, N_Ansatz, N_Joints, ref_n;
  int Used[N_FEs2D];
  int *N_BaseFunct;
  BaseFunct2D *BaseFuncts;
  TBaseFunct2D *bf;
  TFESpace2D *fespace;
  FE2D *UsedElements, LocalUsedElements[N_FEs2D], CurrentElement;
  FE2D TestElement, AnsatzElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1D;
  BaseFunct2D BaseFunctCell;
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
  double **Matrices, *aux, *aux2, *aux3, *aux4;
  double **Matrix;
  double ***LocMatrices, **LocRhs;
  int LocN_BF[N_BaseFuncts2D];
  BaseFunct2D LocBF[N_BaseFuncts2D];
  double *Coeffs[MaxN_QuadPoints_2D];
  int *DOF, ActiveBound, DirichletBound, end, last;
  int *TestDOF, *AnsatzDOF;
  double *Entries,*Entries1,*Entries2,*Entries3, *Entries4, *Entries5;
  int *ColInd, *RowPtr;
  int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
  int *ColInd4, *RowPtr4, *ColInd5, *RowPtr5;
  double *RHS, *MatrixRow;
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

  double *Coefficients1D[MaxN_QuadPoints_2D];
  double *Parameters1D[MaxN_QuadPoints_2D];

  double xi1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D], eta1D[N_BaseFuncts2D][4][MaxN_QuadPoints_1D];
  double**** xietaval_ref1D = new double*** [N_BaseFuncts2D];
  double**** xideriv_ref1D = new double*** [N_BaseFuncts2D];
  double**** etaderiv_ref1D = new double*** [N_BaseFuncts2D];
  double *xyval_ref1D[4][MaxN_QuadPoints_1D];
  double *xderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *yderiv_ref1D[4][MaxN_QuadPoints_1D];
  double *X1D[4], *Y1D[4], *X1D_neigh[4], *Y1D_neigh[4];
  RefTrans2D RefTrans;
  int N_DOF;
  double *Values;
  double*** value_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** xderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double*** yderiv_basefunct_ref1D = new double** [N_BaseFuncts2D];
  double *value_basefunct_ori[6];
  double *xderiv_basefunct_ori[6];
  double *yderiv_basefunct_ori[6];
  double x_pos_ref[6];
  double y_pos_ref[6];
  double x_pos[6];
  double y_pos[6];
  double *value_basefunct_ori_neigh[6];
  double *xderiv_basefunct_ori_neigh[6];
  double *yderiv_basefunct_ori_neigh[6];
  double x_pos_neigh[6];
  double y_pos_neigh[6];
  double dummy2[6];

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
  double *xderiv_Neigh1D, *yderiv_Neigh1D, *xyval_Neigh1D;

  double jump_xyval[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_xderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];
  double jump_yderiv[MaxN_QuadPoints_1D][2*N_BaseFuncts2D];

#ifdef __3D__
  double z0, z1;
#endif

  out=1;
  
  if(TDatabase::ParamDB->SC_VERBOSE > 1)
    OutPut("CIP " << TDatabase::TimeDB->CURRENTTIME << endl);
  weak = TDatabase::ParamDB->WEAK_BC;
  if (weak>=1)
    TDatabase::ParamDB->INTERNAL_DO_NOT_RESPECT_DIRICHLET_BC = 0;


  // ########################################################################
  // store information in local arrays
  // ########################################################################
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  if(n_sqmatrices)
  {
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    for(i=0;i<n_sqmatrices;i++)
    {
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
    }                                             // endfor
  }                                               // endif n_sqmatrices

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

  if(n_rhs)
  {
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    for(i=0;i<n_rhs;i++)
    {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
    }                                             // endfor

    LocRhs = new double* [n_rhs];
    righthand = new double [n_rhs*MaxN_BaseFunctions2D];
    for(i=0;i<n_rhs;i++)
      LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
  }                                               // endif n_rhs

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
  if (out==2)
      OutPut("N_UsedElements:" << N_UsedElements << " " << endl); fflush(0);
  // OutPut("N_BaseFuncts2D:" << N_BaseFuncts2D << " " << endl);
  // OutPut("MaxN_QuadPoints_1D:" << MaxN_QuadPoints_1D << " " << endl);
  // OutPut("MaxN_BaseFunctions2D_Ersatz:" << MaxN_BaseFunctions2D_Ersatz << " " << endl);

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
        
        ref_n = 6;

        for (j=0;j<N_Points1D;j++)                // for all quadrature poin
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
  if (out==2)
      OutPut("basefunct" << endl);

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
    }

  }                                               // endfor m

  for (j=0;j<N_Points1D;j++)                      // arrays for values in reference cell
  {
    xyval_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    xderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
    yderiv_refNeigh1D[j] = new double[MaxN_BaseFunctions2D_Ersatz];
  }

  // ########################################################################
  // Arrays for Parameters
  // ########################################################################

  if (out==2)
      OutPut("parameters" << endl);
  N_Parameters = Parameters->GetN_Parameters();   // get number of parameters of equation
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

    if (out==2)
	OutPut("cell " << i << endl);
    for(n=0;n<n_sqmatrices;n++)
    {
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
      // prepare 1D quadrature formula
      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
      if(out>2){ OutPut("Polynomial degree on cell: " << l << endl);}
      LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1D->GetFormulaData(N_Points1D, weights1D, zeta);

      BoundaryCondition = BoundaryConditions[n];
      BoundaryValue = BoundaryValues[n];
      
      if(out>2)
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

      if(out>2)
      {
        for(r=0;r<N_Edges;r++)
        {
          cell->GetVertex(r)->GetCoords(x0, y0);
          cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
          if(out>2){OutPut("Local edge r: " << r << " Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
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
        // OutPut("Hallo: x_pos_ref[l]: " << x_pos_ref[l] << "value_basefunct_ref1D[BaseFunctCell][l]: " << value_basefunct_ref1D[BaseFunctCell][l] << endl);
      }

      for(r=0;r<N_Edges;r++)
      {                                           // For each edge, get the corresponding neighbour cell.
        neigh=cell->GetJoint(r)->GetNeighbour(cell);
        //#######################################################################//
        // get coefficients on edges
        //only implemented for coeffs that do not depend on the params
        //#######################################################################//

        //  if(N_Parameters>0)                // get parameters of equ.
        // Parameters->GetParameters(N_Points1D, Coll, cell, i, xi1D[BaseFunctCell][r], eta1D[BaseFunctCell][r], X1D[r], Y1D[r], Param1D);

        if(Coeff) Coeff(N_Points1D, X1D[r], Y1D[r], Param, Coefficients1D);
        //#######################################################################//
        // If there is a neighbour to the edge, do...
        if(neigh)
        {                                         // Get the number of this neigbbour cell from the clipboard
          q = neigh->GetClipBoard();
          if(i<q)
          {

            // calculate all needed derivatives of this FE function
                                                  // finite element on neighbour
            CurrEleNeigh = fespaces[n]->GetFE2D(q,neigh);
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
            {                                     //for (k=0; k< MaxN_BaseFunctions2D_Ersatz; k++)
              // OutPut(" xi1D " << xi1D[BaseFunctCell][r][j] << " eta1D " << eta1D[BaseFunctCell][r][j] << " BaseFunct: " << BaseFunctCell << "\n");
              // OutPut(" xi1D " << xi1D[BaseFunctNeigh][neigh_edge][j] << " eta1D " << eta1D[BaseFunctCell][r][j] << " BaseFunctNeigh: " << BaseFunctNeigh << "\n");
              if(out>2){  OutPut("X1D[r][j]: " << X1D[r][j] << " Y1D[r][j]: "  <<  Y1D[r][j] <<  " X1D[neigh_edge][j] " << X1D_neigh[neigh_edge][j] << " Y1D[neigh_edge][j]: " <<  Y1D_neigh[neigh_edge][j] <<  endl);}
            }

            if(X1D_neigh[neigh_edge][0] == X1D[r][0] && Y1D_neigh[neigh_edge][0] == Y1D[r][0] )
            {
              if(out>2)
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
              if(out>2)
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
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    OutPut("Basisfkt: "<< DOF[k] <<"  (x,y)-coordinate: " << x_pos[l] << " " << y_pos[l] << " value: " <<  value_basefunct_ori[l][k] << endl);
                  }
                }
              }
            }                                     //endfor k

            for(k=0;k<N_Neigh; k++)
            {
              if (out>2)
              {
                for(l=0;l<ref_n;l++)
                {
                  if(out>2)
                  {
                    OutPut("Basisfkt neigh: "<< DOF_neigh[k] <<"  (x,y)-coordinate: " << x_pos_neigh[l] << " " << y_pos_neigh[l] << " value: " <<  value_basefunct_ori_neigh[l][k] << endl);
                  }
                }
              }
            }                                     //endfor k

            //Compute the jumps of the basis functions
            //and of their derivatives in the quadrature points on edge r
            //First for the basis functions of cell i

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
                  if(out>2)
                  {
                    OutPut("xyval_Zelle: "<< xyval_ref1D[r][j][k] <<" xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("xderiv_Zelle: "<< xderiv_ref1D[r][j][k] <<  " xderiv_Neigh: " << xderiv_refNeigh1D[j][l] << " jump= " << jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("yderiv_Zelle: "<< yderiv_ref1D[r][j][k] <<  " yderiv_Neigh: " << yderiv_refNeigh1D[j][l] << " jump= " << jump_yderiv[j][k] << " of basefunction: "<< DOF[k]  << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut(endl);
                  }
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

                  if(out>2)
                  {
                    OutPut(" No Neighbour: xyval_Zelle: "<< xyval_ref1D[r][j][k] << " jump= " <<  jump_xyval[j][k] <<  " of basefunction: "<< DOF[k] << " in cell: " << i << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("No Neighbour: x-deriv-jump= " <<  jump_xderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("No Neighbour: y-deriv-jump= " <<  jump_yderiv[j][k] << " of basefunction: "<< DOF[k] << " in cell: " << i <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
                    OutPut(endl);
                  }
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
                  if(out>2)
                  {
                    OutPut("Neighbour!!" << "xyval_Neigh: " << xyval_refNeigh1D[j][l] << " jump= " <<  jump_xyval[j][l+N_]<<  " of basefunction: "<< DOF_neigh[l] << " in cell: " << q << " on edge (local): " << r << " in quadrature point: " << j << endl);
                    OutPut("Neighbour!! " << "x-deriv-jump: " << jump_xderiv[j][l+N_] << " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << endl);
                    OutPut("Neighbour!! y-deriv-jump= " <<  jump_yderiv[j][l+N_]<< " of basefunction: "<< DOF_neigh[l] << " in cell: " << q <<  " on edge: " << r << " in quadrature point: " << j << "\n" << endl);
                    OutPut(endl);
                  }
                }                                 //endfor j
              }                                   //endif
            }                                     //endfor l

            // #################################################################################
            // Compute the edge integrals with the jumps of the basis functions and their derivatives
            // #################################################################################

            //[0][4]: " << Coeffs[0][3]);
            // get vertices of boundary edge
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            if(out>2){  OutPut("Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
            // compute length of the boundary edge
            hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            // tangential normal vector to this boundary (normalized)
            tx = (x1-x0)/hE;
            ty = (y1-y0)/hE;
            tau_par = TDatabase::ParamDB->FACE_SIGMA;
            tau_par = tau_par*hE*hE;
            if (out>2){ OutPut(" tau edge stabilization: " << tau_par << endl)};

            Entries1 = sqmatrices[n]->GetEntries();
            RowPtr1 = sqmatrices[n]->GetRowPtr();
            ColInd1 = sqmatrices[n]->GetKCol();
            fespace = sqmatrices[n]->GetFESpace();
            ActiveBound = fespace->GetActiveBound();
            if(out>2)
            {
              OutPut("ActiveBound von sqmatrix " << n << "  :" << ActiveBound  << endl);
              for(j=0;j<N_Points1D; j++)
              {
                OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<<Coefficients1D[j][1] <<" b2: " << Coefficients1D[j][2] << endl);
              }
            }
            //edge integrals: test function from cell <-> ansatz function from cell
            if(out>2){  OutPut("testfct. cell<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              //OutPut("dof_ii" << dof_ii << endl);
              // Dirichlet node

              if (dof_ii>=ActiveBound)continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  // distinguish different types of the cip method
                  // 0 : Burman, Hansbo 2004, eq. (3)
                  // 1:  Burman, Hansbo 2004, eq. (5), without orthogonal term, ES in Section 3
                  // 2: do nothing
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*fabs(nx*Coefficients1D[j][1] + ny*Coefficients1D[j][2])*
                        (jump_xderiv[j][jj]*nx + jump_yderiv[j][jj]*ny) *
                        (jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj]*Coefficients1D[j][1] + jump_yderiv[j][jj]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii]*Coefficients1D[j][1] + jump_yderiv[j][ii]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }

                  w = weights1D[j]*hE/2;
                  integral += w*integrant;        // integral on the edge
                }

                // update matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from cell <-> ansatz function from neighbour
            if(out>2){  OutPut("testfct. cell<-> ansatzfct. neighbour Integral:" << endl);}
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
                  // distinguish different types of the cip method
                  // 0 : Burman, Hansbo 2004, eq. (3)
                  // 1:  Burman, Hansbo 2004, eq. (5), without orthogonal term, ES in Section 3
                  // 2: do nothing
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*fabs(nx*Coefficients1D[j][1] + ny*Coefficients1D[j][2]) *
                        (jump_xderiv[j][jj+N_]*nx + jump_yderiv[j][jj+N_]*ny) *
                        (jump_xderiv[j][ii]*nx + jump_yderiv[j][ii]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj+N_]*Coefficients1D[j][1] + jump_yderiv[j][jj+N_]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii]*Coefficients1D[j][1] + jump_yderiv[j][ii]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }
                  w = weights1D[j]* hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from cell
            if(out>2){ OutPut("testfct. neighbour<-> ansatzfct. cell Integral:" << endl);}
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

                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral
                {
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*fabs(nx*Coefficients1D[j][1] + ny*Coefficients1D[j][2])*
                        (jump_xderiv[j][jj]*nx + jump_yderiv[j][jj]*ny) *
                        (jump_xderiv[j][ii+N_]*nx + jump_yderiv[j][ii+N_]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj]*Coefficients1D[j][1] + jump_yderiv[j][jj]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii+N_]*Coefficients1D[j][1] + jump_yderiv[j][ii+N_]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }
                  w = weights1D[j]* hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}

            //edge integrals: test function from neighbour <-> ansatz function from neighbour
            if(out>2){     OutPut("testfct. neighbour <-> ansatzfct. neighbour Integral:" << endl);}
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
                  switch (TDatabase::ParamDB->CIP_TYPE)
                  {
                    case 0:
                      integrant = tau_par*fabs(nx*Coefficients1D[j][1] + ny*Coefficients1D[j][2])*
                        (jump_xderiv[j][jj+N_]*nx + jump_yderiv[j][jj+N_]*ny) *
                        (jump_xderiv[j][ii+N_]*nx + jump_yderiv[j][ii+N_]*ny);
                      break;
                    case 1:
                      integrant = tau_par*(jump_xderiv[j][jj+N_]*Coefficients1D[j][1] + jump_yderiv[j][jj+N_]*Coefficients1D[j][2]) *
                        (jump_xderiv[j][ii+N_]*Coefficients1D[j][1] + jump_yderiv[j][ii+N_]*Coefficients1D[j][2]);
                      break;
                    case 2:
                      integrant = 0;
                      break;
                  }

                  w = weights1D[j]* hE/2;
                  integral+= w*integrant;         // integral on the edge
                }
                // update first matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF testfct: " << dof_ii
                        <<  " DOF ansatzfct: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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

        else
        {
          weak = TDatabase::ParamDB->WEAK_BC;
          if(weak==1)
          {
            OutPut("weak" << r << endl);
            joint = cell->GetJoint(r);
            if(joint->GetType() == BoundaryEdge ||
              joint->GetType() == IsoBoundEdge)
            {
              if(joint->GetType() == BoundaryEdge)
              {
                boundedge = (TBoundEdge *)joint;
                BoundComp = boundedge->GetBoundComp();
                boundedge->GetParameters(t0, t1);
              }
              else
              {
                isoboundedge = (TIsoBoundEdge *)joint;
                BoundComp = isoboundedge->GetBoundComp();
                isoboundedge->GetParameters(t0, t1);
              }
              // get id of the boundary component
              OutPut("weak" << r << endl);
              comp=BoundComp->GetID();
              OutPut("weak" << comp << endl);
              // get type of the boundary condition at the beginning
              // and at the end of the current edge
              if (t0 < t1)
              {
                BoundaryCondition(comp, t0+eps, Cond0);
                BoundaryCondition(comp, t1-eps, Cond1);
              }
              else
              {
                BoundaryCondition(comp, t0-eps, Cond0);
                BoundaryCondition(comp, t1+eps, Cond1);
              }
              OutPut("weak" << r << endl);
              // only one boundary condition per edge allowed
              if(Cond0 == Cond1)
              {
                OutPut("weak" << r << endl);
                if(Cond0 != DIRICHLET)
                  continue;
              }
            }
            OutPut("weak" << r << endl);
#ifdef __3D__
            cell->GetVertex(r)->GetCoords(x0, y0, z0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1, z1);
#else
            cell->GetVertex(r)->GetCoords(x0, y0);
            cell->GetVertex((r+1) % N_Edges)->GetCoords(x1, y1);
#endif
            if(out>2){  OutPut("Ecke A " << x0 << " " << y0 << " Ecke B " << x1 << " " << y1 << endl);}
            // compute length of the boundary edge
            hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            // compute normal vector to this boundary (normalized)
            nx = (y1-y0)/hE;
            ny = (x0-x1)/hE;
            // tangential normal vector to this boundary (normalized)
            tx = (x1-x0)/hE;
            ty = (y1-y0)/hE;
            sigma_par = TDatabase::ParamDB->WEAK_BC_SIGMA;
            if(out>2){OutPut("weak_bc: tau edge stabilization: " << tau_par << endl)};

            Entries1 = sqmatrices[n]->GetEntries();
            RowPtr1 = sqmatrices[n]->GetRowPtr();
            ColInd1 = sqmatrices[n]->GetKCol();
            fespace = sqmatrices[n]->GetFESpace();
            ActiveBound = fespace->GetActiveBound();
            if(out>2)
            {
              OutPut("ActiveBound von sqmatrix " << n << "  :" << ActiveBound  << endl);
              for(j=0;j<N_Points1D; j++)
              {
                OutPut("Det["<<j<<"]:" <<  AbsDetjk1D[r][j] <<" b1: "<< Coeffs[j][1] <<" b2: " << Coeffs[j][2] << endl);
              }
            }

            //edge integrals: test function from cell <-> ansatz function from cell
            if(out>2){  OutPut("testfct. cell<-> ansatzfct. cell Integral:" << endl);}
            for (ii=0;ii<N_;ii++)                 //ii - test function
            {
              // look for 'ii'-th row in all matrices
              dof_ii = DOF[ii];
              //OutPut("dof_ii" << dof_ii << endl);
              // Dirichlet node

              if (dof_ii>=ActiveBound)continue;

              // for all dof in the mesh cell
              // jj - ansatz function
              for (jj=0;jj<N_;jj++)
              {
                dof_jj = DOF[jj];
                // initialize the boundary integrals
                integral=0;
                for (j=0;j<N_Points1D;j++)        // compute edge integral; Assumption: N_Points1D cell =  N_Points1D neighbour !!!!
                {
                  integrant = 0;
                                                  // contribution of diffusive term
                  integrant += Coeffs[j][0] * (-xyval_ref1D[r][j][ii]*(xderiv_ref1D[r][j][jj]*nx + yderiv_ref1D[r][j][jj]*ny) - xyval_ref1D[r][j][jj]*(xderiv_ref1D[r][j][ii]*nx + yderiv_ref1D[r][j][ii]*ny));
                                                  // penalty term to achieve coercitivity
                  integrant += sigma_par*Coeffs[j][0]*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj]/hE;
                                                  // contribution of convective term
                  if((Coeffs[j][1]*nx + Coeffs[j][2]*ny) < 0) integrant -= (Coeffs[j][1]*nx + Coeffs[j][2]*ny)*xyval_ref1D[r][j][ii]*xyval_ref1D[r][j][jj];
                  w = weights1D[j]*hE/2;
                  integral += w*integrant;        // integral on the edge
                }

                // update matrix
                found = 0;
                for (ll=RowPtr1[dof_ii];ll < RowPtr1[dof_ii+1]; ll++)
                {
                  if (ColInd1[ll] == dof_jj)
                  {
                    if(out>2)
                    {
                      OutPut("integral: " << integral << " DOF Testfkt: " << dof_ii
                        <<  " DOF Ansatzfkt: " << dof_jj << endl);
                    }

                    Entries1[ll] += integral;
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
            if (out>2){  OutPut(endl);}
          }                                       //endif weak
        }                                         //endfor r (loop over edges)
      }                                           //endfor r (loop over edges)
    }                                             // endfor n (loop over sqmatrices)

    if(weak>=1)
    {
      for(n=0;n<n_rhs;n++)
      {
        fespace = ferhs[n];
        ActiveBound = fespace->GetActiveBound();
        CurrentElement = fespace->GetFE2D(i, cell);

        N_ = N_BaseFunct[CurrentElement];
        if(out>2){OutPut("N_BaseFunct " << N_ << endl);}

        local_rhs = righthand+n*MaxN_BaseFunctions2D;
        RHS = rhs[n];
        // find space for this linear form

        DirichletBound = fespace->GetHangingBound();

        // dof of the rhs nodes connected to this cell
        DOF = RhsGlobalNumbers[n] + RhsBeginIndex[n][i];
        if(out>2)
        {
          for(k=0; k<N_; k++)
          {
            OutPut("DOF[k]: k: " << k << " DOF: " << DOF[k]  << endl);
          }
        }

        BoundaryCondition = BoundaryConditions[n];
        BoundaryValue = BoundaryValues[n];

        N_Joints = cell->GetN_Edges();

        if(out>2)OutPut("N_Joints: " << N_Joints << endl);

        RefTrans = TFEDatabase2D::GetOrig(1, LocalUsedElements,
          Coll, cell, SecondDer,
          N_Points, xi, eta, weights, X, Y, AbsDetjk);
        if(N_Parameters>0)                        // get parameters of equ.
          Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);
                                                  // get coefficients of pde
        if(Coeff) Coeff(N_Points, X, Y, Param, Coeffs);

        for(m=0;m<N_Joints;m++)
        {
          joint = cell->GetJoint(m);
          if(joint->GetType() == BoundaryEdge||
           joint->GetType() == InterfaceJoint ||
            joint->GetType() == IsoBoundEdge)
          {
            if(joint->GetType() == BoundaryEdge||
            joint->GetType() == InterfaceJoint)
            {
              boundedge = (TBoundEdge *)joint;
              BoundComp = boundedge->GetBoundComp();
              boundedge->GetParameters(t0, t1);
            }
            else
            {
              isoboundedge = (TIsoBoundEdge *)joint;
              BoundComp = isoboundedge->GetBoundComp();
              isoboundedge->GetParameters(t0, t1);
            }
            // get id of the boundary component
            comp=BoundComp->GetID();
            // get type of the boundary condition at the beginning
            // and at the end of the current edge
            if (t0 < t1)
            {
              BoundaryCondition(comp, t0+eps, Cond0);
              BoundaryCondition(comp, t1-eps, Cond1);
            }
            else
            {
              BoundaryCondition(comp, t0-eps, Cond0);
              BoundaryCondition(comp, t1+eps, Cond1);
            }
            // only one boundary condition per edge allowed
            if(Cond0 == Cond1)
            {
              if(Cond0 == DIRICHLET)
              {
                // get polynomial degree of fe
                if(out>2){OutPut("Edge " << m << endl);}
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
                  (CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1D = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                if(out>2){OutPut("QuadFormula " << qf1D << endl);}
                qf1D->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
                  ->MakeRefElementData(LineQuadFormula);

                TFEDatabase2D::GetOrigFromRef(RefTrans,N_Points1D, xi1D[BaseFunctCell][m],
                  eta1D[BaseFunctCell][m],
                  X1D[m], Y1D[m], AbsDetjk1D[m]);

                for(j=0;j<N_Points1D;j++)         // get values and derivatives in original cell
                {
                  TFEDatabase2D::GetOrigValues(RefTrans, xi1D[BaseFunctCell][m][j],
                    eta1D[BaseFunctCell][m][j],
                    TFEDatabase2D::GetBaseFunct2D(BaseFunctCell),
                    Coll, (TGridCell *)cell,
                    xietaval_ref1D[BaseFunctCell][m][j],
                    xideriv_ref1D[BaseFunctCell][m][j],
                    etaderiv_ref1D[BaseFunctCell][m][j],
                    xyval_ref1D[m][j],
                    xderiv_ref1D[m][j],
                    yderiv_ref1D[m][j]);
                }

                // get vertices of boundary edge
#ifdef __3D__
                cell->GetVertex(m)->GetCoords(x0, y0, z0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
#else
                cell->GetVertex(m)->GetCoords(x0, y0);
                cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
#endif
                // compute the length of the boundary edge
                hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
                nx = (y1-y0)/hE;
                ny = (x0-x1)/hE;
                // tangential normal vector to this boundary (normalized)
                tx = (x1-x0)/hE;
                ty = (y1-y0)/hE;
                sigma_par = TDatabase::ParamDB->WEAK_BC_SIGMA;

                // compute boundary integral: ansatz function from cell with Dirichlet boundary value
                for(k=0;k<N_;k++)
                {
                  l3 = DOF[k];
                  integral=0;
                  for(l=0;l<N_LinePoints;l++)
                  {
                    if(out>2){OutPut("l:" << l << " normal vector: " << nx <<  " " << ny << endl);};
                    // get quadrature point on the boundary
                    t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
                    if(out>2){OutPut("Quad formula data: Quadrature point[l]:" << zeta[l] <<  " weight: " << LineWeights[l] << endl);};
                    // get value in this quadrature point (in s)
                    BoundaryValue(comp, t, s);
                    if(out>2)
                    {
                      OutPut("Parameter Quadrature Point: " << t << " Boundary Value: " << s << endl);
                      OutPut("Function values Quadrature Point: " << k << " xyval " << xyval_ref1D[m][l][k] << " xderiv " << xderiv_ref1D[m][l][k] << " yderiv " << yderiv_ref1D[m][l][k] << endl);

                    };
                    // multiply value with weights from quadrature formula
                    // and determinant from integral transformation to the
                    // unit edge (-1,1)
                    integrant = 0;
                                                  // contribution of diffusive term
                    integrant += Coeffs[j][0]*(-s*(xderiv_ref1D[m][l][k]*nx + yderiv_ref1D[m][l][k]*ny));
                                                  // penalty term to achieve coercitivity
                    integrant += sigma_par*Coeffs[j][0]*s*xyval_ref1D[m][l][k]/hE;
                    if((Coeffs[l][1]*nx + Coeffs[l][2]*ny) < 0) integrant -= (Coeffs[l][1]*nx + Coeffs[l][2]*ny)*s*xyval_ref1D[m][l][k];
                    integral += hE/2 * LineWeights[l]*integrant;
                  }
                  if(l3 < ActiveBound) RHS[l3] += integral;
                  else
                  {
                    OutPut("Index l3 bigger than ActiveBound!");
                    exit(4711);
                  }
                  if(out>2) OutPut("DOF: " << l3 << " integral: " << integral << " RHS: " << RHS[l3] << endl);
                }
              }                                   // endif
            }                                     // endif (Cond0==Cond1)
            else
            {
              OutPut("different boundary condition on one edge ");
              OutPut("are not allowed!" << endl);
              exit(4711);
            }
          }                                       // endif (boundary joint)
        }                                         // endfor m (N_Joints)
      }                                           // endfor n (rhs)
    }                                             //endif WEAK
  }                                               // endfor i (loop over cells)

  if (out==2)
      OutPut("free memory " << endl);
  if(n_sqmatrices)
  {
    delete GlobalNumbers;
    delete BeginIndex;
  }

  if(n_matrices)
  {
    delete AnsatzGlobalNumbers;
    delete AnsatzBeginIndex;
    delete TestGlobalNumbers;
    delete TestBeginIndex;
  }

  if(N_AllMatrices)
  {
    delete LocMatrices;
    delete Matrices[0];
    delete Matrices;
  }

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
  }
  for(l=0;l<ref_n;l++)
  {
    delete value_basefunct_ori[l];
    delete xderiv_basefunct_ori[l];
    delete yderiv_basefunct_ori[l];
    delete value_basefunct_ori_neigh[l];
    delete xderiv_basefunct_ori_neigh[l];
    delete yderiv_basefunct_ori_neigh[l];
  }

  for (i=0;i<N_BaseFuncts2D;i++)
  {
      for (j=0;j<ref_n;j++)
      {
	  delete value_basefunct_ref1D[i][j];
	  delete xderiv_basefunct_ref1D[i][j];
	  delete yderiv_basefunct_ref1D[i][j];
      }
      delete value_basefunct_ref1D[i];
      delete xderiv_basefunct_ref1D[i];
      delete yderiv_basefunct_ref1D[i];
  }

  for (i=0;i<N_Points1D;i++)
  {
    delete xyval_refNeigh1D[i];
    delete xderiv_refNeigh1D[i];
    delete yderiv_refNeigh1D[i];
  }

  delete aux;
  delete aux2;
  delete aux4;
  delete aux3;

  if(n_rhs)
  {
    delete righthand;
    delete LocRhs;
    delete RhsBeginIndex;
    delete RhsGlobalNumbers;
  }

  for(i=0; i < N_BaseFuncts2D; i++)
  {
    for(j=0; j < 4; j++)
      {
	for(m=0; m < MaxN_QuadPoints_1D; m++)
          {
	    delete xietaval_ref1D[i][j][m];
	    delete xideriv_ref1D[i][j][m];
	    delete etaderiv_ref1D[i][j][m];
	  }
	delete xietaval_ref1D[i][j];
	delete xideriv_ref1D[i][j];
	delete etaderiv_ref1D[i][j];
      }
    delete xietaval_ref1D[i];
    delete xideriv_ref1D[i];
    delete  etaderiv_ref1D[i];
  }

  delete value_basefunct_ref1D;
  delete xderiv_basefunct_ref1D;
  delete yderiv_basefunct_ref1D;

  delete xietaval_ref1D;
  delete xideriv_ref1D;
  delete etaderiv_ref1D;

  if (weak>=1) // ???
  {
    TDatabase::ParamDB->INTERNAL_DO_NOT_RESPECT_DIRICHLET_BC = 1;
    TDatabase::ParamDB->WEAK_BC = 2;
  }
  
  int N_Rows;
  // ####################################################################
  // print the whole matrix -- SECOND
  // ####################################################################
  if(out>2)
  {
    for(k=0;k<n_sqmatrices;k++)
    {
      cout << endl;
      cout << "sqmatrix: " << k << endl;
      RowPtr = sqmatrices[k]->GetRowPtr();
      Entries = sqmatrices[k]->GetEntries();
      ColInd = sqmatrices[k]->GetKCol();
      N_Rows = sqmatrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=RowPtr[i+1];
        for(j=RowPtr[i];j<end;j++)
        {
          // cout << j << endl;
          cout << "Matrix: " << setw(5) << i << setw(5) << ColInd[j] << "   ";
          cout << setw(10) << Entries[j] << endl;
        }
      }
      cout << endl;
    }                                             // endfor k
  }                                               //endif

  /*  for(k=0;k<n_matrices;k++)
    {
      cout << endl;
      cout << "matrix: " << k << endl;
      RowPtr = matrices[k]->GetRowPtr();
      Entries = matrices[k]->GetEntries();
      ColInd = matrices[k]->GetKCol();
      N_Rows = matrices[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
        end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  cout << setw(5) << i << setw(5) << ColInd[j] << "   ";
  cout << setw(10) << Entries[j] << endl;
  }
  }
  cout << endl;
  } // endfor k

  for(k=0;k<n_rhs;k++)
  {
  cout << "rhs: " << k << endl;
  N_Rows = ferhs[k]->GetN_DegreesOfFreedom();
  RHS=rhs[k];
  for(i=0;i<N_Rows;i++)
  cout << setw(5) << i << setw(20) << RHS[i] << endl;
  }
  */

}                                                 // end of Assemble

 
#endif



// void Assemble2D(int n_fespaces, TFESpace2D **fespaces,
//                 int n_sqmatrices, TSquareMatrix2D **sqmatrices,
//                 int n_matrices, TMatrix2D **matrices,
//                 int n_rhs, double **rhs, TFESpace2D **ferhs,
//                 BoundCondFunct2D **BoundaryConditions,
//                 BoundValueFunct2D **BoundaryValues,
//                 LocalAssembling2D& la
// #ifdef __3D__
//                 , TAux2D3D *Aux2D3D
// #endif
//                 , int AssemblePhaseID
// )
// {
//   int N_AllMatrices = n_sqmatrices+n_matrices;
//   int i,j,k,l,l1,l2,l3,n,m, N_LocalUsedElements,ii,jj,ll,ij;
//   int N_Cells, N_Points, N_Parameters, N_, N_Hanging;
//   int N_Test, N_Ansatz, N_Joints, N_Edges;
//   int Used[N_FEs2D];
//   int *N_BaseFunct;
//   BaseFunct2D *BaseFuncts;
//   TFESpace2D *fespace;
//   FE2D LocalUsedElements[N_FEs2D], CurrentElement;
//   FE2D TestElement, AnsatzElement;
//   QuadFormula1D LineQuadFormula;
//   TQuadFormula1D *qf1;
//   TCollection *Coll;
//   TBaseCell *cell;
//   TJoint *joint;
//   TBoundEdge *boundedge;
//   TIsoBoundEdge *isoboundedge;
//   int **GlobalNumbers, **BeginIndex;
//   int **RhsGlobalNumbers, **RhsBeginIndex;
//   int **TestGlobalNumbers, **TestBeginIndex;
//   int **AnsatzGlobalNumbers, **AnsatzBeginIndex;
//   TFE2D *ele;
//   TFEDesc2D *FEDesc_Obj;
//   double *weights, *xi, *eta;
//   double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
//   double AbsDetjk[MaxN_QuadPoints_2D];
//   double *Param[MaxN_QuadPoints_2D];
//   double *local_rhs;
//   double *righthand;
//   double **Matrices, *aux;
//   double **Matrix;
//   double ***LocMatrices, **LocRhs;
//   int LocN_BF[N_BaseFuncts2D];
//   BaseFunct2D LocBF[N_BaseFuncts2D];
//   double *AuxArray[MaxN_QuadPoints_2D];
//   int *DOF, ActiveBound, DirichletBound, end, last;
//   int *TestDOF, *AnsatzDOF;
//   double *Entries,*Entries1,*Entries2,*Entries3;
//   int *ColInd, *RowPtr;
//   int *ColInd1, *RowPtr1,*ColInd2, *RowPtr2, *ColInd3, *RowPtr3;
//   double *RHS, *MatrixRow;
//   double **HangingEntries, **HangingRhs;
//   double *CurrentHangingEntries, *CurrentHangingRhs;
//   int *HangingRowPtr, *HangingColInd;
//   THangingNode *hn, **HangingNodes;
//   HNDesc HNDescr;
//   THNDesc *HNDescr_Obj;
//   double *Coupling, v;
//   TBoundComp *BoundComp;
//   double t0, t1, t, s,integral[2];
//   int comp, dof_ii,dof_jj, found;
//   BoundCond Cond0, Cond1;
//   BoundCondFunct2D *BoundaryCondition;
//   BoundValueFunct2D *BoundaryValue;
//   TNodalFunctional2D *nf;
//   int N_EdgePoints;
//   double *EdgePoints;
//   double PointValues[MaxN_PointsForNodal2D];
//   double FunctionalValues[MaxN_BaseFunctions2D];
//   int *EdgeDOF, N_EdgeDOF;
//   int N_LinePoints;
//   double *LineWeights, *zeta;
//   double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-4;
// #ifdef __3D__
//   double z0, z1;
// #endif
//   double penetration_penalty, friction_parameter;
//   double **JointValues, *JointValue;
//   bool *SecondDer;
//   int *Bounds, NeumannBound, RobinBound;
//   int lr, mr;
//   double RobinScale = 16.0;
//   TVertex *ver0;
//   TRefTrans2D *rt;
// 
//   int N_Rows;
//   int AnsatzActiveBound, AnsatzHangingBound;
// 
//   // ########################################################################
//   // store information in local arrays
//   // ########################################################################
//   BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
//   N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
// 
//   if(n_sqmatrices+n_matrices)
//   {
//     HangingEntries = new double* [n_sqmatrices+n_matrices];
// 
//     for(i=0;i<n_sqmatrices;i++)
//     {
//       j = sqmatrices[i]->GetHangingN_Entries();
//       HangingEntries[i] = new double [j];
//       memset(HangingEntries[i], 0, SizeOfDouble*j);
//     }
// 
//     for(i=0;i<n_matrices;i++)
//     {
//       j = matrices[i]->GetHangingN_Entries();
//       HangingEntries[i+n_sqmatrices] = new double [j];
//       memset(HangingEntries[i+n_sqmatrices], 0, SizeOfDouble*j);
//     }
//   }
// 
//   if(n_sqmatrices)
//   {
//     GlobalNumbers = new int* [n_sqmatrices];
//     BeginIndex = new int* [n_sqmatrices];
//     for(i=0;i<n_sqmatrices;i++)
//     {
//       fespace = sqmatrices[i]->GetFESpace();
//       GlobalNumbers[i] = fespace->GetGlobalNumbers();
//       BeginIndex[i] = fespace->GetBeginIndex();
//     }                                             // endfor
//   }                                               // endif n_sqmatrices
// 
//   if(n_matrices)
//   {
//     TestGlobalNumbers = new int* [n_matrices];
//     AnsatzGlobalNumbers = new int* [n_matrices];
//     TestBeginIndex = new int* [n_matrices];
//     AnsatzBeginIndex = new int* [n_matrices];
//     for(i=0;i<n_matrices;i++)
//     {
//       fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetTestSpace();
//       TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
//       TestBeginIndex[i] = fespace->GetBeginIndex();
// 
//       fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetAnsatzSpace();
//       AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
//       AnsatzBeginIndex[i] = fespace->GetBeginIndex();
//     }                                             // endfor
//   }                                               // endif n_matrices
// 
//   if(n_rhs)
//   {
//     HangingRhs = new double* [n_rhs];
//     RhsBeginIndex = new int* [n_rhs];
//     RhsGlobalNumbers = new int* [n_rhs];
//     for(i=0;i<n_rhs;i++)
//     {
//       fespace = ferhs[i];
//       RhsBeginIndex[i] = fespace->GetBeginIndex();
//       RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
// 
//       j = fespace->GetN_Hanging();
//       HangingRhs[i] = new double [j];
//       memset(HangingRhs[i], 0, SizeOfDouble*j);
//     }                                             // endfor
// 
//     LocRhs = new double* [n_rhs];
//     righthand = new double [n_rhs*MaxN_BaseFunctions2D];
//     for(i=0;i<n_rhs;i++)
//       LocRhs[i] = righthand+i*MaxN_BaseFunctions2D;
//   }                                               // endif n_rhs
// 
//   //N_Parameters = Parameters->GetN_Parameters();
//   N_Parameters = la.GetN_Parameters();
//   
// #ifdef __3D__
//   N_Parameters += 7;                              // (u, ux, uy, uz, nx, ny, nz)
// #endif
// 
//   if(N_Parameters)
//   {
//     aux = new double [MaxN_QuadPoints_2D*N_Parameters];
//     for(j=0;j<MaxN_QuadPoints_2D;j++)
//       Param[j] = aux + j*N_Parameters;
//   }
// 
//   // 40 <= number of terms in bilinear form
//   // DUE NOTE CHANGE BELOW 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
//   aux = new double [MaxN_QuadPoints_2D*40];
//   for(j=0;j<MaxN_QuadPoints_2D;j++)
//     AuxArray[j] = aux + j*40;
// 
//   if(N_AllMatrices)
//   {
//     aux = new double
//       [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
//     Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
//     for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
//       Matrices[j] = aux+j*MaxN_BaseFunctions2D;
// 
//     LocMatrices = new double** [N_AllMatrices];
//     for(i=0;i<N_AllMatrices;i++)
//       LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
//   }                                               // endif N_AllMatrices
// 
//   //SecondDer = DiscreteForm->GetNeeds2ndDerivatives();
//   SecondDer = la.GetNeeds2ndDerivatives();
//   
//   // ########################################################################
//   // loop over all cells
//   // ########################################################################
//   Coll = fespaces[0]->GetCollection();            // all spaces use same Coll
//   N_Cells = Coll->GetN_Cells();
//   for(i=0;i<N_Cells;i++) // set the cell indices
//   {
//     cell = Coll->GetCell(i);
//     cell->SetCellIndex(i);
//   }
// 
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     TDatabase::ParamDB->INTERNAL_CELL = i;
// 
//     // only for multiphase flows
//     if(AssemblePhaseID >= 0)
//     {
//       if(AssemblePhaseID != cell->GetRegionID() )
//         continue;
//     }
// 
//     // ####################################################################
//     // find local used elements on this cell
//     // ####################################################################
// 
//     for(j=0;j<n_fespaces;j++)
//     {
//       CurrentElement = fespaces[j]->GetFE2D(i, cell);
//       LocalUsedElements[j] = CurrentElement;
//       LocN_BF[j] = N_BaseFunct[CurrentElement];
//       LocBF[j] = BaseFuncts[CurrentElement];
//     }
// 
//     N_LocalUsedElements = n_fespaces;
// 
//     // ####################################################################
//     // calculate values on original element
//     // ####################################################################
//  
//     TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
//                            Coll, cell, SecondDer,
//                            N_Points, xi, eta, weights, X, Y, AbsDetjk);
//     
//     //Parameters->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param);
//     la.GetParameters(N_Points, Coll, cell, i, X, Y, Param);
// 
//     if((TDatabase::ParamDB->DISCTYPE == SDFEM)
//       || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
//       || (TDatabase::ParamDB->CELL_MEASURE == 4))
//     {
//       TDatabase::ParamDB->INTERNAL_LOCAL_DOF = i;
//       N_Edges = cell->GetN_Edges();
//       for (ij=0;ij<N_Edges;ij++)
//       {
//         TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
//         TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
//       }
//       if (N_Edges==3)
//         TDatabase::ParamDB->INTERNAL_VERTEX_X[3] = -4711;
//       TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
//     }
// 
// #ifdef __3D__
//     if(Aux2D3D)
//       Aux2D3D->GetGradient(i, N_Points, xi, eta, Param);
// #endif
//     // use DiscreteForm to assemble a few matrices and
//     // right-hand sides at once
// 
//     /*
//     if(DiscreteForm)
//     {
//       DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk,
//         hK, X, Y,
//         LocN_BF, LocBF,
//         Param, AuxArray,
//         cell,
//         N_AllMatrices, n_rhs,
//         LocMatrices, LocRhs);
//     }
//     */
//     la.GetLocalForms(N_Points, weights, AbsDetjk, X, Y, LocN_BF, LocBF,
//                      Param, AuxArray, cell, N_AllMatrices, n_rhs, LocMatrices,
//                      LocRhs);
// 
//     N_Joints = cell->GetN_Joints();
//     // ####################################################################
//     // add local matrices to global matrices (ansatz == test)
//     // ####################################################################
//     for(j=0;j<n_sqmatrices;j++)
//     {
//       // find space for this bilinear form
//       fespace = sqmatrices[j]->GetFESpace();
//       CurrentElement = fespace->GetFE2D(i, cell);
//       N_ = N_BaseFunct[CurrentElement];
// 
//       Matrix = LocMatrices[j];
//       Entries = sqmatrices[j]->GetEntries();
//       RowPtr = sqmatrices[j]->GetRowPtr();
//       ColInd = sqmatrices[j]->GetKCol();
// 
//       CurrentHangingEntries = HangingEntries[j];
//       HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
//       HangingColInd = sqmatrices[j]->GetHangingKCol();
// 
//       ActiveBound = fespace->GetActiveBound();
//       DirichletBound = fespace->GetHangingBound();
//       DOF = GlobalNumbers[j] + BeginIndex[j][i];
//       Bounds = fespace->GetBoundaryNodesBound();
//       NeumannBound = Bounds[0];
//       RobinBound = Bounds[1];
// 
//       /*
//       BoundaryCondition = BoundaryConditions[j];
//       for(m=0;m<N_Joints;m++)
//       {
//       joint = cell->GetJoint(m);
//         if(joint->GetType() == BoundaryEdge ||
//            joint->GetType() == IsoBoundEdge)
//         {
//           if(joint->GetType() == BoundaryEdge)
//           {
//             boundedge = (TBoundEdge *)joint;
//       BoundComp = boundedge->GetBoundComp();
//       boundedge->GetParameters(t0, t1);
//       }
//       else
//       {
//       isoboundedge = (TIsoBoundEdge *)joint;
//       BoundComp = isoboundedge->GetBoundComp();
//       isoboundedge->GetParameters(t0, t1);
//       }
//       // get id of the boundary component
//       comp = BoundComp->GetID();
//       // get type of the boundary condition at the beginning
//       // and at the end of the current edge
//       BoundaryCondition(comp, t0, Cond0);
//       //      cout << "bound1" << endl;
//       if(Cond0 == ROBIN)
//       {
//       #ifdef __2D__
//       //cout << "robin" << endl;
//       // Robin
//       lr = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
// 
//       // get a suitable line quadrature formula
//       LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*lr);
//       qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//       qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
// 
//       TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
//       ->MakeRefElementData(LineQuadFormula);
// 
//       JointValues=TFEDatabase2D::GetJointValues2D(
//       BaseFuncts[CurrentElement], LineQuadFormula, m);
//       // get vertices of boundary edge
//       cell->GetVertex(m)->GetCoords(x0, y0);
//       cell->GetVertex((m+1) % 4)->GetCoords(x1, y1);
//       // compute (half of the) length of the boundary edge
//       hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
//       // cout << "x0: " << x0 << " y0: " << y0 << endl;
//       // cout << "x1: " << x1 << " y1: " << y1 << endl;
//       // compute boundary integral
//       for(n=0;n<N_LinePoints;n++)
//       {
//       // values of test functions in this quadrature point
//       JointValue = JointValues[n];
//       // cout << "Zeta  :" << zeta[n] << endl;
// 
//       // get quadrature point on the boundary
//       for(l=0;l<N_;l++)
//       {
//       MatrixRow = Matrix[l];
//       s = JointValue[l];
//       // multiply value with weights from quadrature formula
//       // and determinant from integral transformation to the
//       // unit edge (-1,1)
//       s *= hE * LineWeights[n];
//       // !! hold the robin boundary values of
//       // !! the function alpha from parameters
//       // s *= alpha;
//       // update rhs for all test functions
//       for(k=0;k<N_;k++)
//       MatrixRow[k] += s*JointValue[k]*RobinScale;
//       } // endfor l
//       } // endfor n
//       #endif
//       } // end Robin
//       } // endif BoundEdge
//       } // endfor m
//       */
// 
//       // add local matrix to global
//       for(m=0;m<N_;m++)
//       {
//         l=DOF[m];
//         MatrixRow = Matrix[m];
//         if(l<ActiveBound)
//         {
//           for(k=0;k<N_;k++)
//           {
//             // DOF[k] is the global index of the k-th local degree of freedom
//             // MatrixRow[k] is the assembled value corresponding to the m-th
//             // local test function and k-th local ansatz function. That means it
//             // corresponds to the l=DOF[m]-th global test function and the 
//             // DOF[k]-th global ansatz function
//              sqmatrices[j]->add(l, DOF[k], MatrixRow[k]);
//           }
//         }                                         // endif l
//         else
//         {
//           if(l<DirichletBound)
//           {
//             // hanging node
//             l -= ActiveBound;
//             end = HangingRowPtr[l+1];
//             for(n=HangingRowPtr[l];n<end;n++)
//             {
//               for(k=0;k<N_;k++)
//               {
//                 if(DOF[k] == HangingColInd[n])
//                 {
//                   CurrentHangingEntries[n] += MatrixRow[k];
//                   break;
//                 }                                 // endif
//               }                                   // endfor k
//             }                                     // endfor n
//           }
//           else
//           {
//             // Dirichlet node
//             sqmatrices[j]->set(l,l,1.0);
//           }
//         }
//       }                                           // endfor m
//     }                                             // endfor j
// 
//     // ####################################################################
//     // add local matrices to global matrices (ansatz != test)
//     // ####################################################################
//     for(j=0;j<n_matrices;j++)
//     {
//       TestElement = ((TFESpace2D *) matrices[j]->GetStructure()->
//         GetTestSpace())->GetFE2D(i, cell);
//       AnsatzElement = ((TFESpace2D *) matrices[j]->GetStructure()->
//         GetAnsatzSpace())->GetFE2D(i, cell);
//       // cout << "non square matrix: " << j << endl;
//       // cout << "TestElement: " << TestElement << endl;
//       // cout << "AnsatzElement: " << AnsatzElement << endl;
// 
//       N_Test = N_BaseFunct[TestElement];
//       N_Ansatz = N_BaseFunct[AnsatzElement];
// 
//       Matrix = Matrices+(j+n_sqmatrices)*MaxN_BaseFunctions2D;
// 
//       Entries = matrices[j]->GetEntries();
//       RowPtr = matrices[j]->GetRowPtr();
//       ColInd = matrices[j]->GetKCol();
// 
//       TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
//       AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];
// 
//       fespace = (TFESpace2D *)(matrices[j]->GetStructure()->GetTestSpace());
//       ActiveBound = fespace->GetActiveBound();
//       DirichletBound = fespace->GetHangingBound();
// 
//       CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
//       HangingRowPtr = matrices[j]->GetHangingRowPtr();
//       HangingColInd = matrices[j]->GetHangingKCol();
// 
//       // add local matrix to global
//       for(m=0;m<N_Test;m++)
//       {
//         l=TestDOF[m];
//         MatrixRow = Matrix[m];
//         // cout << "DOF: " << l << endl;
//         if(l<ActiveBound || l>=DirichletBound)
//         {
//           end=RowPtr[l+1];
//           for(n=RowPtr[l];n<end;n++)
//           {
//             for(k=0;k<N_Ansatz;k++)
//             {
//               if(AnsatzDOF[k] == ColInd[n])
//               {
//                 // cout << m << "   " << k << endl << n << endl;
//                 Entries[n] += MatrixRow[k];
//                 break;
//               }                                   // endif
//             }                                     // endfor k
//           }                                       // endfor n
//         }
//         else
//         {
//           // hanging node
//           l -= ActiveBound;
//           end = HangingRowPtr[l+1];
//           for(n=HangingRowPtr[l];n<end;n++)
//           {
//             // cout << l << " HangingColInd: " << HangingColInd[n] << endl;
//             for(k=0;k<N_Ansatz;k++)
//             {
//               // cout << "AnsatzDOF: " << AnsatzDOF[k] << endl;
//               if(AnsatzDOF[k] == HangingColInd[n])
//               {
//                 CurrentHangingEntries[n] += MatrixRow[k];
//                 break;
//               }                                   // endif
//             }                                     // endfor k
//           }                                       // endfor n
//         }
//       }                                           // endfor m
//     }                                             // endfor j  (n_matrices)
// 
//     // ####################################################################
//     // add local right-hand sides to global right-hand side
//     // ####################################################################
//     for(j=0;j<n_rhs;j++)
//     {
//       fespace = ferhs[j];
//       ActiveBound = fespace->GetActiveBound();
//       CurrentElement = fespace->GetFE2D(i, cell);
// 
//       N_ = N_BaseFunct[CurrentElement];
// 
//       local_rhs = righthand+j*MaxN_BaseFunctions2D;
//       RHS = rhs[j];
//       CurrentHangingRhs = HangingRhs[j];
//       // find space for this linear form
// 
//       ActiveBound = fespace->GetActiveBound();
//       DirichletBound = fespace->GetHangingBound();
// 
//       // dof of the rhs nodes connected to this cell
//       DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
// 
//       // add local right-hand side to the global one
//       for(m=0;m<N_;m++)
//       {
//         if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
//         {
//           l = 0;
//           N_Edges=cell->GetN_Edges();
//           for(jj=0;jj<N_Edges;jj++)                        // loop over all edges of cell
//           {
//             ver0 = cell->GetVertex(jj);
//             t0 =  ver0->GetClipBoard();
//             // vertex not on the boundary
//             if (t0<-1e-8)
//               continue;
//             // component of boundary
//             comp = floor(t0+1e-8);
//             // parameter
//             t0 -= comp;
//             // get boundary condition
//             BoundaryConditions[0](comp, t0, Cond0);
//             // Dirichlet
//             if (Cond0== DIRICHLET)
//             {
//               l = -4711;
//             }
//           }
//           if (l==-4711) 
//           break; // do nothing for this mesh cell
//         }
//     
//         l=DOF[m];
//         //cout << "DOF: " << l << endl;
//         if(l<ActiveBound)
//         {
//           // node l is inner or Neumann node
//           RHS[l] += local_rhs[m];
//           // cout << l << " " << RHS[l] << " " << local_rhs[m]<< " "<<endl;;
//         }                                         // endif l
//         else
//         {
//           if(l<DirichletBound)
//           {
//             // hanging node
//             l -= ActiveBound;
//             CurrentHangingRhs[l] += local_rhs[m];
//           }
//         }
//       }                                           // endfor m
// 
//       BoundaryCondition = BoundaryConditions[j];
//       BoundaryValue = BoundaryValues[j];
//       ele = TFEDatabase2D::GetFE2D(CurrentElement);
//       //if ((ele >= D_P1_2D_Q_A)&&(ele<= D_P3_2D_Q_M))
//       //  continue;
// 
//       nf = ele->GetNodalFunctional2D();
// 
//       if(TDatabase::ParamDB->SUPERCONVERGENCE_ORDER)
//       {
//         /* Superconvergence boundary interpolation */
//         if(nf->GetID() == NF_C_Q_Q2_2D)
//           nf = TFEDatabase2D::GetNodalFunctional2D(NF_S_Q_Q2_2D);
//       }
// 
//       nf->GetPointsForEdge(N_EdgePoints, EdgePoints);
// 
//       FEDesc_Obj = ele->GetFEDesc2D();
//       N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
//       // setting Dirichlet boundary condition
//       N_Joints = cell->GetN_Edges();
// 
//       for(m=0;m<N_Joints;m++)
//       {
//         joint = cell->GetJoint(m);
// 
//         if(joint->GetType() == BoundaryEdge ||
//            joint->GetType() == IsoBoundEdge ||
//            joint->GetType() == InterfaceJoint)
//         {
//           if(joint->GetType() == BoundaryEdge||
//            joint->GetType() == InterfaceJoint)
//           {
//             boundedge = (TBoundEdge *)joint;
//             BoundComp = boundedge->GetBoundComp();
//             boundedge->GetParameters(t0, t1);
//           }
//           else
//           {
//             isoboundedge = (TIsoBoundEdge *)joint;
//             BoundComp = isoboundedge->GetBoundComp();
//             isoboundedge->GetParameters(t0, t1);
//           }
//           // get id of the boundary component
//           comp=BoundComp->GetID();
//           // get type of the boundary condition at the beginning
//           // and at the end of the current edge
//           if (t0 < t1)
//           {
//             BoundaryCondition(comp, t0+eps, Cond0);
//             BoundaryCondition(comp, t1-eps, Cond1);
//           }
//           else
//           {
//             BoundaryCondition(comp, t0-eps, Cond0);
//             BoundaryCondition(comp, t1+eps, Cond1);
//           }
// 
//           // only one boundary condition per edge allowed
//           if(Cond0 == Cond1)
//           {
//             switch(Cond0)
//             {
//               case DIRICHLET:
//                 // if DG
//                 if (N_EdgeDOF==0)
//                   break;
//                 // read boundary values for each quadrature point
//                 for(l=0;l<N_EdgePoints;l++)
//                 {
//                   s = EdgePoints[l];
//                   t = 0.5*(t0*(1-s) + t1*(1+s));
//                   BoundaryValue(comp, t, PointValues[l]);
//                 }                                 // endfor l
//                 // compute boundary values for each dof on the
//                 // boundary edge with the nodal functionals
// 
//                 nf->GetEdgeFunctionals(Coll, cell, m, PointValues,
//                   FunctionalValues);
//                 EdgeDOF = FEDesc_Obj->GetJointDOF(m);
//                 // save boundary values of each dof on the boundary
//                 // edge in the rhs
//                 for(l=0;l<N_EdgeDOF;l++)
//                 {
//                   RHS[DOF[EdgeDOF[l]]] = FunctionalValues[l];
//                 }
//                 break;
// 
//               case NEUMANN:
//                 // get polynomial degree of fe
//                 l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
//                   (CurrentElement);
//                 // get a suitable line quadrature formula
//                 LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
//                 qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//                 qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
//                 TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
//                   ->MakeRefElementData(LineQuadFormula);
//                 JointValues=TFEDatabase2D::GetJointValues2D(
//                   BaseFuncts[CurrentElement], LineQuadFormula, m);
//                 TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
//                   ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
//                 // get vertices of boundary edge
// #ifdef __3D__
//                 cell->GetVertex(m)->GetCoords(x0, y0, z0);
//                 cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1, z1);
// #else
//                 cell->GetVertex(m)->GetCoords(x0, y0);
//                 cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
// #endif
//                 // compute (half of the) length of the boundary edge
//                 hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
//                 // compute boundary integral
//                 for(l=0;l<N_LinePoints;l++)
//                 {
//                   // values of test functions in this quadrature point
//                   JointValue = JointValues[l];
//                   // get quadrature point on the boundary
//                   t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
//                   // get value in this quadrature point (in s)
//                   BoundaryValue(comp, t, s);
//                   // multiply value with weights from quadrature formula
//                   // and determinant from integral transformation to the
//                   // unit edge (-1,1)
//                   s *= hE * LineWeights[l];
//                   // update rhs for all test functions
//                   for(k=0;k<N_;k++)
//                     if((l3 = DOF[k])<ActiveBound)
//                       RHS[l3] += s*JointValue[k];
//                 }
//                 TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])
//                   ->ChangeBF(Coll, cell, N_LinePoints, JointValues);
//                 break;
//               case ROBIN:
// #ifdef __2D__
//                 // get polynomial degree of fe
//                 l = TFEDatabase2D::GetPolynomialDegreeFromFE2D
//                   (CurrentElement);
//                 // get a suitable line quadrature formula
//                 LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
//                 qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
//                 qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
//                 TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)
//                   ->MakeRefElementData(LineQuadFormula);
//                 JointValues=TFEDatabase2D::GetJointValues2D(
//                   BaseFuncts[CurrentElement], LineQuadFormula, m);
//                 // get vertices of boundary edge
//                 cell->GetVertex(m)->GetCoords(x0, y0);
//                 cell->GetVertex((m+1) % N_Joints)->GetCoords(x1, y1);
//                 // compute (half of the) length of the boundary edge
//                 hE = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))/2;
//                 // compute boundary integral
//                 for(l=0;l<N_LinePoints;l++)
//                 {
//                   // values of test functions in this quadrature point
//                   JointValue = JointValues[l];
//                   // get quadrature point on the boundary
//                   t = t0 + 0.5*(t1-t0)*(zeta[l]+1);
//                   // get value in this quadrature point (in s)
//                   BoundaryValue(comp, t, s);
//                   // multiply value with weights from quadrature formula
//                   // and determinant from integral transformation to the
//                   // unit edge (-1,1)
//                   s *= hE * LineWeights[l];
//                   // update rhs for all test functions
//                   for(k=0;k<N_;k++)
//                     if((l3 = DOF[k])<ActiveBound)
//                       RHS[l3] += s*JointValue[k];
//                 }
// #endif
//                 break;
//               case SLIP:
//                 OutPut("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !"<< endl);
//                 exit(4711);
//                 break;
//               case SLIP_FRICTION_PENETRATION_RESISTANCE:
//                 // do nothing here
//                 // everything is done in Assemble2DSlipBC, see below
//                 break;
// 
//             case FREESURF:
//             case INTERFACE:
//               // do nothing here
//               // everything is done in Freesurfint, see Freesurface2D.C
//               break;
//               default :
//                 OutPut("Unknown boundary condition !"<< endl);
//                 exit(4711);
// 
//             }                                     // endswitch Cond0
//           }                                       // endif (Cond0==Cond1)
//           else
//           {
//             OutPut("different boundary condition on one edge ");
//             OutPut("are not allowed!" << endl);
//             exit(4711);
//           }
//         }                                         // endif (boundary joint)
//       }                                           // endfor m (N_Joints)
//     }                                             // endfor j (n_rhs)
//   }                                               // endfor i (N_Cells)
// 
//   // ####################################################################
//   // modify matrix according to coupling
//   // ####################################################################
//   for(j=0;j<n_sqmatrices;j++)
//   {
//     fespace = sqmatrices[j]->GetFESpace();
//     N_ = fespace->GetN_Hanging();
//     HangingNodes = fespace->GetHangingNodes();
// 
//     Entries = sqmatrices[j]->GetEntries();
//     RowPtr = sqmatrices[j]->GetRowPtr();
//     ColInd = sqmatrices[j]->GetKCol();
// 
//     CurrentHangingEntries = HangingEntries[j];
//     HangingRowPtr = sqmatrices[j]->GetHangingRowPtr();
//     HangingColInd = sqmatrices[j]->GetHangingKCol();
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     for(i=0;i<N_;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
//       k = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       end = HangingRowPtr[i+1];
//       for(n=HangingRowPtr[i];n<end;n++)
//       {
//         v = CurrentHangingEntries[n];
//         m = HangingColInd[n];
//         for(l=0;l<k;l++)
//         {
//           l1 = DOF[l];
//           if(l1<ActiveBound)
//           {
//             last=RowPtr[l1+1];
//             for(l2=RowPtr[l1];l2<last;l2++)
//             {
//               if(ColInd[l2] == m)
//               {
//                 Entries[l2] += Coupling[l] * v;
//               }
//             }                                     // endfor l2
//           }                                       // endif
//         }                                         // endfor l
//       }                                           // endfor n
//     }                                             // endfor i
//   }                                               // endfor j
// 
//   for(j=0;j<n_matrices;j++)
//   {
//     // hanging nodes in test space
//     fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetTestSpace());
//     N_ = fespace->GetN_Hanging();
//     HangingNodes = fespace->GetHangingNodes();
// 
//     Entries = matrices[j]->GetEntries();
//     RowPtr = matrices[j]->GetRowPtr();
//     ColInd = matrices[j]->GetKCol();
// 
//     CurrentHangingEntries = HangingEntries[j+n_sqmatrices];
//     HangingRowPtr = matrices[j]->GetHangingRowPtr();
//     HangingColInd = matrices[j]->GetHangingKCol();
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     for(i=0;i<N_;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
//       k = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       end = HangingRowPtr[i+1];
//       for(n=HangingRowPtr[i];n<end;n++)
//       {
//         v = CurrentHangingEntries[n];
//         m = HangingColInd[n];
//         for(l=0;l<k;l++)
//         {
//           l1 = DOF[l];
//           if(l1<ActiveBound)
//           {
//             last=RowPtr[l1+1];
//             for(l2=RowPtr[l1];l2<last;l2++)
//             {
//               if(ColInd[l2] == m)
//               {
//                 Entries[l2] += Coupling[l] * v;
//               }
//             }                                     // endfor l2
//           }                                       // endif
//         }                                         // endfor l
//       }                                           // endfor n
//     }                                             // endfor i
// 
//     // hanging nodes in ansatz space
//     N_Rows =  matrices[j]->GetN_Rows();
//     fespace = (TFESpace2D *) (matrices[j]->GetStructure()->GetAnsatzSpace());
// 
//     N_ = fespace->GetN_Hanging();
//     HangingNodes = fespace->GetHangingNodes();
//     AnsatzActiveBound = fespace->GetActiveBound();
//     AnsatzHangingBound = fespace->GetHangingBound();
//     for(i=0;i<N_Rows;i++)
//     {
//       end = RowPtr[i+1];
//       for(k=RowPtr[i];k<end;k++)
//       {
//         l = ColInd[k];
//         if(l>=AnsatzActiveBound && l<AnsatzHangingBound)
//         {
//           // l is hanging node in ansatz space
//           hn = HangingNodes[l-AnsatzActiveBound];
//           HNDescr = hn->GetType();
//           HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
//           m = HNDescr_Obj->GetN_Nodes();
//           Coupling = HNDescr_Obj->GetCoeff();
//           DOF = hn->GetDOF();
//           v = Entries[k];
//           for(n=0;n<m;n++)
//           {
//             for(l1=RowPtr[i];l1<end;l1++)
//             {
//               if(ColInd[l1] == DOF[n])
//                 Entries[l1] += v*Coupling[n];
//             }
//           }
//           Entries[k] = 0;
//         }                                         // endif l
//       }                                           // endfor k
//     }                                             // endfor i
//   }                                               // endfor j
// 
//   for(j=0;j<n_rhs;j++)
//   {
//     fespace = ferhs[j];
//     N_Hanging = fespace->GetN_Hanging();
//     HangingNodes = fespace->GetHangingNodes();
// 
//     RHS = rhs[j];
//     CurrentHangingRhs = HangingRhs[j];
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     for(i=0;i<N_Hanging;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
//       N_ = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       for(k=0;k<N_;k++)
//       {
//         l = DOF[k];
//         if(l<ActiveBound)
//         {
//           RHS[l] += Coupling[k] * CurrentHangingRhs[i];
//         }
//       }                                           // endfor k
//     }                                             // endfor i
//   }                                               // endfor j
// 
//   // ####################################################################
//   // write coupling into matrix
//   // ####################################################################
//   for(j=0;j<n_sqmatrices;j++)
//   {
//     fespace = sqmatrices[j]->GetFESpace();
//     N_ = fespace->GetN_Hanging();
//     HangingNodes = fespace->GetHangingNodes();
// 
//     Entries = sqmatrices[j]->GetEntries();
//     RowPtr = sqmatrices[j]->GetRowPtr();
//     ColInd = sqmatrices[j]->GetKCol();
// 
//     ActiveBound = fespace->GetActiveBound();
// 
//     n = RowPtr[ActiveBound];
// 
//     for(i=0;i<N_;i++)
//     {
//       hn = HangingNodes[i];
//       HNDescr = hn->GetType();
//       HNDescr_Obj = TFEDatabase2D::GetHNDesc2D(HNDescr);
//       k = HNDescr_Obj->GetN_Nodes();
//       Coupling = HNDescr_Obj->GetCoeff();
//       DOF = hn->GetDOF();
// 
//       Entries[n] = 1.0;
//       n++;
// 
//       for(l=0;l<k;l++)
//       {
//         Entries[n] = - Coupling[l];
//         n++;
//       }                                           // endfor l
// 
//     }                                             // endfor i
//   }
// 
//   if(n_sqmatrices)
//   {
//     delete [] GlobalNumbers;
//     delete [] BeginIndex;
//   }
// 
//   if(n_matrices)
//   {
//     delete [] AnsatzGlobalNumbers;
//     delete [] AnsatzBeginIndex;
//     delete [] TestGlobalNumbers;
//     delete [] TestBeginIndex;
//   }
// 
//   if(n_sqmatrices+n_matrices)
//   {
//     for(i=0;i<n_sqmatrices+n_matrices;i++)
//       delete [] HangingEntries[i];
// 
//     delete [] HangingEntries;
//   }
// 
//   if(n_rhs)
//   {
//     for(i=0;i<n_rhs;i++)
//       delete [] HangingRhs[i];
//     delete [] HangingRhs;
// 
//     delete [] righthand;
//     delete [] LocRhs;
//     delete [] RhsBeginIndex;
//     delete [] RhsGlobalNumbers;
//   }
// 
//   if(N_Parameters)
//   {
//     delete [] Param[0];
//   }
// 
//   if(N_AllMatrices)
//   {
//     delete [] LocMatrices;
//     delete [] Matrices[0];
//     delete [] Matrices;
//   }
// 
//   delete [] AuxArray[0];
// }                                                 // end of Assemble
// 
// 
// 
// 
// 
