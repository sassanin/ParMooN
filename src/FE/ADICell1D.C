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
// @(#)ADICell.C        4.1 07.11.09
// 
// Class:       src for TADICell
// Purpose:     general super class for all ADICells
//              special spaces are implemented in subclasses
//
// Author:      Sashikumaar Ganesan (07.11.09)
//
// History:     start of implementation 07.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADICell1D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <DirectSolver.h>
#include <SquareStructure1D.h>
#include <SquareMatrix1D.h>
#include <FEFunction1D.h>
#include <LineAffin.h>
#include <LinAlg.h>

#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdio.h>
#include <stdlib.h>
// #include <malloc.h>

TADICell1D::TADICell1D(TBaseCell *cell,  TFESpace1D *FESpace1D_internal, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
                       int N_quadPts, double *x, DoubleFunct2D *initial,
                       DoubleFunct2D *exact):
            TADICell(cell, N_quadPts)
{
 int i;

  Initial = initial;
  Exact = exact;

  FESpace1D_Internal = FESpace1D_internal;
  N_V = FESpace1D_internal->GetN_DegreesOfFreedom();
  Collection_Internal = FESpace1D_Internal->GetCollection();

  M_Internal = M_internal;
  A_Internal = A_internal;

 // allocate memory and copy the real domain position values
 X = new double[N_quadPts];
 for(i=0; i<N_quadPts; i++)
  X[i] = x[i];

 ConstructAllInfo();
}

int TADICell1D::ConstructAllInfo()
{
 int i, j, N_Cells_Internal;
 double *val;
 char GString[] = "Grid";
 char VString[] = "V";

  // allocate memory for all QuadPts
  FEFunction1D_Internal = new TFEFunction1D*[N_QuadPts];

  sol = new double[N_QuadPts*N_V];
  rhs = new double[N_QuadPts*N_V];
  solT = new double[N_QuadPts*N_V];
  rhsT = new double[N_QuadPts*N_V];

  memset(sol, 0, N_QuadPts*N_V*SizeOfDouble);
  memset(rhs, 0, N_QuadPts*N_V*SizeOfDouble);
  memset(solT, 0, N_QuadPts*N_V*SizeOfDouble);
  memset(rhsT, 0, N_QuadPts*N_V*SizeOfDouble);

  //OutPut("dof inernal space : "<< setw(10) << N_V << endl);
  for(i=0; i<N_QuadPts; i++)
   {
    FEFunction1D_Internal[i] = new TFEFunction1D(FESpace1D_Internal, VString, VString, sol+i*N_V, N_V);

    FEFunction1D_Internal[i]->Interpolate(0, X[i], Initial);
   }

  GridX_internal = new double[N_V];
  memset(GridX_internal, 0, N_V*SizeOfDouble);
  GridFEFunction1D_Internal = new TFEFunction1D(FESpace1D_Internal, GString, GString, GridX_internal, N_V);
  GridFEFunction1D_Internal->GridToData();

//   //transpose of the solution
//   for(i=0; i<N_QuadPts; i++)
//    {
//     val = sol+i*N_V;
//     for(j=0; j<N_V; j++)
//      {
//       solT[j*N_QuadPts + i ] = val[j];
//      }
//    }

 return 0;
}

void TADICell1D::SolveAllQdPts(double *G, double *QuadPtsRhsT, CoeffFct2D *BilinearCoeffs,
                               BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue,
                               double tau, double *Sol_Loc)
{
 int i, j, k;
 double *B, gamma=0., *Sol_Qpt, *OldSol, *defect, BDValue;
 BoundCond BDType;

  B = new double[N_V];
  defect = new double[N_V];

  for(i=0; i<N_QuadPts; i++)
   {
    OldSol = QuadPtsRhsT+i*N_V;

//      for(j=0; j<N_V; j++)
//      cout<< i << " " <<  j<<  " " << OldSol[j] <<endl;

    Sol_Qpt = sol+i*N_V;
    memcpy(Sol_Qpt, OldSol, N_V*SizeOfDouble);

    // if G is constant then it is enough to assemble A once for all QuadPts
    A_Internal->Reset();

    //working rhs array
    memset(B, 0, N_V*SizeOfDouble);
    Daxpy(N_V, tau*TDatabase::TimeDB->THETA3, rhs, B);

    memset(rhs, 0, N_QuadPts*N_V*SizeOfDouble);
    this->AssembleARhs(X[i], G[i], BilinearCoeffs, BoundaryCondition, BoundValue);

    Daxpy(N_V, tau*TDatabase::TimeDB->THETA4, rhs, B);

   // M = M -(tau*TDatabase::TimeDB->THETA2) A
    MatAdd(M_Internal, A_Internal,  -tau*TDatabase::TimeDB->THETA2);
    // set current factor of steady state matrix
    gamma = -tau*TDatabase::TimeDB->THETA2;

    // defect = M * sol
    // B:= B + defect
    memset(defect, 0, N_V*SizeOfDouble);
    MatVectActive(M_Internal, OldSol, defect);

//     for(j=0; j<N_V; j++)
//      cout<< OldSol[j] << " defect " << i << " " <<  j<<  " " << defect[j] <<endl;

    Daxpy(N_V, 1, defect, B);

  //***************************************************************************************
    // set Dirichlet values
    // starting point

    // starting point
    BoundaryCondition(0, X[i], BDType);
    BoundValue(0, X[i], BDValue);

    if(BDType==DIRICHLET)
     {
      B[0] = BDValue;
      Sol_Qpt[0] = BDValue;
     }

     //endpoint
    BoundaryCondition(1, X[i], BDType);
    BoundValue(1, X[i], BDValue);

    if(BDType==DIRICHLET)
     {
      B[N_V-1] = BDValue;
      Sol_Qpt[N_V-1] = BDValue;
     }
    //***************************************************************************************

    // system matrix
    MatAdd(M_Internal, A_Internal, -gamma + tau*TDatabase::TimeDB->THETA1);
    // set current factor of steady state matrix
    gamma = tau*TDatabase::TimeDB->THETA1;


    DirectSolver(M_Internal, B, Sol_Qpt);

    MatAdd(M_Internal, A_Internal, -gamma);
    // set current factor of steady state matrix
    gamma = 0;

 
   }//for(i=0; i<N_QuadPts; i++)

  memcpy(Sol_Loc, sol, N_QuadPts*N_V*SizeOfDouble);
//   memcpy(Sol_Loc, QuadPtsRhsT, N_QuadPts*N_V*SizeOfDouble);


 delete [] B;
}


/**  Assembling A matrix */
/** mass mat is same for all quad points in this cell */
void TADICell1D::AssembleARhs(double x, double Conv, CoeffFct2D *Bilinear, BoundCondFunct2D *BoundaryCondition, BoundValueFunct2D *BoundValue)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, X[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **aux, **coeff, *Coeff;
 double Space_X[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocRhs[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, BDValue;
 bool Needs2ndDer[1];

 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];
 BoundCond BDType;

  GlobalNumbers = FESpace1D_Internal->GetGlobalNumbers();
  BeginIndex = FESpace1D_Internal->GetBeginIndex();

  RowPtr = A_Internal->GetRowPtr();
  KCol = A_Internal->GetKCol();
  ValuesA = A_Internal->GetEntries(); 

  N_Cells_Internal = Collection_Internal->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    Space_X[i] = x;
    aux[i] = new double[6];
    coeff[i] = new double[6];
   }

  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Collection_Internal->GetCell(i);
    FEId = FESpace1D_Internal->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, X, AbsDetjk);

    Bilinear(N_Points, Space_X, X, aux, coeff);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0];
      g0 = Coeff[1];
      c1 = Coeff[2];
      rhsval = Coeff[3];

      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        // cout<< " uref " << origvaluesD1[j][k] << " orgD0 " << orgD1[k] <<endl;
        LocRhs[k] += rhsval*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];

          val  = c0*ansatz1*test1;
          val += g0*ansatz1*test0;
          val += c1*ansatz0*test0;
          LocMatrixA[k*N_BaseFunct + l] += Mult*val;
        }
       }
     }

    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; // soluble surfactant relation

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct;l++)
        {
         if(KCol[k] == DOF[l])
          {
           ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
   }// for(i=0; i<N_Cells_Internal; i++)

 //update boundary data
 // starting point
  BoundaryCondition(0, x, BDType);
  BoundValue(0, x, BDValue);

  if(BDType==DIRICHLET)
   {
    rhs[0] = BDValue;

    begin = RowPtr[0];
    end = RowPtr[1];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == 0 )
       ValuesA[k] = 1.;
      else
       ValuesA[k] = 0.;
    }
   }
  else if(BDType== NEUMANN)
   {
    rhs[0] += BDValue;
   }
  else
   {
    OutPut("Unknown internal boundary condition !"<< endl); 
    exit(0);
   }

 // end point
  BoundaryCondition(1, x, BDType);
  BoundValue(1, x, BDValue);

  if(BDType==DIRICHLET)
   {
    rhs[N_V-1] = BDValue;

    begin = RowPtr[N_V-1];
    end = RowPtr[N_V];

    for(k=begin;k<end;k++)
     {
      if(KCol[k] == N_V-1 )
       ValuesA[k] = 1.;
      else
       ValuesA[k] = 0.;
     }
   }
  else if(BDType== NEUMANN)
   {
    rhs[N_V-1] += BDValue;
   }
  else
   {
    OutPut("Unknown internal boundary condition !"<< endl); 
    exit(0);
   }


// //print matrix
//   for(j=0;j<N_V;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//      }
//     cout<<endl;
//    }
// exit(0);

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    delete [] aux[i];
    delete [] coeff[i];
   }

  delete [] aux;
  delete [] coeff;
}


TADICell1D::~TADICell1D()
{
 int i;

 for(i=0; i<N_QuadPts; i++)
  if(FEFunction1D_Internal[i]) delete FEFunction1D_Internal[i];

  if(sol) delete [] sol;
  if(rhs) delete [] rhs;
  if(FEFunction1D_Internal) delete [] FEFunction1D_Internal;

//   delete FESpace1D_Internal;
//   delete SqStruct1D_Internal;
//   delete M_Internal;
//   delete A_Internal;

}
