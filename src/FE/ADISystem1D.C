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
// @(#)ADISystem1D.C        4.1 13.11.09
// 
// Class:       TADISystem1D
// Purpose:     class for  ADISystem1D

//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#include <ADISystem1D.h>
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

// 2D
TADISystem1D::TADISystem1D(TFEFunction1D *FEFunction_intl, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
                           TSquareMatrix1D *S_internal, TSquareMatrix1D *K_internal,
                           double x, double y, DoubleFunctND *initial,
                           DoubleFunct2D *bDValue_LMin, DoubleFunct2D *bDVal_LMax,
                           double *Sol,  double *OldSol, double *b, 
                           double *Defect, double *Intlposl, DoubleFunctND *growthAndNuc): TADISystem(Sol, OldSol, b, Defect, growthAndNuc)
{
//  int i;

  FEFunction_Intl = FEFunction_intl;
  FESpace1D_Intl = FEFunction_Intl->GetFESpace1D();
  N_V = FESpace1D_Intl->GetN_DegreesOfFreedom();
  Collection_Intl = FESpace1D_Intl->GetCollection();

  M_Intl = M_internal;
  A_Intl = A_internal;

  S_Intl = S_internal;
  K_Intl = K_internal;

  //  ARowPtr = A_Intl->GetRowPtr();
//  AKCol = A_Intl->GetKCol();
//  AMatValues = A_Intl->GetEntries(); 
  N_MMatValues = M_Intl->GetN_Entries();
   
  Initial = initial;
  BDValue_LMin = bDValue_LMin;
  BDVal_LMax = bDVal_LMax;

  
 // allocate memory and copy the real domain position values
 X = x;
 Y = y;
 Z=0;
 
  IntlPosL = Intlposl;
 
  //set the rhs
  this->ConstructAllInfo();
}


// 3D
TADISystem1D::TADISystem1D(TFEFunction1D *FEFunction_intl, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
                           double x, double y, double z, double *Sol,  double *OldSol, double *b,
                           double *Defect, DoubleFunctND *growthAndNuc): TADISystem(Sol, OldSol, b, Defect, growthAndNuc)
{
 int i;

  Initial = NULL;

  FEFunction_Intl = FEFunction_intl;
  FESpace1D_Intl = FEFunction_Intl->GetFESpace1D();
  N_V = FESpace1D_Intl->GetN_DegreesOfFreedom();
  Collection_Intl = FESpace1D_Intl->GetCollection();

  M_Intl = M_internal;
  A_Intl = A_internal;

 // allocate memory and copy the real domain position values
 X = x;
 Y = y;
 Z = z;
 
 /** set the rhs */
 this->ConstructAllInfo();
}

int TADISystem1D::ConstructAllInfo()
{
 int i, j, N_Cells_Internal;
 double *val;
 char GString[] = "Grid";
 char VString[] = "V";

  rhs = new double[N_V];
  memset(rhs, 0, N_V*SizeOfDouble);

 return 0;
}

// void TADISystem1D::Interpolate(double *Sol, DoubleFunct3D *Exact)
// {
//  int N_Coord = 2;
//  double Coords[3]; // 3+1
//  Coords[0] = X;
//  Coords[1] = Y;
//  
// //   FEFunction_Intl->InterpolateNodalPts(X, Y, Exact, Sol);
// }

void TADISystem1D::Interpolate(double *Sol, DoubleFunctND *Exact)
{
 int N_Coord = 2;
 double Coords[4]; // 3+1
 
 Coords[0] = X;
 Coords[1] = Y;
  
#ifdef __3D__  
 N_Coord++;
 Coords[2] = Z;  
#endif

  FEFunction_Intl->InterpolateNodalPts(N_Coord, Coords, Exact, Sol);
}
// ====================================================================
//  compute defect for the given system (OpenMP specification)
// ====================================================================
void Defect_Scalar(TSquareMatrix *A, double *x, double *b, double *r)
{
 int i, j, k, l, N_DOF, *RowPtr, *KCol;
 int begin, end, index;

 double rhs_i, *Entries;

  N_DOF = A->GetN_Rows();
  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  Entries = A->GetEntries();

  for(i=0;i<N_DOF;i++)
   {
    rhs_i = b[i];
    begin = RowPtr[i];
    end = RowPtr[i+1];

    for(j=begin;j<end;j++)
     {
      index = KCol[j];
      rhs_i -= (Entries[j] * x[index]);
     }
    r[i] = rhs_i;
   }
}


// double GetGrowth(double C, double T)
// {
//   double PbeGC = TDatabase::ParamDB->REACTOR_P27;
//   double G, C_Sat;
//   double T_ref = TDatabase::ParamDB->REACTOR_P23;
//   double C_ref = TDatabase::ParamDB->REACTOR_P25;
// 
//    C_Sat = (1.3045*(T_ref*T - 273.15) + 35.3642)/C_ref;
// 
//     // growth rate based on concentration
//     if(C>C_Sat && (TDatabase::ParamDB->PBE_P0==1 || TDatabase::ParamDB->PBE_P0==4|| TDatabase::ParamDB->PBE_P0==6))
//      { 
//        G = PbeGC*sqrt((C - C_Sat )/C_Sat);defect
//      }
//     else
//      { G = 0.; }
// 
//      if(G<0.) G = 0.;
// 
// //      if(G>0)
// //       cout << " G " << G << endl;  
//  return G;
//  }


void TADISystem1D::SetDirichletBc(BoundCond cond_Lmin, BoundCond cond_Lmax, double BDValue_Lmin, double BDValue_Lmax)
{
 int k;
 int *RowPtr, *KCol, begin, end;
 double *MatValues;

 RowPtr = M_Intl->GetRowPtr();
 KCol = M_Intl->GetKCol();
 MatValues = M_Intl->GetEntries(); 


    if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
     {
       if(cond_Lmin==DIRICHLET)
       {
        begin = RowPtr[0];
        end = RowPtr[1];

        for(k=begin;k<end;k++)
         {
          if(KCol[k] == 0 )
           { MatValues[k] = 1.; }
          else
           { MatValues[k] = 0.; }
         }
       }

      if(cond_Lmax==DIRICHLET)
       {
        begin = RowPtr[N_V-1];
        end = RowPtr[N_V];

        for(k=begin;k<end;k++)
         {
          if(KCol[k] == N_V-1)
           {
            MatValues[k] = 1.;
           }
          else
           {
            MatValues[k] = 0.;
           }
         }
       }
      }

}
// sol will be returned in Sol
void TADISystem1D::SolveAllLevels(double *MatValues_orig, double *SolFromX,  double *AggrRhs, double &C, double C_Sat, 
                                  double &T, CoeffFctND *BilinearCoeffs,
                                  double tau, BoundCond cond_Lmin, BoundCond cond_Lmax, DoubleFunct3D *Exact)
{
 int i, j, k, l;
//  int *RowPtr, *KCol, begin, end;
//  int *ARowPtr, *AKCol, dGDisc=FESpace1D_Intl->IsDGSpace();
 int N_NonLinearSteps = (int)TDatabase::ParamDB->REACTOR_P4;

 double gamma=0., BDValue1, BDValue2, oldC, oldT, tempC, tempT;
 double G, residual_pbe, F, Inn[2], Out[2], Bnuc=0.;
 double *MatValues, limit=1.e-8;
//  double  *AMatValues, *SMatValues, *KMatValues;
 double GC_C = TDatabase::ParamDB->REACTOR_P26;
 double GC_T = TDatabase::ParamDB->REACTOR_P28;
 double rhsfact_C = TDatabase::ParamDB->REACTOR_P6;
 double rhsfact_T = TDatabase::ParamDB->REACTOR_P7;
 double G_DimLessFact = TDatabase::ParamDB->REACTOR_P27;
 
  MatValues = M_Intl->GetEntries(); 

  // PSD at L_{min}
  if(cond_Lmin==DIRICHLET && TDatabase::ParamDB->INTL_DISCTYPE!=DG)
   BDValue1 = SolFromX[0];

//     // PSD at L_{max}
  if(cond_Lmax==DIRICHLET && TDatabase::ParamDB->INTL_DISCTYPE!=DG)
   BDValue2 = SolFromX[N_V -1];

   /** store the original mass mat  values*/
   memcpy(MatValues,  MatValues_orig,  N_MMatValues*SizeOfDouble);

  /** set f as non negative */
//   for(j=0;j<N_V;j++)
//     {  if(SolFromX[j]<0.) SolFromX[j]=0. ; }
  
  
  /** save the t^n values \tilde C^{n-1}=\hat C^n, \tilde T^{n-1}=\hat T^n  before the iteration starts*/
  oldC = C;
  oldT = T;

  for(l=0; l<N_NonLinearSteps; l++)
   {
    A_Intl->Reset();

    if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
     {
      S_Intl->Reset();
      K_Intl->Reset();
     }

    /**growth rate based on concentration defined in Example file*/
    if(GetGrowthAndNuc)
    {
      Inn[0] = C;
      Inn[1] = T;

      GetGrowthAndNuc(2, Inn, Out);

      G = G_DimLessFact*Out[0];
      Bnuc = Out[1];

    }
    else
    {
     G = 0.;
     Bnuc = BDValue1;
    }
//       cout << "G  : " << G <<" Bnuc " << Bnuc <<endl; 
    /** init rhs */
    memset(rhs, 0, N_V*SizeOfDouble);

     /** assemble matrices */
    switch(TDatabase::ParamDB->INTL_DISCTYPE)
     {
       case FD:
        // Finite difference, Il'in-Allen-Southwell scheme
           this->AssembleARhs_FD(G, BilinearCoeffs);
       break;      
       case  GALERKIN:
        // Galerkin
            this->AssembleARhs(G, BilinearCoeffs);
       break;  
       case SDFEM:
        // SUPG
            this->AssembleARhs_SUPG(G, BilinearCoeffs); 
       break;  
       case DG:
        // DG
            this->AssembleARhs_DG(G, Bnuc, BilinearCoeffs, cond_Lmin, cond_Lmax);  
       break;
       default:
            Error("only FD, GALERKIN, SUOG, DG are implemented" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
       break;
     }

   if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
    {
     if(cond_Lmin==DIRICHLET)
      {
       rhs[0] = 0.; // nucleation will be copied from oldsol to B
       sol[0] = Bnuc;
       SolFromX[0] = Bnuc;
       //cout << " Bnuc " << Bnuc << endl;
      }

     if(cond_Lmax==DIRICHLET)
      {
//        cout << " Not Verified Solve Internal all " << endl;
//        exit(0);

       rhs[N_V -1] = BDValue2;
       sol[N_V -1] = BDValue2;
       SolFromX[N_V -1] = BDValue2;
      }
    }

   /** store A and other mat values from previous time step for rhs calculation */
   if(l==0)
    {
      //SUPG, time consistant term
      if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
       {
        MatAdd(M_Intl, S_Intl,  1.);
       }
      
     if(TDatabase::TimeDB->THETA2!=0.)
      {
        // M = M -(tau*TDatabase::TimeDB->THETA2) A
        MatAdd(M_Intl, A_Intl,  -tau*TDatabase::TimeDB->THETA2);
        //SUPG
        if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
         {
          MatAdd(M_Intl, K_Intl,  -tau*TDatabase::TimeDB->THETA2);
         }// if(TDatabase::ParamDB->INTL_DI
      } //if(TDatabase::TimeDB->
      
     this->SetDirichletBc(cond_Lmin, cond_Lmax, Bnuc, BDValue2);
     // defect = M * oldsol
     // B:= B + defect
      //      /** \tilde f^{n-1}=\hat f^n,  , that is, SolFromX at l==0 */
     memset(oldsol, 0, N_V*SizeOfDouble);
     MatVectActive(M_Intl, SolFromX, oldsol);

     /**if aggr need to be considered then Aggr function has to be called within nonlinear loop */
     /**so, aggregation is treated explicitly, since the aggr routines are expensive */
     /** make sure dt is multiplied with Aggr in main prg */
     Daxpy(N_V, 1., AggrRhs,  oldsol);

     //      /**   ( aggregaion,  b\cdot\grad v) */
     if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
       {
        memcpy(sol,  AggrRhs,  N_V*SizeOfDouble);

        if(TDatabase::ParamDB->INTL_DISCTYPE!=DG && cond_Lmin==DIRICHLET)
         {
          sol[0] = Bnuc;
         }

        memset(B, 0, N_V*SizeOfDouble);
        MatVectActive(S_Intl,  sol, B);

       /** add it to old rhs */
        Daxpy(N_V, 1., B,  oldsol);
       }

     /** oldsol now contains : dt*A_{agg} + [M + S^{supg} - dt\theta_2(A+K^{supg})]f^{n-1} + dt*(supg aggr) */
     if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
      {
       if(cond_Lmin==DIRICHLET)
        {
         oldsol[0] = Bnuc;
        }
      if(cond_Lmax==DIRICHLET)
        {
         oldsol[N_V -1] = BDValue2;
        }
      }// !DG
    }// if(l==0)

    /** working rhs array */
    memset(B, 0, N_V*SizeOfDouble);

    /** old rhs */
    Daxpy(N_V, 1., oldsol, B);

    /** usually zero, as no other term except Aggr, so no need dt*[theata3 + theta4] */
    Daxpy(N_V, tau, rhs, B);

    /** reset mass matrix */
    memcpy(MatValues,  MatValues_orig,  N_MMatValues*SizeOfDouble);

    /** system matrix */
    MatAdd(M_Intl, A_Intl, tau*TDatabase::TimeDB->THETA1);

    /** add SUPG, time and space matrices */
    if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
    { 
      MatAdd(M_Intl, S_Intl,  1.);
      MatAdd(M_Intl, K_Intl, tau*TDatabase::TimeDB->THETA1);
//       	        cout << "G ADIS: " << G <<" C : " << C <<" T: " << T <<endl; 
    }

    /** set the Diriclet BC */
    this->SetDirichletBc(cond_Lmin, cond_Lmax, Bnuc, BDValue2);


    // check the convergence of the fixed point linerization
    if(l>0 )
     {
      memset(defect, 0, N_V*SizeOfDouble);
      Defect_Scalar(M_Intl, sol,  B, defect);
      residual_pbe =  Ddot(N_V, defect, defect);

      if(sqrt(residual_pbe)<limit)
       {
        break;
       }
     } // if( l>0 )

    /** set DIRICHLET value */
   if(TDatabase::ParamDB->INTL_DISCTYPE!=DG)
    {
     if(cond_Lmin==DIRICHLET)
      {
       rhs[0] = Bnuc;
       sol[0] = Bnuc;
       B[0]   = Bnuc;
      }

     if(cond_Lmax==DIRICHLET)
      {
       rhs[N_V -1] = BDValue2;
       sol[N_V -1] = BDValue2;   
       B[N_V -1]   = BDValue2;
      }
    }

   DirectSolver(M_Intl, B, sol);

    /** reset the mass matrix for next nonlinear iteration just before the system mat assembling using orig M mat values*/
    memcpy(MatValues,  MatValues_orig,  N_MMatValues*SizeOfDouble);


   //update C
   if((C>0. || T>0.) && (GC_T!=0. || GC_T!=0.) )
    {
     F = this->GetWeightedF();
     if(F<0.) F=0.;

     if(GetGrowthAndNuc)
      {
       Inn[0] = C;
       Inn[1] = T;
       GetGrowthAndNuc(2, Inn, Out);
       G = Out[0]; 
       Bnuc = Out[1];
      }
    
     tempT = oldT - tau*(GC_T*G*F + rhsfact_T*Bnuc);
     if(tempT>0.)// this condition is necessay, if Bnuc is a constant
      T = tempT;

     tempC = oldC - tau*(GC_C*G*F + rhsfact_C*Bnuc);    
     if(tempC>0.) // this condition is necessay, if Bnuc is a constant
      C = tempC;

//      if((tau*GC_C*G*F)>1.e-2)
//       cout <<"C " <<  C << " F " <<  F <<endl;
    }
   } // for(l=0; l<N_NonLinearSteps; l++)

//   if(Exact)
//    FEFunction_Intl->Interpolate(X, Y,  Exact);

// if(TDatabase::TimeDB->CURRENTTIME>0.01)
//     for(j=0;j<N_V;j++)
//     cout << "SolFromX" << j << ": "<< SolFromX[j]<< " sol " << ": "<< sol[j] <<endl;
//  exit(0);

    memcpy(SolFromX, sol, N_V*SizeOfDouble);

//     if(TDatabase::TimeDB->THETA2!=0.)
//      {
//       memcpy(AMatValues,  A_Intl->GetEntries(),  N_AMatValues*SizeOfDouble);
// 
//       if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
//        {
//         memcpy(SMatValues,  S_Intl->GetEntries(),  N_SMatValues*SizeOfDouble);
//         memcpy(KMatValues,  K_Intl->GetEntries(),  N_KMatValues*SizeOfDouble);
//        }
//      }

//    // M = M -(tau*TDatabase::TimeDB->THETA2) A
//     if(TDatabase::TimeDB->THETA2!=0)
//       MatAdd(M_Intl, A_Intl,  -gamma -tau*TDatabase::TimeDB->THETA2);
// 
//       //SUPG
//       if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
//        {
// //       already sub, so needs to be added
//          MatAdd(M_Intl, S_Intl,  1.);
// 
//         if(TDatabase::TimeDB->THETA2!=0)
//           MatAdd(M_Intl, K_Intl,  -gamma -tau*TDatabase::TimeDB->THETA2);
//        }
// 
//     gamma = -tau*TDatabase::TimeDB->THETA2; 
//     // AggrRhs = M * oldsol
//     // B:= B + defect
//     memset(AggrRhs, 0, N_V*SizeOfDouble);
//     MatVectActive(M_Intl, sol, AggrRhs);
//  
// 
//         //restore mass matrix, since M is same for all Intl Pts
//         MatAdd(M_Intl, A_Intl, -gamma);
// 
//         //SUPG
//         if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
//          {
//           MatAdd(M_Intl, S_Intl,  -1.); 
//           MatAdd(M_Intl, K_Intl,  -gamma);
//          } 
//     gamma = 0;
// //   print matrix
//   for(j=0;j<N_V;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
// //       cout << "M(" << j << ", "<< KCol[k] << ") = " << MatValues[k] <<endl;
// //       cout << "A(" << j << ", "<< KCol[k] << ") = " << AMatValues[k] <<endl;
//       cout << "K(" << j << ", "<< KCol[k] << ") = " << KMatValues[k] <<endl;
// //       cout << "S(" << j << ", "<< KCol[k] << ") = " << SMatValues[k] <<endl;
//      }
//      cout <<endl;
//    }
// exit(0);

} // void TADISystem1D::SolveAllLevels(


double TADISystem1D::GetWeightedF()
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double g=0., L, val=0.;
 double *Weights, *zeta, Intl_L[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **origvaluesD0, Mult;
 double *orgD0, test0;
 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  N_Cells_Internal = Coll->GetN_Cells();


  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
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
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      L = Intl_L[j];

//       len +=Mult;

      val=0.;
      for(k=0;k<N_BaseFunct;k++)
       val += sol[DOF[k]]*orgD0[k];

       g +=L*L*Mult*val;

//       cout << " L " << L << " len " << len << endl;

//        if( sol[0]>0.9)
//        cout << " L*L " << L*L << " g " << g << " sol[DOF[k]] " << sol[0] << endl;
      
     } //  for(j=0;j<N_Points;j++)
   }

// exit(0);
 return g;
}


void TADISystem1D::AssembleInitRhs(double C, double T, CoeffFctND *BilinearCoeffs, 
                                   BoundCond cond_Lmin, BoundCond cond_Lmax)
{
//   int i;
//    double gamma, G, Inn[2], Out[2];
// 
//    /** A needed for rhs, other than Euler schemes */
//     A_Intl->Reset();
// 
//     if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
//      {
//       S_Intl->Reset();
//       K_Intl->Reset();
//      }
// 
//    if(TDatabase::TimeDB->THETA2!=0.)
//     {
//     /**growth rate based on concentration defined in Example file*/
//     if(GetGrowthAndNuc)
//     {
//       Inn[0] = C;
//       Inn[1] = T;
// 
//       GetGrowthAndNuc(2, Inn, Out);
//       G = Out[0];
//     }
//     else
//     {G = 0.;}
// 
//      switch(TDatabase::ParamDB->INTL_DISCTYPE)
//       {
//        case FD:
//         // Finite difference, Il'in-Allen-Southwell scheme
//             cout << "Disc type not implemented yet !!!!!!!!!!!!!" << endl;
//             exit(0);
//        break;
//        case  GALERKIN:
//         // Galerkin
//             this->AssembleARhs(G, BilinearCoeffs);
//        break;
//        case SDFEM:
//         // SUPG
//             this->AssembleARhs_SUPG(G, BilinearCoeffs); 
//        break;
//        case DG:
//         // DG
//             this->AssembleARhs_DG(G, BilinearCoeffs, cond_Lmin, cond_Lmax);  
//        break;
//        default:
//             Error("only FD, GALERKIN, SUOG, DG are implemented" << endl);
//             Error("file: " << __FILE__ << " line " << __LINE__ << endl);
//             exit(-1);
//        break;
//       }
// 
//    /** store the mat for the first time step */
//      memcpy(AMatValues,  A_Intl->GetEntries(),  N_AMatValues*SizeOfDouble);
// 
//      if(TDatabase::ParamDB->INTL_DISCTYPE==SUPG)
//       {
//        memcpy(SMatValues,  S_Intl->GetEntries(),  N_SMatValues*SizeOfDouble);
//        memcpy(KMatValues,  K_Intl->GetEntries(),  N_KMatValues*SizeOfDouble);
//       }
//     }
}


// /**  Assembling A matrix Ilen-Southwell finite difference matrix*/
// /** mass mat is same for all quad points in this cell */
void TADISystem1D::AssembleARhs_FD(double Conv, CoeffFctND *Bilinear)
{
  int i, j, N;
  int begin, end, *RowPtr, *KCol; 
  double *ValuesA, a, b, c, h, eps, sigma, q, val1, val2;
  
  if(TDatabase::ParamDB->REACTOR_P3)
   { eps = 1.0/TDatabase::ParamDB->REACTOR_P3;}
  else
   { eps = 1.e-12;}
  
//   tildeps = eps ;
// 
//   sigma = 
//   a = Conv - eps;  
  
  
  
  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 
  N = A_Intl->GetN_Columns();
  
   for(i=0;i<N;i++)
    {
      if(i<(N-1))
       { h = IntlPosL[i+1] - IntlPosL[i]; }
      
      // Ilin-Allen Southwell scheme
      if(fabs(Conv)>0.)
       {
        q = Conv*h/(2.*eps);
        sigma = q/tanh(q);
      
//          cout << " Conv " << Conv/h << " sigma " << eps *sigma/(h*h) << endl;
       }
      else
      {
       sigma = 0.;
      }
       
      begin=RowPtr[i];
      end=RowPtr[i+1];
      
//       val1 = eps *sigma/(h*h);
//       val2 = Conv/h;
      
      
      //BD 
      val1 = Conv/h;
      
//       if(TDatabase::TimeDB->CURRENTTIME>0.01)
//        cout << " i_prev:  " <<  -val1 - val2 << " i:  " << 2*val1<< " i++ " <<  -val1+ val2 <<endl; 
      
      
      for(j=begin;j<end;j++)
       {
//        if( i==(KCol[j]-1) )
//         { ValuesA[j] = -val1 - val2; }
//        else if(i==(KCol[j]))
//         { ValuesA[j] =  2*val1;}	 
//        else if(i==(KCol[j]+1))
//         { ValuesA[j] = -val1 + val2 ; }
//         BD
       if( i==(KCol[j]-1) )
        { ValuesA[j] = -val1; }
       else if(i==(KCol[j]))
       {ValuesA[j] =  val1;}
	
      }   // endfor j
    }    
   
  
  
}


// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TADISystem1D::AssembleARhs(double Conv, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct, N_Coord=2;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocRhs[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, len=0.;
 double **aux, **coeff, *Coeff;
 double **Coords;
 
#ifdef __3D__  
 N_Coord = 3;
#endif
 

 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

  Coords = new double* [N_Coord+1];

  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];

  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    Coords[0][i] = X;
    Coords[1][i] = Y;
#ifdef __3D__  
    Coords[2][i] = Z;
#endif

    aux[i] = new double[6];
    coeff[i] = new double[6];
   }

#ifdef __3D__  
   Intl_L = Coords[3];
#else
   Intl_L = Coords[2];
#endif


  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
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
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

    Bilinear(N_Points, N_Coord, Coords, aux, coeff); 

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
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      //cout<< " c  " << c0 << " g " <<g0 << " c1 " << c1 <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        LocRhs[k] += Mult*rhsval*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];

          val  = Conv*ansatz1*test0;
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;
          //  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA[k*N_BaseFunct + l] += (Mult*val);
        } 
       } 
     }
//      cout<<endl;
    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 
      
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
   
// cout<< " len " << len << endl;
// exit(0);
// //print matrix
//   for(j=0;j<N_V;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//      }
//      
//       cout << "f: " << rhs[j ] <<endl;
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

  
  for(i=0; i<N_Coord+1; i++)  
   delete [] Coords[i]; 
  
 delete [] Coords;  
  
//  printf("    AssembleARhs  \n" );
//   MPI_Finalize();
// exit(0);
} // void TADISystem1D::AssembleARhs(d


// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TADISystem1D::AssembleARhs_SUPG(double Conv, CoeffFctND *Bilinear)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct, N_Coord=2;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixS[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixK[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D]; 
 double LocRhs[MaxN_BaseFunctions1D];
 double **origvaluesD0, **origvaluesD1, **origvaluesD2, Mult;
 double *orgD0, *orgD1, *orgD2, test0, test1, ansatz0, ansatz1, ansatz2, *ValuesA, *ValuesS, *ValuesK;
 double c0, c1, g0, rhsval, val, len=0., bgradv, hE, beta, Pe_K;
 double **aux, **coeff, *Coeff, **Coords;
 
 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];
 TVertex *Vetrex1, *Vetrex2;

  double D_L =  TDatabase::ParamDB->P11;
  double delta0 = TDatabase::ParamDB->P12;
  double delta1 = TDatabase::ParamDB->P13;
  
  if(D_L<1.e-12) D_L = 1.e-12;
  
  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 
  ValuesS = S_Intl->GetEntries();
  ValuesK = K_Intl->GetEntries();
  
  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = TRUE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

#ifdef __3D__  
   N_Coord=3;
#endif 

  Coords = new double* [N_Coord+1]; 

  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];
   
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    Coords[0][i] = X;
    Coords[1][i] = Y;
#ifdef __3D__  
    Coords[2][i] = Z;
#endif

    aux[i] = new double[6];
    coeff[i] = new double[6];
   }

#ifdef __3D__  
   Intl_L = Coords[3];
#else
   Intl_L = Coords[2];
#endif

  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();

    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

//     Bilinear(N_Points, Space_X, Space_Y, Intl_L, aux, coeff);
    Bilinear(N_Points, N_Coord, Coords, aux, coeff);
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);
    origvaluesD2=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D2);

    memset(LocMatrixA, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixS, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixK, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];
    Vetrex1 = Cell->GetVertex(0);
    Vetrex2 = Cell->GetVertex(1);
    hE = fabs(Vetrex1->GetX() - Vetrex2->GetX());
//     cout<< " hE  " << hE  <<endl;

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

      if(fabs(Conv)>0)
       {
        Pe_K = hE*fabs(Conv)/(2.*D_L);

// 	beta based on Lutz book
//         if(Pe_K>1.)
//           beta = delta0 * hE/Conv;
//         else
//           beta = delta1 *hE*hE/D_L ;
// 	
        // beta based on 1d Green's formula
        beta = hE*(1./tanh(Pe_K) - 1./Pe_K)/(2.*fabs(Conv));
       }
      else
      { beta = 0.0; }

//       cout<< " Pe_K  " << Pe_K  <<" beta  " << beta  <<endl; 
//       cout<< " c  " << c0 << " g " <<g0 << " c1 " << c1 <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];
      orgD2 = origvaluesD2[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

	bgradv = Conv*test1;
// 	if(k==0)
// 	  cout<< " bgradv  " << bgradv  <<" bgradv*beta  " << bgradv*beta  <<endl; 
        LocRhs[k] += Mult*rhsval*(test0 + beta*bgradv);

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];
          ansatz2  = orgD2[l];
	  
          val  = Conv*ansatz1*test0;
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;
          //  cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA[k*N_BaseFunct + l] += (Mult*val);

          val  = (-c0*ansatz2 + Conv*ansatz1)*beta*bgradv;
          LocMatrixK[k*N_BaseFunct + l] += (Mult*val);
//            cout<< "S ( " << k <<" , " <<l << " ) " <<bgradv  <<endl;
          val  = ansatz0*beta*bgradv;
          LocMatrixS[k*N_BaseFunct + l] += (Mult*val);

        }
       }
     }
//      cout<<endl;
//  exit(0);
    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 
      
      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
       for(l=0;l<N_BaseFunct;l++)
        {
         if(KCol[k] == DOF[l])
          {
           ValuesA[k] +=LocMatrixA[j*N_BaseFunct + l];
           ValuesK[k] +=LocMatrixK[j*N_BaseFunct + l];
           ValuesS[k] +=LocMatrixS[j*N_BaseFunct + l];
           break;
          }
        } // for(m=0;m<N_BaseFunct_low
      } // for(n=begin;n<end;n++)
     } // for(l=0;l<N_BaseFunct_low
   }// for(i=0; i<N_Cells_Internal; i++)
   
// cout<< " len " << len << endl;
// exit(0);

// // print matrix
//   for(j=0;j<N_V;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
// //       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//       cout << "K(" << j << ", "<< KCol[k] << ") = " << ValuesK[k] <<endl;
// //       cout << "S(" << j << ", "<< KCol[k] << ") = " << ValuesS[k] <<endl;
//      }
//      
// //       cout << "f: " << rhs[j ] <<endl;
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
  
  for(i=0; i<N_Coord+1; i++)  
   delete [] Coords[i]; 
  
 delete [] Coords;    
  
} // void TADISystem1D::AssembleARhs(d



// /**  Assembling A matrix */
// /** mass mat is same for all quad points in this cell */
void TADISystem1D::AssembleARhs_DG(double Conv, double Bnuc, CoeffFctND *Bilinear, BoundCond cond_Lmin, BoundCond cond_Lmax)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct, N_Coord=2;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF, *NeibDOF;
 int TestDOF, begin, end, *RowPtr, *KCol; 
 int m, N_Joints, Neigh_i;
 
 double rec_detjk, Neigh_rec_detjk, Neigh_N_BaseFunct; 
 double *Weights, *zeta, *Intl_L, AbsDetjk[MaxN_QuadPoints_1D];
 double LocMatrixA11[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA12[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA21[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];
 double LocMatrixA22[MaxN_BaseFunctions1D*MaxN_BaseFunctions1D];  
 double LocRhs[MaxN_BaseFunctions1D];
 double BdValues[3];
 double **origvaluesD0, **origvaluesD1, Mult;
 double *orgD0, *orgD1, test0, test1, ansatz0, ansatz1, *ValuesA;
 double c0, c1, g0, rhsval, val, len=0.;
 double **aux, **coeff, *Coeff, **Coords;
 double JointValues[MaxN_BaseFunctions1D], Neigh_JointValues[MaxN_BaseFunctions1D];
 double JointValuesX[MaxN_BaseFunctions1D], Neigh_JointValuesX[MaxN_BaseFunctions1D];
 double Epsilon, Sigma0, Sigma1, h_max, h1, h2;

 bool Needs2ndDer[1], UpdateEdge=FALSE;

 Epsilon = TDatabase::ParamDB->DG_P0;
 Sigma0 = TDatabase::ParamDB->DG_P1;
 Sigma1 = TDatabase::ParamDB->DG_P2;

 TCollection *Coll;
 TBaseCell *Cell, *Neigh=NULL;
 FE1D FEId, Neigh_FEId;
 TFE1D *Element, *Neigh_Element;
 TBaseFunct1D *bf, *Neigh_bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  RowPtr = A_Intl->GetRowPtr();
  KCol = A_Intl->GetKCol();
  ValuesA = A_Intl->GetEntries(); 

  N_Cells_Internal = Coll->GetN_Cells();
  Needs2ndDer[0] = FALSE;

  aux = new double *[MaxN_QuadPoints_1D];
  coeff = new double *[MaxN_QuadPoints_1D];

#ifdef __3D__  
    N_Coord = 3;
#endif  
  
  Coords = new double* [N_Coord+1];
  
  for(i=0; i<N_Coord+1; i++)  
   Coords[i] = new double [MaxN_QuadPoints_1D];
    
  for(i=0; i<MaxN_QuadPoints_1D; i++)
   {
    Coords[0][i] = X;
    Coords[1][i] = Y;
#ifdef __3D__  
    Coords[2][i] = Z;
#endif

    aux[i] = new double[6];
    coeff[i] = new double[6];
   }
   
#ifdef __3D__  
   Intl_L = Coords[3];
#else
   Intl_L = Coords[2];
#endif

  // associate each cell with her number in the collection
  for(i=0;i<N_Cells_Internal;i++)
  {
    Cell = Coll->GetCell(i);
    Cell->SetClipBoard(i);
  }

  for(i=0; i<N_Cells_Internal; i++)
   {
    UpdateEdge = FALSE;
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    bf = Element->GetBaseFunct1D();
    N_BaseFunct = Element->GetN_DOF();
    BaseFunct_ID = Element->GetBaseFunct1D_ID();
    h1 = Cell->GetMeasure();

//     cout << " BaseFunct_ID " << BaseFunct_ID << endl;
    l = bf->GetPolynomialDegree();
    LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(20*l);
    qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
    qf1->GetFormulaData(N_Points, Weights, zeta);

    F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)F_K)->SetCell(Cell);
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);
    rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

//     Bilinear(N_Points, Space_X, Space_Y, Intl_L, aux, coeff);
    Bilinear(N_Points, N_Coord, Coords, aux, coeff); 
    
    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    origvaluesD1=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D1);

    memset(LocMatrixA11, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA12, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA21, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocMatrixA22, 0, N_BaseFunct*N_BaseFunct*SizeOfDouble);
    memset(LocRhs, 0, N_BaseFunct*SizeOfDouble);

    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Coeff = coeff[j];
      c0 = Coeff[0]; // diffusion
      g0 = Coeff[1]; // convection in z direction
      c1 = Coeff[2]; // reaction term
      rhsval = Coeff[3]; //rhs

      // in case of PBE, it is the growth term from concentration
      // in PBE example file define Coeff[1] as < 0, so Conv will be the growth term
      if(g0>=0)
        Conv = g0;

//       cout<< " diffusion  " << c0 << " g " <<g0 << " reaction " << c1 <<endl;   
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      orgD1 = origvaluesD1[j];

      len +=Mult;
      for(k=0;k<N_BaseFunct;k++)
       {
        test0  = orgD0[k];
        test1  = orgD1[k];

        LocRhs[k] += Mult*rhsval*test0;

        for(l=0;l<N_BaseFunct;l++)
         {
          ansatz0  = orgD0[l];
          ansatz1  = orgD1[l];

          val  = Conv*ansatz1*test0;
          val += c0*ansatz1*test1; // difusive term
          val += c1*ansatz0*test0;
//            if(i==0)
//            cout<< "A ( " << k <<" , " <<l << " ) " <<Mult*val <<endl;
          LocMatrixA11[k*N_BaseFunct + l] += (Mult*val);
        }
       }
     }

     N_Joints=Cell->GetN_Edges();
     for(m=0; m<N_Joints; m++)
      {
       Neigh = (Cell->GetJoint(m))->GetNeighbour(Cell);

       if(m==0 && Neigh) // inner joint, since in 1D we have only 2 joints (really vertices)
       {
        UpdateEdge = TRUE;
        // only first joint (really vertices) will be updated, 
        // other joint will be the first joint of the next cell
        h2 = Neigh->GetMeasure();
        h_max = MAX(h1, h2);

        Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
//         Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
        Bilinear(1, N_Coord, Coords, aux, coeff); 
    
        c0 = coeff[0][0]; // diffusion

        //find the current cell basis value at this joint (X_n+)
        bf->GetDerivatives(D0, -1., JointValues);
        bf->GetDerivatives(D1, -1., JointValuesX);

        for(k=0;k<N_BaseFunct;k++)
         JointValuesX[k] *=rec_detjk;

        //find the neib cell basis value at this joint (X_n-)
        Neigh_i=Neigh->GetClipBoard();
        Neigh_FEId = FESpace1D_Intl->GetFE1D(Neigh_i, Neigh);
        Neigh_Element = TFEDatabase2D::GetFE1D(Neigh_FEId);
        Neigh_bf = Neigh_Element->GetBaseFunct1D();
        Neigh_N_BaseFunct = Neigh_Element->GetN_DOF();

        NeibDOF = GlobalNumbers + BeginIndex[Neigh_i];

        F_K = TFEDatabase2D::GetRefTrans1D(LineAffin);
       ((TLineAffin *)F_K)->SetCell(Neigh);
        Neigh_rec_detjk = ((TLineAffin *)F_K)->Getrec_detjk();

        Neigh_bf->GetDerivatives(D0, 1., Neigh_JointValues);
        Neigh_bf->GetDerivatives(D1, 1., Neigh_JointValuesX);

        for(k=0;k<Neigh_N_BaseFunct;k++)
         Neigh_JointValuesX[k] *=Neigh_rec_detjk;

        //(X_n+, X_n+) (test, ansatz)
        for(k=0;k<N_BaseFunct;k++)
         for(l=0;l<N_BaseFunct;l++)
         {
          val = 0.5*c0*JointValuesX[l]*JointValues[k] - 0.5*c0*Epsilon*JointValues[l]*JointValuesX[k]
                + (Sigma0/h_max)*JointValues[l]*JointValues[k] + (Sigma1/h_max)*JointValuesX[l]*JointValuesX[k];

         LocMatrixA11[k*N_BaseFunct + l] += val;
//          cout<<  " val  " << h_max*val  <<endl;
         }

        // e_n
        for(k=0;k<N_BaseFunct;k++) // own test X_n+
         for(l=0;l<Neigh_N_BaseFunct;l++) // neib ansatz X_n-
          {
           val = 0.5*c0*Neigh_JointValuesX[l]*JointValues[k] + 0.5*c0*Epsilon*Neigh_JointValues[l]*JointValuesX[k]
                 - (Sigma0/h_max)*Neigh_JointValues[l]*JointValues[k] - (Sigma1/h_max)*Neigh_JointValuesX[l]*JointValuesX[k];
 
           LocMatrixA21[k*N_BaseFunct + l] += val;
//            cout<<  " val  " << h_max*val  <<endl;
          }

        // d_n
        for(k=0;k<Neigh_N_BaseFunct;k++) // neib test X_n-
         for(l=0;l<N_BaseFunct;l++) // own ansatz X_n+
         {
           val = -0.5*c0*JointValuesX[l]*Neigh_JointValues[k] - 0.5*c0*Epsilon*JointValues[l]*Neigh_JointValuesX[k]
                 - (Sigma0/h_max)*JointValues[l]*Neigh_JointValues[k] - (Sigma1/h_max)*JointValuesX[l]*Neigh_JointValuesX[k];
 
           LocMatrixA12[k*N_BaseFunct + l] += val;
//            cout<<  " val  " << h_max*val  <<endl;
          }

        //(X_n-, X_n-) (test, ansatz)
        for(k=0;k<Neigh_N_BaseFunct;k++)
         for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          val = -0.5*c0*Neigh_JointValuesX[l]*Neigh_JointValues[k] + 0.5*c0*Epsilon*Neigh_JointValues[l]*Neigh_JointValuesX[k]
                + (Sigma0/h_max)*Neigh_JointValues[l]*Neigh_JointValues[k] + (Sigma1/h_max)*Neigh_JointValuesX[l]*Neigh_JointValuesX[k];

          LocMatrixA22[k*N_BaseFunct  + l] += val;
//          cout<<  " val  " << h_max*val  <<endl;
         }
       }//if(neigh)
       else if(!Neigh) // boundary joint
       {
        if(m==0) // L_min 
         {

         if(cond_Lmin==DIRICHLET)
          {
            //find the current cell basis value at this joint (X_n+)
           bf->GetDerivatives(D0, -1., JointValues);
           bf->GetDerivatives(D1, -1., JointValuesX);

           Intl_L[0] = TDatabase::ParamDB->REACTOR_P12;
//            Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
           Bilinear(1, N_Coord, Coords, aux, coeff);
   
           c0 = coeff[0][0]; // diffusion

//            BDValue_LMin(Intl_L[0], 0, BdValues);

           BdValues[0]=Bnuc;
 
           for(k=0;k<N_BaseFunct;k++)
            JointValuesX[k] *=rec_detjk;

           for(k=0;k<N_BaseFunct;k++)
            {
             val =  -Epsilon*c0*JointValuesX[k]*BdValues[0];
             val +=  (Sigma0/h1)*JointValues[k]*BdValues[0];

             LocRhs[k] += val;
//                cout<<  " val    rhsval " << val  <<endl;
             for(l=0;l<N_BaseFunct;l++)
              {
               val = c0*JointValuesX[l]*JointValues[k] - c0*Epsilon*JointValues[l]*JointValuesX[k]
                    + (Sigma0/h1)*JointValues[l]*JointValues[k] + (Sigma1/h1)*JointValuesX[l]*JointValuesX[k];

               LocMatrixA11[k*N_BaseFunct + l] += val;
//                 cout<<  " val L_min  " << val  <<endl;
              }
             }
          } //if(cond_Lmin==DIRICHLET)
         } // if(m==0) 
        else // L_Max
	 {
          if(cond_Lmax==DIRICHLET || cond_Lmax==NEUMANN)
           {
             //find the current cell basis value at this joint (X_n+)
            bf->GetDerivatives(D0, 1., JointValues);
            bf->GetDerivatives(D1, 1., JointValuesX);

            Intl_L[0] = TDatabase::ParamDB->REACTOR_P13;
//             Bilinear(1, Space_X, Space_Y, Intl_L, aux, coeff);
            Bilinear(1, N_Coord, Coords, aux, coeff);
   
            c0 = coeff[0][0]; // diffusion

            BDVal_LMax(Intl_L[0], 0, BdValues);

            for(k=0;k<N_BaseFunct;k++)
             JointValuesX[k] *=rec_detjk;

            for(k=0;k<N_BaseFunct;k++)
             {
              val =  Epsilon*c0*JointValuesX[k]*BdValues[0];
              val +=  (Sigma1/h1)*JointValuesX[k]*BdValues[1];

              LocRhs[k] += val;
//                 cout<<  " val    rhsval " << val  <<endl;
             for(l=0;l<N_BaseFunct;l++)
              {
               val = -c0*JointValuesX[l]*JointValues[k] + c0*Epsilon*JointValues[l]*JointValuesX[k]
                     + (Sigma0/h1)*JointValues[l]*JointValues[k] + (Sigma1/h1)*JointValuesX[l]*JointValuesX[k];

               LocMatrixA11[k*N_BaseFunct + l] += val;
//                cout<<  " val L_max " << val  <<endl;
              }

            }
          }
         }
       }

      }


    // add to global matrices
    for(j=0;j<N_BaseFunct;j++)
     {
      TestDOF = DOF[j];
      rhs[TestDOF] += LocRhs[j]; 

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA11[j*N_BaseFunct + l];
//             cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } // for(m=0;m<N_BaseFunct

       // add edge integrals
      if(UpdateEdge)
       {
        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
         if(KCol[k] == NeibDOF[l])
          {
//             cout << TestDOF << " " << KCol[k] << endl;
           ValuesA[k] +=LocMatrixA21[j*N_BaseFunct + l]; 
           break;
          } 
        } // for(l=0;l<Neigh_N_BaseFunct;l++)
       }

      } // for(n=begin;n<end;n++)
     } //  for(j=0;j<N_BaseFunct;j++)

   // add edge integrals
   if(UpdateEdge)
   {
    for(j=0;j<Neigh_N_BaseFunct;j++)
     {
      TestDOF = NeibDOF[j];

      begin = RowPtr[TestDOF];
      end = RowPtr[TestDOF+1];
      for(k=begin;k<end;k++)
       {
        for(l=0;l<N_BaseFunct;l++)
         {
          if(KCol[k] == DOF[l])
           {
            ValuesA[k] +=LocMatrixA12[j*N_BaseFunct + l];
//             cout << TestDOF << " " << KCol[k] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)

        for(l=0;l<Neigh_N_BaseFunct;l++)
         {
          if(KCol[k] == NeibDOF[l])
           {
            ValuesA[k] +=LocMatrixA22[j*N_BaseFunct + l];
//             cout << TestDOF << " " << KCol[k] <<  " " <<  LocMatrixA22[j*N_BaseFunct + l] << endl;
            break;
           }
          } //  for(l=0;l<N_BaseFunct;l++)
        } // for(k=begin;k<end;k++)
     } // for(j=0;j<Neigh_N_BaseFunct;j++)
    } // if(UpdateEdge)
    

   }// for(i=0; i<N_Cells_Internal; i++)

// cout<< " len " << len << endl;
//  exit(0);
// // print matrix
//   for(j=0;j<N_V;j++)
//    {
//     begin = RowPtr[j];
//     end = RowPtr[j+1];
//     for(k=begin;k<end;k++)
//      {
//       cout << "A(" << j << ", "<< KCol[k] << ") = " << ValuesA[k] <<endl;
//      }
// 
//       cout << "f: " << rhs[j ] <<endl;
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
} // void TADISystem1D::AssembleARhs(d


double TADISystem1D::GetQ3Max(double *currsol)
{
 int i, j, k, l, N_Cells_Internal, N_BaseFunct;
 int N_Points, N_Sets=1, *GlobalNumbers, *BeginIndex, *DOF;
 int TestDOF, begin, end, *RowPtr, *KCol;

 double g=0., L;
 double *Weights, *zeta, Intl_L[MaxN_QuadPoints_1D], AbsDetjk[MaxN_QuadPoints_1D];
 double **origvaluesD0, Mult;
 double *orgD0, test0;
 bool Needs2ndDer[1];

 TCollection *Coll;
 TBaseCell *Cell;
 FE1D FEId;
 TFE1D *Element;
 TBaseFunct1D *bf;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TRefTrans1D *F_K;
 BaseFunct1D BaseFunct_ID, BaseFunct[1];

  Coll = FESpace1D_Intl->GetCollection();
  GlobalNumbers = FESpace1D_Intl->GetGlobalNumbers();
  BeginIndex = FESpace1D_Intl->GetBeginIndex();

  N_Cells_Internal = Coll->GetN_Cells();

//      for(k=0;k<N_V;k++)
//        {
//         cout <<   " currsol " <<  setprecision(8) << currsol[k] << endl;
//        }

  for(i=0; i<N_Cells_Internal; i++)
   {
    Cell = Coll->GetCell(i);
    FEId = FESpace1D_Intl->GetFE1D(i, Cell);
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
    ((TLineAffin *)F_K)->GetOrigFromRef(N_Points, zeta, Intl_L, AbsDetjk);

    BaseFunct[0] = BaseFunct_ID;
    ((TLineAffin *)F_K)->GetOrigValues(N_Sets, BaseFunct, N_Points, zeta,  LineQuadFormula,  Needs2ndDer);

    origvaluesD0=TFEDatabase2D::GetOrigElementValues(BaseFunct_ID, D0);
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Points;j++)
     {
      Mult = Weights[j]*AbsDetjk[j];
      orgD0 = origvaluesD0[j];
      L = Intl_L[j];

//       len +=Mult;

      for(k=0;k<N_BaseFunct;k++)
       {
        if(currsol[DOF[k]]>0.)
	 { g += L*L*L*Mult*currsol[DOF[k]]*orgD0[k]; }
//         cout << " L " << L << " currsol " <<  setprecision(8) << currsol[DOF[k]] << endl;
       }
//        if( currsol[0]>0.9)
//        cout << " L*L " << L*L << " g " << g << " sol[DOF[k]] " << sol[0] << endl;
     } //  for(j=0;j<N_Points;j++)
   }

// exit(0);
 return g;
}


TADISystem1D::~TADISystem1D()
{
  if(rhs) delete [] rhs;
}
