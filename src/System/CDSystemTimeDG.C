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
   
/** ************************************************************************ 
*
* @brief     stores the information of a 2D scalar dG in time discretization
* @author    Sashikumaar Ganesan, 
* @date      24.12.15
* @History    
 ************************************************************************  */


#include <string.h>
#include <CDSystemTimeDG.h>
#include <FEDatabase2D.h>
#include <Database.h>
#include <DirectSolver.h>
#include <stdlib.h>

TCDSystemTimeDG::TCDSystemTimeDG(int n_RowBlockMat, int n_ColBlockMat, TSquareMatrix2D *mat_M, TSquareMatrix2D *mat_A)
{
 int *RowPtr, N_Blocks;
  
  N_RowBlockMat = n_RowBlockMat;
  N_ColBlockMat = n_ColBlockMat;
  N_Blocks = N_RowBlockMat*N_ColBlockMat;
  Mat_M = mat_M;
  Mat_A = mat_A;

  if (Mat_M->GetColOrder() != 1)
   {
    OutPut("dGTime: sort column numbers in Mass and Stiffness mat using TSquareStructure"<<endl);
    OutPut("dGTime: needs it to call umfpack"<<endl);
    exit(0);
   }
  
  N_U = Mat_M->GetN_Columns();  
  N_UActive = Mat_M->GetActiveBound();
  RowPtr = Mat_M->GetRowPtr(); 
  RhsSpace = Mat_M->GetFESpace(); // test and ansatz is same
  
  N_Eqn = N_ColBlockMat*N_U;
  Sys_N_Entries = N_Blocks*RowPtr[N_U];
  
  Sys_rowptr = new int[(N_RowBlockMat*N_U)+1]; 
  Sys_colindex = new int[Sys_N_Entries]; 
  
  Sol_dG = new double [N_Eqn];
  rhs_dG = new double [N_Eqn];  
  memset(Sol_dG, 0, N_Eqn*SizeOfDouble);
  memset(rhs_dG, 0, N_Eqn*SizeOfDouble);

  Sys_structure = new TSquareStructure2D(N_Eqn, Sys_N_Entries, Sys_colindex, Sys_rowptr);
  Sys_Mat = new TSquareMatrix2D(Sys_structure);
  
  //needed only in ALE approach
  CONSERVATIVEALE = FALSE;
  QpMatricsAdded = FALSE;
//      cout <<  " TCDSystemTimeDG  "  << Sys_N_Entries << endl;
//           cout <<  " Sys_Mat  "  << Sys_Mat->GetActiveBound() << endl;
//      exit(0);
}

TCDSystemTimeDG::~TCDSystemTimeDG()
{
  delete [] Sys_rowptr;
  delete [] Sys_colindex; 
  delete [] Sol_dG;
  delete [] rhs_dG;  
}


void TCDSystemTimeDG::AssembleRhs(int N_Rhs, double *T, double *Rhs)
{
  int i, j, k, l, m, N_Cells, N_LocalUsedElements, N_Points;
  int *N_BaseFunct,N_U_LocalDof, *DOF, *BeginIndex, *GlobalNumbers;
  
  double *weights, *xi, *eta, orig_time;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], AbsDetjk[MaxN_QuadPoints_2D];  
  double **uorig, *Orig, Mult, LocRhs[MaxN_BaseFunctions2D];
  double *coeffs[MaxN_QuadPoints_2D], *aux, *rhs_aux;
    
  TCollection *Coll;
  TBaseCell *Me;
  FE2D FEId;
  TFE2D *ele;
  FE2D LocalUsedElements[1];
  BaseFunct2D BaseFunct, *BaseFuncts;
 
  bool *SecondDer;
 
  SecondDer = new bool[2];
  SecondDer[0] = FALSE;
  SecondDer[1] = FALSE;
 
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
 
  memset(Rhs, 0, N_Rhs*N_U*SizeOfDouble);
  
  Coll = RhsSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = RhsSpace->GetBeginIndex();
  GlobalNumbers = RhsSpace->GetGlobalNumbers();
    
  orig_time = TDatabase::TimeDB->CURRENTTIME;
  
  aux = new double [MaxN_QuadPoints_2D*10];
  for(i=0;i<MaxN_QuadPoints_2D;i++)
   coeffs [i] = aux + i*10;
  
  for(i=0;i<N_Cells;i++)
   {
    Me = Coll->GetCell(i);
    FEId = RhsSpace->GetFE2D(i, Me);   
    
    N_LocalUsedElements = 1;
    LocalUsedElements[0]=FEId;   
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                           Coll, Me, SecondDer, 
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);    
    
     BaseFunct = BaseFuncts[FEId];
     N_U_LocalDof = N_BaseFunct[FEId];
      
     uorig = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);   
    
     DOF = GlobalNumbers + BeginIndex[i];
     
     //compute rhs for the given time instances
     for(l=0;l<N_Rhs;l++)
      {
       TDatabase::TimeDB->CURRENTTIME = T[l];
        BilinearCoeffs(N_Points, X, Y, NULL, coeffs);
     
       memset(LocRhs, 0, N_U_LocalDof*SizeOfDouble);
       for(j=0;j<N_Points;j++)
        {
         Orig = uorig[j];
         Mult = weights[j]*AbsDetjk[j]*coeffs[j][4]; // mult * rhs

         // rhs assemble
         for(k=0;k<N_U_LocalDof;k++)
           LocRhs[k] += Mult*Orig[k];

        } // or(j=0;j<N_Poi

       //add local rhs to global rhs
       rhs_aux = Rhs+(l*N_U);  
       for(k=0;k<N_U_LocalDof;k++)
        {
          m = DOF[k];
          rhs_aux[m] += LocRhs[k];
        } // for(k=0;k<N_U_Local
       } //  for(l=0;l<N_Rhs;l
   } // for(i=0;i<N_Cells;i++)
   
  // set the original time
  TDatabase::TimeDB->CURRENTTIME = orig_time;
 
  delete [] aux;
  delete [] SecondDer;
  
//       cout <<  " TCDSystemTimeDG AssembleRhs"    << endl;
//       exit(0);
}


void TCDSystemTimeDG::AssembleSysMat(double *Mu_old, double *Rhs)
{
  
      cout <<  " TCDSystemTimeDG AssembleSysMat"    << endl;
  
      exit(0);
      
}

void TCDSystemTimeDG::AssembleALESysMat_Qp1(double *Mu_old, double *Rhs)
{
  
      cout <<  " TCDSystemTimeDG AssembleALESysMat_Qp1"    << endl;
  
      exit(0);
      
}



void TCDSystemTimeDG::SoveTimedG(double *Sol)
{
  
  cout <<  " TCDSystemTimeDG SoveTimedG: Call appropriate dG"    << endl;
  exit(0); 
}


void TCDSystemTimeDG::SovedGSystem()
{
//     double *Sys_Entries;
//   Sys_Entries = Sys_Mat->GetEntries();    
//        cout <<  " TCDSystemTimeDG  "  << Sys_Mat->GetColOrder() << endl; 
//       exit(0);
  
//   for(int i=Sys_rowptr[6377]; i< Sys_rowptr[6378]; i++)
//       cout << Sys_colindex[i] << " Sys_Entries: " << Sys_Entries[i]  << endl; 
//   
//   
//           cout << Sys_rowptr[6377] << " rhs_dG[i]: " << rhs_dG[6377]  << endl; 
// 	  
// 	  exit(0);
	  
    // call the direct solver
    DirectSolver(Sys_Mat, rhs_dG, Sol_dG);
//             cout << 6377 << " Sol_dG[i]: " << Sol_dG[6377]  << endl; 
// 	    exit(0);
}
