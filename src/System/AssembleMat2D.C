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
* @brief     source file for TAssembleMat2D
* @author    Sashikumaar Ganesan
* @date      16.04.16
* @History   
 ************************************************************************  */
#include <string.h>
#include <AssembleMat2D.h>
#include <SquareMatrix2D.h>
#include <Matrix2D.h>
#include <FESpace2D.h>
#include <AuxParam2D.h>
#include <DiscreteForm2D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <BaseCell.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <IsoBoundEdge.h>
#include <BoundComp2D.h>
#include <stdlib.h>


TAssembleMat2D::TAssembleMat2D(int n_fespaces, TFESpace2D **fespaces,
                         int n_sqmatrices, TSquareMatrix2D **sqmatrices,
                         int n_matrices, TMatrix2D **matrices,
                         int n_rhs, double **rhs, TFESpace2D **ferhs,
                         TDiscreteForm2D *discreteform,
                         BoundCondFunct2D **boundarybonditions,
                         BoundValueFunct2D **boundaryvalues,
                         TAuxParam2D *parameters)
                  :TAssemble(n_fespaces, n_sqmatrices, n_matrices, n_rhs, rhs)
{
 int i, j, N;
 
  TFESpace2D *fespace;

   DiscreteForm = discreteform;
   BoundaryConditions = boundarybonditions;
   BoundaryValues = boundaryvalues;
   AuxParam = parameters; 
   
   N_Parameters = AuxParam->GetN_Parameters();   
   SecondDer = discreteform->GetNeeds2ndDerivatives();
  
  if(n_fespaces)
   { FeSpaces = new TFESpace2D *[n_fespaces]; }
  if(n_rhs)
   { FeRhs = new TFESpace2D *[n_rhs]; }
  if(n_sqmatrices)
   { SqMatrices = new TSquareMatrix2D* [n_sqmatrices]; }
  if(n_matrices)
   { RecMatrices = new TMatrix2D *[n_matrices]; } 
    
    // all spaces use same Coll
    Coll = fespaces[0]->GetCollection();
    N_Cells = Coll->GetN_Cells();  
 
     for(i=0;i<n_fespaces;i++)   
       FeSpaces[i] = fespaces[i];    
    
     for(i=0;i<n_rhs;i++)
     { 
      Rhs[i] = rhs[i];     
      FeRhs[i] = ferhs[i];      
     }
          
    for(i=0;i<N_SqMatrices;i++)
     {
      SqMatrices[i] = sqmatrices[i];
      fespace = sqmatrices[i]->GetFESpace();
      GlobalNumbers[i] = fespace->GetGlobalNumbers();
      BeginIndex[i] = fespace->GetBeginIndex();
     } // endfor 
       
    for(i=0;i<N_Matrices;i++)
     {
      RecMatrices[i] = matrices[i];       
      fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetTestSpace();
       
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();
  
      fespace = (TFESpace2D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
     } // endfor
    
    for(i=0;i<N_Rhs;i++)
     {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
     } // endfor   
} // TAssembleMat2D

TAssembleMat2D::~TAssembleMat2D()
{

delete [] FeSpaces;
delete [] FeRhs;
delete [] SqMatrices; 
delete [] RecMatrices; 
   
}

void TAssembleMat2D::Init()
{
 int i, j, N;
 TFESpace2D *fespace;
   
   for(i=0;i<N_SqMatrices;i++)
     {      
      N = SqMatrices[i]->GetHangingN_Entries();
      if(N>0)
       {
        HangingEntries[i] = new double [N];
        memset(HangingEntries[i], 0, SizeOfDouble*N);
       }
     } // endfor 
     
    for(i=0;i<N_Rhs;i++)
     {
      fespace = FeRhs[i];      
      N = fespace->GetN_Hanging();
      if(N>0)
       {
        HangingRhs[i] = new double [N];
        memset(HangingRhs[i], 0, SizeOfDouble*N);
       }
     } // endfor
     
    if(N_Rhs)
     rhsaux = new double [N_Rhs*MaxN_BaseFunctions2D];
    
    for(i=0;i<N_Rhs;i++)
      LocRhs[i] = rhsaux+i*MaxN_BaseFunctions2D;    
  
    if(N_Parameters)
     {
      paramaux = new double [MaxN_QuadPoints_2D*N_Parameters];
      
      for(j=0;j<MaxN_QuadPoints_2D;j++)
       Param[j] = paramaux + j*N_Parameters;
     }    

  // 40 <= number of terms in bilinear form
  // DUE NOTE CHANGE BELOW 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
    auxarray = new double [MaxN_QuadPoints_2D*40]; 
    for(j=0;j<MaxN_QuadPoints_2D;j++)
     AuxArray[j] = auxarray + j*40;    
 
   if(N_AllMatrices)
    {
      auxmat = new double [N_AllMatrices*MaxN_BaseFunctions2D*MaxN_BaseFunctions2D];
      Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions2D];
      
      for(j=0;j<N_AllMatrices*MaxN_BaseFunctions2D;j++)
        Matrices[j] = auxmat+j*MaxN_BaseFunctions2D;

      for(i=0;i<N_AllMatrices;i++)
        LocMatrices[i] = Matrices+i*MaxN_BaseFunctions2D;
    } // endif N_AllMatrices
} // Init


void TAssembleMat2D::DeAllocate()
{
 int i, j, N;
 TFESpace2D *fespace;  
  
    for(i=0;i<N_SqMatrices;i++)
     {
      N = SqMatrices[i]->GetHangingN_Entries();    
      if(N>0)
        delete []  HangingEntries[i] ;
     } // endfor 
     
    for(i=0;i<N_Rhs;i++)
     {
      fespace = FeRhs[i];
      N = fespace->GetN_Hanging();
      if(N>0)
        delete []  HangingRhs[i] ;
     } // endfor
     
    if(N_Rhs)
     delete []  rhsaux;   
  
    if(N_Parameters)
      delete [] paramaux;  

    delete [] auxarray;    
 
   if(N_AllMatrices)
    {
      delete [] auxmat;
      delete [] Matrices;
    } // endif N_AllMatrices 
  
} // DeAllocate

void TAssembleMat2D::Reset()
{
  int i, N;
  
   // reset to zero
   for(i=0;i<N_SqMatrices;i++)
      SqMatrices[i]->Reset(); 
     
    for(i=0;i<N_Matrices;i++)
      RecMatrices[i]->Reset(); 
   
    for(i=0;i<N_Rhs;i++)
     {   
      N = FeRhs[i]->GetN_DegreesOfFreedom();
      memset(Rhs[i], 0, N*SizeOfDouble);
     }   
} // TAssembleMat2D::Reset()


void TAssembleMat2D::Assemble2D()
{
  int i, j, k, ij;
  int N_Edges, N_VertInCell, *N_BaseFunct, N_LocalUsedElements, N_Vertex;
  int LocN_BF[N_BaseFuncts2D], N_Points;
  
  double hK;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  double *weights, *xi, *eta;
  
  TBaseCell *cell, *neigh;
  TJoint *joint;
  TBoundFace *boundface;
  BaseFunct2D *BaseFuncts;
  FE2D LocalUsedElements[N_FEs2D], CurrentElement;
  BaseFunct2D LocBF[N_BaseFuncts2D];
  RefTrans2D reftrans;
   
  // reset clipboards of all neighbours to -1
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Joints();
    for(j=0;j<k;j++)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if(neigh) neigh->SetClipBoard(-1);
    }
    cell->SetClipBoard(-1);

#ifdef _MPI       
      N_Edges=cell->GetN_Edges();
      for(j=0;j<N_Edges;j++)
       (cell->GetEdge(j))->SetClipBoard(-1);
      
      N_VertInCell = cell->GetN_Vertices();
      for(j=0;j<N_VertInCell;j++)
       (cell->GetVertex(j))->SetClipBoard(-1); 
#endif     
  } // endfor i

  for(i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);  
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  
//========================================================================
// loop over all cells
//========================================================================
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // cout << "Cell No.: " << i << endl;
    TDatabase::ParamDB->INTERNAL_CELL = i;

    // only for multiphase flows
//     if(AssemblePhaseID >= 0)
//      {
//       if(AssemblePhaseID != cell->GetRegionID() )
//        continue;
//      }
     
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
    
    //========================================================================
    // find local used elements on this cell
    //========================================================================
    for(j=0;j<N_FeSpaces;j++)
    {
      CurrentElement = FeSpaces[j]->GetFE2D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = N_FeSpaces;
    //========================================================================
    // calculate values on original element
    //========================================================================

    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                           Coll, cell, SecondDer,
                           N_Points, xi, eta, weights, X, Y, AbsDetjk);
    
//      cout << "AbsDetjk: " << AbsDetjk[0] << endl;
//     OutPut("params " << TDatabase::ParamDB->INTERNAL_LEVEL << endl);
    AuxParam->GetParameters(N_Points, Coll, cell, i, xi, eta, X, Y, Param); 
    
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

// #ifdef __2D__
//     if(Aux2D2D)
//       Aux2D2D->GetGradient(i, N_Points, xi, eta, Param);
// #endif

    // use DiscreteForm to assemble a few matrices and
    // right-hand sides at once 

//     OutPut("local form " << i << endl);
    if(DiscreteForm)
     {
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk,
                                  hK, X, Y,
                                  LocN_BF, LocBF,
                                  Param, AuxArray,
                                  cell,
                                  N_AllMatrices, N_Rhs,
                                  LocMatrices, LocRhs);
     }

    //========================================================================    
    // add local matrices to global matrices (ansatz == test)       
    //========================================================================   
    if(N_SqMatrices)
     this->AddLocalSqMatToGlobal(i, cell, N_BaseFunct);

    //========================================================================       
    // add local matrices to global matrices (ansatz != test)
    //========================================================================      
    if(N_Matrices)
      this->AddLocalRecMatToGlobal(i, cell, N_BaseFunct);   
        
    //========================================================================       
    // add local right-hand sides to global right-hand side
    //========================================================================        
    if(N_Rhs)
      this->AddLocalRhsToGlobal(i, cell, N_BaseFunct, BaseFuncts, reftrans);   
  } // for(i=0;i<N_Cells;i++)
  
   //========================================================================     
   // modify matrix according to coupling 
   //========================================================================   
   this->ModifyMatHang();
  
   //========================================================================     
   // modify matrix according to coupling 
   //========================================================================   
//    this->PrintAllMat();   
   
//   cout << "TAssembleMat2D Assemble2D() done ! " << endl;
//   exit(0);
} // Assemble2D

void TAssembleMat2D::PrintAllMat()
{
 int i, j, k, end, *RowPtr, *ColInd, N_Rows, ActiveBound;
 
 double *Entries, *RHS;
 
 TFESpace2D *fespace;  
  
  for(k=0;k<N_SqMatrices;k++)
  {
    OutPut(endl);
    cout << "sqmatrix: " << k << endl;
    RowPtr = SqMatrices[k]->GetRowPtr();
    Entries = SqMatrices[k]->GetEntries();
    ColInd = SqMatrices[k]->GetKCol();
    N_Rows = SqMatrices[k]->GetN_Rows();
    
    fespace = SqMatrices[k]->GetFESpace();
    ActiveBound = fespace->GetActiveBound();
    for(i=0;i<ActiveBound;i++)
    {
     end=RowPtr[i+1];
      
     OutPut("\033[1;31ma"<<k<<"("<< i+1 << ", " << "j):"<<"\033[0m");
       
     for(j=RowPtr[i];j<end;j++)
      OutPut("\033[1;31m" << ColInd[j]+1 <<"=\033[0m"<<Entries[j] << ", ");
     
    OutPut(" "<< endl);
    }
    cout << endl;
  } // endfor k
  
  for(k=0;k<N_Matrices;k++)
  {
    cout << endl;
    cout << "matrix: " << k << endl;
    RowPtr = RecMatrices[k]->GetRowPtr();
    Entries = RecMatrices[k]->GetEntries();
    ColInd = RecMatrices[k]->GetKCol();
    N_Rows = RecMatrices[k]->GetN_Rows();
    for(i=0;i<N_Rows;i++)
    {
      end=RowPtr[i+1];
      for(j=RowPtr[i];j<end;j++)
      {
        // cout << j << endl;
         OutPut("b"<<k<<"("<< i+1 << "," << ColInd[j]+1 << " )=  ");
         OutPut(Entries[j] << ";" << endl);
              }
    }
    cout << endl;
  } // endfor k

  for(k=0;k<N_Rhs;k++)
  {
    cout << "rhs: " << k << endl;
    N_Rows = FeRhs[k]->GetN_DegreesOfFreedom();
    RHS=Rhs[k];
    for(i=0;i<N_Rows;i++)
      cout << setw(5) << i << setw(20) << RHS[i] << endl;
  }
  
  
} //TAssembleMat2D::PrintAllMat()


void TAssembleMat2D::ModifyMatHang()
{
 int i, j, k, l, m, n, l1, l2, last, N_, end, *DOF, N_Hanging, N_Rows;
 int *HangingRowPtr, *HangingColInd, *RowPtr, *ColInd, ActiveBound, AnsatzActiveBound, AnsatzHangingBound;
  
 double *Coupling, v, *Entries, *RHS;
 double *CurrentHangingEntries, *CurrentHangingRhs;

  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  TFESpace2D *fespace;  

  for(j=0;j<N_SqMatrices;j++)
   {
    // find space for this bilinear form
    fespace = SqMatrices[j]->GetFESpace();  
    N_ = fespace->GetN_Hanging();
    // there are no hanging nodes
    if (N_ == 0)
      continue;
    
    HangingNodes = fespace->GetHangingNodes();  
    
    Entries = SqMatrices[j]->GetEntries();
    RowPtr = SqMatrices[j]->GetRowPtr();
    ColInd = SqMatrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j];
    HangingRowPtr = SqMatrices[j]->GetHangingRowPtr();
    HangingColInd = SqMatrices[j]->GetHangingKCol();

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
   } // for(j=0;j<N_SqMatrices;

 for(j=0;j<N_Matrices;j++)
  {
    // hanging nodes in test space
    fespace = (TFESpace2D *) (RecMatrices[j]->GetStructure()->GetTestSpace());
    N_ = fespace->GetN_Hanging();
    if (N_ == 0)
      continue; // there are no hanging nodes   
        
    HangingNodes = fespace->GetHangingNodes();

    Entries = RecMatrices[j]->GetEntries();
    RowPtr = RecMatrices[j]->GetRowPtr();
    ColInd = RecMatrices[j]->GetKCol();

    CurrentHangingEntries = HangingEntries[j+N_SqMatrices];
    HangingRowPtr = RecMatrices[j]->GetHangingRowPtr();
    HangingColInd = RecMatrices[j]->GetHangingKCol();

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
    N_Rows =  RecMatrices[j]->GetN_Rows();
    fespace = (TFESpace2D *) (RecMatrices[j]->GetStructure()->GetAnsatzSpace());

    N_ = fespace->GetN_Hanging();
    if (N_ == 0)
      continue; // there are no hanging nodes   
            
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
   
 for(j=0;j<N_Rhs;j++)
  {
    fespace = FeRhs[j];
    N_Hanging = fespace->GetN_Hanging();    

    if (N_Hanging == 0)
      continue; // there are no hanging nodes
       
    HangingNodes = fespace->GetHangingNodes();
    RHS = Rhs[j];
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
  for(j=0;j<N_SqMatrices;j++)
  {
    fespace = SqMatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    if (N_ == 0)
      continue; // there are no hanging nodes   
    
    HangingNodes = fespace->GetHangingNodes();

    Entries = SqMatrices[j]->GetEntries();
    RowPtr = SqMatrices[j]->GetRowPtr();
    ColInd = SqMatrices[j]->GetKCol();
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
  }    //   for(j=0;j<N_SqMatrices;j++)  
} // TAssembleMat2D::ModifyMatHang()


void TAssembleMat2D::AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct)
{
 int j, k, m, n, l, N_, l2, end, *RowPtr, *ColInd, *HangingRowPtr, *HangingColInd;
 int *DOF, ActiveBound, DirichletBound, H_N=0;
 
 double **Matrix, *MatrixRow, *CurrentHangingEntries, *Entries;
 
 TFESpace2D *fespace;  
 FE2D CurrentElement;  
 
    for(j=0;j<N_SqMatrices;j++)
    {
      // find space for this bilinear form
      fespace = SqMatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE2D(i, cell);
      N_ = N_BaseFunct[CurrentElement];

      Matrix = Matrices+j*MaxN_BaseFunctions2D;
      Entries = SqMatrices[j]->GetEntries();
      RowPtr = SqMatrices[j]->GetRowPtr();
      ColInd = SqMatrices[j]->GetKCol();

      H_N = SqMatrices[j]->GetHangingN_Entries();
      if(H_N>0)
      {
       CurrentHangingEntries = HangingEntries[j];
       HangingRowPtr = SqMatrices[j]->GetHangingRowPtr();
       HangingColInd = SqMatrices[j]->GetHangingKCol();
      }
      
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = GlobalNumbers[j] + BeginIndex[j][i];

      // add local matrix to global
      for(m=0;m<N_;m++)
      {
        l=DOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          // for all dof
          for(k=0;k<N_;k++)
          {
            // for all columns
            l2 = 0;
            end=RowPtr[l+1];
            for(n=RowPtr[l];n<end;n++)
            {
              // if column of current dof found
              if(DOF[k] == ColInd[n])
              {
                // cout << m << "   " << k << endl << n << endl;
                Entries[n] += MatrixRow[k];
                l2 = 1;
                break;
              } // endif
            } // endfor n
            if(l2 == 0)
	    {
              cout << l << "(assemble3d.c) not found :: " << RowPtr[l]<< " j :: " << j << " end :: " <<
              end << endl;
	      exit(0);
	    }
          } // endfor k
        } // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
	    if(H_N==0) continue;
	    
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
                } // endif
              } // endfor k
            } // endfor n 
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
      } // endfor m  
    } //for(j=0;j<N_SqMatric

//     cout << N_SqMatrices << " TAssembleMat2D AddLocalSqMatToGlobal() done ! " << endl;
//   exit(0);
} //TAssembleMat2D::AddLocalSqMatToGlobal()

void TAssembleMat2D::AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct)
{
 int j, k, l, m, n, end, N_Test, N_Ansatz, *RowPtr, *ColInd, *TestDOF, *AnsatzDOF;

 double **Matrix, *MatrixRow, *Entries;
 
 FE2D TestElement, AnsatzElement;
  
    for(j=0;j<N_Matrices;j++)
    {
      TestElement = ((TFESpace2D *) RecMatrices[j]->GetStructure()->GetTestSpace())->GetFE2D(i, cell);
      AnsatzElement = ((TFESpace2D *) RecMatrices[j]->GetStructure()->GetAnsatzSpace())->GetFE2D(i, cell);

      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+N_SqMatrices)*MaxN_BaseFunctions2D;

      Entries = RecMatrices[j]->GetEntries();
      RowPtr = RecMatrices[j]->GetRowPtr();
      ColInd = RecMatrices[j]->GetKCol();

      TestDOF = TestGlobalNumbers[j] + TestBeginIndex[j][i];
      AnsatzDOF = AnsatzGlobalNumbers[j] + AnsatzBeginIndex[j][i];

      // add local matrix to global
      for(m=0;m<N_Test;m++)
      {
        l=TestDOF[m];
        MatrixRow = Matrix[m];
        // cout << "DOF: " << l << endl;
        for(k=0;k<N_Ansatz;k++)
        {
          end = RowPtr[l+1];
          for(n=RowPtr[l];n<end;n++)
          {
            if(AnsatzDOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += MatrixRow[k];
              break;
            } // endif
          } // endfor n
        } // endfor k
      } // endfor m   
    } // 
  
} //   TAssembleMat2D::AddLocalMatToGlobal


void TAssembleMat2D::AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct2D *BaseFuncts, RefTrans2D reftrans)
{
 int j, k, l, m, l3, N_, jj, l1, dof, ActiveBound, DirichletBound, *DOF, N_H;
 int comp, N_EdgePoints, N_Joints, *JointDOF, N_Edges, N_EdgeDOF, *EdgeDOF;
 const int *TmpFV, *TmpLen, *EdgeVertex, *TmpFE, *ETmpLen;
 int MaxLen, EMaxLen, N_FaceEdges, N_LinePoints, Bd_ID;
  
 double x0, x1, y0, y1, hE, nx, ny, tx, ty, x, y, val, eps=1e-4;
#ifdef __3D__
  double z0, z1;
#endif
 double t0, t1, t, s, *EdgePoints; 
 double *local_rhs, *RHS, *CurrentHangingRhs;
 double *LineWeights, *zeta;
 double PointValues[MaxN_PointsForNodal2D];
 double FunctionalValues[MaxN_BaseFunctions2D];
 double **JointValues, *JointValue;
 
 TFESpace2D *fespace;  
 FE2D CurrentElement;  
 BoundCond Cond0, Cond1;
 BoundCondFunct2D *BoundaryCondition;
 BoundValueFunct2D *BoundaryValue; 
 TNodalFunctional2D *nf;
 TFE2D *ele;
 TFEDesc2D *FEDesc_Obj;
 TJoint *joint;
 TBoundEdge *boundedge;
 TBoundComp *BoundComp;
 TIsoBoundEdge *isoboundedge;
 QuadFormula1D LineQuadFormula;
 TQuadFormula1D *qf1;
 TVertex *Vert; 
     
   for(j=0;j<N_Rhs;j++)
    {
      local_rhs = rhsaux+j*MaxN_BaseFunctions2D;
      RHS = Rhs[j];
      
      fespace = FeRhs[j];
      CurrentElement = fespace->GetFE2D(i, cell);
      N_ = N_BaseFunct[CurrentElement];
   
      N_H = fespace->GetN_Hanging();
      if(N_H>0)
      { CurrentHangingRhs = HangingRhs[j]; }
      
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
      
      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      
      // add local right-hand side to the global one
      for(m=0;m<N_;m++)
      {
	if (TDatabase::ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS)
	  {
	    l = 0;
	    N_Edges=cell->GetN_Edges();
	    for(jj=0;jj<N_Edges;jj++)                        // loop over all edges of cell
	      {
		Vert = cell->GetVertex(jj);
		t0 =  Vert->GetClipBoard();
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
            if(N_H==0) continue;
	    
            // hanging node
            l -= ActiveBound;
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      }                                           // endfor m
   
      ele = TFEDatabase2D::GetFE2D(CurrentElement);
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
                l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(CurrentElement);
                // get a suitable line quadrature formula
                LineQuadFormula = TFEDatabase2D::GetQFLineFromDegree(2*l);
                qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
                qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(BaseFuncts[CurrentElement], LineQuadFormula, m);
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])->ChangeBF(Coll, cell, N_LinePoints, JointValues);
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
                TFEDatabase2D::GetBaseFunct2D(BaseFuncts[CurrentElement])->ChangeBF(Coll, cell, N_LinePoints, JointValues);
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
                TFEDatabase2D::GetBaseFunct2DFromFE2D(CurrentElement)->MakeRefElementData(LineQuadFormula);
                JointValues=TFEDatabase2D::GetJointValues2D(BaseFuncts[CurrentElement], LineQuadFormula, m);
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

    } //for(j=0;j<N_Rhs;j++)
  
/*   cout << "TAssembleMat2D AddLocalRhsToGlobal() done ! " << endl;
  exit(0);*/ 
} //TAssembleMat2D::AddLocalRhsToGlobal


// void TAssembleMat2D::AssembleNavierSlip()
// {
//   int i, j, k, ij, N_, N_H, ActiveBound, DirichletBound, *DOF;
//   int N_Edges, N_VertInCell, *N_BaseFunct, N_LocalUsedElements, N_Vertex;
//   int LocN_BF[N_BaseFuncts2D], N_Points, N_JointDOF;
//   
//   double hK, *local_rhs, *RHS, *CurrentHangingRhs;
//   double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];
//   double Z[MaxN_QuadPoints_2D];
//   double AbsDetjk[MaxN_QuadPoints_2D];
//   double *weights, *xi, *eta, *zeta;
//   
//   TFESpace2D *fespace;   
//   TBaseCell *cell, *neigh;
//   TJoint *joint;
//   TBoundFace *boundface;
//   BaseFunct2D *BaseFuncts;
//   FE2D LocalUsedElements[N_FEs2D], CurrentElement;
//   BaseFunct2D LocBF[N_BaseFuncts2D];
//   RefTrans2D reftrans;
//   BoundCond Cond0, Cond1;
//   BoundCondFunct2D *BoundaryCondition;
//   BoundValueFunct2D *BoundaryValue; 
//   TNodalFunctional2D *nf;
//   TFE2D *ele;
//   TFEDesc2D *FEDesc_Obj;
//   QuadFormula2D FaceQuadFormula;
//   TQuadFormula2D *qf2;
//   TEdge *edge;
//   TVertex *Vert; 
//  
//   if(TDatabase::ParamDB->NSTYPE<4)
//    {
//     ErrMsg("slip with friction not yet inplemented in new format !!!!! " << endl);
//     ErrMsg("For slip with friction bc NSTYPE 4 is necessary !!!!! " << endl);
//     exit(4711);
//    } 
//   
//   
//    BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
//    N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D(); 
// //========================================================================
// // loop over all cells
// //========================================================================
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = Coll->GetCell(i);
//     // cout << "Cell No.: " << i << endl;
//     hK = cell->GetDiameter();    
//     
//     //========================================================================
//     // find local used elements on this cell
//     //========================================================================
//     for(j=0;j<N_FeSpaces;j++)
//     {
//       CurrentElement = FeSpaces[j]->GetFE2D(i, cell);
//       LocalUsedElements[j] = CurrentElement;
//       LocN_BF[j] = N_BaseFunct[CurrentElement];
//       LocBF[j] = BaseFuncts[CurrentElement];
//     }
// 
//     N_LocalUsedElements = N_FeSpaces;
// 
//     //========================================================================
//     // calculate values on original element
//     //========================================================================
//     reftrans=TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
//                                     Coll, cell, SecondDer,
//                                     N_Points, xi, eta, zeta, weights,
//                                     X, Y, Z, AbsDetjk);  
// 
//     //========================================================================    
//     // manipulate global matrices
//     // manipulate global right-hand side
//     //========================================================================
//     for(j=0;j<N_Rhs;j++)
//     {
//       // get fe space
//       fespace = FeRhs[j];
//       // get current finite element
//       CurrentElement = fespace->GetFE2D(i, cell);
//       // get number of basis function of current finite element
//       N_ = N_BaseFunct[CurrentElement];
// 
//       // get information on rhs 
//       local_rhs = rhsaux+j*MaxN_BaseFunctions2D;
//       RHS = Rhs[j];
//       
//       N_H = fespace->GetN_Hanging();
//       if(N_H>0)
//        {
//        CurrentHangingRhs = HangingRhs[j];	 
//        }
// 
//     
//       // find bounds in fe space
//       ActiveBound = fespace->GetActiveBound();
//       DirichletBound = fespace->GetHangingBound();
// 
//       // dof of the rhs nodes connected to this cell
//       DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];
// 
//       // only for faces on the boundary
//       BoundaryCondition = BoundaryConditions[j];
//       BoundaryValue = BoundaryValues[j];
//       ele = TFEDatabase2D::GetFE2D(CurrentElement);
//       nf = ele->GetNodalFunctional2D();
// //       nf->GetPointsForFace(N_EdgePoints, EdgePoints1, EdgePoints2);
// 
//       FEDesc_Obj = ele->GetFEDesc2D();
//       N_JointDOF = FEDesc_Obj->GetN_JointDOF(); 
//       
//       
//       
//       
//       
//     } // for j
//     
//     
//     
//     
//     
//   } //for(i=0;i<N_Cells;i++)
//   
//   
//   
//   
//   cout << " Slip " << TDatabase::ParamDB->NSTYPE << endl;
//   exit(0);
// }
// 
// 
// 
