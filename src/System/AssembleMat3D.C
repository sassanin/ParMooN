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
* @brief     source file for TAssembleMat3D
* @author    Sashikumaar Ganesan
* @date      30.05.15
* @History   
 ************************************************************************  */
#include <string.h>
#include <AssembleMat3D.h>
#include <SquareMatrix3D.h>
#include <Matrix3D.h>
#include <FESpace3D.h>
#include <AuxParam3D.h>
#include <DiscreteForm3D.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <BaseCell.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <BoundFace.h>
#include <BoundComp3D.h>
#include <BDEdge3D.h>
#include <Convolution.h>
#include <TNSE3D_FixPo.h>
#include <stdlib.h>
TAssembleMat3D::TAssembleMat3D(int n_fespaces, TFESpace3D **fespaces,
                         int n_sqmatrices, TSquareMatrix3D **sqmatrices,
                         int n_matrices, TMatrix3D **matrices,
                         int n_rhs, double **rhs, TFESpace3D **ferhs,
                         TDiscreteForm3D *discreteform,
                         BoundCondFunct3D **boundarybonditions,
                         BoundValueFunct3D **boundaryvalues,
                         TAuxParam3D *parameters)
                  :TAssemble(n_fespaces, n_sqmatrices, n_matrices, n_rhs, rhs)
{
 int i, j, N;
 
  TFESpace3D *fespace;

   DiscreteForm = discreteform;
   BoundaryConditions = boundarybonditions;
   BoundaryValues = boundaryvalues;
   AuxParam = parameters; 
   
   N_Parameters = AuxParam->GetN_Parameters();   
   SecondDer = discreteform->GetNeeds2ndDerivatives();
  
//          bool *SecondDer;
//       SecondDer = DiscreteFormARhs->GetNeeds2ndDerivatives();
//       cout << " TAssembleMat3D SecondDer " << SecondDer[0] << endl;
      
      
  if(n_fespaces)
   { FeSpaces = new TFESpace3D *[n_fespaces]; }
  if(n_rhs)
   { FeRhs = new TFESpace3D *[n_rhs]; }
  if(n_sqmatrices)
   { SqMatrices = new TSquareMatrix3D* [n_sqmatrices]; }
  if(n_matrices)
   { RecMatrices = new TMatrix3D *[n_matrices]; } 
    
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
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetTestSpace();
       
      TestGlobalNumbers[i] = fespace->GetGlobalNumbers();
      TestBeginIndex[i] = fespace->GetBeginIndex();
  
      fespace = (TFESpace3D *) matrices[i]->GetStructure()->GetAnsatzSpace();
      AnsatzGlobalNumbers[i] = fespace->GetGlobalNumbers();
      AnsatzBeginIndex[i] = fespace->GetBeginIndex();
     } // endfor
    
    for(i=0;i<N_Rhs;i++)
     {
      fespace = ferhs[i];
      RhsBeginIndex[i] = fespace->GetBeginIndex();
      RhsGlobalNumbers[i] = fespace->GetGlobalNumbers();
     } // endfor   
   
   
} // TAssembleMat3D

TAssembleMat3D::~TAssembleMat3D()
{

   
   
}

void TAssembleMat3D::Init()
{
 int i, j, N;
 TFESpace3D *fespace;
  
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
    { rhsaux = new double [N_Rhs*MaxN_BaseFunctions3D]; }
    
    for(i=0;i<N_Rhs;i++)
      LocRhs[i] = rhsaux+i*MaxN_BaseFunctions3D;    
  
    if(N_Parameters)
     {
      paramaux = new double [MaxN_QuadPoints_3D*N_Parameters];
      
      for(j=0;j<MaxN_QuadPoints_3D;j++)
       Param[j] = paramaux + j*N_Parameters;
     }    

    // 20 <= number of term in bilinear form
    // DUE NOTE CHANGE 20 SINCE THE ENTRY 19 IS USED IN GetLocalForms
    auxarray = new double [MaxN_QuadPoints_3D*20]; 
    for(j=0;j<MaxN_QuadPoints_3D;j++)
     AuxArray[j] = auxarray + j*20;    
 
   if(N_AllMatrices)
    {
      auxmat = new double [N_AllMatrices*MaxN_BaseFunctions3D*MaxN_BaseFunctions3D];
      Matrices = new double* [N_AllMatrices*MaxN_BaseFunctions3D];
      
      for(j=0;j<N_AllMatrices*MaxN_BaseFunctions3D;j++)
        Matrices[j] = auxmat+j*MaxN_BaseFunctions3D;

      for(i=0;i<N_AllMatrices;i++)
        LocMatrices[i] = Matrices+i*MaxN_BaseFunctions3D;
    } // endif N_AllMatrices      
  
} // Init


void TAssembleMat3D::DeAllocate()
{
 int i, j, N;
 TFESpace3D *fespace;  
  
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

void TAssembleMat3D::Reset()
{  
//   cout<<"coming to reset"<<endl;
//   cout<<"N_SqMatrices"<< N_SqMatrices<<endl;
//   cout<<"N_Matrices"<<N_Matrices <<endl;
//   cout<<"N_Rhs"<<N_Rhs <<endl;
  
  
  
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
} // TAssembleMat3D::Reset()


void TAssembleMat3D::Assemble3D()
{
  int i, j, k, ij;
  int N_Edges, N_VertInCell, *N_BaseFunct, N_LocalUsedElements, N_Vertex;
  int LocN_BF[N_BaseFuncts3D], N_Points;
  
  double hK;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *weights, *xi, *eta, *zeta;
  
  TBaseCell *cell, *neigh;
  TJoint *joint;
  TBoundFace *boundface;
  BaseFunct3D *BaseFuncts;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D LocBF[N_BaseFuncts3D];
  RefTrans3D reftrans;
   
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
  
  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  

  double delta, mu;
//========================================================================
// loop over all cells
//========================================================================
  
  double *param;
  for(i=0;i<N_Cells;i++)
  {
//     cout<<"calculating the viscosity:: "<<i<<endl;
    cell = Coll->GetCell(i);
    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
	case 0: // diameter
	    hK = cell->GetDiameter();
	    //OutPut("diameter " << hK <<endl);
	    break;

	case 2: // shortest edge
 	    hK = cell->GetShortestEdge();
	    //OutPut("shortest edge " << hK <<endl);
	    break;
	    
	case 1: // with reference map
	case 3: // measure
 	    hK = cell->GetMeasure();
	    hK = pow(hK,1.0/3.0);
	    //OutPut("measure " << hK <<endl);
	    break;
	    
	default: // diameter
	    hK = cell->GetDiameter();
	    break;
    }

    //========================================================================
    // find local used elements on this cell
    //========================================================================
    for(j=0;j<N_FeSpaces;j++)
    {
      CurrentElement = FeSpaces[j]->GetFE3D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = N_FeSpaces;
    //========================================================================
    // calculate values on original element
    //========================================================================
 
      cout << " SecondDer Assemble 3D: " << SecondDer[0] << endl; 
      
      
    reftrans=TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, zeta, weights,
                                    X, Y, Z, AbsDetjk);    

    AuxParam->GetParameters(N_Points, cell, i, xi, eta, zeta, X, Y, Z, Param);
      cout << " SecondDer Assemble 3D: " << SecondDer[0] << endl; 
    if(TDatabase::ParamDB->DISCTYPE== VMS_PROJECTION && !TDatabase::ParamDB->ASSEMBLEMESHMAT)
     {
      //calculating turbulent viscosity variant
      for(int il=0;il<N_Points;il++)
       {
        param = Param[il]; 
        delta =  CharacteristicFilterWidth(hK);
        mu = TurbulentViscosity3D(delta,&param[3],&param[0], &param[0],&param[12],
                                  &param[13],&param[14], param[21]);
 
         if(TDatabase::ParamDB->viscosity_max < mu)
          {
           TDatabase::ParamDB->viscosity_max = mu;
           }
   
         if(TDatabase::ParamDB->viscosity_min > mu)
          {
           TDatabase::ParamDB->viscosity_min = mu;
          }  
        } // for(int il=0;il<N_Points;il++) for calculating turbulent viscosity
     } // if(TDatabase::ParamDB->DISCTYPE==== VMS_PROJECTION && !TDatabase::ParamDB->ASSEMBLEMESHMAT)
  }

#ifdef _MPI
if(TDatabase::ParamDB->DISCTYPE== VMS_PROJECTION && !TDatabase::ParamDB->ASSEMBLEMESHMAT)
 {
  double send = TDatabase::ParamDB->viscosity_max, recieve;
  MPI_Allreduce(&send, &recieve, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  TDatabase::ParamDB->viscosity_max = recieve;

  send = TDatabase::ParamDB->viscosity_min;
  MPI_Allreduce(&send, &recieve, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  TDatabase::ParamDB->viscosity_min = recieve;
 }
#endif
  

// #ifdef _MPI
// int rank;
// MPI_Comm_rank(MPI_COMM_WORLD,&rank);
// if(rank==0)
// #endif
// {
// /*   if(TDatabase::ParamDB->DISCTYPE ==9)
//   {
//   cout<<"viscosity_max"<<TDatabase::ParamDB->viscosity_max<<endl;
//   cout<<"viscosity_min"<<TDatabase::ParamDB->viscosity_min<<endl;
//   }*/ 
// }  
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // cout << "Cell No.: " << i << endl;

    switch (TDatabase::ParamDB->CELL_MEASURE)
    {
	case 0: // diameter
	    hK = cell->GetDiameter();
	    //OutPut("diameter " << hK <<endl);
	    break;

	case 2: // shortest edge
 	    hK = cell->GetShortestEdge();
	    //OutPut("shortest edge " << hK <<endl);
	    break;
	    
	case 1: // with reference map
	case 3: // measure
 	    hK = cell->GetMeasure();
	    hK = pow(hK,1.0/3.0);
	    //OutPut("measure " << hK <<endl);
	    break;
	    
	default: // diameter
	    hK = cell->GetDiameter();
	    break;
    }
    
    //========================================================================
    // find local used elements on this cell
    //========================================================================
    for(j=0;j<N_FeSpaces;j++)
    {
      CurrentElement = FeSpaces[j]->GetFE3D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = N_FeSpaces;
    //========================================================================
    // calculate values on original element
    //========================================================================

//     OutPut("CELL " << i << endl);
    reftrans=TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, 
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, zeta, weights,
                                    X, Y, Z, AbsDetjk);    
    
//      cout << "AbsDetjk: " << AbsDetjk[0] << endl;
//     OutPut("params " << TDatabase::ParamDB->INTERNAL_LEVEL << endl);
    AuxParam->GetParameters(N_Points, cell, i, xi, eta, zeta, X, Y, Z, Param); 
    
    if ((TDatabase::ParamDB->DISCTYPE == SDFEM)
        || (TDatabase::ParamDB->BULK_REACTION_DISC == SDFEM)
        || (TDatabase::ParamDB->CELL_MEASURE == 4)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 105)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 108)
        || (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE == 109))
    {
      N_Vertex = cell->GetN_Vertices();
      for (ij=0;ij<N_Vertex;ij++)
      {
        TDatabase::ParamDB->INTERNAL_VERTEX_X[ij] = cell->GetVertex(ij)->GetX();
        TDatabase::ParamDB->INTERNAL_VERTEX_Y[ij] = cell->GetVertex(ij)->GetY();
        TDatabase::ParamDB->INTERNAL_VERTEX_Z[ij] = cell->GetVertex(ij)->GetZ();
      }
      if (N_Vertex==4)
        TDatabase::ParamDB->INTERNAL_VERTEX_X[4] = -4711;
      TDatabase::ParamDB->INTERNAL_HK_CONVECTION = -1;
    }    
    
     // use DiscreteForm to assemble a few matrices and 
    // right-hand sides at once

//     OutPut("local form " << i << endl);
    if(DiscreteForm)
      DiscreteForm->GetLocalForms(N_Points, weights, AbsDetjk, 
                                    hK, X, Y, Z,
                                    LocN_BF, LocBF,
                                    Param, AuxArray,
                                    cell,
                                    N_AllMatrices, N_Rhs,
                                    LocMatrices, LocRhs);    
    //========================================================================    
    // add local matrices to global matrices (ansatz == test)       
    //======================================================================== 
//     cout<<"N_SqMatrices"<< N_SqMatrices<<endl; 
    if(N_SqMatrices)
     this->AddLocalSqMatToGlobal(i, cell, N_BaseFunct);

    //========================================================================       
    // add local matrices to global matrices (ansatz != test)
    //========================================================================
//     cout<<"N_Matrices"<< N_Matrices<<endl;
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
   
//   cout << "TAssembleMat3D Assemble3D() done ! " << endl;
//   exit(0);
} // Assemble3D

void TAssembleMat3D::PrintAllMat()
{
 int i, j, k, end, *RowPtr, *ColInd, N_Rows, ActiveBound;
 
 double *Entries, *RHS;
 
 TFESpace3D *fespace;  
  
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
  
  
} //TAssembleMat3D::PrintAllMat()


void TAssembleMat3D::ModifyMatHang()
{
 
 int i, j, k, l, m, n, l1, l2, last, N_, end, *DOF, N_Hanging;
 int *HangingRowPtr, *HangingColInd, *RowPtr, *ColInd, ActiveBound;
  
 double *Coupling, v, *Entries, *RHS;
 double *CurrentHangingEntries, *CurrentHangingRhs;

  THangingNode *hn, **HangingNodes;
  HNDesc HNDescr;
  THNDesc *HNDescr_Obj;
  TFESpace3D *fespace;  

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
      HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
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
            } // endfor l2
          } // endif
        } // endfor l
      } // endfor n
    } // endfor i    
   } // for(j=0;j<N_SqMatrices;
    
    
  for(j=0;j<N_Rhs;j++)
  {
    fespace = FeRhs[j];
    N_Hanging = fespace->GetN_Hanging();
    // there are no hanging nodes
    if (N_Hanging == 0)
      continue;
    HangingNodes = fespace->GetHangingNodes();

    RHS = Rhs[j];
    CurrentHangingRhs = HangingRhs[j];

    ActiveBound = fespace->GetActiveBound();

    for(i=0;i<N_Hanging;i++)
    {
      hn = HangingNodes[i];
      HNDescr = hn->GetType();
      HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
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
      } // endfor k
    } // endfor i
  } // endfor j
  
  
  for(j=0;j<N_SqMatrices;j++)
  {
    fespace = SqMatrices[j]->GetFESpace();
    N_ = fespace->GetN_Hanging();
    // there are no hanging nodes
    if (N_ == 0)
      continue;
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
      HNDescr_Obj = TFEDatabase3D::GetHNDesc3D(HNDescr);
      k = HNDescr_Obj->GetN_Nodes();
      Coupling = HNDescr_Obj->GetCoeff();
      DOF = hn->GetDOF();

      Entries[n] = 1.0;
      n++;

      for(l=0;l<k;l++)
      {
        Entries[n] = - Coupling[l];
        n++;
      } // endfor l

    } // endfor i
  } // endfor j  

//    cout << "TAssembleMat3D ModifyMatHang() done ! " << endl;
//   exit(0);   
} // TAssembleMat3D::ModifyMatHang()


void TAssembleMat3D::AddLocalSqMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct)
{
 int j, k, m, n, l, N_, l2, end, *RowPtr, *ColInd, *HangingRowPtr, *HangingColInd;
 int *DOF, ActiveBound, DirichletBound, H_N=0;
 
 double **Matrix, *MatrixRow, *CurrentHangingEntries, *Entries;
 
 TFESpace3D *fespace;  
 FE3D CurrentElement;  
 
    for(j=0;j<N_SqMatrices;j++)
    {
      // find space for this bilinear form
      fespace = SqMatrices[j]->GetFESpace();
      CurrentElement = fespace->GetFE3D(i, cell);
      N_ = N_BaseFunct[CurrentElement];

      Matrix = Matrices+j*MaxN_BaseFunctions3D;
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
    
    
//     cout << N_SqMatrices << " TAssembleMat3D AddLocalSqMatToGlobal() done ! " << endl;
//   exit(0);
} //TAssembleMat3D::AddLocalSqMatToGlobal()

void TAssembleMat3D::AddLocalRecMatToGlobal(int i, TBaseCell *cell, int *N_BaseFunct)
{
 int j, k, l, m, n, end, N_Test, N_Ansatz, *RowPtr, *ColInd, *TestDOF, *AnsatzDOF;

 double **Matrix, *MatrixRow, *Entries;
 
 FE3D TestElement, AnsatzElement;
  
    for(j=0;j<N_Matrices;j++)
    {
      TestElement = ((TFESpace3D *) RecMatrices[j]->GetStructure()->GetTestSpace())->GetFE3D(i, cell);
      AnsatzElement = ((TFESpace3D *) RecMatrices[j]->GetStructure()->GetAnsatzSpace())->GetFE3D(i, cell);

      // cout << "non square matrix: " << j << endl;
      // cout << "TestElement: " << TestElement << endl;
      // cout << "AnsatzElement: " << AnsatzElement << endl;

      N_Test = N_BaseFunct[TestElement];
      N_Ansatz = N_BaseFunct[AnsatzElement];

      Matrix = Matrices+(j+N_SqMatrices)*MaxN_BaseFunctions3D;

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
  
} //   TAssembleMat3D::AddLocalMatToGlobal


void TAssembleMat3D::AddLocalRhsToGlobal(int i, TBaseCell *cell, int *N_BaseFunct, BaseFunct3D *BaseFuncts, RefTrans3D reftrans)
{
 int j, k, l, m, l3, N_, l1, dof, ActiveBound, DirichletBound, *DOF, N_H;
 int comp, N_EdgeDOF, N_Joints, *JointDOF, *EdgeDOF, N_Edges, N_VertInCell;
 const int *TmpFV, *TmpLen, *EdgeVertex, *TmpFE, *ETmpLen;
 int MaxLen, EMaxLen, N_FaceEdges, N_Points, Bd_ID;
  
 double xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3;
 double nx, ny, nz, t0; 
 double *local_rhs, *RHS, *CurrentHangingRhs;
 double xf, yf, zf, *t, *s, *xi, *eta, *zeta;
 double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D], Z[MaxN_QuadPoints_3D];
 double LinComb[4], *FaceWeights;
 double PointValues[MaxN_PointsForNodal3D];
 double FunctionalValues[MaxN_BaseFunctions3D];
 double AbsDetjk[MaxN_QuadPoints_3D];
 double **JointValues, *JointValue, t1, t2;
  
 bool InnerBoundary, OuterBoundary;
  
 TFESpace3D *fespace;  
 FE3D CurrentElement;  
 BoundCond Cond0, Cond1;
 BoundCondFunct3D *BoundaryCondition;
 BoundValueFunct3D *BoundaryValue; 
 TNodalFunctional3D *nf;
 TFE3D *ele;
 TFEDesc3D *FEDesc_Obj;
 TJoint *joint;
 TBoundFace *boundface;
 QuadFormula2D FaceQuadFormula;
 TQuadFormula2D *qf2;
 TEdge *edge;
 TVertex *Vert; 
 TBoundFace *BdFace;
 TBoundComp3D *BdComp;
 TBDEdge3D * bd_edge;
 
   for(j=0;j<N_Rhs;j++)
    {
      local_rhs = rhsaux+j*MaxN_BaseFunctions3D;
      RHS = Rhs[j];
      
      fespace = FeRhs[j];
      CurrentElement = fespace->GetFE3D(i, cell);
      N_ = N_BaseFunct[CurrentElement];
      ele = TFEDatabase3D::GetFE3D(CurrentElement);
      nf = ele->GetNodalFunctional3D();
      FEDesc_Obj = ele->GetFEDesc3D();
      N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
      N_Joints = cell->GetN_Faces();
   
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
        l=DOF[m];
        // cout << "DOF: " << l << endl;
        if(l<ActiveBound)
        {
          // node l is inner or Neumann node
          RHS[l] += local_rhs[m]; 
        } // endif l
        else
        {
          if(l<DirichletBound)
          {
            // hanging node
            if(N_H==0) continue;

            l -= ActiveBound;
            CurrentHangingRhs[l] += local_rhs[m];
          }
        }
      } // endfor m  
      
      // setting Dirichlet boundary condition
      for(m=0;m<N_Joints;m++)
      {
        joint = cell->GetJoint(m);
        InnerBoundary = false;
        OuterBoundary = false;

        if(joint->GetType() == BoundaryFace || joint->GetType() == IsoBoundFace)
         { OuterBoundary = true; }

        if(InnerBoundary || OuterBoundary)
        {          
          cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
          cell->GetShapeDesc()->GetFaceEdge(TmpFE, ETmpLen, EMaxLen);
          t0 = 1.0/TmpLen[m];
          
          xf = 0; yf = 0; zf = 0;
          for (l1=0;l1<TmpLen[m];l1++)
          {
            cell->GetVertex(TmpFV[m*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
	    
            //LinComb[l1] = t0;
            xf += t0*X[l1];
            yf += t0*Y[l1];
            zf += t0*Z[l1];
          }

          BdFace = (TBoundFace *)joint;
          BdComp    = BdFace->GetBoundComp();
          Bd_ID = BdComp->GetID();
          
          // the face gets the b.c. which is valid at its center
          BoundaryCondition(Bd_ID, xf, yf, zf, Cond0);          
  
          switch(Cond0)
          {
            case DIRICHLET:
 
              if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
              {
                nf->GetPointsForFace(N_Points, t, s);
  
                for(l1=0;l1<N_Points;l1++)
                {
                  switch(TmpLen[m])
                  {
                    case 4:
                      LinComb[0] = (1.-t[l1])*(1.-s[l1]);
                      LinComb[1] = t[l1]*(1.-s[l1]);
                      LinComb[2] = t[l1]*s[l1];
                      LinComb[3] = (1.-t[l1])*s[l1];
                      xf = LinComb[0]*X[0] + LinComb[1]*X[1] + LinComb[2]*X[2] + LinComb[3]*X[3];
                      yf = LinComb[0]*Y[0] + LinComb[1]*Y[1] + LinComb[2]*Y[2] + LinComb[3]*Y[3];
                      zf = LinComb[0]*Z[0] + LinComb[1]*Z[1] + LinComb[2]*Z[2] + LinComb[3]*Z[3];
                    break;

                    case 3:
                      LinComb[0] = 1-t[l1]-s[l1];
                      LinComb[1] = t[l1];
                      LinComb[2] = s[l1];
                      xf = LinComb[0]*X[0] + LinComb[1]*X[1] + LinComb[2]*X[2];
                      yf = LinComb[0]*Y[0] + LinComb[1]*Y[1] + LinComb[2]*Y[2];
                      zf = LinComb[0]*Z[0] + LinComb[1]*Z[1] + LinComb[2]*Z[2];
                    break;
                  }

                  //cout << l1 << ": " << t0 << " " << t1 << " <> ";
                  //cout << xf << " " << yf << " " << zf << endl;
  
                  BoundaryValue(Bd_ID, xf, yf, zf, PointValues[l1]);
                  // cout << " PV: " << PointValues[l1] << endl;
                }
              }
              else
              {
                // no isoparametric
                nf->GetPointsForFace(m, N_Points, xi, eta, zeta);

                TFEDatabase3D::GetOrigFromRef(reftrans, N_Points, xi, eta, zeta, X, Y, Z, AbsDetjk);
  
                for(l1=0;l1<N_Points;l1++)
                {
                  BoundaryValue(Bd_ID, X[l1], Y[l1], Z[l1], PointValues[l1]);
                  // PointValues[l1] = 0;
                  //cout << " PV: " << PointValues[l1] << endl;
                }
              }
              
              nf->GetFaceFunctionals(Coll, cell, m, PointValues, FunctionalValues);
              JointDOF = FEDesc_Obj->GetJointDOF(m);
              N_EdgeDOF = FEDesc_Obj->GetN_JointDOF();
              for(l=0;l<N_EdgeDOF;l++)
              {
                RHS[DOF[JointDOF[l]]] = FunctionalValues[l];
                // cout << " DOF: " << DOF[JointDOF[l]];
                // cout << " FV: " << FunctionalValues[l] << endl;
              }
            break;

            case NEUMANN:
            case ROBIN:  // lhs of robin will be assembled in main program
              // cout << "Neumann condition in Assemble3D" << endl;
              l = TFEDatabase3D::GetPolynomialDegreeFromFE3D(CurrentElement);

              switch(TmpLen[m])
              {
                case 3:
                  // triangular face
                  FaceQuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*l);
//                   FaceQuadFormula = Gauss3Tria;
                break;
                case 4:
                  // quadrilateral face
		    FaceQuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*l);
                break;
              }
              
              qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQuadFormula);
              qf2->GetFormulaData(N_Points, FaceWeights, t, s);
              // generate data on reference mesh cell for the 2d face of 3d cell
	      
// OMP!!!! needs to be atomic 	      
              TFEDatabase3D::GetBaseFunct3DFromFE3D(CurrentElement)->MakeRefElementData(FaceQuadFormula);
	      
              // values of base functions in all quadrature points on face 
              JointValues = TFEDatabase3D::GetJointValues3D(BaseFuncts[CurrentElement], FaceQuadFormula, m);
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])->ChangeBF(Coll, cell, N_Points, JointValues);

              switch(TmpLen[m])
              {
                case 3:
                  xc1 = X[1] - X[0];
                  xc2 = X[2] - X[0];

                  yc1 = Y[1] - Y[0];
                  yc2 = Y[2] - Y[0];

                  zc1 = Z[1] - Z[0];
                  zc2 = Z[2] - Z[0];

                  // normal vector
                  nx = yc1*zc2 - zc1*yc2;
                  ny = zc1*xc2 - xc1*zc2;
                  nz = xc1*yc2 - yc1*xc2;
                  // determinant of reference trafo
                  t2 = sqrt(nx*nx + ny*ny + nz*nz);
                  
                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = t[l];
                    t1 = s[l];
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = 1-t0-t1;
                    LinComb[1] = t0;
                    LinComb[2] = t1;
                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]+LinComb[2]*X[2];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]+LinComb[2]*Y[2];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]+LinComb[2]*Z[2];

                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(Bd_ID, xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*t2;
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                  } // endfor l
                break;
                
                case 4:
                  xc1=(-X[0] + X[1] + X[2] - X[3]) * 0.25;
                  xc2=(-X[0] - X[1] + X[2] + X[3]) * 0.25;
                  xc3=( X[0] - X[1] + X[2] - X[3]) * 0.25;

                  yc1=(-Y[0] + Y[1] + Y[2] - Y[3]) * 0.25;
                  yc2=(-Y[0] - Y[1] + Y[2] + Y[3]) * 0.25;
                  yc3=( Y[0] - Y[1] + Y[2] - Y[3]) * 0.25;

                  zc1=(-Z[0] + Z[1] + Z[2] - Z[3]) * 0.25;
                  zc2=(-Z[0] - Z[1] + Z[2] + Z[3]) * 0.25;
                  zc3=( Z[0] - Z[1] + Z[2] - Z[3]) * 0.25;

                  for(l=0;l<N_Points;l++)
                  {
                    JointValue = JointValues[l];
                    t0 = 0.5*(t[l]+1);
                    t1 = 0.5*(s[l]+1);
                    // cout << "t: " << t0 << " " << t1 << endl;
                    LinComb[0] = (1-t0)*(1-t1);
                    LinComb[1] = t0*(1-t1);
                    LinComb[2] = t0*t1;
                    LinComb[3] = (1-t0)*t1;
                    xf = LinComb[0]*X[0] + LinComb[1]*X[1]+LinComb[2]*X[2] + LinComb[3]*X[3];
                    yf = LinComb[0]*Y[0] + LinComb[1]*Y[1]+LinComb[2]*Y[2] + LinComb[3]*Y[3];
                    zf = LinComb[0]*Z[0] + LinComb[1]*Z[1]+LinComb[2]*Z[2] + LinComb[3]*Z[3];

                    // cout << xf << " " << yf << " " << zf << endl;
                    BoundaryValue(Bd_ID, xf, yf, zf, t0);
                    // cout << "PV: " << t0 << endl;
                    nx = (yc1+s[l]*yc3)*(zc2+t[l]*zc3)-(zc1+s[l]*zc3)*(yc2+t[l]*yc3);
                    ny = (zc1+s[l]*zc3)*(xc2+t[l]*xc3)-(xc1+s[l]*xc3)*(zc2+t[l]*zc3);
                    nz = (xc1+s[l]*xc3)*(yc2+t[l]*yc3)-(yc1+s[l]*yc3)*(xc2+t[l]*xc3);

                    t1 = nx*nx+ny*ny+nz*nz;
                    // cout << t1 << endl;
                    t0 *= FaceWeights[l]*sqrt(t1);
                    for(k=0;k<N_;k++)
                      if((l3 = DOF[k])<ActiveBound)
                        RHS[l3] += t0*JointValue[k];
                  } // endfor l
                break;
              }
              TFEDatabase3D::GetBaseFunct3D(BaseFuncts[CurrentElement])->ChangeBF(Coll, cell, N_Points, JointValues);
            break;
	    
            case SLIP:
                OutPut("Use SLIP_FRICTION_PENETRATION_RESISTANCE boundary condition !"<< endl);
                exit(4711);
                break;
            case SLIP_FRICTION_PENETRATION_RESISTANCE:
                // do nothing here
                // everything is done in Assemble3DSlipBC, see below
             break;
          default: 
           ErrMsg("Unknown Bundary Type ");
	   exit(0);
	  break;
         } // endswitch Cond0

#ifdef _MPI         
        /** edges on this face are already set, so no need of checking in BD edges on this face */
        N_FaceEdges = ETmpLen[m];
        for(l1=0;l1<N_FaceEdges;l1++)
         {
          edge = cell->GetEdge(TmpFE[m*EMaxLen+l1]);
          edge->SetClipBoard(i);
         } // fo 
         
         /** vertices on this face are already set, so no need of checking in BD vert on this face, 30.06.12, sashi */        
         for (l1=0;l1<TmpLen[m];l1++)
          cell->GetVertex(TmpFV[m*MaxLen+l1])->SetClipBoard(i);
#endif      

	} // if(InnerBoundary || OuterBoundary)
      }//  for(m=0;m<N_Joints;m++)
  
#ifdef _MPI   
    N_EdgeDOF = FEDesc_Obj->GetN_EdgeDOF();
    if(N_EdgeDOF> 0) //conforming FE
     {
      cell->GetShapeDesc()->GetEdgeVertex(EdgeVertex);   
      N_Edges=cell->GetN_Edges();
      
      for(m=0;m<N_Edges;m++)
       {
        edge = cell->GetEdge(m);  
    
        if( (edge->GetType()==BDEdge3D || edge->GetType()==IsoEdge3D) && edge->GetClipBoard()==-1)
         {      
          // BD edge is not yet set 
          edge->SetClipBoard(i);  
          EdgeDOF = FEDesc_Obj->GetEdgeDOF(m);   
          
          for(l=0;l<N_EdgeDOF;l++)
           {
            dof = DOF[EdgeDOF[l]];
            fespace->GetDOFPosition(dof, xf, yf, zf);
	    bd_edge = (TBDEdge3D*)edge;	    
   BoundaryCondition(bd_edge->get_Bd_id(),xf, yf, zf, Cond0);
	    
 
            if(Cond0==DIRICHLET) // nothing to do for nin-Dirichlet
             {
           BoundaryValue(bd_edge->get_Bd_id(),xf, yf, zf, t0);
             RHS[dof] = t0;  
             }       
          }  
          /** vertices on this edge are already set, so no need of checking in BD vert on this edge, 30.06.12, sashi */
          cell->GetVertex(EdgeVertex[2*m])->SetClipBoard(i);
          cell->GetVertex(EdgeVertex[2*m+1])->SetClipBoard(i);      
         }
       }//   endfor m
     } 
     
    /** correct the Vert Dirichlet   */
    if( FEDesc_Obj->GetN_VertDOF() > 0) //conforming FE
     {
      N_VertInCell = cell->GetN_Vertices(); 
     
      for(m=0;m<N_VertInCell;m++)
       {  
        Vert = cell->GetVertex(m);
       
        if(Vert->IsBoundVert() && Vert->GetClipBoard()==-1 )
         {
          // BD vert is not yet set 
          Vert->SetClipBoard(i);
          Vert->GetCoords(xf, yf, zf);   
          BoundaryCondition(Vert->get_Bd_id(), xf, yf, zf, Cond0); 
          
          if(Cond0==DIRICHLET) // nothing to do for nin-Dirichlet
           {
            dof = DOF[FEDesc_Obj->GetVertDOF(m)];   
            BoundaryValue(Vert->get_Bd_id(),xf, yf, zf, t0);
            RHS[dof] = t0;  	    
           }
         }
       }//   endfor m    
      } //   if( FEDesc_Obj->GetN_VertDOF() > 0) 
#endif      

    } //for(j=0;j<N_Rhs;j++)
  
//    cout << "TAssembleMat3D AddLocalRhsToGlobal() done ! " << endl;
//   exit(0); 
} //TAssembleMat3D::AddLocalRhsToGlobal


void TAssembleMat3D::AssembleNavierSlip()
{
  int i, j, k, ij, N_, N_H, ActiveBound, DirichletBound, *DOF;
  int N_Edges, N_VertInCell, *N_BaseFunct, N_LocalUsedElements, N_Vertex;
  int LocN_BF[N_BaseFuncts3D], N_Points, N_JointDOF;

  double hK, *local_rhs, *RHS, *CurrentHangingRhs;
  double X[MaxN_QuadPoints_3D], Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double AbsDetjk[MaxN_QuadPoints_3D];
  double *weights, *xi, *eta, *zeta;

  TFESpace3D *fespace;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  TBoundFace *boundface;
  BaseFunct3D *BaseFuncts;
  FE3D LocalUsedElements[N_FEs3D], CurrentElement;
  BaseFunct3D LocBF[N_BaseFuncts3D];
  RefTrans3D reftrans;
  BoundCond Cond0, Cond1;
  BoundCondFunct3D *BoundaryCondition;
  BoundValueFunct3D *BoundaryValue;
  TNodalFunctional3D *nf;
  TFE3D *ele;
  TFEDesc3D *FEDesc_Obj;
  QuadFormula2D FaceQuadFormula;
  TQuadFormula2D *qf2;
  TEdge *edge;
  TVertex *Vert;

  if(TDatabase::ParamDB->NSTYPE<4)
   {
    ErrMsg("slip with friction not yet inplemented in new format !!!!! " << endl);
    ErrMsg("For slip with friction bc NSTYPE 4 is necessary !!!!! " << endl);
    exit(4711);
   }


   BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
   N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
//========================================================================
// loop over all cells
//========================================================================
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    // cout << "Cell No.: " << i << endl;
    hK = cell->GetDiameter();

    //========================================================================
    // find local used elements on this cell
    //========================================================================
    for(j=0;j<N_FeSpaces;j++)
    {
      CurrentElement = FeSpaces[j]->GetFE3D(i, cell);
      LocalUsedElements[j] = CurrentElement;
      LocN_BF[j] = N_BaseFunct[CurrentElement];
      LocBF[j] = BaseFuncts[CurrentElement];
    }

    N_LocalUsedElements = N_FeSpaces;

    //========================================================================
    // calculate values on original element
    //========================================================================
    reftrans=TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements,
                                    Coll, cell, SecondDer,
                                    N_Points, xi, eta, zeta, weights,
                                    X, Y, Z, AbsDetjk);

    //========================================================================
    // manipulate global matrices
    // manipulate global right-hand side
    //========================================================================
    for(j=0;j<N_Rhs;j++)
    {
      // get fe space
      fespace = FeRhs[j];
      // get current finite element
      CurrentElement = fespace->GetFE3D(i, cell);
      // get number of basis function of current finite element
      N_ = N_BaseFunct[CurrentElement];

      // get information on rhs
      local_rhs = rhsaux+j*MaxN_BaseFunctions2D;
      RHS = Rhs[j];

      N_H = fespace->GetN_Hanging();
      if(N_H>0)
       {
       CurrentHangingRhs = HangingRhs[j];
       }


      // find bounds in fe space
      ActiveBound = fespace->GetActiveBound();
      DirichletBound = fespace->GetHangingBound();

      // dof of the rhs nodes connected to this cell
      DOF = RhsGlobalNumbers[j] + RhsBeginIndex[j][i];

      // only for faces on the boundary
      BoundaryCondition = BoundaryConditions[j];
      BoundaryValue = BoundaryValues[j];
      ele = TFEDatabase3D::GetFE3D(CurrentElement);
      nf = ele->GetNodalFunctional3D();
//       nf->GetPointsForFace(N_EdgePoints, EdgePoints1, EdgePoints2);

      FEDesc_Obj = ele->GetFEDesc3D();
      N_JointDOF = FEDesc_Obj->GetN_JointDOF();





    } // for j


    
    
    
  } //for(i=0;i<N_Cells;i++)
  
  
  
  
  cout << " Slip " << TDatabase::ParamDB->NSTYPE << endl;
  exit(0);
}
