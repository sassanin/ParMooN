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
   
// // =======================================================================
// // @(#)StructureNSE3D.C        1.6 09/17/99
// // 
// // Class:       TStructureNSE3D
// //
// // Purpose:     build and store a structure for a square matrix in 3D
// //
// // Author:      Gunar Matthies
// //
// // History:     29.07.2005 start implementation
// //
// // =======================================================================
// 
// #include <DefineParams.h>
// 
// #include <FEDatabase3D.h>
// #include <TetraAffin.h>
// #include <TetraIsoparametric.h>
// #include <HexaAffin.h>
// #include <HexaTrilinear.h>
// #include <HexaIsoparametric.h>
// #include <StructureNSE3D.h>
// #include <MooNMD_Io.h>
// #include <Joint.h>
// #include <string.h>
// #include <stdlib.h>
// 
// /** generate the matrix structure, both space are 3D */
// TStructureNSE3D::TStructureNSE3D(TFESpace3D *Space)
//   : TSquareStructure3D()
// {
//   int i,j,k,l,n,N_, n1, m; 
//   int *GlobalNumbers;
//   int *BeginIndex;
// 
//   int N_Dirichlet;
//   int N_Inner;
//   int N_Hanging;
//   int *BoundNodeBounds, N_NonDiri, N_BoundaryNodeTypes;
//   int HangingBound;
// 
//   int *AuxPtr, *HangingAuxPtr;
//   int *KColAux, *HangingKColAux;
//   int index, oldindex;
//   int N_KCol;
// 
//   int Offset, *DOF, end, begin;
//   THangingNode **HangingNodes;
//   THangingNode *hn;
// 
//   TCollection *Coll;
// 
//   TBaseCell *cell;
//   FE3D CurrentElement;
//   TFE3D *ele;
//   TFEDesc3D *fedesc;
//   int N_InnerDOFs, *InnerDOF;
//   int isinner;
// 
//   FESpace = Space;
// 
//   ActiveBound=Space->GetN_ActiveDegrees();
//   HangingBound=Space->GetHangingBound();
// 
//   N_BoundaryNodeTypes=Space->GetN_DiffBoundaryNodeTypes();
//   BoundNodeBounds=Space->GetN_BoundaryNodes();
//   N_NonDiri = 0;
//   for(i=0;i<N_BoundaryNodeTypes;i++)
//     N_NonDiri += BoundNodeBounds[i];
// 
//   N_Dirichlet=Space->GetN_Dirichlet();
//   N_Inner=Space->GetN_Inner();
//   N_Hanging=Space->GetN_Hanging();
// 
//   if(N_Hanging)
//   {
//     Error("Hanging nodes cannot be managed right now!" << endl);
//     Error("Program aborted!" << endl);
//     exit(-1);
//   }
// 
//   N_Rows = ActiveBound+N_Hanging+N_Dirichlet;
//   N_Columns = N_Rows;
// 
//   // AuxPtr[i] will contain an upper bound for the number of 
//   // matrix entries in row i
//   l=N_Rows+1;
//   AuxPtr=new int[l];
//   memset(AuxPtr, 0, l*sizeof(int));
// 
//   HangingAuxPtr=NULL;
// 
//   GlobalNumbers=Space->GetGlobalNumbers();
// 
//   BeginIndex=Space->GetBeginIndex();
// 
//   Offset=ActiveBound;
// 
//   // loop over all elements 
//   N_=Space->GetN_Cells();
//   Coll = FESpace->GetCollection();
//   for(i=0;i<N_;i++)
//   {
//     cell = Coll->GetCell(i);
//     DOF=GlobalNumbers+BeginIndex[i];
// 
//     CurrentElement = Space->GetFE3D(i, cell);
//     n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();
// 
//     for(j=0;j<n1;j++)
//     {
//       k=DOF[j];
//       if(k<ActiveBound) 
//         AuxPtr[k]+=n1;
//     } // endfor j
//   } // endfor i
// 
//   // add rows for Dirichlet nodes in  space
//   N_Entries = 0;
//   N_KCol = 0;
//   Offset = ActiveBound+N_Hanging;
//   for(i=0,j=Offset;i<N_Dirichlet;i++,j++)
//   {
//     AuxPtr[j]=1;
//   }
// 
//   // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
//   // KColAux array, the column DOF for row i are in the intervall
//   // [ AuxPtr[i],AuxPtr[i+1] )
//   N_=N_Rows;
//   l=AuxPtr[0];
//   AuxPtr[0]=0;
//   // cout << AuxPtr[0] << endl;
//   for(i=0;i<N_;i++)
//   {
//     k=AuxPtr[i+1];
//     AuxPtr[i+1]=AuxPtr[i]+l;
//     l=k;
//     // cout << AuxPtr[i+1] << endl;
//   }
// 
//   /*
//   cout << endl;
//   cout << "AuxPtr array" << endl;
//   for(i=0;i<=N_;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
//   */
// 
//   // cout << "Upper bound: " << AuxPtr[N_Rows] << endl;
// 
//   // get memory for KColAux array, initialize it with -1
//   l=AuxPtr[N_Rows]; // upper bound for number of matrix entries 
//   KColAux=new int[l];
//   memset(KColAux, -1, sizeof(int)*l);
// 
//   N_=Space->GetN_Cells();
//   Coll = FESpace->GetCollection();
//   for(i=0;i<N_;i++)
//   {
//     cell = Coll->GetCell(i);
//     DOF=GlobalNumbers+BeginIndex[i];
// 
//     CurrentElement = Space->GetFE3D(i, cell);
//     ele = TFEDatabase3D::GetFE3D(CurrentElement);
//     n1 = ele->GetN_DOF();
//     fedesc = ele->GetFEDesc3D();
//     N_InnerDOFs = fedesc->GetN_InnerDOF();
//     InnerDOF = fedesc->GetInnerDOF();
// 
//     for(j=0;j<n1;j++)
//     {
//       // does j belongs to cell inner nodes?
//       isinner = 0;
//       for(k=0;k<N_InnerDOFs;k++)
//       {
//         if(j == InnerDOF[k])
//         {
//           isinner = 1;
//           break;
//         }
//       }
//       if(isinner)
//         continue;
// 
//       n=DOF[j];
// 
//       for(k=0;k<n1;k++)
//       {
//         // does k belongs to cell inner nodes?
//         isinner = 0;
//         for(l=0;l<N_InnerDOFs;l++)
//         {
//           if(k == InnerDOF[l])
//           {
//             isinner = 1;
//             break;
//           }
//         }
// 
//         if(isinner)
//           continue;
// 
//         m=DOF[k];
// 
//         if(n<ActiveBound)
//         {
//           // this node is a real node (inner or Neumann)
//           index=AuxPtr[n];
//           l=KColAux[index];
//           // check whether this column is already in this row
//           while(l!=-1 && l!=m)
//           {
//             index++; l=KColAux[index];
//             oldindex++;
//           }
//           if(l==-1)
//           {
//             // this is a new column for this row
//             KColAux[index] = m;
//             index++;
//             N_Entries += 9;
//             N_KCol++;
//           }
//         }
//       } // endfor k
//     } // endfor j
//   } // endfor i
// 
//   HangingKCol=NULL;
//   HangingRowPtr=HangingAuxPtr;
// 
//   // add Dirichlet rows
//   for(i=0;i<N_Dirichlet;i++)
//   {
//     // cout << "index: " << index << endl;
//     j = i+HangingBound;
//     KColAux[AuxPtr[j]]=j;
//     N_Entries += 3;
//     N_KCol++;
//   }
// 
//   /*
//   // check
//   cout << endl;
//   cout << "check" << endl;
//   N_=ActiveBound;
//   N_=N_Rows;
//   for(i=0;i<N_;i++)
//   {
//     cout << "Row: " << setw(4) << i << ": ";
//     index=AuxPtr[i+1];
//     for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
//       cout << setw(4) << KColAux[j];
//     cout << endl;
//   }
//   */
//   
//   /*
//   cout << "Number of matrix entries: ";
//   cout << N_Entries << endl;
//   cout << "Number of KCols: " << N_KCol << endl;
//   cout << endl;
//   */
// 
//   // compress KColAux array to KCol by deleting all -1's
//   // build the RowPtr array
//   N_ = N_Rows;
//   KCol=new int[N_KCol];
//   RowPtr=AuxPtr;
// 
//   index=0;
//   for(i=0;i<N_;i++)
//   {
//     oldindex=index;
//     m=AuxPtr[i+1];
//     // cout << i << " ";
//     // cout << AuxPtr[i] << " " << KColAux[AuxPtr[i]] << endl;
//     for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
//     {
//       KCol[index]=KColAux[j];
//       index++;
//     } // endfor j
//     RowPtr[i]=oldindex;
//     // cout << setw(4) << i << " RowPtr[i]: " << RowPtr[i] << endl;
//   } // endfor i
//   RowPtr[N_] = index;
// 
//   Sort();
// 
//   /*
//   // print out the whole matrix structure
//   cout << endl;
//   N_=N_Rows;
//   for(i=0;i<N_;i++)
//   {
//     cout << RowPtr[i] << " --- " << RowPtr[i+1]-1;
//     cout << " number of entries:" << RowPtr[i+1]-RowPtr[i] << endl;
//     cout << "Row: " << setw(4) << i << ": ";
//     end=RowPtr[i+1];
//     for(j=RowPtr[i];j<end;j++)
//       cout << setw(4) << KCol[j];
//     cout << endl;
//   }
//   */
// 
//   // free KColAux
//   delete KColAux;
// 
//   cout << endl;
//   cout << "Information on the stored StructureNSE3D" << endl;
//   cout << "Number of rows: " << N_Rows << endl;
//   cout << "Number of columns: " << N_Columns << endl;
//   cout << "Number of matrix entries: " << N_Entries << endl;
// 
//   BeginJb = NULL;
//   jb = NULL;
//   N_DOFperJoint = 0;
//   Alpha = NULL;
// 
//   GenerateAlpha();
// } 
// 
// TStructureNSE3D::~TStructureNSE3D()
// {
//   if(HangingKCol) delete HangingKCol;
//   if(HangingRowPtr) delete HangingRowPtr;
// }
// 
// void TStructureNSE3D::GenerateAlpha()
// {
//   int i, j, k, l, m, n;
//   int N_Cells, N_Joints;
//   int N_JointDOF;
//   int N_AllFaces, N_AllFaceDOF;
//   TCollection *coll;
//   TBaseCell *cell, *neigh;
//   TJoint *joint; 
//   FE3D feid, feidNeigh;
//   TFEDesc3D *fedesc, *fedescNeigh;
//   int **JointDOFs, *JointDOF, *JointDOFNeigh;
//   int *GlobalNumbers, *BeginIndex, *DOF, *DOFNeigh;
//   BF3DRefElements refele;
//   RefTrans3D reftrans;
//   TRefTrans3D *rt;
//   TQuadFormula2D *qf;
//   QuadFormula2D QuadFormula;
//   int N_Points;
//   double *Weights, *xi, *eta;
//   double **JointValues, *JointValue;
//   BaseFunct3D *BaseFuncts;
//   double x0, y0, x1, y1;
//   double n1, n2, n3, len;
//   double hE, s;
//   double t11, t12, t13, t21, t22, t23;
//   double PhinE[3*MaxN_BaseFunctions3D];
//   int N_JointsNeigh, NeighNumber, JointNumberNeigh;
//   int N_LocalDOF, N_LocalDOFNeigh;
//   int N_TotalDOF;
//   int JbNumber;
// 
// #ifdef __3D__
//   double z0, z1;
// #endif
// 
//   BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
// 
//   BeginIndex = FESpace->GetBeginIndex();
//   GlobalNumbers = FESpace->GetGlobalNumbers();
//   N_TotalDOF = FESpace->GetN_DegreesOfFreedom();
//   
//   coll = FESpace->GetCollection();
//   N_Cells = coll->GetN_Cells();
// 
//   BeginJb = new int[N_Cells+1];
//   memset(BeginJb, 0, (N_Cells+1)*SizeOfInt);
// 
//   N_AllFaces = 0;
//   BeginJb[0] = 0;
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = coll->GetCell(i);
//     cell->SetClipBoard(i);
//     // OutPut("Cell: " << (long int)cell << endl);
//     N_Joints = cell->GetN_Joints();
//     N_AllFaces += N_Joints;
//     BeginJb[i+1] = N_AllFaces;
// 
//     feid = FESpace->GetFE3D(i, cell);
//     fedesc = TFEDatabase3D::GetFE3D(feid)->GetFEDesc3D();
// 
//     N_JointDOF = fedesc->GetN_JointDOF();
//     if(i==0)
//       N_DOFperJoint = N_JointDOF;
//     else
//       if(N_DOFperJoint != N_JointDOF)
//       {
//         Error("All elements must have the same number of d.o.f. on the joints!" << endl);
//         Error("Programs terminates.");
//         exit(-1);
//       }
//   }
// 
//   // we have three components
//   N_AllFaceDOF = 3*N_DOFperJoint*N_AllFaces;
// 
//   // OutPut("N_AllFaces: " << N_AllFaces << endl);
//   // OutPut("N_AllFaceDOF: " << N_AllFaceDOF << endl);
//   // OutPut(endl);
// 
//   jb = new int[N_AllFaces];
//   Alpha = new double[N_AllFaceDOF];
// 
//   for(i=0;i<N_Cells;i++)
//   {
//     cell = coll->GetCell(i);
//     // OutPut(endl << "Cell: " << i << endl);
// 
//     feid = FESpace->GetFE3D(i, cell);
//     fedesc = TFEDatabase3D::GetFE3D(feid)->GetFEDesc3D();
// 
//     N_LocalDOF = fedesc->GetN_DOF();
//     JointDOFs = fedesc->GetJointDOF();
// 
//     DOF = GlobalNumbers + BeginIndex[i];
// 
//     // find right quadrature rule
//     k = TFEDatabase3D::GetPolynomialDegreeFromFE3D(feid);
//     refele = TFEDatabase3D::GetRefElementFromFE3D(feid);
//     switch(refele)
//     {
//       case BFUnitTetrahedron:
//         QuadFormula = TFEDatabase3D::GetQFTriaFromDegree(2*k);
//       break;
// 
//       case BFUnitHexahedron:
//         QuadFormula = TFEDatabase3D::GetQFQuadFromDegree(2*k);
//       break;
// 
//       default:
//         OutPut("Unknown reference element. Program aborted." << endl);
//         exit(-1);
//     }
//     qf = TFEDatabase3D::GetQuadFormula2D(QuadFormula);
//     qf->GetFormulaData(N_Points, Weights, xi, eta);
// 
//     reftrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(feid);
//     switch(reftrans)
//     {
//       case TetraAffin:
//         rt = (TTetraAffin *)(TFEDatabase3D::GetRefTrans3D(reftrans));
//         ((TTetraAffin *)rt)->SetCell(cell);
//       break;
// 
//       case TetraIsoparametric:
//         rt = (TTetraIsoparametric *)(TFEDatabase3D::GetRefTrans3D(reftrans));
//         ((TTetraIsoparametric *)rt)->SetCell(cell);
//       break;
// 
//       case HexaAffin:
//         rt = (THexaAffin *)(TFEDatabase3D::GetRefTrans3D(reftrans));
//         ((THexaAffin *)rt)->SetCell(cell);
//       break;
// 
//       case HexaTrilinear:
//         rt = (THexaTrilinear *)(TFEDatabase3D::GetRefTrans3D(reftrans));
//         ((THexaTrilinear *)rt)->SetCell(cell);
//       break;
// 
//       case HexaIsoparametric:
//         rt = (THexaIsoparametric *)(TFEDatabase3D::GetRefTrans3D(reftrans));
//         ((THexaIsoparametric *)rt)->SetCell(cell);
//       break;
//       
//       default:
//         OutPut("Wrong reference transformation. Program aborted." << endl);
//         exit(-1);
//     }
// 
//     // calculate scalar base functions on each joint in all points
//     // which are used by the quadrature rule
//     TFEDatabase3D::GetBaseFunct3DFromFE3D(feid)
//       ->MakeRefElementData(QuadFormula);
// 
//     N_Joints = cell->GetN_Joints();
//     for(j=0;j<N_Joints;j++)
//     {
//       // OutPut("Joint: " << j << endl);
//       joint = cell->GetJoint(j);
//       JointDOF = JointDOFs[j];
//       neigh = joint->GetNeighbour(cell);
// 
//       if( (neigh == NULL) || (neigh->GetClipBoard() > i) )
//       {
//         // either boundary or right of definition
//         // => generate information
// 
//         // get function values on the joint
//         JointValues=TFEDatabase3D::GetJointValues3D(
//           BaseFuncts[feid], QuadFormula, j);
// 
//         // change signs if necessary
//         TFEDatabase3D::GetBaseFunct3D(BaseFuncts[feid])
//                 ->ChangeBF((TGridCell *)cell, N_Points, JointValues);
// 
//         // calculate int(phi.nE) over the Face
//         memset(PhinE, 0, 3*N_DOFperJoint*SizeOfDouble);
// 
//         for(k=0;k<N_Points;k++)
//         {
//           // values of functions in this quadrature point for all functions
//           JointValue = JointValues[k];
// 
//           switch(reftrans)
//           {
//             case TetraAffin:
//               ((TTetraAffin *)rt)->GetTangentVectors(j, xi[k], eta[k],
//                                                      t11, t12, t13, t21, t22, t23);
//             break;
//       
//             case TetraIsoparametric:
//               ((TTetraIsoparametric *)rt)->GetTangentVectors(j, xi[k], eta[k],
//                                                      t11, t12, t13, t21, t22, t23);
//             break;
//       
//             case HexaAffin:
//               ((THexaAffin *)rt)->GetTangentVectors(j, xi[k], eta[k],
//                                                      t11, t12, t13, t21, t22, t23);
//             break;
//       
//             case HexaTrilinear:
//               ((THexaTrilinear *)rt)->GetTangentVectors(j, xi[k], eta[k],
//                                                      t11, t12, t13, t21, t22, t23);
//             break;
//       
//             case HexaIsoparametric:
//               ((THexaIsoparametric *)rt)->GetTangentVectors(j, xi[k], eta[k],
//                                                      t11, t12, t13, t21, t22, t23);
//             break;
//           } // end switch
// 
//           n1 = t12*t23 - t22*t13;
//           n2 = t21*t13 - t11*t23;
//           n3 = t11*t22 - t12*t21;
// 
//           len = sqrt(n1*n1+n2*n2+n3*n3);
// 
//           // OutPut("cell: " << i << " joint: " << j << " k: " << k);
//           // OutPut(" t1: " << t11 << " " << t12 << " " << t13); 
//           // OutPut(" t2: " << t21 << " " << t22 << " " << t23); 
//           // OutPut(" n: " << n1 << " " << n2 << " " << n3); 
//           // OutPut(" len: " << len << endl);
// 
//           n1 /= len;
//           n2 /= len;
//           n3 /= len;
//           
//           for(l=0;l<N_DOFperJoint;l++)
//           {
//             m = JointDOF[l];
//             s = Weights[k]*len*JointValue[m];
// 
//             // first component has vector (1,0,0)
//             PhinE[l] += s*n1;
//             // second component has vector (0,1,0)
//             PhinE[l+N_DOFperJoint] += s*n2;
//             // third component has vector (0,0,1)
//             PhinE[l+2*N_DOFperJoint] += s*n3;
//           } // endfor l
//         } // endfor k
// 
//         // reset changes in BF
//         TFEDatabase3D::GetBaseFunct3D(BaseFuncts[feid])
//                 ->ChangeBF((TGridCell *)cell, N_Points, JointValues);
// 
//         // find maximum of flux
//         s = -1;
//         for(k=0;k<3*N_DOFperJoint;k++)
//         {
//           if(fabs(PhinE[k]) > s)
//           {
//             l = k;
//             s = fabs(PhinE[k]);
//           } // endif
//         } // endfor k
// 
//         if(s < 1e-10)
//         {
//           Error("Every flux is zero. Something is wrong. " << endl);
//           exit(-1);
//         }
//         // OutPut("Max flux: " << s << " at local joint DOF: " << l << endl);
// 
//         // store jb for this Face
//         m = BeginJb[i]+j;
//         if(l<N_DOFperJoint)
//         {
//           // maximum in first component
//           jb[m] = JointDOF[l];
//           // OutPut("1 Maximum (" << s << ") at: " << JointDOF[l] << endl);
//         }
//         else
//         {
//           if(l<2*N_DOFperJoint)
//           {
//             // maximum in second component
//             jb[m] = JointDOF[l-N_DOFperJoint]+N_LocalDOF;
//             // OutPut("2 Maximum (" << s << ") at: ");
//             // OutPut(JointDOF[l-N_DOFperJoint]+N_LocalDOF << endl);
//           }
//           else
//           {
//             // maximum in third component
//             jb[m] = JointDOF[l-2*N_DOFperJoint]+2*N_LocalDOF;
//             // OutPut("3 Maximum (" << s << ") at: ");
//             // OutPut(JointDOF[l-2*N_DOFperJoint]+2*N_LocalDOF << endl);
//           }
//         }
// 
//         s = PhinE[l];
//         for(k=0;k<3*N_DOFperJoint;k++)
//         {
//           if(k == l)
//             Alpha[3*N_DOFperJoint*m+k] = s;
//           else
//             Alpha[3*N_DOFperJoint*m+k] = PhinE[k]/s;
//         } // endfor k
// 
//       } // end then
//       else
//       {
//         // neighbour cell has already produced data
//         // => get information from where
// 
//         // get number of neighbour
//         NeighNumber = neigh->GetClipBoard();
//         // OutPut("Neigh: " << NeighNumber << endl);
// 
//         // on which local joint of neigh is cell
//         N_JointsNeigh = neigh->GetN_Joints();
//         for(k=0;k<N_JointsNeigh;k++)
//           if(neigh->GetJoint(k)->GetNeighbour(neigh) == cell)
//             break;
//         JointNumberNeigh = k;
//         // OutPut("JointNumberNeigh: " << JointNumberNeigh << endl);
//         
//         // ask neighbour for global number of jb node
//         feidNeigh = FESpace->GetFE3D(NeighNumber, neigh);
//         fedescNeigh = TFEDatabase3D::GetFE3D(feidNeigh)->GetFEDesc3D();
//         N_LocalDOFNeigh = fedescNeigh->GetN_DOF();
//         JointDOFNeigh = fedescNeigh->GetJointDOF(JointNumberNeigh);
//         
//         // get local number in neigh
//         l = jb[BeginJb[NeighNumber]+JointNumberNeigh];
//         // OutPut("neigh-local-jb: " << l << endl);
// 
//         DOFNeigh = GlobalNumbers + BeginIndex[NeighNumber];
// 
//         m = BeginJb[i]+j;
//         if(l<N_LocalDOFNeigh)
//         {
//           // jb node in first component
//           JbNumber = DOFNeigh[l];
//           for(k=0;k<N_DOFperJoint;k++)
//             if(DOF[JointDOF[k]] == JbNumber) break;
//           if(k==N_DOFperJoint)
//           {
//             Error(i << " " << j << " nothing found 1" << endl);
//           }
//           jb[m] = JointDOF[k];
//         }
//         else
//         {
//           if(l<2*N_LocalDOFNeigh)
//           {
//             // jb node in second component
//             JbNumber = DOFNeigh[l-N_LocalDOFNeigh];
//             for(k=0;k<N_DOFperJoint;k++)
//               if(DOF[JointDOF[k]] == JbNumber) break;
//             if(k==N_DOFperJoint)
//             {
//               Error("nothing found 2" << endl);
//             }
//             JbNumber += N_TotalDOF;
//             jb[m] = JointDOF[k]+N_LocalDOF;
//           }
//           else
//           {
//             // jb node in third component
//             JbNumber = DOFNeigh[l-2*N_LocalDOFNeigh];
//             for(k=0;k<N_DOFperJoint;k++)
//               if(DOF[JointDOF[k]] == JbNumber) break;
//             if(k==N_DOFperJoint)
//             {
//               Error("nothing found 3" << endl);
//             }
//             JbNumber += 2*N_TotalDOF;
//             jb[m] = JointDOF[k]+2*N_LocalDOF;
//           }
//         }
//         // OutPut("global-jb: " << JbNumber << endl);
//         // OutPut("local-jb: " << jb[m] << endl);
//         
//         for(k=0;k<N_DOFperJoint;k++)
//         {
//           n = DOFNeigh[JointDOFNeigh[k]];
//           // OutPut("k: " << k << " " << n << endl);
// 
//           for(l=0;l<N_DOFperJoint;l++)
//           {
//             if( DOF[JointDOF[l]] == n )
//             {
//               // OutPut("l :" << l << endl);
//               break;
//             } // endif
//           } // endfor l
// 
//           // if jb node is found change sign due to different direction of normals
//           if(n != JbNumber)
//             Alpha[3*N_DOFperJoint*m+l] = Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]
//                                                +JointNumberNeigh)+k];
//           else
//           {
//             Alpha[3*N_DOFperJoint*m+l] = -Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]
//                                                 +JointNumberNeigh)+k];
//             // OutPut("Change in first component" << endl);
//           }
// 
//           if(n+N_TotalDOF != JbNumber)
//           {
//             Alpha[3*N_DOFperJoint*m+l+N_DOFperJoint] =
//             Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]+JointNumberNeigh)+k+N_DOFperJoint];
//           }
//           else
//           {
//             Alpha[3*N_DOFperJoint*m+l+N_DOFperJoint] =
//             -Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]+JointNumberNeigh)+k+N_DOFperJoint];
//             // OutPut("Change in second component" << endl);
//           }
// 
//           if(n+2*N_TotalDOF != JbNumber)
//           {
//             Alpha[3*N_DOFperJoint*m+l+2*N_DOFperJoint] =
//             Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]+JointNumberNeigh)+k+2*N_DOFperJoint];
//           }
//           else
//           {
//             Alpha[3*N_DOFperJoint*m+l+2*N_DOFperJoint] =
//             -Alpha[3*N_DOFperJoint*(BeginJb[NeighNumber]+JointNumberNeigh)+k+2*N_DOFperJoint];
//             // OutPut("Change in third component" << endl);
//           }
//         } // endfor k
// 
//       } // endelse
//     } // endfor j
//   } // endfor i
// 
//   /*
//   for(i=0;i<N_Cells;i++)
//   {
//     OutPut("cell: " << i << endl);
//     for(j=BeginJb[i];j<BeginJb[i+1];j++)
//     {
//       OutPut(setw(3) << jb[j]); 
//     }
//     OutPut(endl);
//   }
// 
//   OutPut(endl);
//   for(i=0;i<N_Cells;i++)
//   {
//     OutPut(endl << "cell: " << i << endl);
//     for(j=BeginJb[i];j<BeginJb[i+1];j++)
//     {
//       OutPut("joint: " << j-BeginJb[i] << endl);
//       for(k=0;k<3*N_DOFperJoint;k++)
//       {
//         OutPut(setw(12) << Alpha[3*N_DOFperJoint*j+k]);
//       }
//       OutPut(endl);
//     }
//   }
//   */
// } // end GenerateAlpha
