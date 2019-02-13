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
   
// #include <Constants.h>
// #include <Enumerations.h>
// #include <Database.h>
// #include <FEVectFunct3D.h>
// #include <FEDatabase3D.h>
// #include <TetraAffin.h>
// 
// #include <MooNMD_Io.h>
// 
// #include <stdlib.h>
// 
// void ESGridAssemble(double Mult, double *coeff, 
// 		    double *param, double hK, 
// 		    double **OrigValues, int *N_BaseFuncts,
// 		    double ***LocMatrices, double **LocRhs)
// {
//   double **MatrixA11, **MatrixA22, **MatrixA33, **MatrixA12, **MatrixA13;
//   double **MatrixA21, **MatrixA23, **MatrixA31, **MatrixA32;
//   double *Orig0, *Orig1, *Orig2, c0;
//   double test100, test010, test001, ansatz100, ansatz010, ansatz001;
//   double val, val1;
//   double *A11Row, *A12Row, *A13Row, *A21Row, *A22Row, *A23Row, *A31Row;
//   double *A32Row, *A33Row;
//   
//   int N_;
//   
//   MatrixA11 = LocMatrices[0];
//   MatrixA12 = LocMatrices[1];
//   MatrixA13 = LocMatrices[2];
//   MatrixA21 = LocMatrices[3];
//   MatrixA22 = LocMatrices[4];
//   MatrixA23 = LocMatrices[5];
//   MatrixA31 = LocMatrices[6];
//   MatrixA32 = LocMatrices[7];
//   MatrixA33 = LocMatrices[8];
//   
//   Orig0 = OrigValues[0]; // u_x
//   Orig1 = OrigValues[1]; // u_y
//   Orig2 = OrigValues[2]; // u_z
//   
//   c0 = coeff[0]; // lame constant
//   
//   N_ = N_BaseFuncts[0];
//   
//   for (int i=0;i<N_;++i)
//   {
//     test100 = Orig0[i];
//     test010 = Orig1[i];
//     test001 = Orig2[i];
//     
//     A11Row = MatrixA11[i];
//     A12Row = MatrixA12[i];
//     A13Row = MatrixA13[i];
//     A21Row = MatrixA21[i];
//     A22Row = MatrixA22[i];
//     A23Row = MatrixA23[i];
//     A31Row = MatrixA31[i];
//     A32Row = MatrixA32[i];
//     A33Row = MatrixA33[i];
//     
//     for (int j=0;j<N_;++j)
//     {
//       ansatz100 = Orig0[j];
//       ansatz010 = Orig1[j];
//       ansatz001 = Orig2[j];
//       
// //       val1 =   0.5 *(ansatz100*test100+ansatz010*test010+ansatz001*test001);
// //       
// //       val  = val1 + 0.75*ansatz100*test100;
// //       A11Row[j] +=   Mult*val;
// //       
// //       val  = val1 + 0.75*ansatz010*test010;
// //       A22Row[j] += Mult*val;
// //       
// //       val  = val1 + 0.75*ansatz001*test001;
// //       A33Row[j] += Mult*val;
// //       
// //       val = 2*ansatz010*test100 + ansatz100*test010;
// //       A12Row[j] += Mult*val;
// //       
// //       val = 2*ansatz001*test100 + ansatz100*test001;
// //       A13Row[j] += Mult*val;
// //       
// //       val = 2*ansatz100*test010 + ansatz010*test100;
// //       A21Row[j] += Mult*val;
// //       
// //       val = 2*ansatz001*test010 + ansatz010*test001;
// //       A23Row[j] += Mult*val;
// //       
// //       val = 2*ansatz100*test001 + ansatz001*test100;
// //       A31Row[j] += Mult*val;
// //       
// //       val = 2*ansatz010*test001 + ansatz001*test010;
// //       A32Row[j] += Mult*val;
// 
// 	val = ansatz100*test100 + ansatz010*test010 + ansatz001*test001;
// 	A11Row[j] += Mult*val;
// 	A22Row[j] += Mult*val;
// 	A33Row[j] += Mult*val;
//     }
//   }
// }
// 
// void GridBoundCond (double x, double y, double z, BoundCond &cond)
// {
//   cond = DIRICHLET;
// }
// 
// void GridXBoundValues (double x, double y, double z, double &value)
// {
//   value = 0;
// }
// 
// void GridYBoundValues (double x, double y, double z, double &value)
// {
//   value = 0;
// }
// 
// void GridZBoundValues (double x, double y, double z, double &value)
// {
//   value = 0;
// }
// 
// void GridCoeffs(int n_points, double *X, double *Y, double *Z,
// 		double **parameters, double **coeffs)
// {
//   int i;
//   double *coeff;
// 
//   for (i=0;i<n_points;++i)
//   {
//     coeff = coeffs[i];
//     
//     coeff[0] = 1.0;
//     coeff[1] = 0; 
//     coeff[2] = 0;
//     coeff[3] = 0;
//   }    
// }
// 
// void MovingTimeNSParamsVelo3D(double *in, double *out)
// {
//   if ( TDatabase::ParamDB->PROBLEM_TYPE == 3 )
//   {
//     out[0] = 0;
//     out[1] = 0;
//     out[2] = 0;    
//   }
//   else
//   {
//     out[0] = in[3] - in[6];
//     out[1] = in[4] - in[7];
//     out[2] = in[5] - in[8];
//   }
//   
// }
// 
// void MovingTimeNSParamsVelo3D_TwoPhase(double *in, double *out)
// {
//   int phase = (int) out[0];
//   double ratio = TDatabase::ParamDB->P6;
//   
// //   cout << phase << endl;
//   
//   if ( TDatabase::ParamDB->PROBLEM_TYPE == 3 )
//   {
//     out[0] = 0;
//     out[1] = 0;
//     out[2] = 0;    
//   }
//   else
//   {
//     out[0] = in[3] - in[6];
//     out[1] = in[4] - in[7];
//     out[2] = in[5] - in[8];
//   }
//   
//   if ( phase == 1 )
//   {
//     out[0] *= ratio;
//     out[1] *= ratio;
//     out[2] *= ratio;
//   }
//   
// }
// 
// void SetGridRhs(TFEVectFunct3D *velocity, TFESpace3D *fesp_grid,
// 		int N_BoundFaces, int *CellNumbers, int *JointNumbers,
// 		double dt,
// 		double *rhs1, double *rhs2, double *rhs3)
// {
//   int CellNr, JointNr, *JointDOF_grid, N_DOF_velo;
//   int *N_BaseFuncts, N_BaseFunct, *DOF_velo, *DOF_grid;
//   int *GridGlobalNumbers, *GridBeginIndex;
//   int *VeloGlobalNumbers, *VeloBeginIndex;
//   const int *TmpFV, *TmpLen;
//   int MaxLen;
//   
//   TFESpace3D *fesp_velo;
//   TCollection *Coll;
//   TBaseCell *Cell;
//   BaseFunct3D *BaseFuncts;
//   FE3D FeID_velo, FeID_grid;
//   TBaseFunct3D *BaseFunct_velo;
//   BF3DRefElements RefElement;
//   TFEDesc3D *FEDesc_velo, *FEDesc_grid;
//   TShapeDesc *ShapeDesc;
//   
//   double val1x, val1y, val1z;
//   double val2x, val2y, val2z;
//   double val3x, val3y, val3z, *Values;
//   double n1x, n2x, n3x;
//   double n1y, n2y, n3y;
//   double n1z, n2z, n3z, dot;
//   double uref1[MaxN_BaseFunctions3D], uref2[MaxN_BaseFunctions3D];
//   double uref3[MaxN_BaseFunctions3D];
//   double xi[3], eta[3], zeta[3];
//   double *Coords;
//   double Coords_tetra[][3] = { { 0.0, 0.0, 0.0 },
// 			       { 1.0, 0.0, 0.0 },
// 			       { 0.0, 1.0, 0.0 },
// 			       { 0.0, 0.0, 1.0 } };
//   
//   fesp_velo = velocity->GetFESpace3D();
//   
//   Coll = fesp_grid->GetCollection();
//   
//   GridGlobalNumbers = fesp_grid->GetGlobalNumbers();
//   GridBeginIndex    = fesp_grid->GetBeginIndex();
//   VeloGlobalNumbers = fesp_velo->GetGlobalNumbers();
//   VeloBeginIndex    = fesp_velo->GetBeginIndex();
//   
//   Values = velocity->GetValues();
//   N_DOF_velo = velocity->GetLength();
//   
//   for (int i=0;i<N_BoundFaces;++i)
//   {
//     CellNr = CellNumbers[i];
//     JointNr = JointNumbers[i];
//     
//     Cell = Coll->GetCell(CellNr);
//     ShapeDesc = Cell->GetShapeDesc();
//     ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
//     
//     FeID_velo = fesp_velo->GetFE3D(CellNr, Cell);
//     FeID_grid = fesp_grid->GetFE3D(CellNr, Cell);
//     
//     BaseFunct_velo = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID_velo);
//     FEDesc_grid    = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID_grid);
//     
//     JointDOF_grid = FEDesc_grid->GetJointDOF(JointNr);
//     
//     N_BaseFunct = BaseFunct_velo->GetDimension();
//     
//     RefElement = BaseFunct_velo->GetRefElement();
//     
//     DOF_velo = VeloGlobalNumbers + VeloBeginIndex[CellNr];
//     DOF_grid = GridGlobalNumbers + GridBeginIndex[CellNr];
//     
//     switch (RefElement)
//     {
//       case BFUnitTetrahedron:
// 	xi[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][0];
// 	xi[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][0];
// 	xi[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][0];
// 	
// 	eta[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][1];
// 	eta[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][1];
// 	eta[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][1];
// 	
// 	zeta[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][2];
// 	zeta[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][2];
// 	zeta[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][2];
// 	
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr  ])->SetResetNormal();
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+1])->SetResetNormal();
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+2])->SetResetNormal();
// 	
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr  ])->GetNormal(n1x, n1y, n1z);
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+1])->GetNormal(n2x, n2y, n2z);
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+2])->GetNormal(n3x, n3y, n3z);
// 	break;
// 	
//       case BFUnitHexahedron:
// 	cerr << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
// 	exit(0);
// 	break;
//     }
//     
//     BaseFunct_velo->GetDerivatives(D000, xi[0], eta[0], zeta[0], uref1);
//     BaseFunct_velo->GetDerivatives(D000, xi[1], eta[1], zeta[1], uref2);
//     BaseFunct_velo->GetDerivatives(D000, xi[2], eta[2], zeta[2], uref3);
//  
//     val1x = 0; val2x = 0; val3x = 0;
//     val1y = 0; val2y = 0; val3y = 0;
//     val1z = 0; val2z = 0; val3z = 0;
//     for (int j=0;j<N_BaseFunct;++j)
//     {
//       val1x += Values[             DOF_velo[j]] * uref1[j];
//       val1y += Values[  N_DOF_velo+DOF_velo[j]] * uref1[j];
//       val1z += Values[2*N_DOF_velo+DOF_velo[j]] * uref1[j];
//       
//       val2x += Values[             DOF_velo[j]] * uref2[j];
//       val2y += Values[  N_DOF_velo+DOF_velo[j]] * uref2[j];
//       val2z += Values[2*N_DOF_velo+DOF_velo[j]] * uref2[j];
//       
//       val3x += Values[             DOF_velo[j]] * uref3[j];
//       val3y += Values[  N_DOF_velo+DOF_velo[j]] * uref3[j];
//       val3z += Values[2*N_DOF_velo+DOF_velo[j]] * uref3[j];
//     }
//     
//     if ( (int) TDatabase::ParamDB->P0 )
//     {
//       dot = dt * (val1x*n1x + val1y*n1y + val1z*n1z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[0]]] = dot * n1x;
//       rhs2[DOF_grid[JointDOF_grid[0]]] = dot * n1y;
//       rhs3[DOF_grid[JointDOF_grid[0]]] = dot * n1z;
//       
//       dot = dt * (val2x*n2x + val2y*n2y + val2z*n2z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[1]]] = dot * n2x;
//       rhs2[DOF_grid[JointDOF_grid[1]]] = dot * n2y;
//       rhs3[DOF_grid[JointDOF_grid[1]]] = dot * n2z;
//       
//       dot = dt * (val3x*n3x + val3y*n3y + val3z*n3z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[2]]] = dot * n3x;
//       rhs2[DOF_grid[JointDOF_grid[2]]] = dot * n3y;
//       rhs3[DOF_grid[JointDOF_grid[2]]] = dot * n3z;
//     }
//     else 
//     {
//       rhs1[DOF_grid[JointDOF_grid[0]]] = dt * val1x;
//       rhs2[DOF_grid[JointDOF_grid[0]]] = dt * val1y;
//       rhs3[DOF_grid[JointDOF_grid[0]]] = dt * val1z;
// 
//       rhs1[DOF_grid[JointDOF_grid[1]]] = dt * val2x;
//       rhs2[DOF_grid[JointDOF_grid[1]]] = dt * val2y;
//       rhs3[DOF_grid[JointDOF_grid[1]]] = dt * val2z;
//       
//       rhs1[DOF_grid[JointDOF_grid[2]]] = dt * val3x;
//       rhs2[DOF_grid[JointDOF_grid[2]]] = dt * val3y;
//       rhs3[DOF_grid[JointDOF_grid[2]]] = dt * val3z;
//     }
//       
//   }
// }
// 
// void SetGridRhs(TFEVectFunct3D *velocity, TFESpace3D *fesp_grid,
// 		int N_BoundFaces, int *CellNumbers, int *JointNumbers,
// 		int *GlobalCellNo, double dt,
// 		double *rhs1, double *rhs2, double *rhs3)
// {
//   int CellNr, CellNrLoc, JointNr, *JointDOF_grid, N_DOF_velo;
//   int *N_BaseFuncts, N_BaseFunct, *DOF_velo, *DOF_grid;
//   int *GridGlobalNumbers, *GridBeginIndex;
//   int *VeloGlobalNumbers, *VeloBeginIndex;
//   const int *TmpFV, *TmpLen;
//   int MaxLen;
//   
//   TFESpace3D *fesp_velo;
//   TCollection *Coll;
//   TBaseCell *Cell;
//   BaseFunct3D *BaseFuncts;
//   FE3D FeID_velo, FeID_grid;
//   TBaseFunct3D *BaseFunct_velo;
//   BF3DRefElements RefElement;
//   TFEDesc3D *FEDesc_velo, *FEDesc_grid;
//   TShapeDesc *ShapeDesc;
//   
//   double val1x, val1y, val1z;
//   double val2x, val2y, val2z;
//   double val3x, val3y, val3z, *Values;
//   double n1x, n2x, n3x;
//   double n1y, n2y, n3y;
//   double n1z, n2z, n3z, dot;
//   double uref1[MaxN_BaseFunctions3D], uref2[MaxN_BaseFunctions3D];
//   double uref3[MaxN_BaseFunctions3D];
//   double xi[3], eta[3], zeta[3];
//   double *Coords;
//   double Coords_tetra[][3] = { { 0.0, 0.0, 0.0 },
// 			       { 1.0, 0.0, 0.0 },
// 			       { 0.0, 1.0, 0.0 },
// 			       { 0.0, 0.0, 1.0 } };
//   
//   fesp_velo = velocity->GetFESpace3D();
//   
//   Coll = fesp_velo->GetCollection();
//   
//   GridGlobalNumbers = fesp_grid->GetGlobalNumbers();
//   GridBeginIndex    = fesp_grid->GetBeginIndex();
//   VeloGlobalNumbers = fesp_velo->GetGlobalNumbers();
//   VeloBeginIndex    = fesp_velo->GetBeginIndex();
//   
//   Values = velocity->GetValues();
//   N_DOF_velo = velocity->GetLength();
//   
//   for (int i=0;i<N_BoundFaces;++i)
//   {
//     CellNrLoc = CellNumbers[i];
//     CellNr = GlobalCellNo[CellNrLoc];
//     JointNr = JointNumbers[i];
//     
//     Cell = Coll->GetCell(CellNr);
//     ShapeDesc = Cell->GetShapeDesc();
//     ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
//     
//     FeID_velo = fesp_velo->GetFE3D(CellNr, Cell);
//     FeID_grid = fesp_grid->GetFE3D(CellNrLoc, Cell);
//     
//     BaseFunct_velo = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID_velo);
//     FEDesc_grid    = TFEDatabase3D::GetFEDesc3DFromFE3D(FeID_grid);
//     
//     JointDOF_grid = FEDesc_grid->GetJointDOF(JointNr);
//     
//     N_BaseFunct = BaseFunct_velo->GetDimension();
//     
//     RefElement = BaseFunct_velo->GetRefElement();
//     
//     DOF_velo = VeloGlobalNumbers + VeloBeginIndex[CellNr];
//     DOF_grid = GridGlobalNumbers + GridBeginIndex[CellNrLoc];
//     
//     switch (RefElement)
//     {
//       case BFUnitTetrahedron:
// 	xi[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][0];
// 	xi[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][0];
// 	xi[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][0];
// 	
// 	eta[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][1];
// 	eta[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][1];
// 	eta[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][1];
// 	
// 	zeta[0] = Coords_tetra[TmpFV[MaxLen*JointNr  ]][2];
// 	zeta[1] = Coords_tetra[TmpFV[MaxLen*JointNr+1]][2];
// 	zeta[2] = Coords_tetra[TmpFV[MaxLen*JointNr+2]][2];
// 	
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr  ])->SetResetNormal();
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+1])->SetResetNormal();
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+2])->SetResetNormal();
// 	
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr  ])->GetNormal(n1x, n1y, n1z);
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+1])->GetNormal(n2x, n2y, n2z);
// 	Cell->GetVertex(TmpFV[MaxLen*JointNr+2])->GetNormal(n3x, n3y, n3z);
// 	break;
// 	
//       case BFUnitHexahedron:
// 	cerr << __FILE__ << ":" << __LINE__ << "Nothing implemented for BFUnitHexahedron. Exit!" << endl;
// 	exit(0);
// 	break;
//     }
//     
//     BaseFunct_velo->GetDerivatives(D000, xi[0], eta[0], zeta[0], uref1);
//     BaseFunct_velo->GetDerivatives(D000, xi[1], eta[1], zeta[1], uref2);
//     BaseFunct_velo->GetDerivatives(D000, xi[2], eta[2], zeta[2], uref3);
//  
//     val1x = 0; val2x = 0; val3x = 0;
//     val1y = 0; val2y = 0; val3y = 0;
//     val1z = 0; val2z = 0; val3z = 0;
//     for (int j=0;j<N_BaseFunct;++j)
//     {
//       val1x += Values[             DOF_velo[j]] * uref1[j];
//       val1y += Values[  N_DOF_velo+DOF_velo[j]] * uref1[j];
//       val1z += Values[2*N_DOF_velo+DOF_velo[j]] * uref1[j];
//       
//       val2x += Values[             DOF_velo[j]] * uref2[j];
//       val2y += Values[  N_DOF_velo+DOF_velo[j]] * uref2[j];
//       val2z += Values[2*N_DOF_velo+DOF_velo[j]] * uref2[j];
//       
//       val3x += Values[             DOF_velo[j]] * uref3[j];
//       val3y += Values[  N_DOF_velo+DOF_velo[j]] * uref3[j];
//       val3z += Values[2*N_DOF_velo+DOF_velo[j]] * uref3[j];
//     }
//     
//     if ( (int) TDatabase::ParamDB->P0 )
//     {
//       dot = dt * (val1x*n1x + val1y*n1y + val1z*n1z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[0]]] = dot * n1x;
//       rhs2[DOF_grid[JointDOF_grid[0]]] = dot * n1y;
//       rhs3[DOF_grid[JointDOF_grid[0]]] = dot * n1z;
//       
//       dot = dt * (val2x*n2x + val2y*n2y + val2z*n2z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[1]]] = dot * n2x;
//       rhs2[DOF_grid[JointDOF_grid[1]]] = dot * n2y;
//       rhs3[DOF_grid[JointDOF_grid[1]]] = dot * n2z;
//       
//       dot = dt * (val3x*n3x + val3y*n3y + val3z*n3z);
//   //     dot = 1.0;
//       rhs1[DOF_grid[JointDOF_grid[2]]] = dot * n3x;
//       rhs2[DOF_grid[JointDOF_grid[2]]] = dot * n3y;
//       rhs3[DOF_grid[JointDOF_grid[2]]] = dot * n3z;
//     }
//     else 
//     {
//       rhs1[DOF_grid[JointDOF_grid[0]]] = dt * val1x;
//       rhs2[DOF_grid[JointDOF_grid[0]]] = dt * val1y;
//       rhs3[DOF_grid[JointDOF_grid[0]]] = dt * val1z;
// 
//       rhs1[DOF_grid[JointDOF_grid[1]]] = dt * val2x;
//       rhs2[DOF_grid[JointDOF_grid[1]]] = dt * val2y;
//       rhs3[DOF_grid[JointDOF_grid[1]]] = dt * val2z;
//       
//       rhs1[DOF_grid[JointDOF_grid[2]]] = dt * val3x;
//       rhs2[DOF_grid[JointDOF_grid[2]]] = dt * val3y;
//       rhs3[DOF_grid[JointDOF_grid[2]]] = dt * val3z;
//     }
//       
//   }
// }
// 
// void MapGridVelo(TFESpace3D *grid_space_P, TFEVectFunct3D *g,
// 		  int *GlobalCellNo_P, double *grid_sol_p)
// {
//   int N_Cells_P, N_DOF_P, N_DOF;
//   int *GlobalNumbers_P, *GlobalNumbers;
//   int *BeginIndex_P, *BeginIndex;
//   int *DOF_P, *DOF;
//   TCollection *Coll_P;
//   TFESpace3D *fespace;
//   double *Values;
//   
//   fespace = g->GetFESpace3D();
//   Coll_P = grid_space_P->GetCollection();
//   
//   N_Cells_P = Coll_P->GetN_Cells();
//   N_DOF_P = grid_space_P->GetN_DegreesOfFreedom();
//   N_DOF   = fespace->GetN_DegreesOfFreedom();
//   
//   fespace = g->GetFESpace3D();
//   
//   GlobalNumbers = fespace->GetGlobalNumbers();
//   BeginIndex    = fespace->GetBeginIndex();
//   
//   GlobalNumbers_P = grid_space_P->GetGlobalNumbers();
//   BeginIndex_P    = grid_space_P->GetBeginIndex();
//   
//   Values = g->GetValues();
//   
//   for (int i=0;i<N_Cells_P;++i)
//   {
//     DOF_P = GlobalNumbers_P + BeginIndex_P[i];
//     DOF   = GlobalNumbers   + BeginIndex[GlobalCellNo_P[i]];
//     
//     for (int j=0;j<4;++j) // P1
//     {
//       Values[DOF[j]        ] = grid_sol_p[DOF_P[j]          ];
//       Values[DOF[j]+  N_DOF] = grid_sol_p[DOF_P[j]+  N_DOF_P];
//       Values[DOF[j]+2*N_DOF] = grid_sol_p[DOF_P[j]+2*N_DOF_P];
//     }
//   }
// }
// 
// void SetNormals(TCollection *Coll_P, int N_SurfaceJoints, int *CellNumbers, int *JointNumbers)
// {
//   int CellNr, JointNr;
//   const int *TmpFV, *TmpLen;
//   int MaxLen;
//   double n1, n2, n3;
//   TBaseCell *Cell;
//   TTetraAffin *F_aff;
//   TVertex *Vertex;
//   
// //   // Reset normals
// //   for (int i=0;i<N_SurfaceJoints;++i)
// //   {
// //     CellNr = CellNumbers[i];
// //     JointNr = JointNumbers[i];
// //     
// //     Cell = Coll_P->GetCell(CellNr);
// //     
// //     Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
// //     
// //     for (int j=0;j<TmpLen[JointNr];++j)
// //     {
// //       Vertex = Cell->GetVertex(TmpFV[JointNr*MaxLen+j]);
// //       Vertex->SetResetNormal();
// //     }
// //   }
//   
//   for (int i=0;i<N_SurfaceJoints;++i)
//   {
//     CellNr = CellNumbers[i];
//     JointNr = JointNumbers[i];
//     
//     Cell = Coll_P->GetCell(CellNr);
//     
//     Cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
//     
//     F_aff = (TTetraAffin*) TFEDatabase3D::GetRefTrans3D(TetraAffin);
//     F_aff->SetCell(Cell);
//     
//     F_aff->GetOuterNormal(JointNr, 0., 0., n1, n2, n3);
//     
//     for (int j=0;j<TmpLen[JointNr];++j)
//     {
//        Vertex = Cell->GetVertex(TmpFV[JointNr*MaxLen+j]);
//        Vertex->AddNormal(n1, n2, n3);
//     }
//   }
// }
