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
   
// #include <Database.h>
// #include <MainUtilities.h> // linfb, ave_l2b_quad_points
// #include <FEDatabase2D.h>
// #include <FEFunction2D.h>
// #include <LocalAssembling2D.h>
// #include <ConvDiff.h>
// #include <ConvDiff2D.h> // local assembling routines for 2D convection-diffusion
// #include <Darcy2D.h> // local assembling routines for 2D Darcy problems
// #include <NSE2D_FixPo.h>// local assembling routines for 2D Navier-Stokes
// #include <NSE2D_FixPoSkew.h>// local assembling routines for 2D Navier-Stokes
// #include <NSE2D_FixPoRot.h>// local assembling routines for 2D Navier-Stokes
// #include <MooNMD_Io.h>
// #include <string.h>
// 
// /** @brief a helper function returning a string with for the name of the 
//  *         LocalAssembling2D_type. This returns an empty string in case the type
//  *         is not known. */
// std::string LocalAssembling2D_type_to_string(LocalAssembling2D_type type)
// {
//   switch(type)
//   {
//     ///////////////////////////////////////////////////////////////////////////
//     // CD2D: stationary convection diffusion problems
//     case CD2D_Galerkin:
//       return std::string("CD2D_Galerkin");
//     case CD2D_SUPG:
//       return std::string("CD2D_SUPG");
//     case CD2D_GLS:
//       return std::string("CD2D_GLS");
//     case CD2D_Axiax3D_Galerkin:
//       return std::string("CD2D_Axiax3D_Galerkin");
//     ///////////////////////////////////////////////////////////////////////////
//     // TCD2D: time dependent convection diffusion problems
//     case TCD2D_Mass_Rhs_Galerkin:
//       return std::string("TCD2D_Mass_Rhs_Galerkin");
//     case TCD2D_Stiff_Rhs_Galerkin:
//       return std::string("TCD2D_Stiff_Rhs_Galerkin");
//     case TCD2D_Mass_Rhs_SUPG:
//       return std::string("TCD2D_Mass_Rhs_SUPG");
//     case TCD2D_Stiff_Rhs_SUPG:
//       return std::string("TCD2D_Stiff_Rhs_SUPG");
//     ///////////////////////////////////////////////////////////////////////////
//     // NSE2D: stationary Navier-Stokes problems
//     case NSE2D_Galerkin:
//       return std::string("NSE2D_Galerkin");
//     case NSE2D_Galerkin_Nonlinear:
//       return std::string("NSE2D_Galerkin_Nonlinear");
//     ///////////////////////////////////////////////////////////////////////////
//     // Darcy2D: stationary Darcy problems
//     case Darcy2D_Galerkin:
//       return std::string("Darcy2D_Galerkin");
//     default: return std::string();
//   }
// }
// 
// LocalAssembling2D::LocalAssembling2D(LocalAssembling2D_type type, 
//                                      TFEFunction2D **fefunctions2d,
//                                      CoeffFct2D *coeffs)
//  : FEFunctions2D(fefunctions2d), Coeffs(coeffs)
// {
//   this->name = LocalAssembling2D_type_to_string(type);
//   if(TDatabase::ParamDB->SC_VERBOSE > 1)
//     OutPut("Constructor of LocalAssembling2D: using type " << name << endl);
//   
//   
//   
// //   !!!! Error in Mac Compiler - Sashikumaar !!!!!!!
//   // the values below only matter if you need an existing finite element 
//   // function during your assembly. Change them in such a case
// //   this->N_Parameters = 0;
// //   this->N_ParamFct = 0;
// //   this->ParameterFct = {};
// //   this->N_FEValues = 0;
// //   this->FEValue_FctIndex = {};
// //   this->BeginParameter = {};
// //   
// //   // set all member variables according to the LocalAssembling2D_type
// //   switch(type)
// //   {
// //     ///////////////////////////////////////////////////////////////////////////
// //     // CD2D: stationary convection diffusion problems
// //     case CD2D_Galerkin:
// //       this->N_Terms = 3;
// //       this->Derivatives = { D10, D01, D00 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = false;
// //       this->FESpaceNumber = { 0, 0, 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = BilinearAssembleGalerkin; 
// //       this->Manipulate = NULL;
// //       break;
// //     case CD2D_SUPG:
// //       this->N_Terms = 5;
// //       this->Derivatives = { D10, D01, D00, D20, D02 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = true;
// //       this->FESpaceNumber = { 0, 0, 0, 0, 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = BilinearAssemble_SD; 
// //       if(TDatabase::ParamDB->SDFEM_NORM_B==0)
// //         this->Manipulate = linfb;
// //       else
// //         this->Manipulate = ave_l2b_quad_points;
// //       break;
// //     case CD2D_GLS:
// //       this->N_Terms = 5;
// //       this->Derivatives = { D10, D01, D00, D20, D02 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = true;
// //       this->FESpaceNumber = { 0, 0, 0, 0, 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = BilinearAssemble_GLS; 
// //       if(TDatabase::ParamDB->SDFEM_NORM_B==0)
// //         this->Manipulate = linfb;
// //       else
// //         this->Manipulate = ave_l2b_quad_points;
// //       break;
// //     case CD2D_Axiax3D_Galerkin:
// //       this->N_Terms = 3;
// //       this->Derivatives = { D10, D01, D00 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = false;
// //       this->FESpaceNumber = { 0, 0, 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = BilinearAssemble_Axial3D; 
// //       this->Manipulate = NULL;
// //       break;
// //     ///////////////////////////////////////////////////////////////////////////
// //     // TCD2D: time dependent convection diffusion problems
// //     case TCD2D_Mass_Rhs_Galerkin:
// //       this->N_Terms = 1;
// //       this->Derivatives = { D00 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = false;
// //       this->FESpaceNumber = { 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = MatrixMRhsAssemble; 
// //       this->Manipulate = NULL;
// //       break;
// //     case TCD2D_Stiff_Rhs_Galerkin:
// //       this->N_Terms = 3;
// //       this->Derivatives = { D10, D01, D00 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = false;
// //       this->FESpaceNumber = { 0, 0, 0 };
// //       this->N_Matrices = 1;
// //       this->RowSpace = { 0 };
// //       this->ColumnSpace = { 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = MatrixARhsAssemble; 
// //       this->Manipulate = NULL;
// //       break;
// //     case TCD2D_Mass_Rhs_SUPG:
// //       this->N_Terms = 3;
// //       this->Derivatives = { D10, D01, D00 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[0] = false;
// //       this->FESpaceNumber = { 0, 0, 0 };
// //       this->N_Matrices = 2;
// //       this->RowSpace = { 0, 0 };
// //       this->ColumnSpace = { 0, 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = MatrixMRhsAssemble_SUPG; 
// //       this->Manipulate = NULL;
// //       break;
// //     case TCD2D_Stiff_Rhs_SUPG:
// //       this->N_Terms = 5;
// //       this->Derivatives = { D10, D01, D00, D20, D02 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[1] = true;
// //       this->FESpaceNumber = { 0, 0, 0, 0, 0 }; // number of terms = 5
// //       this->N_Matrices = 2;
// //       this->RowSpace = { 0, 0 };
// //       this->ColumnSpace = { 0, 0 };
// //       this->N_Rhs = 1;
// //       this->RhsSpace = { 0 };
// //       this->AssembleParam = MatricesAKRhsAssemble_SUPG; 
// //       this->Manipulate = NULL;
// //       break;
// //     ///////////////////////////////////////////////////////////////////////////
// //     // NSE2D: stationary Navier-Stokes problems
// //     case NSE2D_Galerkin:
// //     case NSE2D_Galerkin_Nonlinear:
// //       this->set_parameters_for_nse(type);
// //       break;
// //     case Darcy2D_Galerkin:
// //       this->N_Terms = 6;
// //       this->Derivatives = { D00, D00, D10, D01, D10, D01 };
// //       this->Needs2ndDerivatives = new bool[1];
// //       this->Needs2ndDerivatives[1] = false;
// //       this->FESpaceNumber = { 0, 1, 0, 0, 1, 1};
// //       this->N_Matrices = 4;
// //       this->RowSpace = {0, 1, 0, 1};
// //       this->ColumnSpace = { 0, 1, 1, 0};
// //       this->N_Rhs = 2;
// //       this->RhsSpace = { 0, 1 };
// //       this->AssembleParam = BilinearAssembleDarcyGalerkin; 
// //       this->Manipulate = NULL;
// //       break;
// //     default:
// //       ErrMsg("unknown LocalAssembling2D_type " << type << " " << this->name);
// //       throw("unknown LocalAssembling2D_type");
// //   }
// //   
// //   AllOrigValues = new double** [N_Terms];
// //   OrigValues = new double* [N_Terms];
// //   
//   // some consistency checks
//   if(Coeffs == NULL)
//   {
//     ErrMsg("You need to specify a valid function for the coefficients");
//     exit(1);
//   }
//   if(AssembleParam == NULL)
//   {
//     ErrMsg("a local assmebling routine was not set");
//     exit(1);
//   }
// }
// 
// LocalAssembling2D::~LocalAssembling2D()
// {
//   delete [] AllOrigValues;
//   delete [] OrigValues;
//   delete [] Needs2ndDerivatives;
// }
// 
// 
// void LocalAssembling2D::GetLocalForms(int N_Points, double *weights, 
//                                       double *AbsDetjk, double *X, double *Y,
//                                       int *N_BaseFuncts,
//                                       BaseFunct2D *BaseFuncts, 
//                                       double **Parameters, double **AuxArray,
//                                       TBaseCell *Cell, int N_Matrices,
//                                       int N_Rhs,
//                                       double ***LocMatrix, double **LocRhs,
//                                       double factor)
// {
//   int i,j, N_Rows, N_Columns;
//   double **CurrentMatrix, *MatrixRow;
//   double Mult, *Coeff, *Param;
//   const double hK = Cell->Get_hK(TDatabase::ParamDB->CELL_MEASURE);
// 
//   //this->GetParameters(N_Points, NULL, Cell, ??Cell->GetCellIndex(), X, Y, 
//   //                    Parameters);
// 
//   for(i=0; i<N_Matrices; ++i)
//   {
//     CurrentMatrix = LocMatrix[i];
//     N_Rows = N_BaseFuncts[RowSpace[i]];
//     N_Columns = N_BaseFuncts[ColumnSpace[i]];
//     for(j=0;j<N_Rows;j++)
//     {
//       MatrixRow = CurrentMatrix[j];
//       memset(MatrixRow, 0, SizeOfDouble*N_Columns);
//     } // endfor j
//   } // endfor i
// 
//   for(i=0; i<N_Rhs; ++i)
//   {
//     N_Rows = N_BaseFuncts[RhsSpace[i]];
//     memset(LocRhs[i], 0, SizeOfDouble*N_Rows);
//   }
// 
//   // *****************************************************
//   // for 2Phase flow problems (Sashikumaar Ganesan)
//   AuxArray[0][0] = Cell->GetPhase_ID();
//   // *****************************************************
// 
//   if(Coeffs)
//     Coeffs(N_Points, X, Y, Parameters, AuxArray);
// 
//   if(Manipulate)
//     Manipulate(N_Points, AuxArray, Parameters, Cell);
// 
//   for(i=0; i<N_Terms; ++i)
//   {
//     AllOrigValues[i] = 
//       TFEDatabase2D::GetOrigElementValues(BaseFuncts[FESpaceNumber[i]], 
//                                           Derivatives[i]);
//   }
// 
//   for(i=0; i<N_Points; ++i)
//   {
//     Mult = weights[i] * AbsDetjk[i] * factor;
//     Coeff = AuxArray[i];
//     Coeff[19] = AbsDetjk[i];
//     
//     if(TDatabase::ParamDB->Axial3DAxis == 1)
//     {
//       // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)
//       Coeff[20] = Y[i];
//     }
//     else
//     {
//       // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)
//       Coeff[20] = X[i];
//     }
// 
//     Param = Parameters[i];
// 
//     for(j=0; j<N_Terms; j++)
//       OrigValues[j] = AllOrigValues[j][i];
// 
//     AssembleParam(Mult, Coeff, Param, hK, OrigValues, N_BaseFuncts, LocMatrix,
//                   LocRhs);
//   } // end loop over quadrature points 
// }
// 
// 
// void LocalAssembling2D::GetParameters(int n_points, TCollection *Coll,
//                                       TBaseCell *cell, int cellnum,
//                                       double *x, double *y, double **Parameters)
// {
//   int j, k, l, n;
//   double *param, *currparam, s;
//   
//   double *CurrValues, *CurrOrigValues;
//   int *CurrIndex;
// 
//   int *N_BaseFunct = new int[N_FEValues];
//   double **Values = new double* [N_FEValues];
//   double ***orig_values = new double** [N_FEValues];
//   int **Index = new int* [N_FEValues];
//   double Temp[2 + N_FEValues];
//   // collect information
//   for(j=0; j<this->N_FEValues; j++)
//   {
//     TFEFunction2D *fefunction = this->FEFunctions2D[this->FEValue_FctIndex[j]];
//     
//     Values[j] = fefunction->GetValues();
//     
//     TFESpace2D *fespace = fefunction->GetFESpace2D();
//     FE2D FE_Id = fespace->GetFE2D(cellnum, cell);
//     BaseFunct2D BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_Id)->GetBaseFunct2D_ID();
// 
//     N_BaseFunct[j]=TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id)->GetDimension();
//     
//     orig_values[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, 
//                                                          FEValue_MultiIndex[j]);
//     Index[j] = fespace->GetGlobalDOF(cellnum);
//   } // endfor j
// 
//   // loop over all quadrature points
//   if(N_ParamFct != 0)
//   {
//     for(int i=0; i<n_points; ++i)
//     {
//       param = Parameters[i];
// 
//       Temp[0] = x[i];
//       Temp[1] = y[i];
// 
//       // loop to calculate all FE values
//       for(k=2,j=0; j<N_FEValues; j++,k++)
//       {
//         s = 0;
//         n = N_BaseFunct[j];
//         CurrValues = Values[j];
//         CurrOrigValues = orig_values[j][i];
//         CurrIndex = Index[j];
//         for(l=0;l<n;l++)
//           s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
//         Temp[k] = s;
//       }  // endfor j
// 
//       // loop to calculate all parameters
//       for(j=0; j<N_ParamFct; j++)
//       {
//         currparam = param + this->BeginParameter[j];
//         this->ParameterFct[j](Temp, currparam);
//       } // endfor j
//     } // endfor i
//   }
//   
//   delete [] N_BaseFunct;
//   delete [] Values;
//   delete [] orig_values;
//   delete [] Index;
// }
// 
// 
// void LocalAssembling2D::set_parameters_for_nse(LocalAssembling2D_type type)
// {
//   switch(type)
//   {
//     case NSE2D_Galerkin:
//     {
//       switch(TDatabase::ParamDB->NSTYPE)
//       {
//         case 1: // NSE2D_Galerkin, NSTYPE=1,
//         {
//           if(TDatabase::ParamDB->LAPLACETYPE != 0)
//           {
//             ErrMsg("LAPLACETYPE must be set to 0 in case of NSTYPE 1");
//             exit(1);
//           }
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=1, NSE_NONLINEAR_FORM=0
//             {
//               this->N_Terms = 4;
//               this->Derivatives = { D10, D01, D00, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 3;
//               this->RowSpace = { 0, 1, 1 };
//               this->ColumnSpace = { 0, 0, 0 };
//               this->N_Rhs = 2;
//               this->RhsSpace = { 0, 0 };
//               this->AssembleParam = NSType1Galerkin; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=1, NSE_NONLINEAR_FORM=1
//             {
//               this->N_Terms = 4;
//               this->Derivatives = { D10, D01, D00, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 3;
//               this->RowSpace = { 0, 1, 1 };
//               this->ColumnSpace = { 0, 0, 0 };
//               this->N_Rhs = 2;
//               this->RhsSpace = { 0, 0 };
//               this->AssembleParam = NSType1GalerkinSkew; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2:
//             {
//               ErrMsg("Using the rotational form (NSE_NONLINEAR_FORM: 2) is not "
//                      << "possible with NSTYPE: 1. Choose NSTYPE: 3, 4 or 14");
//               exit(1);
//             }
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=1
//         case 2: // NSE2D_Galerkin, NSTYPE=2,
//         {
//           if(TDatabase::ParamDB->LAPLACETYPE != 0)
//           {
//             ErrMsg("LAPLACETYPE must be set to 0 in case of NSTYPE 2");
//             exit(1);
//           }
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=2, NSE_NONLINEAR_FORM=0
//             {
//               this->N_Terms = 4;
//               this->Derivatives = { D10, D01, D00, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 5;
//               this->RowSpace = { 0, 1, 1, 0, 0 };
//               this->ColumnSpace = { 0, 0, 0, 1, 1 };
//               this->N_Rhs = 2;
//               this->RhsSpace = { 0, 0 };
//               this->AssembleParam = NSType2Galerkin; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=2, NSE_NONLINEAR_FORM=1
//             {
//               this->N_Terms = 4;
//               this->Derivatives = { D10, D01, D00, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 5;
//               this->RowSpace = { 0, 1, 1, 0, 0 };
//               this->ColumnSpace = { 0, 0, 0, 1, 1 };
//               this->N_Rhs = 2;
//               this->RhsSpace = { 0, 0 };
//               this->AssembleParam = NSType2GalerkinSkew; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2:
//             {
//               ErrMsg("Using the rotational form (NSE_NONLINEAR_FORM: 2) is not "
//                      << "possible with NSTYPE: 2. Choose NSTYPE: 3, 4 or 14");
//               exit(1);
//             }
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=2
//         case 3: // NSE2D_Galerkin, NSTYPE=3,
//         {
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3Galerkin; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3GalerkinDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3GalerkinSkew; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3GalerkinSkewDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2,
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3GalerkinRot; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 6;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType3GalerkinRotDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=2
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=3
//         case 4: // NSE2D_Galerkin, NSTYPE=4,
//         case 14: // NSE2D_Galerkin, NSTYPE=14,
//         {
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4Galerkin; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4GalerkinDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4GalerkinSkew; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4GalerkinSkewDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2,
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4GalerkinRot; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 4;
//                   this->Derivatives = { D10, D01, D00, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0, 1 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 8;
//                   this->RowSpace = { 0, 0, 0, 0, 1, 1, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0, 0, 0, 1, 1 };
//                   this->N_Rhs = 2;
//                   this->RhsSpace = { 0, 0 };
//                   this->AssembleParam = NSType4GalerkinRotDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=2
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=(1)4
//         default:
//           ErrMsg("unknown NSTYPE " << TDatabase::ParamDB->NSTYPE);
//           exit(1);
//       } // end switch NSTYPE
//       break;
//     } // end case LocalAssembling2D_type=NSE2D_Galerkin
//     case NSE2D_Galerkin_Nonlinear:
//     {
//       switch(TDatabase::ParamDB->NSTYPE)
//       {
//         case 1: // NSE2D_Galerkin, NSTYPE=1,
//         {
//           if(TDatabase::ParamDB->LAPLACETYPE != 0)
//           {
//             ErrMsg("LAPLACETYPE must be set to 0 in case of NSTYPE 1");
//             exit(1);
//           }
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=1, NSE_NONLINEAR_FORM=0
//             {
//               this->N_Terms = 3;
//               this->Derivatives = { D10, D01, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 1;
//               this->RowSpace = { 0 };
//               this->ColumnSpace = { 0 };
//               this->N_Rhs = 0;
//               this->RhsSpace = {};
//               this->AssembleParam = NSType1_2NLGalerkin; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=1, NSE_NONLINEAR_FORM=1
//             {
//               this->N_Terms = 3;
//               this->Derivatives = { D10, D01, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 1;
//               this->RowSpace = { 0 };
//               this->ColumnSpace = { 0 };
//               this->N_Rhs = 0;
//               this->RhsSpace = {};
//               this->AssembleParam = NSType1_2NLGalerkinSkew; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2:
//             {
//               ErrMsg("Using the rotational form (NSE_NONLINEAR_FORM: 2) is not "
//                      << "possible with NSTYPE: 1. Choose NSTYPE: 3, 4 or 14");
//               exit(1);
//             }
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=1
//         case 2: // NSE2D_Galerkin, NSTYPE=2,
//         {
//           if(TDatabase::ParamDB->LAPLACETYPE != 0)
//           {
//             ErrMsg("LAPLACETYPE must be set to 0 in case of NSTYPE 2");
//             exit(1);
//           }
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=2, NSE_NONLINEAR_FORM=0
//             {
//               this->N_Terms = 3;
//               this->Derivatives = { D10, D01, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 1;
//               this->RowSpace = { 0 };
//               this->ColumnSpace = { 0 };
//               this->N_Rhs = 0;
//               this->RhsSpace = {};
//               this->AssembleParam = NSType1_2NLGalerkin; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=2, NSE_NONLINEAR_FORM=1
//             {
//               this->N_Terms = 3;
//               this->Derivatives = { D10, D01, D00 };
//               this->Needs2ndDerivatives = new bool[1];
//               this->Needs2ndDerivatives[0] = false;
//               this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//               this->N_Matrices = 1;
//               this->RowSpace = { 0 };
//               this->ColumnSpace = { 0 };
//               this->N_Rhs = 0;
//               this->RhsSpace = {};
//               this->AssembleParam = NSType1_2NLGalerkinSkew; 
//               this->Manipulate = NULL;
//               
//               this->N_Parameters = 2;
//               this->N_ParamFct = 1;
//               this->ParameterFct =  { NSParamsVelo };
//               this->N_FEValues = 2;
//               this->FEValue_FctIndex = { 0, 1 };
//               this->FEValue_MultiIndex = { D00, D00 };
//               this->BeginParameter = { 0 };
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2:
//             {
//               ErrMsg("Using the rotational form (NSE_NONLINEAR_FORM: 2) is not "
//                      << "possible with NSTYPE: 2. Choose NSTYPE: 3, 4 or 14");
//               exit(1);
//             }
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=2
//         case 3: // NSE2D_Galerkin, NSTYPE=3,
//         {
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkin; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=0
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinSkew; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=1
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinSkewDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2,
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 4;
//                   this->RowSpace = { 0, 0, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinRot; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=3, NSE_NONLINEAR_FORM=2
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 4;
//                   this->RowSpace = { 0, 0, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinRotDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=2
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=3
//         case 4: // NSE2D_Galerkin, NSTYPE=4,
//         case 14: // NSE2D_Galerkin, NSTYPE=14,
//         {
//           switch(TDatabase::ParamDB->NSE_NONLINEAR_FORM)
//           {
//             case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkin; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=0
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=0
//             case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1, 
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinSkew; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=1
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 2;
//                   this->RowSpace = { 0, 0 };
//                   this->ColumnSpace = { 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinSkewDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=1
//             case 2: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2
//             {
//               switch(TDatabase::ParamDB->LAPLACETYPE)
//               {
//                 case 0: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2,
//                         // LAPLACETYPE=0
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 4;
//                   this->RowSpace = { 0, 0, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinRot; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=0
//                 case 1: // NSE2D_Galerkin, NSTYPE=(1)4, NSE_NONLINEAR_FORM=2
//                         // LAPLACETYPE=1
//                 {
//                   this->N_Terms = 3;
//                   this->Derivatives = { D10, D01, D00 };
//                   this->Needs2ndDerivatives = new bool[1];
//                   this->Needs2ndDerivatives[0] = false;
//                   this->FESpaceNumber = { 0, 0, 0 }; // 0: velocity, 1: pressure
//                   this->N_Matrices = 4;
//                   this->RowSpace = { 0, 0, 0, 0 };
//                   this->ColumnSpace = { 0, 0, 0, 0 };
//                   this->N_Rhs = 0;
//                   this->RhsSpace = {};
//                   this->AssembleParam = NSType3_4NLGalerkinRotDD; 
//                   this->Manipulate = NULL;
//                   
//                   this->N_Parameters = 2;
//                   this->N_ParamFct = 1;
//                   this->ParameterFct =  { NSParamsVelo };
//                   this->N_FEValues = 2;
//                   this->FEValue_FctIndex = { 0, 1 };
//                   this->FEValue_MultiIndex = { D00, D00 };
//                   this->BeginParameter = { 0 };
//                   break;
//                 } // end case LAPLACETYPE=1
//                 default:
//                   ErrMsg("unknown LAPLACETYPE " 
//                          << TDatabase::ParamDB->LAPLACETYPE);
//                   exit(1);
//               } // end switch LAPLACETYPE
//               break;
//             } // end case NSE_NONLINEAR_FORM=2
//             default:
//               ErrMsg("unknown NSE_NONLINEAR_FORM "
//                      << TDatabase::ParamDB->NSE_NONLINEAR_FORM);
//               exit(1);
//           } // end switch NSE_NONLINEAR_FORM
//           break;
//         } // end case NSTYPE=(1)4
//         default:
//           ErrMsg("unknown NSTYPE " << TDatabase::ParamDB->NSTYPE);
//           exit(1);
//       } // end switch NSTYPE
//       break;
//     } // end case LocalAssembling2D_type=NSE2D_Galerkin_Nonlinear
//     default:
//       ErrMsg("unknown LocalAssembling2D_type " << type << "  " << this->name);
//       exit(1);
//   } // end switch LocalAssembling2D_type
// }
