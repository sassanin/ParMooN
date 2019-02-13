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
// @(#)AuxParam2D.C        1.2 09/17/99
// 
// Class:       TAuxParam2D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#include <AuxParam2D.h>
#include <FEDatabase2D.h>
#include <MooNMD_Io.h>
// #include <NSE2D_ParamRout.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

/** constructor */
TAuxParam2D::TAuxParam2D(
        int n_fespace2d, int n_fefunction2d, int n_paramfct,
        int n_fevalues,
        TFESpace2D **fespaces2d, TFEFunction2D **fefunctions2d,
        ParamFct **parameterfct,
        int *fevalue_fctindex, MultiIndex2D *fevalue_multiindex,
        int n_parameters, int *beginparameter)
{
  N_FESpace2D = n_fespace2d;
  N_FEFunction2D = n_fefunction2d;
  N_ParamFct = n_paramfct;
  N_FEValues = n_fevalues;

  FESpaces2D = fespaces2d;
  FEFunctions2D = fefunctions2d;
  ParameterFct = parameterfct;

  FEValue_FctIndex = fevalue_fctindex;
  FEValue_MultiIndex = fevalue_multiindex;

  N_Parameters = n_parameters;
  BeginParameter = beginparameter;

  Temp = new double[2 + N_FEValues];

  Values = new double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}

//   !!!! Error in Mac Compiler - Sashikumaar !!!!!!!
// // empty aux class
// TAuxParam2D::TAuxParam2D() : TAuxParam2D(0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, 0, NULL)
// {
// }

// // set aux parameters giving a keyword
// TAuxParam2D::TAuxParam2D(std::string name, TFEFunction2D **fefunctions2d)
// {
//   if (name=="Velocity") // for Navier-Stokes
//   {
//     N_FESpace2D = 0;
//     FESpaces2D = NULL;
//    
//     N_FEFunction2D = 2;
//     FEFunctions2D = fefunctions2d;
// 
//     N_ParamFct = 1;
//     N_FEValues = 2;
//     
//     // for all arrays defined below: see NSE2D_ParamRout.h and NSE2D_FixPo.C
//     ParameterFct = NSFctVelo;
//     FEValue_FctIndex = NSFEFctIndexVelo; // {0,1}
//     FEValue_MultiIndex = NSFEMultiIndexVelo; // {D00,D00}
//     
//     N_Parameters = 2;
//     BeginParameter = NSBeginParamVelo; // {0}
//     
//   } 
//   else
//   {
//     cout << " AuxParam2D:: Constructor: ERROR, name " << name 
//          << " for initialization not imlpemented " << endl;
//     exit(1);
//   }
//   
//   Temp = new double[2 + N_FEValues];
//   
//   Values = new double* [N_FEValues];
//   OrigValues = new double** [N_FEValues];
//   Index = new int* [N_FEValues];
//   N_BaseFunct = new int[N_FEValues];
// 
// 
// }



/** return all parameters at all quadrature points */
void TAuxParam2D::GetParameters(int N_Points, TCollection *Coll,
                              TBaseCell *cell, int cellnum,
                              double *Xi, double *Eta,
                              double *X, double *Y,
                              double **Parameters)
{
  int i, j, k, l, n;
  double xi, *param, *currparam, s;
  double *locvalues;
  TFESpace2D *fespace;
  TFEFunction2D *fefunction;
  FE2D FE_Id;
  BaseFunct2D BaseFunct_Id;
  int *GlobalNumbers, *BeginIndex;

  double *CurrValues, *CurrOrigValues;
  int *CurrIndex;

  // collect information
  for(j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions2D[FEValue_FctIndex[j]];
    
    Values[j] = fefunction->GetValues();
    //  if (N_FEValues==8)
    //  OutPut("aac " << (int) fefunction << " " <<  Values[j][0]<< endl);

    fespace = fefunction->GetFESpace2D();
    FE_Id = fespace->GetFE2D(cellnum, cell);
    BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_Id)->GetBaseFunct2D_ID();

    N_BaseFunct[j]=TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id)->GetDimension();
    
    OrigValues[j] = TFEDatabase2D::GetOrigElementValues(BaseFunct_Id, FEValue_MultiIndex[j]);
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    Index[j] = GlobalNumbers + BeginIndex[cellnum];
  } // endfor j

  // loop over all quadrature points
  for(i=0;i<N_Points;i++)
  {
    param = Parameters[i];

    Temp[0] = X[i];
    Temp[1] = Y[i];

    // loop to calculate all FE values
    for(k=2,j=0;j<N_FEValues;j++,k++)
    {
      s = 0;
      n = N_BaseFunct[j];
      CurrValues = Values[j];
      CurrOrigValues = OrigValues[j][i];
      CurrIndex = Index[j];
      for(l=0;l<n;l++)
        s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(j=0;j<N_ParamFct;j++)
    {
      currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
}

/** return all parameters at all quadrature points on the boundary*/
void TAuxParam2D::GetParameters(int N_Points, TCollection *Coll, TBaseCell *cell, int cellnum,
				double *t, int joint, double **Parameters)
{
  int i,j,k,l,n, N_Cells;
  double xv, yv, xi, eta, eps = 1e-20;
  double s;
  double *param, *currparam, *CurrValues, *CurrOrigValues;
  int *CurrIndex;
  TFESpace2D *fespace;
  TFEFunction2D *fefunction;
  BaseFunct2D BaseFunct_Id;
  FE2D FE_ID;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf, **AllBaseFuncts;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D], uref[MaxN_BaseFunctions2D];
  double uxiref[MaxN_BaseFunctions2D], uetaref[MaxN_BaseFunctions2D];
  
  int *Numbers;
  double u, ux, uy;
  double val;
  int *GlobalNumbers, *BeginIndex;

  double X, Y, absdetjk;

  AllBaseFuncts = new TBaseFunct2D*[N_FEValues];

  for(j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions2D[FEValue_FctIndex[j]];
    
    Values[j] = fefunction->GetValues();

    fespace = fefunction->GetFESpace2D();
    FE_ID = fespace->GetFE2D(cellnum, cell);
    BaseFunct_Id = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D_ID();

    AllBaseFuncts[j] = TFEDatabase2D::GetBaseFunct2D(BaseFunct_Id);
    N_BaseFunct[j] = AllBaseFuncts[j]->GetDimension();
    
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    Index[j] = GlobalNumbers + BeginIndex[cellnum];
  } // endfor j

  FE_Obj = TFEDatabase2D::GetFE2D(FE_ID);
  RefTrans = FE_Obj->GetRefTransID();
  // set cell for reference transformation
  TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
  
  for(i=0;i<N_Points;i++)
  {
    switch(joint)
    {
    case 0: 
      xi=t[i];
      eta=-1;
      break; 
    case 1: 
      xi=1;
      eta=t[i];
      break;
    case 2: 
      xi=-t[i];
      eta=1;
      break;
    case 3: 
      xi=-1;
      eta=-t[i];
      break;
    }//switch
//     cout << "xi eta " << xi << " " << eta << endl;

    param = Parameters[i];

    TFEDatabase2D::GetOrigFromRef(RefTrans, 1, &xi, &eta, &X, &Y, &absdetjk);
    // Temp[0] = X[i];
    // Temp[1] = Y[i];

    // loop to calculate all FE values
    for(k=2,j=0;j<N_FEValues;j++,k++)
    {
      s = 0;
      n = N_BaseFunct[j];
      CurrValues = Values[j];

      // get values and derivatives of basis functions on the
      // reference mesh cell
      bf = AllBaseFuncts[j];
      bf->GetDerivatives(D00, xi, eta, uref);
      bf->GetDerivatives(D10, xi, eta, uxiref);
      bf->GetDerivatives(D01, xi, eta, uetaref);
  
      // compute values on the original mesh cell 
      TFEDatabase2D::GetOrigValues(RefTrans, xi, eta, bf, Coll, (TGridCell *)cell,
				   uref, uxiref, uetaref, 
				   uorig, uxorig, uyorig);
      switch(FEValue_MultiIndex[j])
      {
        case D00:
          CurrOrigValues = uorig;
	break;
        case D10:
          CurrOrigValues = uxorig;
	break;
        case D01:
          CurrOrigValues = uyorig;
	break;
        default:
         cerr << "Second derivatives not added, see AuxParam2D "  << endl;
          exit (-1);
         break;
      } // endswitch
      // cout << "CurrOrigValues" << *uorig << endl;
      // cout << "CurrOrigValuesx" << *uxorig << endl;
      CurrIndex = Index[j];
      for(l=0;l<n;l++)
      { 
	s += CurrValues[CurrIndex[l]]*CurrOrigValues[l];
	// cout << "Parameter   " << CurrOrigValues[l] << endl;
      }
      Temp[k] = s;
    }  // endfor j

    // loop to calculate all parameters
    for(j=0;j<N_ParamFct;j++)
    {
      currparam = param + BeginParameter[j];
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i 
  delete AllBaseFuncts;
}

/** destructor */
TAuxParam2D::~TAuxParam2D()
{
  delete [] Temp;
  delete [] Values;
  delete [] OrigValues;
  delete [] Index;
  delete [] N_BaseFunct;
}
