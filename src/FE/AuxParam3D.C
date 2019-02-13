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
// %W% %G%
// 
// Class:       TAuxParam3D
// Purpose:     store parameter functions and FE functions
//
// Author:      Gunar Matthies (06.08.98)
//
// History:     start of implementation 06.08.98 (Gunar Matthies)
//
// =======================================================================

#include <AuxParam3D.h>
#include <FEDatabase3D.h>
#include <stdlib.h>

/** constructor */
TAuxParam3D::TAuxParam3D(
        int n_fespace3d, int n_fefunction3d, int n_paramfct,
        int n_fevalues,
        TFESpace3D **fespaces3d, TFEFunction3D **fefunctions3d,
        ParamFct **parameterfct,
        int *fevalue_fctindex, MultiIndex3D *fevalue_multiindex,
        int n_parameters, int *beginparameter)
{
  N_FESpace3D = n_fespace3d;
  N_FEFunction3D = n_fefunction3d;
  N_ParamFct = n_paramfct;
  N_FEValues = n_fevalues;

  FESpaces3D = fespaces3d;
  FEFunctions3D = fefunctions3d;
  ParameterFct = parameterfct;

  FEValue_FctIndex = fevalue_fctindex;
  FEValue_MultiIndex = fevalue_multiindex;

  N_Parameters = n_parameters;
  BeginParameter = beginparameter;

  Temp = new double[3 + N_FEValues];

  Values = new double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}


/** set aux parameters giving a keyword*/
TAuxParam3D::TAuxParam3D(std::string name, TFEFunction3D **fefunctions3d)
{
  if(name=="Velocity") // for Navier-Stokes
  {
    N_FESpace3D = 0;
    FESpaces3D = NULL;
    
    N_FEFunction3D = 3;
    FEFunctions3D = fefunctions3d;
    
    N_ParamFct = 1;
    N_FEValues = 3;
    
    //see TNSE3D_FixPo.C
    ParameterFct=new ParamFct*[1];
    ParameterFct[0] = &Velocity_Fct; // see below in this file
    FEValue_FctIndex =new int[3];
    FEValue_FctIndex[0] = 0;
    FEValue_FctIndex[1] = 1;
    FEValue_FctIndex[2] = 2;
    
    FEValue_MultiIndex=new MultiIndex3D[3];
    FEValue_MultiIndex[0]= D000;
    FEValue_MultiIndex[1]= D000;
    FEValue_MultiIndex[2]= D000;    
    
    N_Parameters = 3;
    BeginParameter =new int[1];
    BeginParameter[0]=0;
    
  }
  else
  {
    ErrMsg(" AuxParam3D:: Constructor: ERROR, name " << name 
           << " for initialization not imlpemented ");
    exit(1);
  }
  
  Temp = new double[3 + N_FEValues];

  Values = new double* [N_FEValues];
  OrigValues = new double** [N_FEValues];
  Index = new int* [N_FEValues];
  N_BaseFunct = new int[N_FEValues];
}



/** return all parameters at all quadrature points */
void TAuxParam3D::GetParameters(int N_Points, TBaseCell *cell, int cellnum,
                              double *Xi, double *Eta, double *Zeta,
                              double *X, double *Y, double *Z,
                              double **Parameters)
{
  int i, j, k, l, n;
  double xi, *param, *currparam, s;
  double *locvalues;
  TFESpace3D *fespace;
  TFEFunction3D *fefunction;
  FE3D FE_Id;
  BaseFunct3D BaseFunct_Id;
  int *GlobalNumbers, *BeginIndex;

  double *CurrValues, *CurrOrigValues;
  int *CurrIndex;

   // collect information
  for(j=0;j<N_FEValues;j++)
  {
    fefunction = FEFunctions3D[FEValue_FctIndex[j]];
    Values[j] = fefunction->GetValues();

    fespace = fefunction->GetFESpace3D();
    FE_Id = fespace->GetFE3D(cellnum, cell);
    BaseFunct_Id = TFEDatabase3D::GetFE3D(FE_Id)->GetBaseFunct3D_ID();

    N_BaseFunct[j]=TFEDatabase3D::GetBaseFunct3D(BaseFunct_Id)->GetDimension();
    OrigValues[j] = TFEDatabase3D::GetOrigElementValues
	(BaseFunct_Id, FEValue_MultiIndex[j]);

    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    Index[j] = GlobalNumbers + BeginIndex[cellnum];
  } // endfor j


  // loop over all quadrature points
  for(i=0;i<N_Points;i++)
  {
    // parameters in this quadrature point  
    param = Parameters[i];
    // first three parameters are the coordinates
    Temp[0] = X[i];
    Temp[1] = Y[i];
    Temp[2] = Z[i];

    // loop to calculate all FE values
    for(k=3,j=0;j<N_FEValues;j++,k++)
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
      
      /// change for twophase
      currparam[0] = cell->GetRegionID();
      
      ParameterFct[j](Temp, currparam);
    } // endfor j
  } // endfor i
}

/** destructor */
TAuxParam3D::~TAuxParam3D()
{
  delete [] Temp;
  delete [] Values;
  delete [] OrigValues;
  delete [] Index;
  delete [] N_BaseFunct;
}


void Velocity_Fct(double *inputList, double *outputValues)
{
  outputValues[0] = inputList[3];                // u1old
  outputValues[1] = inputList[4];                // u2old
  outputValues[2] = inputList[5];                // u3old
}
