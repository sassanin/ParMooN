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
   
#include <AuxParam2D3D.h>
#include <FEDatabase3D.h>
#include <TetraIsoparametric.h>
#include <TetraAffin.h>

// CTOR
TAuxParam2D3D::TAuxParam2D3D (int n_fefunctions, TFEFunction3D **fefunctions,
		   int n_fevectfuncts, TFEVectFunct3D **fevectfuncts,
		   int n_parameters, int *fct_index, MultiIndex3D *derivatives,
		   int *isglobal)
{
  int ivar;
  
  mN_FEFunctions = n_fefunctions;
  mFEFunctions = NULL;
  
  if ( n_fefunctions )
  {
    mFEFunctions = new TFEFunction3D* [n_fefunctions];
    
    for (int i=0;i<n_fefunctions;++i)
      mFEFunctions[i] = fefunctions[i];
  }
  
  mN_FEVectFuncts = n_fevectfuncts;
  mFEVectFuncts = NULL;
  
  if ( n_fevectfuncts )
  {
    mFEVectFuncts = new TFEVectFunct3D* [n_fevectfuncts];
    
    for (int i=0;i<n_fevectfuncts;++i)
      mFEVectFuncts[i] = fevectfuncts[i];
  }
  
  mN_TotalParameters = 0;
  for (int i=0;i<n_parameters;++i)
  {
    ivar = fct_index[i]; 
    
    if ( ivar < n_fefunctions ) ++mN_TotalParameters;
    else mN_TotalParameters += 3;
  }
  
  mFct_Index = fct_index;
  mDerivatives = derivatives;
  mN_Parameters = n_parameters;
  mIsGlobal = isglobal;
}

// DTOR
TAuxParam2D3D::~TAuxParam2D3D ()
{
  if (mFEFunctions) delete [] mFEFunctions;
  if (mFEVectFuncts) delete [] mFEVectFuncts;
}

/// Methods
void TAuxParam2D3D::GetParameters(TBaseCell *cell, int cellnr, int jointnr, int globnr,
		      RefTrans3D RefTrans,
		      int n_points, double *xi, double *eta,
		      double **parameters)
{
  int Fct, index;
  double t11, t12, t13, t21, t22, t23, norm;
  TRefTrans3D *F_K;  
  
  // find reftrans
  F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
  
  // find normals
  mProject = true;
  for (int i=0;i<n_points;++i)
  {
    switch (RefTrans)
    {
      case TetraIsoparametric:
	((TTetraIsoparametric*) F_K)->GetTangentVectors(jointnr, xi[i], eta[i],
							  t11, t12, t13,
							  t21, t22, t23);
	break;
	 
      case TetraAffin:
	((TTetraAffin*) F_K)->GetTangentVectors(jointnr, xi[i], eta[i],
						  t11, t12, t13,
						  t21, t22, t23);
	break;
	  
	default:
	  mProject = false;
    }
    
    n1[i] = t12*t23 - t13*t22;
    n2[i] = t13*t21 - t11*t23;
    n3[i] = t11*t22 - t12*t21;
      
    norm = sqrt(n1[i]*n1[i]+n2[i]*n2[i]+n3[i]*n3[i]);
      
    n1[i] /= norm;
    n2[i] /= norm;
    n3[i] /= norm;
  }
  
  mCell = cell;
  mCellNr = cellnr;
  mGlobCellNr = globnr;
  mN_Points = n_points;
  index = 0;
  for (int i=0;i<mN_Parameters;++i)
  {
    Fct = mFct_Index[i];
    
    if ( Fct < mN_FEFunctions )
    {
      GetFunctionsParams(Fct, mDerivatives[i], index, parameters);
      ++index;
    }
    else
    {
      GetVectFunctParams(Fct-mN_FEFunctions, index, parameters);
      index += 3;
    }
  }
}

void TAuxParam2D3D::GetFunctionsParams(int Fct, MultiIndex3D Der, int index, double **parameters)
{
  int N_BaseFuncts, *GlobalNumbers, *BeginIndex, *DOF, CellNr;
  double *values, **origvalues, *orig, *param, val;
  TFESpace3D *fespace;
  FE3D FeID;
  BaseFunct3D BaseFunct;
  
  values = mFEFunctions[Fct]->GetValues();
  fespace = mFEFunctions[Fct]->GetFESpace3D();
  
  CellNr = mCellNr;
  if ( mIsGlobal[Fct] ) CellNr = mGlobCellNr;
  
  FeID = fespace->GetFE3D(CellNr, mCell);
  BaseFunct = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(FeID);
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D(FeID);
  
  origvalues = TFEDatabase3D::GetOrigElementValues(BaseFunct, Der);
  
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex    = fespace->GetBeginIndex();
  DOF = GlobalNumbers + BeginIndex[CellNr];
  
  for (int i=0;i<mN_Points;++i)
  {
    orig = origvalues[i];
    param = parameters[i];
    
    val = 0;
    for (int j=0;j<N_BaseFuncts;++j)
    {
      val += orig[j]*values[DOF[j]];
    }
    
    param[index] = val;
  }
}

void TAuxParam2D3D::GetVectFunctParams(int Fct, int index, double **parameters)
{
  int N_BaseFuncts, *GlobalNumbers, *BeginIndex, *DOF, N_U, CellNr;
  double *values, **origvalues, *orig, *param, valx, valy, valz, d;
  TFESpace3D *fespace;
  FE3D FeID;
  BaseFunct3D BaseFunct;
  
  values = mFEVectFuncts[Fct]->GetValues();
  fespace = mFEVectFuncts[Fct]->GetFESpace3D();
  N_U = mFEVectFuncts[Fct]->GetLength();
  
  CellNr = mCellNr;
  if ( mIsGlobal[Fct+mN_FEFunctions] ) CellNr = mGlobCellNr;
  
  FeID = fespace->GetFE3D(CellNr, mCell);
  BaseFunct = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D(FeID);
  N_BaseFuncts = TFEDatabase3D::GetN_BaseFunctFromFE3D(FeID);
    
  origvalues = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
  
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex    = fespace->GetBeginIndex();
  DOF = GlobalNumbers + BeginIndex[CellNr];
  
  for (int i=0;i<mN_Points;++i)
  {
    orig = origvalues[i];
    param = parameters[i];
    
    valx = valy = valz = 0;
    for (int j=0;j<N_BaseFuncts;++j)
    {
      valx += orig[j]*values[DOF[j]      ];
      valy += orig[j]*values[DOF[j]+  N_U];
      valz += orig[j]*values[DOF[j]+2*N_U];
    }
    
    if (mProject)
    {
      d = n1[i]*valx + n2[i]*valy + n3[i]*valz;
	
      valx -= d*n1[i];
      valy -= d*n2[i];
      valz -= d*n3[i];
    }    
    
    param[index+0] = valx;
    param[index+1] = valy;
    param[index+2] = valz;
  }
}
