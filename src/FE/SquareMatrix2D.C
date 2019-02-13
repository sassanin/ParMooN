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
// @(#)SquareMatrix2D.C        1.2 11/20/98
//
// Class:       TSquareMatrix2D
//
// Purpose:     store a square matrix (ansatz = test space) in 2d
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#include <SquareMatrix2D.h>
#include <string.h>
#include <stdlib.h>

TSquareMatrix2D::TSquareMatrix2D(TSquareStructure2D *squarestructure)
  : TSquareMatrix(squarestructure), structure(squarestructure)
{
}


// set all class variables
void TSquareMatrix2D::SetStructure(TSquareStructure2D *squarestructure)
{
  structure = squarestructure;
  this->TSquareMatrix::SetStructure((TSquareStructure*) squarestructure);
}

// not clear! what is the initial value of structure while calling contructor/
// commented by Sashi 
// TSquareMatrix2D::TSquareMatrix2D(int n) 
//  : structure(new TSquareStructure2D(n)), TSquareMatrix(structure)
// {
// }

TSquareMatrix2D::~TSquareMatrix2D()
{
}


TSquareMatrix2D& TSquareMatrix2D::operator=(const TSquareMatrix2D& rhs)
{
  // compare structures (first pointers, then structures themselves)
  if(structure != rhs.structure && *structure != *rhs.structure)
  {
    OutPut("WARNING: TSquareMatrix2D& operator= Matrices have different "
           << "fe space or structure. Are you sure this is what you want?" << endl);
    structure = rhs.structure;
  }
  // no further tests, make sure you know these two matrices have the same 
  // structure. We cannot compare the two TSquareStructure2D, because during 
  // intitialization every matrix might get its own TSquareStructure2D.
  
  // copy matrix entries.
  for(int i=0; i<structure->GetN_Entries(); i++)
    Entries[i] = rhs.Entries[i];
  
  return *this;
}


TSquareMatrix2D& TSquareMatrix2D::operator*=(double alpha)
{
  int nDOFNonActive = structure->GetN_Rows() - structure->GetActiveBound();
  for (int i=0; i<structure->GetN_Entries() - nDOFNonActive; i++){
      Entries[i] *= alpha;
  }
  return *this;
}


TSquareMatrix2D& TSquareMatrix2D::operator+=(TSquareMatrix2D & rhsMat)
{

  if(*structure != *rhsMat.structure)
  {
    OutPut("ERROR: TSquareMatrix2D::operator+=() The two arguments "
           << "have different structures. Exiting" << endl);
    exit(1);
  }
  //OutPut(" TSquareMatrix2D::operator+= \n");
  double *rhsEntries = rhsMat.Entries;
  
  int nDOFNonActive = structure->GetN_Rows() - structure->GetActiveBound();
  for (int i=0; i<structure->GetN_Entries() - nDOFNonActive; i++){
      Entries[i] += rhsEntries[i];
  }
  return *this;
}

// overloaded operators 
// add to matrices A and B
// note: only active DOF are added
// note: only works for matrices with the same sparsity pattern
TSquareMatrix2D& operator+(const TSquareMatrix2D & A, const TSquareMatrix2D & B)
{
  double *AEntries, *BEntries, *CEntries;
  if (A.GetMatrixStructure() == B.GetMatrixStructure()) 
  {
    // create bew TSquareMatrix2D on heap (otherwise return did not work)
    TSquareMatrix2D *C = new TSquareMatrix2D(A.GetMatrixStructure());
    AEntries = A.GetEntries();
    BEntries = B.GetEntries();
    CEntries = C->GetEntries();
    for (int i=0; i<A.GetActiveBound(); i++) 
    {
      CEntries[i] = AEntries[i] + BEntries[i];
    }
    return *C;
  } else 
  {
    cout << " SquareMatrix2D: ERROR: can't add Matrices "
       << " with different sparse structures " << endl;
    exit(1);
  }
}

// C= A*alpha 
// note: only active DOF are multiplied, others are just copied
TSquareMatrix2D& operator*(const TSquareMatrix2D & A,const double alpha)
{
  double *AEntries, *CEntries;
  TSquareMatrix2D *C = new TSquareMatrix2D(A.GetMatrixStructure());
  AEntries = A.GetEntries();
  CEntries = C->GetEntries();
  // multiply each active entry by alpha and write it into matrix C
  for (int i=0; i<A.GetActiveBound(); i++) 
  {
    CEntries[i] = alpha*AEntries[i];
  }
  // non active entries are just copied
  for (int i=A.GetActiveBound(); i<A.GetN_Entries(); i++) 
  {
    CEntries[i] = AEntries[i];
  }
  return *C;
}
TSquareMatrix2D& operator*(const double alpha,const TSquareMatrix2D & A)
{// just to allow alpha*A as well (in addition to A*alpha)
  return A*alpha;
}



double* operator*(const TSquareMatrix2D & A,const double* x)
{
  double *AEntries;
  int *ARowPtr,*AColIndex;
  AEntries = A.GetEntries();
  ARowPtr = A.GetRowPtr();
  AColIndex = A.GetKCol();

  int nDOFActive = A.GetActiveBound();
  int nDOF = A.GetFESpace()->GetN_DegreesOfFreedom();
  double *y=new double[nDOF];
  double value;
  int index;

  // multiply each active entry by alpha and write it into y
  for(int i=0;i<nDOFActive;i++) {
    value = 0;
    for (int j=ARowPtr[i]; j<ARowPtr[i+1]; j++) {
      index = AColIndex[j];
      value += AEntries[j] * x[index];
    }
    y[i] = value;
  }
  
  for (int i=nDOFActive; i<nDOF; i++) {
    y[i]=x[i];
  }
  
  return y;
}
