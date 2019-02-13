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
// @(#)SquareMatrix.C        1.2 11/20/98
// 
// Class:       TSquareMatrix
//
// Purpose:     store a square matrix (ansatz = test space) in 
//
// Author:      Gunar Matthies
//
// History:     10.08.1998 start implementation
//
// =======================================================================

#include <SquareMatrix.h>
#include <Constants.h>
#include <string.h>
#include <MooNMD_Io.h>
#include <stdlib.h>


TSquareMatrix::TSquareMatrix(TSquareStructure *squarestructure)
 : TMatrix(squarestructure), structure(squarestructure)
{
}

// not clear! what is the initial value of structure while calling contructor/
// commented by Sashi 
// TSquareMatrix::TSquareMatrix(int n)
//  : structure(new TSquareStructure(n)), TMatrix(structure)
// {
// }

TSquareMatrix::~TSquareMatrix()
{
}

void TSquareMatrix::SetStructure(TSquareStructure *structure)
{
  this->structure = structure;
  this->TMatrix::SetStructure((TStructure*)structure);
}

void TSquareMatrix::ResetActive()
{
  memset(Entries, 0, (structure->GetN_Entries() - structure->GetN_Rows()
      + structure->GetActiveBound())*SizeOfDouble);
}

void TSquareMatrix::resetNonActive()
{
  int n_nonactive_rows = structure->GetN_Rows()-structure->GetActiveBound();
  int index_nonactive = structure->GetN_Entries() - n_nonactive_rows;
  memset(Entries + index_nonactive, 0.0, n_nonactive_rows*SizeOfDouble);
}

// find a renumbering of the DOFs from the matrix entries
void TSquareMatrix::ReNumbering(int* &Numbers) const
{
  int i,j,k,l,N_;
  int *Inflow, N_Active;
  int begin, end, beginJ, endJ;
  int Index, State, CurrentNumber;
  double aij, aji;

  N_Active = this->structure->GetActiveBound();
  Inflow = new int[N_Active];
  memset(Inflow, 0, N_Active*SizeOfInt);
  Numbers = new int[N_Active];
  memset(Numbers, -1, N_Active*SizeOfInt);
  int *RowPtr = this->structure->GetRowPtr();
  int *KCol = this->structure->GetKCol();

  for(i=0;i<N_Active;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin+1;j<end;j++)
    {
      k = KCol[j];
      if( (k>i) && (k<N_Active))
      {
        aij = Entries[j];
        beginJ = RowPtr[k];
        endJ = RowPtr[k+1];
        for(l=beginJ+1;l<endJ;l++)
        {
          if( KCol[l] == i )
          {
            aji = Entries[l];
            if(aij>aji) Inflow[k]++;
            if(aji>aij) Inflow[i]++;
          }
        }
      } // endif
    } // endfor j
  } // endfor i

  Index = 0; 
  for(i=0;i<N_Active;i++)
  {
    cout << i << "    " << Inflow[i] << endl;
    if(Inflow[i] == 0)
    {
      Numbers[Index] = i;
      Index++;
    }
  } // endfor i

  State = 0;
  while(State<N_Active)
  {
    CurrentNumber = Numbers[State];
    begin = RowPtr[CurrentNumber];
    end = RowPtr[CurrentNumber+1];
    for(j=begin+1;j<end;j++)
    {
      k = KCol[j];
      if(k<N_Active)
      {
        aij = Entries[j];
        beginJ = RowPtr[k];
        endJ = RowPtr[k+1];
        for(l=beginJ+1;l<endJ;l++)
        {
          if( KCol[l] == CurrentNumber )
          {
            aji = Entries[l];
            if(aij>aji)
            {
              if(Inflow[k] > 0)
              {
                Inflow[k]--;
                if(Inflow[k] == 0)
                {
                  Numbers[Index] = k;
                  Index++;
                }
              }
            }
          }
        } // endfor l
      } // endif
    } // endfor j
    State++;
  } // endwhile

  for(i=0;i<N_Active;i++)
  {
    cout << i << "   " << Numbers[i] << endl;
  }

  delete Inflow;
}

/** write matrix into file */
int TSquareMatrix::Write(const char *filename)
{
  int header[3];

  std::ofstream dat(filename);
  if(!dat)
  {
    cerr << "cannot open file '" << filename << "' for output" << endl;
    return -1;
  }

  header[0] = this->structure->GetN_Rows();
  header[1] = this->structure->GetN_Columns();
  header[2] = this->structure->GetN_Entries();

  dat.write((char *)header, sizeof(int)*3);
  dat.write((char *)this->structure->GetRowPtr(), sizeof(int)*(this->structure->GetN_Rows()+1));
  dat.write((char *)this->structure->GetKCol(), sizeof(int)*this->structure->GetN_Entries());
  dat.write((char *)Entries, sizeof(double)*this->structure->GetN_Entries());

  dat.close();
  
  return 0;
}

void TSquareMatrix::Print()
{
  int begin, end, pos=0;
  
  for (int i=0;i<this->structure->GetN_Rows();++i)
  {
    begin = this->structure->GetRowPtr()[i];
    end   = this->structure->GetRowPtr()[i+1];
    
    for (int j=begin;j<end;++j)
    {
      cout << "a(" << i << "," << this->structure->GetKCol()[pos] << ") = " << Entries[pos] << endl;
      ++pos;
    }
  }
}
