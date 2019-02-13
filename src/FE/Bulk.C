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
   
/****************************************************************************************
 *                                                                                      *
 *                         Bulk.C                                                       *
 *                        -------                                                       *
 *                                                                                      *
 *  common routines for simulation of precipitation process in 2d and 3d cavity         *
 *                                                                                      *
 ***************************************************************************************/
#include <Database.h>
#include <LinAlg.h>
#include <string.h>

#ifdef __2D__
#include <SquareStructure2D.h>
#include <SquareMatrix2D.h>
#include <FEFunction2D.h>
#endif    

#ifdef __3D__
#include <SquareStructure3D.h>
#include <SquareMatrix3D.h>
#include <FEFunction3D.h>
#endif    

// ======================================================================
// lump matrix to diagonal matrix
// the sparsity pattern of the matrix is not condensed
// ======================================================================

#ifdef __2D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix2D *M)
#endif    
#ifdef __3D__
void LumpMassMatrixToDiag_Bulk(TSquareMatrix3D *M)
#endif    
{
  double *Entries, off_diag;
  int *RowPtr, *KCol, i, j, rows, j0, j1, diag;

  RowPtr        = M->GetRowPtr();
  KCol          = M->GetKCol();
  Entries       = M->GetEntries();
  rows          = M->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    off_diag = 0;
    for (j=j0;j<j1;j++)
    {
      // diagonal entry
      if (KCol[j] == i)
        diag = j;
      else
      {
          //if (TDatabase::ParamDB->BULK_REACTION_MASS_LUMPING==1)
	  off_diag +=  Entries[j];
	  Entries[j] = 0;
      }
    }
    Entries[diag] += off_diag;
  }
}

/****************************************************************************************
 *                                                                                       *
 *  Computation of dp_50
 *                                                                                       *
 ****************************************************************************************/

double calculate_dp_50(int N, double *size, double *number)
{
  double a, b, a_0, a_1, dp_0, dp_1;
  double dp_3, val_left, val_right;
  double *Q3;
  double dp_50 = 0.0;

  Q3 = new double[N];
  number[0] = 1e-30;

  // computation of q3 with the trapezoidal rule
  double integral = 0.0, integral_Q3 = 0.0;

  for ( int i=0 ; i<(N-1) ; i++ )
  {
    val_left  = size[i]*size[i]*size[i]*number[i];
    val_right = size[i+1]*size[i+1]*size[i+1]*number[i+1];
    integral = integral + (val_left+val_right)*(size[i+1]-size[i])*0.5;
  }

  for ( int i=0 ; i<N ; i++ )
  {
    dp_3 = size[i]*size[i]*size[i];
    number[i] = number[i] * dp_3 / integral;
  }

  // computation of dp of 0.5 and Q3 (again trapezoidal rule)
  Q3[0] = 0.0;

  for ( int i=0 ; i<(N-1) ; i++ )
  {
    val_left  = number[i];
    val_right = number[i+1];
    integral_Q3 = integral_Q3 + (val_left+val_right) * (size[i+1]-size[i]) * 0.5;
    Q3[i+1] = integral_Q3;
  }

  for ( int i=1 ; i<(N-1) ; i++ )
  {
    if ( (Q3[i]<0.5) && (Q3[i+1]>=0.5) )
    {
      a_0 = Q3[i];
      a_1 = Q3[i+1];
      dp_0 = size[i];
      dp_1 = size[i+1];
      a = (a_1-a_0)/(dp_1-dp_0);
      b = a_0 - dp_0*a;
      dp_50 = (0.5-b)/a;
    }
  }
  delete Q3;
  return dp_50;
}
