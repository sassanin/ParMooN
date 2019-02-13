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
// VMS.C
//
// Purpose:     routines for projection-based VMS
//
// Author:       Volker John  2006/05/18
//
// =======================================================================

#include <Database.h>
#include <MooNMD_Io.h>
#include <Enumerations.h>
#ifdef __2D__
#include <FEDatabase2D.h>
#else
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <FEVectFunct3D.h>
#endif
#include <LinAlg.h>
#include <VMS.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/******************************************************************************/
//
// FULLY IMPLICITE PROJECTION-BASED VMS
//
/******************************************************************************/

#ifdef __2D__

void VMSProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
TSquareMatrix2D **SQMATRICES, TMatrix2D **MATRICES)
{
  TSquareMatrix2D *MatrixA11, *MatrixA22, *MatrixA12, *MatrixA21, *MatrixL;
  TMatrix2D *Matrix_tilde_G11,  *Matrix_tilde_G22;
  TMatrix2D *Matrix_G11,  *Matrix_G22;
  double *val, *Entries_tilde, *Entries, *Entries_A, *Entries_A1, *Entries_L;
  double val1, val2;
  int i, j0, j1, j, k, l0, l1, l, index, m;
  int *RowPtr_tilde, *KCol_tilde, *RowPtr, *KCol, *RowPtr_A, *KCol_A;
  int *RowPtr_L, *KCol_L;

  val = new double[N_U];
  memset(val,0,N_U*SizeOfDouble);

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatrixA11 = SQMATRICES[0];
      MatrixL = SQMATRICES[1];
      Matrix_tilde_G11 = MATRICES[2];
      Matrix_tilde_G22 = MATRICES[3];
      Matrix_G11 = MATRICES[4];
      Matrix_G22 = MATRICES[5];

      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_L     = MatrixL->GetEntries();
      RowPtr_L      = MatrixL->GetRowPtr();
      KCol_L        = MatrixL->GetKCol();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }
        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          val[index] = 0.0;
        }
      }

      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }
        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          val[index] = 0.0;
        }
      }
      break;
    case 3:
    case 4:
      MatrixA11 = SQMATRICES[0];
      MatrixA12 = SQMATRICES[1];
      MatrixA21 = SQMATRICES[2];
      MatrixA22 = SQMATRICES[3];
      MatrixL = SQMATRICES[6];
      Matrix_tilde_G11 = MATRICES[2];
      Matrix_tilde_G22 = MATRICES[3];
      Matrix_G11 = MATRICES[4];
      Matrix_G22 = MATRICES[5];

      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_A1    = MatrixA22->GetEntries();
      Entries_L     = MatrixL->GetEntries();
      RowPtr_L      = MatrixL->GetRowPtr();
      KCol_L        = MatrixL->GetKCol();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          Entries_A1[j] -= val[index]/2;
          val[index] = 0.0;
        }
        // memset(val,0,N_U*SizeOfDouble);
      }

      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_A1    = MatrixA22->GetEntries();
      Entries_L     = MatrixL->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/(val1);
            //                 OutPut(" " <<  Entries[l]);
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          Entries_A1[j] -= val[index];
          //             OutPut(" " << val[index]);
          val[index] = 0.0;
        }
        // memset(val,0,N_U*SizeOfDouble);
      }

      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA12->GetRowPtr();
      KCol_A        = MatrixA12->GetKCol();
      Entries_A     = MatrixA12->GetEntries();
      Entries_L     = MatrixL->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/(val1);
            //                 OutPut(" " <<  Entries[l]);
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          //              OutPut(" " << val[index]);
          val[index] = 0.0;
        }
        // memset(val,0,N_U*SizeOfDouble);
      }

      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();
      RowPtr_A      = MatrixA21->GetRowPtr();
      KCol_A        = MatrixA21->GetKCol();
      Entries_A     = MatrixA21->GetEntries();
      Entries_L     = MatrixL->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/(val1);
            //                 OutPut(" " <<  Entries[l]);
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          //              OutPut(" " << val[index]);
          val[index] = 0.0;
        }
        // memset(val,0,N_U*SizeOfDouble);
      }

      break;
    default:
      OutPut("VMSProjectionUpdateMatrices not implemented !!" << endl);
      exit(4711);
  }
  delete(val);
}


// ======================================================================
// lump matrix to diagonal matrix
// ======================================================================

void LumpMassMatrixToDiag(TSquareMatrix2D *M)
{
  double *Entries;
  int *RowPtr, *KCol, i, j, rows, j0, j1;

  RowPtr        = M->GetRowPtr();
  KCol          = M->GetKCol();
  Entries       = M->GetEntries();
  rows          = M->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for (j=j0;j<j1;j++)
    {
      // diagonal entry
      if (KCol[j] == i)
      {
        RowPtr[i] = i;
        Entries[i] =  Entries[j];
        KCol[i] = i;
        break;
      }
    }
    // for space P00, Q00
    if (fabs(Entries[i])<1e-10)
      Entries[i] = 1.0;
  }
  RowPtr[rows] = rows;
}
#endif                                            // __2D__
#ifdef __3D__

void VMS_ProjectionUpdateMatrices(int N_U,int N_Active,int N_L,
TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES)
{
  TSquareMatrix3D *MatrixA11, *MatrixA12, *MatrixA13, *MatrixA21, *MatrixA22;
  TSquareMatrix3D *MatrixA23, *MatrixA31, *MatrixA32, *MatrixA33, *MatrixL;
  TMatrix3D *Matrix_tilde_G11,  *Matrix_tilde_G22,  *Matrix_tilde_G33;
  TMatrix3D *Matrix_G11,  *Matrix_G22,  *Matrix_G33;
  double *val, *Entries_tilde, *Entries, *Entries_A, *Entries_A1, *Entries_L;
  double *Entries_A2, val1, val2;
  int i, j0, j1, j, k, l0, l1, l, index, m, index1;
  int *RowPtr_tilde, *KCol_tilde, *RowPtr, *KCol, *RowPtr_A, *KCol_A, *RowPtr_L, *KCol_L;

  val = new double[N_U];
  memset(val,0,N_U*SizeOfDouble);

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 1:
      MatrixA11 = SQMATRICES[0];
      MatrixL = SQMATRICES[1];
      Matrix_tilde_G11 = MATRICES[3];
      Matrix_tilde_G22 = MATRICES[4];
      Matrix_tilde_G33 = MATRICES[5];
      Matrix_G11 = MATRICES[6];
      Matrix_G22 = MATRICES[7];
      Matrix_G33 = MATRICES[8];

      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_L     = MatrixL->GetEntries();
      RowPtr_L      = MatrixL->GetRowPtr();
      KCol_L        = MatrixL->GetKCol();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }
        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          val[index] = 0.0;
        }
      }

      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }
        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          val[index] = 0.0;
        }
      }
      RowPtr_tilde  = Matrix_tilde_G33->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G33->GetKCol();
      Entries_tilde = Matrix_tilde_G33->GetEntries();
      RowPtr        = Matrix_G33->GetRowPtr();
      KCol          = Matrix_G33->GetKCol();
      Entries       = Matrix_G33->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr_L[index];
          l1 = RowPtr_L[index+1];
          for (l=l0;l<l1;l++)
          {
            // diagonal entry
            if (KCol_L[l]==index)
            {
              val1 = Entries_L[l];
              break;
            }
          }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          for (l=l0;l<l1;l++)
          {
            val[KCol[l]] += Entries_tilde[j] * Entries[l]/val1;
          }
        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          val[index] = 0.0;
        }
      }
      break;
    case 3:
    case 4:
      MatrixA11 = SQMATRICES[0];
      MatrixA12 = SQMATRICES[1];
      MatrixA13 = SQMATRICES[2];
      MatrixA21 = SQMATRICES[3];
      MatrixA22 = SQMATRICES[4];
      MatrixA23 = SQMATRICES[5];
      MatrixA31 = SQMATRICES[6];
      MatrixA32 = SQMATRICES[7];
      MatrixA33 = SQMATRICES[8];
      MatrixL = SQMATRICES[9];
      Matrix_tilde_G11 = MATRICES[0];
      Matrix_tilde_G22 = MATRICES[1];
      Matrix_tilde_G33 = MATRICES[2];
      Matrix_G11 = MATRICES[3];
      Matrix_G22 = MATRICES[4];
      Matrix_G33 = MATRICES[5];
      // OutPut("u1"<<endl);
      // multiplication of tilde_G11 L^{-1} G11
      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_A1    = MatrixA22->GetEntries();
      Entries_A2    = MatrixA33->GetEntries();
      Entries_L     = MatrixL->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 =  Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index];
          Entries_A1[j] -= val[index]/2;
          Entries_A2[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }

      //OutPut("u2"<<endl);
      // multiplication of tilde_G22 L^{-1} G22
      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_A1    = MatrixA22->GetEntries();
      Entries_A2    = MatrixA33->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G22
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          Entries_A1[j] -= val[index];
          Entries_A2[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }

      //  OutPut("u3"<<endl);
      // multiplication of tilde_G33 L^{-1} G33
      RowPtr_tilde  = Matrix_tilde_G33->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G33->GetKCol();
      Entries_tilde = Matrix_tilde_G33->GetEntries();
      RowPtr        = Matrix_G33->GetRowPtr();
      KCol          = Matrix_G33->GetKCol();
      Entries       = Matrix_G33->GetEntries();
      RowPtr_A      = MatrixA11->GetRowPtr();
      KCol_A        = MatrixA11->GetKCol();
      Entries_A     = MatrixA11->GetEntries();
      Entries_A1    = MatrixA22->GetEntries();
      Entries_A2    = MatrixA33->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G22
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          //l0 = RowPtr_L[index];
          //l1 = RowPtr_L[index+1];
          //for (l=l0;l<l1;l++)
          // {
          // diagonal entry
          //  if (KCol_L[l]==index)
          //  {
          //     val1 = Entries_L[l];
          //     break;
          //  }
          // }
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          Entries_A1[j] -= val[index]/2;
          Entries_A2[j] -= val[index];
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      // OutPut("u4"<<endl);
      // multiplication of tilde_G22 L^{-1} G11
      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA12->GetRowPtr();
      KCol_A        = MatrixA12->GetKCol();
      Entries_A     = MatrixA12->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 =  Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      //OutPut("u5"<<endl);

      // multiplication of tilde_G33 L^{-1} G11
      RowPtr_tilde  = Matrix_tilde_G33->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G33->GetKCol();
      Entries_tilde = Matrix_tilde_G33->GetEntries();
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();
      RowPtr_A      = MatrixA13->GetRowPtr();
      KCol_A        = MatrixA13->GetKCol();
      Entries_A     = MatrixA13->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G11
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      //OutPut("u6"<<endl);

      // multiplication of tilde_G11 L^{-1} G22
      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();
      RowPtr_A      = MatrixA21->GetRowPtr();
      KCol_A        = MatrixA21->GetKCol();
      Entries_A     = MatrixA21->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j]/Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      //OutPut("u7"<<endl);

      // multiplication of tilde_G33 L^{-1} G22
      RowPtr_tilde  = Matrix_tilde_G33->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G33->GetKCol();
      Entries_tilde = Matrix_tilde_G33->GetEntries();
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();
      RowPtr_A      = MatrixA23->GetRowPtr();
      KCol_A        = MatrixA23->GetKCol();
      Entries_A     = MatrixA23->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j] / Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      //OutPut("u8"<<endl);
      // multiplication of tilde_G11 L^{-1} G33
      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();
      RowPtr        = Matrix_G33->GetRowPtr();
      KCol          = Matrix_G33->GetKCol();
      Entries       = Matrix_G33->GetEntries();
      RowPtr_A      = MatrixA31->GetRowPtr();
      KCol_A        = MatrixA31->GetKCol();
      Entries_A     = MatrixA31->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j] / Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      //OutPut("u9"<<endl);
      // multiplication of tilde_G22 L^{-1} G33
      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();
      RowPtr        = Matrix_G33->GetRowPtr();
      KCol          = Matrix_G33->GetKCol();
      Entries       = Matrix_G33->GetEntries();
      RowPtr_A      = MatrixA32->GetRowPtr();
      KCol_A        = MatrixA32->GetKCol();
      Entries_A     = MatrixA32->GetEntries();

      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          // this is the row of Matrix_G22
          index = KCol_tilde[j];
          l0 = RowPtr[index];
          l1 = RowPtr[index+1];
          val1 = Entries_tilde[j] / Entries_L[index];
          for (l=l0;l<l1;l++)
          {
            index1 = KCol[l];
            val[index1] += val1 * Entries[l];
          }

        }
        // add to i-th row of A11
        j0 = RowPtr_A[i];
        j1 = RowPtr_A[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_A[j];
          Entries_A[j] -= val[index]/2;
          val[index] = 0.0;
        }
        //memset(val,0,N_U*SizeOfDouble);
      }
      break;

    default:
      OutPut("VMS_ProjectionUpdateMatrices not implemented !!" << endl);
      exit(4711);
  }

  /*  int end, *ColInd, N_Rows;
    for(k=0;k<8;k++)
    {
      OutPut(endl);
      cout << "sqmatrix: " << k << endl;
      RowPtr = SQMATRICES[k]->GetRowPtr();
      Entries = SQMATRICES[k]->GetEntries();
      ColInd = SQMATRICES[k]->GetKCol();
      N_Rows = SQMATRICES[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  OutPut("ta"<<k<<"("<< i+1 << "," << ColInd[j]+1 << " )=  ");
  OutPut(Entries[j] << ";"<<endl);
  }
  }
  cout << endl;
  } // endfor k
  exit(1);
  */

  delete(val);
}


// ======================================================================
// explicit projection-based VMS
// ======================================================================

void VMS_ProjectionExplUpdateRhs(int N_U,int N_Active,int N_L, TFEVectFunct3D *u,
TSquareMatrix3D **SQMATRICES, TMatrix3D **MATRICES, double *rhs_vms_expl)
{
  double *u1, *u2, *u3;
  TSquareMatrix3D *MatrixL;
  TMatrix3D *Matrix_tilde_G11,  *Matrix_tilde_G22,  *Matrix_tilde_G33;
  TMatrix3D *Matrix_G11,  *Matrix_G22,  *Matrix_G33;
  double *Entries_tilde, *Entries, *Entries_L;
  double val1;
  double *x, *y, *z;
  int i, j0, j1, j, index;
  int *RowPtr_tilde, *KCol_tilde, *RowPtr, *KCol;

  u1 = u->GetComponent(0)->GetValues();
  u2 = u->GetComponent(1)->GetValues();
  u3 = u->GetComponent(2)->GetValues();

  memset(rhs_vms_expl,0,3*N_U*SizeOfDouble);

  x = new double[3*N_L];
  y = x + N_L;
  z = y + N_L;
  memset(x,0,3*N_L*SizeOfDouble);

  switch(TDatabase::ParamDB->NSTYPE)
  {
    case 3:
    case 4:
      MatrixL = SQMATRICES[9];
      Matrix_tilde_G11 = MATRICES[0];
      Matrix_tilde_G22 = MATRICES[1];
      Matrix_tilde_G33 = MATRICES[2];
      Matrix_G11 = MATRICES[3];
      Matrix_G22 = MATRICES[4];
      Matrix_G33 = MATRICES[5];

      // matrix L
      Entries_L = MatrixL->GetEntries();

      // matrix G11
      RowPtr        = Matrix_G11->GetRowPtr();
      KCol          = Matrix_G11->GetKCol();
      Entries       = Matrix_G11->GetEntries();

      // multiplication of G11 u1=:x, G11 u2=:y, G11 u3=:z
      for(i=0;i<N_L;i++)
      {
        // i-th row of G11
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol[j];
          x[i] += Entries[j]*u1[index];
          y[i] += Entries[j]*u2[index];
          z[i] += Entries[j]*u3[index];
        }
        // multiplication of L^{-1} x =:x
        val1 = Entries_L[i];
        x[i] /= val1;
      }

      // matrix tilde_G11
      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();

      // multiplication of tilde_G11 x
      // and add to the right-hand side term 1
      for(i=0;i<N_Active;i++)
      {
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[i] += Entries_tilde[j]*x[index];
        }
      }

      // matrix G22
      RowPtr        = Matrix_G22->GetRowPtr();
      KCol          = Matrix_G22->GetKCol();
      Entries       = Matrix_G22->GetEntries();

      // multiplication of G22 u1 and add to y =:y, G22 u2=:x
      for(i=0;i<N_L;i++)
      {
        // free vector x
        x[i] = 0;
        // i-th row of G22
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol[j];
          y[i] += Entries[j]*u1[index];
          x[i] += Entries[j]*u2[index];
        }
        // multiplication of L^{-1} x =:x, 1/ 2 L^{-1} y =:y
        val1 = Entries_L[i];
        x[i] /= val1;
        y[i] /= (2*val1);
      }

      // multiplication of tilde_G11 y
      // and add to the right-hand side term 2
      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[N_U+i] += Entries_tilde[j]*y[index];
        }
      }

      // matrix tilde_G22
      RowPtr_tilde  = Matrix_tilde_G22->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G22->GetKCol();
      Entries_tilde = Matrix_tilde_G22->GetEntries();

      // multiplication of tilde_G22 y and add to right-hand side term 1,
      // tilde_G22 x and add to right-hand side term 2
      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G22
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[i] += Entries_tilde[j]*y[index];
          rhs_vms_expl[N_U+i] += Entries_tilde[j]*x[index];
        }
      }

      // multiplication of G22 u3=:y
      for(i=0;i<N_L;i++)
      {
        // free vector y
        y[i] = 0;
        // i-th row of  Matrix_G22
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol[j];
          y[i] += Entries[j]*u3[index];
        }
      }

      // matrix G33
      RowPtr        = Matrix_G33->GetRowPtr();
      KCol          = Matrix_G33->GetKCol();
      Entries       = Matrix_G33->GetEntries();

      // multiplication of G33 u2 and add to y =:y,
      // and G33 u3=:x, G33 u1 and add to z=:z
      for(i=0;i<N_L;i++)
      {
        // free vector x
        x[i] = 0;
        // i-th row of  Matrix_G33
        j0 = RowPtr[i];
        j1 = RowPtr[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol[j];
          z[i] += Entries[j]*u1[index];
          y[i] += Entries[j]*u2[index];
          x[i] += Entries[j]*u3[index];
        }
        // multiplication of L^{-1} x =:x, 1/2 L^{-1} y =:y, 1/2 L^{-1} z =:z,
        val1 = Entries_L[i];
        x[i] /= val1;
        y[i] /= (2*val1);
        z[i] /= (2*val1);
      }
      // multiplication of tilde_G22 y
      // and add to right-hand side term 3
      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G33
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[2*N_U+i] += Entries_tilde[j]*y[index];
        }
      }

      // matrix tilde_G33
      RowPtr_tilde  = Matrix_tilde_G33->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G33->GetKCol();
      Entries_tilde = Matrix_tilde_G33->GetEntries();

      // multiplication of tilde_G33 y and add to right-hand side term 2,
      // tilde_G33 x and add to right-hand side term 3,
      // and tilde_G33 z and add to right-hand side term 1
      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G33
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[N_U+i] += Entries_tilde[j]*y[index];
          rhs_vms_expl[2*N_U+i] += Entries_tilde[j]*x[index];
          rhs_vms_expl[i] += Entries_tilde[j]*z[index];
        }
      }

      // Matrix tilde_G11
      RowPtr_tilde  = Matrix_tilde_G11->GetRowPtr();
      KCol_tilde    = Matrix_tilde_G11->GetKCol();
      Entries_tilde = Matrix_tilde_G11->GetEntries();

      // multiplication of tilde_G11 z and add to right-hand side term 3
      for(i=0;i<N_Active;i++)
      {
        // i-th row of  Matrix_tilde_G11
        j0 = RowPtr_tilde[i];
        j1 = RowPtr_tilde[i+1];
        for(j=j0;j<j1;j++)
        {
          index = KCol_tilde[j];
          rhs_vms_expl[2*N_U+i] += Entries_tilde[j]*z[index];
        }
      }

      break;

    default:
      OutPut("VMS_ProjectionExplUpdateMatrices not implemented !!" << endl);
      exit(4711);
  }
  /*  int end, *ColInd, N_Rows;
    for(k=0;k<8;k++)
    {
      OutPut(endl);
      cout << "sqmatrix: " << k << endl;
      RowPtr = SQMATRICES[k]->GetRowPtr();
      Entries = SQMATRICES[k]->GetEntries();
      ColInd = SQMATRICES[k]->GetKCol();
      N_Rows = SQMATRICES[k]->GetN_Rows();
      for(i=0;i<N_Rows;i++)
      {
  end=RowPtr[i+1];
  for(j=RowPtr[i];j<end;j++)
  {
  // cout << j << endl;
  OutPut("ta"<<k<<"("<< i+1 << "," << ColInd[j]+1 << " )=  ");
  OutPut(Entries[j] << ";"<<endl);
  }
  }
  cout << endl;
  } // endfor k
  exit(1);
  */

  delete x;
}


// ======================================================================
// lump matrix to diagonal matrix
// ======================================================================

void LumpMassMatrixToDiag(TSquareMatrix3D *M)
{
  double *Entries;
  int *RowPtr, *KCol, i, j, rows, j0, j1;
  int check;
  RowPtr        = M->GetRowPtr();
  KCol          = M->GetKCol();
  Entries       = M->GetEntries();
  rows          = M->GetN_Rows();

  for (i=0; i<rows; i++)
  {
    check = 0;
    j0 = RowPtr[i];
    j1 = RowPtr[i+1];
    for (j=j0;j<j1;j++)
    {
      // diagonal entry
      if (KCol[j] == i)
      {
        RowPtr[i] = i;
        Entries[i] =  Entries[j];
        KCol[i] = i;
	check=1;
	break;
      }
      if(check)
	break;
//       else
//       {
// 	cout << "error here" << i << endl;
// 	exit(0);
//       }
	
    }
    // for space Q00
    if (fabs(Entries[i])<1e-10)
      Entries[i] = 1.0;
  }
  RowPtr[rows] = rows;
}


// ======================================================================
// compute the VMS projection
// checked 07/03/07
// ======================================================================

void ComputeVMSProjection(TMatrix3D *matG11, TMatrix3D *matG22, TMatrix3D *matG33,
TSquareMatrix3D *MatrixL, TFEFunction3D *u_1,
TFEFunction3D *u_2, TFEFunction3D *u_3,
TFEVectFunct3D *vms_projection_fe)
{
  double *proj, *u1, *u2, *u3, *mass, *val;
  int *LRowPtr, *LKCol, N_LDOF, i, j;

  // get values of the projection
  proj = vms_projection_fe->GetValues();
  // get velocity field
  u1 = u_1->GetValues();
  u2 = u_2->GetValues();
  u3 = u_3->GetValues();
  // get mass matrix
  mass = MatrixL->GetEntries();
  LRowPtr = MatrixL->GetRowPtr();
  LKCol = MatrixL->GetKCol();
  N_LDOF = MatrixL->GetN_Rows();
  // temporary array
  val = new double[N_LDOF];

  // diagonals of tensor
  MatVect1(matG11, u1, proj);
  MatVect1(matG22, u2, proj+3*N_LDOF);
  MatVect1(matG33, u3, proj+5*N_LDOF);
  // off diagonals of tensor
  MatVect1(matG11, u2, proj+N_LDOF);
  MatVect1(matG22, u1, val);
  Daxpy(N_LDOF, 1, val, proj+N_LDOF);
  MatVect1(matG11, u3, proj+2*N_LDOF);
  MatVect1(matG33, u1, val);
  Daxpy(N_LDOF, 1, val, proj+2*N_LDOF);
  MatVect1(matG22, u3, proj+4*N_LDOF);
  MatVect1(matG33, u2, val);
  Daxpy(N_LDOF, 1, val, proj+4*N_LDOF);

  // diagonal matrix
  if (LRowPtr[N_LDOF] == N_LDOF)
  {
    for (i=0;i<6;i++)
    {
      for (j=0;j<N_LDOF;j++)
        proj[j+i*N_LDOF] /= -mass[j];
    }
  }
  else
  {
    OutPut("NOT A DIAGONAL MATRIX " << endl);
    exit(4711);
  }

  delete val;
}

/******************************************************************************/
//
// compute sizes of resolved small scales in L^2(K)
//
/******************************************************************************/
void ComputeSizeOfSmallScales(TMatrix3D *matG11, 
			      TMatrix3D *matG22, 
			      TMatrix3D *matG33,
			      TSquareMatrix3D *MatrixL, 
			      TFEFunction3D *u1,
			      TFEFunction3D *u2, 
			      TFEFunction3D *u3,
			      TFEVectFunct3D *vms_projection_fe, 
			      double *size_small_scales)
{
  int i, j, k, l, N_Cells;
  int N_UsedElements_velo, N_LocalUsedElements_velo, N_Points_velo, N_velo;
  int N_UsedElements_LS, N_LocalUsedElements_LS, N_Points_LS, N_LS;
  int Used[N_FEs3D], *N_BaseFunct, *DOF_velo, *DOF_LS;
  int *GlobalNumbers_velo, *BeginIndex_velo;
  int *GlobalNumbers_LS, *BeginIndex_LS;
  double *u1_values, *u2_values, *u3_values, u1_value, u2_value, u3_value;
  double *vms_proj_11_values, *vms_proj_12_values, *vms_proj_13_values;
  double *vms_proj_22_values, *vms_proj_23_values, *vms_proj_33_values;
  double vms_proj_11_value, vms_proj_12_value, vms_proj_13_value;
  double vms_proj_22_value, vms_proj_23_value, vms_proj_33_value;
  double *Derivatives_u1[MaxN_QuadPoints_3D],*Derivatives_u2[MaxN_QuadPoints_3D];
  double *Derivatives_u3[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_11[MaxN_QuadPoints_3D], *Derivatives_vms_proj_12[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_13[MaxN_QuadPoints_3D], *Derivatives_vms_proj_22[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_23[MaxN_QuadPoints_3D], *Derivatives_vms_proj_33[MaxN_QuadPoints_3D];
  double u1_FEFunctValues[MaxN_BaseFunctions3D], u2_FEFunctValues[MaxN_BaseFunctions3D];
  double u3_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_11_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_12_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_13_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_22_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_23_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_33_FEFunctValues[MaxN_BaseFunctions3D];
  double *weights_velo, *xi_velo, *eta_velo, *zeta_velo;
  double *weights_LS, *xi_LS, *eta_LS, *zeta_LS;
  double X_velo[MaxN_QuadPoints_3D], Y_velo[MaxN_QuadPoints_3D], Z_velo[MaxN_QuadPoints_3D];
  double X_LS[MaxN_QuadPoints_3D], Y_LS[MaxN_QuadPoints_3D], Z_LS[MaxN_QuadPoints_3D];
  double AbsDetjk_velo[MaxN_QuadPoints_3D], AbsDetjk_LS[MaxN_QuadPoints_3D];
  double **OrigFEValues_velo, *Orig_velo;
  double **OrigFEValues_LS, *Orig_LS;
  double hK, K,*aux, *aux1, t, w, val, small_scales = 0, all_scales = 0, all_scales_loc;
  TFEFunction3D *vms_proj_11, *vms_proj_12,*vms_proj_13;
  TFEFunction3D *vms_proj_22, *vms_proj_23,*vms_proj_33;
  MultiIndex3D VeloDerivatives[3] = { D100, D010, D001 };
  MultiIndex3D LargeScaleDerivatives[1] = { D000 };
  TFESpace3D *fespaces[2];
  FE3D LocalUsedElements_velo[N_FEs3D], CurrentElement_velo;
  FE3D LocalUsedElements_LS[N_FEs3D], CurrentElement_LS;
  BaseFunct3D BaseFunct_velo, BaseFunct_LS, *BaseFuncts;
  TBaseCell *cell;
  TCollection *Coll;
  bool SecondDer[2];
  // set quadrature formula to Gauss3 
  // since velocity and large scale spaces have different 
  // quadrature formulas by default
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 2;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<2;i++)
    SecondDer[i] = FALSE;

  // this gives the large scales
  ComputeVMSProjection(matG11, matG22, matG33, MatrixL,
		       u1, u2, u3, vms_projection_fe);

  // velocity
  fespaces[0] = u1->GetFESpace3D();
  GlobalNumbers_velo = fespaces[0]->GetGlobalNumbers();
  BeginIndex_velo = fespaces[0]->GetBeginIndex();

  u1_values = u1->GetValues();
  u2_values = u2->GetValues();
  u3_values = u3->GetValues();

  aux = new double[MaxN_QuadPoints_3D*15];
  aux1 = aux;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u1[j] = aux + j*3;
  aux += 3*MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u2[j] = aux + j*3;
  aux += 3*MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u3[j] = aux + j*3;

  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_11[j] = aux+j;
  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_12[j] = aux+j;
  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_13[j] = aux+j;
  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_22[j] = aux+j;
  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_23[j] = aux+j;
  aux += MaxN_QuadPoints_3D;
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_33[j] = aux+j;

  // large scale tensor
  vms_proj_11 = vms_projection_fe->GetComponent(0);
  vms_proj_12 = vms_projection_fe->GetComponent(1);
  vms_proj_13 = vms_projection_fe->GetComponent(2);
  vms_proj_22 = vms_projection_fe->GetComponent(3);
  vms_proj_23 = vms_projection_fe->GetComponent(4);
  vms_proj_33 = vms_projection_fe->GetComponent(5);

  fespaces[1] = vms_proj_11->GetFESpace3D();
  GlobalNumbers_LS = fespaces[1]->GetGlobalNumbers();
  BeginIndex_LS = fespaces[1]->GetBeginIndex();

  vms_proj_11_values = vms_proj_11->GetValues();
  vms_proj_12_values = vms_proj_12->GetValues();
  vms_proj_13_values = vms_proj_13->GetValues();
  vms_proj_22_values = vms_proj_22->GetValues();
  vms_proj_23_values = vms_proj_23->GetValues();
  vms_proj_33_values = vms_proj_33->GetValues();

  // collection 
  Coll = u1->GetFESpace3D()->GetCollection();
  // number of mesh cells 
  N_Cells = Coll->GetN_Cells();                   
  memset(size_small_scales,0,N_Cells*SizeOfDouble);

  // loop over all mesh cells
  for (i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();
    K = cell->GetMeasure();
    //OutPut("K="<<sqrt(K)<<"  ");
    // ******************************************************
    // VELOCITY
    // find local used elements on this cell
    // ******************************************************
    memset(Used, 0, N_FEs3D*SizeOfInt);
    CurrentElement_velo = fespaces[0]->GetFE3D(i, cell);
    Used[CurrentElement_velo] = 1;

    // compute number of used finite elements
    N_LocalUsedElements_velo = 0;
    memset(LocalUsedElements_velo, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
    {
      LocalUsedElements_velo[j] = (FE3D)k;
      j++;
    }
    // result should be 1
    N_LocalUsedElements_velo = j;

    // *************************************************
    // calculate values on original element
    // *************************************************
    TFEDatabase3D::GetOrig(N_LocalUsedElements_velo, LocalUsedElements_velo,
      Coll, cell, SecondDer,
      N_Points_velo, xi_velo, eta_velo, zeta_velo, weights_velo,
      X_velo, Y_velo, Z_velo, AbsDetjk_velo);
    //OutPut("i: "<< i << "  N_Points_velo: "<< N_Points_velo << endl);
    //OutPut("xi_v: " << xi_velo[0] << "  eta_v: "<< eta_velo[0] << "  zeta_v: "<<
    //         zeta_velo[0] << endl);
    //OutPut("weights_velo: " << weights_velo[0] << endl);
    //OutPut("X_v: " << X_velo[0] << "  Y_v: " << Y_velo[0] <<"  Z_v: " << Z_velo[0] << endl);
    //OutPut("AbsDetjk_v: " << AbsDetjk_velo[0] << endl);
    //OutPut(endl);

    // calculate all needed derivatives of this FE function
    BaseFunct_velo = BaseFuncts[CurrentElement_velo];
    N_velo = N_BaseFunct[CurrentElement_velo];
    // find global dof for mesh cell i
    DOF_velo = GlobalNumbers_velo + BeginIndex_velo[i];
    // get coefficients of vectors
    for(l=0;l<N_velo;l++)
    {
      u1_FEFunctValues[l] = u1_values[DOF_velo[l]];
      u2_FEFunctValues[l] = u2_values[DOF_velo[l]];
      u3_FEFunctValues[l] = u3_values[DOF_velo[l]];
    }
    // compute values and gradient of velocity
    for(k=0;k<3;k++)
    {
      OrigFEValues_velo = TFEDatabase3D::GetOrigElementValues(BaseFunct_velo,
        VeloDerivatives[k]);
      //OutPut("k " << k << " ");
      for(j=0;j<N_Points_velo;j++)
      {
        Orig_velo = OrigFEValues_velo[j];
        u1_value = 0;
        u2_value = 0;
        u3_value = 0;
        for(l=0;l<N_velo;l++)
        {
          u1_value += u1_FEFunctValues[l] * Orig_velo[l];
	  //OutPut(l << " " << u1_FEFunctValues[l] << " " << Orig_velo[l] << " " << u1_value << endl);
          u2_value += u2_FEFunctValues[l] * Orig_velo[l];
          u3_value += u3_FEFunctValues[l] * Orig_velo[l];
        }                                         // endfor l
	//OutPut(endl);
        Derivatives_u1[j][k] = u1_value;
        Derivatives_u2[j][k] = u2_value;
        Derivatives_u3[j][k] = u3_value;
      }                                           // endfor j
    }                                             // endfor k

    // ********************************************
    // LARGE SCALE TENSOR
    // find local used elements on this cell
    // ********************************************
    memset(Used, 0, N_FEs3D*SizeOfInt);
    CurrentElement_LS = fespaces[1]->GetFE3D(i, cell);
    Used[CurrentElement_LS] = 1;

    // compute number of used finite elements
    N_LocalUsedElements_LS = 0;
    memset(LocalUsedElements_LS, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
    {
      LocalUsedElements_LS[j] = (FE3D)k;
      j++;
    }
    // result should be 1
    N_LocalUsedElements_LS = j;

    // ********************************************
    // calculate values on original element
    // ********************************************
    TFEDatabase3D::GetOrig(N_LocalUsedElements_LS, LocalUsedElements_LS,
      Coll, cell, SecondDer,
      N_Points_LS, xi_LS, eta_LS, zeta_LS,
      weights_LS,
      X_LS, Y_LS, Z_LS, AbsDetjk_LS);

    //OutPut("N_Points_LS: "<< N_Points_LS << endl);
    //OutPut("xi_LS: " << xi_LS[0] << "  eta_LS: "<< eta_LS[0] << "  zeta_LS: "<<
    //         zeta_LS[0] << endl);
    //OutPut("weights_LS: " << weights_LS[0] << endl);
    // OutPut("X_LS: " << X_LS[0] << "  Y_LS: " << Y_LS[0] << "  Z_LS: " << Z_LS[0] << " ");
    //OutPut("AbsDetjk_LS: " << AbsDetjk_LS[0] << endl);
    //OutPut("**************" << endl);

    // calculate all needed derivatives of this FE function
    BaseFunct_LS = BaseFuncts[CurrentElement_LS];
    N_LS = N_BaseFunct[CurrentElement_LS];
    // find global dof for mesh cell i
    DOF_LS = GlobalNumbers_LS + BeginIndex_LS[i];
    // get coefficients of vectors
    for(l=0;l<N_LS;l++)
    {
      vms_proj_11_FEFunctValues[l] = vms_proj_11_values[DOF_LS[l]];
      vms_proj_12_FEFunctValues[l] = vms_proj_12_values[DOF_LS[l]];
      vms_proj_13_FEFunctValues[l] = vms_proj_13_values[DOF_LS[l]];
      vms_proj_22_FEFunctValues[l] = vms_proj_22_values[DOF_LS[l]];
      vms_proj_23_FEFunctValues[l] = vms_proj_23_values[DOF_LS[l]];
      vms_proj_33_FEFunctValues[l] = vms_proj_33_values[DOF_LS[l]];
    }

    OrigFEValues_LS = TFEDatabase3D::GetOrigElementValues(BaseFunct_LS,
      LargeScaleDerivatives[0]);

    // compute values of large scale tensor
    for(j=0;j<N_Points_LS;j++)
    {
      Orig_LS = OrigFEValues_LS[j];
      vms_proj_11_value = 0;
      vms_proj_12_value = 0;
      vms_proj_13_value = 0;
      vms_proj_22_value = 0;
      vms_proj_23_value = 0;
      vms_proj_33_value = 0;

      for(l=0;l<N_LS;l++)
      {
        vms_proj_11_value += vms_proj_11_FEFunctValues[l] * Orig_LS[l];
        vms_proj_12_value += vms_proj_12_FEFunctValues[l] * Orig_LS[l];
        vms_proj_13_value += vms_proj_13_FEFunctValues[l] * Orig_LS[l];
        vms_proj_22_value += vms_proj_22_FEFunctValues[l] * Orig_LS[l];
        vms_proj_23_value += vms_proj_23_FEFunctValues[l] * Orig_LS[l];
        vms_proj_33_value += vms_proj_33_FEFunctValues[l] * Orig_LS[l];
	//OutPut("13 " << l << " " << vms_proj_13_FEFunctValues[l] << " " << Orig_LS[l] << endl);
      }                                           // endfor l
      Derivatives_vms_proj_11[j][0] = vms_proj_11_value;
      Derivatives_vms_proj_12[j][0] = vms_proj_12_value/2.0;
      Derivatives_vms_proj_13[j][0] = vms_proj_13_value/2.0;
      Derivatives_vms_proj_22[j][0] = vms_proj_22_value;
      Derivatives_vms_proj_23[j][0] = vms_proj_23_value/2.0;
      Derivatives_vms_proj_33[j][0] = vms_proj_33_value;
    }                                             // endfor j

    /*  ErrorMeth(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives,
                 ExactVal, AuxArray, LocError);*/

    // the L2 norm
    for(j=0;j<N_Points_velo;j++)
    {
      // weight
      w = weights_velo[j]*AbsDetjk_velo[j];

      // term 11
      t = Derivatives_vms_proj_11[j][0] - Derivatives_u1[j][0];
      size_small_scales[i] += w*t*t;
      // terms 12 and 21
      t = Derivatives_vms_proj_12[j][0] - (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0;
      size_small_scales[i] += w*t*t*2;
      // terms 13 and 31
      t = Derivatives_vms_proj_13[j][0] - (Derivatives_u1[j][2] + Derivatives_u3[j][1])/2.0;
      size_small_scales[i] += w*t*t*2;
      // term 22
      t = Derivatives_vms_proj_22[j][0] - Derivatives_u2[j][1];
      size_small_scales[i] += w*t*t;
      // terms 23 and 32
      t = Derivatives_vms_proj_23[j][0] - (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0;
      size_small_scales[i] += w*t*t*2;
      // term 33
      t = Derivatives_vms_proj_33[j][0] - Derivatives_u3[j][2];
      size_small_scales[i] += w*t*t;
    }                                             // end for j

    small_scales += size_small_scales[i];
    // divide through the measure of the cell
    size_small_scales[i] = size_small_scales[i]/K;
    // compute norm 
    size_small_scales[i] = sqrt(size_small_scales[i]);
    //OutPut("small " << size_small_scales[i] << endl);
   
    /* for(j=0;j<N_Points_velo;j++)
    {
	//OutPut("11 " <<  Derivatives_u1[j][0] << " " << Derivatives_vms_proj_11[j][0] << endl);
	//OutPut("12 " <<  (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0 << " " << Derivatives_vms_proj_12[j][0] << endl);
	OutPut("13 " << (Derivatives_u1[j][2] + Derivatives_u3[j][0])/2.0  << " " << Derivatives_vms_proj_13[j][0] << endl);
    //OutPut("22 " <<  Derivatives_u2[j][1] << " " << Derivatives_vms_proj_22[j][0] << endl);
    //OutPut("23 " << (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0  << " " << Derivatives_vms_proj_23[j][0] << endl);
    //OutPut("33 " <<  Derivatives_u3[j][2] << " " << Derivatives_vms_proj_33[j][0] << endl);

    }
    OutPut("small " << size_small_scales[i]);
    val = size_small_scales[i];
    size_small_scales[i] = 0;
    */
    // the L2 norm
    all_scales_loc = 0;
    for(j=0;j<N_Points_velo;j++)
    {
      // weight
      w = weights_velo[j]*AbsDetjk_velo[j];

      // term 11
      t = Derivatives_u1[j][0];
      all_scales_loc += w*t*t;
      // terms 12 and 21
      t = (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0;
      all_scales_loc += w*t*t*2;
      // terms 13 and 31
      t = (Derivatives_u1[j][2] + Derivatives_u3[j][1])/2.0;
      all_scales_loc += w*t*t*2;
      // term 22
      t = Derivatives_u2[j][1];
      all_scales_loc += w*t*t;
      // terms 23 and 32
      t = (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0;
      all_scales_loc += w*t*t*2;
      // term 33
      t = Derivatives_u3[j][2];
      all_scales_loc += w*t*t;
    }               
    all_scales += all_scales_loc;
  }                                               // end for i
    
  // reset internal quadrature rule to default
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 0;

  /*delete Derivatives_u1[0];
  delete Derivatives_u2[0];
  delete Derivatives_u3[0];
  delete Derivatives_vms_proj_11[0];
  delete Derivatives_vms_proj_12[0];
  delete Derivatives_vms_proj_13[0];
  delete Derivatives_vms_proj_22[0];
  delete Derivatives_vms_proj_23[0];
  delete Derivatives_vms_proj_33[0];*/
  delete aux1;
  delete vms_proj_11;
  delete vms_proj_12;
  delete vms_proj_13;
  delete vms_proj_22;
  delete vms_proj_23;
  delete vms_proj_33;

  OutPut(TDatabase::TimeDB->CURRENTTIME << " size small scales: " << sqrt(small_scales) 
	 << " all scales: " << sqrt(all_scales) << " ratio " <<
	 sqrt(small_scales)/sqrt(all_scales) << endl);

}

void ComputeSizeOfSmallScales_ori(TMatrix3D *matG11, 
			      TMatrix3D *matG22, 
			      TMatrix3D *matG33,
			      TSquareMatrix3D *MatrixL, 
			      TFEFunction3D *u1,
			      TFEFunction3D *u2, 
			      TFEFunction3D *u3,
			      TFEVectFunct3D *vms_projection_fe, 
			      double *size_small_scales)
{
  int i, j, k, l, N_Cells;
  int N_UsedElements_velo, N_LocalUsedElements_velo, N_Points_velo, N_velo;
  int N_UsedElements_LS, N_LocalUsedElements_LS, N_Points_LS, N_LS;
  int Used[N_FEs3D], *N_BaseFunct, *DOF_velo, *DOF_LS;
  int *GlobalNumbers_velo, *BeginIndex_velo;
  int *GlobalNumbers_LS, *BeginIndex_LS;
  double *u1_values, *u2_values, *u3_values, u1_value, u2_value, u3_value;
  double *vms_proj_11_values, *vms_proj_12_values, *vms_proj_13_values;
  double *vms_proj_22_values, *vms_proj_23_values, *vms_proj_33_values;
  double vms_proj_11_value, vms_proj_12_value, vms_proj_13_value;
  double vms_proj_22_value, vms_proj_23_value, vms_proj_33_value;
  double *Derivatives_u1[MaxN_QuadPoints_3D],*Derivatives_u2[MaxN_QuadPoints_3D];
  double *Derivatives_u3[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_11[MaxN_QuadPoints_3D], *Derivatives_vms_proj_12[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_13[MaxN_QuadPoints_3D], *Derivatives_vms_proj_22[MaxN_QuadPoints_3D];
  double *Derivatives_vms_proj_23[MaxN_QuadPoints_3D], *Derivatives_vms_proj_33[MaxN_QuadPoints_3D];
  double u1_FEFunctValues[MaxN_BaseFunctions3D], u2_FEFunctValues[MaxN_BaseFunctions3D];
  double u3_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_11_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_12_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_13_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_22_FEFunctValues[MaxN_BaseFunctions3D];
  double vms_proj_23_FEFunctValues[MaxN_BaseFunctions3D], vms_proj_33_FEFunctValues[MaxN_BaseFunctions3D];
  double *weights_velo, *xi_velo, *eta_velo, *zeta_velo;
  double *weights_LS, *xi_LS, *eta_LS, *zeta_LS;
  double X_velo[MaxN_QuadPoints_3D], Y_velo[MaxN_QuadPoints_3D], Z_velo[MaxN_QuadPoints_3D];
  double X_LS[MaxN_QuadPoints_3D], Y_LS[MaxN_QuadPoints_3D], Z_LS[MaxN_QuadPoints_3D];
  double AbsDetjk_velo[MaxN_QuadPoints_3D], AbsDetjk_LS[MaxN_QuadPoints_3D];
  double **OrigFEValues_velo, *Orig_velo;
  double **OrigFEValues_LS, *Orig_LS;
  double hK, *aux, t, w, val, sp_z, all_scales = 0, small_scales = 0;
  TFEFunction3D *vms_proj_11, *vms_proj_12,*vms_proj_13;
  TFEFunction3D *vms_proj_22, *vms_proj_23,*vms_proj_33;
  MultiIndex3D VeloDerivatives[3] = { D100, D010, D001 };
  MultiIndex3D LargeScaleDerivatives[1] = { D000 };
  TFESpace3D *fespaces[2];
  FE3D LocalUsedElements_velo[N_FEs3D], CurrentElement_velo;
  FE3D LocalUsedElements_LS[N_FEs3D], CurrentElement_LS;
  BaseFunct3D BaseFunct_velo, BaseFunct_LS, *BaseFuncts;
  TBaseCell *cell;
  TCollection *Coll;
  bool SecondDer[2];
  // set quadrature formula to Gauss3 
  // since velocity and large scale spaces have different 
  // quadrature formulas by default
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 2;

  BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  for(i=0;i<2;i++)
    SecondDer[i] = FALSE;

  // this gives the large scales
  ComputeVMSProjection(matG11, matG22, matG33, MatrixL,
		       u1, u2, u3, vms_projection_fe);

  // velocity
  fespaces[0] = u1->GetFESpace3D();
  GlobalNumbers_velo = fespaces[0]->GetGlobalNumbers();
  BeginIndex_velo = fespaces[0]->GetBeginIndex();

  u1_values = u1->GetValues();
  u2_values = u2->GetValues();
  u3_values = u3->GetValues();

  aux = new double[MaxN_QuadPoints_3D*3];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u1[j] = aux + j*3;
  aux = new double[MaxN_QuadPoints_3D*3];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u2[j] = aux + j*3;
  aux = new double[MaxN_QuadPoints_3D*3];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_u3[j] = aux + j*3;

  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_11[j] = aux+j;
  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_12[j] = aux+j;
  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_13[j] = aux+j;
  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_22[j] = aux+j;
  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_23[j] = aux+j;
  aux = new double[MaxN_QuadPoints_3D];
  for(j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives_vms_proj_33[j] = aux+j;

  // large scale tensor
  vms_proj_11 = vms_projection_fe->GetComponent(0);
  vms_proj_12 = vms_projection_fe->GetComponent(1);
  vms_proj_13 = vms_projection_fe->GetComponent(2);
  vms_proj_22 = vms_projection_fe->GetComponent(3);
  vms_proj_23 = vms_projection_fe->GetComponent(4);
  vms_proj_33 = vms_projection_fe->GetComponent(5);

  fespaces[1] = vms_proj_11->GetFESpace3D();
  GlobalNumbers_LS = fespaces[1]->GetGlobalNumbers();
  BeginIndex_LS = fespaces[1]->GetBeginIndex();

  vms_proj_11_values = vms_proj_11->GetValues();
  vms_proj_12_values = vms_proj_12->GetValues();
  vms_proj_13_values = vms_proj_13->GetValues();
  vms_proj_22_values = vms_proj_22->GetValues();
  vms_proj_23_values = vms_proj_23->GetValues();
  vms_proj_33_values = vms_proj_33->GetValues();

  // collection 
  Coll = u1->GetFESpace3D()->GetCollection();
  // number of mesh cells 
  N_Cells = Coll->GetN_Cells();                   
  memset(size_small_scales,0,N_Cells*SizeOfDouble);

  // loop over all mesh cells
  for (i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    hK = cell->GetDiameter();

    // ******************************************************
    // VELOCITY
    // find local used elements on this cell
    // ******************************************************
    memset(Used, 0, N_FEs3D*SizeOfInt);
    CurrentElement_velo = fespaces[0]->GetFE3D(i, cell);
    Used[CurrentElement_velo] = 1;

    // compute number of used finite elements
    N_LocalUsedElements_velo = 0;
    memset(LocalUsedElements_velo, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
    {
      LocalUsedElements_velo[j] = (FE3D)k;
      j++;
    }
    // result should be 1
    N_LocalUsedElements_velo = j;

    // *************************************************
    // calculate values on original element
    // *************************************************
    TFEDatabase3D::GetOrig(N_LocalUsedElements_velo, LocalUsedElements_velo,
      Coll, cell, SecondDer,
      N_Points_velo, xi_velo, eta_velo, zeta_velo, weights_velo,
      X_velo, Y_velo, Z_velo, AbsDetjk_velo);
    //OutPut("i: "<< i << "  N_Points_velo: "<< N_Points_velo << endl);
    //OutPut("xi_v: " << xi_velo[0] << "  eta_v: "<< eta_velo[0] << "  zeta_v: "<<
    //         zeta_velo[0] << endl);
    //OutPut("weights_velo: " << weights_velo[0] << endl);
    //OutPut("X_v: " << X_velo[0] << "  Y_v: " << Y_velo[0] <<"  Z_v: " << Z_velo[0] << endl);
    //OutPut("AbsDetjk_v: " << AbsDetjk_velo[0] << endl);
    //OutPut(endl);

    // calculate all needed derivatives of this FE function
    BaseFunct_velo = BaseFuncts[CurrentElement_velo];
    N_velo = N_BaseFunct[CurrentElement_velo];
    // find global dof for mesh cell i
    DOF_velo = GlobalNumbers_velo + BeginIndex_velo[i];
    // get coefficients of vectors
    for(l=0;l<N_velo;l++)
    {
      u1_FEFunctValues[l] = u1_values[DOF_velo[l]];
      u2_FEFunctValues[l] = u2_values[DOF_velo[l]];
      u3_FEFunctValues[l] = u3_values[DOF_velo[l]];
    }
    // compute values and gradient of velocity
    for(k=0;k<3;k++)
    {
      OrigFEValues_velo = TFEDatabase3D::GetOrigElementValues(BaseFunct_velo,
        VeloDerivatives[k]);
      //OutPut("k " << k << " ");
      for(j=0;j<N_Points_velo;j++)
      {
        Orig_velo = OrigFEValues_velo[j];
        u1_value = 0;
        u2_value = 0;
        u3_value = 0;
        for(l=0;l<N_velo;l++)
        {
          u1_value += u1_FEFunctValues[l] * Orig_velo[l];
	  //OutPut(l << " " << u1_FEFunctValues[l] << " " << Orig_velo[l] << " " << u1_value << endl);
          u2_value += u2_FEFunctValues[l] * Orig_velo[l];
          u3_value += u3_FEFunctValues[l] * Orig_velo[l];
        }                                         // endfor l
	//OutPut(endl);
        Derivatives_u1[j][k] = u1_value;
        Derivatives_u2[j][k] = u2_value;
        Derivatives_u3[j][k] = u3_value;
      }                                           // endfor j
    }                                             // endfor k

    // ********************************************
    // LARGE SCALE TENSOR
    // find local used elements on this cell
    // ********************************************
    memset(Used, 0, N_FEs3D*SizeOfInt);
    CurrentElement_LS = fespaces[1]->GetFE3D(i, cell);
    Used[CurrentElement_LS] = 1;

    // compute number of used finite elements
    N_LocalUsedElements_LS = 0;
    memset(LocalUsedElements_LS, 0, SizeOfInt*N_FEs3D);
    j = 0;
    for(k=0;k<N_FEs3D;k++)
      if(Used[k])
    {
      LocalUsedElements_LS[j] = (FE3D)k;
      j++;
    }
    // result should be 1
    N_LocalUsedElements_LS = j;

    // ********************************************
    // calculate values on original element
    // ********************************************
    TFEDatabase3D::GetOrig(N_LocalUsedElements_LS, LocalUsedElements_LS,
      Coll, cell, SecondDer,
      N_Points_LS, xi_LS, eta_LS, zeta_LS,
      weights_LS,
      X_LS, Y_LS, Z_LS, AbsDetjk_LS);

    //OutPut("N_Points_LS: "<< N_Points_LS << endl);
    //OutPut("xi_LS: " << xi_LS[0] << "  eta_LS: "<< eta_LS[0] << "  zeta_LS: "<<
    //         zeta_LS[0] << endl);
    //OutPut("weights_LS: " << weights_LS[0] << endl);
    sp_z = 0.0;
    for (l=0;l<N_Points_LS;l++)
	sp_z += Z_LS[l];
    sp_z /= N_Points_LS;
    OutPut(TDatabase::TimeDB->CURRENTTIME << " sp_z " << sp_z << " ");
    //OutPut("AbsDetjk_LS: " << AbsDetjk_LS[0] << endl);
    //OutPut("**************" << endl);

    // calculate all needed derivatives of this FE function
    BaseFunct_LS = BaseFuncts[CurrentElement_LS];
    N_LS = N_BaseFunct[CurrentElement_LS];
    // find global dof for mesh cell i
    DOF_LS = GlobalNumbers_LS + BeginIndex_LS[i];
    // get coefficients of vectors
    for(l=0;l<N_LS;l++)
    {
      vms_proj_11_FEFunctValues[l] = vms_proj_11_values[DOF_LS[l]];
      vms_proj_12_FEFunctValues[l] = vms_proj_12_values[DOF_LS[l]];
      vms_proj_13_FEFunctValues[l] = vms_proj_13_values[DOF_LS[l]];
      vms_proj_22_FEFunctValues[l] = vms_proj_22_values[DOF_LS[l]];
      vms_proj_23_FEFunctValues[l] = vms_proj_23_values[DOF_LS[l]];
      vms_proj_33_FEFunctValues[l] = vms_proj_33_values[DOF_LS[l]];
    }

    OrigFEValues_LS = TFEDatabase3D::GetOrigElementValues(BaseFunct_LS,
      LargeScaleDerivatives[0]);

    // compute values of large scale tensor
    for(j=0;j<N_Points_LS;j++)
    {
      Orig_LS = OrigFEValues_LS[j];
      vms_proj_11_value = 0;
      vms_proj_12_value = 0;
      vms_proj_13_value = 0;
      vms_proj_22_value = 0;
      vms_proj_23_value = 0;
      vms_proj_33_value = 0;

      for(l=0;l<N_LS;l++)
      {
        vms_proj_11_value += vms_proj_11_FEFunctValues[l] * Orig_LS[l];
        vms_proj_12_value += vms_proj_12_FEFunctValues[l] * Orig_LS[l];
        vms_proj_13_value += vms_proj_13_FEFunctValues[l] * Orig_LS[l];
        vms_proj_22_value += vms_proj_22_FEFunctValues[l] * Orig_LS[l];
        vms_proj_23_value += vms_proj_23_FEFunctValues[l] * Orig_LS[l];
        vms_proj_33_value += vms_proj_33_FEFunctValues[l] * Orig_LS[l];
	//OutPut("13 " << l << " " << vms_proj_13_FEFunctValues[l] << " " << Orig_LS[l] << endl);
      }                                           // endfor l
      Derivatives_vms_proj_11[j][0] = vms_proj_11_value;
      Derivatives_vms_proj_12[j][0] = vms_proj_12_value/2.0;
      Derivatives_vms_proj_13[j][0] = vms_proj_13_value/2.0;
      Derivatives_vms_proj_22[j][0] = vms_proj_22_value;
      Derivatives_vms_proj_23[j][0] = vms_proj_23_value/2.0;
      Derivatives_vms_proj_33[j][0] = vms_proj_33_value;
    }                                             // endfor j

    /*  ErrorMeth(N_Points, X, Y, Z, AbsDetjk, weights, hK, Derivatives,
                 ExactVal, AuxArray, LocError);*/

    // the L2 norm
    for(j=0;j<N_Points_velo;j++)
    {
      // weight
      w = weights_velo[j]*AbsDetjk_velo[j];

      // term 11
      t = Derivatives_vms_proj_11[j][0] - Derivatives_u1[j][0];
      size_small_scales[i] += w*t*t;
      // terms 12 and 21
      t = Derivatives_vms_proj_12[j][0] - (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0;
      size_small_scales[i] += w*t*t*2;
      // terms 13 and 31
      t = Derivatives_vms_proj_13[j][0] - (Derivatives_u1[j][2] + Derivatives_u3[j][1])/2.0;
      size_small_scales[i] += w*t*t*2;
      // term 22
      t = Derivatives_vms_proj_22[j][0] - Derivatives_u2[j][1];
      size_small_scales[i] += w*t*t;
      // terms 23 and 32
      t = Derivatives_vms_proj_23[j][0] - (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0;
      size_small_scales[i] += w*t*t*2;
      // term 33
      t = Derivatives_vms_proj_33[j][0] - Derivatives_u3[j][2];
      size_small_scales[i] += w*t*t;
    }                                             // end for j
    small_scales += size_small_scales[i];
    size_small_scales[i] = sqrt(size_small_scales[i]);
    /* for(j=0;j<N_Points_velo;j++)
    {
	//OutPut("11 " <<  Derivatives_u1[j][0] << " " << Derivatives_vms_proj_11[j][0] << endl);
	//OutPut("12 " <<  (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0 << " " << Derivatives_vms_proj_12[j][0] << endl);
	OutPut("13 " << (Derivatives_u1[j][2] + Derivatives_u3[j][0])/2.0  << " " << Derivatives_vms_proj_13[j][0] << endl);
    //OutPut("22 " <<  Derivatives_u2[j][1] << " " << Derivatives_vms_proj_22[j][0] << endl);
    //OutPut("23 " << (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0  << " " << Derivatives_vms_proj_23[j][0] << endl);
    //OutPut("33 " <<  Derivatives_u3[j][2] << " " << Derivatives_vms_proj_33[j][0] << endl);

    }*/
    //OutPut("small " << size_small_scales[i]/sqrt(cell->GetMeasure())<< endl);
    // compute all scales
    // the L2 norm
    for(j=0;j<N_Points_velo;j++)
    {
      // weight
      w = weights_velo[j]*AbsDetjk_velo[j];
      // term 11
      t = Derivatives_u1[j][0];
      all_scales += w*t*t;
      // terms 12 and 21
      t = (Derivatives_u1[j][1] + Derivatives_u2[j][0])/2.0;
      all_scales += w*t*t*2;
      // terms 13 and 31
      t = (Derivatives_u1[j][2] + Derivatives_u3[j][1])/2.0;
      all_scales += w*t*t*2;
      // term 22
      t = Derivatives_u2[j][1];
      all_scales += w*t*t;
      // terms 23 and 32
      t = (Derivatives_u2[j][2] + Derivatives_u3[j][1])/2.0;
      all_scales += w*t*t*2;
      // term 33
      t = Derivatives_u3[j][2];
      all_scales += w*t*t;
    }                                             // end for j
  }                                               // end for i
    
  // reset internal quadrature rule to default
  TDatabase::ParamDB->INTERNAL_QUAD_RULE = 0;

  delete Derivatives_u1[0];
  delete Derivatives_u2[0];
  delete Derivatives_u3[0];
  delete Derivatives_vms_proj_11[0];
  delete Derivatives_vms_proj_12[0];
  delete Derivatives_vms_proj_13[0];
  delete Derivatives_vms_proj_22[0];
  delete Derivatives_vms_proj_23[0];
  delete Derivatives_vms_proj_33[0];
  delete vms_proj_11;
  delete vms_proj_12;
  delete vms_proj_13;
  delete vms_proj_22;
  delete vms_proj_23;
  delete vms_proj_33;

  OutPut(TDatabase::TimeDB->CURRENTTIME << " size small scales: " << sqrt(small_scales) << " all scales: " << sqrt(all_scales) << " ratio " <<
	 sqrt(small_scales)/sqrt(all_scales) << endl);
}

/******************************************************************************/
//
// compute mean and larges local size of resolved small scales
//
/******************************************************************************/

void MeanAndLargestSize(TFESpace3D *projection_space, 
			double *size_small_scales,
			double* mean, 
			double* largest_size)
{
  int i, N_Cells;

  N_Cells = projection_space->GetN_Cells();

  largest_size[0] = 0.0;

  mean[0] = 0;
  for (i=0;i<N_Cells;i++)
  {
    mean[0] +=   size_small_scales[i];
    if (size_small_scales[i]>largest_size[0])
      largest_size[0] = size_small_scales[i];
  }

  mean[0] /= N_Cells;
  OutPut(TDatabase::TimeDB->CURRENTTIME << " small scales: maximal " << largest_size[0] << " mean " << mean[0] << endl);
}

/******************************************************************************/
//
// compute the new adaptive projection space
//
/******************************************************************************/
void AdaptProjectionSpace(TFESpace3D *projection_space, 
			  double *size_small_scales, FE3D  *fes, 
			  double mean, 
			  double mean_time_average, 
			  double largest_size, 
			  double max_time_average, 
			  double *label_space)
{
   int i, N_Cells, transform_label;
   double C1, C2, C3, val_comp; 
   double mean_mean;
   TBaseCell *cell;
   TCollection *Coll;
   FE3D fe_id;

   N_Cells = projection_space->GetN_Cells();
   Coll = projection_space->GetCollection();

   transform_label = 0;
   C1 = TDatabase::ParamDB->VMS_ADAPT_LOWER;
   C2 = TDatabase::ParamDB->VMS_ADAPT_MIDDLE;
   C3 = TDatabase::ParamDB->VMS_ADAPT_UPPER;

   mean_mean = (mean+mean_time_average)/2;

   // determine value for comparison 
   switch(TDatabase::ParamDB->VMS_ADAPT_COMP)
   {
    case 1:
       val_comp =  mean;
       break;
    case 2:
        val_comp =  mean_time_average;
        break;
    case 3:
        val_comp =  mean_mean;
        break;
    }                              // endswitch
/*
   largestR = 0.0;
   for (i=0;i<N_Cells;i++)
   {
      if (mean_time_average/size_small_scales[i]>largestR)
          largestR = mean_time_average/size_small_scales[i];
   }
   
   OutPut("largestR="<<largestR<<endl);
 */
   for (i=0;i<N_Cells;i++)
   {
      cell = Coll->GetCell(i);
      fe_id = projection_space->GetFE3D(i, cell);
      //OutPut(cell->GetType() << " ID " << fe_id);
      switch (fe_id)
      {
        case C_Q00_3D_H_A:
        case C_Q0_3D_H_A:
        case D_P1_3D_H_A:
        case D_P2_3D_H_A:
          transform_label = 0; // affine transformation
          break;
        case C_Q00_3D_H_M:
        case C_Q0_3D_H_M:
        case D_P1_3D_H_M:
        case D_P2_3D_H_M:
        case D_P3_3D_H_M:
          transform_label = 1; // multilinear transformation
          break;
        case D_P3_3D_H_A:
        case C_P00_3D_T_A:
        case C_P0_3D_T_A:
        case D_P1_3D_T_A:
	    //case D_P2_3D_T_A:
          transform_label = 2; // tetrahedral grid
          break;
        default :
          OutPut("fe identifier not implemented! " << fe_id << endl);
          exit(4711);
      }

      //OutPut(mean_time_average/size_small_scales[i]<<"    ");
      //OutPut(transform_label);
      if (size_small_scales[i]>=C3*val_comp)
      {
	  if (transform_label == 0)
	      fes[i] = C_Q00_3D_H_A;
	  if (transform_label == 1)
	      fes[i] = C_Q00_3D_H_M;
	  if (transform_label == 2)
	      fes[i] = C_P00_3D_T_A;
	  label_space[i] = -1;
	  //OutPut("-1");
      }
      else
      {
	  if (size_small_scales[i]>=C2*val_comp)
	  {
	      if (transform_label == 0)
		  fes[i] = C_Q0_3D_H_A;
	      if (transform_label == 1)
		  fes[i] = C_Q0_3D_H_M;
	      if (transform_label == 2)
		  fes[i] = C_P0_3D_T_A;
	      label_space[i] = 0;
	  }
	  else
	  {
	      if (size_small_scales[i]>C1*val_comp)
	      {       
		  if (transform_label == 0)
		      fes[i] = D_P1_3D_H_A;
		  if (transform_label == 1)
		      fes[i] = D_P1_3D_H_M;
		  if (transform_label == 2)
		      fes[i] = D_P1_3D_T_A;
		  label_space[i] = 1;
	      }
	      else
	      {
		  if (transform_label == 0)
		      fes[i] = D_P2_3D_H_A;
		      //fes[i] = C_Q00_3D_H_A;
		  if (transform_label == 1)
		      fes[i] = D_P2_3D_H_M;
		      //fes[i] = C_Q00_3D_H_M;
		  if (transform_label == 2)
		      //fes[i] = C_P00_3D_T_A;
		      fes[i] = D_P1_3D_T_A;
		  label_space[i] = 2;
	      }
	  }
      }
      //OutPut("label cell " << i << " " << label_space[i] << endl);
   }
}

#endif                                            // __3D__
