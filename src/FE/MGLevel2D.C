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
// @(#)MGLevel2D.C        1.8 05/05/00
//
// Class:       TMGLevel2D
// Purpose:     store all data for one level in a multi grid method
//
// Author:      Gunar Matthies 02.11.1998
//
// History:     02.11.1998 start of implementation
//
// =======================================================================

#include <MGLevel2D.h>
#include <FESpace2D.h>
#include <Database.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <FEFunction2D.h>
// #include <ParFECommunicator2D.h>

#include <stdlib.h>
#include <string.h>

/** constructor */
TMGLevel2D::TMGLevel2D(int level, TSquareMatrix2D *a,
                       double *rhs, double *sol, int n_aux,
                       int *permutation)
{
  int i;
  double *aux;

  Level = level;

  FESpace = a->GetFESpace();

  N_Active = FESpace->GetN_ActiveDegrees();
  HangingNodeBound = FESpace->GetHangingBound();
  N_Dirichlet = FESpace->GetN_Dirichlet();
  N_DOF = FESpace->GetN_DegreesOfFreedom();

/*
  cout << "N_Active: " << N_Active << endl;
  cout << "HangingNodeBound: " << HangingNodeBound << endl;
  cout << "N_Dirichlet: " << N_Dirichlet << endl;
  cout << "N_DOF: " << N_DOF << endl;
*/

  A = a;
  MatrixStructure = a->GetMatrixStructure();
  RowPtr = a->GetRowPtr();
  KCol = a->GetKCol();
  Entries = a->GetEntries();

  Rhs = rhs;
  X = sol;

  N_Aux = n_aux;
  Aux = new double* [N_Aux]; 
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Additional = NULL;

  Permutation = permutation;
}

#ifdef _MPI
/** constructor for parallel */
TMGLevel2D::TMGLevel2D(int level, TSquareMatrix2D *a, double *rhs, double *sol, 
                       TFEFunction2D *c, TParFECommunicator2D *parComm,TFESpace2D *ownScalarSpace, int n_aux,
                       int *permutation)
{
  int i;
  double *aux;
  char UString[] = "u";

  Level = level;

  FESpace = a->GetFESpace();

  N_Active = FESpace->GetN_ActiveDegrees();
  HangingNodeBound = FESpace->GetHangingBound();
  N_Dirichlet = FESpace->GetN_Dirichlet();
  N_DOF = FESpace->GetN_DegreesOfFreedom();

/*
  cout << "N_Active: " << N_Active << endl;
  cout << "HangingNodeBound: " << HangingNodeBound << endl;
  cout << "N_Dirichlet: " << N_Dirichlet << endl;
  cout << "N_DOF: " << N_DOF << endl;
*/

  A = a;
  C = c;
  MatrixStructure = a->GetMatrixStructure();
  RowPtr = a->GetRowPtr();
  KCol = a->GetKCol();
  Entries = a->GetEntries();

  Rhs = rhs;
  X = sol;

  N_Aux = n_aux;
  Aux = new double* [N_Aux]; 
  aux = new double[N_Aux*N_DOF];
  for(i=0;i<N_Aux;i++)
    Aux[i] = aux+i*N_DOF;

  Additional = NULL;

  Permutation = permutation;
  
  
  ParComm = parComm;
  OwnScalarSpace = ownScalarSpace;  
  OwnN_DOF = OwnScalarSpace->GetN_DegreesOfFreedom();
  OwnSolArray = new double[OwnN_DOF];
 
  OwnC = new TFEFunction2D(OwnScalarSpace, UString, UString, OwnSolArray, OwnN_DOF);
}
#endif

/** destructor */
TMGLevel2D::~TMGLevel2D()
{
  delete [] Aux[0];
  delete Aux;

  if(Additional)
    delete Additional;
} // ~TMGLevel2D

/** return i-th auxiliary vector */
double *TMGLevel2D::GetAuxVector(int i)
{
  double *ret;

  if(i<N_Aux)
    ret = Aux[i];
  else
  {
    cerr << "Not enough aux vectors in __FILE__!" << endl;
    exit(-1);
  }
 
  return ret;
} // GetAuxVector

// calculate defect d=f-Ax
void TMGLevel2D::Defect(double *sol, double *f, double *d, double &res)
{
  
  
  ScalarDefect(A, sol, f, d, res);
  
#ifdef _MPI  
  int i, rank, *MasterOfDof, dof;
  double res_global;
  
  MPI_Comm_rank(ParComm->GetComm(), &rank); 
//   MasterOfDof = ParComm->GetMaster();
  
//   for(i=0; i<N_DOF; i++)
//     if(MasterOfDof[i] == rank)
//       res += d[i]*d[i];
//     
//   MPI_Allreduce(&res, &res_global, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
//   res = sqrt(res_global); 
#endif  
    
} // end Defect

// SOR smoother
void TMGLevel2D::SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,k,l,index;
  double s, t, diag;
  double omega;

  omega = Parameters[0];

  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  for(ii=0;ii<N_Active;ii++)
  {
    i = ii;
    // i = Permutation[ii];
    // cout << "row: " << i << "   " << endl;
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        diag = Entries[j];
      }
      else
      {
        s -= Entries[j] * sol[index];
      }
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
    // cout << "sol[i]: " << sol[i] << endl;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

} // SOR

// SSOR smoother
void TMGLevel2D::SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double s, t, diag;
  double omega;

  omega = Parameters[0];
  //OutPut("omega " << omega << endl);
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * sol[index];
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

  // set active nodes
  j = RowPtr[N_Active]-1;
  for(i=N_Active-1;i>=0;i--)
  {
    s = f[i];
    k = RowPtr[i];
    for(;j>=k;j--)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * sol[index];
    } // endfor j
    t = sol[i];
    sol[i] = omega*(s/diag-t) + t;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[HangingNodeBound]-1;
  for(i=HangingNodeBound-1;i>=N_Active;i--)
  {
    s = f[i];
    k = RowPtr[i];
    for(;j>=k;j--)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i
} // SSOR

// Jacobi smoother
void TMGLevel2D::Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double t, s, diag, omega;

  omega = Parameters[0];

  memcpy(aux, sol, N_DOF*SizeOfDouble);

  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index == i)
        diag = Entries[j];
      else
        s -= Entries[j] * aux[index];
    } // endfor j
    t = aux[i];
    sol[i] = (1-omega)*t + omega*s/diag;
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i

} // Jacobi

void TMGLevel2D::Update(double *sol, double *upd)
{
  int i;

  for(i=0;i<N_DOF;i++)
  {
    sol[i] += TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR*upd[i];
  }
}

void TMGLevel2D::Reset(double *vect)
{
  memset(vect, 0, N_DOF*SizeOfDouble);
}

void TMGLevel2D::CorrectNodes(double *vect)
{
  int i, j, k, index;
  double s;

  memset(vect+HangingNodeBound, 0, N_Dirichlet*SizeOfDouble);

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s = 0;
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
      {
        s -= Entries[j] * vect[index];
      }
    } // endfor j
    vect[i] = s;
  } // endfor i
} // CorrectNodes

/** correct defect */
void TMGLevel2D::CorrectDefect(double *vect)
{
  memset(vect+N_Active, 0, SizeOfDouble*(N_DOF-N_Active));
}

// block 2x2 smoother
void TMGLevel2D::Block2x2(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index,maxindex;
  double s1, s2, t, mt, m11, m12, m21, m22, max;
  double omega;

  omega = Parameters[0];

  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
    k = RowPtr[i+1];
    max = 0;
    maxindex = -1;
    m11 = m12 = m21 = m22 = 0;
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == i)
      {
        m11 = Entries[j];
      }
      else
      {
        t = Entries[j];
        if( t > max ) // consider only negative entries
        {
          max = t;
          maxindex = index;
        }
      }
    } // endfor j

    s1 = f[i];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index == maxindex)
        m12 = Entries[j];
      else 
        if(index != i)
          s1 -= Entries[j] * sol[index];
    } // endfor j

    if(maxindex != -1)
    {
      s2 = f[maxindex];
      k=RowPtr[maxindex+1];
      for(j=RowPtr[maxindex];j<k;j++)
      {
        index = KCol[j];
        if(index == maxindex)
          m22 = Entries[j];
        else
          if(index == i)
            m21 = Entries[j];
          else
            s2 -= Entries[j] *sol[index];
      } // endfor j

      t = (m11*m22 - m21*m12);
      sol[i] = (s1*m22 - s2*m12) / t;
      sol[maxindex] = (s2*m11 - s1*m21) / t;
    }
    else
    {
      // no nondiagonal matrix entry found
      sol[i] = s1/m11;
    }
  } // endfor i

  // set hanging nodes 
  j = RowPtr[N_Active];
  for(i=N_Active;i<HangingNodeBound;i++)
  {
    s1 = f[i];
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      if(index != i)
      {
        s1 -= Entries[j] * sol[index];
      }
    } // endfor j
    sol[i] = s1;
  } // endfor i

} // Block2x2

// generate ILU decomposition
void TMGLevel2D::ILUDecomposition()
{
  int i,j,jj,k,l,N_;
  double diag, pivot, update;
  int begin, end, beginJ, endJ, found, index;
  static double beta_ilu = TDatabase::ParamDB->SC_ILU_BETA;

  N_=RowPtr[N_DOF];
  Additional = new double[N_];
  memcpy(Additional, Entries, N_*SizeOfDouble);

  // compute ILU decomposition
  // loop over all rows
  for(i=0;i<N_DOF;i++)
  {
    // pointer to the entries for the columns
    begin = RowPtr[i];
    end = RowPtr[i+1];
    //diag = Additional[begin];
    // find diagonal entry
    diag = -4711;
    for(j=begin;j<end;j++)
    {
      index = KCol[j];
      if (index==i)
      {
        diag = Additional[j];
        break;
      }
    }
    if(fabs(diag)<1e-8)
    {
      cerr << "ILU decomposition failed" << endl;
      return;
    }

    for(j=begin+1;j<end;j++)
    {
      if( (jj=KCol[j]) > i)
      {
        beginJ = RowPtr[jj];
        endJ = RowPtr[jj+1];
        found = 0;
        for(k=beginJ+1;k<endJ;k++)
        {
          if(KCol[k] != i) continue;
          pivot = Additional[k]/diag;
          Additional[k] = pivot;
          found = 1;
          break;
        } // endfor k
        if(!found) continue;
        Additional[beginJ] -= pivot*Additional[j];
        for(k=begin+1;k<end;k++)
        {
          if( KCol[k] <= KCol[j] ) continue;
          update = Additional[k] * pivot;
          for(l=beginJ+1;l<endJ;l++)
          {
            if(KCol[k] == KCol[l])
            {
              Additional[l] -= update;
              update = 0;
              break;
            } // endif
          } // endfor l
          Additional[beginJ] += beta_ilu * fabs(update);
        } // endfor k
      } //endif
    } // endfor j
  } // endfor i

/*
  cout << "Matrix A" << endl;
  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      cout << i << "   " << KCol[j] << "  " << Entries[j] << endl;
  }
  cout << "-----------------" << endl;

  cout << "ILU decomposition" << endl;
  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      cout << i << "   " << KCol[j] << "  " << Additional[j] << endl;
  }
  cout << "-----------------" << endl;
*/

} // ILUDecomposition

// ILU smoother
void TMGLevel2D::ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k, begin, end;
  double diag;

  if (Additional==NULL)
  {
    cout << "do ILU decomposition" << endl;
    ILUDecomposition();
  }

  // do ILU smoothing

  // invert lower
  for(i=1; i<N_DOF; i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      if( (k=KCol[j]) < i)
        aux[i] -= Additional[j]*aux[k];
  }

  // invert upper
  for (i=N_DOF-1; i>=0; i--)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    for(j=begin;j<end;j++)
      if( (k=KCol[j]) > i)
        aux[i] -= Additional[j]*aux[k];
    //aux[i] /= Additional[begin];
    // find diagonal entry
    for(j=begin;j<end;j++)
    {
       if (KCol[j]==i)
      {
        diag = Additional[j];
        break;
      }
    }
    aux[i] /= diag;
  }

  for(i=0;i<N_DOF;i++)
  {
    sol[i] += 1.0*aux[i];
  }
} // ILU

/** solve exact on this level */
void TMGLevel2D::SolveExact(double *u1, double *rhs1)
{  
  double *a, *b;
  int i,j,k,l,index, N_DOF2 = N_DOF*N_DOF, end, begin, m;
  double value;

  // arrays for matrix and one vector
  a = new double[N_DOF2];
  b = new double[N_DOF];
  // initialize
  memset(a, 0, N_DOF2*SizeOfDouble);
  // fill matrix columnwise
  j = RowPtr[0];
  for(i=0;i<N_DOF;i++)
  {
    k = RowPtr[i+1];
    for(;j<k;j++)
    {
      index = KCol[j];
      value = Entries[j];
      a[index*N_DOF+i] = value;
    }
  } // endfor i

  // copy into local data
  memcpy(b, rhs1, N_DOF*SizeOfDouble);
  
  // for (i=0;i<N_DOF;i++)
  // cout<< "b("<<i+1<<") = " << b[i] << ";"<<endl;

  SolveLinearSystemLapack(a, b, N_DOF, N_DOF);

  // for (i=0;i<N_DOF;i++)
  //  cout<< "c("<<i+1<<") = " << b[i] << ";"<<endl;
  //exit(1);
  // copy from local data
  memcpy(u1, b, N_DOF*SizeOfDouble);

  delete a;
  delete b;
}

/** step length control multigrid cycle */
double TMGLevel2D::StepLengthControl(double *u, 
                                     double *uold, 
                                     double *def,
                                     int N_Parameters, 
                                     double *Parameters)
{
  double *x,*y,omega,numerator,nominator;
  int i;

  // allocate auxiliary array
  x = new double[2*N_DOF];
  y = x+N_DOF;

  for (i=0;i<N_DOF;i++)
    x[i] = u[i]-uold[i];
  memset(y,0,N_DOF*SizeOfDouble);
  // if (Ddot(N_DOF-N_Active,x+N_Active,x+N_Active)>0)
  //   cout << "Update in non-active nodes found !!" << 
  //     Ddot(N_DOF-N_Active,x+N_Active,x+N_Active) << 
  //     __FILE__ << endl;
  if (Ddot(N_DOF-HangingNodeBound,x+HangingNodeBound,x+HangingNodeBound)>0)
    cout << "Update in non-active nodes found !!" << 
      Ddot(N_DOF-HangingNodeBound,x+HangingNodeBound,x+HangingNodeBound) << 
      __FILE__ << endl;
  
  // compute matrix times update
  MatVect(A,x,y);

  numerator = Ddot(N_DOF,def,y);
  nominator = Ddot(N_DOF,y,y);

  if (nominator > 0)
    omega = numerator/nominator;
  else
    {
      if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.5;
      cout << "MESSAGE : Step length control failed. Set omega = " << omega << endl;
    }
  if (fabs(omega)<0.1)
    {
       if(N_Parameters>0)
        omega = Parameters[0];
      else
        omega = 0.9;
    }     
  delete x;
  //cout << sqrt(Ddot(N_DOF,def,def)) << " " << numerator << " " << nominator<< " omega " << omega << endl;
  return(omega);
}
