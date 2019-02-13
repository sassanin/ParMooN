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
// Class:       TMGLevel3D
// Purpose:     store all data for one level in a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#include <MGLevel3D.h>
#include <FESpace3D.h>
#include <Database.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <FEFunction3D.h>

#ifdef _MPI 
#include <ParFECommunicator3D.h>
#endif

#include <stdlib.h>
#include <string.h>



#define PreSmooth -1
#define CoarseSmooth 0
#define PostSmooth 1

// #include <omp.h>

double tSor=0.0;
double tD=0.0,tS=0.0;

/** constructor */
TMGLevel3D::TMGLevel3D(int level, TSquareMatrix3D *a,
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
TMGLevel3D::TMGLevel3D(int level, TSquareMatrix3D *a, double *rhs, double *sol, 
                       TParFECommunicator3D *parComm, TParFEMapper3D *parMapper, int n_aux,
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
  
  
  if(TDatabase::ParamDB->MapperType != 2)
  {
    N_InterfaceM = parMapper->GetN_InterfaceM();
    N_Int        = parMapper->GetN_Int_light();
    N_Dept1      = parMapper->GetN_Dept1();
    N_Dept2      = parMapper->GetN_Dept2();
//     N_Dept3      = parMapper->GetN_Dept3();
    
#ifdef _HYBRID
    N_CMaster  = parMapper->GetN_CMaster();
    ptrCMaster = parMapper->GetptrCMaster();
    
    N_CInt  = parMapper->GetN_CInt();
    ptrCInt = parMapper->GetptrCInt();
    
    N_CDept1  = parMapper->GetN_CDept1();
    ptrCDept1 = parMapper->GetptrCDept1();
    
    N_CDept2  = parMapper->GetN_CDept2();
    ptrCDept2 = parMapper->GetptrCDept2();
    
//     N_CDept3  = parComm->GetN_CDept3();
//     ptrCDept3 = parComm->GetptrCDept3();
#endif
  }
  
  ParComm = parComm; 
  ParMapper = parMapper;
  Temp_arr    = new double[N_DOF];

  Reorder_M  = ParMapper->GetReorder_M();
  Reorder_D1 = ParMapper->GetReorder_D1();
  Reorder_D2 = ParMapper->GetReorder_D2();
  Reorder_I  = ParMapper->GetReorder_I();

}
#endif


/** destructor */
TMGLevel3D::~TMGLevel3D()
{
  delete Aux[0];
  delete Aux;

  if(Additional)
    delete Additional;
} // ~TMGLevel3D

/** return i-th auxiliary vector */
double *TMGLevel3D::GetAuxVector(int i)
{
  double *ret;

  if(i<N_Aux)
    ret = Aux[i];
  else
  {
    cerr << "Not enough aux vectors in MGLevel3D.C!" << endl;
    exit(-1);
  }
 
  return ret;
} // GetAuxVector

void TMGLevel3D::Defect(double *sol, double *f, double *d, double &res)
{
  
  double t1,t2;
#ifdef _MPI
  int i, rank, *MasterOfDof, dof,numThreads= TDatabase::ParamDB->OMPNUMTHREADS;
  double res_global=0,res1=0;
  
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank); 
  MasterOfDof = ParComm->GetMaster();
  
  t1 = MPI_Wtime();
  
  ParComm->CommUpdate(sol);
  
#else
  t1 = GetTime();
#endif

  ScalarDefect(A, sol, f, d, res);

#ifdef _MPI
  t2 = MPI_Wtime();
  tS += (t2-t1);
  t1 = MPI_Wtime();
#else
  t2  = GetTime();
  tS += (t2-t1);
  t1  = GetTime();
#endif
  
#ifdef _MPI  
#ifdef _HYBRID
  omp_set_num_threads(numThreads);
#pragma omp parallel default(shared) private(i)
{
#pragma omp for schedule(guided) nowait reduction(+:res1)
#endif
  for(i=0; i<N_DOF; i++){
    if(MasterOfDof[i] == rank)
      res1 += d[i]*d[i];
  }
#ifdef _HYBRID
}
#endif

  MPI_Allreduce(&res1, &res_global, 1, MPI_DOUBLE, MPI_SUM, TDatabase::ParamDB->Comm);
  res = sqrt(res_global);
  //printf("\n.......rank=%d.............res::%lf         res_global::%lf             res1::%lf",rank,res,res_global,res1);
#endif 
#ifdef _MPI
t2 = MPI_Wtime();
tD +=(t2-t1);
#else
t2 = GetTime();
tD +=(t2-t1);
#endif
} // end Defect

// SOR smoother
void TMGLevel3D::SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,k,l,index;
  double s, t, diag;
  double omega;

  omega = Parameters[0];
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  int *master =ParComm->GetMaster();
#endif
  
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i,ii,s,k,j,index,diag,t)
{
#pragma omp for schedule(dynamic) nowait 
#endif
 
  for(ii=0;ii<N_Active;ii++)
  {
#ifdef _MPI
    if(master[ii] != rank)
      continue;
#endif
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
//   j = RowPtr[N_Active];
#ifdef _HYBRID
#pragma omp for schedule(dynamic) nowait 
#endif
  for(i=N_Active;i<HangingNodeBound;i++)
  {
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
    {
      index = KCol[j];
      if(index != i)
        s -= Entries[j] * sol[index];
      else
        diag = Entries[j];
    } // endfor j
    sol[i] = s/diag;
  } // endfor i
#ifdef _HYBRID
}
#endif

} // SOR

// SSOR smoother
void TMGLevel3D::SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double s, t, diag;
  double omega;
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  int *master =ParComm->GetMaster();
#endif

  omega = Parameters[0];
 
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i+1];
    for(j = RowPtr[i];j<k;j++)
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
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i+1];
    for(j=RowPtr[i];j<k;j++)
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
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif    
    s = f[i];
    k = RowPtr[i];
    for(j=RowPtr[i+1]-1;j>=k;j--)
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
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i];
    for(j=RowPtr[i+1]-1;j>=k;j--)
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
void TMGLevel3D::Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k,l,index;
  double t, s, diag, omega;

  omega = Parameters[0];
  
#ifdef _MPI
  int rank;
  MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);
  int *master =ParComm->GetMaster();
#endif

  //SHOULD THIS BE AFTER SETTING DRICHLET????
  memcpy(aux, sol, N_DOF*SizeOfDouble);

//----------------------------------------------------------------------------------------------------  
//   double rhs_norm        = 0.0, sol_norm        = 0.0;
// #ifdef _MPI
//    double rhs_norm_global = 0.0, sol_norm_global = 0.0;
   
//    for(i=0; i<N_DOF; i++)
//     if(master[i] == rank){
//       rhs_norm += f[i]*f[i];
//       sol_norm += sol[i]*sol[i];
//     }
//    MPI_Allreduce(&rhs_norm, &rhs_norm_global, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
//    MPI_Allreduce(&sol_norm, &sol_norm_global, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
   
//    sol_norm = sqrt(sol_norm_global);
//    rhs_norm = sqrt(rhs_norm_global);
   
//    if(rank==0)
// #else
//      sol_norm = sqrt(Ddot((N_DOF)   ,sol,sol));
//      rhs_norm = sqrt(Ddot((N_DOF)   ,f,f));
// #endif
//     {
//        cout << "sol: " << sol_norm << endl;
//        cout << "rhs: " << rhs_norm << endl;
//     }
  
//------------------------------------------------------------------------------------------------------
  // set Dirichlet node
//   memcpy(sol+N_Active, f+N_Active, 
//            (N_DOF-N_Active)*SizeOfDouble);
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

  // set active nodes
  j = RowPtr[0];
  for(i=0;i<N_Active;i++)
  {
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i+1];
    for(j = RowPtr[i];j<k;j++)
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
   // printf("yes\n");
#ifdef _MPI
    if(master[i] != rank)
      continue;
#endif
    s = f[i];
    k = RowPtr[i+1];
    for(j = RowPtr[i];j<k;j++)
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

void TMGLevel3D::Update(double *sol, double *upd)
{
  int i;
#ifdef _HYBRID
#pragma omp parallel default(shared) private(i)
 {
 //  printf("num thrds is %d---\n",omp_get_num_threads());
 #pragma omp for schedule(static) nowait 
#endif
  for(i=0;i<N_DOF;i++)
  {
    sol[i] += TDatabase::ParamDB->SC_GMG_DAMP_FACTOR_SCALAR*upd[i];
  }
 #ifdef _HYBRID
 }
 #endif
}

void TMGLevel3D::Reset(double *vect)
{
  memset(vect, 0, N_DOF*SizeOfDouble);
}

void TMGLevel3D::CorrectNodes(double *vect)
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
void TMGLevel3D::CorrectDefect(double *vect)
{
  memset(vect+N_Active, 0, SizeOfDouble*(N_DOF-N_Active));
  
//   int i, *master;
//   double sum = 0.,global_sum;
// #ifdef _MPI
//   int rank;
//   MPI_Comm_rank(ParComm->GetComm(), &rank);
//   master = ParComm->GetMaster();
//   for(i=0;i<N_DOF;i++){
//     if(master[i] == rank){
//       sum+=vect[i]*vect[i];
//     }
//   }
//   MPI_Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, ParComm->GetComm());
//   sum = global_sum;
//   if(rank==0)
// #else
//   sum = ( Ddot(N_DOF,vect,vect) );
// #endif
//   
//   printf("rhs: %lf\n",sqrt(sum));
}

// block 2x2 smoother
void TMGLevel3D::Block2x2(double *sol, double *f, double *aux,
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

//shamim :: check ILU (currently iterating for slaves also)
// generate ILU decomposition
void TMGLevel3D::ILUDecomposition()
{
  int i,j,jj,k,l,N_;
  double diag, pivot, update;
  int begin, end, beginJ, endJ, found;
  static double beta_ilu = TDatabase::ParamDB->SC_ILU_BETA;

  N_=RowPtr[N_DOF];
  Additional = new double[N_];
  memcpy(Additional, Entries, N_*SizeOfDouble);

  for(i=0;i<N_DOF;i++)
  {
    begin = RowPtr[i];
    end = RowPtr[i+1];
    diag = Additional[begin];
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

//shamim :: check ILU (currently iterating for slaves also)
// ILU smoother
void TMGLevel3D::ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int i,j,k, begin, end;

  if (Additional==NULL)
  {
    if (TDatabase::ParamDB->SC_VERBOSE>1)    
      OutPut("do ILU decomposition" << endl);
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
    aux[i] /= Additional[begin];
  }

  for(i=0;i<N_DOF;i++)
  {
    sol[i] += 1.0*aux[i];
  }
} // ILU

//shamim :: check exact (currently iterating for slaves also)
/** solve exact on this level */
void TMGLevel3D::SolveExact(double *u1, double *rhs1)
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
double TMGLevel3D::StepLengthControl(double *u, 
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

#ifdef _MPI
void TMGLevel3D::SOR_Re(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
  int ii, i,j,k,l,index,rank,loop;
  double s, t, diag;
  double omega;
  int repeat = TDatabase::ParamDB->Par_P6;
  if(repeat <= 0)
    repeat = 1;
  
  omega = Parameters[0];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  for(loop=0;loop<repeat;loop++)
  {
      // set Dirichlet nodes
      memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
	      N_Dirichlet*SizeOfDouble);

    //########################################## MASTERS DOFS ########################################################//
      Reorder = ParMapper->GetReorder_M();
      for(ii=0;ii<N_InterfaceM;ii++)
      {
	i = Reorder[ii];
	if(i >= N_Active)     continue;
	
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
      
      if(loop == (repeat-1))
	ParComm->CommUpdateMS(sol);
    //################################################################################################################//

    //########################################## INDEPENDENT DOFS ####################################################//  
      
      Reorder = ParMapper->GetReorder_I();
      for(ii=0;ii<N_Int;ii++)
      {
	i = Reorder[ii];
	if(i >= N_Active)     continue;
	
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
    //################################################################################################################//

    //########################################## DEPENDENT1 DOFS #####################################################//
      Reorder = ParMapper->GetReorder_D1();
      for(ii=0;ii<N_Dept1;ii++)
      {
	i = Reorder[ii];
	if(i >= N_Active)     continue;
	
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
    //################################################################################################################//

    //########################################## DEPENDENT2 DOFS #####################################################//
      Reorder = ParMapper->GetReorder_D2();
      for(ii=0;ii<N_Dept2;ii++)
      {
	i = Reorder[ii];
	if(i >= N_Active)     continue;
	
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
    //################################################################################################################//

    // //########################################### DEPENDENT3 DOFS ####################################################//
    //   Reorder = ParComm->GetReorder_D3();
    //   for(ii=0;ii<N_Dept3;ii++)
    //   {
    //     i = Reorder[ii];
    //     if(i >= N_Active)     continue;
    //     
    //     s = f[i];
    //     k = RowPtr[i+1];
    //     for(j=RowPtr[i];j<k;j++)
    //     {
    //       index = KCol[j];
    //       if(index == i)
    //       {
    //         diag = Entries[j];
    //       }
    //       else
    //       {
    //         s -= Entries[j] * sol[index];
    //       }
    //     } // endfor j
    //     t = sol[i];
    //     sol[i] = omega*(s/diag-t) + t;
    //     // cout << "sol[i]: " << sol[i] << endl;
    //   } // endfor i
    // //################################################################################################################//

   if(loop == (repeat-1))
      ParComm->CommUpdateH1(sol);

    //############################################# Hanging NODES ####################################################//  
      // set hanging nodes
      int *master = ParComm->GetMaster();
      for(i=N_Active;i<HangingNodeBound;i++)
      {
	if(master[i] != rank)
	  continue;
	s = f[i];
	k = RowPtr[i+1];
	for(j=RowPtr[i];j<k;j++)
	{
	  index = KCol[j];
	  if(index != i)
	    s -= Entries[j] * sol[index];
	  else
	    diag = Entries[j];
	} // endfor j
	sol[i] = s/diag;
      } // endfor i
  }//loop
}
//################################################################################################################//
#endif

//*****************************************************************************************************************//
#ifdef _HYBRID

void TMGLevel3D::SOR_Re_Color(double *sol, double *f, double *aux, int N_Parameters, double *Parameters,int smooth) 
{
//   SOR_Re(sol,f,aux,N_Parameters,Parameters);
//   return;
  int *master = ParComm->GetMaster();
  int itr,ii, i,j,jj,k,l,index,rank,tid,nrows=0,numThreads,end,loop;
  double s, t, diag;
  double omega;
  
  int repeat = TDatabase::ParamDB->Par_P6;
  
  if(repeat <= 0)
    repeat = 1;
  
  omega = Parameters[0];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  omp_set_num_threads(numThreads);
  
  

  if(smooth == -1)
    end = TDatabase::ParamDB->SC_PRE_SMOOTH_SCALAR;
  else if(smooth == 0)
    end = 1;//TDatabase::ParamDB->SC_COARSE_MAXIT_SCALAR;
  else
    end = TDatabase::ParamDB->SC_POST_SMOOTH_SCALAR;
  
#pragma omp parallel default(shared) private(i,ii,s,k,j,jj,tid,index,diag,t,loop,itr)
  {
    // set Dirichlet nodes    
    for(itr=0;itr<end;itr++)
    {
      for(loop=0;loop<repeat;loop++)
      {
#pragma omp master
    memcpy(sol+HangingNodeBound, f+HangingNodeBound, N_Dirichlet*SizeOfDouble);
#pragma omp barrier
      //########################################## MASTERS DOFS ########################################################//
	    if(itr == 0)
	    {
// 	      Reorder = ParMapper->GetReorder_M();
	      for(ii=0;ii<N_CMaster;ii++)
	      {
		#pragma omp for schedule(guided) 
		for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
		{
		  i = Reorder_M[jj];
		  if(i >= N_Active)     continue;
      
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
		} // endfor jj
	      } //end for ii
	    
	      if(loop == (repeat-1))
	      {
	        #pragma omp master
	        {
		  ParComm->CommUpdateMS(sol);
	        }
	        #pragma omp barrier
	      }
	    } //end firstTime
      //########################################## DEPENDENT1 DOFS #####################################################//
// 	    Reorder = ParMapper->GetReorder_D1();
	    for(ii=0;ii<N_CDept1;ii++)
	    {
	      #pragma omp for schedule(guided) 
	      for(jj=ptrCDept1[ii];jj<ptrCDept1[ii+1];jj++)
	      {
		i = Reorder_D1[jj];
		if(i >= N_Active)     continue;
		
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
	      } // endfor jj
	    } //end for ii  
      //################################################################################################################//
	    #pragma omp master
	    {
	      if(loop == (repeat-1))
		ParComm->CommUpdateH1(sol);
	    }    
      //########################################## DEPENDENT2 DOFS #####################################################//
// 	    Reorder = ParMapper->GetReorder_D2();
	    for(ii=0;ii<N_CDept2;ii++)
	    {
	    #pragma omp for schedule(guided) 
	    for(jj=ptrCDept2[ii];jj<ptrCDept2[ii+1];jj++)
	    {
	      i = Reorder_D2[jj];
	      if(i >= N_Active)     continue;
	  
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
	    } // endfor jj
	  } //end for ii
      //################################################################################################################//    

      //########################################## MASTERS DOFS ########################################################//
	    if(itr!=(end-1))
	    {
// 	    Reorder = ParMapper->GetReorder_M();
	    for(ii=0;ii<N_CMaster;ii++)
	    {
	      #pragma omp for schedule(guided) 
	      for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
	      {
		i = Reorder_M[jj];
		if(i >= N_Active)     continue;
	  
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
	      } // endfor jj
	    } //end for ii
	  } //end !lastTime
      //################################################################################################################//
	    if(itr!=(end-1))
	    #pragma omp master
	    {
	      if(loop == (repeat-1))
		ParComm->CommUpdateMS(sol);
	    }    
      //########################################## INDEPENDENT DOFS ####################################################//  
// 	    Reorder = ParMapper->GetReorder_I();
	    for(ii=0;ii<N_CInt;ii++)
	    {
	    #pragma omp for schedule(guided) 
	    for(jj=ptrCInt[ii];jj<ptrCInt[ii+1];jj++)
	    {
	      i = Reorder_I[jj];
	      if(i >= N_Active)     continue;
	    
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
	      } // endfor jj
	      t = sol[i];
	      sol[i] = omega*(s/diag-t) + t;
	      // cout << "sol[i]: " << sol[i] << endl;
	    } // endfor jj
	  } //endfor ii
      //################################################################################################################//      

      //############################################# Hanging NODES ####################################################//  
	    // set hanging nodes
// 	    int *master = ParComm->GetMaster();

	    #pragma omp for schedule(dynamic) nowait 

	    for(i=N_Active;i<HangingNodeBound;i++)
	    {

	    if(master[i] != rank)
	      continue;

	    s = f[i];
	    k = RowPtr[i+1];
	    for(j=RowPtr[i];j<k;j++)
	    {
	      index = KCol[j];
	      if(index != i)
		s -= Entries[j] * sol[index];
	      else
		diag = Entries[j];
	    } // endfor j
	    sol[i] = s/diag;
	  } // endfor i

      }//loop
    }
//################################################################################################################//
  
    
  }
}  

void TMGLevel3D::SOR_Re_Color(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters, bool firstTime, bool lastTime)
{
//     SOR_Re(sol,f,aux,N_Parameters,Parameters);
//   return;
  
  int ii, i,j,jj,k,l,index,rank,tid,nrows=0,numThreads;
  double s, t, diag;
  double omega;
  int *master = ParComm->GetMaster(); 
  omega = Parameters[0];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  omp_set_num_threads(numThreads);
//   omp_set_num_threads(1);
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

#pragma omp parallel default(shared) private(i,ii,s,k,j,jj,tid,index,diag,t)
  {
//########################################## MASTERS DOFS ########################################################//
    if(firstTime)
    {
//       Reorder = ParMapper->GetReorder_M();
      for(ii=0;ii<N_CMaster;ii++)
      {
	#pragma omp for schedule(guided) 
	for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
	{
	  i = Reorder_M[jj];
          if(i >= N_Active)     continue;
    
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
        } // endfor jj
      } //end for ii
      
      #pragma omp master
      {
	ParComm->CommUpdateMS(sol);
      }
      #pragma omp barrier
    } //end firstTime
//################################################################################################################//
    
//########################################## DEPENDENT1 DOFS #####################################################//
//     Reorder = ParMapper->GetReorder_D1();
    for(ii=0;ii<N_CDept1;ii++)
      {
	#pragma omp for schedule(guided) 
	for(jj=ptrCDept1[ii];jj<ptrCDept1[ii+1];jj++)
	{
	  i = Reorder_D1[jj];
	  if(i >= N_Active)     continue;
	  
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
	} // endfor jj
      } //end for ii  
//################################################################################################################//
    #pragma omp master
    {
      ParComm->CommUpdateH1(sol);
    }
//########################################## DEPENDENT2 DOFS #####################################################//
//make a variable to control the fraction of d2 dofs
//     Reorder = ParMapper->GetReorder_D2();
    for(ii=0;ii<N_CDept2;ii++)
    {
      #pragma omp for schedule(guided) 
      for(jj=ptrCDept2[ii];jj<ptrCDept2[ii+1];jj++)
      {
	i = Reorder_D2[jj];
	if(i >= N_Active)     continue;
    
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
      } // endfor jj
    } //end for ii
//################################################################################################################//


//########################################## MASTERS DOFS ########################################################//
    if(!lastTime)
    {
//       Reorder = ParMapper->GetReorder_M();
      for(ii=0;ii<N_CMaster;ii++)
      {
	#pragma omp for schedule(guided) 
	for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
	{
	  i = Reorder_M[jj];
          if(i >= N_Active)     continue;
    
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
        } // endfor jj
      } //end for ii
    } //end !lastTime
//################################################################################################################//
    if(!lastTime)
    #pragma omp master
    {
      ParComm->CommUpdateMS(sol);
    }
    
//########################################## INDEPENDENT DOFS ####################################################//  
//     Reorder = ParMapper->GetReorder_I();
    for(ii=0;ii<N_CInt;ii++)
    {
      #pragma omp for schedule(guided) 
      for(jj=ptrCInt[ii];jj<ptrCInt[ii+1];jj++)
      {
        i = Reorder_I[jj];
        if(i >= N_Active)     continue;
      
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
        } // endfor jj
        t = sol[i];
        sol[i] = omega*(s/diag-t) + t;
        // cout << "sol[i]: " << sol[i] << endl;
      } // endfor jj
    } //endfor ii
//################################################################################################################//

//############################################# Hanging NODES ####################################################//  
    // set hanging nodes
//     int *master = ParComm->GetMaster();

    #pragma omp for schedule(dynamic) nowait 

    for(i=N_Active;i<HangingNodeBound;i++)
    {

      if(master[i] != rank)
        continue;

      s = f[i];
      k = RowPtr[i+1];
      for(j=RowPtr[i];j<k;j++)
      {
        index = KCol[j];
        if(index != i)
          s -= Entries[j] * sol[index];
        else
          diag = Entries[j];
      } // endfor j
      sol[i] = s/diag;
    } // endfor i
  }
//################################################################################################################//
}

void TMGLevel3D::SOR_Re_Color_Coarse(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters)
{
//   SOR_Re(sol,f,aux,N_Parameters,Parameters);
//   return;
  
  int ii, i,j,jj,k,l,index,rank,tid,nrows=0,numThreads;
  double s, t, diag;
  double omega;
  int *master = ParComm->GetMaster();
  omega = Parameters[0];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  numThreads = TDatabase::ParamDB->OMPNUMTHREADS;
  omp_set_num_threads(numThreads);
  // set Dirichlet nodes
  memcpy(sol+HangingNodeBound, f+HangingNodeBound, 
           N_Dirichlet*SizeOfDouble);

#pragma omp parallel default(shared) private(i,ii,s,k,j,jj,tid,index,diag,t)
  {
//########################################## MASTERS DOFS ########################################################//
//       Reorder = ParMapper->GetReorder_M();
      for(ii=0;ii<N_CMaster;ii++)
      {
	#pragma omp for schedule(guided) 
	for(jj=ptrCMaster[ii];jj<ptrCMaster[ii+1];jj++)
	{
	  i = Reorder_M[jj];
          if(i >= N_Active)     continue;
    
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
        } // endfor jj
      } //end for ii
      
      #pragma omp master
      {
	ParComm->CommUpdateMS(sol);
      }
      #pragma omp barrier
//################################################################################################################//
    
//########################################## DEPENDENT1 DOFS #####################################################//
//     Reorder = ParMapper->GetReorder_D1();
    for(ii=0;ii<N_CDept1;ii++)
      {
	#pragma omp for schedule(guided) 
	for(jj=ptrCDept1[ii];jj<ptrCDept1[ii+1];jj++)
	{
	  i = Reorder_D1[jj];
	  if(i >= N_Active)     continue;
	  
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
	} // endfor jj
      } //end for ii  
//################################################################################################################//
    #pragma omp master
    {
      ParComm->CommUpdateH1(sol);
    }
//########################################## DEPENDENT2 DOFS #####################################################//
//     Reorder = ParMapper->GetReorder_D2();
    for(ii=0;ii<N_CDept2;ii++)
    {
      #pragma omp for schedule(guided) 
      for(jj=ptrCDept2[ii];jj<ptrCDept2[ii+1];jj++)
      {
	i = Reorder_D2[jj];
	if(i >= N_Active)     continue;
    
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
      } // endfor jj
    } //end for ii
//################################################################################################################//
    
//########################################## INDEPENDENT DOFS ####################################################//  
//     Reorder = ParMapper->GetReorder_I();
    for(ii=0;ii<N_CInt;ii++)
    {
      #pragma omp for schedule(guided) 
      for(jj=ptrCInt[ii];jj<ptrCInt[ii+1];jj++)
      {
        i = Reorder_I[jj];
        if(i >= N_Active)     continue;
      
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
        } // endfor jj
        t = sol[i];
        sol[i] = omega*(s/diag-t) + t;
        // cout << "sol[i]: " << sol[i] << endl;
      } // endfor jj
    } //endfor ii
//################################################################################################################//

//############################################# Hanging NODES ####################################################//  
    // set hanging nodes
//     int *master = ParComm->GetMaster();

    #pragma omp for schedule(dynamic) nowait 

    for(i=N_Active;i<HangingNodeBound;i++)
    {

      if(master[i] != rank)
        continue;

      s = f[i];
      k = RowPtr[i+1];
      for(j=RowPtr[i];j<k;j++)
      {
        index = KCol[j];
        if(index != i)
          s -= Entries[j] * sol[index];
        else
          diag = Entries[j];
      } // endfor j
      sol[i] = s/diag;
    } // endfor i
  }
//################################################################################################################//

  
}
#endif
