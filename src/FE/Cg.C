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
   
#include <ItMethod.h>
#include <Cg.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** constructor with initialization */
TCg::TCg(
MatVecProc *MatVec,
DefectProc *Defect,
TItMethod *Prec,
int n_aux,
int n_dof,
int scalar)
: TItMethod(MatVec, Defect, Prec, n_aux, n_dof)
{
  int i;
  double *aux;

  matvec = MatVec;
  matvecdefect = Defect;
  prec = Prec;
  N_DOF = n_dof;

  if (scalar)
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SCALAR;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SCALAR;

  }
  else
  {
    res_norm_min = TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SADDLE;
    red_factor= TDatabase::ParamDB->SC_LIN_RED_FACTOR_SADDLE;
    maxit = TDatabase::ParamDB->SC_LIN_MAXIT_SADDLE;
  }

  prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT;
  div_factor = TDatabase::ParamDB->SC_DIV_FACTOR;
  minit = TDatabase::ParamDB->SC_MINIT;

  if (prec_maxit<1)
  {
    OutPut("WARNING: Number of preconditioner iterations too small: "
      << prec_maxit << endl);
    OutPut("         Set number of preconditioner iterations to 1 !!!");
    prec_maxit = TDatabase::ParamDB->SC_AMG_PREC_IT = 1;
  }

  Ap = new double[N_DOF];
  p =  new double[N_DOF];
  r =  new double[N_DOF];
  z =  new double[N_DOF];
  r_last =  new double[N_DOF];

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
}


TCg::~TCg()
{
  delete[] Ap;
  delete[] p;
  delete[] r;
  delete[] z;
  delete[] r_last;

  if (N_Aux>0)
  {
    delete[] AuxArray;
  }
}


int TCg::Iterate (TSquareMatrix **sqmat,
TMatrix **mat, double *sol,
double *rhs)
{
  int maxite = maxit, i, j, k, verbose = TDatabase::ParamDB->SC_VERBOSE, ex_maxit = TDatabase::ParamDB->SC_EX_MAXIT;
  double t1, t2, dnorm, dnorm0, dnormlast, rho,rho_last=1.0,alpha, beta, sp;
  int flexible = TDatabase::ParamDB->SC_FLEXIBLE_KRYLOV_SPACE_SOLVER;

  // 1 - CG-Algorithm out of the Book of Saad
  // 0 - CG-Algo like AMG_solvers and paper of Notay
  int version = 0;

  t1 = GetTime();
  if (verbose>1)
    OutPut("Entering cg" << endl);

  matvecdefect(sqmat, mat, sol, rhs, r);
  /* norm of residual */
  dnorm0 = dnormlast = dnorm = sqrt(Ddot(N_DOF,r,r));

  if (dnorm <= res_norm_min )                     /* stopping criterion fulfilled */
  {
    if ((minit==0)||(dnorm<1e-20))
    {
      OutPut("(no) cg iteration " << 0 << " " << dnorm << endl);
      return(0);
    }
    else
    {
      maxite = minit;
    }
  }

  if (verbose>0)
    OutPut("cg Iteration " << 0 << " " << dnorm << endl);

  //CG-Algorithm out of the Book of Saad
  if(version)
  {
    //Algo wie in Buch
    if (prec==NULL)
    {
      Dcopy(N_DOF, r, z);
    }
    else
    {
      memset(z,0.0,N_DOF*SizeOfDouble);
      for (j=0; j<prec_maxit; j++)
      {
        prec->Iterate(sqmat, mat, z, r);
      }
    }
    Dcopy(N_DOF, z, p);
    for (i=0; i<maxite ; i++)
    {
      //rho
      rho=Ddot(N_DOF, r, z);
      //A*p
      matvec(sqmat, mat, p, Ap);
      //alpha=rho/(Ap,p)
      if(Ddot(N_DOF, p,Ap)!=0)
      {
        alpha=rho/Ddot(N_DOF, p,Ap);
      }
      else
      {
        Error("Error! (Ap,p) = " << Ddot(N_DOF, p,Ap) << endl);
        OutPut("Error in Cg.C !!!" << endl);
        exit(4711);
      }
      //update rho_last
      rho_last=rho;
      //x_j+1=x_j+alpha*p_j
      Daxpy(N_DOF, alpha, p, sol);
      if( flexible)
        //update r_last
        Dcopy(N_DOF, r, r_last);
      //r_j+1=r_j-alpha*Ap_j
      Daxpy(N_DOF, -alpha, Ap, r);
      //checks if residual is small enough already
      dnorm=sqrt(Ddot(N_DOF, r, r));
      if(minit<=i+1 && ex_maxit!=1)
      {
        if (dnorm<dnorm0*red_factor) break;
        if (dnorm<res_norm_min) break;
        if (dnorm>div_factor*dnorm0)
        {
          OutPut("CG iteration diverges !!!\n");
          exit(4711);
        }
      }
      if (verbose>0)
       {
         OutPut("CG Iteration " << i+1 << " dnorm " << dnorm << " " << " dnorm/dnormlast " << dnorm/dnormlast << "\n");
       }

       dnormlast=dnorm;
     //z_j+1=M^-1*r_j+1
      if (prec==NULL)
      {
        Dcopy(N_DOF, r,z);
      }
      else
      {
        memset(z,0.0,N_DOF*SizeOfDouble);
        for (j=0; j<prec_maxit; j++)
        {
          prec->Iterate(sqmat, mat, z, r);
        }
      }
      //update rho
      rho=Ddot(N_DOF, r, z);

      //p[i+1]=z_j+1+beta*p_j
      if(rho_last!=0)
      {
        if (flexible)
        {
          //beta=<z,(r-r_last)>/rho_last
          Dscal(N_DOF, -1.0, r_last);
          Daxpy(N_DOF, 1.0, r, r_last);
          beta=Ddot(N_DOF, r_last, z)/rho_last;
        }
        else
        {
          //beta=<z,r>/rho_last
          beta=rho/rho_last;
        }
      }
      else
      {
        Error("Error! (rho_last) = " << rho_last << endl);
        OutPut("Error in Cg.C !!!" << endl);
        exit(4711);
      }
      Dscal(N_DOF, beta, p);
      Daxpy(N_DOF, 1.0, z, p);

      if (i<minit) continue;
      if (ex_maxit) continue;
    }
    t2 = GetTime();
    if (i==maxit && !ex_maxit)
    {
      OutPut("solver not converged\n");
      OutPut("CG : (maximal) iterations " << i+1 << " residual " << dnorm << " reduction " << dnorm/dnorm0 << "\n");
      return(-1);
    }
    OutPut("CG : iterations " << i+1 << " residual " << dnorm << "\n");
    OutPut("Time: " << t2-t1 << "\n");
    return(i+1);

  }                                               //end of CG-Algorithm out of the Book of Saad

  //CG-Algo like AMG_solvers and paper of Notay
  else
  {
    OutPut("amg" <<  endl);
    //m needed for felxible
    if(flexible)
    {
      OutPut("flexible" <<  endl);
    }
    else
    {
      OutPut("not flexible" <<  endl);
    }
    memset(p,0.0,N_DOF*SizeOfDouble);
    for (i=0; i<maxite ; i++)
    {
      //z_j+1=M^-1*r_j+1
      if (prec==NULL)
      {
        Dcopy(N_DOF, r,z);
      }
      else
      {
        memset(z,0.0,N_DOF*SizeOfDouble);
        for (j=0; j<prec_maxit; j++)
        {
          prec->Iterate(sqmat, mat, z, r);
        }
      }
      if (flexible)
      {
        Dcopy(N_DOF, z, p);
        if(i==0)
        {
          Dcopy(N_DOF, z, r_last);
        }
        else
        {
          sp = Ddot(N_DOF, r_last, Ap);
          if(sp!=0)
          {
            beta=Ddot(N_DOF, z, Ap)/sp;
            Daxpy(N_DOF, -beta, r_last, p);
            Dcopy(N_DOF, p, r_last);
          }
          else
          {
            Error("Error! (r_last,Ap) = " << Ddot(N_DOF, r_last, Ap) << endl);
            OutPut("Error in Cg.C !!!" << endl);
            exit(4711);
          }
        }
      }
      else
      {
        //rho
        rho=Ddot(N_DOF, r, z);
        //p[i+1]=z_j+1+(rho/rho_last)*p_j
        if(rho_last!=0)
        {
          beta=rho/rho_last;
        }
        else
        {
          Error("Error! (rho_last) = " << rho_last << endl);
          OutPut("Error in Cg.C !!!" << endl);
          exit(4711);
        }
        Dscal(N_DOF, beta, p);
        Daxpy(N_DOF, 1.0, z, p);

        //update rho_last
        rho_last=rho;
      }
      //A*p
      matvec(sqmat, mat, p, Ap);

      //alpha=rho/(Ap,p)
      sp = Ddot(N_DOF, p,Ap);
      if(sp!=0)
      {
        if (flexible)
        {
          alpha=Ddot(N_DOF, p, r)/sp;
        }
        else
        {
          alpha=rho/sp;
        }
      }
      else
      {
        Error("Error! (Ap,p) = " << Ddot(N_DOF, p,Ap) << endl);
        OutPut("Error in Cg.C !!!" << endl);
        exit(4711);
      }

      //x_j+1=x_j+alpha*p_j
      Daxpy(N_DOF, alpha, p, sol);
      //r_j+1=r_j-alpha*Ap_j
      Daxpy(N_DOF, -alpha, Ap, r);
      //checks if residual is small enough already
      dnorm=sqrt(Ddot(N_DOF, r, r));
      if(minit<=i+1 && ex_maxit!=1)
      {
        if (dnorm<dnorm0*red_factor) break;
        if (dnorm<res_norm_min) break;
        if (dnorm>div_factor*dnorm0)
        {
          OutPut("CG iteration diverges !!!\n");
          exit(4711);
        }
      }

      if (verbose>0)
      {
        OutPut("CG Iteration " << i+1 << " dnorm " << dnorm << " " << " dnorm/dnormlast " << dnorm/dnormlast << "\n");
      }

      dnormlast=dnorm;
      if (i<minit) continue;
      if (ex_maxit) continue;
    }
    t2 = GetTime();
    if (i==maxit && !ex_maxit)
    {
      OutPut("solver not converged\n");
      OutPut("CG : (maximal) iterations " << i+1 << " residual " << dnorm << " reduction " << dnorm/dnorm0 << "\n");
      return(-1);
    }
    OutPut("CG : iterations " << i+1 << " residual " << dnorm << "\n");
    OutPut("Time: " << t2-t1 << "\n");
    return(i+1);

  }
}
