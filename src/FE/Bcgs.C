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
#include <Bcgs.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <LinAlg.h>
#include <NSE_MultiGrid.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** constructor with initialization */
TBcgs::TBcgs(MatVecProc *MatVec,
DefectProc *Defect,
TItMethod *Prec,
int n_aux, int n_dof,
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
  p =  new double[N_DOF];
  r =  new double[N_DOF];
  y =  new double[N_DOF];
  v = new double[N_DOF];
  t = new double[N_DOF];

  N_Aux = n_aux;
  if (n_aux>0)
  {
    AuxArray = new double* [n_aux];
    aux = new double[n_aux*N_DOF];
    for(i=0;i<n_aux;i++)
      AuxArray[i] = aux+i*N_DOF;
  }
}


TBcgs::~TBcgs()
{
  delete[] p;
  delete[] r;
  delete[] y;
  delete[] v;
  delete[] t;

  if (N_Aux>0)
  {
    delete[] AuxArray;
  }
}


int TBcgs::Iterate (TSquareMatrix **sqmat,
TMatrix **mat, double *sol,
double *rhs)
{

  int maxite = maxit, i, j, verbose = TDatabase::ParamDB->SC_VERBOSE, ex_maxit = TDatabase::ParamDB->SC_EX_MAXIT;
  double t1, t2, t3, t4, dnorm, dnorm0, dnormlast, rho,rho_last=1.0,alpha, beta, omega, delta;

  // 0 - not flexible
  // 1 - flexible
  int flexible = 0;

  t1 = GetTime();
  if (verbose>1)
    OutPut("Entering bcgs" << endl);

  if(flexible)
  {
    OutPut("flexible" <<  endl);
  }
  else
  {
    OutPut("not flexible" <<  endl);
  }

  matvecdefect(sqmat, mat, sol, rhs, r);
  //rhs = tilde_r_0
  Dcopy(N_DOF, r, rhs);
  dnorm=dnorm0=dnormlast=sqrt(Ddot(N_DOF, r,rhs));

  if (dnorm <= res_norm_min )                     /* stopping criterion fulfilled */
  {
    if ((minit==0)||(dnorm<1e-20))
    {
      OutPut("(no) bcgs iteration " << 0 << " " << dnorm << endl);
      return(0);
    }
    else
    {
      maxite = minit;
    }
  }

  if (verbose>0)
    OutPut("bcgs Iteration " << 0 << " " << dnorm << endl);

  for(i=0; i<maxite; i++)
  {
    //rho=(r,r_0)
    rho=Ddot(N_DOF, r,rhs);

    if(fabs(rho)<1.0E-50)
    {
      OutPut("BCGS break down!" << endl);
      return(-1);
    }
    if (i==0)
      Dcopy(N_DOF,r,p);
    else
    {
      //p=r+beta*p-beta*omega*v
      if(omega!=0)
      {
        beta=(rho/rho_last)*(alpha/omega);
      }
      else
      {
        Error("Error! omega = " << omega << endl);
        OutPut("Problem in Cg.C !!!" << endl);
        exit(4711);
      }
      Dscal(N_DOF,beta, p);
      Daxpy(N_DOF, 1.0, r, p);
      Daxpy(N_DOF, -beta*omega,v, p);
    }
    rho_last=rho;

    //y=Kp
    if (prec!=NULL)
    {
      memset(y,0.0,N_DOF*SizeOfDouble);
      for (j=0; j<prec_maxit; j++)
      {
        prec->Iterate(sqmat, mat, y, p);
      }
    }else
    Dcopy(N_DOF,p,y);
    if(flexible)
    {
      //v=Ay
      matvec(sqmat, mat, y, v);
    }
    else
    {
      //v=y
      Dcopy(N_DOF, y, v);
    }
    //delta=(v,r_0)
    delta=Ddot(N_DOF, rhs,v);
    if(fabs(delta)<1.0E-50) delta=1.0E-50;
    //alpha=rho/delta
    alpha=rho/delta;
    //overwrite r instead of using new variable s
    //s=r=r-alpha*v
    Daxpy(N_DOF, -alpha, v, r);
    //update first part of sol(x)
    if(flexible)
    {
      //x=x+alpha*y
      Daxpy(N_DOF, alpha, y, sol);
    }
    else
    {
      //x=x+alpha*p
      Daxpy(N_DOF, alpha, p, sol);
    }
    dnorm=sqrt(Ddot(N_DOF, r, r));
    if ( (dnorm<dnorm0*red_factor || dnorm<res_norm_min) && !ex_maxit )
    {
      if (verbose>0)
      {
        OutPut("BCGS Iteration " << i+1 << " res " << dnorm << " " << dnorm/dnormlast << endl);
      }
      break;
    }
    //y=Ks=Kr
    if(prec!=NULL)
    {
      //no need for y anymore so we reuse it for z
      memset(y,0.0,N_DOF*SizeOfDouble);
      for (j=0; j<prec_maxit; j++)
      {
        prec->Iterate(sqmat, mat, y, r);
      }
    }
    else
    {
      //no need for y anymore so we reuse it
      Dcopy(N_DOF, r, y);
    }
    if(flexible)
    {
      //t=Ay
      matvec(sqmat, mat, y, t);
    }
    else
    {
      //t=y
      Dcopy(N_DOF, y, t);
    }
    //delta=(t,t)
    delta=Ddot(N_DOF, t, t);
    if(fabs(delta)<1.0E-50)
      delta=1.0E-50;
    //omega=(t,s)/delta
    omega=Ddot(N_DOF, t, r)/delta;
    if(flexible)
      //x=x+omega*y
      Daxpy(N_DOF, omega, y, sol);
    else
      //x=x+omega*s
      Daxpy(N_DOF, omega, r, sol);
    //s=r
    //s=r=r-omega*t
    Daxpy(N_DOF, -omega, t, r);
    dnorm=sqrt(Ddot(N_DOF, r, r));
    if (verbose>0)
    {
      OutPut("BCGS Iteration " << i+1 << " dnorm " << dnorm << " " << " dnorm/dnormlast " << dnorm/dnormlast << endl);
    }
    dnormlast=dnorm;
    if (i<minit) continue;
    if (ex_maxit) continue;
    if (dnorm<dnorm0*red_factor) break;
    if (dnorm<res_norm_min) break;
    if (dnorm>div_factor*dnorm0)
    {
      OutPut("BCGS iteration diverges !!!"<<endl);
      exit(4711);
    }
  }

  t2 = GetTime();
  if (i==maxit && !ex_maxit)
  {
    OutPut("solver not converged\n");
    OutPut("BCGS : (maximal) iterations " << i+1 << " residual " << dnorm << " reduction " << dnorm/dnorm0 << endl);
    return(-1);
  }
  OutPut("BCGS : iterations " << i+1 << " residual " << dnorm << endl);
  //OutPut("Time: " << t2-t1 << "\n");
  return(i+1);

}
