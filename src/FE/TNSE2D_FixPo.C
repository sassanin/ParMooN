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
   
// ======================================================================
// @(#)TNSE2D_FixPo.C        1.3 05/05/00
//
// common declaration for all time dependent Navier-Stokes problems
// ======================================================================

#include <Database.h>
#include <Convolution.h>
#include <MooNMD_Io.h>

#include <stdlib.h>

// ======================================================================
// compute turbulent viscosity for LES
// ======================================================================
double TurbulentViscosity(double delta, double* gradU, double* u, double* uConv)
{
  int nu_type = TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE;
  double nu_constant = TDatabase::ParamDB->TURBULENT_VISCOSITY_CONSTANT;
  int nu_tensor =  TDatabase::ParamDB->TURBULENT_VISCOSITY_TENSOR;
  double nu_power, nu_sigma;
  double frobenius_norm_tensor,nu,a,b,c,sigma;

  // compute turbulent viscosity
  switch(nu_type)
  {
      // no turbulent viscosity
      case 0:
	  nu = 0;
	  return(nu);
	  break;
      case 5:
	  nu = nu_constant;
	  return(nu);
	  break;
	  // defect correction  
      case 6:
	  if (TDatabase::ParamDB->DEFECT_CORRECTION_TYPE)	  
	      nu = delta;
	  else
	      nu = 0;  // bwe
	  return(nu);
	  break;
  }

  // compute square of the Frobenius norm of the tensor
  // use deformation tensor
  if (nu_tensor==0)
  {
    // compute (grad(u)+grad(u)^T)/2
    a = gradU[0]+gradU[0];
    b = gradU[1]+gradU[2];
    c = gradU[3]+gradU[3];
    frobenius_norm_tensor = (a*a+ 2*b*b + c*c)/4.0;
  }
  // use grad u
  else
    frobenius_norm_tensor =  gradU[0]*gradU[0] + gradU[1]*gradU[1] +
      gradU[2]*gradU[2] + gradU[3]*gradU[3];

  // compute turbulent viscosity
  switch(nu_type)
  {
    case 1:                      // Smagorinsky
    case 17:                     // Smagorinsky
      nu = nu_constant * delta * delta * sqrt(frobenius_norm_tensor);
      break;
    case 2:                      // p laplacian
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power/2.0);
      nu = nu_constant * delta * delta * pow(frobenius_norm_tensor, nu_power-2.0);
      break;
    case 3:                      // Layton SIAM J. Sci. Comput. 1996
      // nu = nu_0 |ln(h)|^(-0.5*(p-1)) h^sigma \|\nabla u\|^(p-2)
      nu_power = TDatabase::ParamDB->TURBULENT_VISCOSITY_POWER;
      nu_sigma = TDatabase::ParamDB->TURBULENT_VISCOSITY_SIGMA;
      nu = nu_constant * pow(delta,nu_sigma) *
        pow(frobenius_norm_tensor,(nu_power-2)/2.0)
        / pow(fabs(log(delta)),0.5*(nu_power-1));
      break;
    case 4:                      // \|u-g_\delta \ast u\|_2
      nu = nu_constant * delta * sqrt((u[0]-uConv[0])*(u[0]-uConv[0])+
        (u[1]-uConv[1])*(u[1]-uConv[1]));
      //   OutPut(u[0] << " " << uConv[0] << "        "<< endl);
      //    OutPut(u[1] << " " << uConv[1] << endl);
      break;
    default:
      OutPut("This type of turbulent viscosity is not implemented !!!" << endl);
      exit(4711);
      break;
  }
  return(nu);
}

/******************************************************************************/
//
// computation of SUPG parameter following 
// Bazilevs, Calo, Cottrell, Hughes, Reali, Scovazzi
//
/******************************************************************************/

void SUPG_Param2D(double Mult, double* u, double* coeff, double* params)
{
    double x0, x1, x2, y0, y1, y2, g11, g12, g22;
    double d11, d12, d21, d22, delta, nu, tau_c, tau_m;
    double u1, u2, rec_detjk;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;    
    double C_I = TDatabase::ParamDB->DELTA0;

    nu = coeff[0];        
    rec_detjk = coeff[19];
    rec_detjk = 1/rec_detjk;
    u1 = u[0];
    u2 = u[1];
    
    x0 = TDatabase::ParamDB->INTERNAL_VERTEX_X[0];
    y0 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[0];
    x1 = TDatabase::ParamDB->INTERNAL_VERTEX_X[1];
    y1 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[1];
    x2 = TDatabase::ParamDB->INTERNAL_VERTEX_X[2];
    y2 = TDatabase::ParamDB->INTERNAL_VERTEX_Y[2];
    
    // triangle
    if (TDatabase::ParamDB->INTERNAL_VERTEX_X[3]== -4711)
    {
	d11 = (y2-y1) * rec_detjk;  //dxi/dx
	d12 = (x1-x2) * rec_detjk;  //dxi/dy
	d21 = (y0-y1) * rec_detjk;  //deta/dx
	d22 = (x1-x0) * rec_detjk;  //deta/dy
    }
    else
    {
	// quadrilateral
	d11 = (y2-y1) * 0.5 * rec_detjk;  //dxi/dx
	d12 = (x1-x2) * 0.5 * rec_detjk;  //dxi/dy
	d21 = (y0-y1) * 0.5 * rec_detjk;  //deta/dx
	d22 = (x1-x0) * 0.5 * rec_detjk;  //deta/dy
    }
	
    g11 = d11*d11 + d21*d21;
    g12 = d11*d12 + d21*d22;
    g22 = d12*d12 + d22*d22;

    tau_m = g11*g11 + 2*g12*g12 + g22*g22; // G:G
    tau_m *= C_I*nu*nu;
    tau_m +=  4/(time_step*time_step); 
    tau_m += u1 * (g11*u1+g12*u2) + u2*(g12*u1+g22*u2);
    tau_m = 1/sqrt(tau_m); // this is the parameter for the momentum equation
 
    tau_c = (d11+d21)*(d11+d21)+(d12+d22)*(d12+d22);
    tau_c *= tau_m;
    tau_c = 1/tau_c;

    params[0] = tau_m;
    params[1] = tau_c;

/*    delta = (d11*d11+d21*d21)*(d11*d11+d21*d21)+2*(d11*d12+d21*d22)*(d11*d12+d21*d22)+  // G:G
	(d12*d12+d22*d22)*(d12*d12+d22*d22);
    delta *= C_I*nu*nu;         
    delta += 4/(time_step*time_step);  
    delta += u1*u1*(d11*d11+d21*d21)+2*u1*u2*(d11*d12+d21*d22)+u2*u2*(d12*d12+d22*d22);  // uGu
    delta = 1/sqrt(delta); // this is the parameter for the momentum equation
    
    delta *= ((d11+d21)*(d11+d21)+(d12+d22)*(d12+d22));   // gg
    delta = 1/delta; // this is the parameter for the mass equation
*/    
}


// ======================================================================
// Type 1, Standard Galerkin
// Type 1, Coletti
// Type 1, GL00Convolution
// ======================================================================
void TimeNSType1Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// for Type 1 is no SDFEM available
// ======================================================================

// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky
// the nonlinear viscosity is treated implicitly
// ======================================================================
void TimeNSType1Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  
  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 1, Galdi-Layton 98 Model auxiliary problem
// ======================================================================
void TimeNSType1GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM, **AuxMatrix;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, mu2, delta, mu;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}

// ======================================================================
// Type 1, for group fem (the matrices for the convective term)
// ======================================================================
void TimeNSType1GroupFEM(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixCx, **MatrixCy;
  double val;
  double *MatrixRow1, *MatrixRow2;
  double ansatz10, ansatz01;
  double test00;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;

  MatrixCx = LocMatrices[0];
  MatrixCy = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  for(i=0;i<N_U;i++)
  {
    MatrixRow1 = MatrixCx[i];
    MatrixRow2 = MatrixCy[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = ansatz10*test00;
      MatrixRow1[j] += Mult * val;

      val = ansatz01*test00;
      MatrixRow2[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 2, Standard Galerkin
// Type 2, Coletti
// Type 2, GL00Convolution
// ======================================================================
void TimeNSType2Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// for Type 2 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 2, Upwind (only Laplacian in A block)
// ======================================================================
void TimeNSType2Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, Smagorinsky
// ======================================================================
void TimeNSType2Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[2];
  MatrixB2 = LocMatrices[3];
  MatrixB1T = LocMatrices[4];
  MatrixB2T = LocMatrices[5];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 2, GL00AuxProblem
// ======================================================================
void TimeNSType2GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA, **MatrixB1, **MatrixB2, **MatrixM, **AuxMatrix;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow, *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ngu, val1, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA = LocMatrices[0];
  MatrixM = LocMatrices[1];
  MatrixB1 = LocMatrices[3];
  MatrixB2 = LocMatrices[4];
  MatrixB1T = LocMatrices[5];
  MatrixB2T = LocMatrices[6];
  AuxMatrix = LocMatrices[2];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    MatrixMRow = MatrixM[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      MatrixRow[j] += Mult * val;

      val = ansatz00*test00;
      MatrixMRow[j] += Mult * val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v)
// Type 3, Coletti, (grad u, grad v)
// Type 3, GL00Convolution, (grad u, grad v)`
// ======================================================================
void TimeNSType3Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v)
// Type 3, Coletti, D(u):D(v)
// Type 3, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType3GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// for Type 3 is SDFEM not available
// ======================================================================

// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = Mult*c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(2* test10*ansatz10+ test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType3Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType3SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, viscosity, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[6];
  MatrixB2  = LocMatrices[7];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, GL00AuxProblem (grad u, grad v)
// ======================================================================
void TimeNSType3GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *AuxMatrixRow;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double delta, ngu, val1, mu, mu2;
  double u1, u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType3GL00AuxProblemDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, viscosity, val1, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 3, VMSProjection, D(u):D(v)
// ======================================================================
void TimeNSType3VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, viscosity, delta;

  //cout << "vms";

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixL   = LocMatrices[6];
  MatrixB1  = LocMatrices[7];
  MatrixB2  = LocMatrices[8];
  Matrix_tilde_G11  = LocMatrices[9];
  Matrix_tilde_G22  = LocMatrices[10];
  Matrix_G11  = LocMatrices[11];
  Matrix_G22  = LocMatrices[12];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;
  viscosity = c0+mu;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}


// ======================================================================
// Type 4, Standard Galerkin, (grad u, grad v)
// Type 4, Coletti, (grad u, grad v)
// Type 4, GL00Convolution, (grad u, grad v)
// ======================================================================
void TimeNSType4Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v)
// Type 4, Coletti, D(u):D(v)
// Type 4, GL00Convolution, D(u):D(v)
// ======================================================================
void TimeNSType4GalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow, **AuxMatrix;;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

//   if(fabs(u1)>0)
// cout << " u1 " << u1 << " u2 " << u2 <<endl;
// exit(0);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// for Type 4 SDFEM is not available
// ======================================================================

// ======================================================================
// Type 4, Upwind (no convection terms), (grad u, grad v)
// ======================================================================
void TimeNSType4Upwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      // val  = 0;
      // Matrix12Row[j] += Mult * val;

      // val  = 0;
      // Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Upwind (no convection terms), D(u):D(v)
// ======================================================================
void TimeNSType4UpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, (grad u, grad v)
// ======================================================================
void TimeNSType4Smagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, Smagorinsky, D(u):D(v)
// ======================================================================
void TimeNSType4SmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, GL00AuxProblem, (grad u, grad v)
// ======================================================================
void TimeNSType4GL00AuxProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double *AuxMatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, ngu, val1, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, GL00AuxProblem, D(u):D(v)
// ======================================================================
void TimeNSType4GL00AuxProblemDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22, **AuxMatrix;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *AuxMatrixRow;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2;
  double delta, val1, mu, mu2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  AuxMatrix = LocMatrices[6];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    AuxMatrixRow = AuxMatrix[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j

  }                              // endfor i
}


// ======================================================================
// Type 4, VMSProjection, D(u):D(v)
// ======================================================================
void TimeNSType4VMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixL, **Matrix_tilde_G11 ;
  double **Matrix_tilde_G22, **Matrix_G11, **Matrix_G22;
  double *Rhs1, *Rhs2, val;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4;
  int i,j, N_U, N_P, N_L;
  double c0, c1, c2;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixL   = LocMatrices[6];
  MatrixB1 = LocMatrices[7];
  MatrixB2 = LocMatrices[8];
  MatrixB1T = LocMatrices[9];
  MatrixB2T = LocMatrices[10];
  Matrix_tilde_G11  = LocMatrices[11];
  Matrix_tilde_G22  = LocMatrices[12];
  Matrix_G11  = LocMatrices[13];
  Matrix_G22  = LocMatrices[14];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // l

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);
  mu = mu/2.0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (c0+mu)*(2*test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;

      val  = (c0+mu)*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = (c0+mu)*(test10*ansatz10+2*test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;
    }                            // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row  = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }

  for(i=0;i<N_L;i++)
  {
    Matrix11Row = Matrix_G11[i];
    Matrix22Row = Matrix_G22[i];
    test00 = Orig4[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      Matrix11Row[j] -= Mult * ansatz10 * test00;
      Matrix22Row[j] -= Mult * ansatz01 * test00;
    }
  }

  for(i=0;i<N_L;i++)
  {
    test00 = Orig4[i];
    MatrixRow1 = MatrixL[i];
    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig4[j];
      MatrixRow1[j] += Mult * ansatz00 * test00;
    }
  }
}


// ======================================================================
// assemble matrix for auxiliary problem
// ======================================================================

void MatrixAuxiliaryProblem(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *AuxMatrixRow, **AuxMatrix;;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double delta, mu2, val;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  // solution does not need to be convolved
  AuxMatrix = LocMatrices[0];
  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  // delta^2/(4 gamma)
  mu2 = 0.25*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    AuxMatrixRow = AuxMatrix[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = mu2*(test10*ansatz10+test01*ansatz01);
      val += ansatz00*test00;
      AuxMatrixRow[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Assembling routine for all nonlinear matrices
// ======================================================================

// ======================================================================
// Type 1, Standard Galerkin, only nonlinear part
// Type 2, Standard Galerkin, only nonlinear part
// Type 1, Coletti, only nonlinear part
// Type 2, Coletti, only nonlinear part
// Type 1, GL00Convolution, only nonlinear part
// Type 2, GL00Convolution, only nonlinear part
// Type 1, GL00AuxProblem, only nonlinear part
// Type 2, GL00AuxProblem, only nonlinear part
// ======================================================================
void TimeNSType1_2NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, for upwind (only laplacian in A block)
// Type 2, for upwind (only laplacian in A block)
// ======================================================================
void TimeNSType1_2NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_U;
  double c0;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = c0*(test10*ansatz10+test01*ansatz01);

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 1, Smagorinsky, only nonlinear part
// Type 2, Smagorinsky, only nonlinear part
// ======================================================================
void TimeNSType1_2NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;

      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 4, Standard Galerkin, (grad u, grad v), only nonlinear part
// Type 3, Coletti, (grad u, grad v), only nonlinear part
// Type 4, Coletti, (grad u, grad v), only nonlinear part
// Type 3, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 4, GL00Convolution, (grad u, grad v), only nonlinear part
// Type 3, GL00AuxProblem, (grad u, grad v), only nonlinear part
// Type 4, GL00AuxProblem, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLGalerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      Matrix11Row[j] += Mult * val;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// + convection with higher order velocity
// ======================================================================
void TimeNSType3_4NLGalerkin_VMS_1_DD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2, ho_u1, ho_u2;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  ho_u1 = param[0];              // higher order u1old
  ho_u2 = param[1];              // higher order u2old
  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = ((u1+ho_u1)*ansatz10+(u2+ho_u2)*ansatz01)*test00;
      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), (grad u, grad v)
// Type 4, Upwind (no convection term), (grad u, grad v)
// ======================================================================
void TimeNSType3_4NLUpwind(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j, N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = Mult*c0*(test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Upwind (no convection term), D(u):D(v)
// Type 4, Upwind (no convection term), D(u):D(v)
// ======================================================================
void TimeNSType3_4NLUpwindDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test10, test01;
  double *Orig0, *Orig1;
  int i,j,N_U;
  double c0;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  c0 = coeff[0];                 // nu

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = c0*(2*test10*ansatz10+test01*ansatz01);
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01);
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, (grad u, grad v), only nonlinear part
// Type 4, Smagorinsky, (grad u, grad v), only nonlinear part
// ======================================================================
void TimeNSType3_4NLSmagorinsky(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0;
  double u1, u2, mu, delta;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = (c0+mu)*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val *= Mult;
      Matrix11Row[j] += val;
      Matrix22Row[j] += val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Smagorinsky, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLSmagorinskyDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j,N_U;
  double c0, viscosity, delta;
  double u1, u2, mu;
  // cout << "Sma" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);
  mu = mu/2.0;
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Type 3, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// Type 4, VMSProjection, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLVMSProjectionDD(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double **Matrix_tilde_G11,  **Matrix_tilde_G22;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_L;
  double c0, viscosity, delta;
  double u1, u2, mu;
  // cout << "Sma" << endl;
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  Matrix_tilde_G11  = LocMatrices[4];
  Matrix_tilde_G22  = LocMatrices[5];

  N_U = N_BaseFuncts[0];
  N_L = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[6]);
  mu = mu/2.0;
  viscosity = mu+c0;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00;
      val  = viscosity*(2*test10*ansatz10+test01*ansatz01);
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = viscosity*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = viscosity*(test10*ansatz10+2*test01*ansatz01);
      val += val1;
      Matrix22Row[j] += Mult * val;

    }                            // endfor j
  }                              // endfor i

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = Matrix_tilde_G11[i];
    Matrix22Row = Matrix_tilde_G22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_L;j++)
    {
      ansatz00 = Orig3[j];
      Matrix11Row[j] -= Mult * 2*mu * ansatz00 * test10;
      Matrix22Row[j] -= Mult * 2*mu * ansatz00 * test01;
    }
  }
}


// ======================================================================
// ROSENBROCK
// ======================================================================
void TimeNSType1GalerkinRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                              // endfor i
}


void TimeNSType1GalerkinJ(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA;
  double val;
  double *MatrixRow, *MatrixRow1, *MatrixRow2, *MatrixMRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2;
  double u1, u2;

  MatrixA = LocMatrices[0];

  N_U = N_BaseFuncts[0];
  //N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  //Orig3 = OrigValues[3]; // p

  c0 = coeff[0];                 // nu
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    MatrixRow = MatrixA[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      //val += (ansatz00)
      MatrixRow[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


void TimeNSType1GalerkinC(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j,N_U, N_P;
  double c0, c1, c2, c3, c4;
  double u1, u2;

  //cout << "TimeNSType1GalerkinC" << endl;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig2 = OrigValues[2];         // u

  c3 = coeff[3];                 // dot f1
  c4 = coeff[4];                 // dot f2

  for(i=0;i<N_U;i++)
  {
    // test10 = Orig0[i];
    // test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c3;
    Rhs2[i] += Mult*test00*c4;
  }
}


void TimeNSType3GalerkinJ(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0;
  double u1, u2, u1_x, u1_y, u2_x, u2_y;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c0 = coeff[0];                 // nu

  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1_x = param[2];               // u1old
  u2_x = param[3];               // u2old
  u1_y = param[4];               // u1old
  u2_y = param[5];               // u2old

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += u1_x*ansatz00*test00;
      Matrix11Row[j] += Mult * val;

      val  = u1_y*ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = u2_x*ansatz00*test00;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+test01*ansatz01);
      val += (u1*ansatz10+u2*ansatz01)*test00;
      val += u2_y*ansatz00*test00;
      Matrix22Row[j] += Mult * val;
    }                            // endfor j
  }                              // endfor i
}


// ======================================================================
// Assembling routine for right-hand sides only
// ======================================================================

// ======================================================================
// right-hand side ONLY, for NSE
// ======================================================================
void TimeNSRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY for auxiliary problem applied to velocity
// ======================================================================
void TimeNSRHSAuxProblemU(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00;
  double *Orig0;
  int i, N_U;
  double c1, c2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  c1 = param[0];                 // f1
  c2 = param[1];                 // f2

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    //cout <<  Rhs1[i] << " " <<  Rhs2[i] << " ";
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Coletti model
// ======================================================================
void TimeNSRHSColetti(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2;

  double delta, ngu, val1, mu1;
  double D1u1, D2u1, D1u2, D2u2;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // (delta^2)/(2 gamma)
  mu1 = 0.5*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    // LES term
    val1 =  Mult* mu1;
    Rhs1[i] += val1*( test10* (D1u1*D1u1 + D2u1*D2u1)
      +test01* (D1u1*D1u2 + D2u1*D2u2));
    Rhs2[i] += val1*( test10* (D1u2*D1u1 + D2u2*D2u1)
      +test01* (D1u2*D1u2 + D2u2*D2u2));
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Galdi-Layton model with convolution
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSRHSLESModel(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test10, test01;
  double *Orig0, *Orig1;
  int i, N_U;
  double delta, val1, mu1;
  double gdT11, gdT12, gdT22;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y

  gdT11 = param[0];              // gDeltaT11 or 11 - component of auxiliary problem
  gdT12 = param[1];              // gDeltaT12 or 12 - and 21 - component
  gdT22 = param[2];              // gDeltaT22 or 22 - component of auxiliary problem

  // filter width
  delta =  CharacteristicFilterWidth(hK);
  mu1 = 0.5*delta*delta/gamma;
  val1 =  Mult* mu1;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];

    // LES term
    Rhs1[i] += val1*( test10* gdT11 + test01* gdT12);
    Rhs2[i] += val1*( test10* gdT12 + test01* gdT22);
  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSRHSGL00AuxProblemPaper2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2;

  double delta, ngu, val1, mu, mu1;
  double D1u1, D2u1, D1u2, D2u2;
  double AuxProblem11, AuxProblem12, AuxProblem22;
  double AuxProblem11Exact_x, AuxProblem12Exact_x, AuxProblem12Exact_y;
  double AuxProblem22Exact_y;
  double gamma = TDatabase::ParamDB->GAUSSIAN_GAMMA;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  D1u1 = param[0];               // D1u1
  D1u2 = param[1];               // D1u2;
  D2u1 = param[2];               // D2u1
  D2u2 = param[3];               // D2u2;
  AuxProblem11 = param[4];       // 11 - component of auxiliary problem
  AuxProblem12 = param[5];       // 12 - and 21 - component
  AuxProblem22 = param[6];       // 22 - component of auxiliary problem
  AuxProblem11Exact_x = param[7];// 11 - component of auxiliary problem
  AuxProblem12Exact_x = param[8];// 12 - and 21 - component
  AuxProblem12Exact_y = param[9];// 12 - and 21 - component
                                 // 22 - component of auxiliary problem
  AuxProblem22Exact_y = param[10];

  // filter width
  delta =  CharacteristicFilterWidth(hK);

  // turbulent viscosity
  mu = TurbulentViscosity(delta,&param[0],&param[7],&param[9]);

  // (delta^2)/(2 gamma)
  mu1 = 0.5*delta*delta/gamma;

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    // divergence of term with exact solution
    val1 =  Mult*mu1;
    Rhs1[i] -= val1*test00*(AuxProblem11Exact_x+AuxProblem12Exact_y);
    Rhs2[i] -= val1*test00*(AuxProblem12Exact_x+AuxProblem22Exact_y);

    // explicit treatment of artificial viscosity
    //val =  Mult *mu;
    //Rhs1[i] -= val*(D1u1*test10 + D2u1*test01);
    //Rhs2[i] -= val*(D1u2*test10 + D2u2*test01);

    val1 =  Mult*mu1;
    Rhs1[i] += val1*( test10* AuxProblem11 + test01* AuxProblem12);
    Rhs2[i] += val1*( test10* AuxProblem12 + test01* AuxProblem22);

  }                              // endfor i
}


// ======================================================================
// right-hand side for auxiliary problem
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test00;
  double *Orig0;
  int i, N_U;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  D1u1 = param[0];               // D1u1
  D1u2 = param[1];               // D1u2;
  D2u1 = param[2];               // D2u1
  D2u2 = param[3];               // D2u2;

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    val = Mult*test00;

    Rhs1[i] += val*(D1u1*D1u1+D2u1*D2u1);
    Rhs2[i] += val*(D1u1*D1u2+D2u1*D2u2);
    Rhs3[i] += val*(D1u2*D1u2+D2u2*D2u2);

  }                              // endfor i
}


// ======================================================================
// right-hand side for auxiliary problem
// ONLY, Galdi-Layton model with auxiliary problem
// ======================================================================
void TimeNSGL00AuxProblemRHSPaper2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, val;
  double test00;
  double *Orig0;
  int i, N_U;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    val = Mult*test00;

    Rhs1[i] += val*coeff[0];
    Rhs2[i] += val*coeff[1];
    Rhs3[i] += val*coeff[2];

  }                              // endfor i
}


// ======================================================================
// right-hand side ONLY, Smagorinsky Explicit
// ======================================================================
void TimeNSRHSSmagorinskyExplicit(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2, delta;

  double mu;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;

  // turbulent viscosity
  delta =  CharacteristicFilterWidth(hK);
  mu = TurbulentViscosity(delta,&param[2],&param[0],&param[0]);

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;

    val =  Mult * mu;
    Rhs1[i] -= val*(D1u1*test10 + D2u1*test01);
    Rhs2[i] -= val*(D1u2*test10 + D2u2*test01);

  }                              // endfor i
}


// ======================================================================
// right-hand side for additional terms in rhs of small scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_SmallRhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure
  p = param[12];

  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(D1u1*test10+(D1u2+D2u1)*test01/4);
    Rhs1[i] += (u1*D1u1+u2*D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(D2u2*test01+(D1u2+D2u1)*test10/4);
    Rhs2[i] += (u1*D1u2+u2*D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}


// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS
// ======================================================================
void TimeNS_VMS_Large_0_Rhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  // OutPut("l");
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure (small scales)
  p = param[12];
  //OutPut(p);
  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(ho_D1u1*test10+(ho_D1u2+ho_D2u1)*test01/4);
    Rhs1[i] += (ho_u1*ho_D1u1+ho_u2*ho_D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(ho_D2u2*test01+(ho_D1u2+ho_D2u1)*test10/4);
    Rhs2[i] += (ho_u1*ho_D1u2+ho_u2*ho_D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}


// ======================================================================
// right-hand side for additional terms in rhs of large scale systems
// for VMS, Variant 1
// ======================================================================
void TimeNS_VMS_Large_1_Rhs2D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  double u1, u2, ho_u1, ho_u2;
  double D1u1, D2u1, D1u2, D2u2;
  double ho_D1u1, ho_D2u1, ho_D1u2, ho_D2u2;
  double p, c0;
  int i, N_U;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  // large scales
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  // small scales
  ho_u1 = param[6];
  ho_u2 = param[7];
  ho_D1u1 = param[8];            // D1u1
  ho_D1u2 = param[9];            // D1u2;
  ho_D2u1 = param[10];           // D2u1
  ho_D2u2 = param[11];           // D2u2;
  // pressure (small scales)
  p = param[12];
  //OutPut(p);
  c0 = coeff[0];

  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += 2*c0*(ho_D1u1*test10+(ho_D1u2+ho_D2u1)*test01/4);
    Rhs1[i] += (ho_u1*ho_D1u1+ho_u2*ho_D2u1) * test00;
    Rhs1[i] += (u1*ho_D1u1+u2*ho_D2u1) * test00;
    //    Rhs1[i] += (ho_u1*D1u1+ho_u2*D2u1) * test00;
    Rhs1[i] -= p * test10;
    Rhs1[i] *= Mult;
    Rhs2[i] += 2*c0*(ho_D2u2*test01+(ho_D1u2+ho_D2u1)*test10/4);
    Rhs2[i] += (ho_u1*ho_D1u2+ho_u2*ho_D2u2) * test00;
    Rhs2[i] += (u1*ho_D1u2+u2*ho_D2u2) * test00;
    //    Rhs2[i] += (ho_u1*D1u2+ho_u2*D2u2) * test00;
    Rhs2[i] -= p * test01;
    Rhs2[i] *= Mult;
  }                              // endfor i
}

// ======================================================================
// right-hand side ONLY, defect correction type 1, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2, delta;

  double mu;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2

  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  
  //cout << D1u1 << " " << D1u2  << " " << D2u1  << " " << D2u2 <<endl;
  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs1[i] += Mult*hK*(D1u1*test10 + D2u1*test01);
    Rhs2[i] += Mult*test00*c2;
    Rhs2[i] += Mult*hK*(D1u2*test10 + D2u2*test01);
  }                              // endfor i
}
// ======================================================================
// right-hand side ONLY, defect correction type 2, u2
// ======================================================================
void TimeNSRHSDefectCorrectionU2_1(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, val;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i, N_U;
  double c1, c2, delta;
  double dt = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double mu, u1, u2, um11, um12, um21, um22;
  double D1u1, D2u1, D1u2, D2u2;

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u

  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  
  u1 = param[0];
  u2 = param[1];
  D1u1 = param[2];               // D1u1
  D1u2 = param[3];               // D1u2;
  D2u1 = param[4];               // D2u1
  D2u2 = param[5];               // D2u2;
  um11 = param[6];  
  um12 = param[7];
  um21 = param[8];
  um22 = param[9];

  //OutPut( u1<< " " << u2 << " " << um11 << " " << um12  << " " << um21  << " " << um22 <<endl);
  for(i=0;i<N_U;i++)
  {
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs1[i] += Mult*hK*(D1u1*test10 + D2u1*test01);
    Rhs1[i] -= Mult *(u1-2*um11+um21)*test00/(2*dt);
    Rhs2[i] += Mult*test00*c2;
    Rhs2[i] += Mult*hK*(D1u2*test10 + D2u2*test01);
    Rhs2[i] -= Mult *(u2-2*um12+um22)*test00/(2*dt);
  }                              // endfor i
}

// ======================================================================
// Type 4, Standard Galerkin, D(u):D(v) for Axialsymmetric
// ======================================================================
void TimeNSType4GalerkinDD_Axial3D(double Mult, double *coeff,
                                   double *param, double hK,
                                   double **OrigValues, int *N_BaseFuncts,
                                   double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow, **AuxMatrix;;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2;
  double u1, u2, r, x, sign;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // nu
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  u1 = param[0]; // u1old - grid VX   !!check aux parameters
  u2 = param[1]; // u2old - grid VY   !!check aux parameters

  x  = param[2]; // x
  r  = fabs(x);

  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

// exit(0);

  if (x>=0)
     sign = 1;
  else
     sign = -1;
//      sign = 1;    // unles all integral oints are inside it is not feasible
// if(fabs(u2)>1e-2)
// cout << " u1 " << u1<<  " u2 " << u2<< endl;
// cout << "Mult" << Mult<< endl;
// cout << "f1 :" << c1<<  " f2 :" << c2<< endl;
// cout << "c0 :" << c0<<  " Mult :" << Mult<< endl;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i]*sign;
    test01 = Orig1[i];
    test00 = Orig2[i];


     Rhs1[i] += Mult*test00*c1*r;
     Rhs2[i] += Mult*test00*c2*r;


    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r   + 2.*ansatz00*test00/r) ;
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2.*test01*ansatz01)*r;
      val += (u1*ansatz10+u2*ansatz01)*test00*r;
      Matrix22Row[j] += Mult * val;

      val = Mult*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r + ansatz00*test00);
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i
}

// ======================================================================
//  for axial symmetric case
// Type 3, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Standard Galerkin, D(u):D(v), only nonlinear diagonal blocks
// Type 3, Coletti, D(u):D(v), only nonlinear diagonal blocks
// Type 4, Coletti, D(u):D(v), only nonlinear diagonal blocks
// ======================================================================
void TimeNSType3_4NLGalerkinDD_Axial3D(double Mult, double *coeff,
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz10, ansatz01, ansatz00;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0;
  double u1, u2, r, x, sign;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old

  x  = param[2]; // x
  r  = fabs(x);

  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }

  if (x>=0)
     sign = 1;
  else
     sign = -1;
//      sign = 1;    // unles all integral oints are inside it is not feasible
 // cout << "u1" <<u1 << endl;
 // cout << "u2" <<u2 << endl;

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i]*sign;
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j]*sign;
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00*r;
      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r + 2.*ansatz00*test00/r );
      val += val1;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += val1;
      Matrix22Row[j] += Mult * val;

    } // endfor j
  } // endfor i
}


// ======================================================================
// Type 4, Two phase Standard Galerkin, D(u):D(v)
// Type 4, Two phase Coletti, D(u):D(v)
// Type 4, Two phase GL00Convolution, D(u):D(v)
// ======================================================================
// /*
void TimeNSType4GalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixM11, **MatrixM22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double *Rhs1, *Rhs2, val;
  double *AuxMatrixRow, **AuxMatrix;;
  double *Matrix11Row, *Matrix12Row, *Matrix21Row, *Matrix22Row;
  double *MatrixM11Row, *MatrixM22Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i,j, N_U, N_P;
  double c0, c1, c2, c3, c4, c5;
  double u1, u2, x, r;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  MatrixM11 = LocMatrices[4];
  MatrixM22 = LocMatrices[5];
  MatrixB1 = LocMatrices[6];
  MatrixB2 = LocMatrices[7];
  MatrixB1T = LocMatrices[8];
  MatrixB2T = LocMatrices[9];

  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u
  Orig3 = OrigValues[3]; // p

  c0 = coeff[0]; // 1/Re_k
  c1 = coeff[1]; // f1
  c2 = coeff[2]; // f2

  c3 = coeff[3]; // density ratio (inner/outer)

  u1 = param[0]; // u1old - grid VX   !!check aux parameters
  u2 = param[1]; // u2old - grid VY   !!check aux parameters

  x  = param[2]; // x
  
  if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }


  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixM11[i];
    MatrixM22Row  = MatrixM22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1* c3*r;
    Rhs2[i] += Mult*test00*c2* c3*r;

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = c0*((2*test10*ansatz10+test01*ansatz01)*r);
      if(TDatabase::ParamDB->Axial3D==1)
      {
       val += c0*2.*ansatz00*test00/r;
      }
      val += (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      Matrix11Row[j] += Mult *val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      Matrix22Row[j] += Mult * val;

      val = Mult* c3*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    } // endfor j

    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r);
     if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  } // endfor i

  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r);
      if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    } // endfor j

  } // endfor i
}


void TimeNSType3_4NLGalerkinDD_2PhaseAxial3D(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double val, val1;
  double *Matrix11Row, *Matrix22Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_U;
  double c0, c3, c4, c5;
  double u1, u2, x, r;

  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y
  Orig2 = OrigValues[2]; // u

  c0 = coeff[0]; // nu
  c3 = coeff[3]; // density ratio (inner/outer)

  u1 = param[0]; // u1old
  u2 = param[1]; // u2old
  x  = param[2]; // x
  
  if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }

  for(i=0;i<N_U;i++)
  {
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val1 = (u1*ansatz10+u2*ansatz01)*test00*c3*r;
      val  = c0*((2.*test10*ansatz10+test01*ansatz01)*r);
      if(TDatabase::ParamDB->Axial3D==1)
     {
      val+=  c0*2.*(ansatz00*test00/r);
     }
      val += val1;
      Matrix11Row[j] += Mult* val;

      val  = c0*(test10*ansatz10+2*test01*ansatz01)*r;
      val += val1;
      Matrix22Row[j] += Mult* val;

    } // endfor j
  } // endfor i
}


// ======================================================================
// Standard Galerkin
// ======================================================================
void GridAssemble4(double Mult, double *coeff,
                  double *param, double hK,
                  double **OrigValues, int *N_BaseFuncts,
                  double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double val, lambda=0;
  double *MatrixRow11, *MatrixRow12, *MatrixRow21, *MatrixRow22;
  double ansatz10, ansatz01, ansatz00;
  double test10, test01, test00;
  double *Orig0, *Orig1, *Orig2;
  int i,j,k,l, N_U;
  double c0=1., c1, c2, detjk;

  double lame = TDatabase::ParamDB->LameC;

  if(TDatabase::ParamDB->P0)
   { detjk = coeff[19]; } // see DiscreteForm2D.C    
  else
   { detjk = 1; }

//   varying stiffness and divergence term values
//   lame *= detjk;


// cout<< "lame " <<lame<<endl;

  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];

  N_U = N_BaseFuncts[0];

  Orig0 = OrigValues[0]; // u_x
  Orig1 = OrigValues[1]; // u_y

  for(i=0;i<N_U;i++)
  {
    MatrixRow11 = MatrixA11[i];
    MatrixRow12 = MatrixA12[i];
    MatrixRow21 = MatrixA21[i];
    MatrixRow22 = MatrixA22[i];

    test10 = Orig0[i];
    test01 = Orig1[i];

    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

     // New part D(v):D(w)
  
      val = (2*test10*ansatz10+test01*ansatz01) + lame*test10*ansatz10;
      MatrixRow11[j] += Mult * val/detjk;

      val  = test01*ansatz10 + lame*test10*ansatz01;
      MatrixRow12[j] += Mult * val/detjk;

      val = test10*ansatz01 + lame*test01*ansatz10;
      MatrixRow21[j] += Mult * val/detjk;

      val = (test10*ansatz10+2*test01*ansatz01) + lame*test01*ansatz01;
      MatrixRow22[j] += Mult * val/detjk;

   
   // old part gradv:gradv

//       val  = c0*(test10*ansatz10+test01*ansatz01 + lambda*test10*ansatz10);
//       MatrixRow11[j] += Mult * val/detjk;
//       val  = c0*(lambda*test10*ansatz01);
//       MatrixRow12[j] += Mult * val/detjk;
//       val  = c0*(lambda*test01*ansatz10);
//       MatrixRow21[j] += Mult * val/detjk;
//       val  = c0*(test10*ansatz10+test01*ansatz01 + lambda*test01*ansatz01);
//       MatrixRow22[j] += Mult * val/detjk;

    } // endfor j
  } // endfor i
 }


// ========================================================================
// parameter routines
// ========================================================================

// ========================================================================
// parameters: u1old, u2old
// used for : GALERKIN
//            COLETTI (non linear steps)
// ========================================================================
void TimeNSParams2(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
}


// comment: routine can be replaced by  TimeNSParamsVelo_GradVelo !!!
void TimeNSParams2RB(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2
}


// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo_ALE(double *in, double *out)
{
  out[0] = in[2] - in[8];       // u1old - w1old
  out[1] = in[3] - in[9];       // u2old - w2old

  out[2] = in[4];               // D1u1
  out[3] = in[5];               // D1u2
  out[4] = in[6];               // D2u1
  out[5] = in[7];               // D2u2
  
//   cout <<in[8] << " E " << in[9] << endl;
}

// ========================================================================
// parameters: u1old, u2old,
// used for : COLETTI, Smagorinsky
// ========================================================================
void TimeNSParamsVelo_GradVelo_ConvVelo(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
}


// ========================================================================
// parameters: u1old, u2old, Frobenius Norm(grad(u))
// used for : COLETTI, first iteration step, viscosity type = 4
// ========================================================================
void TimeNSParamsVelo_GradVeloNuT4(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // convolution of the solution
  out[6] = in[8];                // g_\delta \ast u1
  out[7] = in[9];                // g_\delta \ast u2
}


// ========================================================================
// used for : GL00Convolution, first assembling
// ========================================================================
void TimeNSParamsGL00Convolution(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  // grad u
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // components of convoluted tensor
  out[6] = in[ 8];               // g_d\ast(D1u1*D1u1+D1u2*D1u2)
  out[7] = in[ 9];               // g_d\ast(D1u1*D1u2+D2u1*D2u2)
  out[8] = in[10];               // g_d\ast(D2u1*D2u1+D2u2*D2u2)

}


// ========================================================================
// used for : GL00Convolution, rhs assembling, viscosity type = 4
// ========================================================================
void TimeNSParamsRHSGL00ConvolutionNuT4(double *in, double *out)
{
  // grad u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  // components of convoluted tensor
  out[4] = in[ 8];               // g_d\ast(D1u1*D1u1+D1u2*D1u2)
  out[5] = in[ 9];               // g_d\ast(D1u1*D1u2+D2u1*D2u2)
  out[6] = in[10];               // g_d\ast(D2u1*D2u1+D2u2*D2u2)

  // velocity
  out[7] = in[2];                // u1
  out[8] = in[3];                // u2

  // convoluted velocity
  out[9] = in[11];               //  g_d\ast u1
  out[10] = in[12];              //  g_d\ast u2

}


// ========================================================================
// used for assembling of term coming from the LES model
// GL00AuxProb and GL00Conv
// parameters used in TimeNSRHSLESModel
// ========================================================================
void TimeNSParamsRHSLES(double *in, double *out)
{
  //  solution of auxiliary problem
  out[0] = in[2];
  out[1] = in[3];
  out[2] = in[4];
}


// ========================================================================
// parameters:
// used for : GL00AuxProblem, viscosity type = 4
// ========================================================================
void TimeNSParamsGL00AuxProblemNuT4(double *in, double *out)
{

  // \nabla u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8];               // sol11
  out[5] = in[ 9];               // sol12 = sol21
  out[6] = in[10];               // sol22

  // solution
  out[7] = in[2];                // u1
  out[8] = in[3];                // u2

  // convolution of the solution
  out[9] = in[11];               // g_\delta \ast u1
  out[10] = in[12];              // g_\delta \ast u2
}


// ========================================================================
// parameters:
// used for : GL00AuxProblem
// ========================================================================
void TimeNSParamsGL00AuxProblemPaper2(double *in, double *out)
{
  // \nabla u
  out[0] = in[4];                // D1u1
  out[1] = in[5];                // D1u2
  out[2] = in[6];                // D2u1
  out[3] = in[7];                // D2u2

  //  solution of auxiliary problem
  out[4] = in[ 8];
  out[5] = in[ 9];
  out[6] = in[10];

  // les term with exact solution
  out[7] = in[11];               // (t11)_x
  out[8] = in[12];               // (t12)_x
  out[9] = in[13];               // (t12)_y
  out[10] = in[14];              // (t22)_y

}


// ========================================================================
// parameters: gradient(u1), gradient(u2)
// ========================================================================
void TimeNSParamsGrad(double *in, double *out)
{
  // cout << "GRAD" << endl;
  out[0] = in[4];                // D10(u1old)
  out[1] = in[5];                // D10(u2old)
  out[2] = in[6];                // D01(u1old)
  out[3] = in[7];                // D01(u2old)
  // cout << in[4] << " E " << in[5] << " " << in[6] << " " << in[7] << endl;
}


// ========================================================================
// used for VMS, assembling of rhs for small scale equation
// ========================================================================
void TimeNSParams_VMS_SmallRhs2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // small scales
  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
  out[8] = in[10];               // D1u1
  out[9] = in[11];               // D1u2
  out[10] = in[12];              // D2u1
  out[11] = in[13];              // D2u2

  // large pressure
  out[12] = in[14];              // p
}


// ========================================================================
// used for VMS, assembling of rhs for large scale equation
// ========================================================================
void TimeNSParams_VMS_LargeRhs2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  // small scales
  out[6] = in[8];                // u1old
  out[7] = in[9];                // u2old
  out[8] = in[10];               // D1u1
  out[9] = in[11];               // D1u2
  out[10] = in[12];              // D2u1
  out[11] = in[13];              // D2u2

  // small pressure
  out[12] = in[14];              // p
}


// ========================================================================
// used for VMS, assembling of matrix, variant 1
// ========================================================================
void TimeNSParams_NLGalerkin_VMS_1_2D(double *in, double *out)
{
  // large scales
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  // small scales
  out[2] = in[4];                // u1current
  out[3] = in[5];                // u2current
}
// ========================================================================
// parameters: u1old, u2old, gradients; values of previous time steps
// defect correction type 2
// ========================================================================
void TimeNSParamsVelo_GradVeloOld2(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old

  out[2] = in[4];                // D1u1
  out[3] = in[5];                // D1u2
  out[4] = in[6];                // D2u1
  out[5] = in[7];                // D2u2

  out[6] = in[8];     //    u_1m11 
  out[7] = in[9];  // u_1m12
  out[8] = in[10]; // u_1m21
  out[9] = in[11]; // u_1m22
}

// ========================================================================
// parameters: x, y, u1old, u2old
// used for : SSMUM
// ========================================================================
void TimeNSParamsVeloPos(double *in, double *out)
{
  out[0] = in[0];                // x
  out[1] = in[1];                // y
  out[2] = in[2];                // u1old
  out[3] = in[3];                // u2old
}

// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams(double *in, double *out)
{
  out[0] = in[2] - in[4];
  out[1] = in[3] - in[5];

// cout << " out[0] " << in[4] << " out[1] " << in[4] <<endl;
}
// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams_Axial3D(double *in, double *out)
{
  out[0] = in[2] - in[4];
  out[1] = in[3] - in[5];
  out[2] = in[0];     // x value for axial symmetric

// cout<< "out[2]  " << out[2]<<endl;
}

// ===============================================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, grad(tau1, tau2, tau3)
// ===============================================================================================
void TimeNSParamsVelo_GradVelo_CST(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9]  = in[11];              // D10 (tau1old)
  out[10] = in[12];              // D10 (tau2old)
  out[11] = in[13];              // D10 (tau3old)
  out[12] = in[14];              // D01 (tau1old)
  out[13] = in[15];             // D01 (tau2old)
  out[14] = in[16];             // D01 (tau3old)
}

// ========================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), grad (tau1, tau2, tau3), D1, D2, D3
// ========================================================================
void TimeNSParamsVelo_GradVelo_CST_DEVSS(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9] = in[11];               // D10(tau1old)
  out[10] = in[12];              // D10(tau2old)
  out[11] = in[13];              // D10(tau3old)
  out[12] = in[14];              // D01(tau1old)
  out[13] = in[15];              // D01(tau2old)
  out[14] = in[16];              // D01(tau3old)
  out[15] = in[17];              // D00 (D1old)
  out[16] = in[18];              // D00 (D2old)
  out[17] = in[19];              // D00 (D3old)

 }

// ======================================================================================
// parameters: tau1old, tau2old, tau3old, gradient(tau1), gradient(tau2), gradient(tau3)
// ======================================================================================
void TimeNSParams_CST(double *in, double *out)
{
  out[0] = in[2];                // tau1old
  out[1] = in[3];                // tau2old
  out[2] = in[4];                // tau3old
  out[3] = in[5];                // tau1old D10
  out[4] = in[6];                // tau2old D10
  out[5] = in[7];                // tau3old D10
  out[6] = in[8];                // tau1old D01
  out[7] = in[9];                // tau2old D01
  out[8] = in[10];               // tau3old D01
}




// ===============================================================================================
// parameters: u1old, u2old, gradient(u1), gradient(u2), tau1, tau2, tau3, grad(tau1, tau2, tau3), u1-w1, u2-w2, r 
// ===============================================================================================
void MovingNSParamsVelo_GradVelo_CST_Axial3D(double *in, double *out)
{
  out[0] = in[2] ;                // u1old 
  out[1] = in[3] ;                // u2old 
  out[2] = in[4];                // D10(u1old)
  out[3] = in[5];                // D10(u2old)
  out[4] = in[6];                // D01(u1old)
  out[5] = in[7];                // D01(u2old)
  out[6] = in[8];                // D00(tau1old)
  out[7] = in[9];                // D00(tau2old)
  out[8] = in[10];               // D00(tau3old)
  out[9]  = in[11];              // D10 (tau1old)
  out[10] = in[12];              // D10 (tau2old)
  out[11] = in[13];              // D10 (tau3old)
  out[12] = in[14];              // D01 (tau1old)
  out[13] = in[15];             // D01 (tau2old)
  out[14] = in[16];             // D01 (tau3old)
  out[15] = in[17] - in[19];     // u1old - w1
  out[16] = in[18] - in[20];     // u2old - w2
  out[17] = in[0];             // x value for axial symmetric
}


// ======================================================================
// parameters u1old, u2old, gridv_x, gridv_y
// ======================================================================
void MovingTNSParams_Axial3D_HeatLine(double *in, double *out)
{
  out[0] = in[2]; 
  out[1] = in[3];
  out[2] = in[4]; 
  out[3] = in[5] - in[7];
  out[4] = in[6] - in[8]; 
  out[5] = in[9] - in[11]; 
  out[6] = in[10] - in[12]; 
  out[7] = in[0];  // x value for axial symmetric

// cout<< "out[0]  " << out[0]<<endl;
}



void TimeNSType2SUPG(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatA, **MatM, **MatK;
  double **MatB1, **MatB2, **MatB1T, **MatB2T;
  double *Rhs1, *Rhs2;
  double *MatRowA, *MatRowM, *MatRowK, *MatRow1, *MatRow2;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  double test00, test10, test01, test20, test02;
  double ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double c0, c1, c2, u1, u2, val;
  int N_U, N_P;

  MatA=LocMatrices[0];
  MatM=LocMatrices[1];
  MatK=LocMatrices[2];
  
  MatB1=LocMatrices[3];
  MatB2=LocMatrices[4];
  MatB1T=LocMatrices[5];
  MatB2T=LocMatrices[6];

  Rhs1=LocRhs[0];
  Rhs2=LocRhs[1];

  N_U=N_BaseFuncts[0];
  N_P=N_BaseFuncts[1];
  

  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
  Orig2=OrigValues[2]; // u 
  Orig3=OrigValues[3]; // u_xx
  Orig4=OrigValues[4]; // u_yy
  Orig5=OrigValues[5]; // p_x 
  Orig6=OrigValues[6]; // p_y
  Orig7=OrigValues[7]; // p

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];

  // stabilization parameters have to defined
  // for initial test its setted to delta0*hK*hK
  double delta=TDatabase::ParamDB->DELTA0*hK*hK;
  double ugrad;

  for(int i=0;i<N_U;i++)
  {
    MatRowA=MatA[i];
    MatRowM=MatM[i];
    MatRowK=MatK[i];
    
    test10=Orig0[i];
    test01=Orig1[i];
    test00=Orig2[i];

    ugrad= delta*(u1*test10 + u2*test01);

    Rhs1[i] += Mult*(test00 + ugrad)*c1;
    Rhs2[i] += Mult*(test00 + ugrad)*c2;

    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];
      ansatz00=Orig2[j];
      ansatz20=Orig3[j];
      ansatz02=Orig4[j];

      // Galerkin terms
      val =c0*(test10*ansatz10+test01*ansatz01); // diffusion term
      val +=(u1*ansatz10 + u2*ansatz01)*test00; // nonlinear (convective term)
      val +=(-c0*(ansatz20+ansatz02) + (u1*ansatz10 + u2*ansatz01))*ugrad; // SUPG terms
      
      MatRowA[j] += Mult*val; // A block
      MatRowM[j] += Mult*ansatz00*test00;

      MatRowK[j] += Mult*ansatz00*ugrad;
    }

    MatRow1=MatB1T[i];
    MatRow2=MatB2T[i];
    
    for(int j=0;j<N_P;j++)
    {
      ansatz10=Orig5[j];
      ansatz01=Orig6[j];
      ansatz00=Orig7[j];

      val = -ansatz00*test10; // Galerkin term
      val += ansatz10*ugrad;
      MatRow1[j] += Mult*val;
      
      val = -ansatz00*test01; // term due to stabilization
      val += ansatz01*ugrad;
      MatRow2[j] += Mult*val;
    }
  }

  for(int i=0;i<N_P; i++)
  {
    MatRow1=MatB1[i];
    MatRow2=MatB2[i];

    test00=Orig7[i];
    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];

      val = -Mult*test00*ansatz10;
      MatRow1[j] += val;
      
      val = -Mult*test00*ansatz01;
      MatRow2[j] += val;
    }
  }
}


void TimeNSType2NLSUPG(double Mult, double *coeff, 
                double *param, double hK, 
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double **MatA, **MatK, **MatB1T, **MatB2T;
  double *MatRowA, *MatRowK, *MatRow1, *MatRow2;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  double test00, test10, test01, ansatz00, ansatz10, ansatz01, ansatz20, ansatz02;
  double *Rhs1, *Rhs2;
  double c0, c1, c2, u1, u2, val, delta, ugrad;
  int N_U, N_P;

  MatA=LocMatrices[0];
  MatK=LocMatrices[1];
  MatB1T=LocMatrices[2];
  MatB2T=LocMatrices[3];


  Rhs1=LocRhs[0];
  Rhs2=LocRhs[1];

  N_U=N_BaseFuncts[0];
  N_P=N_BaseFuncts[1];
  

  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
  Orig2=OrigValues[2]; // u 
  Orig3=OrigValues[3]; // u_xx
  Orig4=OrigValues[4]; // u_yy
  Orig5=OrigValues[5]; // p_x 
  Orig6=OrigValues[6]; // p_y
  Orig7=OrigValues[7]; // p

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];
  
  if(c0<hK)
    delta=TDatabase::ParamDB->DELTA0*hK*hK;
  else 
    delta=TDatabase::ParamDB->DELTA1*hK*hK;

  for(int i=0;i<N_U;i++)
  {
    MatRowA=MatA[i];
    MatRowK=MatK[i];
    
    test10=Orig0[i]; 
    test01=Orig1[i];
    test00=Orig2[i];
    
    ugrad=delta*(u1*test10 + u2*test01);
    
    Rhs1[i] += Mult*ugrad*c1;
    Rhs2[i] += Mult*ugrad*c2;
    
    for(int j=0;j<N_U;j++)
    {
      ansatz10=Orig0[j];
      ansatz01=Orig1[j];
      ansatz00=Orig2[j];
      ansatz20=Orig3[j];
      ansatz02=Orig4[j];
      
      val =c0*(test10*ansatz10 + test01*ansatz01); // diffusion term
      val +=(u1*ansatz10+u2*ansatz01)*test00; // nonlinear term
      val +=(/*-c0*(ansatz20 + ansatz02) +*/ (u1*ansatz10 + u2*ansatz01))*ugrad; // SUPG terms
      
      MatRowA[j] += Mult*val;
      
      MatRowK[j] += Mult*ansatz00*ugrad;
    }
    
    MatRow1=MatB1T[i];
    MatRow2=MatB2T[i];
    for(int j=0;j<N_P;j++)
    {
      ansatz10=Orig5[j];
      ansatz01=Orig6[j];
      ansatz00=Orig7[j];
      
      val = -ansatz00*test10;
      val += ansatz10*ugrad;
      MatRow1[j] += Mult*val;
      
      val = -ansatz00*test01;
      val += ansatz01*ugrad;
      MatRow2[j] += Mult*val;
    }
  }
}

void TimeNSRHSSUPG(double Mult, double *coeff,
                double *param, double hK,
                double **OrigValues, int *N_BaseFuncts,
                double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4;
  double *Orig0, *Orig1, *Orig2;
  double test00, test10, test01;
  double u1, u2, c0, c1, c2, ugrad, delta;
  int N_U;
  
//   Rhs1=LocRhs[0];
//   Rhs2=LocRhs[1];
  Rhs3=LocRhs[0];
  Rhs4=LocRhs[1];

  N_U = N_BaseFuncts[0];
  
  Orig0=OrigValues[0]; // u_x
  Orig1=OrigValues[1]; // u_y 
  Orig2=OrigValues[2]; // u 

  c0=coeff[0];
  c1=coeff[1];
  c2=coeff[2];

  u1=param[0];
  u2=param[1];
  
  if(c0<hK)
    delta=TDatabase::ParamDB->DELTA0*hK*hK;
  else 
    delta=TDatabase::ParamDB->DELTA1*hK*hK;

  for(int i=0;i<N_U;i++)
  {
    test10=Orig0[i]; 
    test01=Orig1[i];
    test00=Orig2[i];
    
    ugrad=delta*(u1*test10 + u2*test01);
    /*
    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;*/
    
    Rhs3[i] += Mult*ugrad*c1;
    Rhs4[i] += Mult*ugrad*c2;
  }   
}

void TimeNSParams4(double *in, double *out)
{
  out[0] = in[2];                // u1old
  out[1] = in[3];                // u2old
  out[2] = in[4];                // u1, previous time
  out[3] = in[5];                // u2, previous time
}





// D(u):D(v) TNSECST Galerkin
void Time_NSEType4_Galerkin_CST_Galerkin(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixMU11, **MatrixMU22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixMS11, **MatrixMS22, **MatrixMS33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, *Rhs7, *Rhs8, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6;
  int i,j,N_U, N_P, N_S;
  double c0, c1, c2, c3, c4, c5, c6, c7;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixMU11 = LocMatrices[4];
  MatrixMU22 = LocMatrices[5];
  
  MatrixG11 = LocMatrices[6];
  MatrixG12 = LocMatrices[7];
  MatrixG21 = LocMatrices[8];
  MatrixG22 = LocMatrices[9];
  MatrixG23 = LocMatrices[10];
  MatrixG32 = LocMatrices[11];
  MatrixG33 = LocMatrices[12];
  
  MatrixMS11 = LocMatrices[13];
  MatrixMS22 = LocMatrices[14];
  MatrixMS33 = LocMatrices[15];
  
  MatrixB1 = LocMatrices[16];
  MatrixB2 = LocMatrices[17];
  
  MatrixB1T = LocMatrices[18];
  MatrixB2T = LocMatrices[19];
  
  MatrixC11 = LocMatrices[20];
  MatrixC12 = LocMatrices[21];
  MatrixC22 = LocMatrices[22];
  MatrixC23 = LocMatrices[23];
  
  MatrixD11 = LocMatrices[24];
  MatrixD12 = LocMatrices[25];
  MatrixD21 = LocMatrices[26];
  MatrixD22 = LocMatrices[27];
  MatrixD31 = LocMatrices[28];
  MatrixD32 = LocMatrices[29];
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];
  Rhs7 = LocRhs[6];
  Rhs8 = LocRhs[7];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixMU11[i];
    MatrixM22Row  = MatrixMU22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*(2*test10*ansatz10+ test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      { 
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;   
      
      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      
    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      val  = c7*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j

    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  

  
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];

    val = c4*test00;
    Rhs3[i] += Mult*val;
    val = 2*c5*test00;
    Rhs4[i] += Mult*val;
    val = c6 *test00;
    Rhs5[i] += Mult*val;
    
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00;
    Rhs6[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00;
    Rhs7[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00;
    Rhs8[i] += Mult*val;
    
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix32Row[j] += Mult * val;
      

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
    MatrixM11Row  = MatrixMS11[i];
    MatrixM22Row  = MatrixMS22[i];
    MatrixM33Row  = MatrixMS33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*test00  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00;
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01) - 2.0*(u1x + u2y)*ansatz00)*test00 ;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00;
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*test00;
      Matrix33Row[j] += Mult * val;
      
      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;
      

    }                            // endfor j

   
  }                              // endfor i

  

  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
  
}


// D(u):D(v) TNSECST Galerkin for two phase Axial3D flows 
void Time_NSEType4_Galerkin_CST_Galerkin_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixMU11, **MatrixMU22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixMS11, **MatrixMS22, **MatrixMS33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, *Rhs7, *Rhs8, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6;
  int i,j,N_U, N_P, N_S;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;
  double u1,u2, u1w1, u2w2, u1x, u1y, u2x, u2y, x, r;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixMU11 = LocMatrices[4];
  MatrixMU22 = LocMatrices[5];
  
  MatrixG11 = LocMatrices[6];
  MatrixG12 = LocMatrices[7];
  MatrixG21 = LocMatrices[8];
  MatrixG22 = LocMatrices[9];
  MatrixG23 = LocMatrices[10];
  MatrixG32 = LocMatrices[11];
  MatrixG33 = LocMatrices[12];
  
  MatrixMS11 = LocMatrices[13];
  MatrixMS22 = LocMatrices[14];
  MatrixMS33 = LocMatrices[15];
  
  MatrixB1 = LocMatrices[16];
  MatrixB2 = LocMatrices[17];
  
  MatrixB1T = LocMatrices[18];
  MatrixB2T = LocMatrices[19];
  
  MatrixC11 = LocMatrices[20];
  MatrixC12 = LocMatrices[21];
  MatrixC22 = LocMatrices[22];
  MatrixC23 = LocMatrices[23];
  
  MatrixD11 = LocMatrices[24];
  MatrixD12 = LocMatrices[25];
  MatrixD21 = LocMatrices[26];
  MatrixD22 = LocMatrices[27];
  MatrixD31 = LocMatrices[28];
  MatrixD32 = LocMatrices[29];
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];
  Rhs7 = LocRhs[6];
  Rhs8 = LocRhs[7];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
   c8 = coeff[8]; // density ratio (inner/outer)
   c9 = coeff[9]; // Fluid type, 1-Newtonian, 2-Oldroyd-B, 3-Giesekus
   c10 = coeff[10];// Tau dot Tau in Giesekus
   
  u1 = param[0];                 // u1old 
  u2 = param[1];                 // u2old 
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold
  u1w1 =  param[15];             // u1old - grid VX   !!check aux parameters
  u2w2 =  param[16];             // u2old - grid VY   !!check aux parameters

  x  = param[17]; // x
  
 if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1*c8*r;
    Rhs2[i] += Mult*test00*c2*c8*r;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixMU11[i];
    MatrixM22Row  = MatrixMU22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*((2*test10*ansatz10+ test01*ansatz01)*r);
        if(TDatabase::ParamDB->Axial3D==1)
      {
       val += c0*2.*ansatz00*test00/r;
      }
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*c8*r;

      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01)*r;
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*c8*r;

      Matrix22Row[j] += Mult * val;   
      
      val = Mult*(ansatz00*test00)*c8*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      
    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      if(c9 == 2 || c9==3)
      {
      val  = c7*test10*ansatz00*r;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00*r;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix23Row[j] += Mult * val;
      }
      else if(c9 == 1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
      }
      else
      {
	cout<<"Invalid fluid type \n";
	exit(4711);
      }
      
    }                            // endfor j

    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r);
       if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];

    val = c4*test00*r;
    Rhs3[i] += Mult*val;
    val = 2*c5*test00*r;
    Rhs4[i] += Mult*val;
    val = c6 *test00*r;
    Rhs5[i] += Mult*val;
    
    if(c9 == 2 || c9==3)
    {
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00*r;
    Rhs6[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00*r;
    Rhs7[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00*r;
    Rhs8[i] += Mult*val;
    }
    else if(c9 == 1)
    {
      val = 0;
      Rhs6[i] += Mult*val;
      Rhs7[i] += Mult*val;
      Rhs8[i] += Mult*val;
    }
     else
    {
	cout<<"Invalid fluid type \n";
	exit(4711);
    }
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      if(c9==2 || c9==3)
      {
      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00*r;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00*r;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix32Row[j] += Mult * val;
      }
      else if(c9==1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix31Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
    {
	cout<<"Invalid fluid type \n";
	exit(4711);
    }
      
    }        // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
    MatrixM11Row  = MatrixMS11[i];
    MatrixM22Row  = MatrixMS22[i];
    MatrixM33Row  = MatrixMS33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      if(c9 == 2 || c9==3)
    {
      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1w1*ansatz10 + u2w2*ansatz01)*test00*r  ;
      if(c9==3)
     {
      val += c10*tau1*ansatz00*test00*r;
     }
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
     Matrix12Row[j] += Mult * val;
     
      val = -2.0*ansatz00*u2x*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
     Matrix21Row[j] += Mult * val;
     
      val = (c3*2.0*ansatz00 + 2.0*(u1w1*ansatz10 + u2w2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00*r ;
      if(c9==3)
     {
      val += c10*(tau1 + tau3)*ansatz00*test00*r;
     }
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
       if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1w1*ansatz10 + u2w2*ansatz01)*test00*r;
      if(c9==3)
     {
      val += c10*tau3*ansatz00*test00*r;
     }
      Matrix33Row[j] += Mult * val;
    }
      
      else if(c9 == 1)
      {
	val  = c3*ansatz00*test00*r;
	Matrix11Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val*2.0;
	Matrix33Row[j] += Mult * val;
	val = 0;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
      
      if(c9 == 2 || c9==3)
    {
      val = Mult*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;
    }
       else if(c9 == 1)
      {
      val = 0;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;
      }
  

    }                            // endfor j

   
  }                              // endfor i

  

  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r);
        if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
  
}



// D(u):D(v) TNSECST Galerkin for impinging droplet Axial3D flows 
void Time_NSEType4_Galerkin_CST_Galerkin_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixMU11, **MatrixMU22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixMS11, **MatrixMS22, **MatrixMS33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, *Rhs7, *Rhs8, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6;
  int i,j,N_U, N_P, N_S;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9;
  double u1,u2, u1w1, u2w2, u1x, u1y, u2x, u2y, x, r;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixMU11 = LocMatrices[4];
  MatrixMU22 = LocMatrices[5];
  
  MatrixG11 = LocMatrices[6];
  MatrixG12 = LocMatrices[7];
  MatrixG21 = LocMatrices[8];
  MatrixG22 = LocMatrices[9];
  MatrixG23 = LocMatrices[10];
  MatrixG32 = LocMatrices[11];
  MatrixG33 = LocMatrices[12];
  
  MatrixMS11 = LocMatrices[13];
  MatrixMS22 = LocMatrices[14];
  MatrixMS33 = LocMatrices[15];
  
  MatrixB1 = LocMatrices[16];
  MatrixB2 = LocMatrices[17];
  
  MatrixB1T = LocMatrices[18];
  MatrixB2T = LocMatrices[19];
  
  MatrixC11 = LocMatrices[20];
  MatrixC12 = LocMatrices[21];
  MatrixC22 = LocMatrices[22];
  MatrixC23 = LocMatrices[23];
  
  MatrixD11 = LocMatrices[24];
  MatrixD12 = LocMatrices[25];
  MatrixD21 = LocMatrices[26];
  MatrixD22 = LocMatrices[27];
  MatrixD31 = LocMatrices[28];
  MatrixD32 = LocMatrices[29];
 
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];
  Rhs7 = LocRhs[6];
  Rhs8 = LocRhs[7];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
   c8 = coeff[8]; // Fluid type, 1-Newtonian, 2-Oldroyd-B, 3-Giesekus
   c9 = coeff[9];// Tau dot Tau in Giesekus
   
  u1 = param[0];                 // u1old 
  u2 = param[1];                 // u2old 
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold
  u1w1 =  param[15];             // u1old - grid VX   !!check aux parameters
  u2w2 =  param[16];             // u2old - grid VY   !!check aux parameters

  x  = param[17]; // x
  
 if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1*r;
    Rhs2[i] += Mult*test00*c2*r;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixMU11[i];
    MatrixM22Row  = MatrixMU22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*((2*test10*ansatz10+ test01*ansatz01)*r);
        if(TDatabase::ParamDB->Axial3D==1)
      {
       val += c0*2.*ansatz00*test00/r;
      }
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*r;

      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10)*r;
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01)*r;
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01)*r;
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*r;

      Matrix22Row[j] += Mult * val;   
      
      val = Mult*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

      
    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      if(c8 == 2 || c8==3)
      {
      val  = c7*test10*ansatz00*r;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00*r;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00*r;
      Matrix23Row[j] += Mult * val;
      }
      else if(c8 == 1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
      }
      else
      {
	cout<<"Invalid fluid type \n";
	exit(4711);
      }
      
    }                            // endfor j

    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*(ansatz00*test10*r);
       if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01*r;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];

    val = c4*test00*r;
    Rhs3[i] += Mult*val;
    val = 2*c5*test00*r;
    Rhs4[i] += Mult*val;
    val = c6 *test00*r;
    Rhs5[i] += Mult*val;
    
    if(c8 == 2 || c8==3)
    {
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00*r;
    Rhs6[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00*r;
    Rhs7[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00*r;
    Rhs8[i] += Mult*val;
    }
    else if(c8 == 1)
    {
      val = 0;
      Rhs6[i] += Mult*val;
      Rhs7[i] += Mult*val;
      Rhs8[i] += Mult*val;
    }
     else
    {
	cout<<"Invalid fluid type \n";
	exit(4711);
    }
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      if(c8==2 || c8==3)
      {
      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00*r;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00*r;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix32Row[j] += Mult * val;
      }
      else if(c8==1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix31Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
    {
	cout<<"Invalid fluid type \n";
	exit(4711);
    }
      
    }        // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
    MatrixM11Row  = MatrixMS11[i];
    MatrixM22Row  = MatrixMS22[i];
    MatrixM33Row  = MatrixMS33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      if(c8 == 2 || c8==3)
    {
      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1w1*ansatz10 + u2w2*ansatz01)*test00*r  ;
      if(c8==3)
     {
      val += c9*tau1*ansatz00*test00*r;
     }
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
     Matrix12Row[j] += Mult * val;
     
      val = -2.0*ansatz00*u2x*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
     Matrix21Row[j] += Mult * val;
     
      val = (c3*2.0*ansatz00 + 2.0*(u1w1*ansatz10 + u2w2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00*r ;
      if(c8==3)
     {
      val += c9*(tau1 + tau3)*ansatz00*test00*r;
     }
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
       if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1w1*ansatz10 + u2w2*ansatz01)*test00*r;
      if(c8==3)
     {
      val += c9*tau3*ansatz00*test00*r;
     }
      Matrix33Row[j] += Mult * val;
    }
      
      else if(c8 == 1)
      {
	val  = c3*ansatz00*test00*r;
	Matrix11Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val*2.0;
	Matrix33Row[j] += Mult * val;
	val = 0;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
      
      if(c8 == 2 || c8==3)
    {
      val = Mult*(ansatz00*test00)*r;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;
    }
       else if(c8 == 1)
      {
      val = 0;
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;
      }
  

    }                            // endfor j

   
  }                              // endfor i

  

  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val = -Mult*(test00*ansatz10*r);
        if(TDatabase::ParamDB->Axial3D==1)
      {
       val += -Mult*(ansatz00*test00);
      }
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01*r;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
  
}



// D(u):D(v) TNSECST DEVSS/SUPG
void Time_NSEType4_Galerkin_CST_SUPG_DEVSS(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA12, **MatrixA21, **MatrixA22;
  double **MatrixMU11, **MatrixMU22;
  double **MatrixB1, **MatrixB2;
  double **MatrixB1T, **MatrixB2T;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixMS11, **MatrixMS22, **MatrixMS33;
  double **MatrixC11, **MatrixC12, **MatrixC22, **MatrixC23;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double **MatrixH11, **MatrixH22, **MatrixH33;
  double **MatrixE11, **MatrixE12, **MatrixE22, **MatrixE23;
  double **MatrixJ11, **MatrixJ21, **MatrixJ22, **MatrixJ32;
  
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5, *Rhs6, *Rhs7, *Rhs8, val;
  
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double *MatrixM11Row, *MatrixM22Row, *MatrixM33Row;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5, *Orig6, *Orig7;
  int i,j,N_U, N_P, N_S, N_D;
  double c0, c1, c2, c3, c4, c5, c6, c7, c8;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y, delta, norm_u, ugrad;
  double delta0 = TDatabase::ParamDB->DELTA0;
    
  MatrixA11 = LocMatrices[0];
  MatrixA12 = LocMatrices[1];
  MatrixA21 = LocMatrices[2];
  MatrixA22 = LocMatrices[3];
  
  MatrixMU11 = LocMatrices[4];
  MatrixMU22 = LocMatrices[5];
  
  MatrixG11 = LocMatrices[6];
  MatrixG12 = LocMatrices[7];
  MatrixG21 = LocMatrices[8];
  MatrixG22 = LocMatrices[9];
  MatrixG23 = LocMatrices[10];
  MatrixG32 = LocMatrices[11];
  MatrixG33 = LocMatrices[12];
  
  MatrixMS11 = LocMatrices[13];
  MatrixMS22 = LocMatrices[14];
  MatrixMS33 = LocMatrices[15];
  
  MatrixH11 = LocMatrices[16];
  MatrixH22 = LocMatrices[17];
  MatrixH33 = LocMatrices[18];
  
  MatrixB1 = LocMatrices[19];
  MatrixB2 = LocMatrices[20];
  
  MatrixB1T = LocMatrices[21];
  MatrixB2T = LocMatrices[22];
  
  MatrixC11 = LocMatrices[23];
  MatrixC12 = LocMatrices[24];
  MatrixC22 = LocMatrices[25];
  MatrixC23 = LocMatrices[26];
  
  MatrixD11 = LocMatrices[27];
  MatrixD12 = LocMatrices[28];
  MatrixD21 = LocMatrices[29];
  MatrixD22 = LocMatrices[30];
  MatrixD31 = LocMatrices[31];
  MatrixD32 = LocMatrices[32];
 
  MatrixE11 = LocMatrices[33];
  MatrixE12 = LocMatrices[34];
  MatrixE22 = LocMatrices[35];
  MatrixE23 = LocMatrices[36];
  
  MatrixJ11 = LocMatrices[37];
  MatrixJ21 = LocMatrices[38];
  MatrixJ22 = LocMatrices[39];
  MatrixJ32 = LocMatrices[40];


  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];
  Rhs6 = LocRhs[5];
  Rhs7 = LocRhs[6];
  Rhs8 = LocRhs[7];

  N_U = N_BaseFuncts[0];
  N_P = N_BaseFuncts[1];
  N_S = N_BaseFuncts[2];
  N_D = N_BaseFuncts[3];

  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // p
  Orig4 = OrigValues[4];         // tau_x
  Orig5 = OrigValues[5];         // tau_y
  Orig6 = OrigValues[6];         // tau
  Orig7 = OrigValues[7];         // deformation

   c0 = coeff[0]; // D(u):D(v) term 
   c1 = coeff[1]; // f1
   c2 = coeff[2]; // f2
   c3 = coeff[3]; // Tau term in CST 
   c4 = coeff[4]; // f3
   c5 = coeff[5]; // f4 
   c6 = coeff[6]; // f5
   c7 = coeff[7]; // div Tau term 
   c8 = coeff[8]; // div Deformation term 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
    
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix12Row = MatrixA12[i];
    Matrix21Row = MatrixA21[i];
    Matrix22Row = MatrixA22[i];
    MatrixM11Row  = MatrixMU11[i];
    MatrixM22Row  = MatrixMU22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*(2*test10*ansatz10+ test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;

      val  = c0*(test01*ansatz10);
      Matrix12Row[j] += Mult * val;

      val  = c0*(test10*ansatz01);
      Matrix21Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      { 
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;   
      
      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += val;

    }                            // endfor j
    
    // Matrix C
    Matrix11Row = MatrixC11[i];
    Matrix12Row = MatrixC12[i];
    Matrix22Row = MatrixC22[i];
    Matrix23Row = MatrixC23[i];
     for(j=0;j<N_S;j++)
    {
      ansatz00 = Orig6[j];

      val  = c7*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c7*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c7*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j

    // Matrix E
    Matrix11Row = MatrixE11[i];
    Matrix12Row = MatrixE12[i];
    Matrix22Row = MatrixE22[i];
    Matrix23Row = MatrixE23[i];
     for(j=0;j<N_D;j++)
    {
      ansatz00 = Orig7[j];

      val  = c8*test10*ansatz00;
      Matrix11Row[j] += Mult * val;

      val  = c8*test01*ansatz00;
      Matrix12Row[j] += Mult * val;

      val  = c8*test10*ansatz00;
      Matrix22Row[j] += Mult * val;

      val  = c8*test01*ansatz00;
      Matrix23Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix BT
    MatrixRow1 = MatrixB1T[i];
    MatrixRow2 = MatrixB2T[i];
    for(j=0;j<N_P;j++)
    {
      ansatz00 = Orig3[j];

      val = -Mult*ansatz00*test10;
      MatrixRow1[j] += val;

      val = -Mult*ansatz00*test01;
      MatrixRow2[j] += val;
    }
  }                              // endfor i

  norm_u = sqrt((u1*u1) + (u2*u2));
  
  if (norm_u < 0.000001)
    delta=0; 
  else
    delta = delta0*hK/norm_u;
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig4[i];
    test01 = Orig5[i];
    test00 = Orig6[i];
    ugrad  = delta * (u1*test10+u2*test01);

    val = c4*(test00+ugrad);
    Rhs3[i] += Mult*val;
    val = 2*c5*(test00+ugrad);
    Rhs4[i] += Mult*val;
    val = c6*(test00+ugrad);
    Rhs5[i] += Mult*val;
    
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*(test00+ugrad);
    Rhs6[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*(test00+ugrad);
    Rhs7[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*(test00+ugrad);
    Rhs8[i] += Mult*val;
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*(test00+ugrad);
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*(test00+ugrad);
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
    MatrixM11Row  = MatrixMS11[i];
    MatrixM22Row  = MatrixMS22[i];
    MatrixM33Row  = MatrixMS33[i];

     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig4[j];
      ansatz01 = Orig5[j];
      ansatz00 = Orig6[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*(test00+ugrad)  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*(test00+ugrad) ;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*(test00+ugrad);
      Matrix33Row[j] += Mult * val;
      
      val = Mult*(ansatz00*test00);
      MatrixM11Row[j] += val;
      MatrixM22Row[j] += 2*val;
      MatrixM33Row[j] += val;

    }                            // endfor j

   
  }                              // endfor i

  
      for(i=0;i<N_D;i++)
  {

    test00 = Orig7[i];
    
    // Matrix J
    Matrix11Row = MatrixJ11[i];
    Matrix21Row = MatrixJ21[i];
    Matrix22Row = MatrixJ22[i];
    Matrix32Row = MatrixJ32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val  = -ansatz10*test00;
      Matrix11Row[j] += Mult * val;

      val  = -ansatz01*test00;
      Matrix21Row[j] += Mult * val;

      val  = -ansatz10*test00;
      Matrix22Row[j] += Mult * val;
      
      val  = -ansatz01*test00;
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix H
    Matrix11Row = MatrixH11[i];
    Matrix22Row = MatrixH22[i];
    Matrix33Row = MatrixH33[i];
    
     for(j=0;j<N_D;j++)
    {
      ansatz00 = Orig7[j];

      val  = ansatz00*test00;
      Matrix11Row[j] += Mult * val;
      
      val = 2.0*ansatz00*test00;
      Matrix22Row[j] += Mult * val;

      val  = ansatz00*test00;
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

  
  // Matrix B
  for(i=0;i<N_P;i++)
  {
    MatrixRow1 = MatrixB1[i];
    MatrixRow2 = MatrixB2[i];

    test00 = Orig3[i];
  
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];

      val = -Mult*test00*ansatz10;
      MatrixRow1[j] += val;

      val = -Mult*test00*ansatz01;
      MatrixRow2[j] += val;
    }                            // endfor j
  }                              // endfor i
  
  
}


// D(u):D(v) TNSECST DEVSS/SUPG non-linear terms only
void Time_NSEType4_Galerkin_CST_SUPG_DEVSS_NLTerms(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_S;
  double c0, c3;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y, delta, norm_u, ugrad;
  double delta0 = TDatabase::ParamDB->DELTA0;
    
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  
  MatrixG11 = LocMatrices[2];
  MatrixG12 = LocMatrices[3];
  MatrixG21 = LocMatrices[4];
  MatrixG22 = LocMatrices[5];
  MatrixG23 = LocMatrices[6];
  MatrixG32 = LocMatrices[7];
  MatrixG33 = LocMatrices[8];
  
  MatrixD11 = LocMatrices[9];
  MatrixD12 = LocMatrices[10];
  MatrixD21 = LocMatrices[11];
  MatrixD22 = LocMatrices[12];
  MatrixD31 = LocMatrices[13];
  MatrixD32 = LocMatrices[14];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // tau_x
  Orig4 = OrigValues[4];         // tau_y
  Orig5 = OrigValues[5];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c3 = coeff[3]; // Tau term in CST 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
   
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*(2*test10*ansatz10+ test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;


      val  = c0*(test10*ansatz10+ 2*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      { 
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;   
    }                            // endfor j

  }                              // endfor i

  norm_u = sqrt((u1*u1) + (u2*u2));
  
  if (norm_u < 0.000001)
    delta=0; 
  else
    delta = delta0*hK/norm_u;
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig3[i];
    test01 = Orig4[i];
    test00 = Orig5[i];
    ugrad  = delta * (u1*test10+u2*test01);

    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*(test00+ugrad);
    Rhs1[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*(test00+ugrad);
    Rhs2[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*(test00+ugrad);
    Rhs3[i] += Mult*val;
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*(test00+ugrad);
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*(test00+ugrad);
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*(test00+ugrad);
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig3[j];
      ansatz01 = Orig4[j];
      ansatz00 = Orig5[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*(test00+ugrad)  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*(test00+ugrad) ;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*(test00+ugrad);
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*(test00+ugrad);
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*(test00+ugrad);
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

  
}


// D(u):D(v) TNSECST LPS non-linear terms only
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_S;
  double c0, c3;
  double u1,u2, u1x, u1y, u2x, u2y;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  
  MatrixG11 = LocMatrices[2];
  MatrixG12 = LocMatrices[3];
  MatrixG21 = LocMatrices[4];
  MatrixG22 = LocMatrices[5];
  MatrixG23 = LocMatrices[6];
  MatrixG32 = LocMatrices[7];
  MatrixG33 = LocMatrices[8];
  
  MatrixD11 = LocMatrices[9];
  MatrixD12 = LocMatrices[10];
  MatrixD21 = LocMatrices[11];
  MatrixD22 = LocMatrices[12];
  MatrixD31 = LocMatrices[13];
  MatrixD32 = LocMatrices[14];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // tau_x
  Orig4 = OrigValues[4];         // tau_y
  Orig5 = OrigValues[5];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c3 = coeff[3]; // Tau term in CST 
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
   
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*(2*test10*ansatz10+ test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      {
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix11Row[j] += Mult * val;


      val  = c0*(test10*ansatz10+ 2*test01*ansatz01);
      if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE==NSE)
      { 
      val += (u1*ansatz10+u2*ansatz01)*test00;
      }
      Matrix22Row[j] += Mult * val;   
    }                            // endfor j

  }                              // endfor i
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig3[i];
    test01 = Orig4[i];
    test00 = Orig5[i];

    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00;
    Rhs1[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00;
    Rhs2[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00;
    Rhs3[i] += Mult*val;
    
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
   
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00;
      Matrix32Row[j] += Mult * val;

    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig3[j];
      ansatz01 = Orig4[j];
      ansatz00 = Orig5[j];

      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1*ansatz10 + u2*ansatz01)*test00  ;
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00;
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1*ansatz10 + u2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00;
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00;
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00;
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1*ansatz10 + u2*ansatz01)*test00;
      Matrix33Row[j] += Mult * val;

    }                            // endfor j

   
  }                              // endfor i

  
}


// D(u):D(v) TNSECST LPS non-linear terms only for two phase axial flows
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_S;
  double c0, c3, c8, c9, c10;
  double u1,u2, u1w1, u2w2, u1x, u1y, u2x, u2y, x, r;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  
  MatrixG11 = LocMatrices[2];
  MatrixG12 = LocMatrices[3];
  MatrixG21 = LocMatrices[4];
  MatrixG22 = LocMatrices[5];
  MatrixG23 = LocMatrices[6];
  MatrixG32 = LocMatrices[7];
  MatrixG33 = LocMatrices[8];
  
  MatrixD11 = LocMatrices[9];
  MatrixD12 = LocMatrices[10];
  MatrixD21 = LocMatrices[11];
  MatrixD22 = LocMatrices[12];
  MatrixD31 = LocMatrices[13];
  MatrixD32 = LocMatrices[14];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // tau_x
  Orig4 = OrigValues[4];         // tau_y
  Orig5 = OrigValues[5];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c3 = coeff[3]; // Tau term in CST 
   c8 = coeff[8]; // density ratio
   c9 = coeff[9]; // Fluid type, 1-Newtonian, 2-Oldroyd-B, 3-Giesekus
   c10 = coeff[10];// Tau dot Tau in Giesekus
  
  u1 = param[0];                 // u1old 
  u2 = param[1];                 // u2old 
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold
  u1w1 =  param[15];             // u1old - grid VX   !!check aux parameters
  u2w2 =  param[16];             // u2old - grid VY   !!check aux parameters
  
  x  = param[17]; // x

 if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
   
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*((2*test10*ansatz10+ test01*ansatz01)*r);
      if(TDatabase::ParamDB->Axial3D==1)
     {
      val+=  c0*2.*(ansatz00*test00/r);
     }
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*c8*r;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01)*r ;
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*c8*r;
      Matrix22Row[j] += Mult * val;   
    }                            // endfor j

  }                              // endfor i
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig3[i];
    test01 = Orig4[i];
    test00 = Orig5[i];

    if(c9==2 || c9==3) 
    {
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00*r;
    Rhs1[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00*r;
    Rhs2[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00*r;
    Rhs3[i] += Mult*val;
    }
    else if(c9==1)
    {
     val = 0; 
     Rhs1[i] += Mult*val;
     Rhs2[i] += Mult*val;
     Rhs3[i] += Mult*val;
    }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
     
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      if(c9==2 || c9==3)
      {
      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00*r;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00*r;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix32Row[j] += Mult * val;
      }
      else if(c9==1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix31Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig3[j];
      ansatz01 = Orig4[j];
      ansatz00 = Orig5[j];

      if(c9==2 || c9==3)
      {
      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1w1*ansatz10 + u2w2*ansatz01)*test00*r  ;
      if(c9==3)
     {
      val += c10*tau1*ansatz00*test00*r;
     }
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1w1*ansatz10 + u2w2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00*r;
      if(c9==3)
     {
      val += c10*(tau1 + tau3)*ansatz00*test00*r;
     }
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
      if(c9==3)
     {
      val += c10*tau2*ansatz00*test00*r;
     }
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1w1*ansatz10 + u2w2*ansatz01)*test00*r;
      if(c9==3)
     {
      val += c10*tau3*ansatz00*test00*r;
     }
      Matrix33Row[j] += Mult * val;
      }
      else if(c9==1)
      {
        val  = c3*ansatz00*test00*r;
	Matrix11Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val*2.0;
	Matrix33Row[j] += Mult * val;
	val = 0;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
       else
      {
	cout<<"Invalid fluid type \n";
	exit(4711);
      }
 
    }                            // endfor j

  }                              // endfor i

}



// D(u):D(v) TNSECST LPS non-linear terms only for impinging droplet 3D-axial flows
void Time_NSEType4_Galerkin_CST_Galerkin_NLTerms_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double **MatrixA11, **MatrixA22;
  double **MatrixG11, **MatrixG12, **MatrixG21, **MatrixG22, **MatrixG23, **MatrixG32, **MatrixG33;
  double **MatrixD11, **MatrixD12, **MatrixD21, **MatrixD22, **MatrixD31, **MatrixD32;
  double *Rhs1, *Rhs2, *Rhs3, val;
  double *Matrix11Row, *Matrix12Row, *Matrix13Row, *Matrix21Row, *Matrix22Row, *Matrix23Row, *Matrix31Row, *Matrix32Row, *Matrix33Row;
  double *MatrixRow1, *MatrixRow2;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3, *Orig4, *Orig5;
  int i,j,N_U, N_P, N_S;
  double c0, c3, c8, c9;
  double u1,u2, u1w1, u2w2, u1x, u1y, u2x, u2y, x, r;
  double tau1, tau2, tau3, tau1x, tau2x, tau3x, tau1y, tau2y, tau3y;
    
  MatrixA11 = LocMatrices[0];
  MatrixA22 = LocMatrices[1];
  
  MatrixG11 = LocMatrices[2];
  MatrixG12 = LocMatrices[3];
  MatrixG21 = LocMatrices[4];
  MatrixG22 = LocMatrices[5];
  MatrixG23 = LocMatrices[6];
  MatrixG32 = LocMatrices[7];
  MatrixG33 = LocMatrices[8];
  
  MatrixD11 = LocMatrices[9];
  MatrixD12 = LocMatrices[10];
  MatrixD21 = LocMatrices[11];
  MatrixD22 = LocMatrices[12];
  MatrixD31 = LocMatrices[13];
  MatrixD32 = LocMatrices[14];
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u_x
  Orig1 = OrigValues[1];         // u_y
  Orig2 = OrigValues[2];         // u
  Orig3 = OrigValues[3];         // tau_x
  Orig4 = OrigValues[4];         // tau_y
  Orig5 = OrigValues[5];         // tau

   c0 = coeff[0]; // D(u):D(v) term 
   c3 = coeff[3]; // Tau term in CST 
   c8 = coeff[8]; // Fluid type, 1-Newtonian, 2-Oldroyd-B, 3-Giesekus
   c9 = coeff[9];// Tau dot Tau in Giesekus
  
  u1 = param[0];                 // u1old 
  u2 = param[1];                 // u2old 
  u1x = param[2];                // u1xold
  u2x = param[3];                // u2xold
  u1y = param[4];                // u1yold
  u2y = param[5];                // u2yold
  tau1 = param[6];               // tau1old
  tau2 = param[7];               // tau2old
  tau3 = param[8];               // tau3old
  tau1x = param[9];              // tau1xold
  tau2x = param[10];             // tau2xold
  tau3x = param[11];             // tau3xold
  tau1y = param[12];             // tau1yold
  tau2y = param[13];             // tau2yold
  tau3y = param[14];             // tau3yold
  u1w1 =  param[15];             // u1old - grid VX   !!check aux parameters
  u2w2 =  param[16];             // u2old - grid VY   !!check aux parameters
  
  x  = param[17]; // x

 if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }

    for(i=0;i<N_U;i++)
  {

    test10 = Orig0[i];
    test01 = Orig1[i];
    test00 = Orig2[i];
   
    // Matrix A
    Matrix11Row = MatrixA11[i];
    Matrix22Row = MatrixA22[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];
      
      val  = c0*((2*test10*ansatz10+ test01*ansatz01)*r);
      if(TDatabase::ParamDB->Axial3D==1)
     {
      val+=  c0*2.*(ansatz00*test00/r);
     }
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = c0*(test10*ansatz10+ 2*test01*ansatz01)*r ;
      val += (u1w1*ansatz10+u2w2*ansatz01)*test00*r;
      Matrix22Row[j] += Mult * val;   
    }                            // endfor j

  }                              // endfor i
  
    for(i=0;i<N_S;i++)
   {

    test10 = Orig3[i];
    test01 = Orig4[i];
    test00 = Orig5[i];

    if(c8==2 || c8==3) 
    {
    val = ((u1*tau1x + u2*tau1y) - 2*(tau1*u1x + tau2*u1y) )*test00*r;
    Rhs1[i] += Mult*val;
    val = (2*(u1*tau2x + u2*tau2y) - 2*(tau1*u2x + tau2*u1x + tau2*u2y + tau3*u1y))*test00*r;
    Rhs2[i] += Mult*val;
    val = ((u1*tau3x + u2*tau3y) - 2*(tau2*u2x + tau3*u2y))*test00*r;
    Rhs3[i] += Mult*val;
    }
    else if(c8==1)
    {
     val = 0; 
     Rhs1[i] += Mult*val;
     Rhs2[i] += Mult*val;
     Rhs3[i] += Mult*val;
    }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
     
    // Matrix D
    Matrix11Row = MatrixD11[i];
    Matrix12Row = MatrixD12[i];
    Matrix21Row = MatrixD21[i];
    Matrix22Row = MatrixD22[i];
    Matrix31Row = MatrixD31[i];
    Matrix32Row = MatrixD32[i];
    
    for(j=0;j<N_U;j++)
    {
      ansatz10 = Orig0[j];
      ansatz01 = Orig1[j];
      ansatz00 = Orig2[j];

      if(c8==2 || c8==3)
      {
      val  = (tau1x*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix11Row[j] += Mult * val;

      val  = tau1y*ansatz00*test00*r;
      Matrix12Row[j] += Mult * val;

      val  = (2*tau2x*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix21Row[j] += Mult * val;

      val  = (2*tau2y*ansatz00 - 2*(tau1*ansatz10 + tau2*ansatz01))*test00*r;
      Matrix22Row[j] += Mult * val;
      
      val  = tau3x*ansatz00*test00*r;
      Matrix31Row[j] += Mult * val;

      val  = (tau3y*ansatz00 - 2*(tau2*ansatz10 + tau3*ansatz01))*test00*r;
      Matrix32Row[j] += Mult * val;
      }
      else if(c8==1)
      {
	val = 0;
	Matrix11Row[j] += Mult * val;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val;
	Matrix31Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
      else
     {
	cout<<"Invalid fluid type \n";
	exit(4711);
     }
    }                            // endfor j
    
    // Matrix G
    Matrix11Row = MatrixG11[i];
    Matrix12Row = MatrixG12[i];
    Matrix21Row = MatrixG21[i];
    Matrix22Row = MatrixG22[i];
    Matrix23Row = MatrixG23[i];
    Matrix32Row = MatrixG32[i];
    Matrix33Row = MatrixG33[i];
    
     for(j=0;j<N_S;j++)
    {
      ansatz10 = Orig3[j];
      ansatz01 = Orig4[j];
      ansatz00 = Orig5[j];

      if(c8==2 || c8==3)
      {
      val  = (c3*ansatz00 - 2.0*ansatz00*u1x + u1w1*ansatz10 + u2w2*ansatz01)*test00*r  ;
      if(c8==3)
     {
      val += c9*tau1*ansatz00*test00*r;
     }
      Matrix11Row[j] += Mult * val;

      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix12Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix21Row[j] += Mult * val;
      
      val = (c3*2.0*ansatz00 + 2.0*(u1w1*ansatz10 + u2w2*ansatz01)- 2.0*(u1x + u2y)*ansatz00)*test00*r;
      if(c8==3)
     {
      val += c9*(tau1 + tau3)*ansatz00*test00*r;
     }
      Matrix22Row[j] += Mult * val;
       
      val = -2.0*ansatz00*u1y*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix23Row[j] += Mult * val;
      
      val = -2.0*ansatz00*u2x*test00*r;
      if(c8==3)
     {
      val += c9*tau2*ansatz00*test00*r;
     }
      Matrix32Row[j] += Mult * val;

      val  = (c3*ansatz00  - 2.0*ansatz00*u2y + u1w1*ansatz10 + u2w2*ansatz01)*test00*r;
      if(c8==3)
     {
      val += c9*tau3*ansatz00*test00*r;
     }
      Matrix33Row[j] += Mult * val;
      }
      else if(c8==1)
      {
        val  = c3*ansatz00*test00*r;
	Matrix11Row[j] += Mult * val;
	Matrix22Row[j] += Mult * val*2.0;
	Matrix33Row[j] += Mult * val;
	val = 0;
	Matrix12Row[j] += Mult * val;
	Matrix21Row[j] += Mult * val;
	Matrix23Row[j] += Mult * val;
	Matrix32Row[j] += Mult * val;
      }
       else
      {
	cout<<"Invalid fluid type \n";
	exit(4711);
      }
 
    }                            // endfor j

  }                              // endfor i

}




// ======================================================================
// right-hand side ONLY, for TNSECST, DEVSS/SUPG
// ======================================================================
void Time_NSEType4_Galerkin_CST_SUPG_RhsOnly(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i, N_U, N_S;
  double c1, c2, c4, c5, c6;
  double u1,u2;
  double delta0 = TDatabase::ParamDB->DELTA0, ugrad, delta, norm_u;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // tau_x
  Orig2 = OrigValues[2];         // tau_y
  Orig3 = OrigValues[3];         // tau
  
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c4 = coeff[4];                // f3
  c5 = coeff[5];                // f4 
  c6 = coeff[6];                // f5
  
  u1 = param[0];                 // u1old
  u2 = param[1];                 // u2old

  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                             
  
  norm_u = sqrt((u1*u1) + (u2*u2));
  
  if (norm_u < 0.000001)
    delta=0; 
  else
    delta = delta0*hK/norm_u;
  
    for(i=0;i<N_S;i++)
   {
    test10 = Orig1[i];
    test01 = Orig2[i];
    test00 = Orig3[i];
    ugrad  = delta * (u1*test10+u2*test01);

    Rhs3[i] += Mult*c4*(test00+ugrad);
    Rhs4[i] += Mult*2*c5*(test00+ugrad);
    Rhs5[i] += Mult*c6*(test00+ugrad);
  
   }
}


// ======================================================================
// right-hand side ONLY, for TNSECST, LPS
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i, N_U, N_S;
  double c1, c2, c4, c5, c6;
  double u1,u2;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // tau_x
  Orig2 = OrigValues[2];         // tau_y
  Orig3 = OrigValues[3];         // tau
  
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c4 = coeff[4];                // f3
  c5 = coeff[5];                // f4 
  c6 = coeff[6];                // f5
  
  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1;
    Rhs2[i] += Mult*test00*c2;
  }                             
  
  
    for(i=0;i<N_S;i++)
   {
    test10 = Orig1[i];
    test01 = Orig2[i];
    test00 = Orig3[i];

    Rhs3[i] += Mult*c4*test00;
    Rhs4[i] += Mult*2*c5*test00;
    Rhs5[i] += Mult*c6*test00;
  
   }
}


// ======================================================================
// right-hand side ONLY, for TNSECST, LPS for 2Phase Axial3D flows
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_2PhaseAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i, N_U, N_S;
  double c1, c2, c4, c5, c6, c8;
  double u1,u2, x, r;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // tau_x
  Orig2 = OrigValues[2];         // tau_y
  Orig3 = OrigValues[3];         // tau
  
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c4 = coeff[4];                // f3
  c5 = coeff[5];                // f4 
  c6 = coeff[6];                // f5
  c8 = coeff[8];                // density ratio
  
  x  = param[17]; // x
  
  if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }
  
  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1*c8*r;
    Rhs2[i] += Mult*test00*c2*c8*r;
  }                             
  
  
    for(i=0;i<N_S;i++)
   {
    test10 = Orig1[i];
    test01 = Orig2[i];
    test00 = Orig3[i];

    Rhs3[i] += Mult*c4*test00*r;
    Rhs4[i] += Mult*2*c5*test00*r;
    Rhs5[i] += Mult*c6*test00*r;
  
   }
}

// ======================================================================
// right-hand side ONLY, for TNSECST, LPS for Impinging droplet Axial3D flows
// ======================================================================
void Time_NSEType4_Galerkin_CST_Galerkin_RhsOnly_ImpDropAxial3D(double Mult, double *coeff,
double *param, double hK,
double **OrigValues, int *N_BaseFuncts,
double ***LocMatrices, double **LocRhs)
{
  double *Rhs1, *Rhs2, *Rhs3, *Rhs4, *Rhs5;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2, *Orig3;
  int i, N_U, N_S;
  double c1, c2, c4, c5, c6;
  double u1,u2, x, r;
  
  Rhs1 = LocRhs[0];
  Rhs2 = LocRhs[1];
  Rhs3 = LocRhs[2];
  Rhs4 = LocRhs[3];
  Rhs5 = LocRhs[4];

  N_U = N_BaseFuncts[0];
  N_S = N_BaseFuncts[2];
  
  Orig0 = OrigValues[0];         // u
  Orig1 = OrigValues[1];         // tau_x
  Orig2 = OrigValues[2];         // tau_y
  Orig3 = OrigValues[3];         // tau
  
  c1 = coeff[1];                 // f1
  c2 = coeff[2];                 // f2
  c4 = coeff[4];                // f3
  c5 = coeff[5];                // f4 
  c6 = coeff[6];                // f5

  
  x  = param[17]; // x
  
  if(TDatabase::ParamDB->Axial3D==1)
  {
    r  = fabs(x);
  if(r<1e-12)
   {
   OutPut("check TNSE2D_FixPo  x value zero !!!!! "<< x <<endl);
   OutPut("Quad formula: Change all integral points as internal points"<<endl);
   }
  }
  else if(TDatabase::ParamDB->Axial3D==0)
  { 
//     OutPut("2D Planar case activated !!!!!!!!! "<<endl);
    r = 1.0;
  }
  else
  {
    OutPut("Invalid value for Axial3D in dat file !!!!! "<<endl);
    exit(4711);
  }
  
  for(i=0;i<N_U;i++)
  {
    test00 = Orig0[i];

    Rhs1[i] += Mult*test00*c1*r;
    Rhs2[i] += Mult*test00*c2*r;
  }                             
  
  
    for(i=0;i<N_S;i++)
   {
    test10 = Orig1[i];
    test01 = Orig2[i];
    test00 = Orig3[i];

    Rhs3[i] += Mult*c4*test00*r;
    Rhs4[i] += Mult*2*c5*test00*r;
    Rhs5[i] += Mult*c6*test00*r;
  
   }
}
