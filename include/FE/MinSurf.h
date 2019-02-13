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
   
// =====================================================================
// for 2D part
// =====================================================================
// part for standard Galerkin
int N_Terms2d = 3;
MultiIndex2D Derivatives2d[3] = { D10, D01, D00 };
int SpacesNumbers2d[3] = { 0, 0, 0 };
int N_Matrices2d = 1;
int RowSpace2d[2] = { 0 };
int ColumnSpace2d[2] = { 0 };
int N_Rhs2d = 1;
int RhsSpace2d[1] = { 0 };

MultiIndex2D AllDerivatives2d[3] = { D00, D10, D01 };

void L2H1Errors2d(int N_Points, double *X, double *Y, double *AbsDetjk, 
                double *Weights, double hK, 
                double **Der, double **Exact,
                double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

void BilinearAssemble2d(double Mult, double *coeff, double hK, 
                    double **OrigValues, int *N_BaseFuncts,
                    double ***LocMatrices, double **LocRhs)
{
  double **Matrix, *Rhs, val, *MatrixRow;
  double ansatz00, ansatz10, ansatz01;
  double test00, test10, test01;
  double *Orig0, *Orig1, *Orig2;
  int i,j, N_;
  double c0, c1, c2, c3, c4; 

  Matrix = LocMatrices[0];
  Rhs = LocRhs[0];

  N_ = N_BaseFuncts[0];

  Orig0 = OrigValues[0];
  Orig1 = OrigValues[1];
  Orig2 = OrigValues[2];

  // coefficients
  c0 = coeff[0];  // eps
  c1 = coeff[1];  // b_1
  c2 = coeff[2];  // b_2
  c3 = coeff[3];  // c
  c4 = coeff[4];  // f

  for(i=0;i<N_;i++)
  {
    MatrixRow = Matrix[i];
    // test function
    test10 = Orig0[i];  // xi derivative
    test01 = Orig1[i];  // eta derivative
    test00 = Orig2[i];  // function 

    // assemble rhs 
    // quad_weigth * test_function * f 
    Rhs[i] += Mult*test00*c4;

    for(j=0;j<N_;j++)
    {
      ansatz10 = Orig0[j]; // xi derivative
      ansatz01 = Orig1[j]; // eta derivative
      ansatz00 = Orig2[j]; // function 

      // assemble viscous term
      // eps (test_x ansatz_x + test_y ansatz_y) 
      val = c0*(test10*ansatz10+test01*ansatz01);
      // assemble convective term
      // (b_1 ansatz_x + b_2 ansatz_y) test
      val += (c1*ansatz10+c2*ansatz01)*test00;
      // assembel reactive term 
      // c  ansatz test
      val += c3*ansatz00*test00;
      
      // quad weigth
      val *= Mult;

      // update matrix entry
      MatrixRow[j] += val;
    } // endfor j
  } // endfor i
}

// exact solution
void Exact2d(double x, double y, double *values)
{
  values[0] = sin(Pi*x)*sin(2*Pi*y);
  values[1] = Pi*cos(Pi*x)*sin(2*Pi*y);
  values[2] = 2*Pi*sin(Pi*x)*cos(2*Pi*y);
  values[3] = 0;
}

void ExactX(double x, double y, double *values)
{
  values[0] = x;
  values[1] = 1;
  values[2] = 0;
  values[3] = 0;
}

void ExactY(double x, double y, double *values)
{
  values[0] = y;
  values[1] = 0;
  values[2] = 1;
  values[3] = 0;
}

void Exact3D(double x, double y, double z, double *values)
{
  values[0] = 4*x*y*z*z;
  values[1] = 4*y*z*z;
  values[2] = 4*x*z*z;
  values[3] = 8*x*y*z;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition2d(int i, double t, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void BoundValue2d(int BdComp, double Param, double &value)
{
  value = 0;
}

void BilinearCoeffs2d(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  static int first = 1;

  double phin;
  double chi = TDatabase::ParamDB->P2 - 1;
  double H0 = TDatabase::ParamDB->P4;
  double a = TDatabase::ParamDB->P6;
  double sigma = TDatabase::ParamDB->P7;
  double g = TDatabase::ParamDB->P8;
  double rho = TDatabase::ParamDB->P9;

  double lambda = a*sqrt(rho*g/sigma);
  double Si = 4*Pi*0.1*chi*chi*H0*H0/(2*sqrt(sigma*g*rho));

  if(first)
  {
    OutPut("lambda: " << lambda << endl);
    OutPut("Si: " << Si << endl);

    first = 0;
  }

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = param[0];
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = lambda*lambda;

    phin = param[2]*param[5] + param[3]*param[6] + param[4]*param[7];
    coeff[4] = lambda*Si/chi*( param[2]*param[2] + param[3]*param[3]
                              +param[4]*param[4])
              +lambda*Si * phin*phin;
  }
}

// for nonlinear magnetisation law
void BilinearCoeffs2dNL(int n_points, double *x, double *y,
        double **parameters, double **coeffs)
{
  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, *param;
  static int first = 1;

  double phin, Im, H, L, arg;

  double chi   = TDatabase::ParamDB->P2 - 1;
  double H0    = TDatabase::ParamDB->P4;
  double Ms    = TDatabase::ParamDB->P5;
  double a     = TDatabase::ParamDB->P6;
  double sigma = TDatabase::ParamDB->P7;
  double g     = TDatabase::ParamDB->P8;
  double rho   = TDatabase::ParamDB->P9;

  double gamma = 3*chi*H0/Ms;

  double lambda = a*sqrt(rho*g/sigma);
  double Si = 4*Pi*0.1*Ms*Ms/(3*chi*sqrt(sigma*g*rho));

  if(first)
  {
    OutPut("lambda: " << lambda << endl);
    OutPut("Si: " << Si << endl);
    OutPut("gamma: " << gamma << endl);

    first = 0;
  }

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    param = parameters[i];

    coeff[0] = param[0];
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = lambda*lambda;

    H = sqrt(param[2]*param[2] + param[3]*param[3] + param[4]*param[4]);
    arg = gamma*H;
    Im = log( sinh(arg) / arg );

    L = 1/tanh(arg) - 1/arg;

    phin = (param[2]*param[5] + param[3]*param[6] + param[4]*param[7])/H;

    coeff[4] = lambda*Si*Im + 1.5*chi*lambda*Si*L*L*phin*phin;
  }
}
// =====================================================================
// for 2D part (END)
// =====================================================================

// Only the fe values are correct, the 3d parameters can be only used
// in the bilinearcoeff subroutine
void MinSurfParams(double *in, double *out)
{
  double n1, n2, n3;
  out[0] = 1/sqrt(1+in[2]*in[2]+in[3]*in[3]); // 1/sqrt(1+|grad u|^2)
}

int N_FESpaceMS = 1;
int N_FEFunctionsMS = 1;
int N_ParamFctMS = 1;
int N_FEValuesMS = 2;
int FEValuesFctIndexMS[2] = { 0, 0 };
MultiIndex2D FEValuesMultiIndexMS[2] = { D10, D01 };
int N_ParametersMS = 1;
int BeginParameterMS[1] = { 0 };
ParamFct *ParameterFctMS[1] = { MinSurfParams };

// generate undisturbed interface at z = height/2
void UnDisturbed(double x, double y, double *values)
{
  static double middle = 0.5*TDatabase::ParamDB->DRIFT_Z;

  values[0] = middle;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}
