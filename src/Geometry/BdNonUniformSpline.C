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
// @(#)BdNonUniformSpline.C        1.2 07/16/99
//
// Class:       BdNonUniformSpline
// Purpose:     spline function as a component of a boundary part
//
// Author:      Volker Behns  18.06.97
//
// =======================================================================

#include <BdNonUniformSpline.h>
#include <math.h>

// Constructor
TBdNonUniformSpline::TBdNonUniformSpline(int id, int N_Spls) : TBoundComp2D(id)
{
  Type = NonUniformSpline2D;
  N_Splines = N_Spls;
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];
}

// Methods
void TBdNonUniformSpline::SetParams (double *params)
{
  Params = params;
  //memcpy(,,);
}

int TBdNonUniformSpline::GetN_Splines ()
{
  return N_Splines;
}

int TBdNonUniformSpline::GetXYofT(double T, double &X, double &Y)
{
  double phi1, phi2, phi3, phi4;
  int i, ISpline;

  Param9[0] = 0;
  for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];
    
  for(i=1;i<=N_Splines;i++) 
  {
    ISpline = (i-1)*10;
    if((T>=Param9[i-1]) && (T<=Param9[i]))
    {
      // further T must be from [0;1] on a subspline
      T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
      break;
    }
  }
  
  phi1 = (2.*T*T - 3.*T)*T + 1.;
  phi2 = (-2.*T + 3.)*T*T;
  phi3 = (T*T - 2.*T + 1.)*T;
  phi4 = (T - 1)*T*T;

  X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
      Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
  Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
      Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

  // cout<<ISpline<<' '<<T<<' '<<X<<' '<<Y<<endl;
  return 0;
}

int TBdNonUniformSpline::GetTofXY(double X, double Y, double &T)
{
  int i,ISpline;
  
  //cout<<"GetTofXY X and Y:"<<X<<' '<<Y<<endl;
  if(ABS(Params[0]-X)<1e-8 && ABS(Params[1]-Y)<1e-8) 
    {
      T = 0.0;
      //cout<<"GetTofXY T:"<<T<<endl;
      return 0; 
    }

  for(i=0;i<N_Splines;i++) // (X,Y) coincides with a spline-point
    {
      ISpline = i*10;  
      if(ABS(Params[ISpline+2]-X)<1e-8 && ABS(Params[ISpline+3]-Y)<1e-8) 
	{
	  T = Params[ISpline+8];
	  //cout<<"GetTofXY T:"<<T<<endl;
	  return 0; 
	}
    }
  for(i=0;i<N_Splines;i++) // (X,Y) lies at the middle of a subspline
    {
      ISpline = i*10;  
      if(ABS(GetLocalXofT(i+1,0.5)-X)<1e-8 && 
	 ABS(GetLocalYofT(i+1,0.5)-Y)<1e-8) 
	{
	  if(i==0)
	    T = 0.5*Params[ISpline+8];
	  else
	    T = 0.5*(Params[(i-1)*10+8]+Params[ISpline+8]);
	  //cout<<"GetTofXY T:"<<T<<endl;
	  return 0; 
	}
    }
  return -1;
}

double TBdNonUniformSpline::GetLocalXofT(int ISpline, double T)
{
  double phi1, phi2, phi3, phi4, X;

  ISpline = (ISpline-1)*10;
    
  phi1 = (2.*T*T - 3.*T)*T + 1.;
  phi2 = (-2.*T + 3.)*T*T;
  phi3 = (T*T - 2.*T + 1.)*T;
  phi4 = (T - 1)*T*T;

  X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
      Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
 
  return X;
}

double TBdNonUniformSpline::GetLocalYofT(int ISpline, double T)
{
  double phi1, phi2, phi3, phi4, Y;

  ISpline = (ISpline-1)*10;
    
  phi1 = (2.*T*T - 3.*T)*T + 1.;
  phi2 = (-2.*T + 3.)*T*T;
  phi3 = (T*T - 2.*T + 1.)*T;
  phi4 = (T - 1)*T*T;

  Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
      Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;
 
  return Y;
}

void TBdNonUniformSpline::GenerateParams1(double *x, double *y, 
                                          double dx0, double dy0, double dx1, double dy1)
{
  double *h, *t;
  double *a, *b, *c;
  double *rhs, *Mx, *My;
  int i, ISpline;
  
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  
  h[0] = 0.; t[0] = 0.;
  for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }
  
  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;
    b[i] = h[i]/(h[i]+h[i+1]);
    c[i] = h[i+1]/(h[i]+h[i+1]);
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);
  
  Solver_3diag(a, b, c, rhs, Mx);
  /*
  cout<<"Mx"<<endl;
  for(i=0;i<=N_Splines;i++)
    {
      cout<<setprecision(15)<<i<<' '<<Mx[i]<<endl;
    }
  */
  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);
  
  Solver_3diag(a, b, c, rhs, My);
  /*
  cout<<"My"<<endl;
  for(i=0;i<=N_Splines;i++)
    {
      cout<<setprecision(15)<<i<<' '<<My[i]<<endl;
    }
  */
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. + 
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1]; 
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. + 
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1]; 
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. + 
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1]; 
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. + 
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1]; 
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;
  }
  /*
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    //cout<<i+1<< " subspline"<<endl;
    cout<<setprecision(15);
    cout<<"  "<<Params[ISpline    ]<<'\t'<<Params[ISpline + 1]<<endl;
    cout<<"  "<<Params[ISpline + 2]<<'\t'<<Params[ISpline + 3]<<endl;
    cout<<"  "<<Params[ISpline + 4]<<'\t'<<Params[ISpline + 5]<<endl;
    cout<<"  "<<Params[ISpline + 6]<<'\t'<<Params[ISpline + 7]<<endl;
    cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  */
  delete h; delete t; delete a; delete b; delete c; delete rhs; delete Mx; delete My;
}

int TBdNonUniformSpline::ReadIn(std::ifstream &dat)
{
  char line[100];
  int i, N_Params;

  N_Params = 5 * N_Splines;

  for(i=0;i<N_Params;i++)
  {
    dat >> Params[2*i] >> Params[2*i+1];
    dat.getline (line, 99);
  }
             
  return 0;
}
void TBdNonUniformSpline::Solver_3diag(double *a, double *b, double *c, double *rhs, double *sol)
{
  double *alpha, *beta, *y;
  int i, N;
  
  N = N_Splines+1;
  alpha = new double[N]; beta = new double[N]; y = new double[N];
  
  alpha[0] = a[0]; y[0] = rhs[0];
  for(i=1;i<N;i++)
  {
    beta[i] = b[i]/alpha[i-1];  
    alpha[i] = a[i]-beta[i]*c[i-1];
    y[i] = rhs[i]-beta[i]*y[i-1];
  }

  sol[N-1] = y[N-1]/alpha[N-1];
  for(i=N-2;i>=0;i--)
    sol[i] = (y[i]-c[i]*sol[i+1])/alpha[i];
    
  delete alpha; delete beta; delete y;
}

/** for BEM */
/**  (x-x_middle)/(t-0.5) */
double TBdNonUniformSpline::AdXofT(int ISpline, double T)
{
  double X;
  ISpline = (ISpline-1)*10;
  
  X = -Params[ISpline    ]*1.5  + Params[ISpline + 2]*1.5 +
      -Params[ISpline + 4]*0.25 - Params[ISpline + 6]*0.25 +
     (-Params[ISpline + 4]      + Params[ISpline + 6])*0.5*(T-0.5)+
     ( Params[ISpline    ]*2.   - Params[ISpline + 2]*2.+
       Params[ISpline + 4]      + Params[ISpline + 6])*(T-0.5)*(T-0.5);
  
  return X;
}

/** for BEM */
/** (y-y_middle)/(t-0.5) */
double TBdNonUniformSpline::AdYofT(int ISpline, double T)
{
  double Y;

  ISpline = (ISpline-1)*10;
  
  Y = -Params[ISpline + 1]*1.5  + Params[ISpline + 3]*1.5 +
      -Params[ISpline + 5]*0.25 - Params[ISpline + 7]*0.25 +
     (-Params[ISpline + 5]      + Params[ISpline + 7])*0.5*(T-0.5)+
     ( Params[ISpline + 1]*2.   - Params[ISpline + 3]*2.+
       Params[ISpline + 5]      + Params[ISpline + 7])*(T-0.5)*(T-0.5);
 
  return Y;
}

/** for BEM */
/** x'(t) */
double TBdNonUniformSpline::GetLocalDXofT(int ISpline, double T)
{
  double dphi1, dphi2, dphi3, dphi4, X;

  ISpline = (ISpline-1)*10;
    
  dphi1 =  6.*T*(T-1);
  dphi2 = -6.*T*(T-1);
  dphi3 = 3.*T*T-4.*T+1;
  dphi4 = T*(3.*T-2);

  X = Params[ISpline    ]*dphi1 + Params[ISpline + 2]*dphi2 +
      Params[ISpline + 4]*dphi3 + Params[ISpline + 6]*dphi4;
 
  return X;
}

/** for BEM */
/** y'(t) */
double TBdNonUniformSpline::GetLocalDYofT(int ISpline, double T)
{
  double dphi1, dphi2, dphi3, dphi4, Y;

  ISpline = (ISpline-1)*10;
    
  dphi1 =  6.*T*(T-1);
  dphi2 = -6.*T*(T-1);
  dphi3 = 3.*T*T-4.*T+1;
  dphi4 = T*(3.*T-2);

  Y = Params[ISpline + 1]*dphi1 + Params[ISpline + 3]*dphi2 +
      Params[ISpline + 5]*dphi3 + Params[ISpline + 7]*dphi4;
 
  return Y;
}
