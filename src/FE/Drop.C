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
   
// free-surface calculation for the drop problem

#include <string.h>
#include <Drop.h>

//#include <LinAlg.h>

//for f the information at all interface points is used
int SchemeA(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a)
{
  double h = 1./(N-1), eps = 1.e-6;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp, sum;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it, i_err_r, i_err_z;
  int K=0; // K=1 for axisymmetrical case, K=0 for plane case
  // for reometer problem
  double alpha = Pi/4.;
  // for adaptation
  double *t; //uniform grid
  double *s; //nonuniform grid
  double *ds, *dsm; // derivatives from s(a,t) on t at vertices and at the midpoints
  if(adap)
    {    
      t = new double[N]; s = new double[N]; ds = new double[N]; dsm = new double[N];
      for(i=0;i<N;i++)
	{
	  t[i] = i*h;
	  s[i] = S(*a,t[i]);
	  ds[i] = dSdt(*a,t[i]);
	  if(i!=0) dsm[i] = dSdt(*a,t[i]-h/2.);
	}
    }
  // ---------------
  f = new double[2*N-1]; f1 = new double[2*N-1]; f2 = new double[2*N-1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];
    
// it is necessary, when the input space coordinates r and z
// are dimentionless on D^1/2  
  
  //for(i=0;i<N;i++)
  //{
  //  r[i]/=L; z[i]/=L;
  //}

  it = 0;
  err_r = 0;
  err_z = 0;
  while(it<5 || err_r>eps || err_z>eps)
  {
    it++;
       
    L = 0.;
    for(i=1;i<N;i++)
      L+=(z[i-1]+z[i])/2.*(r[i]-r[i-1])/h;
    L = pow(h*L,-1./2.)/2.;
    // for reometer problem
    L = 1./r[N-1];

    OutPut("L="<<L<<endl);
  // F = f - K*z'/r + C 
  // f = -2WL/3/hi*f1-WL*f2
  // ----------------------------------------------
  // f[i] the value on the free surface:
  // if i=2j - at the point, i=2j+1 - in the middle
  
  // for the case of magnetic saturation 
    /*
    memset(f1,0,(2*N-1)*SizeOfDouble);
    f2[0] = 1.; f2[2*N-2] = 0.;
    for(i=1;i<2*N-2;i++)
      if(i%2==0)
	if(adap)
	  f2[i] = pow((r[i/2+1]-r[i/2-1])/(s[i/2+1]-s[i/2-1]),2);
	else
	  f2[i] = pow((r[i/2+1]-r[i/2-1])/(2.*h),2);
      else
	if(adap)
	  f2[i] = pow((r[i/2+1]-r[i/2])/(s[i/2+1]-s[i/2]),2);
	else
	  f2[i] = pow((r[i/2+1]-r[i/2])/h,2);
    */
  // taking into considaration the magnetic field
    // for reometer problem
    for(i=0;i<N;i++){
      h_r[i] = 0;
      h_z[i] = 60 - 40*r[i]/r[N-1];
      h_r[i] = (60 - 40*r[i]/r[N-1])*cos(Pi/180.*(90+15));
      h_z[i] = (60 - 40*r[i]/r[N-1])*sin(Pi/180.*(90+15));
    }

    for(i=0;i<2*N-1;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
	// for reometer problem
	Hn = h_z[0]*sin(alpha);
	H = fabs(h_z[0]);
      }
      else if(i==2*N-2)
      {
        Hn = 0;
        H = h_z[N-1];
	// for reometer problem
	Hn = h_z[N-1]*sin(alpha);
	H = fabs(h_z[N-1]);
      }
      else if(i%2==0)
        {
	  if(adap)
	    Hn = -(z[i/2+1]-z[i/2-1])/(s[i/2+1]-s[i/2-1])*h_r[i/2] + 
	      (r[i/2+1]-r[i/2-1])/(s[i/2+1]-s[i/2-1])*h_z[i/2];
	  else
	    Hn = -(z[i/2+1]-z[i/2-1])/(2.*h)*h_r[i/2] + 
	      (r[i/2+1]-r[i/2-1])/(2.*h)*h_z[i/2];
	  H = sqrt(h_r[i/2]*h_r[i/2] + h_z[i/2]*h_z[i/2]);
        }
        else
        {
	  if(adap)
	    Hn = -(z[i/2+1]-z[i/2])/(s[i/2+1]-s[i/2])*(h_r[i/2+1]+h_r[i/2])/2. + 
	      (r[i/2+1]-r[i/2])/(s[i/2+1]-s[i/2])*(h_z[i/2+1]+h_z[i/2])/2.;
	  else
	    Hn = -(z[i/2+1]-z[i/2])/h*(h_r[i/2+1]+h_r[i/2])/2. + 
	      (r[i/2+1]-r[i/2])/h*(h_z[i/2+1]+h_z[i/2])/2.;
	  H = sqrt(pow((h_r[i/2+1]+h_r[i/2])/2.,2) + 
		   pow((h_z[i/2+1]+h_z[i/2])/2.,2));
        }
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      // for nonlinear case      
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    }          
    //OutPut("RHS "<<endl);
    for(i=0;i<2*N-1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i];
	//OutPut(f[i]<<endl);
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
//  ----------------------------------------------  
      
    C = 0.;
    if(adap)
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])*f[2*i-1];
	C = -pow(2.,K)/r[N-1]*(1. + C/pow(r[N-1],K));
      }
    else
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[2*i-1];
	C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));
      }
    // for reometer problem
    C = 0.;
    for(i=1;i<N;i++)
	  C+=(r[i]-r[i-1])*f[2*i-1];
    C = L*(-C+2*cos(alpha));
    OutPut("C="<<C<<endl);

    err = 0.; err_r = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h+
		  (z[i+1]-z[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((r[i-1]-2.*r[i]+r[i+1])/h/h+
		  (z[i+1]-z[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_r)
      { 
        err_r = err;
        i_err_r = i;
      }
    }
    
    err = 0.; err_z = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h-
		  (r[i+1]-r[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((z[i-1]-2.*z[i]+z[i+1])/h/h-
		  (r[i+1]-r[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_z)
      { 
        err_z = err;
        i_err_z = i;
      }
    }
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      if(adap)
	c[i] = 1./dsm[i];
      else
	c[i] = 1.;
 // d
    for(i=1;i<N-1;i++)
      if(adap)
	d[i] = -1.*(1./dsm[i]+1./dsm[i+1]);
      else
	d[i] = -2.;
 // e
    for(i=1;i<N-1;i++)
      if(adap)
	e[i] = 1./dsm[i+1];
      else
	e[i] = 1.;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1.;
    d[N-1] = 1.;
    e[0] = 1.;
    // for reometer problem
    e[0] = 1.;   
    d[0] = -1.;    
    c[N-1] = 0.;
    d[N-1] = 1.;
// b
    if(adap)
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C)*pow(ds[0],2);
    else
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C);
    // for reometer problem
    b[0] = -h*cos(alpha)+h*h/2.*sin(alpha)*(f[0]+C);

    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])*(1.-tau)+
	  h*h*tau*(r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C);
      else
	b[i] = (z[i-1]-2.*z[i]+z[i+1])*(1.-tau)+
	  h*h*tau*(r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);
    b[N-1] = 0.;
    // for reometer problem
    sum = 0.;
    for(i=1;i<N;i++)
      sum+=(r[i-1]+r[i])/2.*(z[i]-z[i-1])/h*h;
    b[N-1] = 0.006/0.01/L+L*sum;

    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);
    
// for r
    c[N-1] = -1.;
    d[0] = 1.;
    d[N-1] = 1.;
    e[0] = 0;
    // for reometer problem
    e[0] = 0.;   
    d[0] = 1.;    
    c[N-1] = -1.;
    d[N-1] = 1.;
// b
    b[0] = 0;
    // for reometer problem
    b[0] = 0.;

    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])*(1.-tau)-
	  h*h*tau*(z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C);
      else
	b[i] = (r[i-1]-2.*r[i]+r[i+1])*(1.-tau)-
	  h*h*tau*(z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);

    if(adap)    
      b[N-1] = -h*h/2.*(f[2*N-2] + K/r[N-1] + C)*pow(ds[N-1],2);
    else
      b[N-1] = -h*h/2.*(f[2*N-2] + K/r[N-1] + C);
    // for reometer problem
    b[N-1] = h*sin(alpha)+h*h/2.*cos(alpha)*(f[2*N-2]+C);
    
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    if(it>10000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
    
    //OutPut("Surface"<<endl);
    //for(i=0;i<N;i++)
    //  OutPut(r[i]<<' '<<z[i]<<endl);
    OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);
  }
  OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);

  //for(i=0;i<N;i++)
  //{
  //  r[i]*=L; z[i]*=L;
  //}

  *LL = L;
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(L),h);
      OutPut("curv = "<<-(f[0]+C)/L<<" a = "<<*a<<" h(0) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }

  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  if(adap) {delete t; delete s; delete ds; delete dsm;}
  return ID;
}
int SchemeA_ax(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a)
{
  double h = 1./(N-1), eps = 1.e-7;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it, i_err_r, i_err_z;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case

  // for adaptation
  double *t; //uniform grid
  double *s; //nonuniform grid
  double *ds, *dsm; // derivatives from s(a,t) on t at vertices and at the midpoints
  if(adap)
    {    
      t = new double[N]; s = new double[N]; ds = new double[N]; dsm = new double[N];
      for(i=0;i<N;i++)
	{
	  t[i] = i*h;
	  s[i] = S(*a,t[i]);
	  ds[i] = dSdt(*a,t[i]);
	  if(i!=0) dsm[i] = dSdt(*a,t[i]-h/2.);
	}
    }
  // ---------------
  f = new double[2*N-1]; f1 = new double[2*N-1]; f2 = new double[2*N-1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];
  
// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3
  
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
  
  it = 0;
  err_r = 0;
  err_z = 0;
  while(it<5 || err_r>eps || err_z>eps)
  {
    it++;
       
     L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);

  // F = f - K*z'/r + C 
  // f = -2WL/3/hi*f1-WL*f2
  // ----------------------------------------------
  // f[i] the value on the free surface:
  // if i=2j - at the point, i=2j+1 - in the middle
  
  // for the case of magnetic saturation 
    /*
    memset(f1,0,(2*N-1)*SizeOfDouble);
    f2[0] = 1.; f2[2*N-2] = 0.;
    for(i=1;i<2*N-2;i++)
      if(i%2==0)
	if(adap)
	  f2[i] = pow((r[i/2+1]-r[i/2-1])/(s[i/2+1]-s[i/2-1]),2);
	else
	  f2[i] = pow((r[i/2+1]-r[i/2-1])/(2.*h),2);
      else
	if(adap)
	  f2[i] = pow((r[i/2+1]-r[i/2])/(s[i/2+1]-s[i/2]),2);
	else
	  f2[i] = pow((r[i/2+1]-r[i/2])/h,2);
    */
  // taking into considaration the magnetic field
    if(it%1==0) OutPut(it<<" Hn"<<endl);
    for(i=0;i<2*N-1;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==2*N-2)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else if(i%2==0)
        {
	  if(adap)
	    Hn = -(z[i/2+1]-z[i/2-1])/(s[i/2+1]-s[i/2-1])*h_r[i/2] + 
	      (r[i/2+1]-r[i/2-1])/(s[i/2+1]-s[i/2-1])*h_z[i/2];
	  else
	    Hn = -(z[i/2+1]-z[i/2-1])/(2.*h)*h_r[i/2] + 
	      (r[i/2+1]-r[i/2-1])/(2.*h)*h_z[i/2];
	  H = sqrt(h_r[i/2]*h_r[i/2] + h_z[i/2]*h_z[i/2]);
	  if(it%1==0) OutPut(z[i/2]<<' '<<Hn<<endl);
        }
        else
        {
	  if(adap)
	    Hn = -(z[i/2+1]-z[i/2])/(s[i/2+1]-s[i/2])*(h_r[i/2+1]+h_r[i/2])/2. + 
	      (r[i/2+1]-r[i/2])/(s[i/2+1]-s[i/2])*(h_z[i/2+1]+h_z[i/2])/2.;
	  else
	    Hn = -(z[i/2+1]-z[i/2])/h*(h_r[i/2+1]+h_r[i/2])/2. + 
	      (r[i/2+1]-r[i/2])/h*(h_z[i/2+1]+h_z[i/2])/2.;
	  H = sqrt(pow((h_r[i/2+1]+h_r[i/2])/2.,2) + 
		   pow((h_z[i/2+1]+h_z[i/2])/2.,2));
	  if(it%1==0) OutPut((z[i/2+1]+z[i/2])/2.<<' '<<Hn<<endl);
        }
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      
      // for nonlinear case      
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    } 
    //if(it%1000==0) OutPut("rhs"<<endl);
    for(i=0;i<2*N-1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i]; 
	//if(it%1000==0) OutPut(i<<' '<<f[i]<<endl);
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
//  ----------------------------------------------  
      
    C = 0.;
    if(adap)
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])*f[2*i-1];
	C = -pow(2.,K)/r[N-1]*(1. + C/pow(r[N-1],K));
      }
    else
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[2*i-1];
	C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));
      }

    err = 0.; err_r = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h+
		  (z[i+1]-z[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((r[i-1]-2.*r[i]+r[i+1])/h/h+
		  (z[i+1]-z[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_r)
      { 
        err_r = err;
        i_err_r = i;
      }
    }
    
    err = 0.; err_z = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h-
		  (r[i+1]-r[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((z[i-1]-2.*z[i]+z[i+1])/h/h-
		  (r[i+1]-r[i-1])/(2.*h)*
		  (f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_z)
      { 
        err_z = err;
        i_err_z = i;
      }
    }
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      if(adap)
	c[i] = 1./dsm[i]/h/h;
      else
	c[i] = 1./h/h;
 // d
    for(i=1;i<N-1;i++)
      if(adap)
	d[i] = -1.*(1./dsm[i]+1./dsm[i+1])/h/h;
      else
	d[i] = -2./h/h;
 // e
    for(i=1;i<N-1;i++)
      if(adap)
	e[i] = 1./dsm[i+1]/h/h;
      else
	e[i] = 1./h/h;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1./h;
    d[N-1] = 1.;
    e[0] = 1./h;
// b
    if(adap)
      b[0] = h/2.*(f[0] - K*(f[0]+C)/2. + C)*pow(ds[0],2);
    else
      b[0] = h/2.*(f[0] - K*(f[0]+C)/2. + C);

    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h*(1.-tau)+
	  tau*(r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C);
      else
	b[i] = (z[i-1]-2.*z[i]+z[i+1])/h/h*(1.-tau)+
	  tau*(r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);
    b[N-1] = 0.;
    
    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);

// for r
    c[N-1] = -1./h;
    d[0] = 1.;
    d[N-1] = 1./h;
    e[0] = 0;
// b
    b[0] = 0;
    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h*(1.-tau)-
	  tau*(z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C);
      else
	b[i] = (r[i-1]-2.*r[i]+r[i+1])/h/h*(1.-tau)-
	  tau*(z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);

    if(adap)    
      b[N-1] = -h/2.*(f[2*N-2] + K/r[N-1] + C)*pow(ds[N-1],2);
    else
      b[N-1] = -h/2.*(f[2*N-2] + K/r[N-1] + C);
    
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    if(it>100000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
  }
  OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);
  
  for(i=0;i<N;i++)
  {
    r[i]*=L; z[i]*=L;
  }
  
  *LL = L;
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(2.*L),h); 
      OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }

  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  if(adap) {delete t; delete s; delete ds; delete dsm;}
  return ID;
}

// the information about field only in the middle of edges and at last points
int SchemeA_ax1(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a)
{
  double h = 1./(N-1), eps = 1.e-8;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it, i_err_r, i_err_z;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case

  // for adaptation
  double *t; //uniform grid
  double *s; //nonuniform grid
  double *ds, *dsm; // derivatives from s(a,t) on t at vertices and at the midpoints
  if(adap)
    {    
      t = new double[N]; s = new double[N]; ds = new double[N]; dsm = new double[N];
      for(i=0;i<N;i++)
	{
	  t[i] = i*h;
	  s[i] = S(*a,t[i]);
	  ds[i] = dSdt(*a,t[i]);
	  if(i!=0) dsm[i] = dSdt(*a,t[i]-h/2.);
	}
    }
  // ---------------
  f = new double[N+1]; f1 = new double[N+1]; f2 = new double[N+1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];
    
// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3
  
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
  
  it = 0;
  err_r = 0;
  err_z = 0;
  //while((W<5.7)&&(it<10) || (W>=5.7)&&(W<5.724)&&(it<50) || 
  //	(W>=5.724)&&(it<5 || err_r>eps || err_z>eps)) 
  //while(it<5 || err_r>eps || err_z>eps) //FEM+BEM
  //while(it<1)
    //while(((W<5.724)&&(it<1)) || ((W>=5.724)&&(it<1))) //BEM+BEM
  //while(it<50)
  while(  (hi==5)  && (it<1) || 
	(hi==15) && ( (W<=4)&&(it<5 || err_r>eps || err_z>eps) || (W>4)&&(it<1) ) ||
	(hi==21) && (it<1)) 
  {
    it++;
       
     L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);
    //OutPut("L "<<it<<' '<<L<<endl);

    for(i=0;i<N+1;i++)
    {
      Hn = h_n[i];
      H = sqrt(h_s[i]*h_s[i] + h_n[i]*h_n[i]);
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      
      // for nonlinear case      
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    } 
   
    //if(W>=5.724 && it==1) OutPut("F "<<it<<endl);
    for(i=0;i<N+1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i]; 
	//if(W>=5.724 && it==1 && i==0) OutPut(z[i]<<' '<<f[i]+C<<endl);
	//if(W>=5.724 && it==1 && i>0 && i!=N) OutPut((z[i-1]+z[i])/2.<<' '<<f[i]+C<<endl);
	//if(W>=5.724 && it==1 && i==N) OutPut(z[N-1]<<' '<<f[N]+C<<endl);
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
//  ----------------------------------------------  
      
    C = 0.;
    if(adap)
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + C/pow(r[N-1],K));
      }
    else
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));
      }
    //OutPut("C/L "<<C/L<<endl);

    err = 0.; err_r = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h+
		  ((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. + 
		  (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((r[i-1]-2.*r[i]+r[i+1])/h/h+
		  ((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		  (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_r)
      { 
        err_r = err;
        i_err_r = i;
      }
    }
    
    err = 0.; err_z = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h-
		  ((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. -
		  (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((z[i-1]-2.*z[i]+z[i+1])/h/h-
		  ((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. -
		  (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_z)
      { 
        err_z = err;
        i_err_z = i;
      }
    }
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      if(adap)
	c[i] = 1./dsm[i];
      else
	c[i] = 1.;
 // d
    for(i=1;i<N-1;i++)
      if(adap)
	d[i] = -1.*(1./dsm[i]+1./dsm[i+1]);
      else
	d[i] = -2.;
 // e
    for(i=1;i<N-1;i++)
      if(adap)
	e[i] = 1./dsm[i+1];
      else
	e[i] = 1.;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1.;
    d[N-1] = 1.;
    e[0] = 1.;
// b
    if(adap)
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C)*pow(ds[0],2);
    else
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C);

    for(i=1;i<N-1;i++)
      {
	if(adap)
	  b[i] = (z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h*(1.-tau)+
	    tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
		 (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
	else
	  b[i] = (z[i-1]-2.*z[i]+z[i+1])/h/h*(1.-tau)+
	    tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
		 (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
	b[i] = b[i]*h*h;
	if(W>=5.724 && it>=15000 && it%10==0 && L>1.4) 
	  {
	    //OutPut(i<<' '<<(z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])*(1.-tau)<<
	    //   ' '<<tau*h*(((r[i]-r[i-1])*f[i]+(r[i+1]-r[i])*f[i+1])/2.)<<' '<< 
	    //	   tau*h*(r[i+1]-r[i-1])/(2.)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C)<<endl);
	    //OutPut(i<<' '<<(r[i+1]-r[i-1])/(2.)<<' '<<(z[i+1]-z[i-1])/(2.*h)<<' '<<
	    //	   r[i]<<' '<<ds[i]<<endl);
	  }
      }
    b[N-1] = 0.;
    
    //if(W>=5.724 && it>=15000 && it%10==0 && L>1.4) 
    /*
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("for z"<<endl);
	for(i=0;i<=N-1;i++)
	  {
	    OutPut(c[i]<<' '<<d[i]<<' '<<e[i]<<' '<<b[i]<<endl);
	  }
      }
    */
    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);

// for r
    c[N-1] = -1.;
    d[0] = 1.;
    d[N-1] = 1.;
    e[0] = 0;
// b
    b[0] = 0;
    for(i=1;i<N-1;i++)
      {
	if(adap)
	  b[i] = (r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h*(1.-tau)-
	    tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		 (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
	else
	  b[i] = (r[i-1]-2.*r[i]+r[i+1])/h/h*(1.-tau)-
	    tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		 (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
	b[i] = b[i]*h*h;
      }

    if(adap)    
      b[N-1] = -h*h/2.*(f[N] + K/r[N-1] + C)*pow(ds[N-1],2);
    else
      b[N-1] = -h*h/2.*(f[N] + K/r[N-1] + C);
    
    //if(W>=5.724 && it>=15000 && it%10==0 && L>1.4) 
    /*
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("for r"<<endl);
	for(i=0;i<=N-1;i++)
	  {
	    OutPut(c[i]<<' '<<d[i]<<' '<<e[i]<<' '<<b[i]<<endl);
	  }
      }
    */
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    //if(W>=5.724 && it>=15000 && it%10==0 && L>1.4) 
    /*
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("Sol "<<it<<endl);    
	for(i=0;i<N;i++)
	  {
	    OutPut(r[i]*L<<' '<<z[i]*L<<' '<<f[i]+C<<endl);
	  }
	OutPut("L="<<L<<" C="<<C<<endl);
      }
    */
    if(it>100000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
  }
  //OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);
  
  L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);

    //if(W>=5.724) OutPut("Surface"<<endl);
    for(i=0;i<N;i++)
      {
	r[i]*=L; z[i]*=L;
	//if(W>=5.724) OutPut(r[i]<<' '<<z[i]<<endl);
      }
  
    /*
    if(W>=5.724)
      {
	OutPut("F"<<endl);
	for(i=0;i<N+1;i++)
	  {
	    if(i==0) OutPut(z[i]<<' '<<f[i]+C<<endl);
	    if(i>0 && i!=N) OutPut((z[i-1]+z[i])/2.<<' '<<f[i]+C<<endl);
	    if(i==N) OutPut(z[N-1]<<' '<<f[N]+C<<endl);
	  }
	OutPut("C="<<C<<' '<<"C/L="<<C/L<<endl);
      }
    */
  *LL = L;
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(2.*L),h); 
      //OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }

  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  if(adap) {delete t; delete s; delete ds; delete dsm;}
  return ID;
}
// the information about field only in the middle of edges and at last points;
// a criterion for interrupting the iterative process is the norm for the difference 
// between an initial surface (r_init,z_init) and a current calculated shape (r,z)
int SchemeA_ax11(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a)
{
  double h = 1./(N-1), eps = 5.e-2;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp, *r_init, *z_init;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case

  // for adaptation
  double *t; //uniform grid
  double *s; //nonuniform grid
  double *ds, *dsm; // derivatives from s(a,t) on t at vertices and at the midpoints
  if(adap)
    {    
      t = new double[N]; s = new double[N]; ds = new double[N]; dsm = new double[N];
      for(i=0;i<N;i++)
	{
	  t[i] = i*h;
	  s[i] = S(*a,t[i]);
	  ds[i] = dSdt(*a,t[i]);
	  if(i!=0) dsm[i] = dSdt(*a,t[i]-h/2.);
	}
    }
  // ---------------
  f = new double[N+1]; f1 = new double[N+1]; f2 = new double[N+1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];
  r_init = new double[N]; z_init = new double[N];
// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3
  
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
    r_init[i] = r[i]; z_init[i] = z[i];
  }
  
  it = 0;
  err = 0.;
  while((err<eps)&&(it<10000)) 
  {
    it++;
   
     L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);
    //OutPut("L "<<it<<' '<<L<<endl);

    for(i=0;i<N+1;i++)
    {
      Hn = h_n[i];
      H = sqrt(h_s[i]*h_s[i] + h_n[i]*h_n[i]);
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      
      // for nonlinear case      
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    } 
   
    //if(W>=5.724 && it==1) OutPut("F "<<it<<endl);
    for(i=0;i<N+1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i]; 
	//if(W>=5.724 && it==1 && i==0) OutPut(z[i]<<' '<<f[i]+C<<endl);
	//if(W>=5.724 && it==1 && i>0 && i!=N) OutPut((z[i-1]+z[i])/2.<<' '<<f[i]+C<<endl);
	//if(W>=5.724 && it==1 && i==N) OutPut(z[N-1]<<' '<<f[N]+C<<endl);
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
//  ----------------------------------------------  
      
    C = 0.;
    if(adap)
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + C/pow(r[N-1],K));
      }
    else
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));
      }
    //OutPut("C/L "<<C/L<<endl);

    err = 0.; 
    for(i=0;i<N;i++)
    {
      err+=pow(fabs(r[i]-r_init[i]),2);
    }
    err_r = sqrt(err);

    err = 0.; 
    for(i=0;i<N;i++)
    {
      err+=pow(fabs(z[i]-z_init[i]),2);
    }
    err_z = sqrt(err);
    err = sqrt(err_r*err_r+err_z*err_z);
    
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      if(adap)
	c[i] = 1./dsm[i];
      else
	c[i] = 1.;
 // d
    for(i=1;i<N-1;i++)
      if(adap)
	d[i] = -1.*(1./dsm[i]+1./dsm[i+1]);
      else
	d[i] = -2.;
 // e
    for(i=1;i<N-1;i++)
      if(adap)
	e[i] = 1./dsm[i+1];
      else
	e[i] = 1.;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1.;
    d[N-1] = 1.;
    e[0] = 1.;
// b
    if(adap)
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C)*pow(ds[0],2);
    else
      b[0] = h*h/2.*(f[0] - K*(f[0]+C)/2. + C);

    for(i=1;i<N-1;i++)
      {
	if(adap)
	  b[i] = (z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h*(1.-tau)+
	    tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
		 (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
	else
	  b[i] = (z[i-1]-2.*z[i]+z[i+1])/h/h*(1.-tau)+
	    tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
		 (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
	b[i] = b[i]*h*h;
      }
    b[N-1] = 0.;
    
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("for z"<<endl);
	for(i=0;i<=N-1;i++)
	  {
	    OutPut(c[i]<<' '<<d[i]<<' '<<e[i]<<' '<<b[i]<<endl);
	  }
      }
    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);

// for r
    c[N-1] = -1.;
    d[0] = 1.;
    d[N-1] = 1.;
    e[0] = 0;
// b
    b[0] = 0;
    for(i=1;i<N-1;i++)
      {
	if(adap)
	  b[i] = (r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h*(1.-tau)-
	    tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		 (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
	else
	  b[i] = (r[i-1]-2.*r[i]+r[i+1])/h/h*(1.-tau)-
	    tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		 (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
	b[i] = b[i]*h*h;
      }

    if(adap)    
      b[N-1] = -h*h/2.*(f[N] + K/r[N-1] + C)*pow(ds[N-1],2);
    else
      b[N-1] = -h*h/2.*(f[N] + K/r[N-1] + C);
    
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("for r"<<endl);
	for(i=0;i<=N-1;i++)
	  {
	    OutPut(c[i]<<' '<<d[i]<<' '<<e[i]<<' '<<b[i]<<endl);
	  }
      }
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    if(W>=5.724 && it%10==0 && L>1.7) 
      {
	OutPut("Sol "<<it<<endl);    
	for(i=0;i<N;i++)
	  {
	    OutPut(r[i]*L<<' '<<z[i]*L<<' '<<f[i]+C<<endl);
	  }
	OutPut("L="<<L<<" C="<<C<<endl);
      }
    //if(it>100000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
  }
  
  L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);

    if(W>=5.724) OutPut("Surface"<<endl);
    for(i=0;i<N;i++)
      {
	r[i]*=L; z[i]*=L;
	if(W>=5.724) OutPut(r[i]<<' '<<z[i]<<endl);
      }
  
    if(W>=5.724)
      {
	OutPut("F"<<endl);
	for(i=0;i<N+1;i++)
	  {
	    if(i==0) OutPut(z[i]<<' '<<f[i]+C<<endl);
	    if(i>0 && i!=N) OutPut((z[i-1]+z[i])/2.<<' '<<f[i]+C<<endl);
	    if(i==N) OutPut(z[N-1]<<' '<<f[N]+C<<endl);
	  }
	OutPut("C="<<C<<' '<<"C/L="<<C/L<<endl);
      }
    OutPut("Number of surface iterations: "<<it<<endl);

  *LL = L;
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(2.*L),h); 
      //OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }

  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  delete r_init; delete z_init;
  if(adap) {delete t; delete s; delete ds; delete dsm;}
  return ID;
}
// the information about field only in the middle of edges and at last points
int SchemeA_ax2(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a)
{
  double h = 1./(N-1), eps = 1.e-8;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it, i_err_r, i_err_z;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case

  // for adaptation
  double *t; //uniform grid
  double *s; //nonuniform grid
  double *ds, *dsm; // derivatives from s(a,t) on t at vertices and at the midpoints
  if(adap)
    {    
      t = new double[N]; s = new double[N]; ds = new double[N]; dsm = new double[N];
      for(i=0;i<N;i++)
	{
	  t[i] = i*h;
	  s[i] = S(*a,t[i]);
	  ds[i] = dSdt(*a,t[i]);
	  if(i!=0) dsm[i] = dSdt(*a,t[i]-h/2.);
	}
    }
  // ---------------
  f = new double[N+1]; f1 = new double[N+1]; f2 = new double[N+1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];
    
// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3
  
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
  
  it = 0;
  err_z = 0;
  err_r = 0;
  while(it<5 || err_r>eps || err_z>eps)
  {
    it++;
       
     L = 0.;
    for(i=1;i<N;i++)
      L+=pow((r[i-1]+r[i])/2.,2)*(z[i]-z[i-1])/h;
    L = pow(-2*Pi*h*L,-1./3.);
    //OutPut("L "<<it<<' '<<L<<endl);

   
    for(i=0;i<N+1;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==N)
      {
        Hn = 0;
        H = h_z[N];
      }
      else 
	{  
	  if(adap)
	    Hn = -(z[i]-z[i-1])/(s[i]-s[i-1])*h_r[i] + 
	      (r[i]-r[i-1])/(s[i]-s[i-1])*h_z[i];
	  else
	    Hn = -(z[i]-z[i-1])/h*h_r[i] + 
	      (r[i]-r[i-1])/h*h_z[i];
	}
	  H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      
      // for nonlinear case      
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    }
    
    if(it%500==0) OutPut("F"<<endl);
    for(i=0;i<N+1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i]; 
	if(it%500==0 && i==0) OutPut(z[i]<<' '<<f[i]+C<<endl);
	if(it%500==0 && i>0 && i!=N) OutPut((z[i-1]+z[i])/2.<<' '<<f[i]+C<<endl);
	if(it%500==0 && i==N) OutPut(z[N-1]<<' '<<f[N]+C<<endl);
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
//  ----------------------------------------------  
      
    C = 0.;
    if(adap)
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + C/pow(r[N-1],K));
      }
    else
      {
	for(i=1;i<N;i++)
	  C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[i];
	C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));
      }
    //if(it%500==0) OutPut("C/L "<<C/L<<' '<<C<<' '<<L<<endl);

    err = 0.; err_r = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h+
		  ((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. + 
		  (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((r[i-1]-2.*r[i]+r[i+1])/h/h+
		  ((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
		  (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_r)
      { 
        err_r = err;
        i_err_r = i;
      }
    }
    
    err = 0.; err_z = 0.;
    for(i=1;i<N-1;i++)
    {
      if(adap)
	err = fabs((z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h-
		  ((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. -
		  (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	err = fabs((z[i-1]-2.*z[i]+z[i+1])/h/h-
		  ((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. -
		  (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_z)
      { 
        err_z = err;
        i_err_z = i;
      }
    }
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      if(adap)
	c[i] = 1./dsm[i]/h/h;
      else
	c[i] = 1./h/h;
 // d
    for(i=1;i<N-1;i++)
      if(adap)
	d[i] = -1.*(1./dsm[i]+1./dsm[i+1])/h/h;
      else
	d[i] = -2./h/h;
 // e
    for(i=1;i<N-1;i++)
      if(adap)
	e[i] = 1./dsm[i+1]/h/h;
      else
	e[i] = 1./h/h;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1./h;
    d[N-1] = 1.;
    e[0] = 1./h;
// b
    if(adap)
      b[0] = h/2.*(f[0] - K*(f[0]+C)/2. + C)*pow(ds[0],2);
    else
      b[0] = h/2.*(f[0] - K*(f[0]+C)/2. + C);

    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (z[i-1]/dsm[i]-z[i]*(1./dsm[i]+1./dsm[i+1])+z[i+1]/dsm[i+1])/h/h*(1.-tau)+
	  tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
	       (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	b[i] = (z[i-1]-2.*z[i]+z[i+1])/h/h*(1.-tau)+
	  tau*(((r[i]-r[i-1])/h*f[i]+(r[i+1]-r[i])/h*f[i+1])/2. +
	       (r[i+1]-r[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
    b[N-1] = 0.;
    
    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);

// for r
    c[N-1] = -1./h;
    d[0] = 1.;
    d[N-1] = 1./h;
    e[0] = 0;
// b
    b[0] = 0;
    for(i=1;i<N-1;i++)
      if(adap)
	b[i] = (r[i-1]/dsm[i]-r[i]*(1./dsm[i]+1./dsm[i+1])+r[i+1]/dsm[i+1])/h/h*(1.-tau)-
	  tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
	       (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i]/ds[i] + C));
      else
	b[i] = (r[i-1]-2.*r[i]+r[i+1])/h/h*(1.-tau)-
	  tau*(((z[i]-z[i-1])/h*f[i]+(z[i+1]-z[i])/h*f[i+1])/2. +
	       (z[i+1]-z[i-1])/(2.*h)*(-K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));

    if(adap)    
      b[N-1] = -h/2.*(f[N] + K/r[N-1] + C)*pow(ds[N-1],2);
    else
      b[N-1] = -h/2.*(f[N] + K/r[N-1] + C);
    
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    if(it>100000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
  }
  OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);
  
  for(i=0;i<N;i++)
  {
    r[i]*=L; z[i]*=L;
  }
  
  *LL = L;
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(2.*L),h); 
      OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }

  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  if(adap) {delete t; delete s; delete ds; delete dsm;}
  return ID;
}
int SchemeB(double *r, double *z, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL)
{
  double h = 1./(N-1), eps = 1.e-6;
  double *f, *f1, *f2, *c, *d, *e, *b, *temp;
  double C, L=*LL, H, Hn;
  double err, err_r, err_z;
  int ID=1, i, it, i_err_r, i_err_z;
  int K=0; // K=1 for axisymmetrical case, K=0 for plane case
  
  f = new double[2*N-1]; f1 = new double[2*N-1]; f2 = new double[2*N-1];
  c = new double[N]; d = new double[N]; e = new double[N];
  b = new double[N]; temp = new double[N];

// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3  

  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }

  it = 0;
  err_r = 0;
  err_z = 0;
  while(it<5 || err_r>eps || err_z>eps)
  {
    it++;

    L = 0.;
    for(i=1;i<N;i++)
      L+=(z[i-1]+z[i])/2.*(r[i]-r[i-1])/h;
    L = pow(h*L,-1./2.)/2.;
    
  // F = f - K*z'/r + C - for axisymmetrical case,
  // f = -2WL/3/hi*f1-WL*f2
  // ----------------------------------------------
  // f[i] the value on the free surface:
  // if i=2j - at the point, i=2j+1 - in the middle
  
  // for the case of magnetic saturation 
/*
    memset(f1,0,(2*N-1)*SizeOfDouble);
    f2[0] = 1.; f2[2*N-2] = 0.;
    for(i=1;i<2*N-2;i++)
      if(i%2==0)
        f2[i] = pow((r[i/2+1]-r[i/2-1])/(2.*h),2);
      else
        f2[i] = pow((r[i/2+1]-r[i/2])/h,2);
*/
  // taking into considaration the magnetic field
    for(i=0;i<2*N-1;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==2*N-2)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else if(i%2==0)
        {
        Hn = -(z[i/2+1]-z[i/2-1])/(2.*h)*h_r[i/2] + 
              (r[i/2+1]-r[i/2-1])/(2.*h)*h_z[i/2];
        H = sqrt(h_r[i/2]*h_r[i/2] + h_z[i/2]*h_z[i/2]);
        }
        else
        {
        Hn = -(z[i/2+1]-z[i/2])/h*(h_r[i/2+1]+h_r[i/2])/2. + 
              (r[i/2+1]-r[i/2])/h*(h_z[i/2+1]+h_z[i/2])/2.;
        H = sqrt(pow((h_r[i/2+1]+h_r[i/2])/2.,2) + 
                 pow((h_z[i/2+1]+h_z[i/2])/2.,2));
        }
      // for linear case
      f1[i] = H*H;
      f2[i] = Hn*Hn;
      // for nonlinear case       
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    }          

    for(i=0;i<2*N-1;i++)
      {
	// for linear case
	f[i] = -W*L*f1[i]-W*L*hi*f2[i];
	// for nonlinear case
	// f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      }
 //  ----------------------------------------------  
      
    C = 0.;
    for(i=1;i<N;i++)
      C+=pow((r[i-1]+r[i])/2.,K)*(r[i]-r[i-1])/h*f[2*i-1];
    C = -pow(2.,K)/r[N-1]*(1. + h*C/pow(r[N-1],K));

    err_r = 0.;
    for(i=1;i<N-1;i++)
    {
      err = fabs((r[i-1]-2.*r[i]+r[i+1])/h/h+
                (z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_r)
      { 
        err_r = err;
        i_err_r = i;
      }
    }
    
    err_z = 0.;
    for(i=1;i<N-1;i++)
    {
      err = fabs((z[i-1]-2.*z[i]+z[i+1])/h/h-
                (r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C));
      if(err>err_z)
      { 
        err_z = err;
        i_err_z = i;
      }
    }
// the common part for r and z    
 // c
    c[0] = 0.;  // must be 0
    for(i=1;i<N-1;i++)
      c[i] = -1./h/h;
 // d
    for(i=1;i<N-1;i++)
      d[i] = 2./h/h + 1./tau;
 // e
    for(i=1;i<N-1;i++)
      e[i] = -1./h/h;
    e[N-1] = 0.;  // must be 0

// for z
    c[N-1] = 0.;
    d[0] = -1./h;
    d[N-1] = 1.;
    e[0] = 1./h;
// b
    b[0] = h/2.*(f[0] - K*(f[0]+C)/2. + C);
    for(i=1;i<N-1;i++)
      b[i] = z[i]/tau-
             (r[i+1]-r[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);
    b[N-1] = 0.;
    
    Solver_3diag(N,c,d,e,b);
    memcpy(temp,b,N*SizeOfDouble);

// for r
    c[N-1] = -1./h;
    d[0] = 1.;
    d[N-1] = 1./h;
    e[0] = 0;
// b
    b[0] = 0;
    for(i=1;i<N-1;i++)
      b[i] = r[i]/tau+
             (z[i+1]-z[i-1])/(2.*h)*(f[2*i] - K*(z[i+1]-z[i-1])/(2.*h)/r[i] + C);
    b[N-1] = -h/2.*(f[2*N-2] + K/r[N-1] + C);
    
    Solver_3diag(N,c,d,e,b);
    memcpy(r,b,N*SizeOfDouble);
    memcpy(z,temp,N*SizeOfDouble);
    
    if(it>10000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
    
//    for(i=0;i<N;i++)
//      OutPut(r[i]<<' '<<z[i]<<endl);
  }
  OutPut("It = "<<it<<", err_r = "<<err_r<<", err_z = "<<err_z<<endl);

  for(i=0;i<N;i++)
  {
    r[i]*=L; z[i]*=L;
  }

  *LL = L;
  delete f; delete f1; delete f2; 
  delete c; delete d; delete e; delete b; delete temp;
  return ID;
}

int SchemeT4(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL, 
        double *F, double *dF)
{
  double h = 1./(N-1), eps = 1.e-6;
  double C, L=*LL, H, Hn;
  double *f, *f1, *f2, *df, *df1, *df2;
  double *beta_old, *r_old, *z_old;
  double err, err_beta;
  int ID=1, i, it, i_err_beta;
  int K=0; // K=1 for axisymmetrical case, K=0 for plane case 
  
  f = new double[N]; f1 = new double[N]; f2 = new double[N];
  df = new double[N]; df1 = new double[N]; df2 = new double[N]; 
  beta_old = new double[N]; r_old = new double[N]; z_old = new double[N];

// it is necessary, when the input space coordinates r and z
// are dimentionless on D^1/2  
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
  
  it = 0;
  err_beta = 0;
  while(it<5 || err_beta>eps)
  {
    it++;
    
    memcpy(beta_old,beta,N*SizeOfDouble);
    memcpy(r_old,r,N*SizeOfDouble);
    memcpy(z_old,z,N*SizeOfDouble);
              
    L = 0.;
    for(i=1;i<N-1;i++)
      L+=cos(beta_old[i])*z[i];
    L = pow(h*L+h*z[0]/2.,-1./2.)/2.;
    
  // F = f - K*z'/r + C 
  // f = -2WL/3/hi*f1-WL*f2
  // ----------------------------------------------
  
  // for the case of magnetic saturation 
    /*
    memset(f1,0,N*SizeOfDouble);
    f2[0] = 1.; f2[N-1] = 0.;
    for(i=1;i<N-1;i++)
      f2[i] = pow(cos(beta_old[i]),2);
    */

  // taking into considaration the magnetic field
    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      }
      // for linear case
       f1[i] = H*H;
       f2[i] = Hn*Hn;
      // for nonlinear case    
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    }          
    for(i=0;i<N;i++)
    {
      // for linear case
       f[i] = -W*L*f1[i]-W*L*hi*f2[i];
      // for nonlinear case
      // f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
    }
 //  ----------------------------------------------  
 
    C = 0.;
    for(i=1;i<N-1;i++)
      C+=h*pow(r[i],K)*cos(beta_old[i])*f[i];
    C = -pow(2.,K)/r[N-1]*(1. + 
        (C+h*pow(r[0],K)*f[0]/2.-
         h*h/12.*(pow(r[N-1],K)*f[N-1]*F[N-1]-
         K*f[0]))/pow(r[N-1],K));
    

    F[0] = f[0] - K*(f[0]+C)/2. + C;
    for(i=1;i<N;i++)
      F[i] = f[i] - K*sin(beta_old[i])/r[i] + C;
      
  // df = -2WL/3/hi*df1-WL*df2
  // ----------------------------------------------
    memset(df1,0,N*SizeOfDouble); // f1 not depend on s

  // for the case of magnetic saturation 
    /*
    df2[0] = 0.;
    for(i=1;i<N;i++)
      df2[i] = -sin(2.*beta_old[i])*F[i];
    */
  // taking into considaration the magnetic field
    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      }
      
      // for linear case
      df2[i] = 2.*Hn*F[i]*
	(-cos(beta_old[i])*h_r[i] -
	 sin(beta_old[i])*h_z[i]);    
      // for nonlinear case
      // arg = H*gamma;
      // df2[i] = pow((1./tanh(arg)-1./arg)/H,2)*2.*Hn*F[i]*
      //         (-cos(beta_old[i])*h_r[i] -
      //          sin(beta_old[i])*h_z[i]);    
    }          

    for(i=0;i<N;i++)
    {
      // for linear case
      df[i] = -W*L*df1[i]-W*L*hi*df2[i];
      // for nonlinear case
      // df[i] = -2.*W*L/3./hi*df1[i]-W*L*df2[i];    
    }
 //  ----------------------------------------------  
 
    dF[0] = 0.;
    for(i=1;i<N;i++)
      dF[i] = df[i] - K*(cos(beta_old[i])*F[i]*r[i]-
              sin(2.*beta_old[i])/2.)/r[i]/r[i];
              
    err_beta = 0.;
    for(i=1;i<N-1;i++)
    {
      err = fabs((beta_old[i+1]-beta_old[i])/h-
                (F[i+1]+F[i])/2. + h/12.*(dF[i+1]-dF[i]));
      if(err>err_beta)
      {
        err_beta = err;
        i_err_beta = i;
      }
    }
    
  // beta
    beta[0] = 0.;
    beta[N-1] = -Pi/2.;
    for(i=N-2;i>0;i--)
      beta[i] = beta[i+1] - h*((F[i]+F[i+1])/2.-h*(dF[i+1]-dF[i])/12.)+
      (1.-tau)*(beta_old[i]-beta_old[i+1]+h*((F[i]+F[i+1])/2.-h*(dF[i+1]-dF[i])/12.));
      
  // r
    r[0] = 0.;
    for(i=1;i<N;i++)
      r[i] = r[i-1] + h*((cos(beta[i-1])+cos(beta[i]))/2.+
             h/12.*(F[i]*sin(beta[i])-F[i-1]*sin(beta[i-1])));

  // z
    z[N-1] = 0.;
    for(i=N-2;i>=0;i--)
      z[i] = z[i+1] - h*((sin(beta[i])+sin(beta[i+1]))/2.-
             h/12.*(F[i+1]*cos(beta[i+1])-F[i]*cos(beta[i])));
 
    if(it>10000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
  }
 
  for(i=0;i<N;i++)
  {
    r[i]*=L; z[i]*=L;
    OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<endl);
  }
 
  *LL = L;
  delete f; delete f1; delete f2; 
  delete df; delete df1; delete df2;
  delete beta_old; delete r_old; delete z_old;

  OutPut("It = "<<it<<", err_beta = "<<err_beta<<endl);

  return ID;
}
int SchemeT4_ax(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL, 
        double *F, double *dF, int adap, double* a, int up)
{
  double hh = 1./(N-1), eps = 1.e-6;
  double C, L=*LL, H, Hn;
  double *f, *f1, *f2, *df, *df1, *df2;
  double *beta_old, *r_old, *z_old;
  double err, err_beta;
  int ID=1, i, it, i_err_beta;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case 

  f = new double[N]; f1 = new double[N]; f2 = new double[N];
  df = new double[N]; df1 = new double[N]; df2 = new double[N]; 
  beta_old = new double[N]; r_old = new double[N]; z_old = new double[N];
  // for adaptation
  double *t; //uniform grid
  double *h; //h[i]=s[i]-s[i-1], s - nonuniform grid

  if(!adap) *a = 100;

  t = new double[N]; h = new double[N];
  for(i=0;i<N;i++)
    {
      t[i] = i*hh;
      if(i!=0) h[i] = S(*a,t[i])-S(*a,t[i-1]);
    }

// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3  
    
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
  
  it = 0;
  err_beta = 0;
  while(it<5 || err_beta>eps) 
  {
    it++;
    
    memcpy(beta_old,beta,N*SizeOfDouble);
    memcpy(r_old,r,N*SizeOfDouble);
    memcpy(z_old,z,N*SizeOfDouble);
              
    L = 0.;
    for(i=1;i<N-1;i++)
      L+=(h[i]+h[i+1])/2.*sin(beta_old[i])*r[i]*r[i];
	//-(h[i]*h[i]-h[i+1]*h[i+1])*(cos(beta_old[i])*r[i]*r[i]*F[i]+
	//r[i]*sin(2.*beta_old[i]))/12.;
    L = pow(-2.*Pi*(L-h[N-1]*r[N-1]*r[N-1]/2.),-1./3.);
    
    if(it%500==0) OutPut("it="<<it<<endl<<"L="<<L<<endl);
  // F = f - K*z'/r + C 
  // ----------------------------------------------
  
  // for the case of magnetic saturation 
    /*            
    memset(f1,0,N*SizeOfDouble);
    f2[0] = 1.; f2[N-1] = 0.;
    for(i=1;i<N-1;i++)
      f2[i] = pow(cos(beta_old[i]),2);
    */
    
  // taking into considaration the magnetic field
    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      }
      // for linear case
       f1[i] = H*H;
       f2[i] = Hn*Hn;
      // for nonlinear case    
      // arg = H*gamma;
      // f1[i] = log(sinh(arg)/arg);    
      // f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
    } 
        
    //if(it%500==0) 
    //OutPut("f"<<endl);
    for(i=0;i<N;i++)
    {
      // for linear case
      f[i] = -W*L*f1[i]-W*L*hi*f2[i]; 
      //if(it%500==0) OutPut(z[i]<<' '<<-W*L*f1[i]<<' '<<-W*L*hi*f2[i]<<' '<<f[i]<<endl);
      // for nonlinear case
      // f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
    }
 //  ----------------------------------------------  
    
    C = 0.;
    for(i=1;i<N-1;i++)
      C+=(h[i]+h[i+1])/2.*pow(r[i],K)*cos(beta_old[i])*f[i];
	//-(h[i]*h[i]-h[i+1]*h[i+1])/12.*(K*cos(beta_old[i])*cos(beta_old[i])*f[i]-
	//pow(r[i],K)*sin(beta_old[i])*F[i]*f[i]+
	//pow(r[i],K)*cos(beta_old[i])*df[i]);
    C = -pow(2.,K)/r[N-1]*(1. + 
			   (C+h[1]/2.*pow(r[0],K)*f[0]
				//-1./12.*(h[N-1]*h[N-1]*pow(r[N-1],K)*f[N-1]*F[N-1]-h[1]*h[1]*K*f[0])
			    )/pow(r[N-1],K));    
    

    if(it%500==0) OutPut("C/L="<<C/L<<endl);

    F[0] = f[0] - K*(f[0]+C)/2. + C; 
    for(i=1;i<N;i++)
      {
	F[i] = f[i] - K*sin(beta_old[i])/r[i] + C;
      }
    
    if(it%500==0) OutPut("F"<<endl);
    for(i=0;i<N;i++)      
      if(it%500==0)
      if(i==0) 
	{OutPut(r[i]<<' '<<z[i]<<' '<<f[i]<<' '<<(f[0]+C)/2.<<' '<<F[i]<<endl);}
      else
	{OutPut(r[i]<<' '<<z[i]<<' '<<f[i]<<' '<<sin(beta_old[i])/r[i]<<' '<<F[i]<<endl);}
      

    if(it%500==0) OutPut("F and beta'"<<endl);
    for(i=0;i<N;i++)      
      if(it%500==0)
      if(i==0) 
	{OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<' '<<(beta_old[1]-beta_old[0])/h[1]<<endl);}
      else
	if(i==N-1)
	  {OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<' '<<
		(beta_old[N-1]-beta_old[N-2])/(h[N-1])<<endl);}
	else
	  {OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<' '<<
		  (beta_old[i+1]-beta_old[i-1])/(h[i+1]+h[i])<<endl);}
    /*
    F[0] = F[4]+(F[3]-F[4])/(S(*a,3)-S(*a,4))*(S(*a,0)-S(*a,4));
    F[1] = F[4]+(F[3]-F[4])/(S(*a,3)-S(*a,4))*(S(*a,1)-S(*a,4));
    F[2] = F[4]+(F[3]-F[4])/(S(*a,3)-S(*a,4))*(S(*a,2)-S(*a,4));
    
    if(it%500==0) OutPut("F"<<endl);
    for(i=0;i<N;i++)      
      if(it%500==0)
      if(i==0) 
	{OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<endl);}
      else
	{OutPut(r[i]<<' '<<z[i]<<' '<<F[i]<<endl);}
    */
    if(it%500==0) OutPut("beta(s)"<<endl);
    for(i=0;i<N;i++)      
      if(it%500==0)
      if(i==0) 
	{OutPut(beta_old[i]<<' '<<S(*a,t[i])<<endl);}
      else
	{OutPut(beta_old[i]<<' '<<S(*a,t[i])<<endl);}
  // df = -2WL/3/hi*df1-WL*df2
  // ----------------------------------------------
    memset(df1,0,N*SizeOfDouble); // f1 not depend on s

  // for the case of magnetic saturation 
    /*           
    df2[0] = 0.;
    for(i=1;i<N;i++)
      df2[i] = -sin(2.*beta_old[i])*F[i];
    */
    
  // taking into considaration the magnetic field
    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H =  Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      }
          
      // for linear case
      df2[i] = 2.*Hn*F[i]*
	(-cos(beta_old[i])*h_r[i] -
	 sin(beta_old[i])*h_z[i]);    

      df2[0] = 0.;
      // for nonlinear case
      // arg = H*gamma;
      // df2[i] = pow((1./tanh(arg)-1./arg)/H,2)*2.*Hn*F[i]*
      //         (-cos(beta_old[i])*h_r[i] -
      //          sin(beta_old[i])*h_z[i]);    
    }          
    
    //if(it%500==0) OutPut("df"<<endl);
    for(i=0;i<N;i++)
    {
      // for linear case
      df[i] = -W*L*df1[i]-W*L*hi*df2[i]; 
      //if(it%500==0) OutPut(z[i]<<' '<<df[i]<<endl);
      // for nonlinear case
      // df[i] = -2.*W*L/3./hi*df1[i]-W*L*df2[i];    
    }
 //  ----------------------------------------------  
 
    dF[0] = 0.;
    for(i=1;i<N;i++)
      dF[i] = df[i] - K*(cos(beta_old[i])*F[i]*r[i]-
              sin(2.*beta_old[i])/2.)/r[i]/r[i];
              
    err_beta = 0.;
    for(i=1;i<N-1;i++)
    {
      if(up)
	{
	  err = fabs((beta_old[i+1]-beta_old[i])/h[i+1]-
		    (F[i+1]+F[i])/2. 
		    //+ h[i+1]/12.*(dF[i+1]-dF[i])
		    );
	}
      else
	err = fabs((beta_old[i]-beta_old[i-1])/h[i]-
		  (F[i-1]+F[i])/2. 
		  //+ h[i]/12.*(dF[i]-dF[i-1])
		  );

      if(err>err_beta)
      {
        err_beta = err;
        i_err_beta = i;
      }
    }
   
    //if(up==1) up = 0; else up = 1; 
  // beta
    beta[0] = 0.;
    beta[N-1] = -Pi/2.;
    if(up)
      {
	for(i=N-2;i>=0;i--)
	  beta[i] = beta[i+1] - h[i+1]*((F[i]+F[i+1])/2.
					//-h[i+1]*(dF[i+1]-dF[i])/12.
					)+
	    (1.-tau)*(beta_old[i]-beta_old[i+1]+h[i+1]*((F[i]+F[i+1])/2.
							//-h[i+1]*(dF[i+1]-dF[i])/12.
							));
      }
    else
      for(i=1;i<=N-1;i++)
	beta[i] = beta[i-1] + h[i]*((F[i-1]+F[i])/2.
				    //-h[i]*(dF[i]-dF[i-1])/12.
				    )+
	  (1.-tau)*(beta_old[i]-beta_old[i-1]-h[i]*((F[i-1]+F[i])/2.
						    //-h[i]*(dF[i]-dF[i-1])/12.
						    ));             
        
    /*   
    if(it==1||it==2||it%100==0||it%100==1)
      { 
	OutPut("it="<<it<<endl);
	for(i=0;i<N;i++)
	  OutPut(beta[i]<<endl);
        OutPut(C<<endl);
      }
    */

    beta[0] = 0.;
  // r
    r[0] = 0.;
    for(i=1;i<N;i++)
      r[i] = r[i-1] + h[i]*((cos(beta[i-1])+cos(beta[i]))/2.
			    //+h[i]/12.*(F[i]*sin(beta[i])-F[i-1]*sin(beta[i-1]))
			    );

  // z
    z[N-1] = 0.;
    for(i=N-2;i>=0;i--)
      z[i] = z[i+1] - h[i+1]*((sin(beta[i])+sin(beta[i+1]))/2.
			      //-h[i+1]/12.*(F[i+1]*cos(beta[i+1])-F[i]*cos(beta[i]))
			      );
 
    //if(it>1000000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
    //if(it%100==0) OutPut(it<<' '<<"curv = "<<-(f[0]+C)/(2.*L)<<endl);
    if(it%1000==0) OutPut("It = "<<it<<", err_beta = "<<err_beta<<endl);
  }

  OutPut("All iterations: "<<it<<endl);
  OutPut("C="<<C<<endl);
  //OutPut("beta'[0]-F[0] = "<<(beta[1]-beta[0])/h[1]-F[0]<<endl);
  /*
  for(i=0;i<N;i++)
  {
    if(i!=0)
      OutPut(r[i]<<' '<<z[i]<<' '<<beta[i]<<' '<<F[i]<<' '<<sin(beta[i])/r[i]<<' '<<f[i]+C<<endl)
    else    
      OutPut(r[i]<<' '<<z[i]<<' '<<beta[i]<<' '<<F[i]<<' '<<(f[0]+C)/2.<<' '<<f[i]+C<<endl);
  }
  */
  for(i=0;i<N;i++)
    {
      r[i]*=L; z[i]*=L;
    }
  *LL = L;
  OutPut("Err "<<err_beta<<endl);

  if (adap) 
    {  
      *a = A(-(f[0]+C)/(2.*L),hh); 
      OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl); 
    }
  /*  
  OutPut("Curv = "<<endl);
  for(i=0;i<N;i++)
    {
      OutPut(i<<' '<<-(f[i]+C)/(2.*L)<<endl);
    }
  */
  delete f; delete f1; delete f2; 
  delete df; delete df1; delete df2;
  delete beta_old; delete r_old; delete z_old;
  delete t; delete h;
  return ID;
}
int SchemeT4_axNL(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL, 
        double *F, double *dF, int adap, double* a)
{
  double hh = 1./(N-1), eps = 1.e-6;
  double C, L=*LL, arg, H, Hn;
  double *f, *f1, *f2, *df, *df1, *df2;
  double *beta_old, *r_old, *z_old;
  double err, err_beta;
  int ID=1, i, it, i_err_beta;
  int K=1; // K=1 for axisymmetrical case, K=0 for plane case 

  f = new double[N]; f1 = new double[N]; f2 = new double[N];
  df = new double[N]; df1 = new double[N]; df2 = new double[N]; 
  beta_old = new double[N]; r_old = new double[N]; z_old = new double[N];
  // for adaptation
  double *t; //uniform grid
  double *h; //h[i]=s[i]-s[i-1], s - nonuniform grid

  if(!adap) *a = 100;

  t = new double[N]; h = new double[N];
  for(i=0;i<N;i++)
    {
      t[i] = i*hh;
      if(i!=0) h[i] = S(*a,t[i])-S(*a,t[i-1]);
    }

// it is necessary, when the input space coordinates r and z
// are dimentionless on V^1/3      
  for(i=0;i<N;i++)
  {
    r[i]/=L; z[i]/=L;
  }
        
  it = 0;
  err_beta = 0;
  while(it<5 || err_beta>eps)
  {
    it++;
    
    memcpy(beta_old,beta,N*SizeOfDouble);
    memcpy(r_old,r,N*SizeOfDouble);
    memcpy(z_old,z,N*SizeOfDouble);
              
    L = 0.;
    for(i=1;i<N-1;i++)
      L+=(h[i]+h[i+1])/2.*sin(beta_old[i])*r[i]*r[i];
    //-(h[i]*h[i]-h[i+1]*h[i+1])*(cos(beta_old[i])*r[i]*r[i]*F[i]+
    //   		     r[i]*sin(2.*beta_old[i]))/12.;
    L = pow(-2.*Pi*(L-h[N-1]*r[N-1]*r[N-1]/2.),-1./3.);

    if(it%500==0) OutPut("it="<<it<<endl<<"L="<<L<<endl);

  // F = f - K*z'/r + C 
  // ----------------------------------------------  
    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H = Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      } 
      arg = H*gamma;
      if(arg>1e-6)
	f1[i] = log(sinh(arg)/arg);    
      else
	f1[i] = 0;
      if(arg>1e-2)   
	f2[i] = pow((1./tanh(arg)-1./arg)*Hn/H,2);    
      else
	f2[i] = pow(hi*Hn,2);
    } 

    if(it%500==0) 
      OutPut("f"<<endl);
    for(i=0;i<N;i++)
    {
      f[i] = -2.*W*L/3./hi*f1[i]-W*L*f2[i];    
      if(it%500==0) 
	OutPut(it<<' '<<z[i]<<' '<<-2.*W*L/3./hi*f1[i]<<' '<<-W*L*f2[i]<<' '<<f[i]<<endl);
    }
 //  ----------------------------------------------  
 
    C = 0.;
    for(i=1;i<N-1;i++)
      C+=(h[i]+h[i+1])/2.*pow(r[i],K)*cos(beta_old[i])*f[i];
    //-(h[i]*h[i]-h[i+1]*h[i+1])/12.*(K*cos(beta_old[i])*cos(beta_old[i])*f[i]-
    //pow(r[i],K)*sin(beta_old[i])*F[i]*f[i]+
    //pow(r[i],K)*cos(beta_old[i])*df[i]);
    C = -pow(2.,K)/r[N-1]*(1. + 
        (C+h[1]/2.*pow(r[0],K)*f[0]
         //-1./12.*(h[N-1]*h[N-1]*pow(r[N-1],K)*f[N-1]*F[N-1]-h[1]*h[1]*K*f[0])
	 )/pow(r[N-1],K));
    
    if(it%500==0) 
      OutPut("C/L="<<C/L<<endl);

    F[0] = f[0] - K*(f[0]+C)/2. + C;
    for(i=1;i<N;i++)
      F[i] = f[i] - K*sin(beta_old[i])/r[i] + C;

    if(it%500==0) 
      OutPut("F"<<endl);
    for(i=0;i<10;i++)
      if(it%500==0)
	if(i==0) 
	  {OutPut(z[i]<<' '<<(f[0]+C)/2.<<' '<<F[i]<<endl);}
	else
	  {OutPut(z[i]<<' '<<sin(beta_old[i])/r[i]<<' '<<F[i]<<endl);}

  // df = -2WL/3/hi*df1-WL*df2
  // ----------------------------------------------
    memset(df1,0,N*SizeOfDouble); // f1 not depend on s

    for(i=0;i<N;i++)
    {
      if(i==0)
      {
        Hn = h_z[0];
        H =  Hn;
      }
      else if(i==N-1)
      {
        Hn = 0;
        H = h_z[N-1];
      }
      else
      {
        Hn = -sin(beta_old[i])*h_r[i] + 
             cos(beta_old[i])*h_z[i];
        H = sqrt(h_r[i]*h_r[i] + h_z[i]*h_z[i]);
      }
      
      arg = H*gamma;
      if(arg>1e-2)
	df2[i] = pow((1./tanh(arg)-1./arg)/H,2)*2.*Hn*F[i]*
	  (-cos(beta_old[i])*h_r[i] -
	   sin(beta_old[i])*h_z[i]);    
      else
	df2[i] = hi*hi*2.*Hn*F[i]*
	  (-cos(beta_old[i])*h_r[i] -
	   sin(beta_old[i])*h_z[i]);    
    }          
    df2[0] = 0.;

    for(i=0;i<N;i++)
    {
      df[i] = -2.*W*L/3./hi*df1[i]-W*L*df2[i];    
    }
 //  ----------------------------------------------  
 
    dF[0] = 0.;
    for(i=1;i<N;i++)
      dF[i] = df[i] - K*(cos(beta_old[i])*F[i]*r[i]-
              sin(2.*beta_old[i])/2.)/r[i]/r[i];
              
    err_beta = 0.;
    for(i=1;i<N-1;i++)
    {
      err = fabs((beta_old[i+1]-beta_old[i])/h[i+1]-
                (F[i+1]+F[i])/2. 
		//+ h[i+1]/12.*(dF[i+1]-dF[i])
		);
      if(err>err_beta)
      {
        err_beta = err;
        i_err_beta = i;
      }
    }
    
  // beta
    beta[0] = 0.;
    beta[N-1] = -Pi/2.;
    for(i=N-2;i>0;i--)
      beta[i] = beta[i+1] - h[i+1]*((F[i]+F[i+1])/2.
				    // -h[i+1]*(dF[i+1]-dF[i])/12.
				    )+
	(1.-tau)*(beta_old[i]-beta_old[i+1]+h[i+1]*((F[i]+F[i+1])/2.
						    // -h[i+1]*(dF[i+1]-dF[i])/12.
						    ));
      
  // r
    r[0] = 0.;
    for(i=1;i<N;i++)
      r[i] = r[i-1] + h[i]*((cos(beta[i-1])+cos(beta[i]))/2.
			    //+h[i]/12.*(F[i]*sin(beta[i])-F[i-1]*sin(beta[i-1]))
			    );

  // z
    z[N-1] = 0.;
    for(i=N-2;i>=0;i--)
      z[i] = z[i+1] - h[i+1]*((sin(beta[i])+sin(beta[i+1]))/2.
			      // -h[i+1]/12.*(F[i+1]*cos(beta[i+1])-F[i]*cos(beta[i]))
			      );
 
    //if(it>1000000) {ID = 0; OutPut("Calculation was stoped"<<endl); break;}
    if(it%1000==0) OutPut("It = "<<it<<", err_beta = "<<err_beta<<endl);
  }
  
  for(i=0;i<N;i++)
  {
    r[i]*=L; z[i]*=L;
  }
  
  *LL = L;
  
  if (adap) 
    {    
      *a = A(-(f[0]+C)/(2.*L),hh);
      OutPut("curv = "<<-(f[0]+C)/(2.*L)<<" a = "<<*a<<" h(1) = "<<S(*a,t[1])-S(*a,t[0])<<endl);
    }
  
  delete f; delete f1; delete f2; 
  delete df; delete df1; delete df2;
  delete beta_old; delete r_old; delete z_old;
  delete t; delete h;
  return ID;
}
void Solver_3diag(int N, double *c, double *d, double *e, double *b)
{
  int i;
  double *aa, *bb;
  aa = new double[N]; bb = new double[N];
  
  aa[0] = -e[0]/d[0];
  bb[0] = b[0]/d[0];
  
  for(i=1;i<N-1;i++)
  {
    aa[i] = -e[i]/(d[i]+aa[i-1]*c[i]);
    bb[i] = (b[i]-bb[i-1]*c[i])/(d[i]+aa[i-1]*c[i]);
  }
  b[N-1] = (b[N-1]-bb[N-2]*c[N-1])/(d[N-1]+aa[N-2]*c[N-1]);

  for(i=N-2;i>=0;i--)
    {
    b[i]=aa[i]*b[i+1]+bb[i];
    }

  delete aa; delete bb; 
}

double S(double a, double t)
{
  return -a+2.*(a+1.)/(1.+pow(1.+2./a,1.-t));
}

double dSdt(double a, double t)
{
  return 2.*(a+1.)*pow(1.+2./a,1.-t)*log(1.+2./a)/pow(1.+pow(1.+2./a,1.-t),2);
}

double A(double curv, double h)
{
  double c, a_old, a_new, dsda;
  // for plane case
  // c = pow(Pi,1./2)*h/curv;
  // for axisymmetrical case
  c = pow(4.*Pi/3.,1./3.)*h/curv;

  a_old = 0.; a_new = c;
 
  while(fabs(a_new-a_old)>1.e-6)
    {
      a_old = a_new; 
      dsda = -1. + 2./(1.+pow(1+2./a_old,1.-h)) +
        4.*(a_old+1)*pow(1.+2./a_old,-h)*(1.-h)/
        (a_old*a_old*pow(1+pow(1+2./a_old,1.-h),2));
      a_new = a_old - (S(a_old,h)-c)/dsda;
    }
  return a_new;
}
