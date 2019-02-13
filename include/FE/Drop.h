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
   
#include <math.h>
#include <Constants.h>
#include <MooNMD_Io.h>

// adap=0 - without adaptation; adap=1 - with adaptation
// a - parameter for grid movement, a=10 corresponds near the uniform grid;
int SchemeA(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);
int SchemeA_ax(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax1(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax11(double *r, double *z, double *h_s, double *h_n, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeA_ax2(double *r, double *z, double *h_r, double *h_z, 
	    int N, double hi, double gamma, double W, double tau, double *LL, 
	    int adap, double *a);   
int SchemeB(double *r, double *z, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL);
        
int SchemeT4(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL,
        double *F, double *dF);
int SchemeT4_ax(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL,
        double *F, double *dF, int adap, double *a, int up);    
int SchemeT4_axNL(double *r, double *z, double *beta, double *h_r, double *h_z, 
        int N, double hi, double gamma, double W, double tau, double *LL, 
        double *F, double *dF, int adap, double* a);   
void Solver_3diag(int N, double *c, double *d, double *e, double *b);

// a generator of a nonuniform grid
double S(double a, double t);

double dSdt(double a, double t);

// a new parameter for the generator of a nonuniform grid
double A(double curv, double h);
