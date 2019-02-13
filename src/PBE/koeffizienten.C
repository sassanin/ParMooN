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
   

#include "koeffizienten.h"
#include "basictools.h"

/*Die folgenden Funktionen sind groesstenteils aus dem eindimensionalen FEPC Algorithmus uebernommen wurden. Die Theorie
 zur Berechnung dieser Koeffizienten ist in einem Paper von Prof. Hackbusch Preprint no.: 38 aus dem Jahr 2007 beschrieben.
 Desweiteren sei auch hier auf die Dokumentation verwiesen*/


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/



/*Umwandeln einer Reellen Zahl x in die naechstgelegene Rationale Zahl der Form r = z/n.

Genau: Ermitteln falls moeglich einer rationalen Zahl der Form r= z/n mit |r-x| <= eps.
Wobei n natuerliche Zahl,z ganze Zahl und |z|,n <= Denom gilt. Falls dies nicht moeglich
ist, wird x zurueckgegeben.

Achtung: Damit der Algorithmus auch Sinn macht muss folgende Ungleichung gelten:
(D*D)/2 > eps, wobei D = 1/Denom ist.

Je groesser Denom ist, desto genauer arbeitet der Algorithmus, dh. desto eher werden korrekte
Werte ermittelt.
Empfehlung: eps = 1e-14 als Default Wert. Denom kann dann bis 1e6 gewaehlt werden (bei groesseren
Werten arbeitet Algorithmus falsch). Nur Denom bestimmt die Genauigkeit. Je groesser desto besser*/

static fepc_real_t
get_rational(fepc_real_t  x, fepc_real_t  eps, long Denom) {
	long  MaxLong;
	fepc_real_t  d, dabs, r;
	long  z, n;

	MaxLong = 2147483647L;

	for (n=1;n<=Denom;n++) {
		d = x*n;
		if (fabs(d) > MaxLong) {
			/*keine passende rat. Zahl gefunden*/
			return x;
		}
		z = (long) floor(d+0.5);
		r = (fepc_real_t)z/(fepc_real_t)n;
		d = x - r;
		dabs = fabs(d);
		if (dabs <= eps) {
			/*passende rationale Zahl gefunden*/
			return r;
		}
	}
	/*falls keine passende rationale Zahl gefunden*/

	return x;
}



/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/


matrix_p
matrix_new(int m,int n) {
	matrix_p  back;
	fepc_real_t  **a;
	int  k;

	back = (matrix_p) malloc(sizeof(matrix_t));
	ASSERT(back != NULL);

	ASSERT( m > 0 );
	ASSERT( n > 0 );

	a = (fepc_real_t**) malloc(sizeof(fepc_real_t*) * m);
	ASSERT(a != NULL);

	a[0] = (fepc_real_t*) malloc(sizeof(fepc_real_t) * m * n);
	ASSERT(a[0] != NULL);

	for (k=0;k<m;k++){
		a[k] = a[0] + k * n;
	}

	back->a = a;
	back->zeilen = m;
	back->spalten = n;

	return back;
}


void
matrix_del(matrix_p matrix) {
	free(*(matrix->a));
	free(matrix->a);
	free(matrix);
}




matrix_p
koeffizienten_xi_1dim(int grad) {
	matrix_p  back, matrix;
	int  q, m, n;
	fepc_real_t  **x;
	fepc_real_t  *a, *b;
	fepc_real_t  rat, reell, faktor, temp;

	ASSERT(grad >= 0);
	a = (fepc_real_t*) malloc(sizeof(fepc_real_t) * (2*grad + 1) );
	ASSERT(a != NULL);
	b = (fepc_real_t*) malloc(sizeof(fepc_real_t) * (2*grad + 1) );
	ASSERT(b != NULL);

	matrix = matrix_new( 2*grad+1, 2*grad+1 );
	x = matrix->a;

	for (n=0;n<=2*grad;n++) {
		a[n] = sqrt(2.*n+3)*sqrt(2.*n+1)/(fepc_real_t)(n+1);
	}
	b[0] = 0.0;
	for (n=1;n<=2*grad;n++) {
		b[n] = (sqrt(2.*n+3)/sqrt(2.*n-1)) * ( (fepc_real_t)n/(fepc_real_t)(n+1) );
	}

	x[0][0] = 1./sqrt(2.);
	for (q=1;q<=2*grad;q++) {
		for (n=0;n<=q;n++) {
			m = q-n;
			if (n<m) {
				x[n][m] = 0.0;
			}
			else {
				temp = a[n-1]/2.0*(x[n-1][m+1]/a[m]+x[n-1][m]);
				if (m>0) {
					temp = temp + a[n-1]/2.0*b[m]/a[m]*x[n-1][m-1];
				}
				if (n>=2) {
					temp = temp - b[n-1]*x[n-2][m];
				}

				/*Korrektur der Werte (siehe Dokumentation)*/
				faktor = pow(2.,((fepc_real_t)n+0.5))/sqrt( (2.*n+1)*(2*m+1) );
				reell = temp;
				rat = reell*faktor;
				rat = get_rational(rat,1e-14,10000);
				reell = rat/faktor;

				x[n][m] = reell;
			}
		}
	}

	back = matrix_new(grad+1,grad+1);
	for (m=0;m<=grad;m++) {
		for (n=0;n<=grad;n++) {
			back->a[m][n] = matrix->a[m][n];
		}
	}

	free(a);
	free(b);
	matrix_del(matrix);

	return back;
}



matrix3_p
koeffizienten_gamma_1dim(int grad) {
	int  a, b, k, n, q;
	fepc_real_t  ***g;
	matrix3_p  back, matrix;
	fepc_real_t  temp, wurzel, rat, reell;
	fepc_real_t  *A, *B;

	ASSERT(grad >= 0);
	A = (fepc_real_t*) malloc(sizeof(fepc_real_t)*(3*grad+1));
	ASSERT(A != NULL);
	B = (fepc_real_t*) malloc(sizeof(fepc_real_t)*(3*grad+1));
	ASSERT(B != NULL);

	for (n=0;n<=3*grad;n++) {
		A[n] = sqrt(2.*n+3)*sqrt(2.*n+1)/(fepc_real_t)(n+1);
	}
	B[0] = 0.0;
	for (n=1;n<=3*grad;n++) {
		B[n] = (sqrt(2.*n+3)/sqrt(2.*n-1)) * ( (fepc_real_t)n/(fepc_real_t)(n+1) );
	}

	matrix = matrix3_new(3*grad+1,3*grad+1,3*grad+1);
	g = matrix->a;

	for(b=0;b<=3*grad;b++) {
		for(k=0;k<=3*grad;k++) {
			g[0][b][k] = 0.0;
		}
	}
	for (k=1;k<=3*grad;k++) {
  	temp = 2.0*sqrt( (2.0*k-1)*(2.0*k+1) );
		g[0][k-1][k] = pow(-1.,k)/temp;
	}
	for (k=1;k<=3*grad;k++) {
		g[0][k][k-1] = g[0][k-1][k];
	}
	g[0][0][0] = 1.0/2.0;


	for(q=1;q<=3*grad;q++) {
		for(a=1;a<=q;a++) {
			for(b=0;b<=(q-a);b++) {
				k = q-a-b;

				if ( (a>b+k+1)||(b>a+k+1)||(k>a+b+1) ) {
					g[a][b][k] = 0.0;
				}
				if (b==0) {
					g[a][0][k] = pow(-1.,a) * g[0][a][k];
				}
				if (k==0) {
					g[a][b][0] = pow(-1.,a) * g[0][b][a];
				}
				if ( ((a>b+k+1)||(b>a+k+1)||(k>a+b+1)||(b==0)||(k==0)) == 0) {
					temp = A[a-1]*(g[a-1][b][k+1]/A[k] + g[a-1][b+1][k]/A[b] + g[a-1][b][k]);
					if (k>0)
						temp = temp + A[a-1]*B[k]*g[a-1][b][k-1]/A[k];
					if (b>0)
						temp = temp + A[a-1]*B[b]*g[a-1][b-1][k]/A[b];
					if (a>1)
						temp = temp - B[a-1]*g[a-2][b][k];
					g[a][b][k] = temp;
					/*Korrektur der Werte (siehe Dokumentation)*/
					wurzel = sqrt( (2.*a+1)*(2*b+1)*(2*k+1) );
					reell = temp;
					rat = reell*wurzel;
					rat = get_rational(rat,1e-14,10000);
					reell = rat/wurzel;
					g[a][b][k] = reell;
				}
			}
		}
	}

	back = matrix3_new(grad+1,grad+1,grad+1);
	for(a=0;a<=grad;a++) {
		for(b=0;b<=grad;b++) {
			for(k=0;k<=grad;k++) {
				back->a[a][b][k] = matrix->a[a][b][k];
			}
		}
	}

	free(A);
	free(B);
	matrix3_del(matrix);

	return back;
}

fepc_real_t
koeffizienten_gamma(int level, int r, fepc_real_t mesh, int alpha, int mu, int kappa, matrix3_p gamma_koef) {
	fepc_real_t  temp;
	int  c;

	/*Berechnen der eindimensionalen Gamma-Koeffizienten*/
	if(r == -1) {
		c = iipow(-1, alpha + mu + kappa);
	}
	else {
		c = 1;
	}

	temp = sqrt( mesh ) / pow(2.0, (fepc_real_t)level / 2.0 );
	return c * temp * gamma_koef->a[ alpha ][ mu ][ kappa ];
}


matrix3_p
matrix3_new(int d1,int d2,int d3) {
	matrix3_p  back;
	fepc_real_t  ***a;
	int  i, j;

	back = (matrix3_p) malloc(sizeof(matrix3_t));
	ASSERT(back != NULL);

	ASSERT( d1 > 0 );
	ASSERT( d2 > 0 );
	ASSERT( d3 > 0 );

	a = (fepc_real_t***) malloc(sizeof(fepc_real_t**) * d1);
	ASSERT(a != NULL);

	a[0] = (fepc_real_t**) malloc(sizeof(fepc_real_t*) * d1 * d2);
	ASSERT(a[0] != NULL);
	for (i=0;i<d1;i++) {
		a[i] = a[0] + i * d2;
	}

	a[0][0] = (fepc_real_t*) malloc(sizeof(fepc_real_t) * d1 * d2 * d3);
	ASSERT(a[0][0] != NULL);
	for(i=0;i<d1;i++) {
		for(j=0;j<d2;j++) {
			a[i][j] = a[0][0] + i*d2*d3 + j*d3;
		}
	}

	back->a = a;
	back->d1 = d1;
	back->d2 = d2;
	back->d3 = d3;

	return back;
}


void
matrix3_del(matrix3_p matrix) {
	free(**(matrix->a));
	free(*(matrix->a));
	free(matrix->a);
	free(matrix);
}


