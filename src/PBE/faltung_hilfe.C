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
   
#include "faltung_hilfe.h"
#include "basictools.h"


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/

static void
lambda(folgen_vektor_p0_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, int alpha, int mu, folge_p0_p back) {
	int  j, lang, alpha_plus_mu = 1 - 2 * ((alpha + mu) % 2);
	fepc_real_t  factor = sqrt(mesh) / pow(2.0, (fepc_real_t)level / 2.0);

	lang = g->vektor->lang + 1;

	back->lang = lang;

	for(j = 1; j < lang - 1; j++) {

		back->glied[j] = factor * 
						( gamma_koef->a[alpha][mu][0] * (g->vektor->glied[j] + g->vektor->glied[j-1]) );
	}

	back->glied[0] = factor * 
					( gamma_koef->a[alpha][mu][0] * g->vektor->glied[0] );

	back->glied[lang-1] = factor * 
						( gamma_koef->a[alpha][mu][0] * g->vektor->glied[lang-2] );
}

static void
LAMBDA( folgen_matrix_p0_p Gamma, matrix_p xi_koef, folge_p0_p back ) {
	folge_p0_p  matrix; fepc_real_t* glied;
	fepc_real_t** xi_a;
	int  i, a, b, imax, lang;
	int  nu;
	fepc_real_t  sum, sub_sum, wert, x;

	matrix = Gamma->matrix;
	glied = back->glied;
	xi_a = xi_koef->a;
	imax = matrix->lang;

	imax = max( imax, matrix->lang );

	lang = imax / 2 + 1;

	back->lang = lang;

	sum = 0.0;

	for(i = 0; i < lang; i++) {
		sum = 0.0; nu = 2 * i;

		x = xi_a[0][0] * xi_a[0][0];
		sub_sum = 0.0;

		for(a = 0; a < 2; a++) {
			for(b = 0; b < 2; b++) {
				wert = folge_glied(nu + a - b, matrix);
				if( wert != 0.0 ) {
					sub_sum += wert;
				}
			}
		}
		sum += x * sub_sum;

		glied[i] = sum;
	}
}

void
faltung_hilfe_restrict(folgen_vektor_p0_p f, matrix_p xi_koef, folgen_vektor_p0_p back) {
	int  lang, nu, j;
	
	lang = (f->vektor->lang - 1) / 2 + 1;

	back->vektor->lang = lang;

	for(j = 0, nu = 0; j < lang - 1; j++) {
		back->vektor->glied[j] = xi_koef->a[0][0] * (f->vektor->glied[nu] + f->vektor->glied[nu + 1]);

		nu += 2;
	}
	back->vektor->glied[lang - 1] = xi_koef->a[0][0] * (f->vektor->glied[nu] + folge_glied(nu+1, f->vektor));
}

void
faltung_hilfe_lambda(folgen_vektor_p0_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, folgen_matrix_p0_p back) {
	lambda( g, level, mesh, gamma_koef, 0, 0, back->matrix );
}

void
faltung_hilfe_LAMBDA(folgen_matrix_p0_p Gamma, matrix_p xi_koef, folgen_matrix_p0_p back) {
	LAMBDA( Gamma, xi_koef, back->matrix );
}

void
faltung_hilfe_Gamma_bauen_1( func_p0_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p back ) {
	static folgen_matrix_p0_p  temp1 = folgen_matrix_p0_new();
	static folgen_matrix_p0_p  temp2 = folgen_matrix_p0_new();
	int  l;

	faltung_hilfe_lambda( g->hierarchie[g->maxlevel], g->maxlevel, mesh, gamma_koef, back->hierarchie[g->maxlevel]);
//	folgen_matrix_p0_print(back->hierarchie[g->maxlevel]);

	for(l=g->maxlevel-1;l>=0;l--) 
	{
		faltung_hilfe_LAMBDA( back->hierarchie[l+1], xi_koef, temp1 );
//		printf("LAMBDA\n");
//		folgen_matrix_p0_print(temp1);
		faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, temp2);
//		printf("lambda\n");
//		folgen_matrix_p0_print(temp2);
		folgen_matrix_p0_add(temp1,temp2,back->hierarchie[l]);

//		folgen_matrix_p0_print(back->hierarchie[l]);
//		printf("%d\n", l);
	}
}

void
faltung_hilfe_Gamma_bauen_2( func_p0_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p back ) {
	static folgen_matrix_p0_p  temp1 = folgen_matrix_p0_new();
	static folgen_matrix_p0_p  temp2 = folgen_matrix_p0_new();
	static folgen_matrix_p0_p  temp3 = folgen_matrix_p0_new();
	int  l, maxlevel = g->maxlevel - 1;

	faltung_hilfe_lambda( g->hierarchie[g->maxlevel],g->maxlevel, mesh, gamma_koef, temp1 );
	faltung_hilfe_LAMBDA( temp1, xi_koef, back->hierarchie[maxlevel] );

	for(l=maxlevel-1;l>=0;l--) 
	{
		faltung_hilfe_lambda( g->hierarchie[l+1], l+1, mesh, gamma_koef, temp2 );
		folgen_matrix_p0_add( back->hierarchie[l+1], temp2, temp3 );
		faltung_hilfe_LAMBDA( temp3, xi_koef, back->hierarchie[l] );
	}
}

static void
lambda(folgen_vektor_p1_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, int alpha, int mu, folge_p1_p back) {
	int  j, lang, alpha_plus_mu = 1 - 2 * ((alpha + mu) % 2);
	fepc_real_t  factor = sqrt(mesh) / pow(2.0, (fepc_real_t)level / 2.0);

	lang = g->vektor0->lang + 1;

	back->lang = lang;

	for(j = 1; j < lang - 1; j++) {

		back->glied[j] = factor * 
						( gamma_koef->a[alpha][mu][0] * (g->vektor0->glied[j] + alpha_plus_mu * g->vektor0->glied[j-1]) +
						  gamma_koef->a[alpha][mu][1] * (g->vektor1->glied[j] - alpha_plus_mu * g->vektor1->glied[j-1]) );
	}

	back->glied[0] = factor * 
					( gamma_koef->a[alpha][mu][0] * g->vektor0->glied[0] +
					  gamma_koef->a[alpha][mu][1] * g->vektor1->glied[0] );

	back->glied[lang-1] = factor * 
						( gamma_koef->a[alpha][mu][0] * alpha_plus_mu * g->vektor0->glied[lang-2] -
						  gamma_koef->a[alpha][mu][1] * alpha_plus_mu * g->vektor1->glied[lang-2] );
}

static void
LAMBDA( folgen_matrix_p1_p Gamma, matrix_p xi_koef, int alpha, int mu, folge_p1_p back ) {
	folge_p1_p  **matrix; fepc_real_t* glied;
	fepc_real_t** xi_a;
	int  i, j, k, a, b, imax, lang;
	int  nu;
	fepc_real_t  sum, sub_sum, wert, x;

	matrix = Gamma->matrix;
	glied = back->glied;
	xi_a = xi_koef->a;
	imax = matrix[0][0]->lang;

	for(k=0;k<=alpha;k++) {
		for(j=0;j<=mu;j++) {
			imax = max( imax, matrix[k][j]->lang );
		}
	}

	lang = imax / 2 + 1;

	back->lang = lang;

	sum = 0.0;
	for(k = 0; k <= alpha; k++) {
		for(j = 0; j <= mu; j++) {
			x = xi_a[alpha][k] * xi_a[mu][j];
			sub_sum = 0.0;

			for(a = 0; a < 2; a++) {
				for(b = 0; b <= a; b++) {
					wert = folge_p1_glied(a - b, matrix[k][j]);
					if( wert != 0.0 ) {
						sub_sum += wert * (1 - 2 * (((k + alpha) * (1 - a) + (j + mu) * (1 - b)) % 2));
					}
				}
			}
			sum += x * sub_sum;
		}
	}
	glied[0] = sum;

	for(i = 1; i < lang; i++) {
		sum = 0.0; nu = 2 * i;

		for(k = 0; k <= alpha; k++) {
			for(j = 0; j <= mu; j++) {
				x = xi_a[alpha][k] * xi_a[mu][j];
				sub_sum = 0.0;

				for(a = 0; a < 2; a++) {
					for(b = 0; b < 2; b++) {
						wert = folge_p1_glied(nu + a - b, matrix[k][j]);
						if( wert != 0.0 ) {
							sub_sum += wert * (1 - 2 * (((k + alpha) * (1 - a) + (j + mu) * (1 - b)) % 2));
						}
					}
				}
				sum += x * sub_sum;
			}
		}
		glied[i] = sum;
	}
}

void
faltung_hilfe_restrict(folgen_vektor_p1_p f, matrix_p xi_koef, folgen_vektor_p1_p back) {
	int  lang, nu, j;
	
	lang = (f->vektor0->lang - 1) / 2 + 1;

	back->vektor0->lang = lang;
	back->vektor1->lang = lang;

	for(j = 0, nu = 0; j < lang - 1; j++) {
		back->vektor0->glied[j] = xi_koef->a[0][0] * (f->vektor0->glied[nu] + f->vektor0->glied[nu + 1]);

		nu += 2;
	}
	back->vektor0->glied[lang - 1] = xi_koef->a[0][0] * (f->vektor0->glied[nu] + folge_p1_glied(nu+1, f->vektor0));

	for(j = 0, nu = 0; j < lang - 1; j++) {
		back->vektor1->glied[j] = xi_koef->a[1][0] * (f->vektor0->glied[nu + 1] - f->vektor0->glied[nu]) +
								  xi_koef->a[1][1] * (f->vektor1->glied[nu + 1] + f->vektor1->glied[nu]);
		
		nu += 2;
	}
	back->vektor1->glied[lang - 1] = xi_koef->a[1][0] * (folge_p1_glied(nu+1, f->vektor0) - f->vektor0->glied[nu]) +
									 xi_koef->a[1][1] * (folge_p1_glied(nu+1, f->vektor1) + f->vektor1->glied[nu]);
}

void
faltung_hilfe_lambda(folgen_vektor_p1_p g, int level, fepc_real_t mesh, matrix3_p gamma_koef, folgen_matrix_p1_p back) {
	int  j, k;

	for(k=0;k<2;k++) {
		for(j=0;j<2;j++) {
			lambda( g, level, mesh, gamma_koef, k, j, back->matrix[k][j] );
		}
	}
}

void
faltung_hilfe_LAMBDA(folgen_matrix_p1_p Gamma, matrix_p xi_koef, folgen_matrix_p1_p back) {
	int  k, j;

	for(k=0; k<2; k++) {
		for(j=0; j<2; j++) {
			LAMBDA( Gamma, xi_koef, k, j, back->matrix[k][j] );
		}
	}
}

void
faltung_hilfe_Gamma_bauen_1( func_p1_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p back ) {
	static folgen_matrix_p1_p  temp1 = folgen_matrix_p1_new();
	static folgen_matrix_p1_p  temp2 = folgen_matrix_p1_new();
	int  l;

	faltung_hilfe_lambda( g->hierarchie[g->maxlevel], g->maxlevel, mesh, gamma_koef, back->hierarchie[g->maxlevel]);
//	folgen_matrix_p1_print(back->hierarchie[g->maxlevel]);

	for(l=g->maxlevel-1;l>=0;l--) 
	{
		faltung_hilfe_LAMBDA( back->hierarchie[l+1], xi_koef, temp1 );
//		printf("LAMBDA\n");
//		folgen_matrix_p1_print(temp1);
		faltung_hilfe_lambda( g->hierarchie[l], l, mesh, gamma_koef, temp2);
//		printf("lambda\n");
//		folgen_matrix_p1_print(temp2);
		folgen_matrix_p1_add(temp1,temp2,back->hierarchie[l]);

//		folgen_matrix_p1_print(back->hierarchie[l]);
//		printf("%d\n", l);
	}
}

void
faltung_hilfe_Gamma_bauen_2( func_p1_p g, fepc_real_t mesh, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p back ) {
	static folgen_matrix_p1_p  temp1 = folgen_matrix_p1_new();
	static folgen_matrix_p1_p  temp2 = folgen_matrix_p1_new();
	static folgen_matrix_p1_p  temp3 = folgen_matrix_p1_new();
	int  l, maxlevel = g->maxlevel - 1;

	faltung_hilfe_lambda( g->hierarchie[g->maxlevel],g->maxlevel, mesh, gamma_koef, temp1 );
	faltung_hilfe_LAMBDA( temp1, xi_koef, back->hierarchie[maxlevel] );

	for(l=maxlevel-1;l>=0;l--) 
	{
		faltung_hilfe_lambda( g->hierarchie[l+1], l+1, mesh, gamma_koef, temp2 );
		folgen_matrix_p1_add( back->hierarchie[l+1], temp2, temp3 );
		faltung_hilfe_LAMBDA( temp3, xi_koef, back->hierarchie[l] );
	}
}
