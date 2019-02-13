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
   
#include "faltung.h"
#include "basictools.h"

static void
faltung_a1( func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p Gamma, func_p0_p back ) {
	static func_p0_p omega = func_p0_new( f, Gamma, f->maxlevel );

	static folgen_vektor_p0_p temp1 = folgen_vektor_p0_new();
	static folgen_vektor_p0_p temp2 = folgen_vektor_p0_new();

	folgen_vektor_p0_faltung( f->hierarchie[back->maxlevel], Gamma->hierarchie[back->maxlevel], omega->hierarchie[back->maxlevel]);
//	folgen_vektor_p0_print_aram(omega->hierarchie[back->maxlevel]);

	for(int l = back->maxlevel - 1; l >= 0; l--) {
		folgen_vektor_p0_faltung( f->hierarchie[l], Gamma->hierarchie[l], temp1 );
		faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef, temp2 );
		folgen_vektor_p0_add( temp1, temp2, omega->hierarchie[l] );
//		folgen_vektor_p0_print_aram(omega->hierarchie[l]);
	}

	func_p0_projekt( omega, w, back );
}

static void
faltung_a2( func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p0_p Gamma, func_p0_p back ) {

	static func_p0_p omega = func_p0_new(f, Gamma, f->maxlevel - 1);
	static folgen_vektor_p0_p temp1 = folgen_vektor_p0_new();
	static folgen_vektor_p0_p temp2 = folgen_vektor_p0_new();

	int min = g->maxlevel-1;

	folgen_vektor_p0_faltung(f->hierarchie[min], Gamma->hierarchie[min], omega->hierarchie[min]);

	for(int l = min - 1; l >= 0; l--) {
		folgen_vektor_p0_faltung( f->hierarchie[l], Gamma->hierarchie[l], temp1 );
		faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef, temp2 );
		folgen_vektor_p0_add( temp1, temp2, omega->hierarchie[l] );
	}

	func_p0_projekt( omega, w, back );
}

/* Berechnung der Algorithmen A1, B1, C1 (Fall l'<=l) */
static void
faltung_sum1( func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func_p0_p back ) {
	static func2_p0_p Gamma = func2_p0_new(w->maxlevel);

	/*Gamma bauen fuer den Fall A1 und B1*/
	faltung_hilfe_Gamma_bauen_1( g, mesh, gamma_koef, xi_koef, Gamma );
//	for(int i = 0; i < Gamma->maxlevel; i++)
//		folgen_matrix_p0_print(Gamma->hierarchie[i]);

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	faltung_a1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma, back );
}

/* Berechnung der Algorithmen A2, B2, C2 (Fall l'>l) */
static void
faltung_sum2( func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func_p0_p back ) {
	func_p0_p  temp1;

	static func2_p0_p Gamma = func2_p0_new(back->maxlevel - 1);

	/*die Funktionen f und g muessen vertauscht werden*/
	temp1 = f; f = g; g = temp1;

	/*Gamma bauen fuer den Fall A2*/
	faltung_hilfe_Gamma_bauen_2( g, mesh, gamma_koef, xi_koef, Gamma );

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	faltung_a2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma, back );
}

void
faltung_fepc(func_p0_p f, func_p0_p g, func_p0_p w, fepc_real_t h, func_p0_p back) {
	static matrix_p xi_koef = koeffizienten_xi_1dim( 1 );
	static matrix3_p gamma_koef = koeffizienten_gamma_1dim( 1 );
	static func_p0_p back1 = func_p0_new(back->maxlevel);
	static func_p0_p back2 = func_p0_new(back->maxlevel);

	/*Berechnung der Faltung*/
	faltung_sum1( f, g, w, h, 1, gamma_koef, xi_koef, back1 );
//	func_p0_print_aram(back1);
	if(back2->maxlevel > 0)
		faltung_sum2( f, g, w, h, 1, gamma_koef, xi_koef, back2 );

//	func_p0_print_aram(back2);
	func_p0_add( back1, back2, back );
//	func_p0_print_aram(back);

	func_p0_grid_zero( back );
}


static void
faltung_a1( func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p Gamma, func_p1_p back ) {
	static func_p1_p omega = func_p1_new( f, Gamma, f->maxlevel );

	static folgen_vektor_p1_p temp1 = folgen_vektor_p1_new();
	static folgen_vektor_p1_p temp2 = folgen_vektor_p1_new();

	folgen_vektor_p1_faltung( f->hierarchie[back->maxlevel], Gamma->hierarchie[back->maxlevel], omega->hierarchie[back->maxlevel]);
//	folgen_vektor_p1_print_aram(omega->hierarchie[back->maxlevel]);

	for(int l = back->maxlevel - 1; l >= 0; l--) {
		folgen_vektor_p1_faltung( f->hierarchie[l], Gamma->hierarchie[l], temp1 );
		faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef, temp2 );
		folgen_vektor_p1_add( temp1, temp2, omega->hierarchie[l] );
//		folgen_vektor_p1_print_aram(omega->hierarchie[l]);
	}

	func_p1_projekt( omega, w, back );
}
   
static void
faltung_a2( func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func2_p1_p Gamma, func_p1_p back ) {

	static func_p1_p omega = func_p1_new(f, Gamma, f->maxlevel - 1);
	static folgen_vektor_p1_p temp1 = folgen_vektor_p1_new();
	static folgen_vektor_p1_p temp2 = folgen_vektor_p1_new();

	int min = g->maxlevel-1;

	folgen_vektor_p1_faltung(f->hierarchie[min], Gamma->hierarchie[min], omega->hierarchie[min]);

	for(int l = min - 1; l >= 0; l--) {
		folgen_vektor_p1_faltung( f->hierarchie[l], Gamma->hierarchie[l], temp1 );
		faltung_hilfe_restrict( omega->hierarchie[l+1], xi_koef, temp2 );
		folgen_vektor_p1_add( temp1, temp2, omega->hierarchie[l] );
	}

	func_p1_projekt( omega, w, back );
}

/* Berechnung der Algorithmen A1, B1, C1 (Fall l'<=l) */
static void
faltung_sum1( func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func_p1_p back ) {
	static func2_p1_p Gamma = func2_p1_new(w->maxlevel);

	/*Gamma bauen fuer den Fall A1 und B1*/
	faltung_hilfe_Gamma_bauen_1( g, mesh, gamma_koef, xi_koef, Gamma );

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	faltung_a1( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma, back );
}

/* Berechnung der Algorithmen A2, B2, C2 (Fall l'>l) */
static void
faltung_sum2( func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t mesh, int maxgrad, matrix3_p gamma_koef, matrix_p xi_koef, func_p1_p back ) {
	func_p1_p  temp1;
	static func2_p1_p Gamma = func2_p1_new(back->maxlevel - 1);;

	/*die Funktionen f und g muessen vertauscht werden*/
	temp1 = f; f = g; g = temp1;

	/*Gamma bauen fuer den Fall A2*/
	faltung_hilfe_Gamma_bauen_2( g, mesh, gamma_koef, xi_koef, Gamma );

	/*Berechnung eines Teil der Faltung (Fall l'<=l)*/
	faltung_a2( f, g, w, mesh, maxgrad, gamma_koef, xi_koef, Gamma, back );
}

void
faltung_fepc(func_p1_p f, func_p1_p g, func_p1_p w, fepc_real_t h, func_p1_p back) {
	static matrix_p xi_koef = koeffizienten_xi_1dim( 3 );
	static matrix3_p gamma_koef = koeffizienten_gamma_1dim( 3 );
	static func_p1_p back1 = func_p1_new(back->maxlevel);
	static func_p1_p back2 = func_p1_new(back->maxlevel);

	/*Berechnung der Faltung*/
	faltung_sum1( f, g, w, h, 1, gamma_koef, xi_koef, back1 );
//	func_p1_print_aram(back1);
	if(back->maxlevel > 0)
		faltung_sum2( f, g, w, h, 1, gamma_koef, xi_koef, back2 );
//	func_p1_print_aram(back2);
	func_p1_add( back1, back2, back );
//	func_p1_print_aram(back);

	func_p1_grid_zero( back );
}
