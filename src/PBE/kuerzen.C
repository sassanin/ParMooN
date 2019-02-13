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
   
#include "kuerzen.h"


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/



/* Bestimmen des Definitionsbereiches von F fuer den Fall B1 und B2, der bei der Faltung von F*Gamma
 von F benoetigt wird. */
static support_p
support_b( func_p f, func_p g, func_p w, func2_p Gamma ) {
	support_p  back, w_sup, Gamma_sup;
	int  l, k, j, dim;
	int  size1, size2, size, min;
	vec_p  n, n1, n2, temp_a, temp_b, temp, a, b, vec_1;


	ASSERT( f->dim == g->dim );
	ASSERT( f->dim == w->dim );
	ASSERT( f->dim == Gamma->dim );

	dim = f->dim;
	vec_1 = vec_one( dim );

	if (g->maxlevel < w->maxlevel) {
		min = g->maxlevel;
	}
	else {
		min = w->maxlevel;
	}

	/*Berechnen des kleinsten gemeinsamen Traegers pro Level von w*/
	w_sup = support_new( min, dim );
	for(l=0;l<=min;l++) {
		n = vec_add( w->hierarchie[l]->grad, vec_1 );
		size = vec_size( n );
		a = vec_copy( w->hierarchie[l]->vektor[0]->start );
		b = vec_add( a, w->hierarchie[l]->vektor[0]->lang );
		for(k=0;k<size;k++) {
			temp_a = w->hierarchie[l]->vektor[k]->start;
			temp_b = vec_add( temp_a, w->hierarchie[l]->vektor[k]->lang );
			temp = vec_min( a, temp_a );
			vec_del( a );
			a = temp;
			temp = vec_max( b, temp_b );
			vec_del( b );
			vec_del( temp_b );
			b = temp;
		}
		vec_del( n );
		vec_del( w_sup->a[l] );
		vec_del( w_sup->b[l] );
		w_sup->a[l] = a;
		w_sup->b[l] = vec_op( 1, b, -1, vec_1 );
		vec_del( b );
	}


	/*Berechnen des kleinsten gemeinsamen Traegers pro Level von Gamma*/
	Gamma_sup = support_new( min, dim );
	for(l=0;l<=min;l++) {
		n1 = vec_add( Gamma->hierarchie[l]->grad1, vec_1 );
		size1 = vec_size( n1 );
		n2 = vec_add( Gamma->hierarchie[l]->grad2, vec_1 );
		size2 = vec_size( n2 );
		a = vec_copy( Gamma->hierarchie[l]->matrix[0][0]->start );
		b = vec_add( a, Gamma->hierarchie[l]->matrix[0][0]->lang );
		for(k=0;k<size1;k++) {
			for(j=0;j<size2;j++) {
				temp_a = Gamma->hierarchie[l]->matrix[k][j]->start;
				temp_b = vec_add( temp_a, Gamma->hierarchie[l]->matrix[k][j]->lang );
				temp = vec_min( a, temp_a );
				vec_del( a );
				a = temp;
				temp = vec_max( b, temp_b );
				vec_del( b );
				vec_del( temp_b );
				b = temp;
			}
		}
		vec_del( n1 );
		vec_del( n2 );
		vec_del( Gamma_sup->a[l] );
		vec_del( Gamma_sup->b[l] );
		Gamma_sup->a[l] = a;
		Gamma_sup->b[l] = vec_op( 1, b, -1, vec_1 );
		vec_del( b );
	}

	back = support_new( min, dim );
	for(l=0;l<=min;l++) {
		vec_del( back->a[l] );
		vec_del( back->b[l] );
		back->a[l] = vec_op( 1, w_sup->a[l], -1, Gamma_sup->b[l] );
		back->b[l] = vec_op( 1, w_sup->b[l], -1, Gamma_sup->a[l] );
	}
	vec_del( vec_1 );
	support_del( w_sup );
	support_del( Gamma_sup );

	return back;
}



/* Diese Funktion berechnet wie sich der Traeger einer Funktion veraendert, falls diese
 nach folgender Methode gebaut wird: f->hierarchie[l] = Prolongation( f->hierarchie[l-1] ).
 Naehere Informationen darueber sind in der Dokumentation zu finden. */
static support_p
support_prolongation( support_p sup ) {
	support_p  back;
	int  l, dim, maxlevel;
	vec_p  *a, *b;
	vec_p  temp, temp_a, temp_b, vec_1;

	maxlevel = sup->maxlevel;
	a = sup->a;
	b = sup->b;
	dim = a[0]->dim;
	back = support_new( maxlevel, dim );
	vec_1 = vec_one( dim );

	vec_del( back->a[ maxlevel ] );
	vec_del( back->b[ maxlevel ] );
	back->a[ maxlevel ] = vec_copy( a[ maxlevel ] );
	back->b[ maxlevel ] = vec_copy( b[ maxlevel ] );

	for(l=maxlevel-1;l>=0;l--) {
		vec_del( back->a[ l ] );
		vec_del( back->b[ l ] );
		temp = vec_op( 1, back->a[ l+1 ], -1, vec_1 );
		temp_a = vec_div( 2, temp );
		vec_del( temp );
		temp_b = vec_div( 2, back->b[ l+1 ] );
		back->a[ l ] = vec_min( a[ l ], temp_a );
		back->b[ l ] = vec_max( b[ l ], temp_b );
		vec_del( temp_a );
		vec_del( temp_b );
	}
	vec_del( vec_1 );
	return back;
}





/* Bestimmen des Definitionsbereiches von F fuer den Fall C1 , der bei der Faltung von F*Gamma
 von F benoetigt wird. */
static support_p
support_c( func_p f, func_p g, func_p w, func2_p Gamma ) {
	support_p  back, w_sup, Gamma_sup;
	int  l, k, j, dim;
	int  size1, size2, min;
	vec_p  n1, n2, temp_a, temp_b, temp, a, b, vec_1;

	ASSERT( f->dim == g->dim );
	ASSERT( f->dim == w->dim );
	ASSERT( f->dim == Gamma->dim );

	dim = f->dim;
	vec_1 = vec_one( dim );

	if( g->maxlevel < (w->maxlevel-1) ) {
		min = g->maxlevel;
	}
	else {
		min = w->maxlevel-1;
	}

	/*Berechnen des kleinsten gemeinsamen Traegers pro Level von w,
	der fuer den Fall C notwendig ist*/
	w_sup = support_omega( w );

	/*Berechnen des kleinsten gemeinsamen Traegers pro Level von Gamma*/
	Gamma_sup = support_new( min, dim );
	for(l=0;l<=min;l++) {
		n1 = vec_add( Gamma->hierarchie[l]->grad1, vec_1 );
		size1 = vec_size( n1 );
		n2 = vec_add( Gamma->hierarchie[l]->grad2, vec_1 );
		size2 = vec_size( n2 );
		a = vec_copy( Gamma->hierarchie[l]->matrix[0][0]->start );
		b = vec_add( a, Gamma->hierarchie[l]->matrix[0][0]->lang );
		for(k=0;k<size1;k++) {
			for(j=0;j<size2;j++) {
				temp_a = Gamma->hierarchie[l]->matrix[k][j]->start;
				temp_b = vec_add( temp_a, Gamma->hierarchie[l]->matrix[k][j]->lang );
				temp = vec_min( a, temp_a );
				vec_del( a );
				a = temp;
				temp = vec_max( b, temp_b );
				vec_del( b );
				vec_del( temp_b );
				b = temp;
			}
		}
		vec_del( n1 );
		vec_del( n2 );
		vec_del( Gamma_sup->a[l] );
		vec_del( Gamma_sup->b[l] );
		Gamma_sup->a[l] = a;
		Gamma_sup->b[l] = vec_op( 1, b, -1, vec_1 );
		vec_del( b );
	}

	back = support_new( min, dim );
	for(l=0;l<=min;l++) {
		vec_del( back->a[l] );
		vec_del( back->b[l] );
		back->a[l] = vec_op( 1, w_sup->a[l], -1, Gamma_sup->b[l] );
		back->b[l] = vec_op( 1, w_sup->b[l], -1, Gamma_sup->a[l] );
	}
	vec_del( vec_1 );
	support_del( w_sup );
	support_del( Gamma_sup );

	return back;
}






/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/

func_p
kuerzen_F_bauen_b( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef ) {
	func_p  back, F;
	int  l, min, dim;
	support_p  sup_1, sup_2;
	folgen_vektor_p  temp1, temp2, temp3;

	dim = f->dim;
	if (g->maxlevel < w->maxlevel) {
		min = g->maxlevel;
	}
	else {
		min = w->maxlevel;
	}

	/*Bestimmen des Definitionsbereiches von F, der bei der Faltung von
	F*Gamma von F benoetigt wird*/
	sup_1 = support_b( f, g, w, Gamma );

	/*Bestimmen des daraus resultierenden Definitionsbereiches von F, der
	dann bei der Konstruktion von F selbst benoetigt wird*/
	sup_2 = support_prolongation( sup_1 );

	/*F bauen fuer den Fall B1 und B2, sowie Einkuerzen auf den fuer die
	Konstruktion notwendigen Bereiche*/
	F = func_new( min, dim );
	for(l=1;l<=min;l++) {
		if ( l > (f->maxlevel+1) ) {
			temp1 = F->hierarchie[l-1];
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( temp1, xi_koef );
		}
		else {
			temp1 = F->hierarchie[l-1];
			temp2 = f->hierarchie[l-1];
			temp3 = folgen_vektor_add( temp1, temp2 );
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( temp3, xi_koef );
			folgen_vektor_del( temp3 );
		}
		/*Einkuerzen*/
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_2->a[l], sup_2->b[l]);
		folgen_vektor_del( temp1 );
		F->hierarchie[l] = temp2;
	}


	/*Einkuerzen von F auf den Definitionsbereich, der bei der
	Faltung von F*Gamma von F tatsaechlich benoetigt wird*/
	back = func_new( min, dim );
	for(l=0;l<=min;l++) {
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_1->a[l], sup_1->b[l] );
		folgen_vektor_del( back->hierarchie[l] );
		back->hierarchie[l] = temp2;
	}

	func_del( F );
	support_del( sup_1 );
	support_del( sup_2 );

	return back;
}



func_p
kuerzen_F_bauen_c1( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef ) {
	func_p  back, F;
	int  l, min, dim;
	support_p  sup_1, sup_2;
	folgen_vektor_p  temp1, temp2, temp3;

	dim = f->dim;
	if( g->maxlevel < (w->maxlevel-1) ) {
		min = g->maxlevel;
	}
	else {
		min = w->maxlevel-1;
	}
	if( min < 0 ) {
		back = func_new( 0, dim );
		return back;
	}


	/*Bestimmen des Definitionsbereiches von F, der bei der Faltung von
	F*Gamma von F benoetigt wird*/
	sup_1 = support_c( f, g, w, Gamma );


	/*Bestimmen des daraus resultierenden Definitionsbereiches von F, der
	dann bei der Konstruktion von F selbst benoetigt wird*/
	sup_2 = support_prolongation( sup_1 );

	/*F bauen fuer den Fall C1, sowie Einkuerzen auf den fuer die
	Konstruktion notwendigen Bereiche*/
	F = func_new( min, dim );
	for(l=0;l<=min;l++) {
		if (l>0) {
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( F->hierarchie[l-1], xi_koef );
		}
		if (l<=f->maxlevel) {
			temp1 = F->hierarchie[l];
			temp2 = folgen_vektor_add( temp1, f->hierarchie[l] );
			folgen_vektor_del( temp1 );
			F->hierarchie[l] = temp2;
		}
		/*Einkuerzen*/
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_2->a[l], sup_2->b[l]);
		folgen_vektor_del( temp1 );
		F->hierarchie[l] = temp2;
	}

	/*Einkuerzen von F auf den Definitionsbereich, der bei der
	Faltung von F*Gamma von F tatsaechlich benoetigt wird*/
	back = func_new( min, dim );
	for(l=0;l<=min;l++) {
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_1->a[l], sup_1->b[l] );
		folgen_vektor_del( back->hierarchie[l] );
		back->hierarchie[l] = temp2;
	}

	func_del( F );
	support_del( sup_1 );
	support_del( sup_2 );

	return back;
}



support_p
support_omega( func_p w ) {
	support_p  back, w_sup;
	int  l, k, dim, size;
	vec_p  n, temp_a, temp_b, temp, a, b, vec_1;


	dim = w->dim;
	vec_1 = vec_one( dim );

	/*Berechnen des kleinsten gemeinsamen Traegers pro Level von w*/
	w_sup = support_new( w->maxlevel, dim );
	for(l=0;l<=w->maxlevel;l++) {
		n = vec_add( w->hierarchie[l]->grad, vec_1 );
		size = vec_size( n );
		a = vec_copy( w->hierarchie[l]->vektor[0]->start );
		b = vec_add( a, w->hierarchie[l]->vektor[0]->lang );
		for(k=0;k<size;k++) {
			temp_a = w->hierarchie[l]->vektor[k]->start;
			temp_b = vec_add( temp_a, w->hierarchie[l]->vektor[k]->lang );
			temp = vec_min( a, temp_a );
			vec_del( a );
			a = temp;
			temp = vec_max( b, temp_b );
			vec_del( b );
			vec_del( temp_b );
			b = temp;
		}
		vec_del( n );
		vec_del( w_sup->a[l] );
		vec_del( w_sup->b[l] );
		w_sup->a[l] = a;
		w_sup->b[l] = vec_op( 1, b, -1, vec_1 );
		vec_del( b );
	}


	/*Bestimmen des daraus resultierenden Definitionsbereiches von F, der
	dann bei der Konstruktion von F selbst benoetigt wird*/
	back = support_prolongation( w_sup );

	vec_del( vec_1 );
	support_del( w_sup );

	return back;
}



folgen_vektor_p
folgen_vektor_support( folgen_vektor_p f, vec_p a, vec_p b ) {
	folgen_vektor_p  back;
	int  k, i, size, size_lang, dim, test, j;
	vec_p  n, grad, r, s, start, lang;
	vec_p  vec_1, temp1, temp2;
	folge_p  h;


	dim = f->grad->dim;
	grad = vec_copy( f->grad );
	ASSERT( dim == a->dim );
	ASSERT( dim == b->dim );

	vec_1 = vec_one( dim );
	n = vec_add( grad, vec_1 );
	size = vec_size( n );
	back = folgen_vektor_new( grad );
	for(k=0;k<size;k++) {
		start = vec_max( a, f->vektor[k]->start );
		temp1 = vec_add( f->vektor[k]->start, f->vektor[k]->lang );
		temp2 = vec_op( 1, temp1, -1, vec_1 );
		vec_del( temp1 );
		temp1 = vec_min( b, temp2 );
		vec_del( temp2 );
		temp2 = vec_op( 1, temp1, -1, start );
		vec_del( temp1 );
		lang = vec_add( temp2, vec_1 );
		vec_del( temp2 );

		test = 0;
		for(j=0;j<dim;j++) {
			if( lang->array[j] <= 0 ) {
				test = test + 1;
			}
		}
		if( test == 0 ) {
			h = folge_new( start, lang );
			size_lang = vec_size( lang );
			for(i=0;i<size_lang;i++) {
				r = entry_one2d( i, lang );
				s = vec_add( r, start );
				vec_del( r );
				h->glied[i] = folge_glied( s , f->vektor[k] );
				vec_del( s );
			}
			folge_del( back->vektor[k] );
			back->vektor[k] = h;
		}
		else {
			vec_del( start );
			vec_del( lang );
		}
	}
	vec_del( vec_1 );
	vec_del( n );

	return back;
}


support_p
support_new( int maxlevel, int dim ) {
	support_p  back;
	vec_p  *a, *b;
	int  k;

	ASSERT( maxlevel >= 0 );
	ASSERT( dim > 0 );

	back = (support_p) malloc(sizeof(support_t));
	ASSERT(back != NULL);

	a = (vec_p*) malloc( sizeof(vec_p) * (maxlevel+1) );
	ASSERT(a != NULL);
	b = (vec_p*) malloc( sizeof(vec_p) * (maxlevel+1) );
	ASSERT(b != NULL);

	for(k=0;k<=maxlevel;k++) {
		a[k] = vec_new( dim );
		b[k] = vec_new( dim );
	}
	back->a = a;
	back->b = b;
	back->maxlevel = maxlevel;
	return back;
}


void
support_del( support_p sup ) {
	int  k;

	for(k=0;k<=sup->maxlevel;k++) {
		vec_del( sup->a[k] );
		vec_del( sup->b[k] );
	}
	free(sup->a);
	free(sup->b);
	free(sup);
}



func_p
kuerzen_F_bauen_c2( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef ) {
	func_p  back, F;
	int  l, min, dim;
	support_p  sup_1, sup_2;
	folgen_vektor_p  temp1, temp2, temp3;

	dim = f->dim;
	if( g->maxlevel < (w->maxlevel-1) ) {
		min = g->maxlevel;
	}
	else {
		min = w->maxlevel-1;
	}
	if( min < 0 ) {
		back = func_new( 0, dim );
		return back;
	}

	/*Bestimmen des Definitionsbereiches von F, der bei der Faltung von
	F*Gamma von F benoetigt wird*/
	sup_1 = support_c( f, g, w, Gamma );

	/*Bestimmen des daraus resultierenden Definitionsbereiches von F, der
	dann bei der Konstruktion von F selbst benoetigt wird*/
	sup_2 = support_prolongation( sup_1 );

	/*F bauen fuer den Fall C2, sowie Einkuerzen auf den fuer die
	Konstruktion notwendigen Bereiche*/
	F = func_new( min, dim );
	for(l=1;l<=min;l++) {
		if ( l > (f->maxlevel+1) ) {
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( F->hierarchie[l-1], xi_koef );
		}
		else {
			temp1 = F->hierarchie[l-1];
			temp2 = f->hierarchie[l-1];
			temp3 =	folgen_vektor_add( temp1, temp2 );
			folgen_vektor_del( F->hierarchie[l] );
			F->hierarchie[l] = faltung_hilfe_prolong( temp3, xi_koef );
			folgen_vektor_del( temp3 );
		}
		/*Einkuerzen*/
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_2->a[l], sup_2->b[l]);
		folgen_vektor_del( temp1 );
		F->hierarchie[l] = temp2;
	}

	/*Einkuerzen von F auf den Definitionsbereich, der bei der
	Faltung von F*Gamma von F tatsaechlich benoetigt wird*/
	back = func_new( min, dim );
	for(l=0;l<=min;l++) {
		temp1 = F->hierarchie[l];
		temp2 = folgen_vektor_support( temp1, sup_1->a[l], sup_1->b[l] );
		folgen_vektor_del( back->hierarchie[l] );
		back->hierarchie[l] = temp2;
	}

	func_del( F );
	support_del( sup_1 );
	support_del( sup_2 );

	return back;
}
