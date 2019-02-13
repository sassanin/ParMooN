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
   
#include "folgen_vektor.h"
#include "basictools.h"

folgen_vektor_p0_p
folgen_vektor_p0_new() {
	folgen_vektor_p0_p  back;
	back = (folgen_vektor_p0_p) malloc(sizeof(folgen_vektor_p0_t));

	back->vektor = folge_new(0);
	ASSERT(back->vektor != NULL);

	return back;
}

void
folgen_vektor_p0_del(folgen_vektor_p0_p f) {
	folge_del(f->vektor);

	free(f);
}

folgen_matrix_p0_p
folgen_matrix_p0_new() {
	folgen_matrix_p0_p  back;
	folge_p0_p  matrix;

	back = (folgen_matrix_p0_p) malloc(sizeof(folgen_matrix_p0_t));
	ASSERT(back != NULL);

	matrix = folge_new(0);

	back->matrix = matrix;

	return back;
}

void
folgen_matrix_p0_del(folgen_matrix_p0_p f) {

	folge_del(f->matrix);

	free(f);
}

void
folgen_vektor_p0_print(folgen_vektor_p0_p f, int info) {
	printf("\n======================================================================");
	printf("\n#Ausgabe eines Folgenvektor der Dimension %d", 1);

	if(info==1) {
		folge_print( f->vektor, 0);
	}

	if(info==2) {
		folge_print( f->vektor, 1);
	}
	printf("\n======================================================================\n\n");
	return;
}

void
folgen_vektor_p0_print_aram(folgen_vektor_p0_p f) {

	printf("\n======================================================================\n");

	folge_print_aram( f->vektor);

	printf("\n");

	return;
}

void
folgen_vektor_p0_faltung(folgen_vektor_p0_p f, folgen_matrix_p0_p Gamma, folgen_vektor_p0_p back) {
	static folge_p0_p temp0 = folge_new(0);
	static folge_p0_p temp1 = folge_new(0);

	folge_faltung( f->vektor, Gamma->matrix, back->vektor );
}

void
folgen_vektor_p0_copy_structure(folgen_vektor_p0_p w, folgen_vektor_p0_p back) {
	back->vektor->lang = w->vektor->lang;
}

void
folgen_vektor_p0_projekt(folgen_vektor_p0_p f,folgen_vektor_p0_p g,folgen_vektor_p0_p back) {
	folge_projekt( f->vektor, g->vektor, back->vektor );
}

void
folgen_vektor_p0_add(folgen_vektor_p0_p f, folgen_vektor_p0_p g, folgen_vektor_p0_p back) {
	 
	folge_add( f->vektor, g->vektor, back->vektor );
}

void
folgen_matrix_p0_print(folgen_matrix_p0_p m)
{
	printf("\n\nooooooooooooooooooooooooooooooooooooo\n");
//	for(int i = 0; i < 2; i++)
//		for(int j = 0; j < 2; j++)
//			folge_print(m->matrix[i][j], 1);

	folge_print(m->matrix, 1);
}

void
folgen_matrix_p0_add(folgen_matrix_p0_p f, folgen_matrix_p0_p g, folgen_matrix_p0_p back) {
	folge_add( f->matrix, g->matrix, back->matrix );
}


folgen_vektor_p1_p
folgen_vektor_p1_new() {
	folgen_vektor_p1_p  back;
	back = (folgen_vektor_p1_p) malloc(sizeof(folgen_vektor_p1_t));

	back->vektor0 = folge_p1_new(0);
	ASSERT(back->vektor0 != NULL);
	back->vektor1 = folge_p1_new(0);
	ASSERT(back->vektor1 != NULL);

	return back;
}

void
folgen_vektor_p1_del(folgen_vektor_p1_p f) {
	folge_p1_del(f->vektor0);
	folge_p1_del(f->vektor1);

	free(f);
}

folgen_matrix_p1_p
folgen_matrix_p1_new() {
	folgen_matrix_p1_p  back;
	folge_p1_p  **matrix;
	int  i, j;


	back = (folgen_matrix_p1_p) malloc(sizeof(folgen_matrix_p1_t));
	ASSERT(back != NULL);
	matrix = (folge_p1_p**) malloc(sizeof(folge_p1_p*) * 2);
	ASSERT(matrix != NULL);

	matrix[0] = (folge_p1_p*) malloc(4 * sizeof(folge_p1_p));
	matrix[1] = matrix[0] + 2;

	for(i=0;i<2;i++) {
		for(j=0;j<2;j++) {
			matrix[i][j] = folge_p1_new(0);
		}
	}
	back->matrix = matrix;

	return back;
}

void
folgen_matrix_p1_del(folgen_matrix_p1_p f) {
	int  i, j;

	for(i=0;i<2;i++) {
		for(j=0;j<2;j++) {
			folge_p1_del(f->matrix[i][j]);
		}
	}
	free(*(f->matrix));
	free(f->matrix);

	free(f);
}

void
folgen_vektor_p1_print(folgen_vektor_p1_p f, int info) {
	printf("\n======================================================================");
	printf("\n#Ausgabe eines Folgenvektor der Dimension %d", 1);

	if(info==1) {
		folge_p1_print( f->vektor0, 0);
		folge_p1_print( f->vektor1, 0);
	}

	if(info==2) {
		folge_p1_print( f->vektor0, 1);
		folge_p1_print( f->vektor1, 1);
	}
	printf("\n======================================================================\n\n");
	return;
}

void
folgen_vektor_p1_print_aram(folgen_vektor_p1_p f) {

	printf("\n======================================================================\n");

	folge_p1_print_aram( f->vektor0);
	folge_p1_print_aram( f->vektor1);

	printf("\n");

	return;
}

void
folgen_vektor_p1_faltung(folgen_vektor_p1_p f, folgen_matrix_p1_p Gamma, folgen_vektor_p1_p back) {
	static folge_p1_p temp0 = folge_p1_new(0);
	static folge_p1_p temp1 = folge_p1_new(0);

	folge_p1_faltung( f->vektor0, Gamma->matrix[0][0], temp0 );
//	printf("0 0\n");
//	folge_p1_print_aram(temp0);

	folge_p1_faltung( f->vektor1, Gamma->matrix[0][1], temp1 );
//	printf("0 1\n");
//	folge_p1_print_aram(temp1);

//	folge_p1_print_aram(temp0);
//	folge_p1_print_aram(temp1);

	folge_p1_add( temp0, temp1, back->vektor0 );

	folge_p1_faltung( f->vektor0, Gamma->matrix[1][0], temp0 );
//	printf("1 0\n");
//	folge_p1_print_aram(temp0);

	folge_p1_faltung( f->vektor1, Gamma->matrix[1][1], temp1 );
//	printf("1 1\n");
//	folge_p1_print_aram(temp1);

	folge_p1_add( temp0, temp1, back->vektor1 );
}

void
folgen_vektor_p1_copy_structure(folgen_vektor_p1_p w, folgen_vektor_p1_p back) {
	back->vektor0->lang = w->vektor0->lang;

	back->vektor1->lang = w->vektor1->lang;
}

void
folgen_vektor_p1_projekt(folgen_vektor_p1_p f,folgen_vektor_p1_p g,folgen_vektor_p1_p back) {
	folge_p1_projekt( f->vektor0, g->vektor0, back->vektor0 );
	folge_p1_projekt( f->vektor1, g->vektor1, back->vektor1 );
}

void
folgen_vektor_p1_add(folgen_vektor_p1_p f, folgen_vektor_p1_p g, folgen_vektor_p1_p back) {
	 
	folge_p1_add( f->vektor0, g->vektor0, back->vektor0 );

	folge_p1_add( f->vektor1, g->vektor1, back->vektor1 );
}

void
folgen_matrix_p1_print(folgen_matrix_p1_p m)
{
	printf("\n\nooooooooooooooooooooooooooooooooooooo\n");
//	for(int i = 0; i < 2; i++)
//		for(int j = 0; j < 2; j++)
//			folge_p1_print(m->matrix[i][j], 1);

	folge_p1_print(m->matrix[0][0], 1);
}

void
folgen_matrix_p1_add(folgen_matrix_p1_p f, folgen_matrix_p1_p g, folgen_matrix_p1_p back) {
	int  i, j;

	for(i=0;i<2;i++) {
		for(j=0;j<2;j++) {
			folge_p1_add( f->matrix[i][j], g->matrix[i][j], back->matrix[i][j] );
		}
	}
}
