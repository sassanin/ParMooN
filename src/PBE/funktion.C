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
   
#include "funktion.h"
#include "basictools.h"

func_p0_p
func_p0_new(int maxlevel) {
	func_p0_p  back;
	folgen_vektor_p0_p  *hierarchie;
	int  k;

	ASSERT(maxlevel >= 0);
	back = (func_p0_p) malloc(sizeof(func_p0_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_vektor_p0_p*) malloc(sizeof(folgen_vektor_p0_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_vektor_p0_new();
	}
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}

func_p0_p
func_p0_new(func_p0_p f, func2_p0_p G, int maxlevel) {
	func_p0_p  back;
	folgen_vektor_p0_p  *hierarchie;
	int  k;

	back = (func_p0_p) malloc(sizeof(func_p0_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_vektor_p0_p*) malloc(sizeof(folgen_vektor_p0_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_vektor_p0_new();
		hierarchie[k]->vektor = folge_new(f->hierarchie[k]->vektor->lang + G->hierarchie[k]->matrix->lang - 1);
	}

	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}

void
func_p0_del(func_p0_p f) {
    if (f && f != NULL) {
        int  k;

	    for(k=0;k<=f->maxlevel;k++) {
		    folgen_vektor_p0_del(f->hierarchie[k]);
	    }
	    free(f->hierarchie);
	    free(f);
    }
}


func2_p0_p
func2_p0_new(int maxlevel) {
	func2_p0_p  back;
	folgen_matrix_p0_p  *hierarchie;
	int  k;

	ASSERT(maxlevel >= 0);
	back = (func2_p0_p) malloc(sizeof(func2_p0_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_matrix_p0_p*) malloc(sizeof(folgen_matrix_p0_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_matrix_p0_new();
	}	
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}


void
func2_p0_del(func2_p0_p f) {
	int  k;

	for(k=0;k<=f->maxlevel;k++) {
		folgen_matrix_p0_del(f->hierarchie[k]);
	}
	free(f->hierarchie);
	free(f);
}



void
func_p0_projekt(func_p0_p f,func_p0_p g,func_p0_p back) {
	int  l;

	for(l=0;l<=g->maxlevel;l++) {
		if (l<=f->maxlevel) {
			folgen_vektor_p0_projekt( f->hierarchie[l], g->hierarchie[l], back->hierarchie[l] );
		}
		else {
			folgen_vektor_p0_copy_structure( g->hierarchie[l], back->hierarchie[l] );
		}
	}
}


void
func_p0_print(func_p0_p f, int info) {
	int  k;

	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	printf("\n#Ausgabe einer Funktion der Dimension %d und des Maxlevel %d ", 1, f->maxlevel);
	if( (info == 1) || (info == 2) || (info == 3) ) {
		info = info -1;
		for(k=0;k<=f->maxlevel;k++) {
			folgen_vektor_p0_print( f->hierarchie[k], info );
		}
	}
	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n\n");
	return;
}

void
func_p0_print_aram(func_p0_p f) {
	int  k;

	printf("\n#Ausgabe einer Funktion der Dimension %d und des Maxlevel %d ", 1, f->maxlevel);
	for(k=0;k<=f->maxlevel;k++)
		folgen_vektor_p0_print_aram( f->hierarchie[k]);

	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n");
	return;
}


void
func_p0_add(func_p0_p f, func_p0_p g, func_p0_p back) {
	for (int k=0;k<=back->maxlevel;k++) {
		folge_add(f->hierarchie[k]->vektor, g->hierarchie[k]->vektor, back->hierarchie[k]->vektor);
	}
}

void
func_p0_clean(func_p0_p f) {
	for(int l=0;l<=f->maxlevel;l++) {
		for(int i=0;i<f->hierarchie[l]->vektor->lang;i++) {
			f->hierarchie[l]->vektor->glied[i] = 0;
		}
	}
}

void
func_p0_grid_zero(func_p0_p f) {
	int  l, i;
	int  b, r;

	for(l=0;l<f->maxlevel;l++) 
	{
		for(i=0;i<f->hierarchie[l]->vektor->lang;i++) 
		{
			r = 2 * i;
			b = f->hierarchie[l+1]->vektor->lang - 1;

			if(!(r>b)) 
			{
				f->hierarchie[l]->vektor->glied[i] = 0;
			}
		}
	}
}

func_p1_p
func_p1_new(int maxlevel) {
	func_p1_p  back;
	folgen_vektor_p1_p  *hierarchie;
	int  k;

	ASSERT(maxlevel >= 0);
	back = (func_p1_p) malloc(sizeof(func_p1_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_vektor_p1_p*) malloc(sizeof(folgen_vektor_p1_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_vektor_p1_new();
	}
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}

func_p1_p
func_p1_new(func_p1_p f, func2_p1_p G, int maxlevel) {
	func_p1_p  back;
	folgen_vektor_p1_p  *hierarchie;
	int  k;

	back = (func_p1_p) malloc(sizeof(func_p1_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_vektor_p1_p*) malloc(sizeof(folgen_vektor_p1_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_vektor_p1_new();
		hierarchie[k]->vektor0 = folge_p1_new(f->hierarchie[k]->vektor0->lang + G->hierarchie[k]->matrix[0][0]->lang - 1);
		hierarchie[k]->vektor1 = folge_p1_new(f->hierarchie[k]->vektor0->lang + G->hierarchie[k]->matrix[0][0]->lang - 1);
	}

	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}



void
func_p1_del(func_p1_p f) {
    if (f && f != NULL) {
        int  k;

	    for(k=0;k<=f->maxlevel;k++) {
		    folgen_vektor_p1_del(f->hierarchie[k]);
	    }
	    free(f->hierarchie);
	    free(f);
    }
}


func2_p1_p
func2_p1_new(int maxlevel) {
	func2_p1_p  back;
	folgen_matrix_p1_p  *hierarchie;
	int  k;

	ASSERT(maxlevel >= 0);
	back = (func2_p1_p) malloc(sizeof(func2_p1_t));
	ASSERT(back != NULL);

	hierarchie = (folgen_matrix_p1_p*) malloc(sizeof(folgen_matrix_p1_p)*(maxlevel+1));
	ASSERT(hierarchie != NULL);

	for(k=0;k<=maxlevel;k++) {
		hierarchie[k] = folgen_matrix_p1_new();
	}
	back->hierarchie = hierarchie;
	back->maxlevel = maxlevel;

	return back;
}


void
func2_p1_del(func2_p1_p f) {
	int  k;

	for(k=0;k<=f->maxlevel;k++) {
		folgen_matrix_p1_del(f->hierarchie[k]);
	}
	free(f->hierarchie);
	free(f);
}



void
func_p1_projekt(func_p1_p f,func_p1_p g,func_p1_p back) {
	int  l;

	for(l=0;l<=g->maxlevel;l++) {
		if (l<=f->maxlevel) {
			folgen_vektor_p1_projekt( f->hierarchie[l], g->hierarchie[l], back->hierarchie[l] );
		}
		else {
			folgen_vektor_p1_copy_structure( g->hierarchie[l], back->hierarchie[l] );
		}
	}
}


void
func_p1_print(func_p1_p f, int info) {
	int  k;

	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
	printf("\n#Ausgabe einer Funktion der Dimension %d und des Maxlevel %d ", 1, f->maxlevel);
	if( (info == 1) || (info == 2) || (info == 3) ) {
		info = info -1;
		for(k=0;k<=f->maxlevel;k++) {
			folgen_vektor_p1_print( f->hierarchie[k], info );
		}
	}
	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n\n");
	return;
}

void
func_p1_print_aram(func_p1_p f) {
	int  k;

	printf("\n#Ausgabe einer Funktion der Dimension %d und des Maxlevel %d ", 1, f->maxlevel);
	for(k=0;k<=f->maxlevel;k++)
		folgen_vektor_p1_print_aram( f->hierarchie[k]);

	printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n\n");
	return;
}


void
func_p1_add(func_p1_p f, func_p1_p g, func_p1_p back) {
	for (int k=0;k<=back->maxlevel;k++) {
		folge_p1_add(f->hierarchie[k]->vektor0, g->hierarchie[k]->vektor0, back->hierarchie[k]->vektor0);
		folge_p1_add(f->hierarchie[k]->vektor1, g->hierarchie[k]->vektor1, back->hierarchie[k]->vektor1);
	}
}

void
func_p1_clean(func_p1_p f) {
	for(int l=0;l<=f->maxlevel;l++) {
		for(int i=0;i<f->hierarchie[l]->vektor0->lang;i++) {
			f->hierarchie[l]->vektor0->glied[i] = 0;
			f->hierarchie[l]->vektor1->glied[i] = 0;
		}
	}
}

void
func_p1_grid_zero(func_p1_p f) {
	int  l, i;
	int  b, r;

	for(l=0;l<f->maxlevel;l++) 
	{
		for(i=0;i<f->hierarchie[l]->vektor0->lang;i++) 
		{
			r = 2 * i;
			b = f->hierarchie[l+1]->vektor0->lang - 1;

			if(!(r>b)) 
			{
				f->hierarchie[l]->vektor0->glied[i] = 0;
				f->hierarchie[l]->vektor1->glied[i] = 0;
			}
		}
	}
}
