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
   
#include "folge.h"
#include "basictools.h"
//#include "timer.h"

int MAX_POSSIBLE_FOLGE_LENGTH = 1;

folge_p0_p
folge_new(int lang) {
	folge_p0_p back;

	back = (folge_p0_p) malloc(sizeof(folge_p0_t));
	ASSERT(back != NULL);

	back->glied = (fepc_real_t*) malloc(sizeof(fepc_real_t)*MAX_POSSIBLE_FOLGE_LENGTH);
	ASSERT(back->glied != NULL);

	for (int k = 0; k < MAX_POSSIBLE_FOLGE_LENGTH; k++) {
		back->glied[k] = 0.0;
	}

	back->lang = lang;

	return back;
}

void
folge_del(folge_p0_p f) {
	free(f->glied);
	free(f);
}




void
folge_faltung(folge_p0_p f, folge_p0_p g, folge_p0_p back) {
/*	extern timer hp_faltung_timer;
	extern timer discrete_faltung_timer;
	hp_faltung_timer.pause_timer();
	discrete_faltung_timer.resume_timer(); */

	back->lang = f->lang + g->lang - 1;
	faltung( f->glied, f->lang, g->glied, g->lang, back->glied);

/*	discrete_faltung_timer.pause_timer();
	hp_faltung_timer.resume_timer(); */
}

void
folge_print(folge_p0_p f, int info) {
	int  k;

	printf("\n\t#:lang");
	printf("\t %d",f->lang);

	if(info == 1) {
		printf("\n\t#:glied");
		for(k=0;k<f->lang;k++) {
			printf("\t%.6lf",f->glied[k]);
		}
	}
	printf("\n------------------------------------------------------------\n");

	return;
}

void
folge_print_aram(folge_p0_p f) {
	int  k;

	for(k=0;k<f->lang;k++)
		printf("\t %.10lf ",f->glied[k]);

	printf("\n");
	return;
}


fepc_real_t
folge_glied( int r, folge_p0_p f ) {
	return (r < 0 || r > f->lang - 1) ? 0 : f->glied[r];
}

void
folge_add(folge_p0_p f, folge_p0_p g, folge_p0_p back) {
	int  k;

	if(!f->lang) {
		folge_copy( g, back );
		return;
	}

	if(!g->lang) {
		folge_copy( f, back );
		return;
	}

	back->lang = max(f->lang, g->lang);

	for(k=0;k<back->lang;k++) {
		back->glied[k] = folge_glied( k, f ) + folge_glied(k, g );
	}
}


void
folge_copy( folge_p0_p f, folge_p0_p back) {
	int  k;

	back->lang = f->lang;
	for(k=0;k<back->lang;k++) {
		back->glied[k] = f->glied[k];
	}
}



void
folge_projekt(folge_p0_p f, folge_p0_p g, folge_p0_p back) {
	back->lang = g->lang;

	for(int k = 0; k < back->lang; k++) {
		back->glied[k] = folge_glied( k, f );
	}
}

folge_p1_p
folge_p1_new(int lang) {
	folge_p1_p back;

	back = (folge_p1_p) malloc(sizeof(folge_p1_t));
	ASSERT(back != NULL);

	back->glied = (fepc_real_t*) malloc(sizeof(fepc_real_t)*MAX_POSSIBLE_FOLGE_LENGTH);
	ASSERT(back->glied != NULL);

	for (int k = 0; k < MAX_POSSIBLE_FOLGE_LENGTH; k++) {
		back->glied[k] = 0.0;
	}

	back->lang = lang;

	return back;
}


void
folge_p1_del(folge_p1_p f) {
	free(f->glied);
	free(f);
}

void
folge_p1_faltung(folge_p1_p f, folge_p1_p g, folge_p1_p back) {

/*	extern timer hp_faltung_timer;
	extern timer discrete_faltung_timer;
	hp_faltung_timer.pause_timer();
	discrete_faltung_timer.resume_timer(); */

	back->lang = f->lang + g->lang - 1;
	faltung( f->glied, f->lang, g->glied, g->lang, back->glied);

/*	discrete_faltung_timer.pause_timer();
	hp_faltung_timer.resume_timer(); */
}

void
folge_p1_print(folge_p1_p f, int info) {
	int  k;

	printf("\n\t#:lang");
	printf("\t %d",f->lang);

	if(info == 1) {
		printf("\n\t#:glied");
		for(k=0;k<f->lang;k++) {
			printf("\t%.6lf",f->glied[k]);
		}
	}
	printf("\n------------------------------------------------------------\n");

	return;
}

void
folge_p1_print_aram(folge_p1_p f) {
	int  k;

	for(k=0;k<f->lang;k++)
		printf("\t %.10lf ",f->glied[k]);

	printf("\n");
	return;
}


fepc_real_t
folge_p1_glied( int r, folge_p1_p f ) {
	return (r < 0 || r > f->lang - 1) ? 0 : f->glied[r];
}

void
folge_p1_add(folge_p1_p f, folge_p1_p g, folge_p1_p back) {
	int  k;

	if(!f->lang) {
		folge_p1_copy( g, back );
		return;
	}

	if(!g->lang) {
		folge_p1_copy( f, back );
		return;
	}

	back->lang = max(f->lang, g->lang);

	for(k=0;k<back->lang;k++) {
		back->glied[k] = folge_p1_glied( k, f ) + folge_p1_glied(k, g );
	}
}


void
folge_p1_copy( folge_p1_p f, folge_p1_p back) {
	int  k;

	back->lang = f->lang;
	for(k=0;k<back->lang;k++) {
		back->glied[k] = f->glied[k];
	}
}

void
folge_p1_projekt(folge_p1_p f, folge_p1_p g, folge_p1_p back) {
	back->lang = g->lang;

	for(int k = 0; k < back->lang; k++) {
		back->glied[k] = folge_p1_glied( k, f );
	}
}

void
faltung(fepc_real_t* a, int n_a, fepc_real_t* b, int n_b, fepc_real_t* c) {
	int n_c = n_a + n_b - 1;
	int minn, maxn;
	for(int i = 0; i < n_c; i++) {
		c[i] = 0; maxn = min(n_a - 1, i); minn = max(0, i - n_b + 1);
		for(int j = minn; j <= maxn; j++)
			c[i] += a[j] * b[i - j];
	}
}
