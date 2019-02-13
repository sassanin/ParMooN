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
   
/*
 * FEPC
 * Copyright (C) 2009 Peter Gerds (gerds@mis.mpg.de)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include <fftw3.h>
#include "fft_faltung.h"


/*******************************************************
 *
 * lokale Funktionen
 *
 *******************************************************/





/*******************************************************
 *
 * globale Funktionen
 *
 *******************************************************/
fepc_real_t*
fft_faltung(fepc_real_t* a, int n_a, fepc_real_t* b, int n_b) {
	int  size_a, size_b, size_c, dim=1;
	int  k, wert, test, temp;
	int  *n;
	int  n_c;
	fepc_real_t  *c;
	fftw_complex  *in, *out,*in_a, *out_a, *in_b, *out_b;
	fftw_plan  p;



	/*Auf Testen von Konsistenz wird verzichtet, da Input bereits auf Konsistenz getestet*/

	/*Faltung ueber Fouriertrafo (Theorie ist in Dokumentation zu finden)*/
	n_c = n_a + n_b - 1;
	n = &n_c;
	size_a = n_a;
	size_b = n_b;
	size_c = n_c;


	/*Initialisieren des Ergebnis Array*/
	c = (fepc_real_t*) malloc(sizeof(fepc_real_t) * size_c);


	/*Berechnen der Fouriertrafo von in_a*/
	in_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	for (k=0;k<size_c;k++) {
		temp = k;
		test = 0;
		if ((temp <0)||(temp >= n_a)) {
			test = test + 1;
		}

		if (test == 0) {
			wert = temp;
			in_a[k][0] = a[wert];
			in_a[k][1] = 0;
		}
		else {
			in_a[k][0] = 0;
			in_a[k][1] = 0;
		}
	}

	p = fftw_plan_dft(1,n,in_a,out_a,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);


	/*Berechnen der Fouriertrafo von in_b*/
	in_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	for (k=0;k<size_c;k++) {
		temp = k;
		test = 0;
		if ((temp <0)||(temp>=n_b)) {
			test = test + 1;
		}
		if (test == 0) {
			wert = temp;
			in_b[k][0] = b[wert];
			in_b[k][1] = 0;
		}
		else {
			in_b[k][0] = 0;
			in_b[k][1] = 0;
		}
	}

	p = fftw_plan_dft(dim,n,in_b,out_b,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);


	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size_c);

	for (k=0;k<size_c;k++) {
		in[k][0] = out_a[k][0]*out_b[k][0] - out_a[k][1]*out_b[k][1];
		in[k][1] = out_a[k][1]*out_b[k][0] + out_a[k][0]*out_b[k][1];
	}

	/*Berechnung der Inversen Fouriertrafo von in*/
	p = fftw_plan_dft(dim,n,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	for (k=0;k<size_c;k++) {
		c[k] = (fepc_real_t) out[k][0]/size_c;
	}

	fftw_free(in);
	fftw_free(in_a);
	fftw_free(in_b);
	fftw_free(out);
	fftw_free(out_a);
	fftw_free(out_b);
	return c;
}

