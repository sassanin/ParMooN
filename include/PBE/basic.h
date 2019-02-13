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
   
#ifndef __BASIC_H
#define __BASIC_H

#include "config.h"

typedef struct {
	int dim;				/* Dimension des Vektors */
	int *array;			/* Array mit Werten */
} vec_t;

typedef vec_t *vec_p;


typedef struct {
	int dim;								/* Dimension des Vektors */
	fepc_real_t *array;			/* Array mit Werten */
} vec_real_t;

typedef vec_real_t *vec_real_p;


int int_round(fepc_real_t d);

/* Diese Routinen dienen der Verwaltung multidimensionaler Arrays. Erlaeuterungen zur Speicherung solcher Arrays,
sowie zur Funktionsweise der Funktionen entry_d2one und entry_one2d sind in der Dokumentation zu finden */

/* Initialisiert einen Vektor der Dimension dim, dessen Eintraege alle mit der Zahl 0 gefuellt sind */
vec_p
vec_new(int dim);

/* Initialisiert einen Vektor der Dimension dim, dessen Eintraege alle mit der Zahl 1.0 gefuellt sind */
vec_real_p
vec_real_new(int dim);

/* Initialisiert einen Vektor der Dimension dim, dessen Eintraege alle mit der Zahl 1 gefuellt sind */
vec_p
vec_one(int dim);

/* Loescht den Vektor s */
void
vec_del(vec_p s);

/* Loescht den Vektor s */
void
vec_real_del(vec_real_p s);

/* Berechnet die Position des Vektors s in einem multidimensionalen Array, dessen Dimension durch den Vektor n bestimmt ist */
int
entry_d2one(vec_p s,vec_p n);

/* Berechnet den zugehoerigen Vektor der Position pos eines multidimensionalen Array, dessen Dimension durch den Vektor n bestimmt ist */
vec_p
entry_one2d(int pos,vec_p n);

/* Berechnet den zugehoerigen Vektor der Position pos eines multidimensionalen Array, dessen Dimension durch den Vektor n bestimmt ist, ohne dabei streng zu sein */
vec_p
entry_one2d_sloppy(int pos,vec_p n);

/* Berechnet den Vektor a*s + b*n */
vec_p
vec_op(int a, vec_p s, int b, vec_p n);

/* Multipliziert den Vektor n mit der ganzen Zahl a */
vec_p
vec_multi(int a, vec_p n);

/* Multipliziert einen reelwertigen Vektor mit einer Konstante */
vec_real_p
vec_real_multi(fepc_real_t factor, vec_real_p vector);

/* Multipliziert einen reelwertigen Vektor mit einer Konstante und speichert das Ergebnis im Eingabevektor */
void
vec_real_multi2(fepc_real_t factor, vec_real_p vector);

/* Dividiert den Vektor n mit der ganzen Zahl a */
vec_p
vec_div(int a, vec_p n);

/* Addiert Vektor s und Vektor n */
vec_p
vec_add(vec_p s, vec_p n);

/* Addiert Vektor s und Vektor n und speichert das Ergebnis in s */
void
vec_add2(vec_p s, vec_p n);

/* Erzeugt eine inhaltliche Kopie des Vektors n */
vec_p
vec_copy(vec_p n);

/* Subtrahiert Vektor s und Vektor n */
vec_real_p
vec_real_substract(vec_real_p s, vec_real_p n);

/* Berechnet das Skalarprodukt der Vektoren r und s */
int
vec_skalar_prod(vec_p r, vec_p s);

/* Berechnet den groessten Vektor n sodass gilt: n<=r und n<=s */
vec_p
vec_min(vec_p r, vec_p s);

/* Berechnet den kleinsten Vektor n sodass gilt: n>=r und n>=s */
vec_p
vec_max(vec_p r, vec_p s);

/* Berechnet r_1 * ... * r_dim */
int
vec_size(vec_p r);

/* Gegeben seien die Vektoren r, n mit n>r>=0 und n > 0. Der Vektor n beschreibt ein multidimensionales Array.
Gesucht sind alle ganzen Zahlen t >= 0 fuer die gilt: r<=s<n mit s = vec_entry_one2d( (int) t, (vec_p) n).
Die Zahlen t werden in aufsteigender Reihenfolge im Ausgabevektor abgespeichert. */
vec_p
vec_r_s_n(vec_p r, vec_p n);

/* Gegeben seien die Vektoren r, n mit n>r>=0 und n > 0. Der Vektor n beschreibt ein multidimensionales Array.
Gesucht sind alle ganzen Zahlen t >= 0 fuer die gilt: 0<=s<=r mit s = vec_entry_one2d( (int) t, (vec_p) n).
Die Zahlen t werden in aufsteigender Reihenfolge im Ausgabevektor abgespeichert. */
vec_p
vec_0_s_r(vec_p r, vec_p n);

/* Sei der Vektor n gegeben. Dann wird eine Zerlegung der Form n = 2*s + r berechnet, wobei 0<=r<=1 gilt. Das zurueckgegebene Array
besitzt die Groesse 2. Es gilt Array[0] = s , Array[1] = r. */
vec_p*
vec_zerlegung(vec_p n);

/**
 * Prints out a vector.
 */
void 
print_vec(vec_p vector);

/**
 * Prints out a real-valued vector 
 */
void 
print_vec_real(vec_real_p vector);

/*
 * Returns the euklidian norm of the vector.
 */
fepc_real_t 
vec_real_norm(vec_real_p vector);

#endif
