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
   
#ifndef __KUERZEN_H
#define __KUERZEN_H

#include "faltung_hilfe.h"

/* Diese Struktur dient der Charakterisierung des Traegers mehrerer Folgenvektoren. Die Struktur dient der Umsetzung der Ergebnisse,
 die in der Dokumentation zum Kapitel "Einkuerzung" gemacht wurden. Ausgehend von einer gegebenen Funktion f des Maxlevel L werden fuer
 0<=l<=L mittels der Arrays a und b Traeger fuer die Folgenvektoren f->hierarchie[l] abgespeichert. Das Intervall [ a[l], b[l] ]
 charakterisiert die Trager der Folgen des Folgenvektors f->hierarchie[l] */

typedef struct {
	vec_p  *a;
	vec_p  *b;
	int  maxlevel;
} support_t;

typedef support_t *  support_p;


/* Algorithmus berechnet die Daten F fuer den Fall B1 und B2 und verkleinert den Traeger
 von F auf den ausschliesslich notwendigen Bereich (siehe Dokumentation). */
func_p
kuerzen_F_bauen_b( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef );


/* Algorithmus berechnet die Daten F fuer den Fall C1 und verkleinert den Traeger
 von F auf den ausschliesslich notwendigen Bereich (siehe Dokumentation). */
func_p
kuerzen_F_bauen_c1( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef );


/* Algorithmus berechnet die Daten F fuer den Fall C2 und verkleinert den Traeger
 von F auf den ausschliesslich notwendigen Bereich (siehe Dokumentation). */
func_p
kuerzen_F_bauen_c2( func_p f, func_p g, func_p w, func2_p Gamma, matrix_p xi_koef );


/* Berechnen des kleinsten gemeinsamen Traegers pro Level von w, der fuer den Fall C notwendig ist */
support_p
support_omega( func_p w );

/* Der Traeger aller Folgen des Folgenvektors f, wird auf das Intervall [a,...,b] eingeschraenkt */
folgen_vektor_p
folgen_vektor_support( folgen_vektor_p f, vec_p a, vec_p b );


/* Initialisierung des Traegers. Man beachte, dass wie bei allen anderen Initalisierungsfunktionen alle Vektoren der Arrays
 ebenfalls initialisiert werden */
support_p
support_new( int maxlevel, int dim );

void
support_del( support_p sup );

#endif
