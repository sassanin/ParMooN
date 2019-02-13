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
   
#ifndef __folgen_vektor_p0_H
#define __folgen_vektor_p0_H

#include "folge.h"

/************************************* Constant ansatz *******************************************/

typedef struct {
	folge_p0_p vektor;
} folgen_vektor_p0_t;

typedef folgen_vektor_p0_t *  folgen_vektor_p0_p;

typedef struct {
	folge_p0_p  matrix;
} folgen_matrix_p0_t;

typedef folgen_matrix_p0_t *  folgen_matrix_p0_p;

folgen_vektor_p0_p
folgen_vektor_p0_new();

void
folgen_vektor_p0_del(folgen_vektor_p0_p f);

/*folgen_matrix wird initialisiert. Beachte: alle Folgen der Matrix werden durch folge_new(r,r) initialisiert, wobei
 r = vec_new(dim) ist. Das heisst alle Folgen der Matrix haben die Laenge 0 (siehe Implementierung von folge_new)*/
folgen_matrix_p0_p
folgen_matrix_p0_new();

void
folgen_matrix_p0_del(folgen_matrix_p0_p f);

/* Gibt die einzelnen Details eines Folgenvektor auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2 sind moegliche Infowerte. 0 wenig Info, 2 viele Infos */
void
folgen_vektor_p0_print(folgen_vektor_p0_p f, int info);

void
folgen_vektor_p0_print_aram(folgen_vektor_p0_p f);

/* Berechnung der erweiterten Version der diskreten Faltung. Faltung von Tupeln von Folgen siehe Dokumentation */
void
folgen_vektor_p0_faltung(folgen_vektor_p0_p f, folgen_matrix_p0_p Gamma, folgen_vektor_p0_p back);

/* Es wird ein Vektor von Folgen erstellt, der die gleiche Struktur wie der Eingabevektor besitzt. Alle Folgen des
 Vektors haben als Folgenglied den Wert 0.0 */
void
folgen_vektor_p0_copy_structure(folgen_vektor_p0_p w, folgen_vektor_p0_p back);

/* Uebertragen der gesamten Folgenwerte des Folgenvektors f auf einen Folgenvektor mit
 der gleichen Gitterstruktur wie der Folgenvektor g */
void
folgen_vektor_p0_projekt(folgen_vektor_p0_p f,folgen_vektor_p0_p g,folgen_vektor_p0_p back);

/* Addition von Folgenvektor f und g */
void
folgen_vektor_p0_add(folgen_vektor_p0_p f, folgen_vektor_p0_p g, folgen_vektor_p0_p back);

/* Addition von Folgenmatritzen f und g */
void
folgen_matrix_p0_add(folgen_matrix_p0_p f, folgen_matrix_p0_p g, folgen_matrix_p0_p back);

void
folgen_matrix_p0_print(folgen_matrix_p0_p m);


/************************************* Linear ansatz *******************************************/

typedef struct {
	folge_p1_p vektor0;
	folge_p1_p vektor1;
} folgen_vektor_p1_t;

typedef folgen_vektor_p1_t *  folgen_vektor_p1_p;

typedef struct {
	folge_p1_p  **matrix;
} folgen_matrix_p1_t;

typedef folgen_matrix_p1_t *  folgen_matrix_p1_p;

folgen_vektor_p1_p
folgen_vektor_p1_new();

void
folgen_vektor_p1_del(folgen_vektor_p1_p f);


/*folgen_matrix wird initialisiert. Beachte: alle Folgen der Matrix werden durch folge_p1_new(r,r) initialisiert, wobei
 r = vec_new(dim) ist. Das heisst alle Folgen der Matrix haben die Laenge 0 (siehe Implementierung von folge_p1_new)*/
folgen_matrix_p1_p
folgen_matrix_p1_new();

void
folgen_matrix_p1_del(folgen_matrix_p1_p f);

/* Gibt die einzelnen Details eines Folgenvektor auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2 sind moegliche Infowerte. 0 wenig Info, 2 viele Infos */
void
folgen_vektor_p1_print(folgen_vektor_p1_p f, int info);

void
folgen_vektor_p1_print_aram(folgen_vektor_p1_p f);

/* Berechnung der erweiterten Version der diskreten Faltung. Faltung von Tupeln von Folgen siehe Dokumentation */
void
folgen_vektor_p1_faltung(folgen_vektor_p1_p f, folgen_matrix_p1_p Gamma, folgen_vektor_p1_p back);

/* Es wird ein Vektor von Folgen erstellt, der die gleiche Struktur wie der Eingabevektor besitzt. Alle Folgen des
 Vektors haben als Folgenglied den Wert 0.0 */
void
folgen_vektor_p1_copy_structure(folgen_vektor_p1_p w, folgen_vektor_p1_p back);

/* Uebertragen der gesamten Folgenwerte des Folgenvektors f auf einen Folgenvektor mit
 der gleichen Gitterstruktur wie der Folgenvektor g */
void
folgen_vektor_p1_projekt(folgen_vektor_p1_p f,folgen_vektor_p1_p g,folgen_vektor_p1_p back);

/* Addition von Folgenvektor f und g */
void
folgen_vektor_p1_add(folgen_vektor_p1_p f, folgen_vektor_p1_p g, folgen_vektor_p1_p back);

/* Addition von Folgenmatritzen f und g */
void
folgen_matrix_p1_add(folgen_matrix_p1_p f, folgen_matrix_p1_p g, folgen_matrix_p1_p back);

void
folgen_matrix_p1_print(folgen_matrix_p1_p m);


#endif

