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
   
#ifndef __FUNKTION_H
#define __FUNKTION_H

#include "folgen_vektor.h"


/************************************* Constant ansatz *******************************************/

typedef struct {
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_vektor_p0_p  *hierarchie;				/* Array von Zeigern auf die einzelnen Folgenvektoren */
} func_p0_t;

typedef func_p0_t *  func_p0_p;

typedef struct {
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_matrix_p0_p  *hierarchie;   		/* Array von Zeigern auf die einzelnen Folgenmatritzen */
} func2_p0_t;

typedef func2_p0_t *  func2_p0_p;

/* Funktion wird initialisiert. Beachte: alle Vektoren von Folgen der Funktion werden durch folgen_vektor_p0_new( vec_new(dim) )
initialisiert. f = func_p0_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Vektoren von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Vektoren von Folgen ist vom grad=vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func_p0_p
func_p0_new(int maxlevel);

func_p0_p
func_p0_new(func_p0_p f, func2_p0_p G, int maxlevel);

void
func_p0_del(func_p0_p f);

/* Funktion2 wird initialisiert. Beachte: alle Matritzen von Folgen der Funktion werden durch
folgen_matrix_p0_new( vec_new(dim), vec_new(dim) ) initialisiert. f = func2_p0_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Matritzen von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Matritzen von Folgen ist vom grad1 = vec_new(dim) und grad2 = vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func2_p0_p
func2_p0_new(int maxlevel);

void
func2_p0_del(func2_p0_p f);


/* Uebertragen der gesamten Folgenwerte der Funtkion f auf eine Funktion mit
 der gleichen Gitterstruktur wie die Funktion g */
void
func_p0_projekt(func_p0_p f,func_p0_p g,func_p0_p back);

/* Gibt die einzelnen Details einer Funktion auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2, 3 sind moegliche Infowerte. 0 wenig Info, 3 viele Infos */
void
func_p0_print(func_p0_p f, int info);

void
func_p0_print_aram(func_p0_p f);


/* Es werden zwei Funktionen addiert */
void
func_p0_add(func_p0_p f, func_p0_p g, func_p0_p back);


/* Entsprechend der Gitterstruktur der hp-Funktion werden gegebenenfalls Folgenwerte gleich Null gesetzt. Falls das zu
 einem Folgenwert zugehoerige Intervall weiter verfeinert wird und somit durch Folgenwerte hoeherer
 Level repraesentiert wird, so wird dieser Folgenwert gleich Null gesetzt. */
void
func_p0_grid_zero(func_p0_p f);

void
func_p0_clean(func_p0_p f);


/************************************* Linear ansatz *******************************************/

typedef struct {
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_vektor_p1_p  *hierarchie;				/* Array von Zeigern auf die einzelnen Folgenvektoren */
} func_p1_t;

typedef func_p1_t *  func_p1_p;


typedef struct {
	int  maxlevel;											/* hoechste Level in Funktionendarstellung */
	folgen_matrix_p1_p  *hierarchie;   		/* Array von Zeigern auf die einzelnen Folgenmatritzen */
} func2_p1_t;

typedef func2_p1_t *  func2_p1_p;



/* Funktion wird initialisiert. Beachte: alle Vektoren von Folgen der Funktion werden durch folgen_vektor_p1_new( vec_new(dim) )
initialisiert. f = func_p1_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Vektoren von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Vektoren von Folgen ist vom grad=vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func_p1_p
func_p1_new(int maxlevel);

func_p1_p
func_p1_new(func_p1_p f, func2_p1_p G, int maxlevel);

void
func_p1_del(func_p1_p f);

/* Funktion2 wird initialisiert. Beachte: alle Matritzen von Folgen der Funktion werden durch
folgen_matrix_p1_new( vec_new(dim), vec_new(dim) ) initialisiert. f = func2_p1_new(maxlevel) bedeutet:
f besteht aus maxlevel+1 Matritzen von Folgen f->hierarchie[k] fuer k=0,...,maxlevel.
Jeder dieser Matritzen von Folgen ist vom grad1 = vec_new(dim) und grad2 = vec_new(dim) ,dh. besteht aus einer trivialen Folge.
Somit ist die neu initialisierte Funktionen volldefiniert und vollwertig. */
func2_p1_p
func2_p1_new(int maxlevel);

void
func2_p1_del(func2_p1_p f);


/* Uebertragen der gesamten Folgenwerte der Funtkion f auf eine Funktion mit
 der gleichen Gitterstruktur wie die Funktion g */
void
func_p1_projekt(func_p1_p f,func_p1_p g,func_p1_p back);

/* Gibt die einzelnen Details einer Funktion auf dem Bildschirm wieder. Je groesser der Wert info ist, umso mehr
 Informationen werden wiedergegeben. 0, 1, 2, 3 sind moegliche Infowerte. 0 wenig Info, 3 viele Infos */
void
func_p1_print(func_p1_p f, int info);

void
func_p1_print_aram(func_p1_p f);


/* Es werden zwei Funktionen addiert */
void
func_p1_add(func_p1_p f, func_p1_p g, func_p1_p back);


/* Entsprechend der Gitterstruktur der hp-Funktion werden gegebenenfalls Folgenwerte gleich Null gesetzt. Falls das zu
 einem Folgenwert zugehoerige Intervall weiter verfeinert wird und somit durch Folgenwerte hoeherer
 Level repraesentiert wird, so wird dieser Folgenwert gleich Null gesetzt. */
void
func_p1_grid_zero(func_p1_p f);

void
func_p1_clean(func_p1_p f);

#endif
