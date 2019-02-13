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
   
#ifndef __FEPCEASY
#define __FEPCEASY

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "faltung.h"

typedef struct {
    fepc_real_t start;
    fepc_real_t end;
} interval_t; 

typedef interval_t * interval_p;

#ifndef INT_STEPS
#define INT_STEPS 10
#endif

#define SQRT_3 1.7320508075688772935274463415058723669428
#define DIV_1_6 0.16666666666666666666666666666666666666666666

typedef fepc_real_t (*Funcimpl_step) (int, int, int, fepc_real_t);

typedef fepc_real_t (*Funcimpl) (fepc_real_t);

/**
 * Sets the gridstructure in the needed form out of a common list of intervals.
 */
void set_gridstructure(func_p0_p function, interval_p* intervals, fepc_real_t stepping);
void set_gridstructure(func_p1_p function, interval_p* intervals, fepc_real_t stepping);

/**
 * Returns the h_l value corrensponding to the step.
 */
inline fepc_real_t get_h_l(int step, fepc_real_t& stepping) {
    return pow(.5, step)*stepping;
}

/**
 * Returns true if the vector is in the latter interval of the folge.
 */
inline bool_t is_in_latter_interval(int& v, folge_p0_p folge) {
	return v >= folge->lang / 2. ? FEPC_FALSE : FEPC_TRUE;
}

fepc_real_t get_value_at_step(func_p0_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count);

/**
 * Returns true if the vector is in the latter interval of the folge.
 */
inline bool_t is_in_latter_interval(int& v, folge_p1_p folge) {
	return v >= folge->lang / 2. ? FEPC_FALSE : FEPC_TRUE;
}

fepc_real_t get_value_at_step(func_p1_p function, fepc_real_t& x, int& step, fepc_real_t& stepping, int* count);

/**
 * Returns the value for phi_l (basis function) for the given vector, point x and step.
 */
fepc_real_t phi_l(int& step, int& v, int p, fepc_real_t& x, fepc_real_t& stepping);

#endif
