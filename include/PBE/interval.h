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
   
#ifndef __INTERVALS
#define __INTERVALS

#include "basic.h"

typedef struct {
    int dimension;
    fepc_real_t* start;
    fepc_real_t* end;
} interval_t; 

typedef interval_t * interval_p;

/**
 * Creates a new interval with the given dimension.
 */
interval_p 
interval_new(int dimension);

void
interval_del(interval_p interval);

void
intervals_del(interval_p* intervals, int interval_count);

/**
 * Prints out the range of the given interval.
 */
void 
print_interval(interval_p interval);

/**
 * Prints out the ranges of all given intervals.
 */
void 
print_intervals(interval_p * intervals, int count);

interval_p 
interval_clone(interval_p interval);

interval_p *
intervals_clone(interval_p *interval, int count);

#endif
