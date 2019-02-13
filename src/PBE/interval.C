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
   
#include "interval.h"


interval_p 
interval_new(int dimension) {
    interval_p result = (interval_p) malloc(sizeof(interval_t));
    
    result->dimension = dimension;
    result->start = (fepc_real_t*) malloc(sizeof(fepc_real_t)*dimension);
    result->end = (fepc_real_t*) malloc(sizeof(fepc_real_t)*dimension);
    
    int n;
    
    for (n = 0; n < dimension; n++) {
        result->start[n] = 0.;
    }
    return result;    
}

void interval_del(interval_p interval) {
    free(interval->start);
    free(interval->end);
    free(interval);
}

void intervals_del(interval_p* intervals, int interval_count) {
    int n;
    
    for (n = 0; n < interval_count; n++) {
        interval_del(intervals[n]);
    }
    free(intervals);
}

void print_interval(interval_p interval) {
    int n;
    
    for (n = 0; n < interval->dimension; n++) {
        printf("%f\t<\t%f\n", interval->start[n], interval->end[n]);
    } 
}

void print_intervals(interval_p * intervals, int count) {
    int n;
    
    for (n = 0; n < count; n++) {
        printf("\nInterval %i\n", n+1);
        print_interval(intervals[n]);
    }
}

interval_p interval_clone(interval_p interval) {
    int n, dimension;
    
    dimension = interval->dimension;
    interval_p result = interval_new(dimension);
    for (n = 0; n < dimension; n++) {
        result->start[n] = interval->start[n];
        result->end[n] = interval->end[n];
    }
    return result;
}

interval_p * intervals_clone(interval_p* intervals, int count) {
    interval_p* result = (interval_p*) malloc(sizeof(interval_p)*count);
    
    int n;
    
    for (n = 0; n < count; n++) {
        result[n] = interval_clone(intervals[n]);
    }
    return result;
}
