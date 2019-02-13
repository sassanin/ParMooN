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
   
#include "discont.h"
#include "basictools.h"


linear_function_p
linear_function_new(fepc_real_t y_0, fepc_real_t slope) {
    linear_function_p result = (linear_function_p) malloc(sizeof(linear_function_t));
    
    result->y_0 = y_0;
    result->slope = slope;
    return result;
}


linear_function_p
linear_function_new_points(fepc_real_t y_0, fepc_real_t y_1, fepc_real_t x_0, fepc_real_t x_1) {
	fepc_real_t slope = (y_1-y_0) / (x_1-x_0);

    return linear_function_new(y_0-slope*x_0, slope);
}


linear_function_set_p
linear_function_set_new(int count) {
    linear_function_set_p result = (linear_function_set_p) malloc(sizeof(linear_function_set_t));
    
    result->count = count;
    result->functions = (linear_function_p*) malloc(sizeof(linear_function_p)*count);
    return result;
}


discont_function_p
discont_function_new(int stepcount) {
    discont_function_p result = (discont_function_p) malloc(sizeof(discont_function_t));
    
    result->stepcount = stepcount;
    result->intervals = (interval_p*) malloc(sizeof(interval_p)*stepcount);
    result->function_sets = (linear_function_set_p*) malloc(sizeof(linear_function_set_p)*stepcount);
    return result;
}

void
discont_function_del(discont_function_p function) {
    int n;
    
    for (n = 0; n < function->stepcount; n++) 
	{
        linear_function_set_del(function->function_sets[n]);
    }
    free(function->intervals);
    free(function->function_sets);
    free(function);
}

void
linear_function_set_del(linear_function_set_p function_set) {
    int n;
    

    for (n = 0; n < function_set->count; n++) {
    	if (function_set->functions[n] != NULL) {
    		free(function_set->functions[n]);
    	}
    }
    free(function_set->functions);
    free(function_set);
}

void
discont_function_setup_points(discont_function_p function, int step, fepc_real_t start, fepc_real_t end, fepc_real_t * y1, fepc_real_t * y2, fepc_real_t stepping) {
    function->intervals[step]->start = start;
    function->intervals[step]->end = end;
    
    fepc_real_t h_l = get_h_l(step, stepping);
    
    int count, n;
    
    count = iround((end-start)/h_l); // this has to be equal to the array sizes of y1 and y2

    if (y1 != NULL && y2 != NULL) {
        for (n = 0; n < count; n++) {
			fepc_real_t slope = (y2[n]-y1[n]) / h_l;

            function->function_sets[step]->functions[n]->slope = slope;
			function->function_sets[step]->functions[n]->y_0 = y1[n]-slope*(start+n*h_l);
        }
    } else {
    	for (n = 0; n < count; n++) {
            function->function_sets[step]->functions[n] = NULL;
        }
    }
}


void
discont_function_print(discont_function_p function) {
    int n;
    
    printf("Discontinuous function with %i levels\n=====================================\n", function->stepcount);
    
    for (n = 0; n < function->stepcount; n++) {

        printf("Level %i\n", n);
        linear_function_set_print(function->function_sets[n]);
    }
}


void
linear_function_set_print(linear_function_set_p function_set) {
    int n;
    
    printf("Linear function set\n=====\n");
    for (n = 0; n < function_set->count; n++) {
        printf("   %f*x + %f\n", function_set->functions[n]->slope, function_set->functions[n]->y_0);
    }
}

/*
func_p1_p 
convert_discont_function(discont_function_p function, fepc_real_t stepping) {
    func_p1_p result = func_p1_new(function->stepcount-1);
    
    set_gridstructure(result, function->intervals, stepping); // set the correct intervals
    add_folgenentries_discont(result, function, stepping); // add folgenentries
    return result;
} */


void
convert_func(func_p1_p function, discont_function_p result, fepc_real_t stepping, int* count) {
    fepc_real_t h_l, temp_x1, temp_x2, temp_y1, temp_y2, slope;
    
    int n, k, stepcount;
    
    fepc_real_t x;
    
    for (n = 0; n < result->stepcount; n++) {
        h_l = get_h_l(n, stepping);
        stepcount = iround( (result->intervals[n]->end - result->intervals[n]->start) / h_l );

		for (k = 0; k < stepcount; k++) {
            temp_x1 = (ONE_THIRD+k)*h_l;
            temp_x2 = (TWO_THIRD+k)*h_l;
            x = temp_x1;
            temp_y1 = get_value_at_step(function, x, n, stepping, count);
            x = temp_x2;
            temp_y2 = get_value_at_step(function, x, n, stepping, count);
            slope = 3*(temp_y2-temp_y1)/h_l;
            result->function_sets[n]->functions[k]->slope = slope;
			result->function_sets[n]->functions[k]->y_0 = temp_y1 - slope*temp_x1;
        }
    }
}


fepc_real_t
integrate_coeff_discont(discont_function_p function, int v, int position, int p, int level, fepc_real_t stepping) {
    // integrate from v[0]*h_l till (v[0]+1)*h_l

    fepc_real_t h_l, slope, y_0;
    h_l = get_h_l(level, stepping);

    slope = function->function_sets[level]->functions[position]->slope;
    y_0 = function->function_sets[level]->functions[position]->y_0;

    if (p == 0) { // legendre(0) = sqrt(1/h_l)
        return sqrt(h_l)*(y_0 + h_l*slope*(v+0.5));
    } else { // p == 1 --> legendre(1) = sqrt(12)(x-(v+0.5)*h_l)/(h_l^1.5)
        return h_l*sqrt(h_l)*slope/SQRT_12;
    }
}


void
add_folgenentries_discont(func_p1_p function, discont_function_p discont_function, fepc_real_t stepping, int* count) {

    int level, n, v, lang;

	level = function->maxlevel;
	lang = function->hierarchie[level]->vektor0->lang;

	for (n = 0; n < count[level]; n++) {
		v = n % lang;

		function->hierarchie[level]->vektor0->glied[v] = integrate_coeff_discont(discont_function, v, v, 0, level, stepping);
		function->hierarchie[level]->vektor1->glied[v] = integrate_coeff_discont(discont_function, v, v, 1, level, stepping);
    }

	for (level = 0; level < function->maxlevel; level++) {
		lang = function->hierarchie[level]->vektor0->lang;

		for (n = 0; n < count[level]; n++) {
			v = n % lang;

			if (!is_in_latter_interval(v, function->hierarchie[level+1]->vektor0)) {
				function->hierarchie[level]->vektor0->glied[v] = integrate_coeff_discont(discont_function, v, v, 0, level, stepping);
				function->hierarchie[level]->vektor1->glied[v] = integrate_coeff_discont(discont_function, v, v, 1, level, stepping);
			}
        }
	}
}
