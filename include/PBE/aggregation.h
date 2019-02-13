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
   
#ifndef MPIMIS_AGGREGATION
#define MPIMIS_AGGREGATION

#include "basictools.h"
#include "piecewise_lin_cont.h"
#include "constants.h"
#include "timer.h"

#include <stdlib.h>
#include <math.h>

typedef class Aggregation aggregation;
typedef aggregation* paggregation;

enum GRID_TYPE {UNIFORM, REFINED};

double u0_brown(double x);
double v0_brown(double y);
double u1_brown(double x);
double v1_brown(double y);

double u0_shear(double x);
double v0_shear(double y);
double u1_shear(double x);
double v1_shear(double y);

class Aggregation {
public:
	// Constructor for uniform grid
	Aggregation(double _xmin, double _xmax, int _na, int _nx, int _ny, int _nz, double* _grid, double* params) {
				xmin = _xmin, xmax = _xmax, na = _na - 1, nx = _nx, ny = _ny, nz = _nz;
				space_dim = iround(params[0]);
				INCLUDE_BROWNIAN_KERNEL = iround(params[1]);
				POLYNOM_DEGREE = iround(params[2]);
				BROWNIAN_FACTOR_TYPE = iround(params[3]);
				BROWNIAN_FACTOR = params[4];
				SHEAR_INDUCED_FACTOR_TYPE = iround(params[5]);
				SHEAR_INDUCED_FACTOR = params[6];
				ETA_C =  params[7];

				a = allocate_vector(na + 1);
				memcpy(a, _grid, (na + 1) * sizeof(double));
				scale_vector(a, na + 1, cube(xmax));

				stepping = a[na] - a[na - 1];

				create_hierarchical_structure();

				generate_piecewise_constant_template(pcf);
				generate_piecewise_constant_template(pcg);
				generate_piecewise_constant_template(pcw);

				generate_piecewise_linear_template(plf);
				generate_piecewise_linear_template(plg);
				generate_piecewise_linear_template(plw);

				generate_piecewise_linear_continuous_template(psd);
				generate_piecewise_linear_continuous_template(rhs);
				generate_piecewise_linear_continuous_template(fc);
				generate_piecewise_linear_continuous_template(gc);
				generate_piecewise_linear_continuous_template(wb0);
				generate_piecewise_linear_continuous_template(wb1);
				generate_piecewise_linear_continuous_template(ws0);
				generate_piecewise_linear_continuous_template(ws1);
				generate_piecewise_linear_continuous_template(wbps);

				generate_func_p0_p_template(pf0);
				generate_func_p0_p_template(pg0);
				generate_func_p0_p_template(pw0);
				generate_func_p0_p_template(convolution_result0);

				u0_brown_v = allocate_vector(na);
				fill_vector_from_function(u0_brown, a, na, u0_brown_v, DS_L2_CONSTANT);
				u1_brown_v = allocate_vector(na);
				fill_vector_from_function(u1_brown, a, na, u1_brown_v, DS_L2_CONSTANT);
				v0_brown_v = allocate_vector(na);
				fill_vector_from_function(v0_brown, a, na, v0_brown_v, DS_L2_CONSTANT);
				v1_brown_v = allocate_vector(na);
				fill_vector_from_function(v1_brown, a, na, v1_brown_v, DS_L2_CONSTANT);
				u0_shear_v = allocate_vector(na);
				fill_vector_from_function(u0_shear, a, na, u0_shear_v, DS_L2_CONSTANT);
				u1_shear_v = allocate_vector(na);
				fill_vector_from_function(u1_shear, a, na, u1_shear_v, DS_L2_CONSTANT);
				v0_shear_v = allocate_vector(na);
				fill_vector_from_function(v0_shear, a, na, v0_shear_v, DS_L2_CONSTANT);
				v1_shear_v = allocate_vector(na);
				fill_vector_from_function(v1_shear, a, na, v1_shear_v, DS_L2_CONSTANT);

				generate_func_p1_p_template(pf1);
				generate_func_p1_p_template(pg1);
				generate_func_p1_p_template(pw1);
				generate_func_p1_p_template(convolution_result1);

				diag = new double[na + 1];
				subdiag = new double[na];
				fill_diagonal_part_of_Gramm_matrix(na, a, diag);
				fill_subdiagonal_part_of_Gramm_matrix(na, a, subdiag);

				fill_l2m_and_m2l_values();
			};

	double* a; 

private:
	GRID_TYPE gt;
	double xmin;
	double xmax;
	int na;
	
	int L; // number of levels in the hierarchy
	double* length; // array of size L. Describes the length of each interval
	int* count; // array of size L. Indicates the number of intervals in the corresponding level
	int* count2; // array of size L. Indicates the number of intervals in the corresponding level
	double** nonzero_value_start;
	double** nonzero_value_end;
	int* nl; // array of size L. Indicates the first grid index of the corresponding level

	double* diag;
	double* subdiag;

	int INCLUDE_BROWNIAN_KERNEL;
	int POLYNOM_DEGREE;
	int BROWNIAN_FACTOR_TYPE;
	double BROWNIAN_FACTOR;
	int SHEAR_INDUCED_FACTOR_TYPE;
	double SHEAR_INDUCED_FACTOR;
	double ETA_C;

public:
	piecewise_constant pcf, pcg, pcw;
	piecewise_linear plf, plg, plw;

	piecewise_linear_continuous psd, rhs, fc, gc;
	piecewise_linear_continuous wb0, wb1, ws0, ws1, wbps;

	interval_p* intervals;
	double* hl;
	double* sqrt_hl;

	func_p0_p pf0, pg0, pw0, convolution_result0;
	func_p1_p pf1, pg1, pw1, convolution_result1;

	double* u0_brown_v;
	double* v0_brown_v;
	double* u1_brown_v;
	double* v1_brown_v;
	double* u0_shear_v;
	double* v0_shear_v;
	double* u1_shear_v;
	double* v1_shear_v;
	double* l2m_start_value_k;
	double* l2m_start_value_b;
	double* l2m_end_value_k;
	double* l2m_end_value_b;
	double* m2l_start_value_k;
	double* m2l_start_value_b;
	double* m2l_end_value_k;
	double* m2l_end_value_b;


	double stepping;
	int nx, ny, nz, space_dim;

public:
	void apply_aggregation(double* input, double* output, double* grad_v, double* temp);

private:
	void generate_piecewise_linear_template(piecewise_linear& pl);
	void generate_piecewise_constant_template(piecewise_constant& pl);
	void generate_piecewise_linear_continuous_template(piecewise_linear_continuous& plc);
	void generate_func_p0_p_template(func_p0_p& fp);
	void generate_func_p1_p_template(func_p1_p& fp);

	void convert_func_p0_to_piecewise_constant(func_p0_p& function, piecewise_constant& pl);
	void convert_pl2func1(piecewise_linear pl, func_p1_p); // l is the level
	void convert_func12pl(func_p1_p fncp1, piecewise_linear plin);
	int create_hierarchical_structure();

	void fill_l2m_and_m2l_values();
	void transform_length2mass(double* input, double* output);
	void transform_mass2length(double* input, double* output);

	// f, g, w are the hidden variables for this
	int rank1_aggregation(piecewise_linear_continuous& r1w);

	int convolution();
};

#endif // MPIMIS_AGGREGATION
