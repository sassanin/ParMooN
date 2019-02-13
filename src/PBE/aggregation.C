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
   
#include <aggregation.h>
#include "exported.h"
//#include "MooNMD_Io.h"

static double factor_shear = 0;
static double factor_brown = 0;

double u0_brown(double x) { return 1.0; };
double v0_brown(double y) { return 1.0; };
double u1_brown(double x) { return cubert(x); };
double v1_brown(double y) { return y == 0 ? 0 : 1.0 / (cubert(y)); };

double u0_shear(double x) { return x; };
double v0_shear(double y) { return 1.0; };
double u1_shear(double x) { return 3 * cubert(sqr(x)); };
double v1_shear(double y) { return cubert(y); };

void Aggregation::fill_l2m_and_m2l_values()
{
	double h;
	l2m_start_value_k = new double[na];
	l2m_end_value_k = new double[na];
	l2m_start_value_b = new double[na];
	l2m_end_value_b = new double[na];
	m2l_start_value_k = new double[na];
	m2l_end_value_k = new double[na];
	m2l_start_value_b = new double[na];
	m2l_end_value_b = new double[na];

	for(int i = 0; i < na; i++) {
		h = a[i+1] - a[i];
		l2m_start_value_k[i] = (pow(a[i+1],2./3.) - pow(a[i], 2./3.)) * 1.5 / h - 6 / cube(h) * h * 
							( (pow(a[i+1], 5./3.) - pow(a[i], 5./3.)) * 0.6 - 
								  0.75 * (a[i+1] + a[i]) * ( pow(a[i+1], 2./3.) - pow(a[i], 2./3.) ) );
		l2m_end_value_k[i] = (pow(a[i+1], 2./3.) - pow(a[i], 2./3.)) * 1.5 / h + 6 / cube(h) * h * 
							( (pow(a[i+1], 5./3.) - pow(a[i], 5./3.)) * 0.6 - 
							   0.75 * (a[i+1] + a[i]) * ( pow(a[i+1], 2./3.) - pow(a[i], 2./3.) ) );
		l2m_start_value_b[i] = (pow(a[i+1], 1./3.) - pow(a[i], 1./3.)) * 3 / h - 6 / cube(h) * h * 
							( (pow(a[i+1], 4./3.) - pow(a[i], 4./3.)) * 0.75 - 
							   1.5 * (a[i+1] + a[i]) * ( pow(a[i+1], 1./3.) - pow(a[i], 1./3.) ) );
		l2m_end_value_b[i] = (pow(a[i+1], 1./3.) - pow(a[i], 1./3.)) * 3 / h + 6 / cube(h) * h * 
							( (pow(a[i+1], 4./3.) - pow(a[i], 4./3.)) * 0.75 - 
							   1.5 * (a[i+1] + a[i]) * ( pow(a[i+1], 1./3.) - pow(a[i], 1./3.) ) );

		h = cubert(a[i+1]) - cubert(a[i]);
		m2l_start_value_k[i] = (sqr(a[i+1]) - sqr(a[i])) / (6 * h) - 6 / cube(h) * h * 
							( (pow(a[i+1], 7./3.) - pow(a[i], 7./3.)) / 7. - 
							  (cubert(a[i+1]) + cubert(a[i])) * ( sqr(a[i+1]) - sqr(a[i])) / 12. );
		m2l_end_value_k[i] = (pow(a[i+1], 2.) - pow(a[i], 2.)) / (6 * h) + 6 / cube(h) * h * 
							( (pow(a[i+1], 7./3.) - pow(a[i], 7./3.)) * 1./7. - 
							  (cubert(a[i+1]) + cubert(a[i])) * ( pow(a[i+1], 2.) - pow(a[i], 2.)) / 12. );
		m2l_start_value_b[i] = (pow(a[i+1], 1.) - pow(a[i], 1.)) / (3 * h) - 6 / cube(h) * h * 
							( (pow(a[i+1], 4./3.) - pow(a[i], 4./3.)) * 1./4. - 
							  (cubert(a[i+1]) + cubert(a[i])) * (a[i+1] - a[i]) / 6. );
		m2l_end_value_b[i] = (pow(a[i+1], 1.) - pow(a[i], 1.)) / (3 * h) + 6 / cube(h) * h * 
							( (pow(a[i+1], 4./3.) - pow(a[i], 4./3.)) * 1./4. - 
							  (cubert(a[i+1]) + cubert(a[i])) * (a[i+1] - a[i]) / 6. );
	}
}

void Aggregation::transform_length2mass(double* input, double* output)
{
	// plf - projection of kx / 3 * cubert(x^2)
	// plg - projection of b / 3 * cubert(x^2)
	double k, b;
	for(int i = 0; i < na; i++)
	{
		k = (input[i+1] - input[i]) / (cubert(a[i+1]) - cubert(a[i]));
		b = input[i] - k * cubert(a[i]);
		k /= 3; b /= 3;

		plf.start_value[i] = l2m_start_value_k[i] * k;
		plf.end_value[i]   = l2m_end_value_k[i] * k;
		plg.start_value[i] = l2m_start_value_b[i] * b;
		plg.end_value[i]   = l2m_end_value_b[i] * b;

	}
	output[0] = 0;
	output[na] = plf.end_value[na-1] + plg.end_value[na-1];
	if(output[na] < 0) output[na] = 0;
	for(int i = 1; i < na; i++)
	{
		output[i] = 0.5 * ( (plf.start_value[i] + plg.start_value[i]) + (plf.end_value[i-1] + plg.end_value[i-1]) );

		if(output[i] < 0) output[i] = 0;
	}
}

void Aggregation::transform_mass2length(double* input, double* output)
{
	// plf - projection of k * l^3 * 3 * l^2
	// plg - projection of b * 3 * l^2
	double k, b;
	for(int i = 0; i < na; i++)
	{
		k = (input[i+1] - input[i]) / (a[i+1] - a[i]);
		b = input[i] - k * a[i];
		k *= 3; b *= 3;
		plf.start_value[i] = m2l_start_value_k[i] * k;
		plf.end_value[i]   = m2l_end_value_k[i] * k;
		plg.start_value[i] = m2l_start_value_b[i] * b;
		plg.end_value[i]   = m2l_end_value_b[i] * b;
	}
	output[0] = 0;
	output[na] = plf.end_value[na-1] + plg.end_value[na-1];
	if(output[na] < 0) output[na] = 0;
	for(int i = 1; i < na; i++)
	{
		output[i] = 0.5 * ( (plf.start_value[i] + plg.start_value[i]) + (plf.end_value[i-1] + plg.end_value[i-1]) );
	}
}

void Aggregation::apply_aggregation(double* input, double* output, double* grad_v, double* temp) {
	int Nr = nx * ny * nz, i, j;
	double input_mass = 0.0, rhs_mass = 0.0;

	for(i = 0; i < Nr; i++)
	{
		for(j = 0; j <= na; j++) {
			psd.coeff[j] = input[i + j * Nr]; // a[j] > 0 ? input[i + j * Nr] / (3 * cubert(sqr(a[j]))) : 0;
		}
		transform_length2mass(psd.coeff, psd.coeff);

/*		if(psd.mass() > 2) {
			printf("Mass violation!\n");
			exit(0);
		} */

// 		input_mass += psd.mass();
// 		printf("Input mass is %g\n", input_mass);

		factor_brown = BROWNIAN_FACTOR * 
					  (BROWNIAN_FACTOR_TYPE == 0 ? 1 : 
					   2 * k_B * temp[i] / (3 * ETA_C));
// 					  printf("eta_c is %g\n", ETA_C);
		factor_shear = SHEAR_INDUCED_FACTOR * 
					  (SHEAR_INDUCED_FACTOR_TYPE  == 0 ? 1 : 
					   SQRT2 * euclidean_norm(grad_v + space_dim * space_dim * i, space_dim * space_dim));
 
// 					  /* 
// 					  if(fabs(factor_brown)>1e-6)
// 					  printf("brown, shear is %g  %g\n",factor_brown, factor_shear);
// printf("shear is %g  %g\n",factor_shear, 
// euclidean_norm(grad_v + space_dim * space_dim * i, space_dim * space_dim) );
// */
		if(POLYNOM_DEGREE == 0)
		{
			if(INCLUDE_BROWNIAN_KERNEL) 
			{
				psd.multiplication_with_function(u0_brown_v, pcf);
				psd.multiplication_with_function(v0_brown_v, pcg);
				rank1_aggregation(wb0);

				psd.multiplication_with_function(u1_brown_v, pcf);
				psd.multiplication_with_function(v1_brown_v, pcg);
				rank1_aggregation(wb1);
			}

			psd.multiplication_with_function(u0_shear_v, pcf);
			psd.multiplication_with_function(v0_shear_v, pcg);
			rank1_aggregation(ws0);

			psd.multiplication_with_function(u1_shear_v, pcf);
			psd.multiplication_with_function(v1_shear_v, pcg);
			rank1_aggregation(ws1);
 		}
		else if(POLYNOM_DEGREE == 1)
		{
			if(INCLUDE_BROWNIAN_KERNEL) 
			{
				psd.multiplication_with_function(u0_brown_v, plf);
				psd.multiplication_with_function(v0_brown_v, plg);
				rank1_aggregation(wb0);

				psd.multiplication_with_function(u1_brown_v, plf);
				psd.multiplication_with_function(v1_brown_v, plg);
				rank1_aggregation(wb1);
			}

			psd.multiplication_with_function(u0_shear_v, plf);
			psd.multiplication_with_function(v0_shear_v, plg);
			rank1_aggregation(ws0);

			psd.multiplication_with_function(u1_shear_v, plf);
			psd.multiplication_with_function(v1_shear_v, plg);
			rank1_aggregation(ws1);
		}

		for(j = 0; j <= na; j++)
		{
			wbps.coeff[j] = factor_shear * (ws0.coeff[j] + ws1.coeff[j]) +
							(INCLUDE_BROWNIAN_KERNEL == 0 ? 0 : factor_brown * (wb0.coeff[j] + wb1.coeff[j]));
							
//    if (fabs((wb0.coeff[j] + wb1.coeff[j]))>1.e-3)
//    {
//  printf("INCLUDE_BROWNIAN_KERNEL, shear is %d, %d, %f  %f\n", POLYNOM_DEGREE,  INCLUDE_BROWNIAN_KERNEL, wb0.coeff[j], wb1.coeff[j]);
// //  exit(0);
//    }
	}
//		wbps.ProjectToMassConservedShift();
//		rhs_mass += wbps.mass();

/*		for(j = 0; j <= na; j++)
			output[i + j * Nr] = wbps.coeff[j] * 3 * cubert(sqr(a[j])); */
		transform_mass2length(wbps.coeff, wbps.coeff);
		for(j = 0; j <= na; j++)
			output[i + j * Nr] = wbps.coeff[j];
	}
// 	fprint_vector("output1.dat", output, Nr * (na+1));
// exit(0);

// 	OutPut("Total mass of the psd was " << input_mass << endl << "Total mass of the right-hand side was " << rhs_mass << endl);
	
/*	global_timer.stop_timer();
	projection_timer.stop_timer();
	hp_faltung_timer.stop_timer();
	discrete_faltung_timer.stop_timer();

	printf("%d aggregations with %d gridpoints were done in %g seconds.\n", Nr, na, global_timer.get_time());
	printf("hp-convolution without discrete convolution took %g seconds or %g %%.\n", hp_faltung_timer.get_time(), 100.0 * hp_faltung_timer.get_time() / global_timer.get_time());
	printf("Discrete convolution took %g seconds or %g %%.\n", discrete_faltung_timer.get_time(), 100.0 * discrete_faltung_timer.get_time() / global_timer.get_time());
	printf("Projections and convertions took %g seconds or %g %%.\n", projection_timer.get_time(), 100.0 * projection_timer.get_time() / global_timer.get_time()); */
}

int aggregation::create_hierarchical_structure()
{
	double h_min = a[1] - a[0];
	double h_max = a[na] - a[na - 1];
	L = iround(log(h_max / h_min, 2.0)) + 1;
	length = new double[L];
	hl = new double[L];
	sqrt_hl = new double[L];
	count = new int[L];
	count2 = new int[L];
	nl = new int[L];
	nonzero_value_start = new double*[L];
	nonzero_value_end = new double*[L];
	intervals = (interval_p*) malloc(sizeof(interval_p)*L);

	extern int MAX_POSSIBLE_FOLGE_LENGTH;

	double h = h_min;
	int i, l = L-1;
	for(i = 0; i < na; i++)
	{
		if(iround((a[i + 1] - a[i]) / h) > 1)
		{
			length[l] = a[i] - a[0];
			nl[l-1] = i;
			count[l] = iround(length[l] / h);
			count2[l] = (l < (L - 1) ? iround((length[l] - length[l + 1]) / h)  : count[l]);
			nonzero_value_start[l] = new double[count[l]];
			nonzero_value_end[l] = new double[count[l]];

			if(2 * count[l] > MAX_POSSIBLE_FOLGE_LENGTH)
				MAX_POSSIBLE_FOLGE_LENGTH = 2 * count[l] + 1;

			i--; l--; h *= 2.0;
		}
	}
	if(L == 1)
		MAX_POSSIBLE_FOLGE_LENGTH = 2 * na;

	hl[L-1] = a[1] - a[0];
	nl[L-1] = 0;
	length[0] = a[na] - a[0];
	count[0] = iround(length[0] / h_max);
	count2[0] = (0 < (L - 1) ? iround((length[0] - length[1]) / h)  : count[0]);
	nonzero_value_start[0] = new double[count[0]];
	nonzero_value_end[0] = new double[count[0]];
	assert(l == 0);

	for(int l = 0; l < L; l++)
	{
		intervals[l] = new interval_t;

		intervals[l]->start = a[0];
		intervals[l]->end = a[0] + length[l];

		hl[l] = hl[L-1] * pow(2., L-l-1);
		sqrt_hl[l] = sqrt(hl[l]);
	}

	return 0;
}

int aggregation::rank1_aggregation(piecewise_linear_continuous& r1w)
{
	int n = pcw.n, i;
	convolution();

	double Integral_of_g;
	double Integral_of_f;

	if(POLYNOM_DEGREE == 0)
	{
		Integral_of_g = pcg.integral();
		Integral_of_f = pcf.integral();

		for(i = 0; i < n; i++)
			pcw.value[i] -= pcf.value[i] * Integral_of_g + pcg.value[i] * Integral_of_f;

		r1w.IsContinuousProjectionOf(pcw, diag, subdiag);

//		r1w.IsContinuousProjectionOf(pcw, diag, subdiag);
//		fc.IsContinuousProjectionOf(pcf, diag, subdiag);
//		gc.IsContinuousProjectionOf(pcg, diag, subdiag);
	}
	else if(POLYNOM_DEGREE == 1)
	{
		Integral_of_g = plg.integral();
		Integral_of_f = plf.integral();

		for(i = 0; i <= n; i++) {
			if(plw.start_value[i] < 0) plw.start_value[i] = 0;
			if(plw.end_value[i] < 0) plw.end_value[i] = 0;

			plw.start_value[i] -= plf.start_value[i] * Integral_of_g + plg.start_value[i] * Integral_of_f;
			plw.end_value[i] -= plf.end_value[i] * Integral_of_g + plg.end_value[i] * Integral_of_f;
		}
		r1w.IsContinuousProjectionOf(plw, diag, subdiag);

//		r1w.IsContinuousProjectionOf(plw, diag, subdiag);
//		fc.IsContinuousProjectionOf(plf, diag, subdiag);
//		gc.IsContinuousProjectionOf(plg, diag, subdiag);
	}

//	double Integral_of_g = gc.integral();
//	double Integral_of_f = fc.integral();
//	for(i = 0; i <= n; i++)
//	{
//		r1w.coeff[i] = r1w.coeff[i] - fc.coeff[i] * Integral_of_g - gc.coeff[i] * Integral_of_f;
//	}

	return 0;
}

int aggregation::convolution()
{
	int l, i, start, end;
	if(POLYNOM_DEGREE == 0)
	{
		for(l = 0; l < L; l++)
		{
			end = count[l] - 1;
			start = count[l] - count2[l];
			for(i = start - 1; i>=0; i--) {
				pf0->hierarchie[l]->vektor->glied[i] = 0;
				pg0->hierarchie[l]->vektor->glied[i] = 0;
			}

			for(i = end; i >= start; i--) {
				pf0->hierarchie[l]->vektor->glied[i] = sqrt_hl[l] * pcf.value[nl[l] + i - start];
				pg0->hierarchie[l]->vektor->glied[i] = sqrt_hl[l] * pcg.value[nl[l] + i - start];
			}
		}

	//	projection_timer.pause_timer();
	//	hp_faltung_timer.resume_timer();

		faltung_fepc(pf0, pg0, pw0, stepping, convolution_result0);

	//	hp_faltung_timer.pause_timer();
	//	projection_timer.resume_timer();

		convert_func_p0_to_piecewise_constant(convolution_result0, pcw);
	//	print_vector(pcw.value, na);
	}
	else if(POLYNOM_DEGREE == 1)
	{
		convert_pl2func1(plf, pf1);
		convert_pl2func1(plg, pg1);

		faltung_fepc(pf1, pg1, pw1, stepping, convolution_result1);

		convert_func12pl(convolution_result1, plw);
	}

	return 0;
}

void aggregation::convert_pl2func1(piecewise_linear plin, func_p1_p fncp1)
{
	int i, start, end;
	double dif, y, h_l, sqrt_h_l;

	for(int l = 0; l < L; l++)
	{
		end = count[l] - 1; start = count[l] - count2[l];
		h_l = hl[l]; /*get_h_l(l, stepping);*/ sqrt_h_l = sqrt_hl[l];

		for(i = start - 1; i>=0; i--) {
			nonzero_value_start[l][i] = 0;
			nonzero_value_end[l][i] = 0;
		}

		for(i = start; i <= end; i++) {
			nonzero_value_start[l][i] = plin.start_value[nl[l] + (i - start)];
			nonzero_value_end[l][i] = plin.end_value[nl[l] + (i - start)];

			dif = (nonzero_value_end[l][i] - nonzero_value_start[l][i]);
			y = nonzero_value_start[l][i] - dif * (a[0] + i * h_l) / h_l;

			if (l == L-1 || !is_in_latter_interval(i, pf1->hierarchie[l+1]->vektor0)) {
				fncp1->hierarchie[l]->vektor0->glied[i] = sqrt_h_l * (y + dif * (i + 0.5));
				fncp1->hierarchie[l]->vektor1->glied[i] = sqrt_h_l * dif / SQRT_12;
			}
		}
	}
}
void aggregation::convert_func12pl(func_p1_p fncp1, piecewise_linear plin)
{
	double x1, x2, slope, y1, y2, h_l, y_0;
	for (int l = 0; l < L; l++) {
		h_l = hl[l];

		for (int i = count[l] - count2[l]; i < count[l]; i++) {
			x1 = (ONE_THIRD+i)*h_l;
			y1 = get_value_at_step(fncp1, x1, l, stepping, count);
			x2 = (TWO_THIRD+i)*h_l;
			y2 = get_value_at_step(fncp1, x2, l, stepping, count);
			slope = 3 * (y2 - y1) / h_l;
			y_0 = y1 - slope * x1;

			plin.start_value[nl[l] + i - (count[l] - count2[l])] = y_0 + slope * a[nl[l] + i - (count[l] - count2[l])];

			plin.end_value[nl[l] + i - (count[l] - count2[l])] = y_0 + slope * a[nl[l] + i - (count[l] - count2[l]) + 1];
		}
	}
}

void aggregation::generate_piecewise_constant_template(piecewise_constant& pl)
{
	pl.n = na;
	pl.grid = a;
	pl.value = new double[na];
}
void aggregation::generate_piecewise_linear_template(piecewise_linear& pl)
{
	pl.n = na;
	pl.grid = a;
	pl.start_value = new double[na + 1];
	pl.end_value = new double[na + 1];
}
void aggregation::generate_piecewise_linear_continuous_template(piecewise_linear_continuous& plc)
{
	plc.n = na;
	plc.grid = a;
	plc.coeff = new double[na + 1];
}
void aggregation::generate_func_p0_p_template(func_p0_p& fp)
{
	int l;
	fp = func_p0_new(L-1);

    for (l = 0; l < L; l++) {
        folge_del(fp->hierarchie[l]->vektor);
		fp->hierarchie[l]->vektor = folge_new(count[l]);
    }
}
void aggregation::generate_func_p1_p_template(func_p1_p& fp)
{
	fp = func_p1_new(L-1);
	set_gridstructure(fp, intervals, stepping); // set the correct intervals
}
void aggregation::convert_func_p0_to_piecewise_constant(func_p0_p& function, piecewise_constant& pl)
{
	for(int l = 0; l < L; l++)
	{
		double hl = stepping * pow(0.5, l);
		for(int i = count[l] - count2[l]; i < count[l]; i++)
		{
			pl.value[nl[l] + i - (count[l] - count2[l])] = function->hierarchie[l]->vektor->glied[i] * sqrt(1./hl);
		}
	}
}
