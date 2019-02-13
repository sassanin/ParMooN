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
   
#include "pbe.h"
#include <assert.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>

double entry(int r, int c, void* data)
{
	void** param;
	double *grid, s;
	int n;

	s = 0;
	param = (void**)data;
	grid = (double*)param[1];
	n = *((int*)param[2]);

	if(c != 0 && r != 0)
		s += project_pcw_lnr_fnc2D((function_2D*)param[0], grid[c - 1], grid[c], grid[r - 1], grid[r], 0, 0);
	if(c != 0 && r != (n - 1))
		s += project_pcw_lnr_fnc2D((function_2D*)param[0], grid[c - 1], grid[c], grid[r + 1], grid[r], 0, 0);
	if(c != (n - 1) && r != 0)
		s += project_pcw_lnr_fnc2D((function_2D*)param[0], grid[c + 1], grid[c], grid[r - 1], grid[r], 0, 0);
	if(c != (n - 1) && r != (n - 1))
		s += project_pcw_lnr_fnc2D((function_2D*)param[0], grid[c + 1], grid[c], grid[r + 1], grid[r], 0, 0);

	return s;
}

pintegraloperator new_integraloperator(int nx, int nr)
{
	pintegraloperator m;

	m = (pintegraloperator)malloc(sizeof(integraloperator));
	m->nx = nx;
	m->nr = nr;
	m->xgrid = allocate_vector(nx);
	m->s = new_sparsematrix(nr * nr, nr* nr, CMAX);

	return m;
}

void del_integraloperator(pintegraloperator m)
{
	if(m->rk != 0)
		del_rkmatrix(m->rk);
	if(m->hmatr != 0)
		del_hmatrix(m->hmatr);
	if(m->s != 0)
		del_sparsematrix(m->s);
	free_vector(m->xgrid);
	m->nx = 0;
	m->nr = 0;
}

int load_integraloperator(pintegraloperator mo, double* internal_grid, KernelType kt, 
						  rkdata* lowr, function_2D kl, function_2D kh)
{
	int ret;

	generate_mass_matrix(mo->s);

	memcpy(mo->xgrid, internal_grid, mo->nx * sizeof(double));

	mo->kt = kt;
	if(mo->kt & kt_hmatrix)
	{
		ret = load_hmatrix(mo, kh);
		if(ret != 0)
			return ret;
	}
	return 0;
}

int apply_integraloperator(pintegraloperator mo, double* input, double* output)
{
	int i;
	double* testoutput;

	testoutput = allocate_vector(mo->nx * (mo->nr * mo->nr));
	clear_vector(testoutput, mo->nx * (mo->nr * mo->nr));
	clear_vector(output, mo->nx * (mo->nr * mo->nr));

	if(mo->kt & kt_hmatrix)
	{
		assert(mo->hmatr != 0);
		assert(mo->hmatr->rows == mo->hmatr->cols);

		trihmmm(mo->hmatr, "L", input, mo->nr * mo->nr, output);

		for(i = 0; i < mo->nx * mo->nr * mo->nr; i++)
			testoutput[i] += output[i];
	}
	matrix_times_sparsematrix(mo->s, testoutput, mo->nx, output);

	free_vector(testoutput);

	return 0;
}

int load_hmatrix(pintegraloperator mo, function_2D k)
{
	void** param;

	mo->hmatr = new_hmatrix(mo->nx, mo->nx);
	param = (void**)malloc(3 * sizeof(void*));
	param[0] = k;
	param[1] = mo->xgrid;
	param[2] = &(mo->nx);

	create_hmatrix(mo->hmatr, entry, (void*)param, 0, 0);

	return 0;
}

void generate_mass_matrix(psparsematrix s)
{
	int n = (int)sqrt((double)s->cols);
	double h = 1.0 / (double)(n - 1);
	int c = 0, i;
	int left, up, right, bottom;
	for(i = 0; i < s->cols; i++)
	{
		left = i < n ? 1 : 0;
		up = i % n == 0 ? 1 : 0;
		right = i + n >= n * n ? 1 : 0;
		bottom = i % n == n - 1 ? 1 : 0;

		s->nzindices[c] = (left || up ? -1 : i - n - 1);
		s->data[c++] = h * h / 36;

		s->nzindices[c] = (left ? -1 : i - n);
		s->data[c++] = (up || bottom ? h * h / 18 : h * h / 9);

		s->nzindices[c] = (left || bottom ? -1 : i - n + 1);
		s->data[c++] = h * h / 36;

		s->nzindices[c] = (up ? -1 : i - 1);
		s->data[c++] = (left || right ? h * h / 18 : h * h / 9);

		s->nzindices[c] = i;
		s->data[c++] = h * h * (left || right ? 1 : 2) * (up || bottom ? 1 : 2) / 9;

		s->nzindices[c] = (bottom ? -1 : i + 1);
		s->data[c++] = (left || right ? h * h / 18 : h * h / 9);

		s->nzindices[c] = (up || right ? -1 : i + n - 1);
		s->data[c++] = h * h / 36;

		s->nzindices[c] = (right ? -1 : i + n);
		s->data[c++] = (up || bottom ? h * h / 18 : h * h / 9);

		s->nzindices[c] = (right || bottom ? -1 : i + n + 1);
		s->data[c++] = h * h / 36;
	}
}

