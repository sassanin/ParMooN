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
   
#include "attrition.h"
#include "ifunctions.h"

pattrition new_attrition(int nx, int nr, function_3D3D* velocity)
{
	pattrition a = (pattrition)malloc(sizeof(attrition));
	a->nx = nx;
	a->nr = nr;
	a->s = new_sparsematrix(nr, nr, CMAX);
	a->velocity = velocity;

	return a;
}

double get_x(int ix, int nx)
{
	return (double)ix / (double)nx;
}

double* get_r(int ir, int nr)
{
	double* r = allocate_vector(2);
	r[0] = (double)((int)(ir / nr)) / (double)nr;
	r[1] = (double)(ir - nr * (int)(ir / nr)) / (double)nr;

	return r;
}

int is_on_boundary(double* r)
{
	if(r[0] == 0.0 || r[0] == 1.0 || r[1] == 0.0 || r[1] == 1.0)
		return 1;
	return 0;
}

int apply_attrition(pattrition mo, double* input, double* output, double t)
{
	int ir, ix;
	double x, y, s = 0;
	double* r;
	for(ir = 0; ir < mo->nr * mo->nr; ir++)
	{
		r = get_r(ir, mo->nr);
		for(ix = 0; ix < mo->nx; ix++)
		{
			y = get_x(ix, mo->nx);
			s += nu_frag(y, r, mo->velocity, t) * 2.25 / 
				(pow(L_min, -2.25) - pow(L_max(y, r, mo->velocity, t), -2.25)) * 
				input[ir * mo->nx + ix];
		}
		double db = b(r, mo->velocity, t);
		for(ix = 0; ix < mo->nx; ix++)
		{
			x = get_x(ix, mo->nx);
			output[ir * mo->nx + ix] = pow(x, -3.25) * s + L_star(x, r, mo->velocity, t);
			output[ir * mo->nx + ix] *= db;
		}
		free_vector(r);
	}

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
