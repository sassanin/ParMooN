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
   
#include "basic.h"

int int_round(fepc_real_t d) {
	return (int)(d + 0.5);
}

vec_p
vec_new(int dim) {
	vec_p  back;
	int  k;

	back = (vec_p) malloc(sizeof(vec_t));
	ASSERT(back != NULL);
	ASSERT(dim > 0);
	back->array = (int*) malloc(sizeof(int)*dim);
	ASSERT(back->array != NULL);
	for(k=0;k<dim;k++) {
		back->array[k] = 0;
	}
	back->dim = dim;

	return back;
}

vec_real_p
vec_real_new(int dim) {
	vec_real_p  back;
	int  k;

	back = (vec_real_p) malloc(sizeof(vec_real_t));
	ASSERT(back != NULL);
	ASSERT(dim > 0);
	back->array = (fepc_real_t*) malloc(sizeof(fepc_real_t)*dim);
	ASSERT(back->array != NULL);
	for(k=0;k<dim;k++) {
		back->array[k] = 1.0;
	}
	back->dim = dim;

	return back;
}


vec_p
vec_one(int dim) {
	vec_p  back;
	int  k;

	back = (vec_p) malloc(sizeof(vec_t));
	ASSERT(back != NULL);
	ASSERT(dim > 0);
	back->array = (int*) malloc(sizeof(int)*dim);
	ASSERT(back->array != NULL);
	for(k=0;k<dim;k++) {
		back->array[k] = 1;
	}
	back->dim = dim;

	return back;
}


void
vec_del(vec_p s) {
	free(s->array);
	free(s);
}

void
vec_real_del(vec_real_p s) {
	free(s->array);
	free(s);
}

int
entry_d2one(vec_p s,vec_p n) {
	int  i, l, sum, prod, test;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT(s->dim == n->dim);

/*	test = 0;
	for(i=0;i<s->dim;i++) {
		if ((s->array[i]>=n->array[i]) || (s->array[i]<0)) {
			test = test + 1;
		}
		if ((n->array[i]<=0)) {
			test = test + 1;
		}
	}
	ASSERT(test == 0);*/

	/*Berechnung*/
	sum = s->array[s->dim-1];

	for(i=0;i<s->dim-1;i++) {
		prod = 1;
		for(l=i+1;l<s->dim;l++) {
			prod = prod * n->array[l];
		}
		sum = sum + s->array[i]*prod;
	}

	return sum;
}


vec_p
entry_one2d(int pos,vec_p n) {
	int  k, size, test;
	vec_p  s;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	for(k=0;k<n->dim;k++) {
		ASSERT( n->array[k] > 0 );
	}


	size = vec_size( n );
	ASSERT(pos >= 0);
	ASSERT(pos < size);

	/*Berechnung*/
	s = vec_new( n->dim );

	for(k=s->dim-1;k>=1;k--) {
		s->array[k] = pos % n->array[k];
		pos = (pos - s->array[k]) / n->array[k];
	}
	s->array[0] = pos;

	return s;
}

vec_p
entry_one2d_sloppy(int pos,vec_p n) {
	int  k, size, test;
	vec_p  s;

	size = vec_size( n );
	ASSERT(pos >= 0);

	/*Berechnung*/
	s = vec_new( n->dim );

	for(k=s->dim-1;k>=1;k--) {
	    if (n->array[k] == 0) {
	        s->array[k] = 0;
	    } else {
	        s->array[k] = pos % n->array[k];
		    pos = (pos - s->array[k]) / n->array[k];
	    }
	}
	s->array[0] = pos;

	return s;
}

vec_p
vec_add(vec_p s, vec_p n) {
	int  k;
	vec_p  back;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == n->dim );

	back = vec_new( s->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] = n->array[k] + s->array[k];
	}
	return back;
}

void
vec_add2(vec_p s, vec_p n) {
	int  k;
	
	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == n->dim );

	for(k=0; k<s->dim;k++) {
		s->array[k] += n->array[k];
	}
}

vec_real_p
vec_real_substract(vec_real_p s, vec_real_p n) {
	int  k;
	vec_real_p  back;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == n->dim );

	back = vec_real_new( s->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] = s->array[k] - n->array[k];
	}
	return back;
}


vec_p
vec_multi(int a, vec_p n) {
	int  k;
	vec_p  back;

	back = vec_new( n->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] =  a * n->array[k];
	}
	return back;
}

vec_p
vec_div(int a, vec_p n) {
	int  k;
	vec_p  back;

	back = vec_new( n->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] =  ( n->array[k] ) / a;
	}
	return back;
}


vec_real_p
vec_real_multi(fepc_real_t factor, vec_real_p vector) {
    int n;
    
    vec_real_p back = vec_real_new(vector->dim);
    
    for (n = 0; n < vector->dim; n++) {
        back->array[n] = vector->array[n]*factor;
    }
    return back;
}


void
vec_real_multi2(fepc_real_t factor, vec_real_p vector){
    int n;
    
    for (n = 0; n < vector->dim; n++) {
        vector->array[n] *= factor;
    }
}

vec_p
vec_copy(vec_p n) {
	int  k;
	vec_p  back;

	back = vec_new( n->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] = n->array[k];
	}
	return back;
}

int
vec_skalar_prod(vec_p r, vec_p s) {
	int  k, sum;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == r->dim );

	sum = 0;
	for(k=0;k<s->dim;k++) {
		sum = sum + r->array[k] * s->array[k];
	}
	return sum;
}

vec_p
vec_min(vec_p r, vec_p s) {
	int  k;
	vec_p  back;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == r->dim );

	back = vec_new( s->dim );
	for(k=0; k<back->dim;k++) {
		if ( s->array[k] < r->array[k] ) {
			back->array[k] = s->array[k];
		}
		else {
			back->array[k] = r->array[k];
		}
	}
	return back;
}

vec_p
vec_max(vec_p r, vec_p s) {
	int  k;
	vec_p  back;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == r->dim );

	back = vec_new( s->dim );
	for(k=0; k<back->dim;k++) {
		if ( s->array[k] > r->array[k] ) {
			back->array[k] = s->array[k];
		}
		else {
			back->array[k] = r->array[k];
		}
	}
	return back;
}


int
vec_size(vec_p r) {
	int  k, prod;

	prod = 1;
	for(k=0;k<r->dim;k++) {
		prod = prod * r->array[k];
	}
	return prod;
}






vec_p
vec_r_s_n(vec_p r, vec_p n) {
	int  k, d, dim;
	vec_p  back, constant, prod, modulo;
	int  *b, *c, *p , *m;
	int  q, size, sum, i, temp;

	ASSERT( r->dim == n->dim );
	dim = r->dim;
	for(k=0;k<dim;k++) {
		ASSERT( n->array[k] > 0 );
	}


	/*Berechnungen der Hilfsvektoren*/
	prod = vec_new(dim);
	p = prod->array;
	p[dim-1] = 1;
	for(k=dim-2;k>=0;k--) {
		p[k] = n->array[k+1] * p[k+1];
	}

	constant = vec_new(dim);
	c = constant->array;
	for(k=0;k<dim;k++) {
		c[k] = (n->array[k] - r->array[k]);
	}

	modulo = vec_new(dim);
	m = modulo->array;
	m[dim-1] = 1;
	for(k=dim-2;k>=0;k--) {
		m[k] = m[k+1] * (c[k+1]);
	}

	/*Berechnung des Ausgabearrays*/
	size = vec_size(constant);
	for(k=0;k<dim;k++) {
		ASSERT( constant->array[k] > 0 );
	}
	back = vec_new(size);
	b = back->array;

	q = entry_d2one( r , n);
	for(k=0;k<size;k++) {
		sum = q;
		for(d=dim-1;d>=0;d--) {
			i = k / m[d];
			temp = ( i % c[d] ) * p[d];
			sum = sum + temp;
		}
		b[k] = sum;
	}

	vec_del( constant );
	vec_del( prod );
	vec_del( modulo );
	return back;
}






vec_p
vec_op(int a, vec_p s, int b, vec_p n) {
	int  k;
	vec_p  back;

	/*Testen auf Konsistenz (siehe Dokumentation)*/
	ASSERT( s->dim == n->dim );

	back = vec_new( s->dim );
	for(k=0; k<back->dim;k++) {
		back->array[k] = a*s->array[k] + b*n->array[k];
	}
	return back;
}



vec_p*
vec_zerlegung(vec_p n) {
	int  k, dim;
	vec_p  *back;
	vec_p  s , r;

	dim = n->dim;
	r = vec_new(dim);
	s = vec_new(dim);

	for(k=0;k<dim;k++) {
		if( n->array[k] >= 0 ) {
			s->array[k] = n->array[k] / 2;
			r->array[k] = n->array[k]%2;
		}
		else {
			s->array[k] = (n->array[k]-1) / 2;
			r->array[k] = -(n->array[k]%2);
		}
	}
	back = (vec_p*) malloc(sizeof(vec_p) * 2);
	ASSERT(back != NULL);

	back[0] = s;
	back[1] = r;

	return back;
}



vec_p
vec_0_s_r(vec_p r, vec_p n) {
	int  k, d, dim;
	vec_p  back, constant, prod, modulo;
	int  *b, *c, *p , *m;
	int  size, sum, i, temp;

	ASSERT( r->dim == n->dim );
	dim = r->dim;
	for(k=0;k<dim;k++) {
		ASSERT( n->array[k] > 0 );
	}

	/*Berechnungen der Hilfsvektoren*/
	prod = vec_new(dim);
	p = prod->array;
	p[dim-1] = 1;
	for(k=dim-2;k>=0;k--) {
		p[k] = n->array[k+1] * p[k+1];
	}

	constant = vec_new(dim);
	c = constant->array;
	for(k=0;k<dim;k++) {
		c[k] = r->array[k] + 1;
	}

	modulo = vec_new(dim);
	m = modulo->array;
	m[dim-1] = 1;
	for(k=dim-2;k>=0;k--) {
		m[k] = m[k+1] * (c[k+1]);
	}

	/*Berechnung des Ausgabearrays*/
	size = vec_size(constant);
	for(k=0;k<dim;k++) {
		ASSERT( constant->array[k] > 0 );
	}
	back = vec_new(size);
	b = back->array;


	for(k=0;k<size;k++) {
		sum = 0;
		for(d=dim-1;d>=0;d--) {
			i = k / m[d];
			temp = ( i % c[d] ) * p[d];
			sum = sum + temp;
		}
		b[k] = sum;
	}

	vec_del( constant );
	vec_del( prod );
	vec_del( modulo );
	return back;
}



void 
print_vec(vec_p vector) {
    int n;
    
    for (n = 0; n < vector->dim; n++) {
        printf("\t%i\n", vector->array[n]);
    }
}



void 
print_vec_real(vec_real_p vector) {
    int n;
    
    for (n = 0; n < vector->dim; n++) {
        printf("\t%f\n", vector->array[n]);
    }
}

fepc_real_t
vec_real_norm(vec_real_p vector) {
    int n;
    
    fepc_real_t result = 0.;
    
    for (n = 0; n < vector->dim; n++) {
        result += vector->array[n]*vector->array[n];
    }
    return result;
}
