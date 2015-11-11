/* Copyright (C) 2009-2011, Stanford University
 * This file is part of S4
 * Written by Victor Liu (vkl@stanford.edu)
 *
 * S4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * S4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "sort.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int isqrt(int u){
	if(u <= 0){ return 0; }
	int a = 2;
	int b = u/a;
	while(((a-b) > 1) || ((b-a) > 1)){
		a = (a+b)/2;
		b = u/a;
	}
	return (a<b)?a:b;
}

static double Gcmp_d(const void *i, const void *j, void *arg){
	// Comparing lengths of two vectors:
	//  Vector u is defined by
	//   u[0] = G[2*i+0] * Lk[0] + G[2*i+1] * Lk[2]
	//   u[1] = G[2*i+0] * Lk[1] + G[2*i+1] * Lk[3]
	//  Similarly for v and index j.
	//  The length is then:
	//   |u|^2 = w^T [ G[2*i+0]*G[2*i+0] ]
	//               [ G[2*i+0]*G[2*i+1] ]
	//               [ G[2*i+1]*G[2*i+1] ]
	//  where
	//   w = [    Lk[0]*Lk[0]+Lk[1]*Lk[1]  ]
	//       [ 2*(Lk[0]*Lk[2]+Lk[1]*Lk[3]) ]
	//       [    Lk[2]*Lk[2]+Lk[3]*Lk[3]  ]
	const int *Gi = (int*)i;
	const int *Gj = (int*)j;

	const double *w = (const double*)arg;
	int dv[3] = {
		Gi[0]*Gi[0] - Gj[0]*Gj[0],
		Gi[0]*Gi[1] - Gj[0]*Gj[1],
		Gi[1]*Gi[1] - Gj[1]*Gj[1]
	};
	return dv[0]*w[0] + dv[1]*w[1] + dv[2]*w[2];
}

static int Gcmp(const void *i, const void *j, void *arg){
	double d = Gcmp_d(i, j, arg);
	if(d > 0){ return 1; }
	if(d < 0){ return -1; }
	return 0;
}
static int G_same(const int Gi[2], const int Gj[2], const double Lk[4]){
	double Lkprod[3] = {
		Lk[0]*Lk[0]+Lk[1]*Lk[1],
		2.*(Lk[0]*Lk[2]+Lk[1]*Lk[3]),
		Lk[2]*Lk[2]+Lk[3]*Lk[3]
	};
	const double ilen = hypot(
		Gi[0]*Lk[0]+Gi[1]*Lk[2],
		Gi[0]*Lk[1]+Gi[1]*Lk[3]
	);
	const double jlen = hypot(
		Gj[0]*Lk[0]+Gj[1]*Lk[2],
		Gj[0]*Lk[1]+Gj[1]*Lk[3]
	);
	const double maxlen = (ilen > jlen) ? ilen : jlen;
	double d = fabs(Gcmp_d(&Gi[0], &Gj[0], &Lkprod[0]));
	return d < 2*DBL_EPSILON*maxlen;
}

static int Gsel_parallelogramic(unsigned int *NG, const double Lk[4], int *G){
	int M, NGroot = isqrt(*NG);
	int i, j;
	double Lkprod[3] = {
		Lk[0]*Lk[0]+Lk[1]*Lk[1],
		2.*(Lk[0]*Lk[2]+Lk[1]*Lk[3]),
		Lk[2]*Lk[2]+Lk[3]*Lk[3]
	};

	if(0 == (NGroot % 2)){ NGroot--; }
	M = NGroot/2;

	for(j = 0; j < NGroot; ++j){
		for(i = 0; i < NGroot; ++i){
			G[2*(i+j*NGroot)+0] = i-M;
			G[2*(i+j*NGroot)+1] = j-M;
		}
	}
	*NG = NGroot*NGroot;
	sort(G, *NG, 2*sizeof(int), &Gcmp, &Lkprod[0]);
	return 0;
}

static int Gsel_circular(unsigned int *NG, const double Lk[4], int *G){
	// NG * |u x v| is approximately the area in k-space we will need cover
	// with a circular disc. (u and v are the 2 shortest lattice vectors)
	// From the area, we can find the radius (and round it up). Then, we
	// can find the minimum extends in each of the two lattice directions.

	const double u = hypot(Lk[0],Lk[1]);
	const double v = hypot(Lk[2],Lk[3]);
	const double u2 = Lk[0]*Lk[0] + Lk[1]*Lk[1];
	const double v2 = Lk[2]*Lk[2] + Lk[3]*Lk[3];
	const double uv = Lk[0]*Lk[2] + Lk[1]*Lk[3];
	double Lkprod[3] = { u2, 2*uv, v2 };
	const double uxv = fabs(Lk[0]*Lk[3] - Lk[1]*Lk[2]);

	const double circ_area = (double)(*NG) * uxv;
	const double circ_radius = sqrt(circ_area/M_PI) + u+v;

	const int u_extent = 1+(int)(circ_radius/(u*sqrt(1.-uv*uv/(u2*v2))));
	const int v_extent = 1+(int)(circ_radius/(v*sqrt(1.-uv*uv/(u2*v2))));
	const int uext21 = 2*u_extent+1;
	const int vext21 = 2*v_extent+1;
	int *Gtemp = (int*)malloc(sizeof(int)*2*uext21*vext21);
	int i, j;

	for(i = 0; i < uext21; ++i){
		for(j = 0; j < vext21; ++j){
			Gtemp[2*(i+j*uext21)+0] = i-u_extent;
			Gtemp[2*(i+j*uext21)+1] = j-v_extent;
		}
	}
	sort(Gtemp, uext21*vext21, 2*sizeof(int), &Gcmp, &Lkprod[0]);
	j = uext21*vext21;
	if((int)*NG < j){ j = *NG; }
	for(i = j; i > 0; --i){
		if(!G_same(&Gtemp[2*(i-1)], &Gtemp[2*(i-0)], Lk)){
			break;
		}
	}
	*NG = i;
	for(i = 0; i < (int)*NG; ++i){
		G[2*i+0] = Gtemp[2*i+0];
		G[2*i+1] = Gtemp[2*i+1];
	}
	free(Gtemp);

	return 0;
}

int G_select(const int method, unsigned int *NG, const double Lk[4], int *G){
	if(NULL == NG){ return -2; }
	if(*NG < 1){ return -2; }
	if(NULL == Lk){ return -3; }
	if(NULL == G){ return -4; }

	if(1 == *NG){
		G[0] = 0;
		G[1] = 0;
		return 0;
	}

	if(0 == method){
		return Gsel_circular(NG, Lk, G);
	}else if(1 == method){
		return Gsel_parallelogramic(NG, Lk, G);
	}
	return -1;
}
