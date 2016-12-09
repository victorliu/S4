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

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include "pattern.h"
#include <stdio.h>
#include <float.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_2_PI
#define M_2_PI 0.63661977236758134308
#endif

double Jinc(double x){
	static const double a1[] = {
		 0.1171875,
		-0.1441955566406250,
		 0.6765925884246826,
		-6.883914268109947,
		 1.215978918765359e2,
		-3.302272294480852e3,
		 1.276412726461746e5,
		-6.656367718817688e6,
		 4.502786003050393e8,
		-3.833857520742790e10,
		 4.011838599133198e12,
		-5.060568503314727e14,
		 7.572616461117958e16,
		-1.326257285320556e19};
	static const double b1[] = {
		-0.1025390625,
		 0.2775764465332031,
		-1.993531733751297,
		 2.724882731126854e1,
		-6.038440767050702e2,
		 1.971837591223663e4,
		-8.902978767070678e5,
		 5.310411010968522e7,
		-4.043620325107754e9,
		 3.827011346598605e11,
		-4.406481417852278e13,
		 6.065091351222699e15,
		-9.833883876590679e17,
		 1.855045211579828e20};
	x *= 2*M_PI;
	/* return 2*J_1(x)/x; */
	if(x < 1e-9){ /* J_1(x) ~ x(1-x^2/8)/2 */
		return 1.-x*x/8.;
	}

	if(x <= 0.0){
		return 1.0;
	}
	if(x <= 12.0){
		double x2 = x*x;
		double j1 = 1.0;
		double r = 1.0;
		int k;
		for(k=1;k<=30;k++){
			r *= -0.25*x2/(k*(k+1));
			j1 += r;
			if (fabs(r) < fabs(j1)*1e-15) break;
		}
		return j1;
	}else{
		double cu = sqrt(M_2_PI/x);
		double t2 = x-0.75*M_PI;
		double p1 = 1.0;
		double q1 = 0.375/x;
		int k;
		int kz;
		if (x >= 50.0) kz = 8;          /* Can be changed to 10 */
		else if (x >= 35.0) kz = 10;    /*  "       "        12 */
		else kz = 12;                   /*  "       "        14 */
		for(k=0;k<kz;k++){
			p1 += a1[k]*pow(x,-2*k-2);
			q1 += b1[k]*pow(x,-2*k-3);
		}
		return 2.0*cu*(p1*cos(t2)-q1*sin(t2))/x;
	}
}
static double my_j0(double x){ /* spherical bessel of first kind, order 0 */
	if(fabs(x) < 1e-9){
		x = 1.-x*x/6.;
	}else{
		x = sin(x)/x;
	}
	return x;
}
double Sinc(double x){
	return my_j0(M_PI*x);
}
static double shape_area(const shape *s){
	switch(s->type){
	case CIRCLE:
		return M_PI * s->vtab.circle.radius*s->vtab.circle.radius;
	case ELLIPSE:
		return M_PI * s->vtab.ellipse.halfwidth[0]*s->vtab.ellipse.halfwidth[1];
	case RECTANGLE:
		return 4.0 * s->vtab.ellipse.halfwidth[0]*s->vtab.ellipse.halfwidth[1];
	case POLYGON:
		{
			int i;
			double a = 0;
			for(i = 0; i < s->vtab.polygon.n_vertices; ++i){
				int j = (i+1)%s->vtab.polygon.n_vertices;
				a += (s->vtab.polygon.vertex[2*i+0]*s->vtab.polygon.vertex[2*j+1] - s->vtab.polygon.vertex[2*i+1]*s->vtab.polygon.vertex[2*j+0]);
			}
			return 0.5*a;
		}
	default:
		return 0;
	}
}

static void shape_get_interior_point(const shape *s, double x[2]){
	switch(s->type){
	case CIRCLE:
	case ELLIPSE:
	case RECTANGLE:
		x[0] = s->center[0]; x[1] = s->center[1];
		return;
	case POLYGON:
		x[0] = s->center[0] + (1./3.) * (s->vtab.polygon.vertex[0] + s->vtab.polygon.vertex[2] + s->vtab.polygon.vertex[4]);
		x[1] = s->center[1] + (1./3.) * (s->vtab.polygon.vertex[1] + s->vtab.polygon.vertex[3] + s->vtab.polygon.vertex[5]);
		return;
	default:
		x[0] = s->center[0]; x[1] = s->center[1];
		return;
	}
}

static int shape_contains_point(const shape *s, const double x_[2]){
	double x[2];
	const double ca = cos(s->angle);
	const double sa = sin(s->angle);
	x[0] = (x_[0] - s->center[0]) * ca + (x_[1] - s->center[1]) * sa;
	x[1] = (x_[0] - s->center[1]) *-sa + (x_[1] - s->center[1]) * ca;
	switch(s->type){
	case CIRCLE:
		return x[0]*x[0] + x[1]*x[1] <= s->vtab.circle.radius*s->vtab.circle.radius;
	case ELLIPSE:
		{
			double r = s->vtab.ellipse.halfwidth[0]*s->vtab.ellipse.halfwidth[0] - s->vtab.ellipse.halfwidth[1]*s->vtab.ellipse.halfwidth[1];
			double f, L, d; /* f = focus length, L is twice the semi-major axis length */
			if(r >= 0){
				f = sqrt(r);
				L = 2 * s->vtab.ellipse.halfwidth[0];
				d = sqrt((x[0]-f)*(x[0]-f) + x[1]*x[1]) + sqrt((x[0]+f)*(x[0]+f) + x[1]*x[1]);
			}else{
				f = sqrt(-r);
				L = 2 * s->vtab.ellipse.halfwidth[1];
				d = sqrt((x[1]-f)*(x[1]-f) + x[0]*x[0]) + sqrt((x[1]+f)*(x[1]+f) + x[0]*x[0]);
			}
			return d < L;
		}
	case RECTANGLE:
		return (fabs(x[0]) < s->vtab.rectangle.halfwidth[0]) && (fabs(x[1]) < s->vtab.rectangle.halfwidth[1]);
	case POLYGON:
		{ /* From http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html */
			int i, j;
			int c = 0;
			for(i = 0, j = s->vtab.polygon.n_vertices-1; i < s->vtab.polygon.n_vertices; j = i++){
				double vix = s->vtab.polygon.vertex[2*i+0];
				double viy = s->vtab.polygon.vertex[2*i+1];
				double vjx = s->vtab.polygon.vertex[2*j+0];
				double vjy = s->vtab.polygon.vertex[2*j+1];
				if ( ((viy>x[1]) != (vjy>x[1]))
				&& (x[0] < (vjx-vix) * (x[1]-viy) / (vjy-viy) + vix) ){ c = !c; }
			}
			return c;
		}
	default:
		return 0;
	}
}

int shape_get_normal(const shape *s, const double x[2], double n[2]){
	double r[2];
	if(NULL == s){ return -1; }
	if(NULL == x){ return -2; }
	if(NULL == n){ return -3; }
	r[0] = x[0] - s->center[0];
	r[1] = x[1] - s->center[1];

	switch(s->type){
	case CIRCLE:
		n[0] = r[0];
		n[1] = r[1];
		break;
	case ELLIPSE:
		/* The ellipse is defined by the equation
		 *   Norm_2(B.{x,y})^2 == 1
		 * where
		 *   B = diag{1/halfwidth[0],1/halfwidth[1]}.[ cos(angle) sin(angle) ]
		 *                                           [-sin(angle) cos(angle) ]
		 * let ca = cos(angle), sa = sin(angle), ilx2 = 1/lx^2, ily2 = 1/ly^2
		 * Let A = B^T B = [ ca*ca*ilx2 + sa*sa*ily2  ca*sa*(ilx2-ily2) ]
		 *                 [ same as other      ca*ca*ily2 + sa*sa*ilx2 ]
		 * The gradient is then A.{x,y}
		 */
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double ilx2 = 1./(s->vtab.ellipse.halfwidth[0]*s->vtab.ellipse.halfwidth[0]);
			const double ily2 = 1./(s->vtab.ellipse.halfwidth[1]*s->vtab.ellipse.halfwidth[1]);
			const double a = ca*ca*ilx2 + sa*sa*ily2;
			const double b = ca*ca*ily2 + sa*sa*ilx2;
			const double c = ca*sa*(ilx2-ily2);
			n[0] = a*r[0] + c*r[1];
			n[1] = c*r[0] + b*r[1];
		}
		break;
	case RECTANGLE:
		/* The rectangle is defined by the equation
		 *   max(abs(rx)-halfwidth[0], abs(ry)-halfwidth[1]) == 0
		 * where
		 *   rx = proj onto { ca,sa} of r
		 *   ry = proj onto {-sa,ca} of r
		 */
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double rx = ca*r[0] + sa*r[1];
			const double ry =-sa*r[0] + ca*r[1];
			if((fabs(rx)-s->vtab.rectangle.halfwidth[0]) > (fabs(ry)-s->vtab.rectangle.halfwidth[1])){
				double sgn = (rx > 0. ? 1. : -.1);
				n[0] = sgn*ca; n[1] = sgn*sa;
			}else{
				double sgn = (ry > 0. ? 1. : -.1);
				n[0] = sgn*(-sa); n[1] = sgn*ca;
			}
		}
		break;
	case POLYGON:
		{
			int i, j;
			double maxdist = -1;
			for(j = 0, i = s->vtab.polygon.n_vertices-1; j < s->vtab.polygon.n_vertices; i = j++){
				/* compute distance from r to segment */
				double v[2], pr[2];
				v[0] = s->vtab.polygon.vertex[2*j+0] - s->vtab.polygon.vertex[2*i+0];
				v[1] = s->vtab.polygon.vertex[2*j+1] - s->vtab.polygon.vertex[2*i+1];

				pr[0] = r[0] - s->vtab.polygon.vertex[2*i+0];
				pr[1] = r[1] - s->vtab.polygon.vertex[2*i+1];
				{
					const double v2 = v[0]*v[0] + v[1]*v[1];
					const double prj = (pr[0]*v[0] + pr[1]*v[1])/v2;
					const double voff[2] = {pr[0] - prj*v[0], pr[1] - prj*v[1]};
					const double dist = hypot(voff[0], voff[1]);
					if(dist > maxdist){
						maxdist = dist;
						n[0] = v[1];
						n[1] = -v[0];
					}
				}
			}
		}
		break;
	default:
		n[1] = n[0] = 0;
		break;
	}
	{ /* normalize */
		double ln = hypot(n[0], n[1]);
		if(ln > 0){
			n[0] /= ln;
			n[1] /= ln;
		}
	}
	return 0;
}

#include "intersection.h"

/* Returns intersection area between shape s and the rectangle with
 * bottom left corner p0 and dimensions dp.
 */
double shape_get_intersection_area_quad(const shape *s, const double *p0, const double *duv){
	/* We will use the strategy of decomposing polygonal shapes into triangles.
	 * For the ellipse, we will apply a coordinate transformation to make it a circle.
	 */

	double p0p[2] = {p0[0] - s->center[0], p0[1] - s->center[1]};
	double upv[2] = {duv[0]+duv[2], duv[1]+duv[3]};
	double umv[2] = {duv[0]-duv[2], duv[1]-duv[3]};
	double lupv = hypot(upv[0], upv[1]);
	double lumv = hypot(umv[0], umv[1]);
	switch(s->type){
	case CIRCLE:
		{
			double a;
			double tri_org[2], tri_u[2], tri_v[2];
			if(lupv < lumv){ /* split by u+v */
				tri_org[0] = p0p[0];
				tri_org[1] = p0p[1];
				tri_u[0] = duv[0]; tri_u[1] = duv[1];
				tri_v[0] = upv[0]; tri_v[1] = upv[1];
				a = intersection_area_circle_triangle(s->vtab.circle.radius, tri_org, tri_u, tri_v);
				tri_u[0] = upv[0]; tri_u[1] = upv[1];
				tri_v[0] = duv[2]; tri_v[1] = duv[3];
			}else{
				tri_org[0] = duv[2] + p0p[0];
				tri_org[1] = duv[3] + p0p[1];
				tri_u[0] = -duv[2]; tri_u[1] = -duv[3];
				tri_v[0] = umv[0]; tri_v[1] = umv[1];
				a = intersection_area_circle_triangle(s->vtab.circle.radius, tri_org, tri_u, tri_v);
				tri_u[0] = umv[0]; tri_u[1] = umv[1];
				tri_v[0] = duv[0]; tri_v[1] = duv[1];
			}
			return a + intersection_area_circle_triangle(s->vtab.circle.radius, tri_org, tri_u, tri_v);
		}
	case ELLIPSE:
		{
			double a;
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double ratio = (s->vtab.ellipse.halfwidth[0] / s->vtab.ellipse.halfwidth[1]);
			double tri_org[2], tri_u[2], tri_v[2];
			if(lupv < lumv){ /* split by u+v */
				tri_org[0] = (p0p[0]* ca + p0p[1]*sa);
				tri_org[1] = (p0p[0]*-sa + p0p[1]*ca);
				tri_u[0] = duv[0]* ca + duv[1]*sa;
				tri_u[1] = duv[0]*-sa + duv[1]*ca;
				tri_v[0] = upv[0]* ca + upv[1]*sa;
				tri_v[1] = upv[0]*-sa + upv[1]*ca;
				tri_org[1] *= ratio; tri_u[1] *= ratio; tri_v[1] *= ratio;
				a = intersection_area_circle_triangle(s->vtab.ellipse.halfwidth[0], tri_org, tri_u, tri_v);
				tri_u[0] = upv[0]* ca + upv[1]*sa;
				tri_u[1] = upv[0]*-sa + upv[1]*ca;
				tri_v[0] = duv[2]* ca + duv[3]*sa;
				tri_v[1] = duv[2]*-sa + duv[3]*ca;
			}else{
				tri_org[0] = ((p0p[0]+duv[2])* ca + (p0p[1]+duv[3])*sa);
				tri_org[1] = ((p0p[0]+duv[2])*-sa + (p0p[1]+duv[3])*ca);
				tri_u[0] = -(duv[2]* ca + duv[3]*sa);
				tri_u[1] = -(duv[2]*-sa + duv[3]*ca);
				tri_v[0] = umv[0]* ca + umv[1]*sa;
				tri_v[1] = umv[0]*-sa + umv[1]*ca;
				tri_org[1] *= ratio; tri_u[1] *= ratio; tri_v[1] *= ratio;
				a = intersection_area_circle_triangle(s->vtab.ellipse.halfwidth[0], tri_org, tri_u, tri_v);
				tri_u[0] = umv[0]* ca + umv[1]*sa;
				tri_u[1] = umv[0]*-sa + umv[1]*ca;
				tri_v[0] = duv[0]* ca + duv[1]*sa;
				tri_v[1] = duv[0]*-sa + duv[1]*ca;
			}
			tri_org[1] *= ratio; tri_u[1] *= ratio; tri_v[1] *= ratio;
			return a + intersection_area_circle_triangle(s->vtab.ellipse.halfwidth[0], tri_org, tri_u, tri_v);
		}
	case RECTANGLE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double P[8], Q[8], Pi[16], u[2], v[2];
			int ret, nPi = 8;

			P[2*0+0] = p0p[0];
			P[2*0+1] = p0p[1];
			P[2*1+0] = P[2*0+0] + duv[0];
			P[2*1+1] = P[2*0+1] + duv[1];
			P[2*2+0] = P[2*1+0] + duv[2];
			P[2*2+1] = P[2*1+1] + duv[3];
			P[2*3+0] = P[2*0+0] + duv[2];
			P[2*3+1] = P[2*0+1] + duv[3];
			u[0] = s->vtab.rectangle.halfwidth[0] * ca;
			u[1] = s->vtab.rectangle.halfwidth[0] * sa;
			v[0] = s->vtab.rectangle.halfwidth[1] *-sa;
			v[1] = s->vtab.rectangle.halfwidth[1] * ca;
			Q[2*0+0] = -u[0]-v[0];
			Q[2*0+1] = -u[1]-v[1];
			Q[2*1+0] = u[0]-v[0];
			Q[2*1+1] = u[1]-v[1];
			Q[2*2+0] = u[0]+v[0];
			Q[2*2+1] = u[1]+v[1];
			Q[2*3+0] = v[0]-u[0];
			Q[2*3+1] = v[1]-u[1];
			/*{int i;for(i = 0;i<4;++i){fprintf(stderr, " {%f,%f},\n", P[2*i+0], P[2*i+1]);}}*/
			/*{int i;for(i = 0;i<4;++i){fprintf(stderr, " {%f,%f},\n", Q[2*i+0], Q[2*i+1]);}}*/
			ret = convex_polygon_intersection(4,P,4,Q,&nPi,Pi);
			if(1 == ret){ /* pixel completely in rectangle */
				return duv[0]*duv[3]-duv[1]*duv[2];
			}
			/*{int i;for(i = 0;i<nPi;++i){fprintf(stderr, " %f %f\n", Pi[2*i+0], Pi[2*i+1]);}}*/
			return polygon_area(nPi,Pi);
		}
	case POLYGON:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double a = 0;
			double P[8], Q[6], Pi[14], u[2], v[2];
			int i, j, nPi;
			int *tri = (int*)malloc(sizeof(int)*3*(s->vtab.polygon.n_vertices-2));

			P[2*0+0] = p0p[0]*ca + p0p[1]*sa;
			P[2*0+1] = p0p[0]*-sa + p0p[1]*ca;
			u[0] = duv[0]*ca + duv[1]*sa;
			u[1] = duv[0]*-sa+ duv[1]*ca;
			v[0] = duv[2]*ca + duv[3]*sa;
			v[1] = duv[2]*-sa+ duv[3]*ca;
			P[2*1+0] = P[2*0+0] + u[0];
			P[2*1+1] = P[2*0+1] + u[1];
			P[2*2+0] = P[2*1+0] + v[0];
			P[2*2+1] = P[2*1+1] + v[1];
			P[2*3+0] = P[2*0+0] + v[0];
			P[2*3+1] = P[2*0+1] + v[1];

			polygon_triangulate(s->vtab.polygon.n_vertices, s->vtab.polygon.vertex, tri);
			/*
			for(i = 0; i < s->vtab.polygon.n_vertices-2; ++i){
				fprintf(stderr, " (%d %d %d)", tri[3*i+0], tri[3*i+1], tri[3*i+2]);
			}fprintf(stderr, "\n");
			*/
			for(i = 0; i < s->vtab.polygon.n_vertices-2; ++i){
				for(j = 0; j < 3; ++j){
					Q[2*j+0] = s->vtab.polygon.vertex[2*tri[3*i+j]+0];
					Q[2*j+1] = s->vtab.polygon.vertex[2*tri[3*i+j]+1];
				}
				nPi = 7;
				convex_polygon_intersection(4,P,3,Q,&nPi,Pi);
				a += polygon_area(nPi,Pi);
			}
			free(tri);
			return a;
		}
	default:
		break;
	}
	return 0;
}

/* Returns intersection area between shape s and the rectangle with
 * bottom left corner p0 and dimensions dp.
 */
double shape_get_intersection_area(const shape *s, const double p0[2], const double dp[2]){
	/* We will use the strategy of decomposing polygonal shapes into triangles.
	 * For the ellipse, we will apply a coordinate transformation to make it a circle.
	 */

	double p0p[2] = {p0[0] - s->center[0], p0[1] - s->center[1]};
	switch(s->type){
	case CIRCLE:
		{
			double a;
			double tri_org[2], tri_u[2], tri_v[2];
			tri_org[0] = p0[0] - s->center[0];
			tri_org[1] = p0[1] - s->center[1];
			tri_u[0] = dp[0]; tri_u[1] = 0;
			tri_v[0] = dp[0]; tri_v[1] = dp[1];
			a = intersection_area_circle_triangle(s->vtab.circle.radius, tri_org, tri_u, tri_v);
			tri_u[0] = 0; tri_u[1] = dp[1];
			return a + intersection_area_circle_triangle(s->vtab.circle.radius, tri_org, tri_v, tri_u);
		}
	case ELLIPSE:
		{
			double a;
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double ratio = (s->vtab.ellipse.halfwidth[0] / s->vtab.ellipse.halfwidth[1]);
			double tri_org[2], tri_u[2], tri_v[2];
			tri_org[0] = p0p[0]* ca + p0p[1]*sa;
			tri_org[1] = (p0p[0]*-sa + p0p[1]*ca) * ratio;
			tri_u[0] = dp[0]*ca; tri_u[1] = dp[0]*-sa;
			tri_v[0] = tri_u[0] + dp[1]*sa; tri_v[1] = tri_u[1] + dp[1]*ca;
			tri_u[1] *= ratio; tri_v[1] *= ratio;
			a = intersection_area_circle_triangle(s->vtab.ellipse.halfwidth[0], tri_org, tri_u, tri_v);
			tri_u[0] = dp[1]*sa; tri_u[1] = dp[1]*ca * ratio;
			return a + intersection_area_circle_triangle(s->vtab.ellipse.halfwidth[0], tri_org, tri_v, tri_u);
		}
	case RECTANGLE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double P[8], Q[8], Pi[16], u[2], v[2];
			int nPi = 8;
			P[2*0+0] = p0[0] - s->center[0];
			P[2*0+1] = p0[1] - s->center[1];
			P[2*1+0] = P[2*0+0] + dp[0];
			P[2*1+1] = P[2*0+1];
			P[2*2+0] = P[2*1+0];
			P[2*2+1] = P[2*1+1] + dp[1];
			P[2*3+0] = P[2*0+0];
			P[2*3+1] = P[2*0+1] + dp[1];
			u[0] = s->vtab.rectangle.halfwidth[0] * ca;
			u[1] = s->vtab.rectangle.halfwidth[0] * sa;
			v[0] = s->vtab.rectangle.halfwidth[1] *-sa;
			v[1] = s->vtab.rectangle.halfwidth[1] * ca;
			Q[2*0+0] = -u[0]-v[0];
			Q[2*0+1] = -u[1]-v[1];
			Q[2*1+0] = u[0]-v[0];
			Q[2*1+1] = u[1]-v[1];
			Q[2*2+0] = u[0]+v[0];
			Q[2*2+1] = u[1]+v[1];
			Q[2*3+0] = v[0]-u[0];
			Q[2*3+1] = v[1]-u[1];
			convex_polygon_intersection(4,P,4,Q,&nPi,Pi);
			return polygon_area(nPi,Pi);
		}
	case POLYGON:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double a = 0;
			double P[8], Q[6], Pi[14], u[2], v[2];
			int i, j, nPi;
			int *tri = (int*)malloc(sizeof(int)*3*(s->vtab.polygon.n_vertices-2));

			P[2*0+0] = p0p[0]*ca + p0p[1]*sa;
			P[2*0+1] = p0p[0]*-sa + p0p[1]*ca;
			u[0] = dp[0]*ca; u[1] = dp[0]*-sa;
			v[0] = dp[1]*sa; v[1] = dp[1]* ca;
			P[2*1+0] = P[2*0+0] + u[0];
			P[2*1+1] = P[2*0+1] + u[1];
			P[2*2+0] = P[2*1+0] + v[0];
			P[2*2+1] = P[2*1+1] + v[1];
			P[2*3+0] = P[2*0+0] + v[0];
			P[2*3+1] = P[2*0+1] + v[1];

			polygon_triangulate(s->vtab.polygon.n_vertices, s->vtab.polygon.vertex, tri);
			/*
			for(i = 0; i < s->vtab.polygon.n_vertices-2; ++i){
				fprintf(stderr, " (%d %d %d)", tri[3*i+0], tri[3*i+1], tri[3*i+2]);
			}fprintf(stderr, "\n");
			*/
			for(i = 0; i < s->vtab.polygon.n_vertices-2; ++i){
				for(j = 0; j < 3; ++j){
					Q[2*j+0] = s->vtab.polygon.vertex[2*tri[3*i+j]+0];
					Q[2*j+1] = s->vtab.polygon.vertex[2*tri[3*i+j]+1];
				}
				nPi = 7;
				convex_polygon_intersection(4,P,3,Q,&nPi,Pi);
				a += polygon_area(nPi,Pi);
			}
			free(tri);
			return a;
		}
	default:
		break;
	}
	return 0;
}

/* segment is from p0 to p0+d
 * returns number of intersections
 */
int shape_get_tangent_cross_segment(const shape *s, const double p0[2], const double d0[2], double *cross){
	int i, c = 0;
	double p0p[2] = {p0[0] - s->center[0], p0[1] - s->center[1]};
	*cross = 0;
	switch(s->type){
	case CIRCLE:
		{
			double isect[4];
			c = intersection_circle_segment(s->vtab.circle.radius, p0p, d0, isect, NULL);
			for(i = 0; i < c; ++i){
				/* Each intersection (x,y) implies a tangential normal vector (-y,x)/r */
				*cross -= (d0[0]*isect[2*i+0] + d0[1]*isect[2*i+1]);
			}
			*cross /= s->vtab.circle.radius;
			break;
		}
	case ELLIPSE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double ratio = (s->vtab.ellipse.halfwidth[0] / s->vtab.ellipse.halfwidth[1]);
			const double iratio = 1./ratio;

			double p1[2], u[2], isect[4], u1;
			p1[0] = p0p[0]* ca + p0p[1]*sa;
			p1[1] = (p0p[0]*-sa + p0p[1]*ca) * ratio;
			u1 = d0[0]*-sa;
			u[0] = d0[0]*ca; u[1] = u1 * ratio;

			c = intersection_circle_segment(s->vtab.ellipse.halfwidth[0], p1, u, isect, NULL);
			for(i = 0; i < c; ++i){
				double L;
				isect[2*c+1] *= iratio;
				L = hypot(isect[2*c+0],isect[2*c+1]);
				if(0 == L){ L = 1; }
				*cross += (u[0]*isect[2*c+0] + u1*isect[2*c+1]) / L;
			}
			break;
		}
	case RECTANGLE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double P[8], u[2], v[2];
			u[0] = s->vtab.rectangle.halfwidth[0] * ca;
			u[1] = s->vtab.rectangle.halfwidth[0] * sa;
			v[0] = s->vtab.rectangle.halfwidth[1] *-sa;
			v[1] = s->vtab.rectangle.halfwidth[1] * ca;
			P[2*0+0] = -u[0]-v[0];
			P[2*0+1] = -u[1]-v[1];
			P[2*1+0] = u[0]-v[0];
			P[2*1+1] = u[1]-v[1];
			P[2*2+0] = u[0]+v[0];
			P[2*2+1] = u[1]+v[1];
			P[2*3+0] = v[0]-u[0];
			P[2*3+1] = v[1]-u[1];

			c = intersection_polygon_segment(4, P, p0p, d0, NULL, cross, NULL, NULL);
			break;
		}
	case POLYGON:
		{
			c = intersection_polygon_segment(
				s->vtab.polygon.n_vertices,
				s->vtab.polygon.vertex,
				p0p, d0, NULL, cross, NULL, NULL);
			break;
		}
	default:
		break;
	}
	return c;
}

/* cross0 is the cross product weighted by the t's
 * cross1 is the cross product weighted by (1-t)'s
 */
int shape_get_tangent_cross_segment_tri(const shape *s, const double p0[2], const double d0[2], double *cross0, double *cross1){
	int i, c = 0;
	double p0p[2] = {p0[0] - s->center[0], p0[1] - s->center[1]};
	*cross0 = 0;
	*cross1 = 0;
	switch(s->type){
	case CIRCLE:
		{
			double isect[4];
			double t[2];
			c = intersection_circle_segment(s->vtab.circle.radius, p0p, d0, isect, t);
			for(i = 0; i < c; ++i){
				/* Each intersection (x,y) implies a tangential normal vector (-y,x)/r */
				double cross = (d0[0]*isect[2*i+0] + d0[1]*isect[2*i+1]);
				*cross0 -= cross*t[i];
				*cross1 -= cross*(1-t[i]);
			}
			*cross0 /= s->vtab.circle.radius;
			*cross1 /= s->vtab.circle.radius;
			break;
		}
	case ELLIPSE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			const double ratio = (s->vtab.ellipse.halfwidth[0] / s->vtab.ellipse.halfwidth[1]);
			const double iratio = 1./ratio;

			double p1[2], u[2], isect[4], u1, t[2];
			p1[0] = p0p[0]* ca + p0p[1]*sa;
			p1[1] = (p0p[0]*-sa + p0p[1]*ca) * ratio;
			u1 = d0[0]*-sa;
			u[0] = d0[0]*ca; u[1] = u1 * ratio;

			c = intersection_circle_segment(s->vtab.ellipse.halfwidth[0], p1, u, isect, t);
			for(i = 0; i < c; ++i){
				double L;
				isect[2*c+1] *= iratio;
				L = hypot(isect[2*c+0],isect[2*c+1]);
				if(0 == L){ L = 1; }
				*cross0 += t[i]*(u[0]*isect[2*c+0] + u1*isect[2*c+1]) / L;
				*cross1 += (1-t[i])*(u[0]*isect[2*c+0] + u1*isect[2*c+1]) / L;
			}
			break;
		}
	case RECTANGLE:
		{
			const double ca = cos(s->angle);
			const double sa = sin(s->angle);
			double dummy;
			double P[8], u[2], v[2];
			u[0] = s->vtab.rectangle.halfwidth[0] * ca;
			u[1] = s->vtab.rectangle.halfwidth[0] * sa;
			v[0] = s->vtab.rectangle.halfwidth[1] *-sa;
			v[1] = s->vtab.rectangle.halfwidth[1] * ca;
			P[2*0+0] = -u[0]-v[0];
			P[2*0+1] = -u[1]-v[1];
			P[2*1+0] = u[0]-v[0];
			P[2*1+1] = u[1]-v[1];
			P[2*2+0] = u[0]+v[0];
			P[2*2+1] = u[1]+v[1];
			P[2*3+0] = v[0]-u[0];
			P[2*3+1] = v[1]-u[1];

			c = intersection_polygon_segment(4, P, p0p, d0, NULL, &dummy, cross0, cross1);
			break;
		}
	case POLYGON:
		{
			double dummy;
			c = intersection_polygon_segment(
				s->vtab.polygon.n_vertices,
				s->vtab.polygon.vertex,
				p0p, d0, NULL, &dummy, cross0, cross1);
			break;
		}
	default:
		break;
	}
	return c;
}

int shape_valid(const shape *s){
	switch(s->type){
	case CIRCLE:
		if(s->vtab.circle.radius < 0){ return 0; }
		break;
	case ELLIPSE:
		if(s->vtab.ellipse.halfwidth[0] < 0){ return 0; }
		if(s->vtab.ellipse.halfwidth[1] < 0){ return 0; }
		break;
	case RECTANGLE:
		if(s->vtab.rectangle.halfwidth[0] < 0){ return 0; }
		if(s->vtab.rectangle.halfwidth[1] < 0){ return 0; }
		break;
	case POLYGON:
		if(s->vtab.polygon.n_vertices < 3){ return 0; }
		if(NULL == s->vtab.polygon.vertex){ return 0; }
		if(0 != polygon_check(s->vtab.polygon.n_vertices, s->vtab.polygon.vertex)){ return 0; }
		break;
	default:
		return 0;
	}
	return 1;
}

int shapes_intersect(const shape *s0, const shape *s1){
	switch(s0->type){
	case CIRCLE:
		switch(s1->type){
		case CIRCLE:
			{ /* circle-circle */
				/*return hypot(s1->center[0]-s0->center[0], s1->center[1]-s0->center[1]) <= s0->vtab.circle.radius + s1->vtab.circle.radius;*/
			}
			break;
		case ELLIPSE:
			break;
		case RECTANGLE:
			break;
		case POLYGON:
			break;
		default:
			return 0;
		}
		break;
	case ELLIPSE:
		break;
	case RECTANGLE:
		break;
	case POLYGON:
		break;
	default:
		return 0;
	}
	return 0;
}

int pattern_get_containment_tree(
	int nshapes,
	shape *shapes,
	int *parent
){
	double *area;
	int i, j, k;
	if(0 == nshapes){ return 0; }
	if(nshapes < 0){ return -1; }
	if(NULL == shapes){ return -2; }
	if(NULL == parent){ return -3; }

	/* We will assume for now that the non-self-intersection criterion is met. */
	if(nshapes > 1){
		for(i = 0; i < nshapes; ++i){
			if(!shape_valid(&shapes[i])){ return i+1; }
			for(j = i+1; i < nshapes; ++i){
				if(shapes_intersect(&shapes[i], &shapes[j])){
					return nshapes+i+1;
				}
			}
		}
	}

	area = (double*)malloc(sizeof(double)*nshapes);
	for(i = 0; i < nshapes; ++i){
		parent[i] = -1;
		area[i] = shape_area(&shapes[i]);
	}

	/* Sort by area */
	for(k = 1; k < nshapes; ++k){
		for(j = 0; j < k; ++j){
			if(area[j] < area[k]){
				shape ts = shapes[j];
				double t = area[j]; area[j] = area[k]; area[k] = t;
				shapes[j] = shapes[k];
				shapes[k] = ts;
			}
		}
	}

	/* The following is O(n^2), we might be able to do better with sorting
	 * by area but n is usually small.
	 */
	for(i = 1; i < nshapes; ++i){
		double p[2];
		shape_get_interior_point(&shapes[i], p);

		for(j = i-1; j >= 0; --j){
			if(shape_contains_point(&shapes[j], p)){
				parent[i] = j;
				break;
			}
		}
	}
	free(area);

	return 0;
}
int Pattern_GetContainmentTree(
	Pattern *p
){
	return pattern_get_containment_tree(p->nshapes, p->shapes, p->parent);
}

int pattern_get_shape(
	int nshapes, /* number of shapes in the array 'shapes' */
	const shape *shapes, /* array of shapes, ordered by decreasing area */
	const int *parent, /* from Pattern_GetContainmentTree */
	const double x[2], /* location in the plane to get the tag */
	int *shape_index, /* returned shape index */
	double n[2]
){
	int i, found = 0;
	if(0 == nshapes){ return 1; }
	if(nshapes < 0){ return -1; }
	if(NULL == shapes){ return -2; }
	if(NULL == parent){ return -3; }
	if(NULL == x){ return -4; }
	if(NULL == shape_index){ return -5; }

	for(i = nshapes-1; i >= 0; --i){
		if(shape_contains_point(&shapes[i], x)){
			*shape_index = i;
			found = 1;
			break;
		}
	}
	if(!found){ return 1; }
	if(NULL != n){
		shape_get_normal(&shapes[i], x, n);
	}
	return 0;
}
int Pattern_GetShape(
	const Pattern *p,
	const double x[2],
	int *shape_index,
	double n[2]
){
	return pattern_get_shape(p->nshapes, p->shapes, p->parent, x, shape_index, n);
}

/* returns 0 on success
 * returns -n if n-th argument is invalid
 */
int pattern_get_fourier_transform(
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double *value,
	const double k_[2], /* remember, this is k/2pi */
	int ndim,
	double unit_cell_size,
	double f[2]
){
	int i;
	const int DC = (0 == k_[0] && 0 == k_[1]) ? 1 : 0;
	double inv_size;
	double k[2];
	
	if(nshapes < 0){ return -1; }
	if(nshapes > 0 && NULL == shapes){ return -2; }
	if(NULL == parent){ return -3; }
	if(NULL == value){ return -4; }
	if(NULL == k_){ return -5; }
	if(ndim < 1 || ndim > 2){ return -6; }
	if(unit_cell_size <= 0){ return -7; }
	if(NULL == f){ return -8; }
	if(1 == ndim){
		if(k_[1] != 0){ return -6-5; }
		for(i = 0; i < nshapes; ++i){
			const shape *s = &shapes[i];
			if(RECTANGLE != s->type){
				return -6-2;
			}
		}
	}

	
	if(DC){
		f[0] = value[0]; f[1] = value[1];
	}else{
		f[0] = 0; f[1] = 0;
	}
	if(0 == nshapes){
		return 0;
	}

	inv_size = 1./unit_cell_size;

	
	for(i = 0; i < nshapes; ++i){
		const shape *s = &shapes[i];

		double dval[2] = {value[2*(i+1)+0]-value[2*(parent[i]+1)+0], value[2*(i+1)+1]-value[2*(parent[i]+1)+1]};

		double phase_angle = -2*M_PI*(k_[0]*s->center[0] + k_[1]*s->center[1]); /* phase = exp(i*phase_angle); */
		double z[2] = {0,0};
		double area;

		const double ca = cos(s->angle);
		const double sa = sin(s->angle);
		k[0] = k_[0] * ca + k_[1] * sa;
		k[1] = k_[0] *-sa + k_[1] * ca;

		/* Each shape should set z to be the Fourier component, but without dval, and without area. */
		switch(s->type){
		case CIRCLE:
			area = M_PI*s->vtab.circle.radius*s->vtab.circle.radius;
			z[0] =Jinc(s->vtab.circle.radius*hypot(k[0],k[1]));
			break;
		case ELLIPSE:
			area = M_PI*s->vtab.ellipse.halfwidth[0]*s->vtab.ellipse.halfwidth[1];
			if(s->vtab.ellipse.halfwidth[0] >= s->vtab.ellipse.halfwidth[1]){
				double r = s->vtab.ellipse.halfwidth[1] /  s->vtab.ellipse.halfwidth[0] * k[1];
				z[0] = Jinc(s->vtab.ellipse.halfwidth[0]*hypot(k[0],r));
			}else{
				double r = s->vtab.ellipse.halfwidth[0] /  s->vtab.ellipse.halfwidth[1] * k[0];
				z[0] =Jinc(s->vtab.ellipse.halfwidth[1]*hypot(r,k[1]));
			}
			break;
		case RECTANGLE:
			if(1 == ndim){
				area = 2*s->vtab.rectangle.halfwidth[0];
				z[0] = Sinc(2*k[0]*s->vtab.rectangle.halfwidth[0]);
			}else{
				area = 4*s->vtab.rectangle.halfwidth[0]*s->vtab.rectangle.halfwidth[1];
				z[0] = Sinc(2*k[0]*s->vtab.rectangle.halfwidth[0])*Sinc(2*k[1]*s->vtab.rectangle.halfwidth[1]);
			}
			break;
		case POLYGON:
			{
				area = polygon_area(s->vtab.polygon.n_vertices, s->vtab.polygon.vertex);
				if(DC){
					z[0] = 1;
					z[1] = 0;
				}else{
					/* For k != 0,
					 * S(k) = i/|k|^2 * Sum_{i=0,n-1} z.((v_{i+1}-v_{i}) x k) j0(k.(v_{i+1}-v_{i})/2) e^{ik.(v_{i+1}+v_{i})/2}
					 */
					int p,q;
					double num, pa;
					double rc[2], u[2];
					for(p=s->vtab.polygon.n_vertices-1,q=0; q < s->vtab.polygon.n_vertices; p = q++){
						u[0] = s->vtab.polygon.vertex[2*q+0]-s->vtab.polygon.vertex[2*p+0];
						u[1] = s->vtab.polygon.vertex[2*q+1]-s->vtab.polygon.vertex[2*p+1];
						rc[0] = 0.5*(s->vtab.polygon.vertex[2*q+0]+s->vtab.polygon.vertex[2*p+0]);
						rc[1] = 0.5*(s->vtab.polygon.vertex[2*q+1]+s->vtab.polygon.vertex[2*p+1]);

						num = (u[0]*k[1]-u[1]*k[0]) * Sinc(k[0]*u[0]+k[1]*u[1]);
						pa = -2*M_PI*(k[0]*rc[0]+k[1]*rc[1]);

						// Multiplication by i means we mess up the order here
						z[0] += num * sin(pa);
						z[1] -= num * cos(pa);
					}
					// Our k lacks a 2pi factor
					//z[0] /= 2*M_PI*(k[0]*k[0]+k[1]*k[1])*area;
					//z[1] /= 2*M_PI*(k[0]*k[0]+k[1]*k[1])*area;
					area = 1;
					z[0] /= 2*M_PI*(k[0]*k[0]+k[1]*k[1]);
					z[1] /= 2*M_PI*(k[0]*k[0]+k[1]*k[1]);
				}
			}
			break;
		default:
			area = 0;
			break;
		}
		//f += dval*z*phase/area;
		//f += dval*(z[0]+i*z[1])*(cos(phase_angle)+i*sin(phase_angle)) / area;
		//f += dval*( z[0]*cos-z[1]*sin + i*(z[1]*cos+z[0]*sin) ) / area;
		{
			double cpa = cos(phase_angle);
			double spa = sin(phase_angle);
			double t[2] = {area, area};
			if(DC){
				t[0] *= cpa;
				t[1] *= spa;
			}else{
				t[0] *= ( z[0]*cpa-z[1]*spa );
				t[1] *= ( z[1]*cpa+z[0]*spa );
			}
			f[0] += inv_size*(t[0]*dval[0]-t[1]*dval[1]);
			f[1] += inv_size*(t[0]*dval[1]+t[1]*dval[0]);
		}
	}
	return 0;
}
int Pattern_GetFourierTransform(
	const Pattern *p,
	const double *value,
	const double k[2],
	int ndim,
	double unit_cell_size,
	double f[2]
){
	return pattern_get_fourier_transform(p->nshapes, p->shapes, p->parent, value, k, ndim, unit_cell_size, f);
}

int pattern_discretize_cell(
	int nshapes, // number of shapes in the array 'shapes'
	const shape *shapes, // array of shapes, ordered by decreasing area
	const int *parent, // from Pattern_GetContainmentTree
	const double L[4],
	int nu, int nv, // size of grid; number of steps in x and y
	int iu, int iv,
	double *value
){
	int k;
	const double det = L[0]*L[3] - L[1]*L[2];

	const int dim = (0 == L[1] && 0 == L[2] && 0 == L[3]) ? 1 : 2;

	if(nshapes < 0){ return -1; }
	if(nshapes != 0 && NULL == shapes){ return -2; }
	if(NULL == parent){ return -3; }
	if(2 == dim && 0 == det){ return -4; }
	if(nu < 1){ return -5; }
	if(nv < 1){ return -6; }
	if(iu < 0 || iu >= nu){ return -7; }
	if(iv < 0 || iv >= nv){ return -8; }
	if(NULL == value){ return -9; }

	for(k = 0; k <= nshapes; ++k){
		value[k] = 0;
	}
	value[0] = 1;

	if(1 == dim){
		const double du = L[0] / (double)nu;
		const double p0 = L[0]*((double)iu/(double)nu - 0.5); // left edge
		for(k = 0; k < nshapes; ++k){
			double a = 0;
			const shape *s = &shapes[k];
			if(RECTANGLE != s->type){ return -2; }
			{
				const double l = p0 - s->center[0];
				const double h = s->vtab.rectangle.halfwidth[0];
				// Rectangle centered at origin with halfwidth h
				// Pixel left edge at l
				if(l+du <= -h || l >= h){ continue; } // pixel outside
				if(l >= -h && l+du <= h){ // pixel entirely contained
					a = 1.;
				}else if(l+du >= h){ // pixel intersects right edge of rect
					a = (h-l) / du;
				}else{ // pixel intersects right edge of rect
					a = 1. + (h+l) / du;
				}
				if(a > 0){
					value[k+1] += a;
					value[parent[k]+1] -= a;
				}
			}
		}
	}else{
		const double duv[4] = {
			L[0]/(double)nu, L[1]/(double)nu,
			L[2]/(double)nv, L[3]/(double)nv
			};
		const double inv_pixel_area = 1./(duv[0]*duv[3]-duv[1]*duv[2]);
		const double nuv[2] = {
			(double)iu/(double)nu - 0.5,
			(double)iv/(double)nv - 0.5
		};
		const double p0[2] = { // bottom left of pixel
			L[0]*nuv[0] + L[2]*nuv[1],
			L[1]*nuv[0] + L[3]*nuv[1]
		};

		for(k = 0; k < nshapes; ++k){
			double a = shape_get_intersection_area_quad(&shapes[k], p0, duv) * inv_pixel_area;
			if(a > 0){
				value[k+1] += a;
				value[parent[k]+1] -= a;
			}
		}
	}
	/*
	// Dumb sampling
	for(int i = 0; i < nx; ++i){
		for(int j = 0; j < ny; ++j){
			double p0[2] = {-0.5 + ((double)i+0.5)/nx, -0.5 + ((double)j+0.5)/ny}; // bottom left of pixel
			int index;
			int ret = Pattern_GetShape(nshapes, shapes, parent, p0, &index);
			if(ret == 0){
				grid[i+j*nx] = value[index + 1];
			}else if(ret == 1){
				grid[i+j*nx] = value[0];
			}
		}
	}
	*/

	return 0;
}
int Pattern_DiscretizeCell(
	const Pattern *p,
	const double L[4],
	int nu, int nv,
	int iu, int iv,
	double *value
){
	return pattern_discretize_cell(p->nshapes, p->shapes, p->parent, L, nu, nv, iu, iv, value);
}


/*
static void fprint_double(FILE *fp, double x){
	static const unsigned int buflen = 64;
	static const char fmt[] = "%.16g";
	static const char expsymb[] = "*^";
	static const unsigned int lexpsymb = sizeof(expsymb)/sizeof(char) - 1;
	int len;
	char buf[buflen];
	len = snprintf(buf, buflen, fmt, x);

	char *target = &buf[0];
	int freetarget = 0;
	if(len > (int)(buflen - lexpsymb)){
		target = (char*)malloc(sizeof(char)*(len+lexpsymb));
		len = snprintf(target, len+lexpsymb, fmt, x);
		freetarget = 1;
	}
	char *exppos = strchr(target, 'e');
	if(NULL != exppos){
		int expoff = exppos - target;
		int shift = lexpsymb -1;
		for(; len > expoff; --len){
			target[len+shift] = target[len];
		}
		for(len = 0; len < (int)lexpsymb; ++len){
			target[expoff+len] = expsymb[len];
		}
	}
	fputs(buf, fp);
	if(freetarget){ free(target); }
}

static void fprint_spmatrix(FILE *fp, int n, const double *M, const int *Mcol, int Mnd){
	int i, j;
	fprintf(fp, "SparseArray[{\n");
	for(i = 0; i < n; ++i){ // r += A*x
		fprintf(fp, "{%d,%d} -> ", i, i);
		fprint_double(fp, M[(Mnd+1)*i+0]);
		fprintf(fp, "\n");
		for(j = 0; j < Mnd; ++j){
			fprintf(fp, "{%d,%d} -> ", i, Mcol[Mnd*i+j]);
			fprint_double(fp, M[(Mnd+1)*i+j+1]);
			fprintf(fp, ",\n");
		}
	}
	fprintf(fp, "},2,0]\n");
}
*/
// the least-squares system needs a damping factor to remove the null-space.
// The null space is dimension 2*g where g is the genus of the surface.
// since our surface is always isomorphic to a 2-torus, the genus should be zero.
#define DAMPING_FACTOR 1e-8

#ifndef HAVE_LIBCHOLMOD

static double norm(int n, const double *x){
#ifdef USE_BLAS
	extern double dnrm2_(const int *n, const double *x, const int *incx);
	const int ione = 1;
	return dnrm2_(&n, x, &ione);
#else
	double ssq = 0, scale = 0;
	while(n --> 0){
		if(0 != *x){
			double temp = fabs(*x);
			if(scale < temp){
				double r = scale/temp;
				ssq = ssq*r*r + 1.;
				scale = temp;
			}else{
				double r = temp/scale;
				ssq += r*r;
			}
		}
		++x;
	}
	return scale*sqrt(ssq);
#endif
}
static double dot(int n, const double *x, const double *y){
#ifdef USE_BLAS
	extern double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
	const int ione = 1;
	return ddot_(&n, x, &ione, y, &ione);
#else
	double z = 0;
	while(n --> 0){
		z += (*x)*(*y);
		++x; ++y;
	}
	return z;
#endif
}

static void cg(int n, const double *M, const int *Mcol, int Mnd, const double *b, double *x, double anorm){
	const double tol = 1e-10*n;
	const int max_iter = 2*n;
	const double damp = DAMPING_FACTOR*anorm;
	int i, j, iter;
	double resid = 0, alpha = 0, beta = 0, rho = 0, rho_1 = 0, normb = 0;
	double *work = (double*)malloc(sizeof(double) * 4*n);
	double *p = work;
	double *z = p + n;
	double *q = z + n;
	double *r = q + n;

	normb = norm(n, b);
	for(i = 0; i < n; ++i){ // r += A*x
		r[i] = (damp+M[(Mnd+1)*i+0])*x[i];
		for(j = 0; j < Mnd; ++j){
			r[i] += M[(Mnd+1)*i+j+1]*x[Mcol[Mnd*i+j]];
		}
	}
	for(i = 0; i < n; ++i){ r[i] = b[i] - r[i]; } // r = b-A*x
//for(i = 0; i < n; ++i){ fprintf(stderr, " %e\n", r[i]); }
	if(0 == normb){ normb = 1; }
//fprintf(stderr, "normb = %e\n", normb);
	resid = norm(n, r);
	if(resid < tol*anorm || (resid = resid / normb) <= tol){
		//tol = resid;
		//max_iter = 0;
		free(work);
		return;
	}
//fprintf(stderr, "resid = %e\n", resid);
	for(iter = 0; iter < max_iter; ++iter){
		for(i = 0; i < n; ++i){ // z = M.solve(r);
			z[i] = r[i]/(M[(Mnd+1)*i+0]+damp);
		}
		rho = dot(n, r, z);
//fprintf(stderr, "rho = %e\n", rho);
		if(iter == 0){
			for(i = 0; i < n; ++i){ p[i] = z[i]; }
		}else{
			beta = rho / rho_1;
			for(i = 0; i < n; ++i){ p[i] = z[i] + beta*p[i]; }
		}

		for(i = 0; i < n; ++i){ // q += A*p
			q[i] = (damp+M[(Mnd+1)*i+0])*p[i];
			for(j = 0; j < Mnd; ++j){
				q[i] += M[(Mnd+1)*i+j+1]*p[Mcol[Mnd*i+j]];
			}
		}
		alpha = rho / dot(n, p, q);
//fprintf(stderr, "alpha = %e\n", alpha);
		for(i = 0; i < n; ++i){ x[i] += alpha*p[i]; }
		if(0 == (iter % 50)){
			for(i = 0; i < n; ++i){ // r = A*x
				r[i] = (damp+M[(Mnd+1)*i+0])*x[i];
				for(j = 0; j < Mnd; ++j){
					r[i] += M[(Mnd+1)*i+j+1]*x[Mcol[Mnd*i+j]];
				}
			}
			for(i = 0; i < n; ++i){ r[i] = b[i] - r[i]; } // r = b-A*x
		}else{
			for(i = 0; i < n; ++i){ r[i] -= alpha*q[i]; }
		}

		if((resid = norm(n, r) / normb) <= tol){
			//tol = resid;
			iter = max_iter;
			break;
		}
//fprintf(stderr, "resid = %e\n", resid);
		rho_1 = rho;
	}
	//tol = resid;
	free(work);
}
#endif // !HAVE_LIBCHOLMOD

#ifdef HAVE_LIBCHOLMOD
#include <cholmod.h>
#endif

static void sparse_linsolve_d(int n, const double *M, const int *Mcol, int Mnd, const double *b, double *x, double anorm){
#ifdef HAVE_LIBCHOLMOD
	//const double damp = DAMPING_FACTOR*anorm;
	int i, j, nnz;
	double *Ax;
	cholmod_sparse *A;
	cholmod_dense *xs, *bs;
	cholmod_factor *L;
	cholmod_common c;
	cholmod_start(&c);

	if(anorm){} // anorm is not used except with CG

	A = cholmod_allocate_sparse(n,n, n*Mnd, 0, 1, 1, CHOLMOD_REAL, &c);
	Ax = (double*)(A->x);
	nnz = 0;
	for(i = 0; i < n; ++i){
		switch(A->itype){
		case CHOLMOD_INT:
			((int*)A->p)[i] = nnz;
			((int*)A->i)[nnz] = i;
			break;
		case CHOLMOD_INTLONG:
			((UF_long*)A->p)[i] = nnz;
			((int*)A->i)[nnz] = i;
			break;
		case CHOLMOD_LONG:
			((UF_long*)A->p)[i] = nnz;
			((UF_long*)A->i)[nnz] = i;
			break;
		}
		Ax[nnz] = M[(Mnd+1)*i+0] ;//+ damp;
		nnz++;
		for(j = 0; j < Mnd; ++j){
			if(0 == M[(Mnd+1)*i+j+1]){ break; }
			if(Mcol[Mnd*i+j] <= i){
				Ax[nnz] = M[(Mnd+1)*i+j+1];
				switch(A->itype){
				case CHOLMOD_INT:
					((int*)A->i)[nnz] = Mcol[Mnd*i+j];
					break;
				case CHOLMOD_INTLONG:
					((int*)A->i)[nnz] = Mcol[Mnd*i+j];
					break;
				case CHOLMOD_LONG:
					((UF_long*)A->i)[nnz] = Mcol[Mnd*i+j];
					break;
				}
				nnz++;
			}
		}
	}
	switch(A->itype){
	case CHOLMOD_INT:
		((int*)A->p)[i] = nnz;
		break;
	case CHOLMOD_INTLONG:
		((UF_long*)A->p)[i] = nnz;
		break;
	case CHOLMOD_LONG:
		((UF_long*)A->p)[i] = nnz;
		break;
	}

	bs = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, &c);
	memcpy(bs->x, b, sizeof(double)*n);

	L = cholmod_analyze(A, &c);
	cholmod_factorize(A, L, &c);
	xs = cholmod_solve(CHOLMOD_A, L, bs, &c);
	memcpy(x, xs->x, sizeof(double)*n);
	cholmod_free_factor(&L, &c);
	cholmod_free_sparse(&A, &c);
	cholmod_free_dense(&xs, &c);
	cholmod_free_dense(&bs, &c);
	cholmod_finish(&c);
#else
	cg(n, M, Mcol, Mnd, b, x, anorm);
#endif
}

static int pattern_generate_flow_field_constraints(
	int type,
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double p0[2],
	const double dp[2],
	double *value
){
	int si, c = 0;
	*value = 0;

	if(type){} // currently unused
	if(parent){} // parent information is not currently used

	for(si = 0; si < nshapes; ++si){
		double t;
		int ct;
		ct = shape_get_tangent_cross_segment(&shapes[si], p0, dp, &t);
		if(ct > 0){
			*value += t;
			c += ct;
		}
	}
	return c;
}
static int pattern_generate_flow_field_constraints_tri(
	int type,
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double p0[2],
	const double dp[2],
	double *value0,
	double *value1
){
	int si, c = 0;
	*value0 = 0;
	*value1 = 0;

	if(type){} // type is currently not used
	if(parent){} // parent info is currently not used

	for(si = 0; si < nshapes; ++si){
		double t0, t1;
		int ct;
		ct = shape_get_tangent_cross_segment_tri(&shapes[si], p0, dp, &t0, &t1);
		if(ct > 0){
			*value0 += t0;
			*value1 += t1;
			c += ct;
		}
	}
	return c;
}

static double pattern_generate_flow_field_rect(
	int nshapes,
	const shape *shapes,
	const int *parent,
	int type,
	const double L[4],
	int nu, int nv,
	double *field
){
	const int N = nu*nv;
	const int N2 = 2*N;
	int i, j;

	double *work = (double*)malloc(sizeof(double)*N2*(5+2));
	double *x = work;
	double *b = x+N2;
	double *M = b+N2;
	int *iwork = (int*)malloc(sizeof(int)*N2*4);
	int *Mcol = iwork;

	// The 1-form Laplacian M has exactly 5 non-zero entries per row on a
	// rectangular grid. The constraint matrix Z is diagonal, with one
	// entry per edge with a shape intersection.

	// Sparse matrix storage format:
	// M is arranged in blocks of 5 elements, for each row,
	// in ascending row order.
	// The first element in each row is the diagonal entry.
	// The remaining entries are for off-diagonal entries.
	// Mcol gives the columns of these off-diagonal entries, so
	// it is arranged in blocks of 4 elements, for each row, in
	// the same order as M.

	// === Assemble the M matrix and RHS b ===
	for(i = 0; i < N2; ++i){ b[i] = 0; }
	for(i = 0; i < N2; ++i){ x[i] = 0; }
{
	// We index the edges by the cell order of `field', and within that,
	// the bottom x-edge comes first, then the left y-edge.
	const double dL[4] = {L[0]/nu,L[1]/nu, L[2]/nv,L[3]/nv};
	/*
	const double _dxdy = (double)nu*nv;
	const double dy_dx3 = (double)nu*nu*nu/(double)nv;
	const double dx_dy3 = (double)nv*nv*nv/(double)nu;
	*/
	const double _dxdy = 1;
	const double dy_dx3 = (double)nu*nu/(double)(nv*nv);
	const double dx_dy3 = (double)nv*nv/(double)(nu*nu);

	for(j = 0; j < nv; ++j){
		int jp1 = (j+1)%nv;
		int jm1 = (j+nv-1)%nv;
		for(i = 0; i < nu; ++i){
			int ip1 = (i+1)%nu;
			int im1 = (i+nu-1)%nu;
			int k = i+j*nu;

			// The u-edge
			M[5*(2*k+0)+0] = 2*(dy_dx3 + _dxdy);
			Mcol[4*(2*k+0)+0] = 2*(i  +jm1*nu)+0; M[5*(2*k+0)+1] = -_dxdy;
			Mcol[4*(2*k+0)+1] = 2*(im1+j  *nu)+0; M[5*(2*k+0)+2] = -dy_dx3;
			Mcol[4*(2*k+0)+2] = 2*(ip1+j  *nu)+0; M[5*(2*k+0)+3] = -dy_dx3;
			Mcol[4*(2*k+0)+3] = 2*(i  +jp1*nu)+0; M[5*(2*k+0)+4] = -_dxdy;

			// The v-edge
			M[5*(2*k+1)+0] = 2*(dx_dy3 + _dxdy);
			Mcol[4*(2*k+1)+0] = 2*(i  +jm1*nu)+1; M[5*(2*k+1)+1] = -dx_dy3;
			Mcol[4*(2*k+1)+1] = 2*(im1+j  *nu)+1; M[5*(2*k+1)+2] = -_dxdy;
			Mcol[4*(2*k+1)+2] = 2*(ip1+j  *nu)+1; M[5*(2*k+1)+3] = -_dxdy;
			Mcol[4*(2*k+1)+3] = 2*(i  +jp1*nu)+1; M[5*(2*k+1)+4] = -dx_dy3;
		}
	}
	/*
	for(j = 0; j < nv; ++j){
		for(i = 0; i < nu; ++i){
			int k, l;
			for(k = 0; k < 2; ++k){
				int row = 2*(i+j*nu)+k;
				printf("%d\t%d\t%f\n", row+1, row+1, M[5*row+0]);
				for(l = 0; l < 4; ++l){
					printf("%d\t%d\t%f\n", row+1, 1+Mcol[4*row+l], M[5*row+l+1]);
				}
			}
		}
	}
	for(i = 0; i < N2; ++i){ fprintf(stderr, " %f\n", b[i]); }
	*/
	// === Assemble the constraints and move them to RHS ===
{
	int found_multiple_xsects = 0;
	if(0 == type){
		for(j = 0; j < nv; ++j){
			for(i = 0; i < nu; ++i){
				int k = i+j*nu;
				int ei; // edge index;
				int si, sj;

				for(sj = -1; sj <= 1; ++sj){ for(si = -1; si <= 1; ++si){
				const double p0[2] = {
					(si-0.5 + (double)i/nu) * L[0] + (sj-0.5 + (double)j/nv) * L[2],
					(si-0.5 + (double)i/nu) * L[1] + (sj-0.5 + (double)j/nv) * L[3]
					}; // bottom left of pixel
				double cross;
				int c;

				for(ei = 0; ei < 2; ++ei){ const int nei = 1^ei;
					double pshifted[2];
					int row = 2*k+ei;
					pshifted[0] = p0[0]+0.5*dL[2*ei+0]-0.5*dL[2*nei+0];
					pshifted[1] = p0[1]+0.5*dL[2*ei+1]-0.5*dL[2*nei+1];
					c = 0; cross = 0;

					c = pattern_generate_flow_field_constraints(
						type,
						nshapes, shapes, parent,
						pshifted, &dL[2*nei],
						&cross);
//if(0 == si && 0 == sj) fprintf(stderr, "%d\t%d\t%d\t%d\t%f\t%f\n", i, j, 0 == ei ? c:0, 1 == ei ? c:0, 0 == ei ? cross:0, 1 == ei ? cross:0);
					if(ei){ cross = -cross; }
					if(c > 0){
						x[row] = 1;
						b[row] += cross;
						if(c > 1){
							found_multiple_xsects = 1;
						}
					}
				}
				}}
			}
//fprintf(stderr, "\n");
		}
	}else{ // assume type == 1
		// We must first change the periodic boundary conditions to natural boundaries.
	}
	for(i = 0; i < N2; ++i){ // apply the LS constraint modification to M
		M[5*i+0] += x[i];
		x[i] = b[i];
	}
}

	/*
	for(j = 0; j < nv; ++j){
		for(i = 0; i < nu; ++i){
			int k, l;
			for(k = 0; k < 2; ++k){
				int row = 2*(i+j*nu)+k;
				printf("%d\t%d\t%f\n", row+1, row+1, M[5*row+0]);
				for(l = 0; l < 4; ++l){
					printf("%d\t%d\t%f\n", row+1, 1+Mcol[4*row+l], M[5*row+l+1]);
				}
			}
		}
	}
	for(i = 0; i < N2; ++i){ fprintf(stderr, " %f\n", b[i]); }
	*/

	//fprint_spmatrix(stderr, N2, M, Mcol, 4);
	sparse_linsolve_d(N2, M, Mcol, 4, b, x, _dxdy);
{
	// === Reconstruct the vector field ===
	double max_field = 0;
	for(j = 0; j < nv; ++j){
		int jp1 = (j+1)%nv;
		for(i = 0; i < nu; ++i){
			int ip1 = (i+1)%nu;
			double af;
			field[2*(i+j*nu)+0] = 0.5*(x[2*(i+j*nu)+0]+x[2*(i+jp1*nu)+0]);
			field[2*(i+j*nu)+1] = 0.5*(x[2*(i+j*nu)+1]+x[2*(ip1+j*nu)+1]);
			af = hypot(field[2*(i+j*nu)+0],field[2*(i+j*nu)+1]);
			if(af > max_field){ max_field = af; }
		}
	}
	free(iwork);
	free(work);
	return max_field;

}}
}

static double pattern_generate_flow_field_tri(
	int nshapes,
	const shape *shapes,
	const int *parent,
	int type,
	const double L[4],
	int nu, int nv,
	double *field
){
	const int N = nu*nv;
	const int N3 = 3*N;
	int i, j;

	double *work = (double*)malloc(sizeof(double)*N3*(11+2));
	double *x = work;
	double *b = x+N3;
	double *M = b+N3;
	int *iwork = (int*)malloc(sizeof(int)*N3*10);
	int *Mcol = iwork;

	// The fundamental parallelogram is like
	//
	//          +-----------+        +-----------+         |
	//         / \         /          \         / \        |
	//        /| |\       /           |\       /|  \       |
	//       v     d     /              v     d     \      |
	//      /       \   /                \   /       \     |
	//     /         \ /                  \ /         \    |
	//    +-----u--->-+                    +-----u--->-+   |
	//      sd == 1                          sd == 0
	//
	// There are 3 edges per unit cell instead of 2 as in the rectangular
	// case. We store the u-directed edge, then the v-directed edge, then
	// the diagonal d. The diagonal can be in one of two configurations
	// depending on the angle between u and v. The orientation of the
	// diagonals is always upwards toward v (v-u or v+u).

	// The 1-form Laplacian M has exactly 11 non-zero entries per row on a
	// regular triangular grid. The constraint matrix Z is diagonal, with
	// one entry per edge with a shape intersection.

	// pre-compute Hodge star components
	// Note that due to the regularity of the grid, primal area = dual area
	// sd is 1 if d = v-u, 0 if d = v+u
	double dL[6] = {
		L[0] / nu, L[1] / nu,
		L[2] / nv, L[3] / nv,
		0, 0
	};
	const int sd = (dL[0]*dL[2] + dL[1]*dL[3]) > 0;

	const double cross = dL[0]*dL[3] - dL[1]*dL[2];
	const double cross2 = cross*cross;
	double r_u, r_v, r_d;
	if(sd){
		dL[4] = dL[2]-dL[0];
		dL[5] = dL[3]-dL[1];
		{
		double Cu = 0.5*(dL[2]*dL[2]+dL[3]*dL[3])/cross2 * -(dL[0]*dL[4] + dL[1]*dL[5]);
		double Cv = 0.5*(dL[0]*dL[0]+dL[1]*dL[1])/cross2 *  (dL[2]*dL[4] + dL[3]*dL[5]);
		r_u = 2*hypot((Cu-0.5)*dL[0] +  Cv     *dL[2], (Cu-0.5)*dL[1] +  Cv     *dL[3]);
		r_v = 2*hypot( Cu     *dL[0] + (Cv-0.5)*dL[2],  Cu     *dL[1] + (Cv-0.5)*dL[3]);
		r_d = 2*hypot((Cu-0.5)*dL[0] + (Cv-0.5)*dL[2], (Cu-0.5)*dL[1] + (Cv-0.5)*dL[3]);
		}
	}else{
		dL[4] = dL[2]+dL[0];
		dL[5] = dL[3]+dL[1];
		{
		double Cu = 0.5*(dL[4]*dL[4]+dL[5]*dL[5])/cross2 * -(dL[0]*dL[2] + dL[1]*dL[3]);
		double Cv = 0.5*(dL[0]*dL[0]+dL[1]*dL[1])/cross2 *  (dL[4]*dL[2] + dL[5]*dL[3]);
		r_u = 2*hypot((Cu-0.5)*dL[0] +  Cv     *dL[4], (Cu-0.5)*dL[1] +  Cv     *dL[5]);
		r_v = 2*hypot( Cu     *dL[0] + (Cv-0.5)*dL[4],  Cu     *dL[1] + (Cv-0.5)*dL[5]);
		r_d = 2*hypot((Cu-0.5)*dL[0] + (Cv-0.5)*dL[4], (Cu-0.5)*dL[1] + (Cv-0.5)*dL[5]);
		}
	}
/*
	const double r_uu = r_u*r_u;
	const double r_vv = r_v*r_v;
	const double r_dd = r_d*r_d;
	const double r_uv = r_u*r_v;
	const double r_du = r_u*r_d;
	const double r_dv = r_v*r_d;
*/
{
	const double _r_uu = r_u*r_u;
	const double _r_vv = r_v*r_v;
	const double _r_dd = r_d*r_d;
	const double _r_uv = r_u*r_v;
	const double _r_du = r_u*r_d;
	const double _r_dv = r_v*r_d;
	double _r_scale = _r_uu;
	if(_r_vv > _r_scale){ _r_scale = _r_vv; }
	if(_r_dd > _r_scale){ _r_scale = _r_dd; }
	_r_scale = 1./_r_scale;
{
	const double r_uu = _r_uu * _r_scale;
	const double r_vv = _r_vv * _r_scale;
	const double r_dd = _r_dd * _r_scale;
	const double r_uv = _r_uv * _r_scale;
	const double r_du = _r_du * _r_scale;
	const double r_dv = _r_dv * _r_scale;

	// Sparse matrix storage format:
	// M is arranged in blocks of 11 elements, for each row,
	// in ascending row order.
	// The first element in each row is the diagonal entry.
	// The remaining entries are for off-diagonal entries.
	// Mcol gives the columns of these off-diagonal entries, so
	// it is arranged in blocks of 10 elements, for each row, in
	// the same order as M.

	for(i = 0; i < N3; ++i){ b[i] = 0; }
	for(i = 0; i < N3; ++i){ x[i] = 0; }
	memset(M, 0, sizeof(double)*N3*11);
{
	// Neighborhood information
	const signed char neigh_sd[] = {
		// if sd
		//   u edge
		//    The curl(curl(.)) part is
		//    [ (u_this + d_this - v_this) - (-d_below - u_this + v_belowright) ] / area(uvd)
		//    The grad(div(.)) part is
		//    r_u / area(dual) [ (r_u*(u_this-u_left) + r_v*(v_this-v_below) + r_d*(d_left-d_below))
		//                     - (r_u*(u_right-u_this) + r_v*(v_right-v_belowright) + r_d*(d_this-d_belowright)) ]
		//     Curl-loop 1-rings
		1,  1, -1,  1*-1, // v edge, belowright, pos (neg curl side)
		1,  0,  0, -1*+1, // v edge, this cell, neg (pos curl side)
		2,  0,  0,  1*+1, // d edge, this cell, pos (pos curl side)
		2,  0, -1, -1*-1, // d edge, below, neg (neg curl side)
		//     The div stars
		1,  1, -1, -1*-1, // v, belowright, neg (sink side, neg)
		1,  0,  0,  1*+1, // v edge, this, pos (source side, pos)
		2,  0,  0,  1*-1, // d, this, pos (sink side, neg)
		2,  0, -1, -1*+1, // d edge, below, neg (source side, pos)

		0, -1,  0, -1*+1, // u edge, left, neg (source side, pos)
		2, -1,  0,  1*+1, // d edge, left, pos (source side, pos)
		1,  0, -1, -1*+1, // v edge, below, neg (source side, pos)
		0,  1,  0,  1*-1, // u, right, pos (sink side, neg)
		1,  1,  0,  1*-1, // v, right, pos (sink side, neg)
		2,  1, -1, -1*-1, // d, belowright, neg (sink side, neg)

		//   v edge
		//    The curl(curl(.)) part is
		//    [ (v_this - u_aboveleft - d_left) - (u_this + d_this - v_this) ] /area(uvd)
		//    The grad(div(.)) part is
		//    r_v / area(dual) [ (r_u*(u_this-u_left) + r_v*(v_this-v_below) + r_d*(d_left-d_below))
		//                     - (r_u*(u_above-u_aboveleft) + r_v*(v_above-v_this) + r_d*(d_aboveleft-d_this)) ]
		//     Curl-loop 1-rings
		0, -1,  1, -1*+1, // u, aboveleft, neg (pos curl side)
		0,  0,  0,  1*-1, // u, this, pos (neg curl side)
		2,  0,  0,  1*-1, // d, this, pos (neg curl side)
		2, -1,  0, -1*+1, // d, left, neg (pos curl side)
		//     The div stars
		0, -1,  1, -1*-1, // u, aboveleft, neg (sink side, neg)
		0,  0,  0,  1*+1, // u, this, pos (source side, pos)
		2,  0,  0, -1*-1, // d, this, neg (sink side, neg)
		2, -1,  0,  1*+1, // d, left, pos (source side, pos)

		0, -1,  0, -1*+1, // u, left, neg (source side, pos)
		1,  0, -1, -1*+1, // v, below, neg (source side, pos)
		2,  0, -1, -1*+1, // d, below, neg (source side, pos)
		0,  0,  1,  1*-1, // u, above, pos (sink side, neg)
		1,  0,  1,  1*-1, // v, above, pos (sink side, neg)
		2, -1,  1,  1*-1, // d, aboveleft, pos (sink side, neg)

		//   The d-edge
		//    The curl(curl(.)) part is
		//    [ (u_this + d_this - v_this) - (-d_this + v_right - u_above) ] /area(uvd)
		//    The grad(div(.)) part is
		//    r_d / area(dual) [ (r_u*(u_right-u_this) + r_v*(v_right-v_belowright) + r_d*(d_this-d_belowright))
		//                     - (r_u*(u_above-u_aboveleft) + r_v*(v_above-v_this) + r_d*(d_aboveleft-d_this) ) ]
		//     Curl-loop 1-rings
		0,  0,  1, -1*-1, // u, above, neg (neg curl side)
		0,  0,  0,  1*+1, // u, this, pos (pos curl side)
		1,  0,  0, -1*+1, // v, this, neg (pos curl side)
		1,  1,  0,  1*-1, // v, right, pos (neg curl side)
		//     The div stars
		0,  0,  1,  1*-1, // u, above, pos (sink side, neg)
		0,  0,  0, -1*+1, // u, this, neg (source side, pos)
		1,  0,  0, -1*-1, // v, this, neg (sink side, neg)
		1,  1,  0,  1*+1, // v, right, pos (source side, pos)

		0,  1,  0,  1*+1, // u, right, pos (source side, pos)
		1,  1, -1, -1*+1, // v, belowright, neg (source side, pos)
		2,  1, -1, -1*+1, // d, belowright, neg (source side, pos)
		0, -1,  1, -1*-1, // u, aboveleft, neg (sink side, neg)
		1,  0,  1,  1*-1, // v, above, pos (sink side, neg)
		2, -1,  1,  1*-1, // d, aboveleft, pos (sink side, neg)

		// if !sd
		//   u edge
		//    The curl(curl(.)) part is
		//    [ (u_this + v_right - d_this) - (d_below - u_this - v_below) ] / area(uvd)
		//    The grad(div(.)) part is
		//    r_u / area(dual) [ (r_u*(u_this-u_left) + r_v*(v_this-v_below) + r_d*(d_this-d_belowleft))
		//                     - (r_u*(u_right-u_this) + r_v*(v_right-v_belowright) + r_d*(d_right-d_below)) ]
		//     Curl-loop 1-rings
		1,  1,  0,  1*+1, // v, right, pos (pos curl side)
		1,  0, -1, -1*-1, // v, below, neg (neg curl side)
		2,  0, -1,  1*-1, // d, below, pos (neg curl side)
		2,  0,  0, -1*+1, // d, this, neg (pos curl side)
		//     The div stars
		1,  1,  0,  1*-1, // v, right, pos (sink side, neg)
		1,  0, -1, -1*+1, // v, below, neg (source side, pos)
		2,  0, -1, -1*-1, // d, below, neg (sink side, neg)
		2,  0,  0,  1*+1, // d, this, pos (source side, pos)

		0, -1,  0, -1*+1, // u, left, neg (source side, pos)
		1,  0,  0,  1*+1, // v, this, pos (source side, pos)
		2, -1, -1, -1*+1, // d, belowleft, neg (source side, pos)
		0,  1,  0,  1*-1, // u, right, pos (sink side, neg)
		1,  1, -1, -1*-1, // v, belowright, neg (sink side, neg)
		2,  1,  0,  1*-1, // d, right, pos (sink side, neg)

		//   v edge
		//    The curl(curl(.)) part is
		//    [ (v_this - d_left + u_left) - (-u_above + d_this - v_this) ] /area(uvd)
		//    The grad(div(.)) part is
		//    r_v / area(dual) [ (r_u*(u_this-u_left) + r_v*(v_this-v_below) + r_d*(d_this-d_belowleft))
		//                     - (r_u*(u_above-u_aboveleft) + r_v*(v_above-v_this) + r_d*(d_above-d_left)) ]
		//     Curl-loop 1-rings
		0,  0,  1, -1*-1, // u, above, neg (neg curl side)
		0, -1,  0,  1*+1, // u, left, pos (pos curl side)
		2, -1,  0, -1*+1, // d, left, neg (pos curl side)
		2,  0,  0,  1*-1, // d, this, pos (neg curl side)
		//     The div stars
		0,  0,  1,  1*-1, // u, above, pos (sink side neg)
		0, -1,  0, -1*+1, // u, left, neg (source side pos)
		2, -1,  0, -1*-1, // d, left, neg (sink side neg)
		2,  0,  0,  1*+1, // d, this, pos (source side pos)

		0,  0,  0,  1*+1, // u, this, pos (source side pos)
		1,  0, -1, -1*+1, // v, below, neg (source side pos)
		2, -1, -1, -1*+1, // d, belowleft, neg (source side pos)
		0, -1,  1, -1*-1, // u, aboveleft, neg (sink side neg)
		1,  0,  1,  1*-1, // v, above, pos (sink side neg)
		2,  0,  1,  1*-1, // d, above, pos (sink side neg)

		//   The d-edge
		//    The curl(curl(.)) part is
		//    [ (-u_above + d_this - v_this) - (u_this + v_right - d_this) ] /area(uvd)
		//    The grad(div(.)) part is
		//    r_d / area(dual) [ (r_u*(u_this-u_left) + r_v*(v_this-v_below) + r_d*(d_this-d_belowleft))
		//                     - (r_u*(u_aboveright-u_above) + r_v*(v_aboveright-v_right) + r_d*(d_aboveright-d_this) ) ]
		//     Curl-loop 1-rings
		0,  0,  1, -1*+1, // u, above, neg (pos curl side)
		0,  0,  0,  1*-1, // u, this, pos (neg curl side)
		1,  1,  0,  1*-1, // v, right, pos (neg curl side)
		1,  0,  0, -1*+1, // v, this, neg (pos curl side)
		//     The div stars
		0,  0,  1, -1*-1, // u, above, neg (sink side neg)
		0,  0,  0,  1*+1, // u, this, pos (source side pos)
		1,  1,  0, -1*-1, // v, right, neg (sink side neg)
		1,  0,  0,  1*+1, // v, this, pos (source side pos)

		0, -1,  0, -1*+1, // u, left, neg (source side pos)
		1,  0, -1, -1*+1, // v, below, neg (source side pos)
		2, -1, -1, -1*+1, // d, belowleft, neg (source side pos)
		0,  1,  1,  1*-1, // u, aboveright, pos (sink side neg)
		1,  1,  1,  1*-1, // v, aboveright, pos (sink side neg)
		2,  1,  1,  1*-1  // d, aboveright, pos (sink side neg)
	};
	const signed char iconstraint_sd[] = { // only the last column in this table is unique
		// if sd
		//    u-edge
		1,  1, -1, -1, // v edge, belowright
		1,  0,  0, -1, // v edge, this cell
		2,  0,  0, -1, // d edge, this cell
		2,  0, -1, -1, // d edge, below
		//    v-edge
		0, -1,  1,  1, // u, aboveleft
		0,  0,  0,  1, // u, this
		2,  0,  0, -1, // d, this
		2, -1,  0, -1, // d, left
		//    d-edge
		0,  0,  1,  1, // u, above
		0,  0,  0,  1, // u, this
		1,  0,  0,  1, // v, this
		1,  1,  0,  1, // v, right
		// if !sd
		//    u-edge
		1,  1,  0, -1, // v, right
		1,  0, -1, -1, // v, below
		2,  0, -1, -1, // d, below
		2,  0,  0, -1, // d, this
		//    v-edge
		0,  0,  1,  1, // u, above
		0, -1,  0,  1, // u, left
		2, -1,  0,  1, // d, left
		2,  0,  0,  1, // d, this
		//    d-edge
		0,  0,  1,  1, // u, above
		0,  0,  0,  1, // u, this
		1,  1,  0, -1, // v, right
		1,  0,  0, -1, // v, this
	};
	const signed char *neigh = (sd ? neigh_sd : neigh_sd+3*14*4);
	const signed char *iconstraint = (sd ? iconstraint_sd : iconstraint_sd+3*4*4);
	const double hodge[9] = {
		r_uu, r_uv, r_du,
		r_uv, r_vv, r_dv,
		r_du, r_dv, r_dd
	};

	for(j = 0; j < nv; ++j){
		int jp1 = (j+1)%nv;
		int jm1 = (j+nv-1)%nv;
		int j3[3] = {jm1,j,jp1};
		for(i = 0; i < nu; ++i){
			int ip1 = (i+1)%nu;
			int im1 = (i+nu-1)%nu;
			int i3[3] = {im1,i,ip1};
			int k = i+j*nu;
			int ei;

			for(ei = 0; ei < 3; ++ei){
				int ci;
				M[11*(3*k+ei)+0] = 2*(1 + hodge[ei+ei*3]);
				// Add curl contribution
				for(ci = 0; ci < 4; ++ci){
					const signed char *edge = &neigh[(ci+ei*14)*4];
					Mcol[10*(3*k+ei)+ci] = 3*(i3[1+edge[1]]+j3[1+edge[2]]*nu)+edge[0];
					M[11*(3*k+ei)+ci+1] = neigh[3+(ci+ei*14)*4];
				}
				for(ci = 0; ci < 10; ++ci){
					const signed char *edge = &neigh[(4+ci+ei*14)*4];
					Mcol[10*(3*k+ei)+ci] = 3*(i3[1+edge[1]]+j3[1+edge[2]]*nu)+edge[0];
					M[11*(3*k+ei)+ci+1] += edge[3] * hodge[edge[0]+ei*3];
				}
			}
		}
	}

	if(0 == type){
		for(j = 0; j < nv; ++j){
			int jp1 = (j+1)%nv;
			int jm1 = (j+nv-1)%nv;
			int j3[3] = {jm1,j,jp1};
			for(i = 0; i < nu; ++i){
				int ip1 = (i+1)%nu;
				int im1 = (i+nu-1)%nu;
				int i3[3] = {im1,i,ip1};

				int si, sj;
				for(sj = -1; sj <= 1; ++sj){ for(si = -1; si <= 1; ++si){
				int ei; // edge index;
				double p0[2] = {
					(si-0.5 + (double)i/nu) * L[0] + (sj-0.5 + (double)j/nv) * L[2],
					(si-0.5 + (double)i/nu) * L[1] + (sj-0.5 + (double)j/nv) * L[3]
					}; // bottom left of pixel
				double cross[6];
				int icross = 0;
				for(ei = 0; ei < 3; ++ei){
					if(sd && 2 == ei){ // for d, the last edge, we can change p0 without worry
						p0[0] += dL[0];
						p0[1] += dL[1];
					}
					{
					int c = pattern_generate_flow_field_constraints_tri(
						type,
						nshapes, shapes, parent,
						p0, &dL[2*ei],
						&cross[2*ei+0], &cross[2*ei+1]);
					if(c > 0){
						icross |= (1<<ei);
					}
					}
				}
				// We store the alterations to the diagonal of M in x for now
				for(ei = 0; ei < 3; ++ei){
					if(icross & (1 << ei)){
						int k;
						for(k = 0; k < 4; ++k){
							const signed char *edge = &iconstraint[(k+ei*4)*4];
							int row = 3*(i3[1+edge[1]]+j3[1+edge[2]]*nu)+edge[0];
							int ti = (k&1);
							x[row] = 1;
							b[row] += edge[3] * cross[2*ei+ti];
						}
					}
				}
				}}
			}
		}
	}else{ // assume type == 1
		// We must first change the periodic boundary conditions to natural boundaries.
	}
	for(i = 0; i < N3; ++i){ // apply the LS constraint modification to M
		M[11*i+0] += x[i];
		x[i] = b[i];
	}

	/*
	printf("A=SparseArray[{");
	for(j = 0; j < nv; ++j){
		for(i = 0; i < nu; ++i){
			int k, l;
			for(k = 0; k < 3; ++k){
				int row = 3*(i+j*nu)+k;
				//printf("%d\t%d\t%f\n", row+1, row+1, M[11*row+0]);
				printf("{%d,%d}->%f,\n", row+1, row+1, M[11*row+0]);
				for(l = 0; l < 10; ++l){
					//printf("%d\t%d\t%f\n", row+1, 1+Mcol[10*row+l], M[11*row+l+1]);
					printf("{%d,%d}->%f,\n", row+1, 1+Mcol[10*row+l], M[11*row+l+1]);
				}
			}
		}
	}
	printf("{_,_}->0}];\n");
	for(i = 0; i < N3; ++i){ fprintf(stderr, " %f\n", b[i]); }
	*/
	//fprintf(stderr, "Starting CG\n");
	sparse_linsolve_d(N3, M, Mcol, 10, b, x, M[0]);
	//fprintf(stderr, "Ending CG\n");
	//for(i = 0; i < N3; ++i){ x[i] = b[i]; }
{
	// === Reconstruct the vector field ===
	double max_field = 0;
	for(j = 0; j < nv; ++j){
		int jp1 = (j+1)%nv;
		for(i = 0; i < nu; ++i){
			int ip1 = (i+1)%nu;
			double af;
			// Use proper Whitney 1-form interpolation
			if(sd){
				double cuv = x[3*(i+j*nu)+0]+x[3*(i+j*nu)+1]+x[3*(ip1+j*nu)+1]+x[3*(i+jp1*nu)+0];
				double cd  = x[3*(i+j*nu)+2];
				field[2*(i+j*nu)+0] = -0.5*cuv*(-dL[5]) + cd*(-dL[1]-dL[3]);
				field[2*(i+j*nu)+1] = -0.5*cuv*( dL[4]) + cd*( dL[0]+dL[2]);
			}else{
				double cuv = x[3*(i+j*nu)+1]-x[3*(i+j*nu)+0]+x[3*(ip1+j*nu)+1]-x[3*(i+jp1*nu)+0];
				double cd  = x[3*(i+j*nu)+2];
				field[2*(i+j*nu)+0] = 0.5*cuv*(-dL[5]) + cd*(-(dL[1]-dL[3]));
				field[2*(i+j*nu)+1] = 0.5*cuv*( dL[4]) + cd*( (dL[0]-dL[2]));
			}
			af = hypot(field[2*(i+j*nu)+0],field[2*(i+j*nu)+1]);
			if(af > max_field){ max_field = af; }
		}
	}

	free(iwork);
	free(work);
	return max_field;
}}}}
}
int pattern_generate_flow_field(
	int nshapes,
	const shape *shapes,
	const int *parent,
	int type,
	const double L[4],
	int nu, int nv,
	double *field
){
	int i, j;
	if(nshapes < 0){ return -1; }
	if(nshapes != 0 && NULL == shapes){ return -2; }
	if(NULL == parent){ return -3; }
	if(0 != type && 1 != type){ return -4; }
	if(NULL == L){ return -5; }
	if(nu < 1){ return -6; }
	if(nv < 1){ return -7; }
	if(NULL == field){ return -8; }


	if(0 == nshapes){
		for(j = 0; j < nv; ++j){
			for(i = 0; i < nu; ++i){
				field[2*(i+j*nu)+0] = 0;
				field[2*(i+j*nu)+1] = 0;
			}
		}
		return 0;
	}
{
	// Determine if we can get away with a simpler Laplacian formulation
	int lattice_type = 0;

	double u_len = hypot(L[0],L[1]);
	double v_len = hypot(L[2],L[3]);
	double uv = L[0]*L[2] + L[1]*L[3];
	double uv_max = ((u_len > v_len) ? u_len : v_len);
	if(fabs(uv / uv_max) < DBL_EPSILON){
		lattice_type |= 1;
	}
{
	double max_field = 1;
	if(0 != (lattice_type&1)){ // square/rectangular lattice
		// The u and v direction Laplacians are decoupled
		max_field = pattern_generate_flow_field_rect(
			nshapes, shapes, parent,
			type, L, nu, nv, field);
	}else{
		// Use the natural induced triangular mesh
		// This will create problems when u/nu is very different in
		// length from v/nv.
		max_field = pattern_generate_flow_field_tri(
			nshapes, shapes, parent,
			type, L, nu, nv, field);
	}

	max_field = 1./max_field;
	for(j = 0; j < nv; ++j){
		for(i = 0; i < nu; ++i){
			field[2*(i+j*nu)+0] *= max_field;
			field[2*(i+j*nu)+1] *= max_field;
		}
	}
	return 0;
}}
}
int Pattern_GenerateFlowField(
	const Pattern *p,
	int type,
	const double L[4],
	int nu, int nv,
	double *value
){
	return pattern_generate_flow_field(p->nshapes, p->shapes, p->parent, type, L, nu, nv, value);
}
