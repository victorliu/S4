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
#include "intersection.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double polygon_area(int n, const double *v){
	int p, q;
	double area = 0;
	for(p=n-1, q=0; q < n; p = q++){
		area += v[2*p+0]*v[2*q+1] - v[2*q+0]*v[2*p+1];
	}
	return 0.5*area;
}

inline double LeftTurn(const double a[2], const double b[2], const double c[2]){
	extern double orient2d(double*, double*, double*);
	return orient2d((double*)a, (double*)b, (double*)c);
}

inline int dsign(double x){
	if(0 == x){ return 0; }
	else if(x > 0){ return 1; }
	else{ return -1; }
}

int SegmentsIntersect(const double a[2], const double b[2], const double c[2], const double d[2], double *x, double *t){
	double c0, c1, c2, c3;
	c0 = LeftTurn(a,b,c);
	c1 = LeftTurn(a,b,d);
	if(0 == c0 && 0 == c1){ return 0; } // collinear -> no intersection
	if(dsign(c0) != dsign(c1)){
		c2 = LeftTurn(c, d, a);
		c3 = LeftTurn(c, d, b);
		if(dsign(c2) != dsign(c3)){
			if(NULL != x){
				double r;
				c1 = c0-c1; c3 = c2-c3;
				if(fabs(c1) > fabs(c3)){
					r = c0/c1;
					x[0] = c[0] + r*(d[0]-c[0]);
					x[1] = c[1] + r*(d[1]-c[1]);
				}else{
					r = c2/c3;
					x[0] = a[0] + r*(b[0]-a[0]);
					x[1] = a[1] + r*(b[1]-a[1]);
				}
				if(NULL != t){ *t = r; }
			}
			return 1;
		}else{ return 0; }
	}else{ return 0; }
}

int convex_polygon_intersection(
	int n, // n >= 3
	const double *P,
	int m, // m >= 3
	const double *Q,
	int *ni, // on input, size of Pi, on output, numer of points in Pi
	double *Pi // output intersection polygon
){
	int i, j;
	if(n < 3){ return -1; }
	if(NULL == P){ return -2; }
	if(m < 3){ return -3; }
	if(NULL == Q){ return -4; }
	if(NULL == ni){ return -5; }
	if(NULL == Pi){ return -6; }

	// Implementation of:
	// "A new linear algorithm for intersecting convex polygons"
	// Joseph O'Rourke, Chi-Bin Chien, Thomas Olson, and David Naddor
	// Computer Graphics and Image Processing 19, pp. 384-391 (1982)
	
	const int nPi = *ni; *ni = 0;
	
	
	int ip = 1, iq = 1;
	int ipp = 0, iqp = 0; // prev of ip and iq
	char inside = ' ';
	// record first intersection
	int first_xsected = 0;
	int ipf = n, iqf = m;
	int first_iter = 0;
	
	int Pi_full = 0;
	int iter;

	// First, a bounding box check
	{
		int iP_x_min = 0, iP_x_max = 0, iP_y_min = 0, iP_y_max = 0;
		int iQ_x_min = 0, iQ_x_max = 0, iQ_y_min = 0, iQ_y_max = 0;
		for(i = 1; i < n; ++i){
			if(P[2*i+0] < P[2*iP_x_min+0]){ iP_x_min = i; }
			if(P[2*i+1] < P[2*iP_y_min+1]){ iP_y_min = i; }
			if(P[2*i+0] > P[2*iP_x_max+0]){ iP_x_max = i; }
			if(P[2*i+1] > P[2*iP_y_max+1]){ iP_y_max = i; }
		}
		for(i = 1; i < m; ++i){
			if(Q[2*i+0] < Q[2*iQ_x_min+0]){ iQ_x_min = i; }
			if(Q[2*i+1] < Q[2*iQ_y_min+1]){ iQ_y_min = i; }
			if(Q[2*i+0] > Q[2*iQ_x_max+0]){ iQ_x_max = i; }
			if(Q[2*i+1] > Q[2*iQ_y_max+1]){ iQ_y_max = i; }
		}
		if(
			( Q[2*iQ_x_min+0] > P[2*iP_x_max+0] ) ||
			( P[2*iP_x_min+0] > Q[2*iQ_x_max+0] ) ||
			( Q[2*iQ_y_min+1] > P[2*iP_y_max+1] ) ||
			( P[2*iP_y_min+1] > Q[2*iQ_y_max+1] )
		){ return 0; }
	}

	for(iter = 0; iter <= 2*(m+n); ++iter){
//fprintf(stderr, "iter %d, ip = %d, iq = %d, inside = %c\n", iter, ip, iq, inside);
		double xp[2];
		if(SegmentsIntersect(&P[2*ipp],&P[2*ip],&Q[2*iqp],&Q[2*iq],xp, NULL)){
//fprintf(stderr, " xsect! %f,%f %f,%f %f,%f %f,%f\n", P[2*ipp+0],P[2*ipp+1],P[2*ip+0],P[2*ip+1],Q[2*iqp+0],Q[2*iqp+1],Q[2*iq+0],Q[2*iq+1]);
			if(first_xsected && first_iter != iter-1){ // if the first intersection was NOT found during the previous iteration
				if(ip == ipf && iq == iqf){ break; } // if this intersection is the same as the first intersection
			}
			if(*ni >= nPi){ Pi_full = 1; }
			if(!Pi_full){ Pi[2*(*ni)+0] = xp[0]; Pi[2*(*ni)+1] = xp[1]; (*ni)++; }
//fprintf(stderr, "  Adding %f,%f\n", Pi[2*((*ni)-1)+0], Pi[2*((*ni)-1)+1]);
			if(LeftTurn(&Q[2*iqp],&Q[2*iq],&P[2*ip]) >= 0){
				inside = 'P';
			}else{ inside = 'Q'; }
			
			if(!first_xsected){
				first_xsected = 1;
				ipf = ip; iqf = iq;
				first_iter = iter;
			}
		}
		xp[0] = P[2*ip+0] + (Q[2*iq+0] - P[2*ipp+0]);
		xp[1] = P[2*ip+1] + (Q[2*iq+1] - P[2*ipp+1]);
		if(LeftTurn(&Q[2*iqp],&Q[2*iq],xp)/*Cross(Q[2*iq]-Q[2*iqp],P[2*ip]-P[2*ipp])*/ >= 0){
			if(LeftTurn(&Q[2*iqp],&Q[2*iq],&P[2*ip]) >= 0){ // advance Q
				if(inside == 'Q'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = Q[2*iq+0]; Pi[2*(*ni)+1] = Q[2*iq+1]; (*ni)++; }
				}
				iqp = iq;
				iq = (iq+1)%m;
			}else{ // advance P
				if(inside == 'P'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = P[2*ip+0]; Pi[2*(*ni)+1] = P[2*ip+1]; (*ni)++; }
				}
				ipp = ip;
				ip = (ip+1)%n;
			}
		}else{
			if(LeftTurn(&P[2*ipp],&P[2*ip],&Q[2*iq]) >= 0){ // advance P
				if(inside == 'P'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = P[2*ip+0]; Pi[2*(*ni)+1] = P[2*ip+1]; (*ni)++; }
				}
				ipp = ip;
				ip = (ip+1)%n;
			}else{ // advance Q
				if(inside == 'Q'){
					if(*ni >= nPi){ Pi_full = 1; }
					if(!Pi_full){ Pi[2*(*ni)+0] = Q[2*iq+0]; Pi[2*(*ni)+1] = Q[2*iq+1]; (*ni)++; }
				}
				iqp = iq;
				iq = (iq+1)%m;
			}
		}
	}
	// At this point, either P in Q, Q in P, or they don't intersect
	if(*ni == 0){
		int flag = 1;
		for(j = 0; j < n; ++j){ // really we only need to check j == 0, but due to degeneracy, it is safest to check all
			for(i = 0; i < m; ++i){
				if(LeftTurn(&Q[2*i],&Q[2*((i+1)%m)], &P[2*j]) < 0){
					flag = 0; j = n+1; break;
				}
			}
		}
		if(flag){ // P in Q
			if(*ni+n >= nPi){ Pi_full = 1; }
			if(!Pi_full){
				for(i = 0; i < n; ++i){
					Pi[2*(*ni)+0] = P[2*i+0]; Pi[2*(*ni)+1] = P[2*i+1]; (*ni)++;
				}
				return 1;
			}
		}else{
			flag = 1;
			for(j = 0; j < m; ++j){ // really we only need to check j == 0, but due to degeneracy, it is safest to check all
				for(i = 0; i < n; ++i){
					if(LeftTurn(&P[2*i],&P[2*((i+1)%n)],&Q[2*j]) < 0){
						flag = 0; j = m+1; break;
					}
				}
			}
			if(flag){ // Q in P
				if(*ni+m >= nPi){ Pi_full = 1; }
				if(!Pi_full){
					for(i = 0; i < m; ++i){
						Pi[2*(*ni)+0] = Q[2*i+0]; Pi[2*(*ni)+1] = Q[2*i+1]; (*ni)++;
					}
					return 2;
				}
			}
		}
	}
	if(Pi_full){
		return -10;
	}else{
		return 0;
	}
}

// returns 1 if inside or on boundary, 0 otherwise
int TriangleContainsPoint(
	const double org[2], // triangle vertices are {org,org+u,org+v}, in CCW orientation
	const double u[2],
	const double v[2],
	const double p[2] // query point
){
	if(LeftTurn(org,u,p) >= 0){
		if(LeftTurn(org,v,p) <= 0){
			// treat org as origin, we have points u and v, just need p-org
			double x[2] = {p[0] - org[0], p[1] - org[1]};
			if(LeftTurn(u,v,x) >= 0){
				return 1;
			}
		}
	}
	return 0;
}

// circle radius r, chord length s
// r > 0, s > 0, 2*r > s
double CircularSectorArea(double r, double s){
	if(s <= 0){ return 0; }

	// area = area of circular wedge - area of triangle part
	// area = theta*r*r - s/2*sqrt(r^2 - (s/2)^2)
	s *= 0.5;
	// area = asin(s/r)*r*r - s*sqrt(r^2 - s^2)
	// area = r*s*[ asin(s/r)/(s/r) - sqrt(1-(s/r)^2) ]
	double x = s/r;
	if(x < (1./32.)){ // use taylor expansion for stuff in brackets
		static const double c2 = 2./3.;
		static const double c4 = 1./5.;
		static const double c6 = 3./28.;
		static const double c8 = 3./72.;
		x *= x;
		return r*s*(c2 + (c4 + (c6 + c8*x)*x)*x)*x;
	}else{
		return r*s*(asin(x)/x - sqrt((1.+x)*(1.-x)));
	}
}

double intersection_area_circle_triangle(
	double radius,
	// triangle vertices: {org, org+u, org+v} are in CCW orientation
	const double tri_org[2],
	const double tri_u[2],
	const double tri_v[2]
){
	const double origin[2] = {0,0};
	int i, j;
	const double iradius = 1./radius;
	const double iradius2 = 1./(radius*radius);
	int inside = 0; // bitfield of which triangle vertices are in circle, 4th bit is if circle center is in triangle
	
	double vert[6];
	vert[2*0+0] = tri_org[0];
	vert[2*0+1] = tri_org[1];
	vert[2*1+0] = (tri_org[0] + tri_u[0]);
	vert[2*1+1] = (tri_org[1] + tri_u[1]);
	vert[2*2+0] = (tri_org[0] + tri_v[0]);
	vert[2*2+1] = (tri_org[1] + tri_v[1]);
	
	double vert_org[6];
	vert_org[2*0+0] = 0;
	vert_org[2*0+1] = 0;
	vert_org[2*1+0] = tri_u[0];
	vert_org[2*1+1] = tri_u[1];
	vert_org[2*2+0] = tri_v[0];
	vert_org[2*2+1] = tri_v[1];
	
	double tri_r[3]; // distance from circle center of each vertex of triangle, normalized to radius
	for(i = 0; i < 3; ++i){
		tri_r[i] = hypot(vert[2*i+0], vert[2*i+1]) * iradius;
		if(tri_r[i] <= 1.0){ inside |= (1<<i); }
	}
	if(TriangleContainsPoint(tri_org, tri_u, tri_v, origin)){ inside |= (1<<3); }
	if((inside & 0x7) == 0x7){ // all triangle points in circle
		return 0.5*fabs(LeftTurn(origin,tri_u,tri_v));
	}

	double seg[6];
	seg[2*0+0] = tri_u[0];
	seg[2*0+1] = tri_u[1];
	seg[2*1+0] = (tri_v[0] - tri_u[0]);
	seg[2*1+1] = (tri_v[1] - tri_u[1]);
	seg[2*2+0] = -tri_v[0];
	seg[2*2+1] = -tri_v[1];
	
	double side_length[3]; // normalized to radius
	for(i = 0; i < 3; ++i){
		side_length[i] = hypot(seg[2*i+0], seg[2*i+1]) * iradius;
	}
	
	double seg_dot[3]; // p0.(p1-p0)/r^2
	for(i = 0; i < 3; ++i){
		seg_dot[i] = (vert[2*i+0]*seg[2*i+0] + vert[2*i+1]*seg[2*i+1]) * iradius2;
	}
	
	// Get intersections of each segment with each circle
	// segment 0 is org to u, segment 1 is u to v, segment 2 is v to org
	double xp[12]; // intersection points
	int nxp[3]; // number of intersections with each segment
	int nx = 0;
	for(i = 0; i < 3; ++i){
		int ip1 = (i+1)%3;
		int in0 = (inside & (1<<i)) ? 1 : 0;
		int in1 = (inside & (1<<ip1)) ? 1 : 0;
		if(in0 && in1){
			nxp[i] = 0;
		}else{
			// line:   x = p0 + t(p1-p0)
			// circle: x^2 = r^2
			// (p0 + t(p1-p0))^2 = r^2
			// t^2(p1-p0)^2 + 2t(p0).(p1-p0) + (p0)^2 - r^2 = 0
			// t^2 * side_length^2 + 2*t*seg_dot + tri_r^2 - 1 = 0
			// t = -seg_dot/side_length^2 +/- sqrt(seg_dot^2/side_length^4 - (tri_r^2-1)/side_length^2)
			double isl2 = 1./(side_length[i]*side_length[i]);
			double disc = (seg_dot[i]*seg_dot[i]*isl2 - (tri_r[i]*tri_r[i]-1)) * isl2;
			double t0 = -seg_dot[i]*isl2;
			double t, t2, tdist, t2dist;
			if(in0 != in1){
				// get the one intersection point
				nxp[i] = 1;
				nx += 1;

				if(disc < 0){ disc = 0; }
				disc = sqrt(disc);
				t = t0+disc;
				t2 = t0-disc;
				if(t > 0.5){ tdist = fabs(t-1.); }else{ tdist = fabs(t); }
				if(t2 > 0.5){ t2dist = fabs(t2-1.); }else{ t2dist = fabs(t2); }
				if(t2dist < tdist){ t = t2; }
				if(t < 0){ t = 0; }
				if(t > 1){ t = 1; }
				xp[2*2*i+0] = vert_org[2*i+0] + t*seg[2*i+0];
				xp[2*2*i+1] = vert_org[2*i+1] + t*seg[2*i+1];
			}else{
				// possibly 0 or 2 intersection points; we count 1 degenerate point as none
				if(disc <= 0){
					nxp[i] = 0;
				}else{
					disc = sqrt(disc);
					nxp[i] = 0;
					t = t0-disc;
					t2 = t0+disc;
					if(0 < t && t < 1 && 0 < t2 && t2 < 1){
						xp[2*(2*i+0)+0] = vert_org[2*i+0] + t*seg[2*i+0];
						xp[2*(2*i+0)+1] = vert_org[2*i+1] + t*seg[2*i+1];
						xp[2*(2*i+1)+0] = vert_org[2*i+0] + t*seg[2*i+0];
						xp[2*(2*i+1)+1] = vert_org[2*i+1] + t*seg[2*i+1];
						nxp[i] += 2;
						nx += 2;
					}
				}
			}
		}
	}
/*	
	printf("tri: (%f,%f) (%f,%f) (%f,%f)\n",
		vert[0], vert[1], vert[2], vert[3], vert[4], vert[5]);
	printf("rad: %f\n", radius);
	printf("inside = %d, nx = %d\n", inside, nx);
	printf("xp:");
	for(i = 0; i < 3; ++i){
		for(j = 0; j < nxp[i]; ++j){
			printf(" (%d,%d,{%f,%f})", i, j, tri_org[0]+xp[2*(2*i+j)+0], tri_org[1]+xp[2*(2*i+j)+1]);
		}
	}printf("\n");
*/
	if(0 == nx){ // either no intersection area, or triangle entirely in circle, or circle in triangle
		if((inside & 0x7) == 0x7){ // all triangle points in circle
			// we already dealt with this above; getting here would be an error.
			return -1;
		}else{ // either no intersection area, or circle in triangle
			if(inside & 0x8){ // triangle contains circle center, intersection area is either circle area or triangle area
				return M_PI*radius*radius;
			}else{
				return 0;
			}
		}
	}else if(2 == nx){
		// Either the 2 intersections are a single side or on two different sides
		if(nxp[0] < 2 && nxp[1] < 2 && nxp[2] < 2){ // on different sides
			// The area is determined by tracing
		}else{
			for(i = 0; i < 3; ++i){
				if(nxp[i] > 1){ break; }
			}
			// Either the circle is mostly inside with a wedge poking out a side
			// or the circle is mostly outside with a wedge poking inside
			double sector_area = CircularSectorArea(radius, hypot(xp[2*(2*i+1)+0]-xp[2*2*i+0],xp[2*(2*i+1)+1]-xp[2*2*i+1]));
			if(inside & (1 << 3)){
				// Area of circle minus a wedge
				return M_PI*radius*radius - sector_area;
			}else{
				return sector_area;
			}
		}
	}else if(4 == nx){
		// The area is determined by tracing
	}else if(6 == nx){
		// The area is determined by tracing
	}else{
		// There is no way we can get here
		return -1;
	}
	
	// At this point we expect to just trace out the intersection shape
	// The vertices of the intersection shape is either a triangle vertex
	// or a intersection point on a triangle edge.
	int vtype[6]; // 1 = triangle vertex, 0 = intersection point
	double vp[12];
	int nv = 0; // number of actual vertices
	
	for(i = 0; i < 3; ++i){
		if(inside & (1 << i)){
			vp[2*nv+0] = vert_org[2*i+0];
			vp[2*nv+1] = vert_org[2*i+1];
			vtype[nv++] = 1;
		}
		for(j = 0; j < nxp[i]; ++j){
			vp[2*nv+0] = xp[2*(2*i+j)+0];
			vp[2*nv+1] = xp[2*(2*i+j)+1];
			vtype[nv++] = 0;
		}
	}

	if(nv < 3){ // this should not be possible
		return -1;
	}
	
	// All neighboring points in v which are intersection points should have circular caps added
	double area = polygon_area(nv, vp);
	for(i = 0; i < nv; ++i){
		int im1 = i-1; if(im1 < 0){ im1 = nv-1; }
		if((0 == vtype[im1]) && (0 == vtype[i])){
			area += CircularSectorArea(radius, hypot(vp[2*i+0]-vp[2*im1+0],vp[2*i+1]-vp[2*im1+1]));
		}
	}
	return area;
}
