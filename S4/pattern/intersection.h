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

#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

// This contains a bunch of routines for computing intersections between
// shapes, and various other computational geometry functions.

// Computes sqrt(x*x + y*y) in a stable way
inline double pythag2(double x, double y);

// Assumed no self-crossings, no holes, CCW orientation.
// v is length 2*n, xy-pairs.
double polygon_area(int n, const double *v);

// returns 0 if valid, 1 if self intersections, 2 if orientation incorrect
// currently does nothing (returns 0)
int polygon_check(int n, const double *v);

// Circle assumed at origin
double intersection_area_circle_triangle(
	double radius,
	// triangle vertices: {org, org+u, org+v} are in CCW orientation
	const double tri_org[2],
	const double tri_u[2],
	const double tri_v[2]
);

// Computes the intersection of (assumed) convex polygons P and Q.
// |P| == n, |Q| == m, returns Pi = intersection(P,Q), |Pi| == *ni
// On entry, *ni should be the size of Pi.
// If Pi is too small, returns -10, and *ni is the needed size.
// If P is in Q, returns 1.
// If Q is in P, returns 2.
// Otherwise, on success, returns 0.
// Returns -n if there is an error with the n-th argument.
int convex_polygon_intersection(
	int n, // n >= 3
	const double *P,
	int m, // m >= 3
	const double *Q,
	int *ni, // on input, size of Pi (length/2), on output, numer of points in Pi
	double *Pi // output intersection polygon
);

int polygon_triangulate(
	int n, const double *v, // the polygon, v is length 2*n
	int *t // length 3*(n-2), stores the triangle as triples of vertex indices into v
);

// return number of points in isect
// cross0 and cross1 only referenced when both are non-NULL and
// when cross is non-NULL
int intersection_polygon_segment(
	int n, // n >= 3
	const double *P,
	const double seg0[2],
	const double segd[2],
	double *isect, // length 2*n, intersection xy-pairs. Can be NULL.
	double *cross, // sum of cross products between tangent vector at intersection and segment, can be NULL
	double *cross0, // weighted sum of cross products, weighted by t-parameter of intersection
	double *cross1 // same as cross0 but for (1-t)
);

// returns number of intersections
int intersection_circle_segment(
	double radius,
	const double seg0[2],
	const double segd[2],
	double isect[4], // intersection xy-pairs
	double t[2] // can be NULL; the t-parameter of the segment's intersection(s)
);

#endif // _INTERSECTION_H_
