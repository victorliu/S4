/* Copyright (C) 2009-2011, Stanford University
 * This file is part of S4r
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

#ifndef S4R_INTERSECTION_H_INCLUDED
#define S4R_INTERSECTION_H_INCLUDED

/* Returns area of polygon with n vertices stored in v */
double polygon_area(int n, const double *v);

/* Circle assumed at origin */
double intersection_area_circle_triangle(
	double radius,
	/* triangle vertices: {org, org+u, org+v} are in CCW orientation */
	const double tri_org[2],
	const double tri_u[2],
	const double tri_v[2]
);

/* Computes the intersection of (assumed) convex polygons P and Q.
 * |P| == n, |Q| == m, returns Pi = intersection(P,Q), |Pi| == *ni
 * On entry, *ni should be the size of Pi.
 * If Pi is too small, returns -10, and *ni is the needed size.
 * If P is in Q, returns 1.
 * If Q is in P, returns 2.
 * Otherwise, on success, returns 0.
 * Returns -n if there is an error with the n-th argument.
 * Note that ni = n+m is sufficient.
 */
int convex_polygon_intersection(
	int n, const double *P, /* n >= 3, number of vertices of P */
	int m, const double *Q, /* m >= 3, number of vertices of Q */
	int *ni, /* in: capacity of Pi (points), out: numer of points in Pi */
	double *Pi /* output intersection polygon */
);

#endif /* S4R_INTERSECTION_H_INCLUDED */
