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

#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

typedef struct Interpolator_* Interpolator;
typedef enum Interpolator_type_{
	Interpolator_LINEAR,
	Interpolator_CUBIC_SPLINE,
	Interpolator_CUBIC_HERMITE_SPLINE
} Interpolator_type;

Interpolator Interpolator_New(int n, int ny, double *xy, Interpolator_type type);
void Interpolator_Destroy(Interpolator I);
double* Interpolator_Get(const Interpolator I, double x, int *ny);

#endif /* _INTERPOLATOR_H_ */
