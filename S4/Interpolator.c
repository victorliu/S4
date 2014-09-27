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

#include "Interpolator.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>

struct Interpolator_{
	double *xy;
	int ny;
	int n;
	Interpolator_type type;
	double *result;
};

Interpolator Interpolator_New(int n, int ny, double *xy, Interpolator_type type){
	Interpolator I = (Interpolator)malloc(sizeof(struct Interpolator_));
	I->n = n;
	I->ny = ny;
	I->result = (double*)malloc(sizeof(double)*(ny + (1+ny)*n));
	I->xy = I->result + ny;
	memcpy(I->xy, xy, sizeof(double)*(1+ny)*n);
	I->type = type;
	return I;
}

void Interpolator_Destroy(Interpolator I){
	if(NULL == I){ return; }
	if(NULL != I->result){ free(I->result); }
	free(I);
}

double* Interpolator_Get(const Interpolator I, double x, int *ny){
	if(NULL == I){ return NULL; }
	if(I->n < 1){ return NULL; }
	const int ld = 1+I->ny;
	int i, j;
	
	*ny = I->ny;
	switch(I->type){
	case Interpolator_CUBIC_SPLINE:
		{
			return I->result;
		}
	case Interpolator_CUBIC_HERMITE_SPLINE:
		{
			// Hermite spline interpolation with Kochanek-Bartels tangents
			/*
			static const double tension = 0, bias = 0, continuity = 0;
			static const double dpp = (1-tension)*(1+bias)*(1+continuity)/2;
			static const double dpm = (1-tension)*(1+bias)*(1-continuity)/2;
			static const double dmp = (1-tension)*(1-bias)*(1+continuity)/2;
			static const double dmm = (1-tension)*(1-bias)*(1-continuity)/2;
			*/
			static const double dpp = 0.5;
			static const double dpm = 0.5;
			static const double dmp = 0.5;
			static const double dmm = 0.5;

			int im1;
			for(im1 = 0, i = 1; i < I->n; im1 = i++){
				if((x < I->xy[i*ld]) || (i == I->n-1)){// x coordinate lies between i-1 and i
					int im2 = (im1 > 0 ? im1-1 : 0);
					int ip1 = (i < I->n-1 ? i+1 : i);
					int j;
					const double t = (x - I->xy[ld*im1+0]) / (I->xy[ld*i+0] - I->xy[ld*im1+0]);
					const double h00 = (2*t - 3)*t*t + 1;
					const double h10 = ((t - 2)*t + 1) * t;
					const double h01 = (-2*t + 3)*t*t;
					const double h11 = (t-1)*t*t;
					for(j = 0; j < I->ny; ++j){
						double m0 = dpp*(I->xy[ld*im1+j+1] - I->xy[ld*im2+j+1])/(I->xy[ld*im1+0] - I->xy[ld*im2+0]) + dmm*(I->xy[ld*i+j+1] - I->xy[ld*im1+j+1])/(I->xy[ld*i+0] - I->xy[ld*im1+0]);
						double m1 = dpm*(I->xy[ld*i+j+1] - I->xy[ld*im1+j+1])/(I->xy[ld*i+0] - I->xy[ld*im1+0]) + dmp*(I->xy[ld*ip1+j+1] - I->xy[ld*i+j+1])/(I->xy[ld*ip1+0] - I->xy[ld*i+0]);

						I->result[j] = h00 * I->xy[ld*im1+j+1] + h01 * I->xy[ld*i+j+1] + (I->xy[ld*i+0] - I->xy[ld*im1+0])*(h10*m0 + h11*m1);
					}
					break;
				}
			}
			return I->result;
		}
	default: // Interpolator_LINEAR
		{
			// piecewise linear interpolation
			for(i = 1; i < I->n; ++i){
				if((x < I->xy[i*ld]) || (i == I->n-1)){
					double t = (x - I->xy[(i-1)*ld]) / (I->xy[i*ld] - I->xy[(i-1)*ld]);
					for(j = 0; j < I->ny; ++j){
						I->result[j] = t * I->xy[1+j+i*ld] + (1-t) * I->xy[1+j+(i-1)*ld];
					}
					break;
				}
			}
			return I->result;
		}
	}
	return NULL;
}
