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

#include <config.h>

#include <cstdio>
#include <cmath>
#include <S4.h>
#include "fmm.h"

#include <limits>

extern "C" double Jinc(double x);

double GetLanczosSmoothingOrder(const S4_Simulation *S){
	double Lkuv = hypot(
		S->G[2*(S->n_G-1)+0]*S->Lk[0] + S->G[2*(S->n_G-1)+1]*S->Lk[2],
		S->G[2*(S->n_G-1)+0]*S->Lk[1] + S->G[2*(S->n_G-1)+1]*S->Lk[3]
	);
	double Lku = hypot(S->Lk[0],S->Lk[1]);
	double Lkv = hypot(S->Lk[2],S->Lk[3]);
	if(Lkv < Lku){ Lku = Lkv; } // Lku is the smaller length of the two Lk vectors
	return Lkuv + Lku;
}

/*
// non-overshooting smoothing (most of the time at least):
double GetLanczosSmoothingFactor(double mp1, int power, double f[2]){
	double fr = hypot(f[0],f[1]);
	double x = fr/mp1;
	double s = Jinc(x);
	return 0.5*(s*s + (1.-x))*s;
}
*/

double GetLanczosSmoothingFactor(double mp1, int power, double f[2]){
	double fr = hypot(f[0],f[1]);
	//fprintf(stderr, "%f\t%f\t%f\t%f\t%f\n", f[0], f[1], fr, mp1, 2*mp1+1-fr);
	double j = Jinc(fr/(2*mp1+1));
	return pow(j, power);
}

