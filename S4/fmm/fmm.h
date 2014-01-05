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

#ifndef _S4_FMM_H_
#define _S4_FMM_H_

#include <S4.h>
#include <complex>

typedef int (*FMMGetEpsilon)(const Simulation *S, const Layer *L, int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);

int FMMGetEpsilon_Kottke(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);
int FMMGetEpsilon_ClosedForm(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);
int FMMGetEpsilon_FFT(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);
int FMMGetEpsilon_Experimental(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);

// The following need to be called after FFT or ClosedForm
int FMMGetEpsilon_PolBasisNV(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);
int FMMGetEpsilon_PolBasisVL(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);
int FMMGetEpsilon_PolBasisJones(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv);

double GetLanczosSmoothingOrder(const Simulation *S);
double GetLanczosSmoothingFactor(double order, int power, double f[2]);

#endif // _S4_FMM_H_
