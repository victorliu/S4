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

#include <complex>
#include <cstddef>

typedef struct tag_fft_plan *fft_plan;

std::complex<double> *fft_alloc_complex(size_t n);
void fft_free(void *p);

fft_plan fft_plan_dft_2d(
	int n[2],
	std::complex<double> *in, std::complex<double> *out,
	int sign
);
void fft_plan_exec(const fft_plan plan);
void fft_plan_destroy(fft_plan plan);

int fft_next_fast_size(int n);

extern "C" void fft_init();
extern "C" void fft_destroy();
