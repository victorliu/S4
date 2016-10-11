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

#include <cmath>
#include <S4.h>
#include "../RNP/TBLAS.h"
#ifdef HAVE_BLAS
# include "../RNP/TBLAS_ext.h"
#endif
#include "../RNP/LinearSolve.h"
#ifdef HAVE_LAPACK
# include "../RNP/LinearSolve_lapack.h"
#endif
#include "fmm.h"

#include <limits>
//#include <kiss_fft.h>
//#include <tools/kiss_fftnd.h>
#include <cstdio>
#include "fft_iface.h"
#include <cstring>

static void vfield_to_Jones_vector(double v[2], std::complex<double> j[2]){
	double tx = v[0];
	double ty = v[1];
	double abs_t = hypot(tx,ty);
	/*
	double r[2] = {-0.5+((double)ai[0]+0.5)/ngrid[0], -0.5+((double)ai[1]+0.5)/ngrid[1]};
	abs_t = pythag2(r[0],r[1]);
	tx = -r[1]/abs_t; ty = r[0]/abs_t;
	abs_t = 1;
	*/

	double theta;
	double eta;
	if(abs_t < std::numeric_limits<double>::epsilon()){
		tx = 1;
		ty = 0;
		theta = 0;
		eta = 0.25*M_PI;
	}else{
		theta = atan2(ty,tx);
		tx /= abs_t;
		ty /= abs_t;
		eta = abs_t;
		if(eta > 1){ eta = 1; }
		eta = 0.125*M_PI*(1+cos(M_PI*eta));
	}

	// \hat{parallel} = exp(i*theta)/abs_t * [t t^\perp] [cos eta; i sin eta]
	// = exp(i*theta)/abs_t * [ tx*cos(eta) - i*ty*sin(eta) ]
	//                        [ ty*cos(eta) + i*tx*sin(eta) ]
	double ceta = cos(eta);
	double seta = sin(eta);
	double cth = cos(theta);
	double sth = sin(theta);
	std::complex<double> phase = std::complex<double>(cth,sth);
	j[0] = phase * std::complex<double>(cth*ceta, -sth*seta);
	j[1] = phase * std::complex<double>(sth*ceta,  cth*seta);
}

// We assume that coming into this function, Epsilon_inv is already filled in
// and so we do not modify it. We assume that Epsilon2 is of the form
//   [ Epsilon     0    ]
//   [    0     Epsilon ]
// where the diagonal blocks are identical (for isotropic materials).
// This matrix is used to effectively compute (Dy,-Dx) from (Ey,-Ex).
// We need to form the coordinate transformation matrix
//   F = [  uy ux* ]
//       [ -ux uy* ]
// where u is the Jones vector field, and each of the 4 blocks are the
// Fourier transformed matrices.
int FMMGetEpsilon_PolBasisJones(const S4_Simulation *S, const S4_Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv){
	double mp1 = 0;
	int pwr = S->options.lanczos_smoothing_power;
	if(S->options.use_Lanczos_smoothing){
		mp1 = GetLanczosSmoothingOrder(S);
		S4_TRACE("I   Lanczos smoothing order = %f\n", mp1);
		mp1 *= S->options.lanczos_smoothing_width;
	}

	if(Epsilon_inv){} // prevent unused parameter warning

	const int n2 = 2*n;
	const int nn = n*n;
	const double unit_cell_size = Simulation_GetUnitCellSize(S);
	const int *G = S->G;
	const int ndim = (0 == S->Lr[2] && 0 == S->Lr[3]) ? 1 : 2;
	double *ivalues = (double*)S4_malloc(sizeof(double)*(2+10)*(L->pattern.nshapes+1));
	double *values = ivalues + 2*(L->pattern.nshapes+1);

	// Get all the dielectric tensors
	//bool have_tensor = false;
	for(int i = -1; i < L->pattern.nshapes; ++i){
		const S4_Material *M;
		if(-1 == i){
			M = &S->material[L->material];
		}else{
			M = &S->material[L->pattern.shapes[i].tag];
		}
		if(0 == M->type){
			std::complex<double> eps_temp(M->eps.s[0], M->eps.s[1]);
			//eps_temp = Simulation_GetEpsilonByIndex(S, L->pattern.shapes[i].tag);
			values[2*(i+1)+0] = eps_temp.real();
			values[2*(i+1)+1] = eps_temp.imag();
			eps_temp = 1./eps_temp;
			ivalues[2*(i+1)+0] = eps_temp.real();
			ivalues[2*(i+1)+1] = eps_temp.imag();
		}else{
			//have_tensor = true;
		}
	}

	// P will in the end store the LU decomposition of the F matrix
	std::complex<double> *P = Simulation_GetCachedField(S, L);
	std::complex<double> *work = NULL;
	std::complex<double> *Eta = NULL;
	size_t *ipiv = NULL;
	if(NULL == P){
		// We need to compute the vector field

		// Make vector fields
		// Determine size of the vector field grid
		int ngrid[2] = {1,1};
		for(int i = 0; i < 2; ++i){ // choose grid size
			for(int j = 0; j < n; ++j){
				if(abs(G[2*j+i]) > ngrid[i]){ ngrid[i] = abs(G[2*j+i]); }
			}
			if(ngrid[i] < 1){ ngrid[i] = 1; }
			ngrid[i] *= S->options.resolution;
			ngrid[i] = fft_next_fast_size(ngrid[i]);
		}
		const int ng2 = ngrid[0]*ngrid[1];

		work = (std::complex<double>*)S4_malloc(
			sizeof(std::complex<double>)*(5*nn + 4*ng2) +
			sizeof(size_t)*(2*n));
		P = work;
		ipiv = (size_t*)(P + 4*nn);
		Eta = (std::complex<double>*)(ipiv + 2*n);
		std::complex<double> *Ffrom = Eta + nn; // Fourier source
		std::complex<double> *Fto = Ffrom + ng2; // Fourier dest
		std::complex<double> *par = Fto + ng2; // real space parallel vector

		// Generate the vector field
		const double ing2 = 1./(double)ng2;
		int ii[2];

		double *vfield = (double*)S4_malloc(sizeof(double)*2*ng2);
		if(0 == S->Lr[2] && 0 == S->Lr[3]){ // 1D, generate the trivial field
			double nv[2] = {-S->Lr[1], S->Lr[0]};
			double nva = hypot(nv[0],nv[1]);
			nv[0] /= nva; nv[1] /= nva;
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					vfield[2*(ii[0]+ii[1]*ngrid[0])+0] = nv[0];
					vfield[2*(ii[0]+ii[1]*ngrid[0])+1] = nv[1];
				}
			}
		}else{
			int error = 0;
			S4_VERB(1, "Generating polarization vector field of size %d x %d\n", ngrid[0], ngrid[1]);
			error = Pattern_GenerateFlowField(&L->pattern, 0, S->Lr, ngrid[0], ngrid[1], vfield);

			if(0 != error){
				S4_TRACE("< Simulation_ComputeLayerBands (failed; Pattern_GenerateFlowField returned %d) [omega=%f]\n", error, S->omega[0]);
				if(NULL != vfield){ S4_free(vfield); }
				if(NULL != work){ S4_free(work); }
				if(NULL != ivalues){ S4_free(ivalues); }
				return error;
			}
		}

		// Normalize the field to max length and transform into Jones vector
		{
			double scale = 0;
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					double a = hypot(
						vfield[2*(ii[0]+ii[1]*ngrid[0])+0],
						vfield[2*(ii[0]+ii[1]*ngrid[0])+1]);
					if(a > scale){
						scale = a;
					}
				}
			}
			scale = 1./scale;
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					vfield[2*(ii[0]+ii[1]*ngrid[0])+0] *= scale;
					vfield[2*(ii[0]+ii[1]*ngrid[0])+1] *= scale;
					vfield_to_Jones_vector(&vfield[2*(ii[0]+ii[1]*ngrid[0])], &par[2*(ii[0]+ii[1]*ngrid[0])]);
				}
			}
		}

		if(NULL != S->options.vector_field_dump_filename_prefix){
			const char *layer_name = NULL != L->name ? L->name : "";
			const size_t prefix_len = strlen(S->options.vector_field_dump_filename_prefix);
			char *filename = (char*)malloc(sizeof(char) * (prefix_len + strlen(layer_name) + 1));
			strcpy(filename, S->options.vector_field_dump_filename_prefix);
			strcpy(filename+prefix_len, layer_name);
			FILE *fp = fopen(filename, "wb");
			if(NULL != fp){
				for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
					for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
						fprintf(fp, "%d\t%d\t%f\t%f\t%f\t%f\n", ii[0], ii[1], par[2*(ii[0]+ii[1]*ngrid[0])+0].real(), par[2*(ii[0]+ii[1]*ngrid[0])+0].imag(), par[2*(ii[0]+ii[1]*ngrid[0])+1].real(), par[2*(ii[0]+ii[1]*ngrid[0])+1].imag());
					} fprintf(fp, "\n");
				}
				fclose(fp);
			}
			free(filename);
		}

		/*
		for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
			for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
				std::cout << ii[0] << "\t" << ii[1] << "\t"
					<< par[2*(ii[0]+ii[1]*ngrid[0])+0].real() << "\t"
					<< par[2*(ii[0]+ii[1]*ngrid[0])+0].imag() << "\t"
					<< par[2*(ii[0]+ii[1]*ngrid[0])+1].real() << "\t"
					<< par[2*(ii[0]+ii[1]*ngrid[0])+1].imag() << std::endl;
			} std::cout << std::endl;
		}
		*/

		fft_plan plan = fft_plan_dft_2d(ngrid, Ffrom, Fto, 1);

		// We fill in the quarter blocks of F in Fortran order

		// 00 block: uy
		{
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				const int si1 = ii[1] >= ngrid[1]/2 ? ii[1]-ngrid[1]/2 : ii[1]+ngrid[1]/2;
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					const int si0 = ii[0] >= ngrid[0]/2 ? ii[0]-ngrid[0]/2 : ii[0]+ngrid[0]/2;
					Ffrom[si1+si0*ngrid[1]] = par[2*(ii[0]+ii[1]*ngrid[0])+1];
				}
			}
			fft_plan_exec(plan);
			for(int j = 0; j < n; ++j){
				for(int i = 0; i < n; ++i){
					int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
					if(f[0] < 0){ f[0] += ngrid[0]; }
					if(f[1] < 0){ f[1] += ngrid[1]; }
					P[0+i+(0+j)*n2] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
				}
			}
		}
		// 10 block: -ux
		{
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				const int si1 = ii[1] >= ngrid[1]/2 ? ii[1]-ngrid[1]/2 : ii[1]+ngrid[1]/2;
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					const int si0 = ii[0] >= ngrid[0]/2 ? ii[0]-ngrid[0]/2 : ii[0]+ngrid[0]/2;
					Ffrom[si1+si0*ngrid[1]] = -par[2*(ii[0]+ii[1]*ngrid[0])+0];
				}
			}
			fft_plan_exec(plan);
			for(int j = 0; j < n; ++j){
				for(int i = 0; i < n; ++i){
					int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
					if(f[0] < 0){ f[0] += ngrid[0]; }
					if(f[1] < 0){ f[1] += ngrid[1]; }
					P[n+i+(0+j)*n2] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
				}
			}
		}
		// 01 block: ux*
		{
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				const int si1 = ii[1] >= ngrid[1]/2 ? ii[1]-ngrid[1]/2 : ii[1]+ngrid[1]/2;
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					const int si0 = ii[0] >= ngrid[0]/2 ? ii[0]-ngrid[0]/2 : ii[0]+ngrid[0]/2;
					Ffrom[si1+si0*ngrid[1]] = std::conj(par[2*(ii[0]+ii[1]*ngrid[0])+0]);
				}
			}
			fft_plan_exec(plan);
			for(int j = 0; j < n; ++j){
				for(int i = 0; i < n; ++i){
					int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
					if(f[0] < 0){ f[0] += ngrid[0]; }
					if(f[1] < 0){ f[1] += ngrid[1]; }
					P[0+i+(n+j)*n2] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
				}
			}
		}
		// 11 block: uy*
		{
			for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
				const int si1 = ii[1] >= ngrid[1]/2 ? ii[1]-ngrid[1]/2 : ii[1]+ngrid[1]/2;
				for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
					const int si0 = ii[0] >= ngrid[0]/2 ? ii[0]-ngrid[0]/2 : ii[0]+ngrid[0]/2;
					Ffrom[si1+si0*ngrid[1]] = std::conj(par[2*(ii[0]+ii[1]*ngrid[0])+1]);
				}
			}
			fft_plan_exec(plan);
			for(int j = 0; j < n; ++j){
				for(int i = 0; i < n; ++i){
					int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
					if(f[0] < 0){ f[0] += ngrid[0]; }
					if(f[1] < 0){ f[1] += ngrid[1]; }
					P[n+i+(n+j)*n2] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
				}
			}
		}

		fft_plan_destroy(plan);
		if(NULL != vfield){ S4_free(vfield); }

		// do LU decomposition on P
		RNP::TLASupport::LUDecomposition(n2,n2, P,n2, ipiv);

		// Add to cache (assume that ipiv is immediately after P
		Simulation_AddFieldToCache((S4_Simulation*)S, L, S->n_G, P, 4*nn+2*n);
	}else{
		// P contains the cached version
		//work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*2*nn);
		work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*nn);
		Eta = work;
		ipiv = (size_t*)(P + 4*nn);
	}

	// Generate the Fourier matrix of epsilon^{-1} into first block of Epsilon2
	for(int j = 0; j < n; ++j){
		for(int i = 0; i < n; ++i){
			int dG[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
			double f[2] = {
				dG[0] * S->Lk[0] + dG[1] * S->Lk[2],
				dG[0] * S->Lk[1] + dG[1] * S->Lk[3]
				};
			double ft[2];
			Pattern_GetFourierTransform(&L->pattern, ivalues, f, ndim, unit_cell_size, ft);
			if(S->options.use_Lanczos_smoothing){
				double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
				ft[0] *= sigma;
				ft[1] *= sigma;
			}
			Eta[i+j*n] = std::complex<double>(ft[0],ft[1]);
		}
	}
	RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., &Epsilon2[n+n*n2],n2);
	RNP::LinearSolve<'N'>(n,n, Eta,n, &Epsilon2[n+n*n2],n2, NULL, NULL);

	// P stores the LU factors of F at this point.
	// We want:
	// Epsilon2 <- P L U Epsilon2 inv(U) inv(L) inv(P)
	RNP::TBLAS::MultTrM<'L','U','N','N'>(n2,n2, 1.0, P,n2, Epsilon2,n2);
	RNP::TBLAS::SolveTrM<'R','U','N','N'>(n2,n2, 1.0, P,n2, Epsilon2,n2);
	RNP::TBLAS::MultTrM<'L','L','N','U'>(n2,n2, 1.0, P,n2, Epsilon2,n2);
	RNP::TBLAS::SolveTrM<'R','L','N','U'>(n2,n2, 1.0, P,n2, Epsilon2,n2);
	RNP::TLASupport::ApplyPermutations<'L','N'>(n2,n2, Epsilon2,n2, ipiv);
	RNP::TLASupport::ApplyPermutations<'R','I'>(n2,n2, Epsilon2,n2, ipiv);

	if(NULL != work){ S4_free(work); }

	S4_free(ivalues);

	return 0;
}
