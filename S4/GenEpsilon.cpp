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
#include "RNP/TBLAS.h"
#ifdef HAVE_BLAS
# include "RNP/TBLAS_ext.h"
#endif
#include "RNP/LinearSolve.h"
#ifdef HAVE_LAPACK
# include "RNP/LinearSolve_lapack.h"
#endif
#include "Patterning.hpp"
#include "fmm/fmm.h"

#include <limits>


int GenEpsilon(const S4_Simulation *S, const S4_Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv){
	const int n2 = 2*n;
	const int *G = S->G;
	double *ivalues = (double*)S4_malloc(sizeof(double)*(2+10)*(S->n_materials));
	double *values = ivalues + 2*(S->n_materials);

	S4_TRACE("I  Closed-form epsilon\n");

	// Get all the dielectric tensors
	bool have_tensor = false;
	for(int i = 0; i < S->n_materials; ++i){
		const S4_Material *M;
		M = &S->material[i];
		if(0 == M->type){
			std::complex<double> eps_temp(M->eps.s[0], M->eps.s[1]);
			//eps_temp = Simulation_GetEpsilonByIndex(S, L->pattern.shapes[i].tag);
			values[2*(i+1)+0] = eps_temp.real();
			values[2*(i+1)+1] = eps_temp.imag();
			eps_temp = 1./eps_temp;
			ivalues[2*(i+1)+0] = eps_temp.real();
			ivalues[2*(i+1)+1] = eps_temp.imag();
		}else{
			have_tensor = true;
		}
	}

	double mp1 = 0;
	int pwr = S->options.lanczos_smoothing_power;
	if(S->options.use_Lanczos_smoothing){
		mp1 = GetLanczosSmoothingOrder(S);
		S4_TRACE("I   Lanczos smoothing order = %f\n", mp1);
		mp1 *= S->options.lanczos_smoothing_width;
	}

	const double unit_cell_size = Simulation_GetUnitCellSize(S);

	if(!have_tensor){
		L->pat->SetTagToValueMap((std::complex<double>*)values, 1);
		// Make Epsilon in upper left block of Epsilon2
		for(int j = 0; j < n; ++j){
			for(int i = 0; i < n; ++i){
				int dG[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
				std::complex<double> ft;
				L->pat->FourierSeries(S->Lk, 1, dG, &ft);
				if(S->options.use_Lanczos_smoothing){
					double f[2] = {
						dG[0] * S->Lk[0] + dG[1] * S->Lk[2],
						dG[0] * S->Lk[1] + dG[1] * S->Lk[3]
					};
					double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
					ft *= sigma;
				}
				Epsilon2[i+j*n2] = ft;
			}
		}
		S4_TRACE("I  Epsilon(0,0) = %f,%f [omega=%f]\n", Epsilon2[0].real(), Epsilon2[0].imag(), S->omega[0]);

		if(0 == S->Lr[2] && 0 == S->Lr[3]){ // 1D proper FFF rule
			L->pat->SetTagToValueMap((std::complex<double>*)ivalues, 1);
			for(int j = 0; j < n; ++j){
				for(int i = 0; i < n; ++i){
					int dG[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
					std::complex<double> ft;
					L->pat->FourierSeries(S->Lk, 1, dG, &ft);
					if(S->options.use_Lanczos_smoothing){
						double f[2] = {
							dG[0] * S->Lk[0] + dG[1] * S->Lk[2],
							dG[0] * S->Lk[1] + dG[1] * S->Lk[3]
						};
						double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
						ft *= sigma;
					}
					Epsilon_inv[i+j*n] = ft;
				}
			}
			RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., &Epsilon2[n+n*n2],n2);
			RNP::LinearSolve<'N'>(n,n, Epsilon_inv,n, &Epsilon2[n+n*n2],n2, NULL, NULL);
			RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., Epsilon_inv,n);
			RNP::TBLAS::CopyMatrix<'A'>(n,n,&Epsilon2[0+0*n2],n2, &Epsilon2[n+0*n2],n2); // use lower block for temp storage; will be cleaned up later
			RNP::LinearSolve<'N'>(n,n, &Epsilon2[n+0*n2],n2, Epsilon_inv,n, NULL, NULL);
		}else{
			// Upper block of diagonal of Epsilon2 is already Epsilon
			RNP::TBLAS::CopyMatrix<'A'>(n,n,&Epsilon2[0+0*n2],n2, &Epsilon2[n+n*n2],n2);
			RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., Epsilon_inv,n);
			RNP::LinearSolve<'N'>(n,n, &Epsilon2[0+0*n2],n2, Epsilon_inv,n, NULL, NULL);
			RNP::TBLAS::CopyMatrix<'A'>(n,n,&Epsilon2[n+n*n2],n2, &Epsilon2[0+0*n2],n2);
		}
		RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,0., &Epsilon2[n+0*n2],n2);
		RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,0., &Epsilon2[0+n*n2],n2);
		// Epsilon2 has Epsilon's on its diagonal
	}else{ // have tensor dielectric
		const int ldv = 2*S->n_materials;
		for(int i = 0; i < S->n_materials; ++i){
			const S4_Material *M;
			M = &S->material[i];
			if(0 == M->type){
				std::complex<double> eps_temp(M->eps.s[0], M->eps.s[1]);
				values[0*ldv+2*(i+1)+0] = eps_temp.real();
				values[0*ldv+2*(i+1)+1] = eps_temp.imag();
				values[1*ldv+2*(i+1)+0] = 0;
				values[1*ldv+2*(i+1)+1] = 0;
				values[2*ldv+2*(i+1)+0] = 0;
				values[2*ldv+2*(i+1)+1] = 0;
				values[3*ldv+2*(i+1)+0] = eps_temp.real();
				values[3*ldv+2*(i+1)+1] = eps_temp.imag();
				values[4*ldv+2*(i+1)+0] = eps_temp.real();
				values[4*ldv+2*(i+1)+1] = eps_temp.imag();
			}else{
				std::complex<double> eps_temp(M->eps.s[0], M->eps.s[1]);
				// We must transpose the values array here, as well as transpose the tensor
				/*
				values[10*(i+1)+0] = M->eps.abcde[0];
				values[10*(i+1)+1] = M->eps.abcde[1];
				values[10*(i+1)+2] = M->eps.abcde[4];
				values[10*(i+1)+3] = M->eps.abcde[5];
				values[10*(i+1)+4] = M->eps.abcde[2];
				values[10*(i+1)+5] = M->eps.abcde[3];
				values[10*(i+1)+6] = M->eps.abcde[6];
				values[10*(i+1)+7] = M->eps.abcde[7];
				values[10*(i+1)+8] = M->eps.abcde[8];
				values[10*(i+1)+9] = M->eps.abcde[9];
				*/
				values[0*ldv+2*(i+1)+0] = M->eps.abcde[0];
				values[0*ldv+2*(i+1)+1] = M->eps.abcde[1];
				values[1*ldv+2*(i+1)+0] = M->eps.abcde[4];
				values[1*ldv+2*(i+1)+1] = M->eps.abcde[5];
				values[2*ldv+2*(i+1)+0] = M->eps.abcde[2];
				values[2*ldv+2*(i+1)+1] = M->eps.abcde[3];
				values[3*ldv+2*(i+1)+0] = M->eps.abcde[6];
				values[3*ldv+2*(i+1)+1] = M->eps.abcde[7];
				values[4*ldv+2*(i+1)+0] = M->eps.abcde[8];
				values[4*ldv+2*(i+1)+1] = M->eps.abcde[9];
			}
		}

		for(int k = -1; k < 4; ++k){
			L->pat->SetTagToValueMap((std::complex<double>*)&values[(k+1)*ldv], 1);
			if(-1 == k){
				for(int j = 0; j < n; ++j){
					for(int i = 0; i < n; ++i){
						int dG[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
						std::complex<double> ft;
						L->pat->FourierSeries(S->Lk, 1, dG, &ft);
						if(S->options.use_Lanczos_smoothing){
							double f[2] = {
								dG[0] * S->Lk[0] + dG[1] * S->Lk[2],
								dG[0] * S->Lk[1] + dG[1] * S->Lk[3]
							};
							double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
							ft *= sigma;
						}
						Epsilon2[i+j*n2] = ft;
					}
				}
				RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., Epsilon_inv,n);
				RNP::LinearSolve<'N'>(n,n, &Epsilon2[0+0*n2],n2, Epsilon_inv,n, NULL, NULL);
			}else{
				const int ib = k&1 ? n : 0;
				const int jb = k&2 ? n : 0;
				for(int j = 0; j < n; ++j){
					for(int i = 0; i < n; ++i){
						int dG[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
						std::complex<double> ft;
						L->pat->FourierSeries(S->Lk, 1, dG, &ft);
						if(S->options.use_Lanczos_smoothing){
							double f[2] = {
								dG[0] * S->Lk[0] + dG[1] * S->Lk[2],
								dG[0] * S->Lk[1] + dG[1] * S->Lk[3]
							};
							double sigma = GetLanczosSmoothingFactor(mp1, pwr, f);
							ft *= sigma;
						}
						Epsilon2[ib+i+(jb+j)*n2] = ft;
					}
				}
			}
		}
	}

	S4_free(ivalues);

	return 0;
}
