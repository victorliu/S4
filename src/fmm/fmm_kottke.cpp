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

#include <cstdlib>
#include <cstring>
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


// M and RMR are row-major
static void sym2x2rot(const std::complex<double> M[4], const double nvec[2], std::complex<double> RMR[4]){
	// Computes R^T M R where R = [ nvec[0] -nvec[1] ]
	//                            [ nvec[1]  nvec[0] ]
	std::complex<double> MR[4] = {
		nvec[0]*M[0] + nvec[1]*M[1], nvec[0]*M[1] - nvec[1]*M[0],
		nvec[0]*M[2] + nvec[1]*M[3], nvec[0]*M[3] - nvec[1]*M[2]
	};
	RMR[0] = nvec[0]*MR[0] + nvec[1]*MR[2]; RMR[1] = nvec[0]*MR[1] + nvec[1]*MR[3];
	RMR[2] = nvec[0]*MR[2] - nvec[1]*MR[0]; RMR[3] = nvec[0]*MR[3] - nvec[1]*MR[1];
}

int FMMGetEpsilon_Kottke(const Simulation *S, const Layer *L, const int n, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv){
	const int n2 = 2*n;
	const int *G = S->solution->G;

	// Make grid
	// Determine size of the grid
	int ngrid[2] = {1,1};
	for(int i = 0; i < 2; ++i){ // choose grid size
		for(int j = 0; j < n; ++j){
			if(abs(G[2*j+i]) > ngrid[i]){ ngrid[i] = abs(G[2*j+i]); }
		}
		if(ngrid[i] < 1){ ngrid[i] = 1; }
		ngrid[i] *= S->options.resolution;
		ngrid[i] = fft_next_fast_size(ngrid[i]);
	}
	S4_TRACE("I  Subpixel smoothing on %d x %d grid\n", ngrid[0], ngrid[1]);
	const int ng2 = ngrid[0]*ngrid[1];
	const double ing2 = 1./(double)ng2;
	// The grid needs to hold 5 matrix elements: xx,xy,yx,yy,zz
	// We actually make 5 different grids to facilitate the fft routines
	
	//std::complex<double> *work = (std::complex<double>*)S4_malloc(sizeof(std::complex<double>)*(6*ng2));
	std::complex<double> *work = fft_alloc_complex(6*ng2);
	std::complex<double>*fxx = work;
	std::complex<double>*fxy = fxx + ng2;
	std::complex<double>*fyx = fxy + ng2;
	std::complex<double>*fyy = fyx + ng2;
	std::complex<double>*fzz = fyy + ng2;
	std::complex<double>*Fto = fzz + ng2;
//memset(work, 0, sizeof(std::complex<double>) * 6*ng2);
	double *discval = (double*)S4_malloc(sizeof(double)*(L->pattern.nshapes+1));

	fft_plan plans[5];
	for(int i = 0; i <= 4; ++i){
		plans[i] = fft_plan_dft_2d(ngrid, fxx+i*ng2, Fto, 1);
	}

	int ii[2];
	for(ii[0] = 0; ii[0] < ngrid[0]; ++ii[0]){
		const int si0 = ii[0] >= ngrid[0]/2 ? ii[0]-ngrid[0]/2 : ii[0]+ngrid[0]/2;
		for(ii[1] = 0; ii[1] < ngrid[1]; ++ii[1]){
			const int si1 = ii[1] >= ngrid[1]/2 ? ii[1]-ngrid[1]/2 : ii[1]+ngrid[1]/2;
			Pattern_DiscretizeCell(&L->pattern, S->Lr, ngrid[0], ngrid[1], ii[0], ii[1], discval);
			int nnz = 0;
			int imat[2] = {-1,-1};
			for(int i = 0; i <= L->pattern.nshapes; ++i){
				if(fabs(discval[i]) > 2*std::numeric_limits<double>::epsilon()){
					if(0 == nnz){ imat[0] = i; }
					else if(1 == nnz){ imat[1] = i; }
					++nnz;
				}
			}
//S4_TRACE("I   %d,%d nnz = %d\n", ii[0], ii[1], nnz);

			if(nnz < 2){ // just one material
//fprintf(stderr, "%d\t%d\t0\t0\n", ii[0], ii[1]);
				const Material *M;
				if(0 == imat[0]){
					M = Simulation_GetMaterialByName(S, L->material, NULL);
				}else{
					M = Simulation_GetMaterialByIndex(S, L->pattern.shapes[imat[0]-1].tag);
				}
				if(0 == M->type){
					std::complex<double> eps_scalar(M->eps.s[0], M->eps.s[1]);
					fxx[si1+si0*ngrid[1]] = eps_scalar;
					fxy[si1+si0*ngrid[1]] = 0;
					fyx[si1+si0*ngrid[1]] = 0;
					fyy[si1+si0*ngrid[1]] = eps_scalar;
					fzz[si1+si0*ngrid[1]] = eps_scalar;
				}else{
					fxx[si1+si0*ngrid[1]] = std::complex<double>(M->eps.abcde[0],M->eps.abcde[1]);
					fyy[si1+si0*ngrid[1]] = std::complex<double>(M->eps.abcde[6],M->eps.abcde[7]);
					fxy[si1+si0*ngrid[1]] = std::complex<double>(M->eps.abcde[2],M->eps.abcde[3]);
					fyx[si1+si0*ngrid[1]] = std::complex<double>(M->eps.abcde[4],M->eps.abcde[5]);
					fzz[si1+si0*ngrid[1]] = std::complex<double>(M->eps.abcde[8],M->eps.abcde[9]);
				}
			}else{
				if(imat[1] > imat[0]){
					std::swap(imat[0],imat[1]);
				}
				// imat[0] is the more-contained shape
				
				double nvec[2];
				const double nxvec[2] = { ((double)ii[0]+0.5)/(double)ngrid[0]-0.5, ((double)ii[1]+0.5)/(double)ngrid[1]-0.5 };
				const double xvec[2] = { // center of current parallelogramic pixel
					S->Lr[0]*nxvec[0] + S->Lr[2]*nxvec[1],
					S->Lr[1]*nxvec[0] + S->Lr[3]*nxvec[1]
				};
				shape_get_normal(&(L->pattern.shapes[imat[0]-1]), xvec, nvec);

				if(2 == nnz && imat[1] != imat[0] && !(0 == nvec[0] && 0 == nvec[1])){ // use Kottke averaging
//fprintf(stderr, "%d\t%d\t%f\t%f\n", ii[0], ii[1], nvec[0], nvec[1]);
					const double fill = discval[imat[0]];

					// Get the two tensors
					std::complex<double> abcde[2][5];
					for(int i = 0; i < 2; ++i){
						const Material *M;
						if(0 == imat[i]){
							M = Simulation_GetMaterialByName(S, L->material, NULL);
						}else{
							M = Simulation_GetMaterialByIndex(S, L->pattern.shapes[imat[i]-1].tag);
						}
						if(0 == M->type){
							std::complex<double> eps_scalar(M->eps.s[0], M->eps.s[1]);
							abcde[i][0] = eps_scalar;
							abcde[i][1] = 0;
							abcde[i][2] = 0;
							abcde[i][3] = eps_scalar;
							abcde[i][4] = eps_scalar;
						}else{
							for(int j = 0; j < 4; ++j){
								abcde[i][j] = std::complex<double>(M->eps.abcde[2*j+0],M->eps.abcde[2*j+1]);
							}
							abcde[i][4] = std::complex<double>(M->eps.abcde[8],M->eps.abcde[9]);
						}
					}
					// The zz component is just directly obtained by averaging
					fzz[si1+si0*ngrid[1]] = discval[imat[0]] * abcde[0][4] + discval[imat[1]] * abcde[1][4];

					// Compute rotated tensors
					// abcd1 = Rot^T abcd1 Rot
					// Rot = [ nvec[0] -nvec[1] ]
					//       [ nvec[1]  nvec[0] ]
					std::complex<double> abcd[2][4];
//fprintf(stderr, " nvec = {%f,%f}, fill = %f\n", nvec[0], nvec[1], fill);
//fprintf(stderr, " abcde[0] = {%f,%f,%f,%f}\n", abcde[0][0].real(), abcde[0][1].real(), abcde[0][2].real(), abcde[0][3].real());
//fprintf(stderr, " abcde[1]  = {%f,%f,%f,%f}\n", abcde[1][0].real(), abcde[1][1].real(), abcde[1][2].real(), abcde[1][3].real());
					sym2x2rot(abcde[0], nvec, abcd[0]);
					sym2x2rot(abcde[1], nvec, abcd[1]);

					// Compute the average tau tensor into abcde[0][0-3]
					// tau(e__) = [ -1/e11         e12/e11      ]
					//            [ e21/e11   e22 - e21 e12/e11 ]
					abcde[0][0] = fill * (-1./abcd[0][0]) + (1.-fill) * (-1./abcd[1][0]);
					abcde[0][1] = fill * (abcd[0][1]/abcd[0][0]) + (1.-fill) * (abcd[1][1]/abcd[1][0]);
					abcde[0][2] = fill * (abcd[0][2]/abcd[0][0]) + (1.-fill) * (abcd[1][2]/abcd[1][0]);
					abcde[0][3] = fill * (abcd[0][3]-abcd[0][1]*abcd[0][2]/abcd[0][0]) + (1.-fill) * (abcd[1][3]-abcd[1][1]*abcd[1][2]/abcd[1][0]);
					
					// Invert the tau transform into abcd[1]
					// invtau(t__) = [ -1/t11         -t12/t11      ]
					//               [ -t21/t11   t22 - t21 t12/t11 ]
					abcd[1][0] = -1./abcde[0][0];
					abcd[1][1] = -abcde[0][1]/abcde[0][0];
					abcd[1][2] = -abcde[0][2]/abcde[0][0];
					abcd[1][3] = abcde[0][3] - abcde[0][1]*abcde[0][2]/abcde[0][0];
//fprintf(stderr, " abcd[1]  = {%f,%f,%f,%f}\n", abcd[1][0].real(), abcd[1][1].real(), abcd[1][2].real(), abcd[1][3].real());
				
					// Unrotate abcd[1] into abcd[0]
					nvec[1] = -nvec[1];
					sym2x2rot(abcd[1], nvec, abcd[0]);
//fprintf(stderr, " abcd[0] = {%f,%f,%f,%f}\n\n", abcd[0][0].real(), abcd[0][1].real(), abcd[0][2].real(), abcd[0][3].real());

					fxx[si1+si0*ngrid[1]] = abcd[0][0];
					fxy[si1+si0*ngrid[1]] = abcd[0][1];
					fyx[si1+si0*ngrid[1]] = abcd[0][2];
					fyy[si1+si0*ngrid[1]] = abcd[0][3];
				}else{ // too many, just use the area weighting
//fprintf(stderr, "%d\t%d\t3\n", ii[0], ii[1]);
					fxx[si1+si0*ngrid[1]] = 0;
					fxy[si1+si0*ngrid[1]] = 0;
					fyx[si1+si0*ngrid[1]] = 0;
					fyy[si1+si0*ngrid[1]] = 0;
					fzz[si1+si0*ngrid[1]] = 0;
					for(int i = 0; i <= L->pattern.nshapes; ++i){
						if(0 == discval[i]){ continue; }
						int j = i-1;
						const Material *M;
						if(-1 == j){
							M = Simulation_GetMaterialByName(S, L->material, NULL);
						}else{
							M = Simulation_GetMaterialByIndex(S, L->pattern.shapes[j].tag);
						}
						if(0 == M->type){
							std::complex<double> eps_scalar(M->eps.s[0], M->eps.s[1]);
							fxx[si1+si0*ngrid[0]] += discval[i]*eps_scalar;
							fyy[si1+si0*ngrid[0]] += discval[i]*eps_scalar;
							fzz[si1+si0*ngrid[0]] += discval[i]*eps_scalar;
						}else{
							std::complex<double> ea(M->eps.abcde[0],M->eps.abcde[1]);
							std::complex<double> eb(M->eps.abcde[2],M->eps.abcde[3]);
							std::complex<double> ec(M->eps.abcde[4],M->eps.abcde[5]);
							std::complex<double> ed(M->eps.abcde[6],M->eps.abcde[7]);
							fxx[si1+si0*ngrid[1]] += discval[i]*ea;
							fxy[si1+si0*ngrid[1]] += discval[i]*eb;
							fyx[si1+si0*ngrid[1]] += discval[i]*ec;
							fyy[si1+si0*ngrid[1]] += discval[i]*ed;
							fzz[si1+si0*ngrid[1]] += discval[i]*std::complex<double>(M->eps.abcde[8],M->eps.abcde[9]);
						}
					}
				}
			}
/*
fprintf(stderr, "%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ii[0], ii[1],
fxx[si1+si0*ngrid[1]].real(), fxx[si1+si0*ngrid[1]].imag(),
fxy[si1+si0*ngrid[1]].real(), fxy[si1+si0*ngrid[1]].imag(),
fyx[si1+si0*ngrid[1]].real(), fyx[si1+si0*ngrid[1]].imag(),
fyy[si1+si0*ngrid[1]].real(), fyy[si1+si0*ngrid[1]].imag(),
fzz[si1+si0*ngrid[1]].real(), fzz[si1+si0*ngrid[1]].imag()
);//*/
		}
//fprintf(stderr, "\n");
	}

	// Make Epsilon_inv first
	{
//fprintf(stderr, "fzz1[0] = %f+I %f\n", fzz[0].real(), fzz[0].imag());
//memset(Fto, 0, sizeof(std::complex<double>)*ng2);
		fft_plan_exec(plans[4]);
//fprintf(stderr, "fzz2[0] = %f+I %f\n", fzz[0].real(), fzz[0].imag());
//fprintf(stderr, "Fto[0] = %f+I %f\n", Fto[0].real(), Fto[0].imag());
		for(int j = 0; j < n; ++j){
			for(int i = 0; i < n; ++i){
				int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
				if(f[0] < 0){ f[0] += ngrid[0]; }
				if(f[1] < 0){ f[1] += ngrid[1]; }
				Epsilon2[i+j*n] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
			}
		}
	}
//fprintf(stderr, "Epsilon2[0] = %f+I %f\n", Epsilon2[0].real(), Epsilon2[0].imag());
	// Epsilon_inv needs inverting
	RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., Epsilon_inv,n);
	int solve_info;
	RNP::LinearSolve<'N'>(n,n, Epsilon2,n, Epsilon_inv,n, &solve_info, NULL);
//fprintf(stderr, "Epsilon_inv[0] = %f+I %f\n", Epsilon_inv[0].real(), Epsilon_inv[0].imag());

	// We fill in the quarter blocks of F in Fortran order
	for(int w = 0; w < 4; ++w){
		int Ecol = (w&1 ? n : 0);
		int Erow = (w&2 ? n : 0);
//memset(Fto, 0, sizeof(std::complex<double>)*ng2);
		fft_plan_exec(plans[w]);
		
//fprintf(stderr, "Fto(%d)[0] = %f+I %f\n", w, Fto[0].real(), Fto[0].imag());
		for(int j = 0; j < n; ++j){
			for(int i = 0; i < n; ++i){
				int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
				if(f[0] < 0){ f[0] += ngrid[0]; }
				if(f[1] < 0){ f[1] += ngrid[1]; }
				Epsilon2[Erow+i+(Ecol+j)*n2] = ing2 * Fto[f[1]+f[0]*ngrid[1]];
			}
		}
	}
//fprintf(stderr, "Epsilon2[0] = %f+I %f\n", Epsilon2[0].real(), Epsilon2[0].imag());
	for(int i = 0; i <= 4; ++i){
		fft_plan_destroy(plans[i]);
	}
	
	S4_free(discval);
	//S4_free(work);
	fft_free(work);
	
	return 0;
}
