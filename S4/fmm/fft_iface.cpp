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

#include "fft_iface.h"
#include <cstdlib>

#ifdef HAVE_LIBFFTW3
#include <fftw3.h>
#else
#include <kiss_fft.h>
#include <tools/kiss_fftnd.h>
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
static pthread_mutex_t mutex;
#endif

int fft_next_fast_size(int n){
    while(1){
        int m=n;
        while( (m%2) == 0 ) m/=2;
        while( (m%3) == 0 ) m/=3;
        while( (m%5) == 0 ) m/=5;
        if(m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

std::complex<double> *fft_alloc_complex(size_t n){
#ifdef HAVE_LIBFFTW3
	return (std::complex<double>*)(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
#else
	return (std::complex<double>*)KISS_FFT_MALLOC(sizeof(std::complex<double>) * n);
#endif
}

void fft_free(void *p){
#ifdef HAVE_LIBFFTW3
	fftw_free(p);
#else
	KISS_FFT_FREE(p);
#endif
}

struct tag_fft_plan{
#ifdef HAVE_LIBFFTW3
	fftw_plan plan;
#else
	kiss_fftnd_cfg cfg;
	std::complex<double> *in, *out;
#endif
};

fft_plan fft_plan_dft_2d(
	int n[2],
	std::complex<double> *in, std::complex<double> *out,
	int sign
){
	fft_plan plan = NULL;
#ifdef HAVE_LIBFFTW3
# ifdef HAVE_LIBPTHREAD
	pthread_mutex_lock(&mutex);
# endif
	fftw_plan p;
	p = fftw_plan_dft(2, n, (fftw_complex*)in, (fftw_complex*)out, sign, FFTW_ESTIMATE);
# ifdef HAVE_LIBPTHREAD
	pthread_mutex_unlock(&mutex);
# endif
	if(NULL != p){
		plan = (fft_plan)malloc(sizeof(tag_fft_plan));
		plan->plan = p;
	}
#else
	kiss_fftnd_cfg cfg;
	cfg = kiss_fftnd_alloc(n, 2, sign, NULL, NULL);
	if(NULL != cfg){
		plan = (fft_plan)malloc(sizeof(tag_fft_plan));
		plan->cfg = cfg;
		plan->in = in;
		plan->out = out;
	}
#endif
	return plan;
}

void fft_plan_exec(const fft_plan plan){
#ifdef HAVE_LIBFFTW3
	fftw_execute(plan->plan);
#else
	kiss_fftnd(plan->cfg, (const kiss_fft_cpx *)plan->in, (kiss_fft_cpx *)plan->out);
#endif
}

void fft_plan_destroy(fft_plan plan){
	if(NULL == plan){ return; }
#ifdef HAVE_LIBFFTW3
# ifdef HAVE_LIBPTHREAD
	pthread_mutex_lock(&mutex);
# endif
	fftw_destroy_plan(plan->plan);
# ifdef HAVE_LIBPTHREAD
	pthread_mutex_unlock(&mutex);
# endif
#else
	free(plan->cfg);
#endif
	free(plan);
}

void fft_init(){
#ifdef HAVE_LIBPTHREAD
	if(pthread_mutex_init(&mutex, NULL)){
        printf("Unable to initialize a mutex for FFT module\n");
    }
#endif
}

void fft_destroy(){
#ifdef HAVE_LIBFFTW3
	fftw_cleanup();
#endif
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_destroy(&mutex);
#endif
}

