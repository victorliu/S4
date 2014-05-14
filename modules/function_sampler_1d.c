/* Copyright (C) 2009-2013, Stanford University
 * Written by Victor Liu (vkl@stanford.edu)
 */

#include "function_sampler_1d.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#ifdef DEBUG
#include <stdio.h>
#endif

/* Refinement is determined by considering consecutive triplets of samples.
 * When normalizing the three sample points to a unit square, if the angle
 * deviation between the segments exceeds the threshold, then the intervals
 * on either side of the middle sample will be subdivided.
 *
 * However, priority of the intervals is determined by a different metric.
 * Each interval is assigned a badness parameter, which is larger of the two
 * triangle areas formed by the two triplets that contain the interval.
 * Intervals with larger badness are suggested first for refinement.
 */

typedef struct sample_{
	double x, y;
	double badness;
	int id;
} sample;

struct function_sampler_1d_{
	function_sampler_1d_options opts;
	
	/* List of samples, nsamp is number of samples, and samp is a linked
	 * list of samples.
	 */
	int ns, ns_alloc;
	sample *s;
	int nbad;
	
	int ntmp_alloc;
	double *tmp;
	
	/* Min and max x and y values encountered so far */
	double x0, x1;
	double y0, y1;
};

#ifdef DEBUG
void function_sampler_1d_dump_state(const function_sampler_1d sampler, FILE *fp){
	int i;
	fprintf(fp, "{\n");
	fprintf(fp, "  samples = %d, nbad = %d\n", sampler->ns, sampler->nbad);
	for(i = 0; i < sampler->ns; ++i){
		fprintf(fp, "    %g\t%g\t%d\t%g\n", sampler->s[i].x, sampler->s[i].y, sampler->s[i].id, sampler->s[i].badness);
	}
	fprintf(fp, "}\n");
}
#endif /* DEBUG */

void function_sampler_1d_options_defaults(
	function_sampler_1d_options *opts
){
	if(NULL == opts){ return; }
	opts->min_dy_abs = 0;
	opts->min_dy_rel = 1e-3;
	opts->max_curvature = 0.1736; /* sin 10 degrees */
	opts->min_dx = 1e-6;
	opts->range_bias = 0;
}

function_sampler_1d function_sampler_1d_new(
	const function_sampler_1d_options *options
){
	function_sampler_1d sampler = (function_sampler_1d)malloc(sizeof(struct function_sampler_1d_));
	function_sampler_1d_options *opts = &(sampler->opts);
	if(NULL != options){
		opts->min_dy_rel    = options->min_dy_rel;
		opts->min_dy_abs    = options->min_dy_abs;
		opts->max_curvature = options->max_curvature;
		opts->min_dx        = options->min_dx;
		opts->range_bias    = options->range_bias;
		if(opts->min_dy_rel < 0){
			opts->min_dy_rel = 0;
		}
		if(opts->min_dy_abs < 0){
			opts->min_dy_abs = 0;
		}
		if(opts->max_curvature < 0){
			opts->max_curvature = DBL_EPSILON;
		}
		if(opts->max_curvature > 1){
			opts->max_curvature = 1;
		}
		if(opts->min_dx < DBL_EPSILON){
			opts->min_dx = DBL_EPSILON;
		}
		if(opts->range_bias < 0 || opts->range_bias > 2){
			opts->range_bias = 0;
		}
	}else{
		function_sampler_1d_options_defaults(opts);
	}
	
	sampler->x0 = DBL_MAX;
	sampler->x1 = -DBL_MAX;
	sampler->y0 = DBL_MAX;
	sampler->y1 = -DBL_MAX;
	sampler->ns = 0;
	sampler->ns_alloc = 0;
	sampler->s = NULL;
	sampler->nbad = 0;
	sampler->ntmp_alloc = 0;
	sampler->tmp = NULL;
	return sampler;
}

function_sampler_1d_options *function_sampler_1d_get_options(
	const function_sampler_1d sampler
){
	if(NULL == sampler){ return NULL; }
	return &(sampler->opts);
}

void function_sampler_1d_destroy(function_sampler_1d sampler){
	if(NULL == sampler){ return; }
	free(sampler->s);
	free(sampler->tmp);
	free(sampler);
}

void function_sampler_1d_clear(const function_sampler_1d sampler){
	if(NULL == sampler){ return; }
	sampler->ns = 0;
	sampler->nbad = 0;
	sampler->x0 = DBL_MAX;
	sampler->x1 = -DBL_MAX;
	sampler->y0 = DBL_MAX;
	sampler->y1 = -DBL_MAX;
}

int function_sampler_1d_is_done(const function_sampler_1d sampler){
	if(NULL == sampler){ return 1; }
	return (0 == sampler->nbad);
}

int function_sampler_1d_num_samples(const function_sampler_1d sampler){
	if(NULL == sampler){ return 0; }
	return sampler->ns;
}

void function_sampler_1d_get(
	function_sampler_1d sampler, int i,
	double *x, double *y, int *id
){
	*x = sampler->s[i].x;
	*y = sampler->s[i].y;
	*id = sampler->s[i].id;
}
void function_sampler_1d_get_min(
	function_sampler_1d sampler,
	double *x, double *y, int *id
){
	int i;
	if(NULL == sampler){ return; }
	if(sampler->ns <= 0){ return; }
	
	*x = sampler->s[0].x;
	*y = sampler->s[0].y;
	*id = sampler->s[0].id;
	for(i = 1; i < sampler->ns; ++i){
		if(sampler->s[i].y < *y){
			*x = sampler->s[i].x;
			*y = sampler->s[i].y;
			*id = sampler->s[i].id;
		}
	}
}
void function_sampler_1d_get_max(
	function_sampler_1d sampler,
	double *x, double *y, int *id
){
	int i;
	if(NULL == sampler){ return; }
	if(sampler->ns <= 0){ return; }
	
	*x = sampler->s[0].x;
	*y = sampler->s[0].y;
	*id = sampler->s[0].id;
	for(i = 1; i < sampler->ns; ++i){
		if(sampler->s[i].y > *y){
			*x = sampler->s[i].x;
			*y = sampler->s[i].y;
			*id = sampler->s[i].id;
		}
	}
}

int  function_sampler_1d_num_refine(const function_sampler_1d sampler){
	if(NULL == sampler){ return 0; }
	return sampler->nbad;
}

int function_sampler_1d_get_refine(
	const function_sampler_1d sampler,
	int nx, double *x
){
	int i;
	int nret = 0;
	double best = 0;
	if(NULL == sampler){ return -1; }
	if(nx < 1){ return 0; }
	if(NULL == x){ return -3; }
	
	if(1 == nx){
		double worst = 0;
		int found = 0;
		for(i = 0; i+1 < sampler->ns; ++i){
			if(sampler->s[i].badness > worst){
				found = 1;
				*x = 0.5*sampler->s[i].x + 0.5*sampler->s[i+1].x;
				worst = sampler->s[i].badness;
			}
		}
		return found;
	}

	sampler->tmp = (double*)realloc(sampler->tmp, sizeof(double) * nx);	
	
	/* Look through each interval */
	for(i = 0; i+1 < sampler->ns; ++i){
		/* no priority */
		/*
		if(sampler->s[i].badness > 0){
			double xnew = 0.5*sampler->s[i].x + 0.5*sampler->s[i+1].x;
			*x = xnew;
			++x;
			++nret;
			if(nret >= nx){ return nret; }
		}
		*/
//printf("badness[%d] = %g\n", i, sampler->s[i].badness);
		if(sampler->s[i].badness > 0 && (nret < nx || sampler->s[i].badness > best)){
			int j;
			double xnew = 0.5*sampler->s[i].x + 0.5*sampler->s[i+1].x;
//printf("xnew = %g\n", xnew);
			if(nret >= nx){ nret--; }
			
			x[nret] = xnew;
			sampler->tmp[nret] = sampler->s[i].badness;
			++nret;
			
			for(j = nret-1; j > 0; --j){
				if(sampler->tmp[j] > sampler->tmp[j-1]){
					double dt;
					dt = sampler->tmp[j];
					sampler->tmp[j] = sampler->tmp[j-1];
					sampler->tmp[j-1] = dt;
					dt = x[j]; x[j] = x[j-1]; x[j-1] = dt;
				}
			}
			best = sampler->tmp[nret-1];
		}
	}
//printf("nret = %d\n", nret);
	return nret;
}


// returns 0 if no refinement is needed
//         1 if the segment before p should be subdivided
//         2                after
//         3                both
static void update(const function_sampler_1d sampler, int is){
//printf("REFINE x=%g\n", p->x);
	double yp, y0, yn;
	double xp, x0, xn;
	xp = sampler->s[is-1].x;
	yp = sampler->s[is-1].y;
	x0 = sampler->s[is].x;
	y0 = sampler->s[is].y;
	xn = sampler->s[is+1].x;
	yn = sampler->s[is+1].y;
	
//printf("xp,x0,xn = %g, %g, %g\n", xp, x0, xn);
//printf("yp,y0,yn = %g, %g, %g\n", yp, y0, yn);
//printf("   xn-x0=%e, x0-xp=%e, min_dx=%e\n", xn-x0, x0-xp,sampler->opts.min_dx);
	if(xn-x0 < sampler->opts.min_dx && x0-xp < sampler->opts.min_dx){
		return;
	}
	if(fabs(y0-yp) < sampler->opts.min_dy_abs && fabs(yn-y0) < sampler->opts.min_dy_abs){
		return;
	}
	const double min_dy = (sampler->y1 - sampler->y0) * sampler->opts.min_dy_rel;
//printf("   yn-y0=%e, y0-yp=%e, min_dy=%e\n", yn-y0, y0-yp,min_dy);
	if(fabs(y0-yp) < min_dy && fabs(yn-y0) < min_dy){
		return;
	}
	
	double local_y_max = yp;
	if(y0 > local_y_max){ local_y_max = y0; }
	if(yn > local_y_max){ local_y_max = yn; }
	double local_y_min = yp;
	if(y0 < local_y_min){ local_y_min = y0; }
	if(yn < local_y_min){ local_y_min = yn; }
	double dx0 = (x0-xp)/(xn-xp);
	double dx1 = (xn-x0)/(xn-xp);
	double dy0 = (y0-yp)/(local_y_max-local_y_min);
	double dy1 = (yn-y0)/(local_y_max-local_y_min);
	double il0 = 1./sqrt(dx0*dx0 + dy0*dy0);
	double il1 = 1./sqrt(dx1*dx1 + dy1*dy1);
	double sinq = (dx0*dy1 - dy0*dx1) * il0 * il1;
	if(1 == sampler->opts.range_bias && sinq < 0){ return; }
	if(2 == sampler->opts.range_bias && sinq > 0){ return; }
	if(fabs(sinq) > sampler->opts.max_curvature){
		double adx0 = (x0-xp);
		double adx1 = (xn-x0);
		double ady0 = (y0-yp);
		double ady1 = (yn-y0);
		double area = fabs(adx0*ady1 - adx1*ady0);
		if(x0-xp > sampler->opts.min_dx){
			double newbad = area;
			if(sampler->s[is-1].badness <= 0){ sampler->nbad++; }
			if(newbad > sampler->s[is-1].badness){
				sampler->s[is-1].badness = newbad;
			}
		}
		if(xn-x0 > sampler->opts.min_dx){
			double newbad = area;
			if(sampler->s[is].badness <= 0){ sampler->nbad++; }
			if(newbad > sampler->s[is].badness){
				sampler->s[is].badness = newbad;
			}
		}
	}
}

int function_sampler_1d_add(
	function_sampler_1d sampler,
	double x, double y, int id
){
	int i;
	int ipos; /* position in the list where the new sample should go */

	if(NULL == sampler){ return -1; }

//function_sampler_1d_dump_state(sampler, stdout);
	
	if(x < sampler->x0){ sampler->x0 = x; }
	if(x > sampler->x1){ sampler->x1 = x; }
	if(y < sampler->y0){ sampler->y0 = y; }
	if(y > sampler->y1){ sampler->y1 = y; }
	
	/* Locate where x sits */
	for(ipos = 0; ipos < sampler->ns; ++ipos){
		if(x < sampler->s[ipos].x){
			break;
		}
	}
	
	/* At this point, we know where the new sample should be located */
	/* First add the actual sample */
	if(sampler->ns_alloc <= sampler->ns){ /* Allocated if needed */
		if(0 == sampler->ns_alloc){ sampler->ns_alloc = 64; }
		else{ sampler->ns_alloc *= 2; }
		sampler->s = (sample*)realloc(sampler->s, sizeof(sample) * sampler->ns_alloc);
	}
	{ /* Add new sample at end, swap into place */
		double badness = 0;
		for(i = ipos; i < sampler->ns; ++i){
			double dt; int it;
			dt = sampler->s[i].x ; sampler->s[i].x  = x ; x  = dt;
			dt = sampler->s[i].y ; sampler->s[i].y  = y ; y  = dt;
			it = sampler->s[i].id; sampler->s[i].id = id; id = it;
			dt = sampler->s[i].badness; sampler->s[i].badness = badness; badness = dt;
		}
		sampler->s[i].x = x;
		sampler->s[i].y = y;
		sampler->s[i].id = id;
		sampler->s[i].badness = 0;
		sampler->ns++;
	}
	
	/* Now update the badness */
	if(sampler->ns < 3){
		if(1 == sampler->ns){
			sampler->y0 = sampler->s[0].y;
			sampler->y1 = sampler->s[0].y;
		}
		return 0;
	}
	if(0 == ipos){
		update(sampler, 1);
	}else if(ipos+1 == sampler->ns){
		if(sampler->s[ipos-1].badness > 0){ sampler->nbad--; }
		sampler->s[ipos-1].badness = 0;
		update(sampler, ipos-1);
	}else{ /* internal sample */
		if(sampler->s[ipos-1].badness > 0){ sampler->nbad--; }
		sampler->s[ipos-1].badness = 0;
		if(ipos > 1){ update(sampler, ipos-1); }
		update(sampler, ipos);
		if(ipos+2 < sampler->ns){ update(sampler, ipos+1); }
	}
	
	return 0;
}
