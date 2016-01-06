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

#include "SpectrumSampler.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct data_point_{
	double x, y;
	struct data_point_ *prev, *next;
} data_point;

struct SpectrumSampler_{
	data_point *value;
	union{
		data_point *active;
		struct{
			data_point **active;
			int n_active;
			int n_alloc;
			double *buf;
		} active_list;
	} active_set;
	int done;
	int state;
	double x0, x1;
	double y_max, y_min;
	SpectrumSampler_Options options;
};

SpectrumSampler SpectrumSampler_New(double x0, double x1, const SpectrumSampler_Options *options){
	int i;
	SpectrumSampler sampler = (SpectrumSampler)malloc(sizeof(struct SpectrumSampler_));
	if(NULL != options){
		sampler->options.initial_num_points = options->initial_num_points;
		sampler->options.range_threshold = options->range_threshold;
		sampler->options.max_bend = options->max_bend;
		sampler->options.min_dx = options->min_dx;
		sampler->options.parallelize = options->parallelize;
		if(sampler->options.initial_num_points < 3){
			sampler->options.initial_num_points = 3;
		}
		if(sampler->options.range_threshold < 0){
			sampler->options.range_threshold = 0;
		}
		if(sampler->options.max_bend < 0){
			sampler->options.max_bend = 0;
		}
		if(sampler->options.max_bend > 1){
			sampler->options.max_bend = 1;
		}
		if(sampler->options.min_dx < DBL_EPSILON){
			sampler->options.min_dx = DBL_EPSILON;
		}
	}else{
		sampler->options.initial_num_points = 33;
		sampler->options.range_threshold = 0.001;
		sampler->options.max_bend = cos(10.0*M_PI/180.0);
		sampler->options.min_dx = 1e-6;
		sampler->options.parallelize = 0;
	}

	sampler->state = 0;
	sampler->x0 = x0;
	sampler->x1 = x1;
	if(sampler->x0 > sampler->x1){
		double t = sampler->x1;
		sampler->x1 = sampler->x0;
		sampler->x0 = t;
	}
	sampler->done = 0;
	if(sampler->options.parallelize){
		const double dx = (sampler->x1 - sampler->x0) / (double)(sampler->options.initial_num_points-1);
		sampler->active_set.active_list.n_active = sampler->options.initial_num_points;
		sampler->active_set.active_list.n_alloc  = sampler->options.initial_num_points;
		sampler->active_set.active_list.active = (data_point**)malloc(sizeof(data_point*) * sampler->active_set.active_list.n_alloc);
		sampler->active_set.active_list.buf = NULL;

		sampler->value = NULL;
		data_point *tail = NULL;
		for(i = 0; i < sampler->active_set.active_list.n_active; ++i){
			// make new node
			data_point *newnode = (data_point*)malloc(sizeof(data_point));
			if(NULL == sampler->value){ sampler->value = newnode; }
			newnode->x = sampler->x0 + (double)i*dx;
			newnode->y = 0;
			newnode->prev = tail;
			newnode->next = NULL;
			if(NULL != tail){ tail->next = newnode; }
			tail = newnode;
			sampler->active_set.active_list.active[i] = tail;
		}
	}else{
		sampler->value = (data_point*)malloc(sizeof(data_point));
		sampler->active_set.active = sampler->value;
		sampler->value->prev = NULL;
		sampler->value->next = NULL;

		sampler->active_set.active->x = sampler->x0;
	}
	sampler->y_max = -DBL_MAX;
	sampler->y_min = DBL_MAX;

	sampler->options.min_dx *= (sampler->x1 - sampler->x0);

	return sampler;
}

void SpectrumSampler_Destroy(SpectrumSampler sampler){
	if(NULL == sampler){ return; }
	data_point *p = sampler->value;
	while(NULL != p){
		data_point *t = p;
		p = p->next;
		free(t);
	}
	if(sampler->options.parallelize){
		if(NULL != sampler->active_set.active_list.buf){
			free(sampler->active_set.active_list.buf);
		}
		if(NULL != sampler->active_set.active_list.active){
			free(sampler->active_set.active_list.active);
		}
	}
	free(sampler);
}

int SpectrumSampler_IsDone(const SpectrumSampler sampler){
	if(NULL == sampler){ return 1; }
	return sampler->done;
}
int SpectrumSampler_IsParallelized(const SpectrumSampler sampler){
	if(NULL == sampler){ return 0; }
	return sampler->options.parallelize;
}
double SpectrumSampler_GetFrequency(const SpectrumSampler sampler){
	if(NULL == sampler){ return 0; }
	if(sampler->done){ return 0; }
	if(sampler->options.parallelize){ return 0; }
	return sampler->active_set.active->x;
}
int SpectrumSampler_GetNumPoints(const SpectrumSampler sampler){
	if(NULL == sampler){ return 0; }
	return sampler->state;
}


SpectrumSampler_Enumerator SpectrumSampler_GetPointEnumerator(const SpectrumSampler sampler){
	return (SpectrumSampler_Enumerator)&sampler->value;
}
int SpectrumSampler_Enumerator_Get(SpectrumSampler_Enumerator e, double pt[2]){
	data_point **pp = (data_point**)e;
	data_point *p = *pp;
	if(NULL == p){ return 0; }
	pt[0] = p->x;
	pt[1] = p->y;
	*pp = p->next;
	return NULL != p->next;
}

void SpectrumSampler_GetPoints(const SpectrumSampler sampler, double *pt){
	if(NULL == sampler){ return; }
	data_point *p = sampler->value;
	while(NULL != p){
		*pt = p->x; ++pt;
		*pt = p->y; ++pt;
	}
}

// returns 0 if no refinement is needed
//         1 if the segment after  p should be subdivided
//        -1                before
static int needs_refinement(const SpectrumSampler sampler, const data_point *p){
	double yp, y0, yn;
	double xp, x0, xn;
	yp = p->prev->y; y0 = p->y; yn = p->next->y;
	xp = p->prev->x; x0 = p->x; xn = p->next->x;
//printf("   xn-x0=%e, x0-xp=%e, min_dx=%e\n", xn-x0, x0-xp,sampler->options.min_dx);
	if(xn-x0 < sampler->options.min_dx && x0-xp < sampler->options.min_dx){
		return 0;
	}
	const double min_dy = (sampler->y_max - sampler->y_min) * sampler->options.range_threshold;
//printf("   yn-y0=%e, y0-yp=%e, min_dy=%e\n", yn-y0, y0-yp,min_dy);
	if(fabs(y0-yp) < min_dy && fabs(yn-y0) < min_dy){
		return 0;
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

	double cosq = (dx0*dx1 + dy0*dy1) / sqrt((dx0*dx0 + dy0*dy0) * (dx1*dx1 + dy1*dy1));

	if(cosq < sampler->options.max_bend || dx1 > 3*dx0 || dx0 > 3*dx1){
		if(x0-xp < sampler->options.min_dx){
			return 1;
		}
		if(xn-x0 < sampler->options.min_dx){
			return -1;
		}
		if(x0-xp > xn-x0){
			return -1;
		}else{
			return 1;
		}
	}
	return 0;
}

// only for sampler->options.parallelize == 0
static void do_refinement(const SpectrumSampler sampler, data_point *p, int which){
	if(which < 0){
		p = p->prev;
	}
	data_point *newpt = (data_point*)malloc(sizeof(data_point));
	newpt->prev = p;
	newpt->next = p->next;
	newpt->x = 0.5*(p->x + p->next->x);
	p->next->prev = newpt;
	p->next = newpt;
	sampler->active_set.active = newpt;
}

int SpectrumSampler_SubmitResult(SpectrumSampler sampler, double y){
	if(NULL == sampler){ return -1; }
	if(sampler->options.parallelize){ return -1; }
	if(y > sampler->y_max){ sampler->y_max = y; }
	if(y < sampler->y_min){ sampler->y_min = y; }

	int found = 0;
	sampler->active_set.active->y = y;
	if(0 <= sampler->state && sampler->state < sampler->options.initial_num_points){
		sampler->active_set.active->next = (data_point*)malloc(sizeof(data_point));
		sampler->active_set.active->next->prev = sampler->active_set.active;
		sampler->active_set.active->next->next = NULL;
		sampler->active_set.active = sampler->active_set.active->next;
		sampler->state++;
		sampler->active_set.active->x = sampler->x0 + ((double)sampler->state / (double)(sampler->options.initial_num_points-1)) * (sampler->x1-sampler->x0);
		if(sampler->state >= sampler->options.initial_num_points){
			// Perform a scan
			sampler->active_set.active = sampler->value;
			sampler->active_set.active = sampler->active_set.active->next;
			while(NULL != sampler->active_set.active->next->next){
				int which = needs_refinement(sampler, sampler->active_set.active);
				if(0 != which){
					do_refinement(sampler, sampler->active_set.active, which);
					found = 1;
					break;
				}
				sampler->active_set.active = sampler->active_set.active->next;
			}
		}else{
			found = 1;
		}
	}else{
		if(sampler->active_set.active->prev->prev != NULL){
			sampler->active_set.active = sampler->active_set.active->prev;
		}
		while(NULL != sampler->active_set.active->next->next){
			int which = needs_refinement(sampler, sampler->active_set.active);
			if(0 != which){
				do_refinement(sampler, sampler->active_set.active, which);
				found = 1;
				break;
			}
			sampler->active_set.active = sampler->active_set.active->next;
		}
	}
	if(!found){
		sampler->done = 1;
	}
	return sampler->done;
}








int SpectrumSampler_GetFrequencies(const SpectrumSampler sampler, double **freqs){
	int i;
	if(NULL == sampler){ return 0; }
	if(NULL == freqs){ return 0; }
	if(!sampler->options.parallelize){ return 0; }

	sampler->active_set.active_list.buf = (double*)realloc(sampler->active_set.active_list.buf, sizeof(double) * (sampler->active_set.active_list.n_active));
	for(i = 0; i < sampler->active_set.active_list.n_active; ++i){
		sampler->active_set.active_list.buf[i] = sampler->active_set.active_list.active[i]->x;
	}
	*freqs = sampler->active_set.active_list.buf;
	return sampler->active_set.active_list.n_active;
}
int SpectrumSampler_GetSubmissionBuffer(const SpectrumSampler sampler, double **y){
	if(NULL == sampler){ return 0; }
	if(NULL == y){ return 0; }
	if(NULL == sampler->active_set.active_list.buf){ return 0; }
	*y = sampler->active_set.active_list.buf;
	return sampler->active_set.active_list.n_active;
}
// only for sampler->options.parallelize == 1
static data_point* do_refinement_par(data_point *p, int which){
	if(which < 0){
		p = p->prev;
	}
	data_point *newpt = (data_point*)malloc(sizeof(data_point));
	newpt->prev = p;
	newpt->next = p->next;
	newpt->x = 0.5*(p->x + p->next->x);
	p->next->prev = newpt;
	p->next = newpt;
	return newpt;
}
// pass in the same freqs pointer from above; it will get freed.
int SpectrumSampler_SubmitResults(SpectrumSampler sampler){
	int i;
	if(NULL == sampler){ return -1; }
	if(!sampler->options.parallelize){ return -1; }

	for(i = 0; i < sampler->active_set.active_list.n_active; ++i){
		double y = sampler->active_set.active_list.buf[i];
		if(y > sampler->y_max){ sampler->y_max = y; }
		if(y < sampler->y_min){ sampler->y_min = y; }
		sampler->active_set.active_list.active[i]->y = y;
	}
	sampler->state += sampler->active_set.active_list.n_active;
	sampler->active_set.active_list.n_active = 0;

	// find new active set
	data_point *t = sampler->active_set.active_list.active[0];
	if(NULL == t->prev){ // cannot start at left-most point
		t = t->next;
	}else if(NULL != t->prev->prev){ // try to move one point to the left
		t = t->prev;
	}
	while(NULL != t && NULL != t->next){
		const int which = needs_refinement(sampler, t);
		if(0 != which){
//printf("refining at %e,%e,%e  %e,%e,%e\n", t->prev->x, t->x, t->next->x, t->prev->y, t->y, t->next->y);
			if(sampler->active_set.active_list.n_active >= sampler->active_set.active_list.n_alloc){
				sampler->active_set.active_list.n_alloc *= 2;
				sampler->active_set.active_list.active = (data_point**)realloc(sampler->active_set.active_list.active, sizeof(data_point*) * sampler->active_set.active_list.n_alloc);
			}
			sampler->active_set.active_list.active[sampler->active_set.active_list.n_active] = do_refinement_par(t, which);
			++sampler->active_set.active_list.n_active;
			if(1 == which){
				t = t->next;
			}
			t = t->next;
		}
		if(NULL != t){
			t = t->next;
		}
	}

	if(0 == sampler->active_set.active_list.n_active){ sampler->done = 1; }

	return sampler->done;
}

