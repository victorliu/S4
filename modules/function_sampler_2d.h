/* Copyright (C) 2009-2013, Stanford University
 * Written by Victor Liu (vkl@stanford.edu)
 */

#ifndef FUNCTION_SAMPLER_2D_H_INCLUDED
#define FUNCTION_SAMPLER_2D_H_INCLUDED

typedef struct function_sampler_2d_ * function_sampler_2d;

typedef struct{
	double min_dz_abs, min_dz_rel;
	double max_gaussian_curvature;
	double max_principal_curvature;
	double min_dxy;
	int range_bias;
} function_sampler_2d_options;

void function_sampler_2d_options_defaults(function_sampler_2d_options *options);
function_sampler_2d function_sampler_2d_new(
	const function_sampler_2d_options *options,
	int ninit, const double *xy, const double *z, const int *id
);

function_sampler_2d_options *function_sampler_2d_get_options(
	const function_sampler_2d T
);
void function_sampler_2d_destroy(function_sampler_2d sampler);
int function_sampler_2d_is_done(const function_sampler_2d sampler);
int function_sampler_2d_num_refine(const function_sampler_2d sampler);
int function_sampler_2d_get_refine(
	const function_sampler_2d sampler,
	int nxy, double *xy
);

int function_sampler_2d_size(const function_sampler_2d sampler);
int function_sampler_2d_add(
	function_sampler_2d sampler,
	const double *xy, double z, int id
);

/* Returns the current number of samples */
int  function_sampler_2d_num_samples(const function_sampler_2d sampler);
/* Gets the i-th sample. No parameter checking is performed. */
void function_sampler_2d_get(
	const function_sampler_2d sampler, int i,
	double *xy, double *z, int *id
);
void function_sampler_2d_get_min(
	const function_sampler_2d sampler,
	double *xy, double *z, int *id
);
void function_sampler_2d_get_max(
	const function_sampler_2d sampler,
	double *xy, double *z, int *id
);


#endif /* FUNCTION_SAMPLER_2D_H_INCLUDED */
