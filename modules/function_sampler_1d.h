/* Copyright (C) 2009-2013, Stanford University
 * Written by Victor Liu (vkl@stanford.edu)
 */

#ifndef FUNCTION_SAMPLER_1D_H_INCLUDED
#define FUNCTION_SAMPLER_1D_H_INCLUDED

typedef struct function_sampler_1d_ * function_sampler_1d;

typedef struct{
	double min_dy_abs, min_dy_rel;
	/* The absolute and relative y-spacing below which further refinement is
	 * not needed. min_dy_rel is relative to the maximum y-range observed.
	 */
	
	double max_curvature;
	/* Cosine of the maximum deviation angle between line segments. Angles
	 * greater than this angle will cause refinement. Note that the cosine
	 * of the angle should be specified here.
	 */
	
	double min_dx;
	/* The minimum absolute x-spacing below which further refinement is not needed.
	 * This is an absolute measure, not relative to the total x-range.
	 */

	int range_bias;
	/* Whether to bias sampling towards the top or bottom of the range.
	 * When 0, no bias is applied, and both peaks and troughs are resolved.
	 * When 1, only troughs are resolved, and when 2, only peaks are resolved.
	 */
} function_sampler_1d_options;

void function_sampler_1d_options_defaults(function_sampler_1d_options *options);
function_sampler_1d function_sampler_1d_new(
	const function_sampler_1d_options *options
);
function_sampler_1d_options *function_sampler_1d_get_options(
	const function_sampler_1d sampler
);
void function_sampler_1d_destroy(function_sampler_1d sampler);
/* Clears all added points */
void function_sampler_1d_clear(const function_sampler_1d sampler);
/* Returns nonzero if no further refinement is needed. */
int  function_sampler_1d_is_done(const function_sampler_1d sampler);
/* Returns the current number of refinement positions. */
int  function_sampler_1d_num_refine(const function_sampler_1d sampler);
/* Gets the first nx of the current refinement positions.
 * These are prioritized towards the most desirable sample positions.
 */
int  function_sampler_1d_get_refine(
	const function_sampler_1d sampler,
	int nx, double *x
);
/* Add a sample located at x with value y, with an optional id tag. */
int  function_sampler_1d_add(
	function_sampler_1d sampler,
	double x, double y, int id
);

/* Returns the current number of samples */
int  function_sampler_1d_num_samples(const function_sampler_1d sampler);
/* Gets the i-th sample. No parameter checking is performed. */
void function_sampler_1d_get(
	const function_sampler_1d sampler, int i,
	double *x, double *y, int *id
);
void function_sampler_1d_get_min(
	const function_sampler_1d sampler,
	double *x, double *y, int *id
);
void function_sampler_1d_get_max(
	const function_sampler_1d sampler,
	double *x, double *y, int *id
);

#endif /* FUNCTION_SAMPLER_1D_H_INCLUDED */
