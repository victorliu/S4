/* Copyright (C) 2009-2013, Stanford University
 * Written by Victor Liu (vkl@stanford.edu)
 */

#include "function_sampler_2d.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
//#ifdef DEBUG
#include <stdio.h>
//#endif

static void update_edge(const function_sampler_2d T, int h);
void DT_eval_bad(const function_sampler_2d T);

extern double orient2d(const double*, const double*, const double*);
extern double incircle(const double*, const double*, const double*, const double*);

/* Vertex data structure */
typedef struct{
	double r[3]; /* x and y coordinates */
	int h; /* index of an outgoing halfedge */
	/* If a point lies on the boundary of the convex hull, then the halfedge
	 * is guaranteed to be "rewound" in the sense that it is on the boundary.
	 */
	int id;
} Vert;

/* Halfedge data structure */
typedef struct{
	int next; /* index of the next halfedge in the ring */
	int flip; /* index of the opposite halfedge across an edge (-1 if none) */
	int from; /* index of origin vertex */
	int face; /* index of the parent face */
} Half;

typedef struct{
	int h;
	double badness;
} Face;

struct function_sampler_2d_{
	function_sampler_2d_options opts;
	
		/* vector of vertices */
	int nv, nv_alloc;
	Vert *v;
	
	/* vector of halfedges */
	int nh, nh_alloc;
	Half *h;
	
	/* vector of faces */
	int nf, nf_alloc;
	Face *f;
	/* The halfedge of a face is the one (of three) with the smallest index */

	int hlast;

	int nbad;
	int ntmp_alloc;
	double *tmp;
	
	int bad_valid;

	/* Min and max values encountered so far */
	double minxy[2], maxxy[2];
	double z0, z1;
};

#define NEXT(H) (T->h[H].next)
#define FLIP(H) (T->h[H].flip)
#define FROM(H) (T->h[H].from)
#define FACE(H) (T->h[H].face)
#define SETHALF(H,NXT,OPP,ORG,TRI) do{ \
	T->h[H].next = NXT; \
	T->h[H].flip = OPP; \
	T->h[H].from = ORG; \
	T->h[H].face = TRI; \
}while(0)

void DT_check(const function_sampler_2d T){
	int i;
	for(i = 0; i < T->nv; ++i){
		assert(0 <= T->v[i].h && T->v[i].h < T->nh);
	}
	for(i = 0; i < T->nf; ++i){
		assert(0 <= T->f[i].h && T->f[i].h < T->nh);
	}
	for(i = 0; i < T->nh; ++i){
		assert(0 <= T->h[i].next && T->h[i].next < T->nh);
		assert(-1 == T->h[i].flip || (0 <= T->h[i].flip && T->h[i].flip < T->nh));
		assert(0 <= T->h[i].from && T->h[i].from < T->nv);
		assert(0 <= T->h[i].face && T->h[i].face < T->nf);
		assert(NEXT(NEXT(NEXT(i))) == i);
		if(-1 != FLIP(i)){ assert(FLIP(FLIP(i)) == i); }
	}
}

void DT_dump(const function_sampler_2d T){
	int i;
	printf("{\nnxt : opp : from : to   : face\n");
	for(i = 0; i < T->nh; ++i){
		printf("%d\t%3d : %3d : %3d : %3d : %3d\n",
			i, NEXT(i), FLIP(i), FROM(i), FROM(NEXT(i)), FACE(i)
		);
	}
	printf("}\n");
}

void function_sampler_2d_options_defaults(function_sampler_2d_options *opts){
	opts->min_dz_abs = 0;
	opts->min_dz_rel = 0.001;
	opts->max_gaussian_curvature = DBL_MAX;
	opts->max_principal_curvature = sin(20*3.141592654/180.);
	opts->min_dxy = 1e-6;
	opts->range_bias = 0;
}

function_sampler_2d function_sampler_2d_new(
	const function_sampler_2d_options *options,
	int ninit, const double *xy, const double *z, const int *id
){
	int i;
	function_sampler_2d T;
	function_sampler_2d_options *opts;
	if(ninit < 3){ return NULL; }
	T = (function_sampler_2d)malloc(sizeof(struct function_sampler_2d_));
	if(NULL == T){ return NULL; }
	opts = &(T->opts);
	if(NULL != options){
		opts->min_dz_abs = options->min_dz_abs;
		opts->min_dz_rel = options->min_dz_rel;
		opts->max_principal_curvature = options->max_principal_curvature;
		opts->max_gaussian_curvature  = options->max_gaussian_curvature;
		opts->min_dxy         = options->min_dxy;
		if(opts->min_dz_abs < 0){
			opts->min_dz_abs = 0;
		}
		if(opts->min_dz_rel < 0){
			opts->min_dz_rel = 0;
		}
		if(opts->max_principal_curvature < 0){
			opts->max_principal_curvature = DBL_EPSILON;
		}
		if(opts->max_principal_curvature > 1){
			opts->max_principal_curvature = 1;
		}
		if(opts->min_dxy < DBL_EPSILON){
			opts->min_dxy = DBL_EPSILON;
		}
	}else{
		function_sampler_2d_options_defaults(opts);
	}
	/* Create a new initial DT with the following layout:
	 *                c
	 *                +
	 *              .'|
	 *            .'  |
	 *       e2 .'    |
	 *        .'h2  h1|e1
	 *      .'    f0  |
	 *    .'    h0    |
	 *   +------------+
	 *  a      e0      b
	 */
	/* Initialize vertex list */
	T->nv = 3;
	T->nv_alloc = 4;
	T->v = (Vert*)malloc(sizeof(Vert) * T->nv_alloc);
	/* Initialize halfedge list */
	T->nh = 3;
	T->nh_alloc = 4;
	T->h = (Half*)malloc(sizeof(Half) * T->nh_alloc);
	/* Initialize face list */
	T->nf = 1;
	T->nf_alloc = 4;
	T->f = (Face*)malloc(sizeof(Face) * T->nf_alloc);
	/* Add the vertices */
	T->v[0].r[0] = xy[0];
	T->v[0].r[1] = xy[1];
	T->v[0].r[2] = z[0];
	T->v[0].h = 0;
	T->v[0].id = (NULL != id ? id[0] : 0);
	T->minxy[0] = xy[0];
	T->minxy[1] = xy[1];
	T->maxxy[0] = xy[0];
	T->maxxy[1] = xy[1];
	
	T->v[1].r[0] = xy[2];
	T->v[1].r[1] = xy[3];
	T->v[1].r[2] = z[1];
	T->v[1].h = 1;
	T->v[1].id = (NULL != id ? id[1] : 0);
	if(xy[2] < T->minxy[0]){ T->minxy[0] = xy[2]; }
	if(xy[3] < T->minxy[1]){ T->minxy[1] = xy[3]; }
	if(xy[2] > T->maxxy[0]){ T->maxxy[0] = xy[2]; }
	if(xy[3] > T->maxxy[1]){ T->maxxy[1] = xy[3]; }
	
	T->v[2].r[0] = xy[4];
	T->v[2].r[1] = xy[5];
	T->v[2].r[2] = z[2];
	T->v[2].h = 2;
	T->v[2].id = (NULL != id ? id[2] : 0);
	if(xy[4] < T->minxy[0]){ T->minxy[0] = xy[4]; }
	if(xy[5] < T->minxy[1]){ T->minxy[1] = xy[5]; }
	if(xy[4] > T->maxxy[0]){ T->maxxy[0] = xy[4]; }
	if(xy[5] > T->maxxy[1]){ T->maxxy[1] = xy[5]; }

	/* Swap vertices if negative orientation */
	if(orient2d(T->v[0].r, T->v[1].r, T->v[2].r) < 0){
		T->v[1].r[0] = xy[4];
		T->v[1].r[1] = xy[5];
		T->v[1].r[2] = z[2];
		T->v[1].id = id[2];
		T->v[2].r[0] = xy[2];
		T->v[2].r[1] = xy[3];
		T->v[2].r[2] = z[1];
		T->v[2].id = id[1];
	}

	/* Create the new halfedges */
	SETHALF(0, 1, -1, 0, 0);
	SETHALF(1, 2, -1, 1, 0);
	SETHALF(2, 0, -1, 2, 0);
	T->f[0].h = 0;
	T->f[0].badness = 0;
	
	T->hlast = 0;
	
	T->bad_valid = 0;
	
	T->nbad = 0;
	T->ntmp_alloc = 0;
	T->tmp = NULL;
	T->z0 = z[0];
	T->z1 = z[0];
	
	for(i = 3; i < ninit; ++i){
		function_sampler_2d_add(T, &xy[2*i], z[i], (NULL != id ? id[i] : 0));
	}
	
	return T;
}

function_sampler_2d_options *function_sampler_2d_get_options(
	const function_sampler_2d T
){
	if(NULL == T){ return NULL; }
	return &(T->opts);
}

void function_sampler_2d_destroy(function_sampler_2d T){
	if(NULL == T){ return; }
	free(T->h);
	free(T->f);
	free(T->v);
	free(T);
}

int function_sampler_2d_is_done(const function_sampler_2d T){
	if(NULL == T){ return 1; }
	if(!T->bad_valid){ DT_eval_bad(T); }
	return (0 == T->nbad);
}

int function_sampler_2d_num_samples(const function_sampler_2d T){
	if(NULL == T){ return 0; }
	return T->nv;
}

void function_sampler_2d_get(
	function_sampler_2d T, int i,
	double *xy, double *z, int *id
){
	xy[0] = T->v[i].r[0];
	xy[1] = T->v[i].r[1];
	*z    = T->v[i].r[2];
	*id   = T->v[i].id;
}
void function_sampler_2d_get_min(
	function_sampler_2d T,
	double *xy, double *z, int *id
){
	int i;
	if(NULL == T){ return; }
	if(T->nv <= 0){ return; }
	
	xy[0] = T->v[0].r[0];
	xy[1] = T->v[0].r[1];
	*z    = T->v[0].r[2];
	*id   = T->v[0].id;
	for(i = 1; i < T->nv; ++i){
		if(T->v[i].r[2] < *z){
			xy[0] = T->v[i].r[0];
			xy[1] = T->v[i].r[1];
			*z    = T->v[i].r[2];
			*id   = T->v[i].id;
		}
	}
}
void function_sampler_2d_get_max(
	function_sampler_2d T,
	double *xy, double *z, int *id
){
	int i;
	if(NULL == T){ return; }
	if(T->nv <= 0){ return; }
	
	xy[0] = T->v[0].r[0];
	xy[1] = T->v[0].r[1];
	*z    = T->v[0].r[2];
	*id   = T->v[0].id;
	for(i = 1; i < T->nv; ++i){
		if(T->v[i].r[2] > *z){
			xy[0] = T->v[i].r[0];
			xy[1] = T->v[i].r[1];
			*z    = T->v[i].r[2];
			*id   = T->v[i].id;
		}
	}
}

void DT_eval_bad(const function_sampler_2d T){
	int i;
	for(i = 0; i < T->nf; ++i){
		T->f[i].badness = 0;
	}
	T->nbad = 0;
	for(i = 0; i < T->nh; ++i){
		if(i < FLIP(i)){ continue; }
		update_edge(T, i);
	}
	T->bad_valid = 1;
}

int function_sampler_2d_size(const function_sampler_2d T){
	if(NULL == T){ return 0; }
	return T->nv;
}
int function_sampler_2d_num_refine(const function_sampler_2d T){
	if(NULL == T){ return 0; }
	if(!T->bad_valid){ DT_eval_bad(T); }
	return T->nbad;
}

void geom_circum_tri2d(
        const double a[2],
        const double b[2],
        const double c[2],
        double circumcenter[2],
        double *xi,
        double *eta
){
  double xba, yba, xca, yca;
  double balength, calength;
  double denominator;
  double xcirca, ycirca;

  /* Use coordinates relative to point `a' of the triangle. */
  xba = b[0] - a[0];
  yba = b[1] - a[1];
  xca = c[0] - a[0];
  yca = c[1] - a[1];
  /* Squares of lengths of the edges incident to `a'. */
  balength = xba * xba + yba * yba;
  calength = xca * xca + yca * yca;

  /* Calculate the denominator of the formulae. */
#ifdef EXACT
  /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html */
  /* to ensure a correctly signed (and reasonably accurate) result, */
  /* avoiding any possibility of division by zero. */
  denominator = 0.5 / orient2d(b, c, a);
#else
  /* Take your chances with floating-point roundoff. */
  denominator = 0.5 / (xba * yca - yba * xca);
#endif

  /* Calculate offset (from `a') of circumcenter. */
  xcirca = (yca * balength - yba * calength) * denominator;
  ycirca = (xba * calength - xca * balength) * denominator;
  circumcenter[0] = xcirca;
  circumcenter[1] = ycirca;

  if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a */
    /* coordinate system with a xi-axis directed from `a' to `b' and */
    /* an eta-axis directed from `a' to `c'. The values for xi and eta */
    /* are computed by Cramer's Rule for solving systems of linear */
    /* equations. */
    *xi = (xcirca * yca - ycirca * xca) * (2.0 * denominator);
    *eta = (ycirca * xba - xcirca * yba) * (2.0 * denominator);
  }
}

static inline void get_center(const double *a, const double *b, const double *c, double *cc){
	static const double third = (1./3.);
	cc[0] = (a[0] + b[0] + c[0]) * third;
	cc[1] = (a[1] + b[1] + c[1]) * third;
/*
	double avg[2] = {
		(a[0] + b[0] + c[0]) / 3.,
		(a[1] + b[1] + c[1]) / 3.
	};
	double xi, eta;
	geom_circum_tri2d(a, b, c, cc, &xi, &eta);
	if(xi <= 0 || eta <= 0 || xi+eta >= 1){
		cc[0] = avg[0];
		cc[1] = avg[1];
	}else{
		cc[0] += a[0];
		cc[1] += a[1];
	}
*/
}

int function_sampler_2d_get_refine(
	const function_sampler_2d T,
	int nxy, double *xy
){
	int i;
	int nret = 0;
	double best = 0;
	if(NULL == T){ return -1; }
	if(nxy < 1){ return 0; }
	if(NULL == xy){ return -3; }
	
	if(!T->bad_valid){ DT_eval_bad(T); }
	
	if(1 == nxy){
		double worst = 0;
		int found = 0;
		for(i = 0; i < T->nf; ++i){
			if(T->f[i].badness > worst){
				const int h = T->f[i].h;
				const int ia = FROM(h);
				const int ib = FROM(NEXT(h));
				const int ic = FROM(NEXT(NEXT(h)));
				found = 1;
				get_center(T->v[ia].r, T->v[ib].r, T->v[ic].r, xy);
				worst = T->f[i].badness;
			}
		}
		return found;
	}
	
	T->tmp = (double*)realloc(T->tmp, sizeof(double) * nxy);
	
	/* Look through each interval */
	for(i = 0; i+1 < T->nf; ++i){
		/* no priority */
		/*
		if(T->f[i].badness > 0){
			const int h = T->f[i].h;
			const int ia = FROM(h);
			const int ib = FROM(NEXT(h));
			const int ic = FROM(NEXT(NEXT(h)));
			xy[0] = (T->v[ia].r[0] + T->v[ib].r[0] + T->v[ic].r[0]) / 3.;
			xy[1] = (T->v[ia].r[1] + T->v[ib].r[1] + T->v[ic].r[1]) / 3.;
			xy += 2;
			++nret;
			if(nret >= nx){ return nret; }
		}
		*/
		if(T->f[i].badness > 0 && (nret < nxy || T->f[i].badness > best)){
			int j;
			const int h = T->f[i].h;
			const int ia = FROM(h);
			const int ib = FROM(NEXT(h));
			const int ic = FROM(NEXT(NEXT(h)));
			/*
			const double xynew[2] = {
				(T->v[ia].r[0] + T->v[ib].r[0] + T->v[ic].r[0]) / 3.,
				(T->v[ia].r[1] + T->v[ib].r[1] + T->v[ic].r[1]) / 3.
			};*/
			double xynew[2];
			get_center(T->v[ia].r, T->v[ib].r, T->v[ic].r, xynew);
			if(nret >= nxy){ nret--; }
			
			xy[2*nret+0] = xynew[0];
			xy[2*nret+1] = xynew[1];
			T->tmp[nret] = T->f[i].badness;
			++nret;
			
			for(j = nret-1; j > 0; --j){
				if(T->tmp[j] > T->tmp[j-1]){
					double dt;
					dt = T->tmp[j];
					T->tmp[j] = T->tmp[j-1];
					T->tmp[j-1] = dt;
					dt = xy[2*j+0]; xy[2*j+0] = xy[2*(j-1)+0]; xy[2*(j-1)+0] = dt;
					dt = xy[2*j+1]; xy[2*j+1] = xy[2*(j-1)+1]; xy[2*(j-1)+1] = dt;
				}
			}
			best = T->tmp[nret-1];
		}
	}
	return nret;
}

static double geom_norm3d(const double v[3]){
        double a[3] = {fabs(v[0]),fabs(v[1]),fabs(v[2])};
        double w = a[0];
        if(a[1] > w){ w = a[1]; }
        if(a[2] > w){ w = a[2]; }
        if(0 == w){
                return a[0] + a[1] + a[2];
        }else{
                a[0] /= w; a[1] /= w; a[2] /= w;
                w *= sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
                return w;
        }
}
static double geom_normalize3d(double v[3]){
	double n = geom_norm3d(v);
	v[0] /= n; v[1] /= n; v[2] /= n;
	return n;
}
static double geom_dot3d(const double a[3], const double b[3]){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static double geom_cross2d(const double a[2], const double b[2]){
	return a[0]*b[1]-a[1]*b[0];
}
static void geom_cross3d(const double a[3], const double b[3], double result[3]){
	result[0] = a[1]*b[2] - a[2]*b[1];
	result[1] = a[2]*b[0] - a[0]*b[2];
	result[2] = a[0]*b[1] - a[1]*b[0];
}
// returns 0 if no refinement is needed
//         1 if the segment before p should be subdivided
//         2                after
//         3                both
static void update_edge(const function_sampler_2d T, int h){
	if(-1 == FLIP(h)){ return; }
	//assert(0 <= h && h < T->nh);
	//assert(0 <= FLIP(h) && FLIP(h) < T->nh);
	
//printf("update_edge on h=%d, from %d to %d\n", h, FROM(h), FROM(NEXT(h)));

	const int f0 = FACE(h);
	const int f1 = FACE(FLIP(h));
//printf("f0,f1 = %d, %d, nf = %d\n", f0, f1, T->nf);
	assert(0 <= f0 && f0 < T->nf);
	assert(0 <= f1 && f1 < T->nf);
	const int ia = FROM(h);
	const int ib = FROM(NEXT(h));
	const int ic = FROM(NEXT(NEXT(h)));
	const int id = FROM(NEXT(NEXT(FLIP(h))));
	
//printf("id = %d/%d, (%g, %g, %g)\n", id, T->nv, T->v[id].r[0], T->v[id].r[1], T->v[id].r[2]);
	assert(0 <= ia && ia < T->nv);
	assert(0 <= ib && ib < T->nv);
	assert(0 <= ic && ic < T->nv);
	assert(0 <= id && id < T->nv);
	const double ab[3] = {
		T->v[ib].r[0] - T->v[ia].r[0],
		T->v[ib].r[1] - T->v[ia].r[1],
		T->v[ib].r[2] - T->v[ia].r[2]
	};
//printf(" ab = %g, %g, %g\n", ab[0], ab[1], ab[2]);
	const double ab2 = hypot(ab[0],ab[1]);
//printf(" ab2 = %g\n", ab2);
	const double cd[2] = {
		T->v[id].r[0] - T->v[ic].r[0],
		T->v[id].r[1] - T->v[ic].r[1]
	};
//printf(" cd = %g, %g\n", cd[0], cd[1]);
	const double nrm2 = ab2 / fabs(geom_cross2d(cd, ab));
//printf(" nrm2 = %g\n", nrm2);
	double ac[3] = {
		T->v[ic].r[0] - T->v[ia].r[0],
		T->v[ic].r[1] - T->v[ia].r[1],
		T->v[ic].r[2] - T->v[ia].r[2]
	};
	double ad[3] = {
		T->v[id].r[0] - T->v[ia].r[0],
		T->v[id].r[1] - T->v[ia].r[1],
		T->v[id].r[2] - T->v[ia].r[2]
	};
	double n0[3], n1[3], n01[3], sinq;
	double zrange[2] = { T->v[ia].r[2], T->v[ia].r[2] };
	double minxy[2] = { T->v[ia].r[0], T->v[ia].r[1] };
	double maxxy[2] = { T->v[ia].r[0], T->v[ia].r[1] };
	//double zab = 0.5*(T->v[ia].r[2] + T->v[ib].r[2]);
	//double zcd = 0.5*(T->v[ic].r[2] + T->v[id].r[2]);
	const double d[2] = {
		geom_cross2d(ab, ac) / ab2,
		geom_cross2d(ad, ab) / ab2
	};
//printf(" d = %g,%g\n", d[0], d[1]);
	if(d[0] < T->opts.min_dxy && d[1] < T->opts.min_dxy){
//printf("  d's too small, returning\n");
		return;
	}
	if(T->v[ib].r[0] < minxy[0]){ minxy[0] = T->v[ib].r[0]; }
	if(T->v[ib].r[0] > maxxy[0]){ maxxy[0] = T->v[ib].r[0]; }
	if(T->v[ib].r[1] < minxy[1]){ minxy[1] = T->v[ib].r[1]; }
	if(T->v[ib].r[1] > maxxy[1]){ maxxy[1] = T->v[ib].r[1]; }
	if(T->v[ib].r[2] < zrange[0]){ zrange[0] = T->v[ib].r[2]; }
	if(T->v[ib].r[2] > zrange[1]){ zrange[1] = T->v[ib].r[2]; }
	if(T->v[ic].r[0] < minxy[0]){ minxy[0] = T->v[ic].r[0]; }
	if(T->v[ic].r[0] > maxxy[0]){ maxxy[0] = T->v[ic].r[0]; }
	if(T->v[ic].r[1] < minxy[1]){ minxy[1] = T->v[ic].r[1]; }
	if(T->v[ic].r[1] > maxxy[1]){ maxxy[1] = T->v[ic].r[1]; }
	if(T->v[ic].r[2] < zrange[0]){ zrange[0] = T->v[ic].r[2]; }
	if(T->v[ic].r[2] > zrange[1]){ zrange[1] = T->v[ic].r[2]; }
	if(T->v[id].r[0] < minxy[0]){ minxy[0] = T->v[id].r[0]; }
	if(T->v[id].r[0] > maxxy[0]){ maxxy[0] = T->v[id].r[0]; }
	if(T->v[id].r[1] < minxy[1]){ minxy[1] = T->v[id].r[1]; }
	if(T->v[id].r[1] > maxxy[1]){ maxxy[1] = T->v[id].r[1]; }
	if(T->v[id].r[2] < zrange[0]){ zrange[0] = T->v[id].r[2]; }
	if(T->v[id].r[2] > zrange[1]){ zrange[1] = T->v[id].r[2]; }
	if(maxxy[0] - minxy[0] < T->opts.min_dxy && maxxy[1] - minxy[1] < T->opts.min_dxy){
		return;
	}
	if(zrange[1] - zrange[0] < T->opts.min_dz_abs){ return; }
	if(zrange[1] - zrange[0] < (T->z1 - T->z0) * T->opts.min_dz_rel){
//printf("  z range too small, returning\n");
		return;
	}
//printf(" ac = %g, %g, %g\n", ac[0], ac[1], ac[2]);
//printf(" ad = %g, %g, %g\n", ad[0], ad[1], ad[2]);
//printf(" nrm2 = %g\n", nrm2);
	//ac[0] *= nrm2; ac[1] *= nrm2;
	//ad[0] *= nrm2; ad[1] *= nrm2;
//printf(" ab = %g, %g, %g\n", ab[0], ab[1], ab[2]);
//printf(" ac = %g, %g, %g\n", ac[0], ac[1], ac[2]);
//printf(" ad = %g, %g, %g\n", ad[0], ad[1], ad[2]);
	geom_cross3d(ab, ac, n0);
	geom_cross3d(ad, ab, n1);
//printf(" n0 = %g, %g, %g\n", n0[0], n0[1], n0[2]);
//printf(" n1 = %g, %g, %g\n", n1[0], n1[1], n1[2]);
	geom_normalize3d(n0);
	geom_normalize3d(n1);
	geom_cross3d(n0, n1, n01);
	sinq = geom_norm3d(n01);
//printf("  sinq = %g\n", sinq);
	if(sinq > T->opts.max_principal_curvature){
		double t[3], vol;
		geom_cross3d(ac, ab, t);
		vol = fabs(geom_dot3d(ad, t));
		if(d[0] > T->opts.min_dxy){
			if(T->f[f0].badness <= 0){ T->nbad++; }
			if(vol > T->f[f0].badness){
				T->f[f0].badness = vol;
//printf("     marking badness = %g, nbad = %d\n", vol, T->nbad);
			}
		}
		if(d[1] > T->opts.min_dxy){
			if(T->f[f1].badness <= 0){ T->nbad++; }
			if(vol > T->f[f1].badness){
				T->f[f1].badness = vol;
//printf("     marking badness = %g, nbad = %d\n", vol, T->nbad);
			}
		}
	}
}

/* Returns a halfedge of a face in ih (smallest index ih)
 * Returns 0 if inside a triangle, 1 if on exterior
 */
int DT_locate(const function_sampler_2d T, const double r[2], int *ih){
	while(1){
		const int h0 = *ih;
		const int h1 = NEXT(h0);
		const int h2 = NEXT(h1);
		int nih[3];
		int nf = 0;
		if(orient2d(T->v[FROM(h0)].r, T->v[FROM(h1)].r, r) < 0){
			nih[nf] = FLIP(*ih);
			if(-1 == nih[nf]){ return 1; }
			nf++;
		}
		if(orient2d(T->v[FROM(h1)].r, T->v[FROM(h2)].r, r) < 0){
			nih[nf] = FLIP(h1);
			if(-1 == nih[nf]){ *ih = h1; return 1; }
			nf++;
		}
		if(orient2d(T->v[FROM(h2)].r, T->v[FROM(h0)].r, r) < 0){
			nih[nf] = FLIP(h2);
			if(-1 == nih[nf]){ *ih = h2; return 1; }
			nf++;
		}
		if(0 == nf){
			if(h1 < *ih){ *ih = h1; }
			if(h2 < *ih){ *ih = h2; }
			return 0;
		}else if(1 == nf){
			*ih = nih[0];
		}else{
			*ih = nih[rand()%nf];
		}
	}
}

void DT_flip(function_sampler_2d T, int h){
	/* assert(0 <= h && h < T->nh && -1 != FLIP(h)); */
	const int h0 = h;
	const int h1 = FLIP(h0), h2 = NEXT(h0), h3 = NEXT(h2);
	const int h4 = NEXT(h1), h5 = NEXT(h4);
	const int ia = FROM(h0), ib = FROM(h2);
	const int ic = FROM(h3), id = FROM(h5);
	const int f0 = FACE(h0), f1 = FACE(h1);
	/* Flips an edge in its parent quadrilateral
	 *
	 *     c                  b          c                  b
	 *      +----------------+            +----------------+
	 *      |       h2     .'|            |`.     h2       |
	 *      |   f0       .'  |            |  `.      f1    |
	 *      |          .'    |            |    `.          |
	 *      |h3    h0.'      |            |h3    `.h1    h5|
	 *      |      .' h1   h5|    ===>    |      h0`.      |
	 *      |    .'          |            |          `.    |
	 *      |  .'      f1    |            |   f0       `.  |
	 *      |.'     h4       |            |       h4     `.|
	 *      +----------------+            +----------------+
	 *     a                  d          a                  d
	 */
	NEXT(h0) = h3; FROM(h0) = id;
	NEXT(h1) = h5; FROM(h1) = ic;
	NEXT(h2) = h1; FACE(h2) = f1;
	NEXT(h3) = h4;
	NEXT(h4) = h0; FACE(h4) = f0;
	NEXT(h5) = h2;
	T->v[ia].h = h4; T->v[ib].h = h2;
	T->v[ic].h = h3; T->v[id].h = h5;
	
	/* Re-set the face pointers */
	T->f[f0].h = h0;
	if(h3 < h0){ T->f[f0].h = h3; }
	if(h4 < T->f[f0].h){ T->f[f0].h = h4; }
	
	T->f[f1].h = h1;
	if(h2 < h1){ T->f[f1].h = h2; }
	if(h5 < T->f[f1].h){ T->f[f1].h = h5; }
}

/* h should be a halfedge emanating from the newly added vertex.
 * We assume that h is "rewound" so that if the new vertex is on
 * the boundary, h is on the boundary.
 */
void DT_fixup(function_sampler_2d T, int h){
	const int h0 = h;
	const int ip = FROM(h);
	h = NEXT(h);
	do{
		const int ht = FLIP(h);
		if(-1 != ht){
			const int ia = FROM(ht);
			const int ib = FROM(NEXT(ht));
			const int ic = FROM(NEXT(NEXT(ht)));
			if(/*orient2d(T->v[ia].r, T->v[ic].r, T->v[ip].r) < 0 &&*/
				incircle(T->v[ia].r, T->v[ib].r, T->v[ic].r, T->v[ip].r) > 0
			){
				DT_flip(T, h);
				h = NEXT(NEXT(h));
			}else{ /* advance */
				h = FLIP(NEXT(h));
				if(-1 == h || h == h0){ break; }
				h = NEXT(h);
			}
		}else{
			h = FLIP(NEXT(h));
			if(-1 == h || h == h0){ break; }
			h = NEXT(h);
		}
	}while(1);
}

/* Allocate more halfedges */
static void DT_addhalfs(function_sampler_2d T, int n){ /* Only performs allocation */
	if(T->nh+n >= T->nh_alloc){
		while(T->nh+n >= T->nh_alloc){ T->nh_alloc *= 2; }
		T->h = (Half*)realloc(T->h, sizeof(Half) * T->nh_alloc);
	}
}
static void DT_addface(function_sampler_2d T, int h){
	if(T->nf >= T->nf_alloc){
		T->nf_alloc *= 2;
		T->f = (Face*)realloc(T->f, sizeof(Face) * T->nf_alloc);
	}
	T->f[T->nf].h = h;
	T->f[T->nf].badness = 0;
	T->nf++;
}

int function_sampler_2d_add(
	function_sampler_2d T,
	const double *xy, double z, int id
){
	int h0 = T->hlast;
	int outside;
	const int ip = T->nv;
	
	if(NULL == T){ return -1; }

//function_sampler_1d_dump_state(sampler, stdout);
	
	if(z < T->z0){ T->z0 = z; }
	if(z > T->z1){ T->z1 = z; }
	if(xy[0] < T->minxy[0]){ T->minxy[0] = xy[0]; }
	if(xy[1] < T->minxy[1]){ T->minxy[1] = xy[1]; }
	if(xy[0] > T->maxxy[0]){ T->maxxy[0] = xy[0]; }
	if(xy[1] > T->maxxy[1]){ T->maxxy[1] = xy[1]; }
	
	/* Add the vertex*/
	if(T->nv >= T->nv_alloc){
		T->nv_alloc *= 2;
		T->v = (Vert*)realloc(T->v, sizeof(Vert) * T->nv_alloc);
	}
	T->v[T->nv].r[0] = xy[0];
	T->v[T->nv].r[1] = xy[1];
	T->v[T->nv].r[2] = z;
	T->v[T->nv].id = id;
	T->nv++;
	
	outside = DT_locate(T, xy, &h0);
	T->hlast = h0;
	if(outside){
		int hp = T->nh+1; /* Pointer to current prev boundary edge wrt p */
		int hn = T->nh+2; /* Pointer to current next boundary edge wrt p */
		int ia, ib;
		/* Find all visible edges, crawl out from h0 */
		int ht; /* temporary halfedge pointer that sits on the boundary */
		
		/* Add first triangle */
		ia = FROM(h0);
		ib = FROM(NEXT(h0));
		FLIP(h0) = T->nh+0;
		DT_addhalfs(T, 3);
		SETHALF(T->nh+0, hp   , h0, ib, T->nf);
		SETHALF(T->nh+1, hn   , -1, ia, T->nf);
		SETHALF(T->nh+2, T->nh, -1, ip, T->nf);
		DT_addface(T, T->nh);
		T->v[T->nv-1].h = hn;
		T->nh += 3;
		
		//update_edge(T, h0);
		
		ht = h0;
		do{ /* loop until not visible */
			/* advance ht forwards */
			ht = NEXT(ht);
			while(-1 != FLIP(ht)){
				ht = NEXT(FLIP(ht));
			}
			ia = FROM(ht);
			ib = FROM(NEXT(ht));
			if(orient2d(T->v[ia].r, T->v[ib].r, xy) < 0){ /* if visible */
				FLIP(ht) = T->nh+0;
				FLIP(hn) = T->nh+1;
				DT_addhalfs(T, 3);
				SETHALF(T->nh+0, T->nh+1, ht, ib, T->nf);
				SETHALF(T->nh+1, T->nh+2, hn, ia, T->nf);
				SETHALF(T->nh+2, T->nh+0, -1, ip, T->nf);
				DT_addface(T, T->nh);
				hn = T->nh+2;
				
				T->v[T->nv-1].h = hn;
				T->nh += 3;
				
				//update_edge(T, ht);
				//update_edge(T, hn);
			}else{ break; }
		}while(1);
		
		ht = h0;
		do{ /* loop until not visible */
			/* advance ht backwards */
			ht = NEXT(NEXT(ht));
			while(-1 != FLIP(ht)){
				ht = NEXT(NEXT(FLIP(ht)));
			}
			ia = FROM(ht);
			ib = FROM(NEXT(ht));
			if(orient2d(T->v[ia].r, T->v[ib].r, xy) < 0){ /* if visible */
				FLIP(ht) = T->nh+0;
				FLIP(hp) = T->nh+2;
				DT_addhalfs(T, 3);
				SETHALF(T->nh+0, T->nh+1, ht, ib, T->nf);
				SETHALF(T->nh+1, T->nh+2, -1, ia, T->nf);
				SETHALF(T->nh+2, T->nh+0, hp, ip, T->nf);
				DT_addface(T, T->nh);
				
				hp = T->nh+1;
				T->nh += 3;
				
				//update_edge(T, ht);
				//update_edge(T, hp);
			}else{ break; }
		}while(1);
	}else{ /* inside convex hull */
		const int h1 = NEXT(h0), h2 = NEXT(h1);
		const int ia = FROM(h0), ib = FROM(h1), ic = FROM(h2);
		const int f0 = FACE(h0);
		/*                    c
		 *                    +
		 *                   /|\
		 *                  / | \
		 *                 /  |  \
		 *                /  4|5  \
		 *               / fl1|fl0 \
		 *              /h2   |   h1\
		 *             /      +p     \
		 *            /  1 .-' `-. 2  \
		 *           /  .-'0 f0  3`-.  \
		 *          /.-'     h0      `-.\
		 *        a+---------------------+b
		 */
		/* connect p to vertices */
		DT_addhalfs(T, 6);
		SETHALF(T->nh+0, h0     , T->nh+1, ip, f0     );
		SETHALF(T->nh+1, T->nh+4, T->nh+0, ia, T->nf+1);
		SETHALF(T->nh+2, h1     , T->nh+3, ip, T->nf+0);
		SETHALF(T->nh+3, T->nh+0, T->nh+2, ib, f0     );
		SETHALF(T->nh+4, h2     , T->nh+5, ip, T->nf+1);
		SETHALF(T->nh+5, T->nh+2, T->nh+4, ic, T->nf+0);
		NEXT(h0) = T->nh+3;
		NEXT(h1) = T->nh+5; FACE(h1) = T->nf+0;
		NEXT(h2) = T->nh+1; FACE(h2) = T->nf+1;
		T->v[T->nv-1].h = T->nh;
		T->f[f0].h = h0;
		DT_addface(T, h1);
		DT_addface(T, h2);
		
		if(T->f[f0].badness > 0){ T->nbad--; }
		T->f[f0].badness = 0;
		
		//update_edge(T, h0);
		//update_edge(T, h1);
		//update_edge(T, h2);
		//update_edge(T, T->nh+0);
		//update_edge(T, T->nh+2);
		//update_edge(T, T->nh+4);
		
		T->nh += 6;
	}
	DT_fixup(T, T->v[T->nv-1].h);
	
	//DT_check(T);
	T->bad_valid = 0;
	return 0;
}
