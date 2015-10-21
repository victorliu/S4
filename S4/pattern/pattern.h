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

#ifndef _PATTERN_H_
#define _PATTERN_H_

/* A pattern is defined as a set of 1D geometrical shapes (e.g. circles)
 * in the 2D plane. All shapes must be closed, simply connected contours.
 *
 * Available shapes:
 *   circle
 *   ellipse
 *   rectangle
 *   polygon
 *
 * A valid pattern must not have any shapes that intersect with each other.
 * Due to this requirement, each shape is either contained entirely within
 * another shape, or has no containing shape. The set of shapes then forms
 * a forest (a set of trees), where the parent-child relationship within a
 * tree is the immediate container-containee relationship.
 *
 * Currently, the no-intersection criterion is not checked due its
 * complexity, so the input shapes must be sanitized elsewhere.
 */
typedef enum shape_type_{
	CIRCLE,
	ELLIPSE,
	RECTANGLE,
	POLYGON
} shape_type;

typedef struct shape_{
	/* `type' determines which field in the union is used. */
	shape_type type;

	/* Location of the center of the shape. For circles, ellipses, and
	 * rectangles, the center is the geometric center (center of mass).
	 */
	double center[2]; /* x and y coordinates */

	/* The angle by which the shape is rotated, CCW, in radians.
	 * Obviously, the rotation is applied first, then the translation
	 * given in `center' is applied.
	 */
	double angle;

	/* The proper choice of field to use is determined by `type' above.
	 */
	union{
		struct{
			double radius;
		} circle;
		struct{
			double halfwidth[2]; /* x and y semimajor axes */
		} ellipse;
		struct{
			double halfwidth[2]; /* x and y halfwidths */
		} rectangle;
		struct{
			int n_vertices;
			/* Separately allocated, length 2*n_vertices.
			 * A list of xy-pairs in CCW orientation.
			 */
			double *vertex;
		} polygon;
	} vtab;
	int tag; /* available for user defined purposes */
} shape;

/* Convenience structure to pack up all the pattern information.
 * Fields `nshapes' and `shapes' are to be filled in, and `parent'
 * is to be computed using pattern_get_containment_tree below.
 */
typedef struct Pattern_{
	int nshapes;
	shape *shapes;
	int *parent;
} Pattern;

/* Computes the containment relationship of all the shapes.
 *  The shapes must not intersect one another, so that any shape is either
 * not contained within any other shape, or has a unique immediately
 * containing shape. The containing shape of shapes[i] is shapes[parent[i]]
 * unless parent[i] == -1, in which case shapes[i] is not contained within
 * any other shape.
 *
 * Arguments:
 *    nshapes  IN      Number of shapes; length of `shapes'.
 *    shapes   IN/OUT  Array of `shape' structures. On exit, the array
 *                     is sorted in order of decreasing area of shapes.
 *    parent   OUT     Length `nshapes'. Gives the immediate containing
 *                     shape of each shape, or -1 if none.
 * Return values:
 *                           0: If no shapes intersect each other and containment tree
 *                              returned is valid.
 *          0 < n <= nshapes  : If shapes[n-1] is invalid.
 *    nshapes < n <= 2*nshapes: If shapes[n-1] intersects with another shape.
 *                          -n: If n-th argument is invalid.
 */
int pattern_get_containment_tree(
	int nshapes,
	shape *shapes,
	int *parent
);
/* Convenience version of the above.
 * p->nshapes and p->shapes must be filled in.
 */
int Pattern_GetContainmentTree(
	Pattern *p
);

/* Determines the child-most (smallest area) shape which contains a given
 * point. Optionally returns an approximate outward normal vector to the
 * shape.
 *
 * Arguments:
 *    nshapes      IN   Number of shapes; length of `shapes'.
 *    shapes       IN   Array of `shape' structures in order of decreasing
 *                      area. (as output from pattern_get_containment_tree).
 *    parent       IN   Immediate containing shape relationship for each
 *                      shape. (as output from pattern_get_containment_tree).
 *    x            IN   The point at which to determine the smallest-area
 *                      shape. (x and y coordinates)
 *    shape_index  OUT  The returned shape as an index into `shapes'.
 *    n            OUT  Optional. If non-NULL, the normal vector.
 * Return values:
 *    0: If *shape_index contains the sought shape index.
 *    1: If `x' is not in any shape; *shape_index not referenced.
 *   -n: If n-th argument is invalid.
 */
int pattern_get_shape(
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double x[2],
	int *shape_index,
	double n[2]
);
/* Convenience version of the above. */
int Pattern_GetShape(
	const Pattern *p,
	const double x[2],
	int *shape_index,
	double n[2]
);

/* Returns an approximate outward normal vector to the shape.
 *
 * Arguments:
 *    shape        IN   `shape' structure whose normal to compute
 *    x            IN   The point at which to determine the normal
 *    n            OUT  The normal vector.
 * Return values:
 *   -n: If n-th argument is invalid.
 */
int shape_get_normal(const shape *s, const double x[2], double n[2]);

/* Returns the Fourier transform of the pattern, given values for the interior
 * of each shape and the background.
 *  The Fourier transform convention is
 *          F(kx,ky) = Integral f(x,y) exp(-i*(kx*x + ky*y)) dx dy
 *
 * Arguments:
 *    nshapes      IN   Number of shapes; length of `shapes'.
 *    shapes       IN   Array of `shape' structures in order of decreasing
 *                      area. (as output from pattern_get_containment_tree).
 *    parent       IN   Immediate containing shape relationship for each
 *                      shape. (as output from pattern_get_containment_tree).
 *    value        IN   Length `2*(nshapes+1)'. Real and imaginary parts of
 *                      of the interior value of each shape. The background
 *                      value is `value[0]+i*value[1]', and the value of
 *                      shape[i] is `value[2*(i+1)+0]+i*value[2*(i+1)+1]'.
 *    f            IN   Point in reciprocal space to obtain the F.T., but
 *                      coordinates should be divided by 2*pi (e.g. On a unit
 *                      square real space lattice, the point (1,0) should be
 *                      used instead of (2*pi,0).
 *  unit_cell_size IN Area of the unit cell in real space.
 *    FT           OUT  Real and imaginary part of the F.T. at `2*pi*f'.
 * Return values:
 *    0: If successful.
 *   -n: If n-th argument is invalid.
 */
int pattern_get_fourier_transform(
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double *value,
	const double f[2],
	int ndim,
	double unit_cell_size,
	double FT[2]
);
/* Convenience version of the above. */
int Pattern_GetFourierTransform(
	const Pattern *p,
	const double *value,
	const double f[2],
	int ndim,
	double unit_cell_size,
	double FT[2]
);

/* Returns an area weighting of each shape within one cell of a uniform
 * discretization of the origin-centered unit square.
 *  The area-fraction of each shape within the rectangle
 * [x0,x0+dx] x [y0,y0+dy] is returned in `value'. The discretization has
 * `nx' cells along the x-direction, and `ny' cells along the y-direction.
 * The total region is the unit square [-0.5,0.5] x [-0.5,0.5], so
 *            x0 = -0.5, 0.5 + dx, 0.5 + 2*dx, ..., 0.5 - dx
 * where dx = 1/nx. Similarly for y0.
 *
 * Arguments:
 *    nshapes      IN   Number of shapes; length of `shapes'.
 *    shapes       IN   Array of `shape' structures in order of decreasing
 *                      area. (as output from pattern_get_containment_tree).
 *    parent       IN   Immediate containing shape relationship for each
 *                      shape. (as output from pattern_get_containment_tree).
 *    L            IN   Lattice basis vectors, (L[0],L[1]) are x- and y-
 *                      coordinates of the first (u) basis vector, and
 *                      (L[2],L[3]) are the coordinates of the second (v)
 *                      basis vector.
 *    nu, nv       IN   Number of cells in u and v directions. nu > 0, nv > 0.
 *    iu, iv       IN   u and v index of the cell at which the area fraction
 *                      is to be computed. 0 <= iu < nu, 0 <= iv < nv.
 *    value        OUT  Length `nshapes+1'. value[0] is the area fraction of
 *                      the background (no shape). The area fraction of
 *                      shape[i] is stored in value[i+1].
 * Return values:
 *    0: If successful.
 *   -n: If n-th argument is invalid.
 */
int pattern_discretize_cell(
	int nshapes,
	const shape *shapes,
	const int *parent,
	const double L[4],
	int nu, int nv,
	int iu, int iv,
	double *value
);
/* Convenience version of the above. */
int Pattern_DiscretizeCell(
	const Pattern *p,
	const double L[4],
	int nu, int nv,
	int iu, int iv,
	double *value
);

/* Returns a vector field which is tangential or normal to the shapes.
 *  The field is generated as the minimizer of the Dirichlet 1-form energy
 * on the specified uniform rectangular grid subject to the constraints
 * that the resulting vector field be equal to the tangent or normal vectors
 * of the parameterizations of the shape contours, with periodic boundary
 * conditions over the origin-centered unit square.
 *
 * Arguments:
 *    nshapes      IN   Number of shapes; length of `shapes'.
 *    shapes       IN   Array of `shape' structures in order of decreasing
 *                      area. (as output from pattern_get_containment_tree).
 *    parent       IN   Immediate containing shape relationship for each
 *                      shape. (as output from pattern_get_containment_tree).
 *    type         IN   0 for a field tangential to shapes, 1 for normal using
 *                      shape centers as sources.
 *    L            IN   Lattice basis vectors, (L[0],L[1]) are x- and y-
 *                      coordinates of the first (u) basis vector, and
 *                      (L[2],L[3]) are the coordinates of the second (v)
 *                      basis vector.
 *    nu, nv       IN   Number of cells in u and v directions. nu > 0, nv > 0.
 *    field        OUT  Length 2*nu*nv, array of xy-vectors. The field vector
 *                      located at uv-index (i,j) is
 *                      { field[2*(i+j*nu)+0], field[2*(i+j*nu)+1] }
 * Return values:
 *    0: If successful.
 *    1: If the discretization is not fine enough to properly apply
 *       constraints, although a best-attempt field was generated.
 *   -n: If the n-th argument is invalid.
 *
 * Notes:
 *   The details of implementation depend on the lattice type.
 *    In the case of an orthogonal lattice, we use a rectangular quad mesh.
 *   We use an interpolation scheme on the dual mesh to compute the
 *   constraints, which is piecewise constant over each edge. The Laplacian
 *   is then just a generalized (4,-1,-1,-1,-1) stencil.
 *    In the case of a non-orthogonal lattice, we split the parallelogramic
 *   discretization cells along the short diagonal into two triangles. We
 *   assume that the triangulation is well centered and use a circumcentric
 *   dual mesh with a diagonal Hodge star to construct the Laplacian. The
 *   constraints are enforced using the primal Whitney 1-form interpolator.
 *    In both cases, the resulting field is computed by the proper
 *   interpolation scheme at the centers of the parallelogramic discretization
 *   cells. The constraints are enforced using weighted least squares with
 *   uniform weights per edge (not per edge-shape intersection).
 *    To deal with shapes extending beyond the origin-centered fundamental
 *   parallelogram of the lattice, edge-shape intersections are computed using
 *   each shape and its 8 nearest neighbor unit cell copies.
 */
int pattern_generate_flow_field(
	int nshapes,
	const shape *shapes,
	const int *parent,
	int type,
	const double L[4],
	int nu, int nv,
	double *field
);
/* Convenience version of the above. */
int Pattern_GenerateFlowField(
	const Pattern *p,
	int type,
	const double L[4],
	int nu, int nv,
	double *value
);

#endif /* _PATTERN_H_ */
