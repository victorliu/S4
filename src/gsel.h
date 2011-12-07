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

// Purpose:
//   Selects the G vectors for a given lattice Lk, and a truncation order NG
//   The returned number of vectors in G is always less than or equal to NG.
//
// Arguments:
//   method - (INPUT) Type of truncation to use.
//            0 = circular truncation
//            1 = parallelogramic truncation
//   NG     - (IN/OUT) On input, the upper limit of the number of G vectors
//            to return (capacity of G). On exit, the number of G vectors
//            returned in G
//   Lk     - (INPUT) Primitive lattice vectors. (Lk[0],Lk[1]) are the
//            cartesian coordinates of the first lattice vector, etc.
//   G      - (OUTPUT) - Size 2*(*NG). Stores the returned G vectors as
//            successive pairs of integer coefficients of the lattice basis.
// Returns
//   -n if the n-th argument is invalid
//   0 on success
int G_select(const int method, unsigned int *NG, const double Lk[4], int *G);
