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

#ifndef _CONVERT_H_
#define _CONVERT_H_

/* Converts photonics quantities.
 * Supported units:
 *   Lengths: um, nm, m, cm, mm
 *   Energies: eV, J
 *   Frequencies: THz, GHz, Hz, rad/s
 * Relations used:
 *   c = f*l, E=h*f
 */
 
/* Returns 0 on success,
          -n if n-th argument is invalid,
		   1 if from_units not supported,
		   2 if to_units not supported
 */
int convert_units(double *value, const char *from_units, const char *to_units);

#endif /* _CONVERT_H_ */
