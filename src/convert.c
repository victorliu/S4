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

#include <string.h>

int convert_units_from_um(double *value, const char *to_units){
	if(0 == strcmp(to_units, "um")){
		return 0;
	// Lengths
	}else if(0 == strcmp(to_units, "nm")){
		*value *= 1e3;
	}else if(0 == strcmp(to_units, "m")){
		*value *= 1e-6;
	}else if(0 == strcmp(to_units, "cm")){
		*value *= 1e-4;
	}else if(0 == strcmp(to_units, "mm")){
		*value *= 1e-3;
	// Energies
	}else if(0 == strcmp(to_units, "eV")){
		*value = 1.2398412230660532 / *value;
	}else if(0 == strcmp(to_units, "J")){
		*value = 1.9864455003959036e-19  / *value;
	// Frequencies
	}else if(0 == strcmp(to_units, "THz")){
		*value = 299.792458 / *value;
	}else if(0 == strcmp(to_units, "GHz")){
		*value = 299792.458 / *value;
	}else if(0 == strcmp(to_units, "Hz")){
		*value = 2.99792458e14 / *value;
	}else if(0 == strcmp(to_units, "rad/s")){
		*value = 1.8836515673088532e15 / *value;
	}else{
		return 2;
	}
	return 0;
}

int convert_units_from_eV(double *value, const char *to_units){
	if(0 == strcmp(to_units, "eV")){
		return 0;
	// Lengths
	}else if(0 == strcmp(to_units, "um")){
		*value = 1.2398412230660532 / *value;
	}else if(0 == strcmp(to_units, "nm")){
		*value = 1.2398412230660532e3 / *value;
	}else if(0 == strcmp(to_units, "m")){
		*value = 1.2398412230660532e-6 / *value;
	}else if(0 == strcmp(to_units, "cm")){
		*value = 1.2398412230660532e-4 / *value;
	}else if(0 == strcmp(to_units, "mm")){
		*value = 1.2398412230660532e-3 / *value;
	// Energies
	}else if(0 == strcmp(to_units, "J")){
		*value *= 1.60217733e-19;
	// Frequencies
	}else if(0 == strcmp(to_units, "THz")){
		*value *= 241.7990726737019;
	}else if(0 == strcmp(to_units, "GHz")){
		*value *= 241799.0726737019;
	}else if(0 == strcmp(to_units, "Hz")){
		*value *= 2.4179907267370188e14;
	}else if(0 == strcmp(to_units, "rad/s")){
		*value *= 1.5192683807130525e15;
	}else{
		return 2;
	}
	return 0;
}


int convert_units_from_THz(double *value, const char *to_units){
	if(0 == strcmp(to_units, "THz")){
		return 0;
	// Lengths
	}else if(0 == strcmp(to_units, "um")){
		*value = 2.99792458e2 / *value;
	}else if(0 == strcmp(to_units, "nm")){
		*value = 2.99792458e5 / *value;
	}else if(0 == strcmp(to_units, "m")){
		*value = 2.99792458e-4 / *value;
	}else if(0 == strcmp(to_units, "cm")){
		*value = 2.99792458e-2 / *value;
	}else if(0 == strcmp(to_units, "mm")){
		*value = 2.99792458e-1 / *value;
	// Energies
	}else if(0 == strcmp(to_units, "eV")){
		*value *= 0.004135665157613982;
	}else if(0 == strcmp(to_units, "J")){
		*value *= 6.6260689599999995e-22;
	// Frequencies
	}else if(0 == strcmp(to_units, "GHz")){
		*value *= 1e3;
	}else if(0 == strcmp(to_units, "Hz")){
		*value *= 1e12;
	}else if(0 == strcmp(to_units, "rad/s")){
		*value *= 6.283185307179586e12;
	}else{
		return 2;
	}
	return 0;
}

int convert_units(double *value, const char *from_units, const char *to_units){
	if(NULL == value){ return -1; }
	if(NULL == from_units){ return -2; }
	if(NULL == to_units){ return -3; }
	
	if(0 == strcmp(from_units, to_units)){
		return 0;
	// Lengths
	}else if(0 == strcmp(from_units, "um")){
		return convert_units_from_um(value, to_units);
	}else if(0 == strcmp(from_units, "nm")){
		*value *= 1e-3;
		return convert_units_from_um(value, to_units);
	}else if(0 == strcmp(from_units, "m")){
		*value *= 1e6;
		return convert_units_from_um(value, to_units);
	}else if(0 == strcmp(from_units, "cm")){
		*value *= 1e4;
		return convert_units_from_um(value, to_units);
	}else if(0 == strcmp(from_units, "mm")){
		*value *= 1e3;
		return convert_units_from_um(value, to_units);
	// Energies
	}else if(0 == strcmp(from_units, "eV")){
		return convert_units_from_eV(value, to_units);
	}else if(0 == strcmp(from_units, "J")){
		*value *= 6.241506363094027e18;
		return convert_units_from_eV(value, to_units);
	// Frequencies
	}else if(0 == strcmp(from_units, "THz")){
		return convert_units_from_THz(value, to_units);
	}else if(0 == strcmp(from_units, "GHz")){
		*value *= 1e-3;
		return convert_units_from_THz(value, to_units);
	}else if(0 == strcmp(from_units, "Hz")){
		*value *= 1e-12;
		return convert_units_from_THz(value, to_units);
	}else if(0 == strcmp(from_units, "rad/s")){
		*value *= 1.5915494309189535e-13;
		return convert_units_from_THz(value, to_units);
	}else{
		return 1;
	}
}
