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

#include <stdlib.h>
#include <stdio.h>

#ifdef WIN32
# include <malloc.h>
void * _aligned_malloc(size_t size, size_t alignment);
void _aligned_free(void *ptr);
#else
#include <inttypes.h>
typedef uintptr_t malloc_aligned_ULONG_PTR;
#endif

// size : the size of allocated memory
//        The actual size of allocation will be greater than this size.
// alignment : the alignment boundary
void *malloc_aligned(size_t size, size_t alignment){
#ifdef WIN32
	return (void*)_aligned_malloc(size, alignment);
#else
	void *pa, *ptr;

	//pa=malloc(((size+alignment-1)&~(alignment-1))+sizeof(void *)+alignment-1);

	pa = malloc((size+alignment-1)+sizeof(void*));
	if(!pa){ return NULL; }

	ptr = (void*)( ((malloc_aligned_ULONG_PTR)pa+sizeof(void*)+alignment-1)&~(alignment-1) );
	*((void **)ptr-1) = pa;
	
	return ptr;
#endif
}

void free_aligned(void *ptr){
#ifdef WIN32
	_aligned_free(ptr);
#else
	if(ptr){
		free(*((void **)ptr-1));
	}
#endif
}
