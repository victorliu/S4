% S4
% Victor Liu (vkl@stanford.edu)
% Sep 17, 2011
<style type="text/css">
@import url(s4.css);
</style>

[S4 Home](index.html) | [Download](download.html) | [FAQ](faq.html) | [Lua API](s4_lua_api.html) | [Developer information](dev_info.html) | [Changelog](changelog.html)

# License and Copyright

S4 is copyright (c) 2009-2011, Stanford University.

S4 is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this library; if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the GNU web site:

->[http://www.gnu.org/copyleft/gpl.html](http://www.gnu.org/copyleft/gpl.html)<-

The files `cubature.c` and `cubature.h` were taken from Stephen G. Johnson's Cubature package. It is copyright (c) 2005-2010 by Steven G. Johnson.

The files `sort.c` and `sort.h` were taken from Rusty Russell's modification of the original libc routines. See the [CCAN page](http://ccan.ozlabs.org/info/asort.html) for details.

The file `predicates.c` is public domain, written by Jonathan Richard Shewchuk. See the [Quake project page](http://www.cs.cmu.edu/~quake/robust.html) for details.

Distributed with S4 is the Kiss FFT library. It is copyright (c) 2003-2011 by Mark Borgerding.

As a clarification, we should note that Lua script (`.lua`) files, written by the user (i.e. not containing code distributed with S4) and loaded at runtime by the S4 software, are not considered derived works of S4 and do not fall thereby under the restrictions of the GNU General Public License.

In addition, all of the example Lua code in these pages, as well as the example Lua files in the `examples/` directory, may be freely used, modified, and redistributed, without any restrictions. (The warranty disclaimer still applies, of course.)

# Referencing

We kindly ask you to reference the S4 package and its authors in any publication for which you used S4. (You are not legally required to do so; it is up to your common sense to decide whether you want to comply with this request or not.) We are working on putting out a paper on the package, but in the meantime the preferred citation is something like:

> Simulations were performed with the Fourier Modal Method (FMM) [ref FMM], using a freely available software package [ref S4].

For the FMM reference, you might use, for example, Lifeng Li, "New formulation of the Fourier modal method for crossed surface-relief gratings," J. Opt. Soc. Am. A **14**, p. 2758-2767 (1997). For referencing S4, currently refer to this website and its authors Victor Liu and Shanhui Fan.
