% S4 Developer information
% Victor Liu (vkl@stanford.edu)
% Sep 17, 2011
<style type="text/css">
@import url(s4.css);
</style>

[S4 Home](index.html) | [Download](download.html) | [FAQ](faq.html) | [Lua API](s4_lua_api.html) | [Developer information](dev_info.html) | [Changelog](changelog.html)

# Developer information

## Compiling and installation

The usual `./configure; make; make install` sequence should just work. Without installing, the executable is in the `src/` directory.
There is also a provided `Makefile.custom` and `src/config.h.custom` (which needs to be renamed into config.h) which can be tailored to the particular installation environment, allowing one to bypass the configure process.

### Dependencies
* Lua 5.2, both headers and libraries.
* BLAS, or equivalent, library only.
* LAPACK (optional, 3.2 or later), library only.
* Pthreads (optional), both headers and libraries.
* FFTW3 (optional, 3.x), both headers and libraries.
* CHOLMOD (optional), both headers and libraries.

#### BLAS
The performance of S4 is directly dependent on the performance of the BLAS
library that it is linked against. Without a BLAS library, the default RNP
routines are used, which are suboptimal. It is recommended that either
vendor specific routines (Intel MKL or AMD AMCL) are used, or tuned
libraries such as ATLAS or GotoBLAS are used. Vendor specific routines have
incompatible calling/naming conventions than the standard Fortran interface,
so some porting effort is required. Note that GotoBLAS should be compiled
without threading if threading is enabled within S4; otherwise results may
not be correct.

#### LAPACK
If a precompiled LAPACK library is used, ensure ZGEEV is thread safe for
large matrices. Specifically, bug0061 in the LAPACK errata means that the
maximum stack space reserved by the Fortran compiler used to build LAPACK
must be sufficient to hold certain temporaries on the stack instead of in
static (shared) memory. The best way to avoid this with GFortran is to use
the -frecursive flag. A simple test of thread safety is provided in the
testing/Lapack_thread_test/ directory.

If no LAPACK library is available, the default RNP routines are used, which
are suboptimal due to lack of blocked routines. However, they are definitely
thread-safe.

#### CHOLMOD
Although not required, it is recommended when polarization decomposition
settings are enabled, since without it, a conjugate gradient solver is used
instead. The conjugate gradient solver does not always converge for fine
discretizations.

#### FFTW
This library is optional; the Kiss FFT library is used instead when FFTW
is not present. When available the version must be 3.x. This does not have
a great effect on simulation speed since the FFT is required only a handful
of times, and tends to be rather small in problem size.

### Porting

The code is developed in Windows under MinGW+MSYS and runs primarily in UNIX-like environments. For other platforms, the following list details the main problem areas.

1. In main.c near the top dealing with threads (this is not a problem if Pthreads are disabled).
2. In pattern/predicates.c, the settings may need modification to enable robust computations.
3. In numalloc.c, a more efficient aligned allocator may be supplied instead of the default.

## Program structure

S4 is built as a library on top of the Lua programming language.

S4 has 3 main layers:

=main.c=
	Specifies the Lua interface, which passes arguments on to the S4 layer.
=S4.cpp=
	The glue between the user specified information and the RCWA
	specification. The heavy lifting here is mainly determining which
	Fourier orders to use, and generating the Fourier dielectric matrices.
	The header S4.h is C-compatible and allows binding of a custom
	interface to the internals, or use programmatically as a component.
=rcwa.cpp=
	The real RCWA and S-matrix algorithms are contained here. The main
	functions are solving the layer band structure, composing S-matrices
	for stacks of layers, and solving for solutions given input fields.
	Also contained within here are functions for manipulating solution
	vectors and computing various quantities of interest.

Auxilliary functions:

=pattern.c/intersection.c=
	A self-contained library for representing layer patterning, as well as
	generating proper discretizations, Fourier transforms of the patterns,
	and quasiharmonic conformal flow fields.
=fmm=
	A set of different formulations to generate the Fourier components of
	the dielectric function.
=kiss_fft=
	A simple FFT library for performing numerical Fourier transforms of the
	conformal vector fields.
=RNP=
	A library of numerical linear algebra functions, and a thin layer over
	the BLAS and LAPACK libraries.
=SpectrumSampler/Interpolator=
	These are addtional convenience modules to perform commonly needed tasks.

## Areas of improvement

* **Robustness and efficiency of G-vector selection.** The current
  implementation is likely unnecessarily slow. Robustness of the algorithm
  for a very wide variety of lattices remains to be investigated.
* **Conjugate Gradient implementation in `pattern.c`.** The current CG
  implementation is a textbook implementation, with hard coded regularization.
* **Workspace pre-allocation for repeated calls.** Currently, many workspace
  parameters in `rcwa.cpp` are ignored, defaulting to automatic allocation.
  For repeated small calls, the allocation overhead may be substantial.
  Profiling is needed.
* **A full programming API in `S4.h`.** The current Lua interface reaches
  into the `Simulation` structure for certain functions. It is best to hide
  the details of the `Simulation` structure behind an opaque pointer and fix
  an API.
* **Parallelization over layer band structure computations.** For a large
  number of non-trivial layers, a rather simple parallelization would be to
  compute the bands of the layers in parallel. It is unclear whether this
  requires an invasive modification of S4.cpp.
* **Parallelization over S-matrix assembly.** The construction of the
  S-matrix should be parallelizable by a divide-and-conquer algorithm.
  Furthermore, it ought to be possible to cache intermediate results to
  speed up computing many S-matrices.
* **Graceful handling of diffraction threshold frequencies.** If a frequency
  happens to cause a zero eigenvalue in a layer's band structure, some
  internal matrices become singular, and no solution is obtainable. This
  requires a thorough theoretical investigation.

## Coding conventions

The code is mostly C-styled C++. The need to use C++ is mainly for the complex number type. There should never be non-trivial objects (only Plain Old Data structs), and certainly no inheritance, polymorphism, or templates (except in RNP).

Indentation should be 4 spaces per tab. Always use tabs instead of spaces. The actual tab width should never matter for readability, meaning tabs can only exist contiguously starting from the left-most column, and comment blocks should never sit at the end of lines of code. Lines should be kept to around 72 or 80 characters in length if possible, especially for comment blocks.

Functions generally are in Lapack-style, where there are a large number of well defined inputs, and a number of outputs returned by pointers. Functions return an integer code, with negative values corresponding to invalid parameters. When appropriate, workspaces can be passed in to reduce the number of dynamic allocations. Also, when convenient, workspace querying should be supported.

The code should compile cleanly with all warnings enabled. The only exemptions are external libraries like Kiss FFT or the geometric predicates sources.
