% S4 Known issues
% Victor Liu (vkl@stanford.edu)
% Sep 17, 2011
<style type="text/css">
@import url(s4.css);
</style>

[S4 Home](index.html) | [Download](download.html) | [FAQ](faq.html) | [Lua API](s4_lua_api.html) | [Developer information](dev_info.html) | [Changelog](changelog.html)

# Known issues

* Vector field generation in pattern.c fails for large discretization sizes without CHOLMOD.

  The iteration limit of the conjugate gradient solver may need to scale up
  faster at large grid sizes.

* Crash on unknown material name in GetEpsilon.

  There is no stringent checking to make sure all materials referenced in
  layers are defined. This needs to be fixed.

* Incorrect results when using multithreading (`SolveInParallel`).

  This is most likely an indication of a poorly compiled LAPACK library.
  Make sure to compile LAPACK with no static variables. This means using
  the `-frecursive` flag with gfortran, or else manually patching up the
  functions which get called (mainly `zgeev` and its callees).

* Vector field generation in pattern.c does not correctly handle shapes
  which have portions extremely far outside the unit cell.
  
  **Workaround**: When defining shapes, do not allow the extent of the
  shapes to extend beyond the 8 neighbors of the real space lattice's
  (origin centered) fundamental parallelogram.
  
* Vector field generation in pattern.c currently does not support generating
  everywhere-normal fields with natural boundary conditions.
  
  **Proposal**: Implementation of natural boundaries and interpolated source
  placement. Constraint computation should be similar to the
  tangential case, except cross products are taken with the pixel
  edges themselves.
  
* Vector field generation in pattern.c currently fails silently or produces
  poor output on non-orthogonal lattices with substantially different basis
  vector length to sampling ratios.
  
  **Proposal**: The problem is partially mitigated by the use of circular
  G-vector truncation; the length-to-sampling ratios will usually be
  approximately equal. Correctly supporting the general case requires
  generating a non-trivial implicit mesh for the fundamental parallelogram.
  The resulting sparse matrix system would then become substantially larger
  and the overall complexity of the field generation would be much greater.
  No attempt to fix will be made any time soon.
  
* Dipole excitation source is not implemented in main.c and S4.cpp.

  **Proposal**: Dipole sources require knowledge of (px,py,pz) instead of (hx,hy)
  in the planewave case. It may be possible to express the dipole
  as a set of tangential H field values. Further planning is needed.

