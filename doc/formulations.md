% S4 Formulations
% Victor Liu (vkl@stanford.edu)
% August 31, 2011
<style type="text/css">
@import url(s4.css);
</style>

[S4 Home](index.html) | [Download](download.html) | [FAQ](faq.html) | [Lua API](s4_lua_api.html) | [Developer information](dev_info.html) | [Changelog](changelog.html)

# Fourier Modal Method formulations

There has been extensive literature on the best way to generate the Fourier series coefficients for the in-plane dielectric profiles of each layer.
S4 implements a number of different formulations.
The following functions determine which formulation is selected:

* [UseDiscretizedEpsilon](s4_lua_api.html#S4_Simulation_UseDiscretizedEpsilon)
* [UsePolarizationDecomposition](s4_lua_api.html#S4_Simulation_UsePolarizationDecomposition)
* [UseSubpixelSmoothing](s4_lua_api.html#S4_Simulation_UseSubpixelSmoothing)
* [UseJonesVectorBasis](s4_lua_api.html#S4_Simulation_UseJonesVectorBasis)
* [UseNormalVectorBasis](s4_lua_api.html#S4_Simulation_UseNormalVectorBasis)

In addition, the following functions control accuracy and the lattice truncation:

* [SetResolution](s4_lua_api.html#S4_Simulation_SetResolution)
* [SetLatticeTruncation](s4_lua_api.html#S4_Simulation_SetLatticeTruncation)

To simplify the choice for users, the table below summarizes the recommended settings.
It is recommended to always use circular truncation unless there is a good reason to do otherwise.
Speed indicates the speed of the Fourier coefficient generation, which is usually not the dominant part of the simulation time.

Options         | Can handle Anisotropic? | Recommended resolution | Speed   | Accuracy
----------------|-------------------------|------------------------|---------|----------
none            | yes                     | -                      |  fast   |   poor
Disc            | yes                     | 8                      |  medium |   poor
Subpixel        | yes                     | 4                      |  medium |  medium
Pol             | no                      | 8                      |  slow   |   good
Pol+Normal      | no                      | 8                      |  slow   |   good
Pol+Jones       | no                      | 8                      |  slow   |   good
Disc+Pol        | no*                     | 4                      |  slow   |  medium
Disc+Pol+Normal | no*                     | 4                      |  slow   |  medium
Disc+Pol+Jones  | no*                     | 4                      |  slow   |  medium

*: The formulation does not strictly work correctly for anisotropic media
   however it may still work.
   Proper support for anisotropic materials is in principle possible. There
   are currently no plans for implementing generation of the proper basis
   fields for this feature.
