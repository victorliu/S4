Units & Coordinates
===================

|S4| solves the linear Maxwell's equations, which are scale-invariant.
Therefore, |S4| uses normalized units, so the conversion between the numbers in |S4| and physically relevant numbers can sometimes be confusing.
Here we show how to perform these conversions and provide some examples.

Length, time, and frequency
---------------------------

The speed of light and the vacuum permittivity and vacuum permeability are all normalized to unity in |S4|.
Because of this, time is measured in units of length and frequency is measured in units of inverse length.
Due to the scale invariant nature of the linear Maxwell's equations, there is no intrinsic length unit.
Instead, one can imagine that all lengths specified in |S4| are multiples of a common reference length unit.

When the lattice vectors are set, this determines a length scale.
Suppose we have a square lattice with a periodicity of 680nm. It might be logical then to choose 1 micron as the base length unit, and to specify all lengths in microns. The lattice would then be set to the vectors (0.68, 0) and (0, 0.68).
Since frequency is in units of inverse length, then ``SetFrequency(1)`` corresponds to a wavelength of 1um or a physical frequency of (c/1um) = 300 THz, and ``SetFrequency(1.1)`` corresponds to a wavelength of (1/1.1) = 909nm, etc.

For numerical stability, it is best to choose a length scale in which most dimensions are close to unity. For typical nanophotonic structures, using (implicit) microns is preferable to using (implicit) nanometers.

Amplitudes and powers
---------------------

Typically |S4| is used to compute transmission or reflection spectra, so the absolute units of incident power or amplitude are irrelevant. However, sometimes, it is necessary to be able to translate these figures into physical units.

Under a normally-incident planewave excitation, the incident power is 1 if the amplitude is 1, regardless of unit cell size.
At off-normal incidence, the incident power is reduced by a factor of cos(phi), where phi is the polar angle (angle from normal incidence).

The specified amplitudes are for the electric field.
If the incident power is considered to have units of Watts per unit-cell-area, then the electric field has units of Volts per square-root-of-unit-cell-area.

.. |S4| replace:: S\ :sup:`4`
