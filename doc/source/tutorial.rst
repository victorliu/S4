Tutorial
========

Basic Usage
-----------

This section assumes that the Lua frontend is used.
The program may be run with a Lua script as an argument, in which case the script is run and then the program terminates, or without an argument, in which case it enters interactive mode with a Lua shell.

The basic flow of a script is listed below.
To find more information about these functions and for a full listing, see the :ref:`lua-api-label`.

#. Obtain a new simulation object::

	S = S4.NewSimulation()

   S now contains a simulation object with a blank specification, and no solutions.

#. Set the lattice and number of basis functions::

	S:SetLattice({ux,uy}, {vx,vy})
	S:SetNumG(100)

   This sets the periodicity based on the lattice vectors (ux,uy) and (vx,vy).
   The number of basis functions is set to 100.

#. Define all materials::

	S:AddMaterial('Vacuum', 1)
	S:AddMaterial('Silicon', {12, 0.01})

   This adds a material called "Vacuum" with a dielectric constant of 1, and a material called "Silicon", with a dielectric constant of 12+0.01i.

#. Add all layers::

	S:AddLayer('top', 0, 'Vacuum')
	S:AddLayer('slab', 0.5, 'Silicon')
	S:AddLayerCopy('bottom', 0, 'top')

   Here we are adding a semi-infinite Vacuum layer (the thickness is largely irrelevant), a middle slab layer made of silicon with a thickness of 0.5, and finally a bottom semi-infinite layer that is a copy of the top layer.
   It is always preferable to copy layers if possible since it reduces the computational cost.

#. Add patterning to layers::

	S:SetLayerPatternCircle('slab',   -- layer to pattern
	                        'Vacuum', -- material inside circle
	                        {0,0},    -- center of circle
	                        0.2)      -- radius of circle

#. Specify the excitation mechanism::

	S:SetExcitationPlanewave(
		{0, 0}, -- phi in [0,180), theta in [0,360)
		{1, 0}, -- s-polarization amplitude and phase in degrees
		{0, 0}) -- p-polarization

   This sets the excitation to a normally incident planewave with amplitude 1.

#. Specify the operating frequency::

	S:SetFrequency(0.4)

#. Obtain desired output::

	transmission = S:GetPowerFlux('bottom')
	print(transmission)

   This obtains the forward power flux in the bottom layer, which is also the transmitted power.

.. _fmm-formulations-label:

Fourier Modal Method formulations 
---------------------------------

There has been extensive literature on the best way to generate the Fourier series coefficients for the in-plane dielectric profiles of each layer. S4 implements a number of different formulations. The following functions determine which formulation is selected:

* UseDiscretizedEpsilon
* UsePolarizationDecomposition
* UseSubpixelSmoothing
* UseJonesVectorBasis
* UseNormalVectorBasis

In addition, the following functions control accuracy and the lattice truncation:

* SetResolution
* SetLatticeTruncation

To simplify the choice for users, the table below summarizes the recommended settings. It is recommended to always use circular truncation unless there is a good reason to do otherwise. Speed indicates the speed of the Fourier coefficient generation, which is usually not the dominant part of the simulation time.

+-----------------+-------------------------+------------------------+--------+----------+
| Options         | Can handle Anisotropic? | Recommended resolution | Speed  | Accuracy |
+=================+=========================+========================+========+==========+
| none            | yes                     | N/A                    | fast   | poor     |
+-----------------+-------------------------+------------------------+--------+----------+
| Disc            | yes                     | 8                      | medium | poor     |
+-----------------+-------------------------+------------------------+--------+----------+
| Subpixel        | yes                     | 4                      | medium | medium   |
+-----------------+-------------------------+------------------------+--------+----------+
| Pol             | no                      | 8                      | slow   | good     |
+-----------------+-------------------------+------------------------+--------+----------+
| Pol+Normal      | no                      | 8                      | slow   | good     |
+-----------------+-------------------------+------------------------+--------+----------+
| Pol+Jones       | no                      | 8                      | slow   | good     |
+-----------------+-------------------------+------------------------+--------+----------+
| Disc+Pol        | no [#f1]_               | 4                      | slow   | medium   |
+-----------------+-------------------------+------------------------+--------+----------+
| Disc+Pol+Normal | no [#f1]_               | 4                      | slow   | medium   |
+-----------------+-------------------------+------------------------+--------+----------+
| Disc+Pol+Jones  | no [#f1]_               | 4                      | slow   | medium   |
+-----------------+-------------------------+------------------------+--------+----------+

.. rubric:: Footnotes

.. [#f1] The formulation does not strictly work correctly for anisotropic media however it may still work. Proper support for anisotropic materials is in principle possible. There are currently no plans for implementing generation of the proper basis fields for this feature.

Examples
--------

The source distribution of S4 includes numerous fully working didactic examples as well as examples replicating published results.
You can find these examples in the ``examples/`` directory of the source distribution.

.. |S4| replace:: S\ :sup:`4`
