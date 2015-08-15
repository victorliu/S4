.. default-domain:: py

.. _python-api-label:

Python API reference
====================

|S4| can be built as a `Python <http://python.org>`_ extension, in addition to the original `Lua <http://www.lua.org>`_ interface.
The current Python interface is not as fully featured as the Lua interface, but it should ultimately achieve feature parity.
Much auxiliary functionality, such as numerical integration, is not included here since `Numpy <http://www.numpy.org/>`_ and `Scipy <http://www.scipy.org/>`_ can easily be used instead.

S4 module
---------

.. module:: S4

All top level functions of |S4| are located in the ``S4`` library.

.. function:: New(Lattice, NumBasis)

   Returns a new `Simulation object`_.
   
   Usage::
   
	S = S4.New(Lattice=((1,0),(0,1)), NumBasis=100)
   
   Arguments:
   
	Lattice
		Sets the lattice vectors for the structure. This can be a single real number in the case of a 1D periodic lattice, or a pair of vectors for a 2D periodic lattice.
	NumBasis
		Sets the maximum number of in-plane (x and y) Fourier expansion orders to use.
		All fields and eigenmodes of the system use the same Fourier basis of the same dimension.
		The computation time is roughly proportional to the cube of this number, and the memory usage is roughly proportional to the square.
		This is equivalent to the "NumG" setting in the Lua interface.

   Return values:

	S
		A new `Simulation object`_.


Simulation object
-----------------

.. class:: Simulation

The Simulation object is the primary object which computes solutions to systems.
When a new Simulation object is requested from :func:`New`, all settings are in a blank state, with no materals, layers, or excitation.
When solutions are requested, only a minimal set of internal computations are performed in order to satisfy the request.

Parameter specification
^^^^^^^^^^^^^^^^^^^^^^^

.. method:: Simulation.SetMaterial(Name, Epsilon)

	Updates an existing material with a new dielectric constant or adds a material if none exists.

   Usage::

	S.SetMaterial(Name = 'Silicon, Epsilon = 12+0.01j)
	S.SetMaterial(Name = 'Silicon', Epsilon = (
		(12+0.01j, 0, 0),
		(0, 12+0.01j, 0),
		(0, 0, 12+0.01j)
		))

   Arguments:

	Name
		(string) The name of the material to update, or the name of a new material if no material by that name exists.
	Epsilon
		(number, complex number, or 3x3 tensor) The relative permittivity of the material. The imaginary part should generally be positive for lossy materials.

   Return values:

	None

.. method:: Simulation.AddLayer(Name, Thickness, Material)

	Adds a new unpatterned layer with a specified thickness and material.

   Usage::

	S.AddLayer(Name = 'slab', Thickness = 0.6, Material = 'Silicon')

   Arguments:

	Name
		(string) The name of the layer. Each layer should have a unique name if it is to be referenced later.
	Thickness
		(number) The thickness of the layer.
	Material
		(string) The name of the material which comprises the layer. With patterning, this is the default (background) material of the layer.

   Return values:

	None

.. method:: Simulation.AddLayerCopy(Name, Thickness, Layer)

	Adds a new layer with a specified thickness, but identical patterning as another existing layer.
	Note that this merely creates a reference to the copied layer; further patterning of the copied layer also affects the new layer. Additionally, a copy of a copy cannot be made.

   Usage::

	S.AddLayerCopy(Name = 'slab2', Thickness = 0.5, Layer = 'slab')

   Arguments:

	Name
		(string) The name of the new layer, different from the layer being copied.
	Thickness
		(number) The thickness of the new layer.
	Layer
		(string) The name of the layer which whose pattern is to be copied. That layer cannot itself be a copy of a layer.

   Return values:

	None

.. method:: Simulation.SetLayerThickness(Layer, Thickness)

	Updates an existing layer with a new thickness.
	Previously cached layer eigenmodes are preserved, making this function the preferred way to update a layer's thickness.

   Usage::

	S.SetLayerThickness(Layer = 'slab', thickness)

   Arguments:

	Layer
		(string) The name of the layer to update.
	Thickness
		(number) The new thickness of the layer.

   Return values:

	None

.. method:: Simulation.RemoveLayerRegions(Layer)

	Removes all layer regions from an existing layer.

   Usage::

	S.RemoveLayerRegions(Layer = 'slab')

   Arguments:

	Layer
		(string) The name of the layer to modify.

   Return values:

	None

.. method:: Simulation.SetRegionCircle(Layer, Material, Center, Radius)

	Adds a (filled) circle of a specified material to an existing non-copy layer.
	The circle should not intersect any other patterning shapes, but may contain or be contained within other shapes.

   Usage::

	S.SetRegionCircle(
		Layer = 'slab',
		Material = 'Vacuum',
		Center = (0,0),
		Radius = 0.2
	)

   Arguments:

	Layer
		(string) The name of the layer to pattern. This layer cannot be a copy of another layer.
	Material
		(string) The name of the material which fills the interior of the circle.
	Center
		(pair of numbers) x- and y-coordinates of the center of the circle relative to the center of the unit cell (the origin).
	Radius
		(number) Radius of the circle.

   Return values:

	None

.. method:: Simulation.SetRegionEllipse(Layer, Material, Center, Angle, Halfwidths)

	Adds a (filled) ellipse of a specified material to an existing non-copy layer.
	The ellipse should not intersect any other patterning shapes, but may contain or be contained within other shapes.

   Usage::

	S.SetRegionEllipse(
		Layer = 'slab',
		Material = 'Vacuum',
		Center = (0,0),
		Angle = 30,            # in degrees
		Halfwidths = (0.3,0.4) # semi-axis lengths
	)

   Arguments:

	Layer
		(string) The name of the layer to pattern. This layer cannot be a copy of another layer.
	Material
		(string) The name of the material which fills the interior of the ellipse.
	Center
		(numeric table, length 2) x- and y-coordinates of the center of the ellipse relative to the center of the unit cell (the origin).
	Angle
		(number) The angle (in degrees) by which the x-axis of the ellipse should be rotated (CCW).
	Halfwidths
		(pair of numbers) The lengths of the semi-major axes of the ellipse. For an angle of 0, the first length is the semi-major axis in the x-direction, and the second length is the semi-major axis in the y-direction.

   Return values:

	None

.. method:: Simulation.SetRegionRectangle(Layer, Material, Center, Angle, Halfwidths)

	Adds a (filled) rectangle of a specified material to an existing non-copy layer.
	The rectangle should not intersect any other patterning shapes, but may contain or be contained within other shapes.

   Usage::

	S.SetRegionRectangle(
		Layer = 'slab',
		Material = 'Vacuum',
		Center = (0,0),
		Angle = 30,            # in degrees
		Halfwidths = (0.3,0.4)
	)

   Arguments:

	Layer
		(string) The name of the layer to pattern. This layer cannot be a copy of another layer.
	Material
		(string) The name of the material which fills the interior of the rectangle.
	Center
		(numeric table, length 2) x- and y-coordinates of the center of the rectangle relative to the center of the unit cell (the origin).
	Angle
		(number) The angle (in degrees) by which the x-axis of the rectangle should be rotated (CCW).
	Halfwidths
		(pair of numbers) The half-widths of the rectangle. For an angle of 0, the first length is half the width of the rectangle in the x-direction, and the second length is half the height in the y-direction.

   Return values:

	None


.. method:: Simulation.SetRegionPolygon(Layer, Material, Center, Angle, Vertices)

	Adds a (filled) polygon of a specified material to an existing non-copy layer.
	The polygon should not self-intersect nor intersect any other patterning shapes, but may contain or be contained within other shapes. The polygon must also be specified with positive orientation (the vertices circle CCW about an interior point).

   Usage::

	S.SetRegionPolygon(
		Layer = 'slab',
		Material = 'Vacuum',
		Center = (0,0),
		Angle = 10,            # in degrees
		Vertices = (
			(0,0),
			(0.2,0),
			(0.2,0.2),
			(0.1,0.2),
			(0.1,0.1),
			(0,0.1)
		)
	)

   Arguments:

	Layer
		(string) The name of the layer to pattern. This layer cannot be a copy of another layer.
	Material
		(string) The name of the material which fills the interior of the polygon.
	Center
		(numeric table, length 2) x- and y-coordinates of the center of the polygon relative to the center of the unit cell (the origin).
	Angle
		(number) The angle (in degrees) by which the polygon should be rotated (CCW).
	Vertices
		(tuple of pairs) The x- and y-coordinates of the vertices of the (unrotated) polygon, one after another, in counter-clockwise order.

   Return values:

	None


.. method:: Simulation.SetExcitationPlanewave(IncidenceAngles, sAmplitude=0, pAmplitude=0, Order=0])

	Sets the excitation to be a planewave incident upon the front (first layer specified) of the structure.
	If both tilt angles are specified to be zero, then the planewave is normally incident with the electric field polarized along the x-axis for the p-polarization.
	The phase of each polarization is defined at the origin (z = 0).

   Usage::

	S.SetExcitationPlanewave(
		IncidenceAngles=(
			10, # polar angle in [0,180)
			30  # azimuthal angle in [0,360)
		),
		sAmplitude = 0.707+0.707j,
		pAmplitude = 0.707-0.707j,
		Order = 0
	)

   Arguments:

	IncidenceAngles
		(pair of numbers) Of the form (phi,theta) with angles in degrees.
		``phi`` and ``theta`` give the spherical coordinate angles of the planewave k-vector.
		For zero angles, the k-vector is assumed to be (0, 0, kz), while the electric field is assumed to be (E0, 0, 0), and the magnetic field is in (0, H0, 0). The angle ``phi`` specifies first the angle by which the E,H,k frame should be rotated (CW) about the y-axis, and the angle ``theta`` specifies next the angle by which the E,H,k frame should be rotated (CCW) about the z-axis.
		Note the different directions of rotations for each angle.
	sAmplitude
		(complex number) The electric field amplitude of the s-polarizations of the planewave.
	pAmplitude
		(complex number) The electric field amplitude of the p-polarizations of the planewave.
	Order
		(integer) An optional positive integer specifying which order (mode index) to excite. Defaults to 0. Refer to :func:`GetBasisSet` for details.

   Return values:

	None

.. method:: Simulation.SetFrequency(freq)

	Sets the operating frequency of the system (and excitation).

   Usage::

	S.SetFrequency(1.2)

   Arguments:

	freq
		(complex number) The frequency of the excitation. This is not the angular frequency (the angular frequency is 2pi times of this).
		If a complex number is specified (typically for mode solving), the imaginary part should be negative for a physical (decaying in time) system.

   Return values:

	None

Outputs requiring no solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: Simulation.GetReciprocalLattice()

	Retrieves the reciprocal lattice basis vectors.
	The vectors lack the scaling by 2pi (multiply them by 2pi to obtain the true reciprocal lattice basis vectors).

   Usage::

	(Gu,Gv) = S.GetReciprocalLattice()

   Arguments:

	None.

   Return values:

	Gu, Gv
		The first and second reciprocal lattice basis vectors.
		Their relative geometric orientation is the same as the lattice specified in :func:`New`.
		Each vector is a tuple of length 2, holding the x- and y-coordinates of the vector.

.. method:: Simulation.GetEpsilon(x, y, z)

	Retrieves the dielectric constant at a particular point in the system by reconstructing the Fourier series using the G-vectors of the system.

	Note that this reconstruction is not representative of the actual dielectric constant profile used in simulations (such a notion is not meaningful). The reconstruction is created using the closed-form Fourier series coefficients of the specified patterning, summed over the terms comprising the G-vector list obtained from lattice truncation. This function exists to provide an intuitive sense for the spatial resolution of a particular G-vector truncation order.

   Usage::

	eps = S.GetEpsilon(0.1, 0.2, 0.3)

   Arguments:

	x, y, z
		(number) The coordinates of the point at which to retrieve the dielectric constant.

   Return values:

	eps
		The (usually complex) relative dielectric constant at the specified point.

.. method:: Simulation.OutputLayerPatternPostscript(Layer[, Filename])

	Outputs a list of PostScript commands to render the exact layer pattern description from the specified patterning commands. Assumes letter-sized paper.

   Usage::

	S.OutputLayerPatternPostscript(Layer = 'slab', Filename = 'out.ps')

   Arguments:

	Layer
		(string) Name of the layer whose pattern description should be output.
	Filename
		(string, optional) Filename to which the description should be output. If this argument is not provided, standard output is used.

   Return values:

	None.

Outputs requiring solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. method:: Simulation.OutputStructurePOVRay([Filename])

	Outputs a `POV-Ray <http://www.povray.org>`_ script that will render one unit cell of the structure in 3D. Materials named ``air`` or ``vacuum`` (case insensitive) will be completely transparent.
	Note that the output is not guaranteed to be correct or parsable by POV-Ray; it is only meant as a starting point to produce publication-ready figures.

   Usage::

	S.OutputStructurePOVRay(Filename = 'out.pov')

   Arguments:

	Filename
		(string, optional) Filename to which the structure should be output. If this argument is not provided, standard output is used.

   Return values:

	None.

.. method:: Simulation.GetBasisSet()

	Returns a tuple of reciprocal lattice coordinates of the Fourier series orders used.
	The coordinates are in the reciprocal lattice basis, and so they are integers.

   Usage::

	Glist = S.GetBasisSet()

   Arguments:

	None.

   Return values

	Glist
		A tuples of tuples of length 2 holding the pairs of integer recprical lattice coordinates.


.. method:: Simulation.GetAmplitudes(Layer, zOffset=0)

    Returns the raw mode amplitudes within a particular layer.
    For uniform (unpatterned) layers, the modes are simply the diffracted orders, and the indexing of the returned amplitudes corresponds to the value obtained from GetDiffractedOrder.
    The first value is guaranteed to be the straight transmitted or specularly reflected diffraction order.
    For patterned layers, there is typically no meaningful information in these amplitudes.

   Usage::

    (forw,back) = S.GetAmplitudes(Layer = 'substrate', zOffset = 0)

   Arguments:

    Layer
        (string) The name of the layer in which to obtain mode amplitudes.
    zOffset
        (number) The z-offset at which to obtain the mode amplitudes.

   Return values:

    forw,back
        Tuples of length 2*NumBasis containing the complex amplitudes of each forward and backward mode.
        In each tuple, the first NumBasis entries are for one of the polarizations, and the second
        NumBasis entries are for the other polarization.

.. method:: Simulation.GetPowerFlux(Layer, zOffset=0)

	Returns the integral of the power flux density over a unit cell surface normal to the z-direction.
	In other words, the z-component of the integrated Poynting flux is returned.

   Usage::

	(forw,back) = S.GetPowerFlux(Layer = 'substrate', zOffset = 0)

   Arguments:

	Layer
		(string) The name of the layer in which the integration surface lies.
	zOffset
		(number, optional) The z-offset of the integration surface from the beginning of the layer. This only matters for lossy layers.

   Return values:

	forw
		The forward component of the complex Poynting vector. Note that the result is not time averaged (no factor of 0.5 multiplied in). The forward component is defined as (E_total^* x H_forw + H_total^* x E_forw) / 2, where E_forw and H_forw are the fields reconstructed from only the forward propagating waveguide modes of the layer.
	back
		The backward component of the complex Poynting vector. Note that the result is not time averaged (no factor of 0.5 multiplied in). An analogous definition of the backward component of the Poynting vector follows from above.

.. method:: Simulation.GetPowerFluxByOrder(Layer, zOffset=0)

    Returns the integral of the Poynting flux density over a unit cell surface normal to the z-direction for each Fourier series order.
    In other words, the z-component of the Poynting flux for each order is returned.

   Usage::

    P = S.GetPowerFluxByOrder(Layer = 'substrate', zOffset = 0)

   Arguments:

    Layer
        (string) The name of the layer in which the integration surface lies.
    zOffset
        (number) The z-offset of the integration surface from the beginning of the layer. This only matters for lossy layers.

   Return values:

    P
        A tuple with length equal to the number of Fourier series orders used.
        Each entry of the is a pair (tuple) of forward and backward complex powers.
        These power quantities are described in the section for :func:`GetPowerFlux`.


.. method:: Simulation:GetStressTensorIntegral(Layer, zOffset=0)

    Returns the integral of the electromagnetic stress tensor over a unit cell surface normal to the z-direction.

   Usage::

    (Tx, Ty, Tz) = S.GetStressTensorIntegral(layer, offset)

   Arguments:

    Layer
        (string) The name of the layer in which the integration surface lies.
    zOffset
        (number) The z-offset of the integration surface from the beginning of the layer.

   Return values:

    Tx,Ty,Tz
        The real and imaginary parts of the x-, y-, and z-components of the stress tensor integrated over the specified surface, assuming a unit normal vector in the +z direction. Note that the result is not time averaged (no factor of 0.5 multiplied in).

.. method:: Simulation.GetLayerVolumeIntegral(Layer, Quantity)

    Returns the volume integral of a particular density over a unit cell throughout the entire thickness of a layer.

   Usage::

    I = S.GetLayerVolumeIntegral(Layer = 'slab', Quantity = 'U')

   Arguments:

	Layer
		(string) The name of the layer in which to integrate over.
	Quantity
		(string) The quantity to integrate. Currently, the choices are:
		
		U
			(epsilon*\|E\|^2 + \|H\|^2)
		E
			(epsilon*\|E\|^2)
		H
			(\|H\|^2)
		e
			(\|E\|^2)

   Return values:

    I
        The complex integral of the density throughout the volume of the layer's unit cell. Note that the result is not time averaged (no factor of 0.5 multiplied in).


.. method:: Simulation:GetLayerZIntegral(Layer, xy)

	Returns the line integral along z (depth direction) of the squared magnitudes of electric and magnetic field components (\|Ex\|^2, \|Ey\|^2, etc.) throughout the entire thickness of a layer.

   Usage::

    E,H = S.GetLayerZIntegral(Layer = 'slab', xy = (0.1, 0.2))

   Arguments:

    Layer
        (string) The name of the layer in which to integrate through.
    xy
        (pair of numbers) The in-plane coordinates at which to integrate.

   Return values:

    E,H
        Tuples of length 3 containing the x-, y-, and z-components of the integrated squared magnitude of the E or H field throughout the thickness of the layer. Note that the result is not time averaged (no factor of 0.5 multiplied in).


.. method:: Simulation.GetFields(x, y, z)

	Returns the electric and magnetic field at a particular point within the structure.

   Usage::

	(E,H) = S.GetFields(0.1, 0.2, 0.3)

   Arguments:

	x, y, z
		(number) The coordinates of the point at which to obtain the field.

   Return values:

	E,H
		Tuples of length 3 of the complex electric or magnetic field at the specified point. Note that the result is not time averaged (no factor of 0.5 multiplied in).
	
.. method:: Simulation.GetFieldsOnGrid(z, NumSamples, Format, BaseFilename)

	Returns the electric and magnetic fields on a regular grid over the unit cell (fundamental parallelogram) at a particular z coordinate.
	It is more efficient to use this function than :func:`GetFields`.

   Usage::

	E,H = S.GetFieldsOnGrid(z = 0.2, NumSamples=(100,100), Format = 'Array')
	S.GetFieldsOnGrid(z = 0.2, NumSamples=(100,100), Format = 'FileWrite', BaseFilename = 'field')
	S.GetFieldsOnGrid(z = 0.2, NumSamples=(100,100), Format = 'FileAppend', BaseFilename = 'field')

   Arguments:

	z
		(number) The z-coordinate of the plane on which to obtain the field.
	NumSamples
		(pair of integers) The number of sample points to use in each lattice vector direction.
	Format
		(string) Specifies the format of the output. Current choices are:
		
		Array
			Returns a pair of (tuple) arrays of dimension `nu` by `nv`, each element is a tuple of length 3, containing complex numbers of the E or H field components.
		FileWrite
			Outputs the field data to files, overwriting the files.
		FileAppend
			Outputs the field data to files, appending to the files. This is useful of volume fields are needed.
		
	BaseFilename
		(string) The base filename for file output. The outputs are named ``BaseFilename + '.E'`` and ``basefilename + '.H'``.

   Return values:

	E,H
		Only returned if format is 'Array'.
		Arrays (tuples) of dimension `nu` by `nv`, each element is a tuple of length 3, containing complex numbers of the E or H field components.

.. method:: Simulation.GetSMatrixDeterminant()

    Returns the determinant of the S-matrix (scattering matrix) of the entire structure.
    The determinant is an analytic function in the complex frequency plane and has poles at the complex modal frequencies of the system.

   Usage::

    (mant, base, expo) = S.GetSMatrixDeterminant()

   Arguments:

    None.

   Return values:

    mant
        The determinant typically causes overflow or underflow, so it is returned as a mantissa multiplying a base raised to an exponent. The value of the determinant is mant*base^expo.
    base
        The base of the determinant representation (see above).
    expo
        The exponent of the determinant representation (see above).


Options
^^^^^^^

.. method:: Simulation.SetOptions(...)

	Sets various options for a `Simulation object`_. The options are described below, and any option not specified is left unchanged.

   Usage::
   
	S.SetOptions( # these are the defaults
		Verbosity = 0,
		LatticeTruncation = 'Circular',
		DiscretizedEpsilon = False,
		DiscretizationResolution = 8,
		PolarizationDecomposition = False,
		PolarizationBasis = 'Default',
		LanczosSmoothing = False,
		SubpixelSmoothing = False,
		ConserveMemory = False
	)

   Arguments:
   
	Verbosity
		(integer) The larger this value, the more status output is generated. Valid values are in the range of 0-9, inclusive. A value of 0 disables all status output.
	LatticeTruncation
		(string) Sets the type of lattice truncation to use when selecting G-vectors. Can be one of the following values:
		
		Circular
			This is the default. The G-vectors are selected to have shortest length (by l2 norm).
		Parallelogramic
			Chooses the G-vectors within a parallelogram aligned with the reciprocal lattice basis. The number chosen will always be a perfect square of an odd number.

	DiscretizedEpsilon
		(boolean) Enables or disables the use of discretization in generating the Fourier coefficients of the in-plane epsilon profiles, instead of using values from closed-form equations. When enabled, the coefficients are obtained by FFT.
	DiscretizationResolution
		(integer) This option only has an effect when ``DiscretizedEpsilon`` or ``SubpixelSmoothing`` are used.
		This function sets the resolution of the FFT grid and vector field generated by ``PolarizationDecomposition``.
		The resolution is multiplied by the largest G-vector extent (integer lattice coordinate), and should be at least 2 to satisfy the Nyquist limit. It is best to use a number with small integer factors in order for the FFT to be computed efficiently. The size of each dimension of the FFT is obviously proportional to this value. The default is 8.
		See the :ref:`fmm-formulations-label` for details.
	PolarizationDecomposition
		(boolean) Enables or disables the use of proper in-plane Fourier factorization rules by decomposing fields into a polarization basis which conforms to the material boundaries.
		The polarization basis field is generated automatically by computing a quasi-harmonic vector field everywhere tangent to the layer pattern boundaries.
		This option is not guaranteed to work in the presence of tensor dielectric constants.
		Enabling this feature typically improves convergence with respect to the number of G-vectors.
		See the :ref:`fmm-formulations-label` for details.
	PolarizationBasis
		(string) Sets the method by which the polarization decomposition basis is generated.
		This option only has an effect when ``PolarizationDecomposition`` is set.
		See the :ref:`fmm-formulations-label` for details.
		Valid choices are:
		
		Default
			Uses a smooth tangent vector field with respect to layer patterning.
		Normal
			Uses a unit normal vector field with respect to layer patterning.
		Jones
			Uses a complex-valued Jones polarization vector field.
		
	LanczosSmoothing
		(boolean or dict) A boolean value enables or disables smoothing of the Fourier series representations of the layer dielectric constants using the Lanczos sigma factor (box filtering). This reduces the Gibbs phenomenon ringing in the real space reconstruction.
		Otherwise, specify a dictionary with keys ``Power`` (positive integer) and/or ``Width`` (positive number) to change the properties of the smoothing function.
	SubpixelSmoothing
		(boolean) Enables or disables the use of second-order accurate epsilon averaging rules within a pixel.
		The average epsilon within a pixel is computed using the fill factor of each material and the interface direction.
		Enabling this feature may improve convergence with respect to the number of G-vectors.
		See the :ref:`fmm-formulations-label` for details.
	ConserveMemory
		(boolean) Setting this option will prevent storage of certain intermediate results. This will save approximately 30% memory for non-trivial layers.
		The drawback is slower computation of any output quantities that require solutions.
		
   Return values:
   
	None.

Miscellaneous
^^^^^^^^^^^^^

.. method:: Simulation.Clone()

    Duplicates an existing `Simulation object`_, copying all materials, layers, and excitation information.
    No partial solution information is copied.

   Usage::

    S2 = S.Clone()

   Arguments:

    None.

   Return values:

    A copy of the `Simulation object`_.


.. |S4| replace:: S\ :sup:`4`
