# S4v2 Lua API

This file describes the new Lua interface to S4.
Not all functionality of the original API is exposed yet.

## Top level functions

### NewSimulation

Creates a new simulation object with a specified lattice and basis set.
The lattice may be specified as a single number (a 1D lattice), or a pair of 2D lattice vectors.
The basis set may be specified by a single number (the number of basis functions) or an explicit list
of basis waves (single integers for 1D, or pairs of integers for 2D).

Example:

    S = S4v2.NewSimulation{ lattice = 0.500, bases = 7 }
    S = S4v2.NewSimulation{ lattice = {{1, 0}, {0, 1}}, bases = 5 }
    S = S4v2.NewSimulation{ lattice = 0.5, bases = {0, -1, 1, -2, 2, -3, 3} }
    S = S4v2.NewSimulation{ lattice = {{1, 0}, {0, 1}}, bases = { {0,0}, {1,0}, {-1,0}, {0,1}, {0,-1} } }

## Simulation object methods:

* AddMaterial
* AddLayer
* Clone
* SetLattice/GetLattice
* SetBases/GetBases
* SetFrequency/GetFrequency
* GetMaterial/GetLayer
* ExcitationPlanewave
* GetEpsilon

## Material object methods:

* SetName/GetName
* SetEpsilon/GetEpsilon

## Layer object methods:

* SetName/GetName
* SetThickness/GetThickness
* ClearRegions
* SetRegion
* GetPowerFlux
* GetWaves

### Simulation: AddMaterial

Add a new material to the simulation material library and return a material ID handle that may be used to refer to it later.
A material may optionally be referred to by name as well. The dielectric constant may be specified as a scalar real or complex number,
or as a tensor quantity. Currently, only z-decoupled tensors are supported (xz, yz, zx, zy components must be zero).

Example:

    M = S:AddMaterial{ name = 'Glass', epsilon = 1.5^2 }
    M = S:AddMaterial{ name = 'Glass', epsilon = { 1.5^2, 0 } }
    mGlass = S:AddMaterial{ epsilon = {
        { { 1.5^2, 0 }, {0,0}, {0,0} },
        { {0,0}, { 1.5^2, 0 }, {0,0} },
        { {0,0}, {0,0}, { 1.5^2, 0 } }
    }

### Simulation: AddLayer

Add a new layer to the simulation. Layers are added in order of increasing z-coordinate.
The first and last layers are typically specified to have zero thickness (but not required).
A layer ID handle is returned that may be used to refer to it later.
The material may be specified by either a material ID handle or by a material's name if it was specified.

Alternatively, if a layer is identical in composition (all attributes except thickness), then it may be specified as a copy of a previously defined layer.
In this case, the previously defined layer may be referred to by a layer ID handle or its layer name.

Example:

    L  = S:AddLayer{ name = 'Substrate', material = mGlass, thickness = 0.200 }
    L2 = S:AddLayer{ name = 'Superstrate', thickness = 0.100, copy = 'Substrate' }

### Simulation: Clone

Create a copy of a simulation object.

Example:

    S2 = S:Clone()

### Simulation: SetLattice/GetLattice

Set or get the real-space lattice.

Example:

    S:SetLattice{ {1, 0}, {0, 1} }
    L = S:GetLattice()

### Simulation: SetFrequency/GetFrequency

Set or get the frequency of operation.

Example:

    S:SetFrequency(1 / wavelength)
    wavelength = 1 / S:GetFrequency()

### Simulation: GetMaterial/GetLayer

Get a material or layer handle ID by name.

Example:

    M = S:GetMaterial('Glass')
    L = S:GetLayer('Substrate')

### Simulation: ExcitationPlanewave

Set a planewave excitation source.
A planewave is specified by its incident direction vector (k, not necessarily normalized),
its electric field polarization vector (u, orthogonal to k, not necessarily normalized),
and its complex amplitudes (cu and cv) along u and v (the cross product of u with v is in the k direction).
The complex amplitudes are specified by real and imaginary parts.

Example:

    S:ExcitationPlanewave{ k = {0, 0, 1}, u = {0, 1, 0}, cu = {1, 0}, cv = {0, 0} }

### Simulation: GetEpsilon

Returns the dielectric constant at a point in space.
Currently, the dielectric constant is assumed to be a scalar complex number.

Example:

    eps_real, eps_imag = S:GetEpsilon{x, y, z}

### Material: SetName/GetName

Set or get the name of an existing material.

Example:

    M:SetName('Glass')
    print('Glass' == M:GetName())

### Material: SetEpsilon/GetEpsilon

Sets or gets the dielectric constant of an existing material.
The dielectric constant may be specified by a scalar real or complex number, or as a 3x3 real or complex tensor.
The dielectric constant is always returned as a 3x3 complex tensor.

Example:

    M:SetEpsilon(1.5^2)
    eps = M:GetEpsilon()
    for row = 1,3 do
        for col = 1,3 do
            print(row, col, eps[row][col][1], eps[row][col][2])
        end
    end

### Layer: SetName/GetName

Set or get the name of an existing layer.

Example:

    L:SetName('Substrate')
    print('Substrate' == L:GetName())

### Layer: SetThickness/GetThickness

Set or get the thickness of an existing layer.

Example:

    L:SetThickness(0.250)
    print(0.250 == L:GetThickness())

### Layer: ClearRegions

Clear all patterned regions of a layer.
The layer is reset to a homogeneous region of its specified material.

Example:

    L:ClearRegions()

### Layer: SetRegion

Adds a patterned region to a layer. The region is specified by a shape and a material.
A variety of shapes are available: circle, ellipse, rectangle, or polygon.
The shape may optionally be centered away from the origin or rotated about its center by some angle.

Example:

    L:SetRegion{ material = mGlass, shape = 'circle', center = {0, 0}, radius = 0.2 }
    L:SetRegion{ material = mGlass, shape = 'rectangle', center = {0, 0}, halfwidths = {0.1, 0.2}, angle_degrees = 45 }
    L:SetRegion{ material = mGlass, shape = 'ellipse', center = {0, 0}, halfwidths = {0.1, 0.2}, angle_radians = math.atan(0.2) }
    L:SetRegion{ material = mGlass, shape = 'polygon', vertices = { {0,0}, {1,0}, {0,1} }, center = {-0.2, -0.2}, angle_degrees = 30 }

### Layer: GetPowerFlux

Returns the forward and backwards power flux within the layer.
This method solves for the scattered fields of the specified structure given the specified excitation.
The forward and backwards powers are returned by real and imaginary parts in the case that there is reactive power in a layer.

Example:

    forw_real, back_real, forw_imag, back_imag = L:GetPowerFlux()

### Layer: GetWaves

Returns the forwards and backwards propagating waves within the layer.
This method solves for the scattered fields of the specified structure given the specified excitation.
The forwards and backwards waves have physical only for unpatterned layers, in which they correspond to the Rayleigh (planewave) modes.
Similar to `ExcitationPlanewave`, each wave is specified by its direction (k), its electric field polarization vector direction (u, orthogonal to k), and its complex amplitudes along u and v.
Different to `ExcitationPlanewave`, the k-vector has a 4-th elementh, the imaginary part of kz which will be nonzero for evanescent waves.

Example:

    waves = L:GetWaves()
    for i,wave = ipairs(waves)
        print(i, wave.k[1], wave.k[2], wave.k[3], wave.k[4], wave.u[1], wave.u[2], wave.u[3], wave.cu[1], wave.cu[2], wave.cv[1], wave.cv[2])
    end
