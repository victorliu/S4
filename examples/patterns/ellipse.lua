-- Run as
--   S4 ellipse.lua >ellipse_eps_cell.txt 2>ellipse_eps_real.txt

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(80)

-- Material definition
S:AddMaterial('Silicon', {12,0}) -- real and imag parts
S:AddMaterial('Vacuum', {1,0})

-- Structure definition
S:AddLayer('AirAbove',  -- layer name
           0,           -- thickness
           'Vacuum')    -- background material
S:AddLayer('Slab', 0.5, 'Silicon')
S:SetLayerPatternEllipse('Slab',     -- which layer to alter
                         'Vacuum',   -- material in ellipse
	                     {0.2,0},    -- center
	                     45,         -- tilt angle (degrees)
	                     {0.3, 0.2}) -- half-widths (semi-axes)
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

-- Dump the epsilon reconstruction in the unit cell to stdout
S:OutputLayerPatternRealization('Slab', 64, 64)

-- Dump the epsilon reconstruction in real space to stderr
for x=-1.5,1.5,0.02 do
	for y=-1.5,1.5,0.02 do
		er,ei = S:GetEpsilon({x,y,0.25}) -- returns real and imag parts
		io.stderr:write(x .. '\t' .. y .. '\t' .. er .. '\t' .. ei .. '\n')
	end
	io.stderr:write('\n')
end
