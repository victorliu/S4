pcall(loadstring(S4.arg))

if not NG then NG = 10 end

nset = 3
Svec = {}

-- Bottom pane of Fig. 12 in
-- Shanhui Fan and J. D. Joannopoulos,
-- "Analysis of guided resonances in photonic crystal slabs",
-- Phys. Rev. B, Vol. 65, 235112

--function make_sims()
	local S = S4.NewSimulation()
	S:SetLattice({1,0}, {0,1})
	S:SetNumG(NG)
	S:AddMaterial("Silicon", {12,0}) -- real and imag parts
	S:AddMaterial("Vacuum", {1,0})

	S:AddLayer('AirAbove', 0 , 'Vacuum')
	S:AddLayer('Slab', 0.5, 'Silicon')
	S:SetLayerPatternCircle('Slab', 'Vacuum', {0,0}, 0.2)
	S:AddLayerCopy('AirBelow', 0, 'AirAbove')

	S:SetExcitationPlanewave(
		{0,0}, -- incidence angles
		{1,0}, -- s-polarization amplitude and phase (in degrees)
		{0,0}) -- p-polarization amplitude and phase


	if form == 'ref' then
		S:UseDiscretizedEpsilon(false)
	elseif form == 'fft' then
		S:UseDiscretizedEpsilon(true)
		S:SetResolution(8)
	elseif form == 'sps' then
		S:UseDiscretizedEpsilon(true)
		S:UseSubpixelSmoothing()
		S:SetResolution(8)
	elseif form == 'nv' then
		S:UsePolarizationDecomposition()
		S:UseNormalVectorBasis()
		S:SetResolution(8)
	elseif form == 'cpx' then
		S:UsePolarizationDecomposition()
		S:UseJonesVectorBasis()
		S:SetResolution(8)
	elseif form == 'pol' then
		S:UsePolarizationDecomposition()
		S:SetResolution(8)
	end


	local Svec = {}
	for i = 1,nset do
		Svec[i] = S:Clone()
	end
--end

function solve_freqs(fvec)
	Svec_short = {}
	for i,f in ipairs(fvec) do
		Svec[i]:SetFrequency(f)
		Svec_short[i] = Svec[i]
	end
	S4.SolveInParallel('AirBelow', unpack(Svec_short))
	
	local tvec = {}
	for i,f in ipairs(fvec) do
		tvec[i] = Svec[i]:GetPoyntingFlux('AirBelow', 0)
	end

	return tvec
end

function solve_freq(freq)
	return solve_freqs({freq})[1]
end

function solve_set(fvec)
	local tvec = {}
	for i,v in SetsOfN(fvec,nset) do
		local nv = #v
		local t = solve_freqs(v);
		for j = 1,nset do
			tvec[i-nv+j] = t[j]
		end
	end
	return tvec
end


function iterN(a,i,N)
	if a[i+1] then
		local v = {}
		for j = 1,N do
			if a[i+1] then
				v[j] = a[i+1]
				i = i+1
			else
				break;
			end
		end
		return i,v
	end
end
function SetsOfN(a,N)
	return function (a,i) return iterN(a,i,N) end, a, 0
end
--[[
for freq=0.25,0.6,0.003 do
--	print(freq .. '\t' .. solve_freq(freq))
end
]]


sampler = S4.NewSpectrumSampler(0.5, 0.55, {
	InitialNumPoints = 33,
	RangeThreshold = 0.001,
	MaxBend = math.cos(math.rad(10)),
	MinimumSpacing = 1e-6,
	Parallelize = true
	})
while not sampler:IsDone() do
	local xvec = sampler:GetFrequencies()
	local y = solve_set(xvec)
	
	for i,x in ipairs(xvec) do
		print(x,y[i])
	end
	sampler:SubmitResults(y)
end

--print(solve_freq(0.54392194747925))
