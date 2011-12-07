function f(x)
	x = x/0.1;
	return math.exp(-x*x)
end

function f2(x)
	return math.sin(1/x)
end

--[[
sampler = S4.NewSpectrumSampler(-1,1, {
	InitialNumPoints = 33,
	RangeThreshold = 0.001,
	MaxBend = math.cos(math.rad(10)),
	MinimumSpacing = 1e-6
	})
while not sampler:IsDone() do
	x = sampler:GetFrequency()
	y = f2(x)
	--print(x,y)
	sampler:SubmitResult(y)
end

spectrum = sampler:GetSpectrum()
for i,xy in ipairs(spectrum) do
	print(xy[1],xy[2])
end


]]


sampler = S4.NewSpectrumSampler(-1,1, {
	InitialNumPoints = 33,
	RangeThreshold = 0.001,
	MaxBend = math.cos(math.rad(10)),
	MinimumSpacing = 1e-6,
	Parallelize = true
	})
while not sampler:IsDone() do
	local xvec = sampler:GetFrequencies()
	local y = {}
	for i,x in ipairs(xvec) do
		y[i] = f2(x)
		--print(x,y[i])
	end
	--print(x,y)
	sampler:SubmitResults(y)
end

spectrum = sampler:GetSpectrum()
for i,xy in ipairs(spectrum) do
	print(xy[1],xy[2])
end
