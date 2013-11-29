FunctionSampler1D = require('FunctionSampler1D')

function f(x)
	return math.sin(1/(x*x+0.05))
end

sampler = FunctionSampler1D.New()

sampler:Add(-4, f(-4))
sampler:Add( 0, f( 0))
sampler:Add( 4, f( 4))

while not sampler:IsDone() do
	local xvals = sampler:GetNext()
	for ix,x in ipairs(xvals) do
		local y = f(x)
		sampler:Add(x, y)
		print(x, y)
	end
end

