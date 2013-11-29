FunctionSampler2D = require('FunctionSampler2D')

function f(x,y)
	local x2 = x*x
	if 0 == x and 0 == y then
		return 0
	else
		return x2*y/(x2*x2+y*y)
	end
end

sampler = FunctionSampler2D.New{
	{-2, -2, f(-2,-2)},
	{ 2, -2, f( 2,-2)},
	{ 2,  2, f( 2, 2)},
	{-2,  2, f(-2, 2)}
}
sampler:SetParameters{
	MaxPrincipalCurvature = 30
}

print(-2, -2, f(-2,-2))
print( 2, -2, f( 2,-2))
print( 2,  2, f( 2, 2))
print(-2,  2, f(-2, 2))

for x = -1,1,0.2 do
	for y = -1,1,0.2 do
		local z = f(x,y)
		sampler:Add(x, y, z)
		print(x, y, z)
	end
end

while not sampler:IsDone() do
	local xyvals = sampler:GetNext(10)
	for ixy,xy in ipairs(xyvals) do
		local x = xy[1]
		local y = xy[2]
		local z = f(x,y)
		sampler:Add(x, y, z)
		print(x, y, z)
	end
end

