import math
import FunctionSampler1D

sampler = FunctionSampler1D.new()

def f(x):
	return math.sin(1/(x*x+0.05))


sampler.add(-4, f(-4))
sampler.add(0, f(0))
sampler.add(4, f(4))

while not sampler.is_done():
	xvals = sampler.get_next()
	for x in xvals:
		y = f(x)
		sampler.add(x, y)

print('#', len(sampler))
for t in sampler:
	x, y, id = t;
	print(x, y)