import FunctionSampler2D
import numpy as np

def f(x,y):
	x2 = x*x
	if 0 == x and 0 == y:
		return 0
	else:
		return x2*y/(x2*x2+y*y)

sampler = FunctionSampler2D.new(
	(-2, -2, f(-2,-2)),
	( 2, -2, f( 2,-2)),
	( 2,  2, f( 2, 2)),
	(-2,  2, f(-2, 2))
)
sampler.options.max_principal_curvature = 30
sampler.options.xy_tol = 1e-3

print(-2, -2, f(-2,-2))
print( 2, -2, f( 2,-2))
print( 2,  2, f( 2, 2))
print(-2,  2, f(-2, 2))

for x in np.arange(-1,1, 0.2):
	for y in np.arange(-1,1, 0.2):
		z = f(x,y)
		sampler.add(x, y, z)
		#print(x, y, z)

while not sampler.is_done():
	xyvals = sampler.get_next(10)
	for xy in xyvals:
		x,y = xy
		z = f(x,y)
		sampler.add(x, y, z)
		#print(x, y, z)

for t in sampler:
	x,y,z,id = t
	print(x,y,z)
