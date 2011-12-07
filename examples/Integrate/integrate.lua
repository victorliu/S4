function f(x,y)
	return math.cos(x*x)*math.sin(y*y)
end


print(S4.Integrate(f, {0,1}, {0,1}, {
	MaxEval = 1000,
	AbsoluteError = 0,
	RelativeError = 1e-5
	}))
