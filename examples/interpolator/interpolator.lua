interpolator = S4.NewInterpolator('cubic hermite spline', {
	{3.0, {14.2, 32}}, -- x, and list of y values
	{5.4, {4.6, 10}},
	{5.7, {42.7, 20}},
	{8.0, {35.2, 40}}
	})

for x = 3, 8, 0.1 do
	y1, y2 = interpolator:Get(x)
	print(x, y1, y2)
end
