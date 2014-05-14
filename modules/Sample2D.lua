FunctionSampler1D = require('FunctionSampler1D')
FunctionSampler2D = require('FunctionSampler2D')

function f1(x,y)
	local z = math.sin(x*y)
	print(x,y,z)
	return z
end
function flist(xylist)
	local zlist = {}
	for i,v in ipairs(xylist) do
		zlist[i] = f1(unpack(v))
	end
	return zlist
end

-- func(xylist) returns list of z
function Sample2D(xrange, yrange, func, npar, num_edge_samples, num_int_samples, sampler1)
	npar = npar or 1 -- default run 1 at a time
	num_edge_samples = num_edge_samples or 8
	num_int_samples = num_int_samples or 4
	sampler1 = sampler1 or FunctionSampler1D.New()

	local dx = xrange[2] - xrange[1]
	local dy = yrange[2] - yrange[1]

	local fnormalize = function(xy) return (xy[1]-xrange[1])/dx, (xy[2]-yrange[1])/dy end
	local funnormalize = function(uv) return uv[1]*dx+xrange[1], uv[2]*dy+yrange[1] end
	local fpaste = function(A,B)
		local ret = {}
		for i = 1,#A do
			ret[i] = { A[i][1], A[i][2], B[i] }
		end
		return ret
	end

	-- Generate initial corner samples
	local xy0list = {
		{ xrange[1], yrange[1] },
		{ xrange[2], yrange[1] },
		{ xrange[2], yrange[2] },
		{ xrange[1], yrange[2] }
	}
	-- Determine if we can afford to also throw in some interior points
	local n_interior = math.floor(math.sqrt(npar-4))
	for i = 1,n_interior do
		local ti = (i-0.5) / n_interior
		for j = 1,n_interior do
			local tj = (j-0.5) / n_interior
			table.insert(xy0list, { funnormalize{ti,tj} })
		end
	end

	local z0list = func(xy0list)
	local sampler = FunctionSampler2D.New(fpaste(xy0list, z0list))

	-- Sample the edges
	-- eitab entries: 1 index of starting point in xylist
	--                2 index of ending point in xylist
	--                3 dimension of variation (1 = x, 2 = y)
	local eitab = { {1,2,1},{2,3,2},{3,4,1},{4,1,2} }
	for iedge = 1,4 do
		sampler1:Clear()
		sampler1:Add(
			0, --xylist[eitab[iedge][1]][eitab[iedge][3]],
			z0list[eitab[iedge][1]]
		)
		sampler1:Add(
			1, --xylist[eitab[iedge][2]][eitab[iedge][3]],
			z0list[eitab[iedge][2]]
		)
		-- Add initial uniformly spaced edge samples
		xylist = {}
		tlist = {}
		for isamp = 1,num_edge_samples do
			xylist[isamp] = { xy0list[eitab[iedge][1]][1], xy0list[eitab[iedge][1]][2] }
			local t = isamp/(num_edge_samples+1)
			tlist[isamp] = t
			xylist[isamp][eitab[iedge][3]] = (
				(1-t) * xy0list[eitab[iedge][1]][eitab[iedge][3]] + t* xy0list[eitab[iedge][2]][eitab[iedge][3]]
			)
		end
		zlist = func(xylist)
		for i = 1,#xylist do
			sampler1:Add(tlist[i], zlist[i])
			local xx,yy = unpack(xylist[i])
			sampler:Add(xx, yy, zlist[i])
		end
		-- Adaptively sample edge
		while not sampler1:IsDone() do
			local tlist = sampler1:GetNext(npar)
			xylist = {}
			for isamp,t in ipairs(tlist) do
				xylist[isamp] = { xy0list[eitab[iedge][1]][1], xy0list[eitab[iedge][1]][2] }
				xylist[isamp][eitab[iedge][3]] = (
					(1-t) * xy0list[eitab[iedge][1]][eitab[iedge][3]] + t* xy0list[eitab[iedge][2]][eitab[iedge][3]]
				)
			end
			zlist = func(xylist)
			for i = 1,#xylist do
				sampler1:Add(tlist[i], zlist[i])
				local xx, yy = unpack(xylist[i])
				sampler:Add(xx, yy, zlist[i])
			end
		end
	end

	-- Sample the interior
	xylist = {}
	for i = 1,num_int_samples do
		local ti = i / (num_int_samples+1)
		for j = 1,num_int_samples do
			local tj = j / (num_int_samples+1)
			table.insert(xylist, { funnormalize{ti,tj} })
		end
	end
	local zlist = func(xylist)
	for i = 1,#xylist do
		local xx, yy = unpack(xylist[i])
		sampler:Add(xx, yy, zlist[i])
	end

	-- Adaptively sample interior
	while not sampler:IsDone() do
		local xylist = sampler:GetNext(npar)
		zlist = func(xylist)
		for i = 1,#xylist do
			local xx, yy = unpack(xylist[i])
			sampler:Add(xx, yy, zlist[i])
		end
	end
end

Sample2D({0,5}, {0,5}, flist)
