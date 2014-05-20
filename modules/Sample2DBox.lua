inspect = require('inspect')
function Sample2DBox(argtable)
	if not argtable['Xlim'] then
		error('Must provide Xlim table: Xlim = { xmin, xmax }')
	end
	if not argtable['Ylim'] then
		error('Must provide Ylim table: Ylim = { ymin, ymax }')
	end
	if not argtable['UpdateFunc'] then
		error('Must provide UpdateFunc')
	end
	if not argtable['ComputeFunc'] then
		error('Must provide ComputeFunc')
	end

	local FunctionSampler1D = require('FunctionSampler1D')
	local FunctionSampler2D = require('FunctionSampler2D')

	local Xlim = argtable['Xlim']
	local Ylim = argtable['Ylim']
	local Ncores = argtable['Ncores'] or 1
	local Fupdate = argtable['UpdateFunc']
	local Fcompute = argtable['ComputeFunc']
	local Sampler1D = argtable['Sampler1D'] or FunctionSampler1D.New()
	local Nedge_initial = argtable['NumUniformEdgeSamples'] or 16
	local Nedge_max = argtable['MaxEdgeSamples'] or 128
	local Ninterior = argtable['NumUniformInteriorSamples'] or 100

	local dx = Xlim[2] - Xlim[1]
	local dy = Ylim[2] - Ylim[1]

	---------------------------------------------------------------------------
	---------- STAGE 0: Generate a list of cloned Simulation objects ----------
	---------------------------------------------------------------------------

	local Slist = { argtable['Simulation'] }
	if nil == Slist[0] then
		Slist = nil
	else
		for i = 2,Ncores do
			Slist[i] = Slist[0]:Clone()
		end
	end
	---------------------------------------------------------------------------
	-------- STAGE 1: Perform uniform and adative sampling along edges --------
	---------------------------------------------------------------------------

	local corners = { -- A list of the 4 corner vertices
		{ Xlim[1], Ylim[1] },
		{ Xlim[2], Ylim[1] },
		{ Xlim[2], Ylim[2] },
		{ Xlim[1], Ylim[2] }
	}
	local edge_endpoints = { -- indices into `corners' of the endpoints of 4 edges
		{ 1, 2 }, { 2, 3 }, { 3, 4 }, { 4,  1 }
	};
	local edge_varying_dim = { 1, 2, 1, 2 } -- which dimension of the edge is varying

	local edge_sample_points = {} -- list of size 4, each a list of edge samples
	local edge_sample_values = {}

	-- Perform uniform sampling on each edge
	-- Here, each edge includes the beginning vertex but excludes the ending vertex
	for iedge = 1,4 do
		edge_sample_points[iedge] = {}
		edge_sample_values[iedge] = {}
		local icorner1 = edge_endpoints[iedge][1] -- index of starting vertex of edge
		local icorner2 = edge_endpoints[iedge][2] -- index of ending vertex of edge
		local ivardim = edge_varying_dim[iedge] -- 1 or 2 for which dimension of the edge is varying

		-- Reset the 1D sampler.
		-- The 1D sampler will have normalized x values (the t parameter below)
		Sampler1D:Clear()

		-- A local function to translate a list of parameter t values to xy coords
		local t2xy_and_update = function(tlist)
			local xylist = {}
			for i = 1,#tlist do
				local t = tlist[i]
				xylist[i] = { corners[icorner1][1], corners[icorner1][2] }
				xylist[i][ivardim] = (
					(1-t) * corners[icorner1][ivardim] + t* corners[icorner2][ivardim]
				)
				if Slist then
					Fupdate(Slist[i], xylist[i][1], xylist[i][2])
				end
			end
			return xylist
		end

		-- Do the sampling in batches of `Ncores' points at a time
		local n_edge_batches = math.ceil(Nedge_initial/Ncores)
		for ibatch = 1,n_edge_batches do
			local tbatch = {} -- list of parameter values [0-1] along edge for current batch
			for icore = 1,Ncores do
				tbatch[icore] = (icore-1 + (ibatch-1)*Ncores) / (Ncores*n_edge_batches)
			end
			-- Translate from parameter values along edge to actual points along edge and update
			local xylist = t2xy_and_update(tbatch)
			local zlist = Fcompute(Slist, xylist)

			for i = 1,Ncores do
				table.insert(edge_sample_points[iedge], xylist[i])
				table.insert(edge_sample_values[iedge], zlist[i])
				Sampler1D:Add(tbatch[i], zlist[i])
			end
		end

		-- Finish up with adaptive sampling along the edge
		while not Sampler1D:IsDone() and #edge_sample_points[iedge] < Nedge_max do
			local tbatch = Sampler1D:GetNext(Ncores) -- Note: tbatch can have fewer than Ncores items
			local xylist = t2xy_and_update(tbatch)
			local zlist = Fcompute(Slist, xylist)
			for i = 1,Ncores do
				table.insert(edge_sample_points[iedge], xylist[i])
				table.insert(edge_sample_values[iedge], zlist[i])
				Sampler1D:Add(tbatch[i], zlist[i])
			end
		end
	end

	---------------------------------------------------------------------------
	----------- STAGE 2: Initialize a new FunctionSampler2D object ------------
	---------------------------------------------------------------------------
	local xy2uv = function(xy) return { (xy[1]-Xlim[1])/dx, (xy[2]-Ylim[1])/dy } end
	local uv2xy = function(uv) return { uv[1]*dx+Xlim[1], uv[2]*dy+Ylim[1] } end

	-- Add the four corners
	local uvlist = {
		xy2uv(edge_sample_points[1][1]),
		xy2uv(edge_sample_points[2][1]),
		xy2uv(edge_sample_points[3][1]),
		xy2uv(edge_sample_points[4][1])
	}
	local Sampler2D = FunctionSampler2D.New{
		{ uvlist[1][1], uvlist[1][2], edge_sample_values[1][1] },
		{ uvlist[2][1], uvlist[2][2], edge_sample_values[2][1] },
		{ uvlist[3][1], uvlist[3][2], edge_sample_values[3][1] },
		{ uvlist[4][1], uvlist[4][2], edge_sample_values[4][1] }
	}

	-- Add the rest of the edge points
	for iedge = 1,4 do
		for i = 2,#edge_sample_points[iedge] do
			local uv = xy2uv(edge_sample_points[iedge][i])
			Sampler2D:Add(
				uv[1], uv[2], edge_sample_values[iedge][i]
			)
		end
	end
	---------------------------------------------------------------------------
	----------------- STAGE 3: Uniformly sample the interior ------------------
	---------------------------------------------------------------------------

	-- Produce a list of uv locations
	local uvlist = {}
	local xylist = {}
	local n_int_side = math.ceil(math.sqrt(Ninterior))
	for i = 1,n_int_side do
		local ti = i/(n_int_side+1)
		for j = 1,n_int_side do
			local tj = j/(n_int_side+1)
			local uv = { ti, tj }
			table.insert(uvlist, uv)
			table.insert(xylist, uv2xy(uv))
		end
	end
	-- Compute at the uv locations in chunks of size Ncores
	for i = 1,#xylist,Ncores do
		local xybatch = {}
		local uvbatch = {}
		local Sbatch = {}
		for j = 1,Ncores do
			if i+j-1 <= #xylist then
				table.insert(uvbatch, uvlist[i+j-1])
				table.insert(xybatch, xylist[i+j-1])
				if Slist then
					Sbatch[j] = Slist[j]
					Fupdate(Sbatch[j], xybatch[j][1], xybatch[j][2])
				end
			end
		end
		local zbatch = Fcompute(Sbatch, xybatch)
		for j = 1,#xybatch do
			Sampler2D:Add(uvbatch[j][1], uvbatch[j][2], zbatch[j])
		end
	end

	---------------------------------------------------------------------------
	----------------- STAGE 4: Adaptively sample the interior -----------------
	---------------------------------------------------------------------------
	while not Sampler2D:IsDone() do
		local uvbatch = Sampler2D:GetNext(Ncores)
		local xybatch = {}
		local Sbatch = {}
		for i = 1,#uvbatch do
			xybatch[i] = uv2xy(uvbatch[i])
			if Slist then
				Sbatch[i] = Slist[i]
				Fupdate(Sbatch[i], xybatch[i][1], xybatch[i][2])
			end
		end
		zbatch = Fcompute(Sbatch, xybatch)
		for i = 1,#xybatch do
			Sampler2D:Add(uvbatch[i][1], uvbatch[i][2], zbatch[i])
		end
	end
end

function f(x,y)
	local z = math.sin(x*y)
	print(x,y,z)
	return z
end
function update(Slist)
end
function compute(Slist, xylist)
	local ret = {}
	for i = 1,#xylist do
		ret[i] = f(xylist[i][1], xylist[i][2])
	end
	return ret
end

Sample2DBox{
	Xlim = { 0, 5 },
	Ylim = { 0, 5 },
	Ncores = 4,
	UpdateFunc = update,
	ComputeFunc = compute
}
