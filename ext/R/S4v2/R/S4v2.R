# Object oriented structure as described here:
# https://rstudio-pubs-static.s3.amazonaws.com/150296_904158e070594471864e89c10c0d14f9.html
S4_Simulation <- function(lattice, bases){
	# argument checking
	if(is.complex(lattice) || !(
		is.numeric(lattice) &&
		(1 == length(lattice) || (is.matrix(lattice) && 2 == nrow(lattice) && 2 == ncol(lattice)))
	)){
		stop("lattice must be numeric and either a real scalar or 2x2 matrix")
	}
	if(!(
		((1 == length(bases)) && (bases > 0)) ||
		(length(bases) > 1) ||
		(is.matrix(bases) && (2 == nrow(bases)))
	)){
		stop("bases must be a positive integer, an integer vector, or a 2xN matrix")
	}
	storage.mode(lattice) <- "double"
	storage.mode(bases) <- "integer"
	S <- .Call("R_S4v2_Simulation_New", lattice, bases)


	add_material <- function(name = NULL, epsilon = 1){
		if(!(is.null(name) | is.character(name))){
			stop("name must be null or a string")
		}
		if(!(
			(is.numeric(epsilon) || is.complex(epsilon)) &&
			(1 == length(epsilon) || (is.matrix(epsilon) && 3 == nrow(epsilon) && 3 == ncol(epsilon)))
		)){
			stop("epsilon must be real or complex, and either scalar or a 3x3 matrix")
		}
		S <- S
		M <- .Call("R_S4v2_Simulation_AddMaterial", S, name, epsilon)

		# get_epsilon <- function(){}
		# get_name <- function(){}

		structure(class = "S4_Material", environment())
	}

	add_layer <- function(name = NULL, thickness = 0, copy = NULL, material = NULL){
		if(!(is.null(name) | is.character(name))){
			stop("name must be null or a string")
		}
		if(!(
			is.numeric(thickness) && thickness >= 0
		)){
			stop("thickness must be a non-negative number")
		}
		if((is.null(material) && is.null(copy)) || (!is.null(material) && !is.null(copy))){
			stop("Only one of material or copy may be non-NULL")
		}
		if(!is.null(material) && class(material) != "S4_Material"){
			stop("material must be a valid S4_Material")
		}
		if(!is.null(copy) && class(copy) != "S4_Layer"){
			stop("copy must be a valid S4_Layer")
		}

		storage.mode(thickness) <- "double"
		S <- S
		L <- .Call("R_S4v2_Simulation_AddLayer", S, name, thickness, copy$L, material$M)

		set_region <- function(
			material,
			shape = "interval",
			center = c(0,0),
			halfwidths = NULL,
			radius = NULL,
			angle_frac = 0,
			vertices = NULL
		){
			if(!is.null(material) && class(material) != "S4_Material"){
				stop("material must be a valid S4_Material")
			}
			if(
				"rectangle" == shape ||
				"interval" == shape ||
				"ellipse" == shape
			){
				if(!(
					!is.null(halfwidths) &&
					is.null(radius) && is.null(vertices)
				)){
					stop("Must specify halfwidths for shape")
				}
			}
			if("circle" == shape){
				if(!(
					!is.null(radius) &&
					is.null(halfwidths) && is.null(vertices)
				)){
					stop("Must specify radius for 'circle' shape")
				}
				halfwidths = c(radius, radius)
			}
			if("polygon" == shape){
				if(!(
					!is.null(vertices) &&
					is.null(halfwidths) && is.null(radius)
				)){
					stop("Must specify vertices for 'polygon' shape")
				}
			}
			storage.mode(center) <- "double"
			storage.mode(halfwidths) <- "double"
			storage.mode(angle_frac) <- "double"
			storage.mode(vertices) <- "double"
			if("polygon" == shape){
				.Call("R_S4v2_Layer_SetRegionVertices", S, L, material$M, shape, center, vertices, angle_frac);
			}else{
				.Call("R_S4v2_Layer_SetRegionHalfwidths", S, L, material$M, shape, center, halfwidths, angle_frac);
			}
			invisible()
		}

		get_power_flux <- function(offset = 0){
			storage.mode(offset) <- "double"
			.Call("R_S4v2_Layer_GetPowerFlux", S, L, offset);
		}
		get_waves <- function(){
			waves <- .Call("R_S4v2_Layer_GetWaves", S, L);
			rownames(waves) <- c("kx", "ky", "kzr", "kzi", "ux", "uy", "uz", "cur", "cui", "cvr", "cvi")
			return(waves)
		}

		structure(class = "S4_Layer", environment())
	}

	set_frequency <- function(freq){
		if(!(
			(is.numeric(freq) && freq > 0) ||
			(is.complex(freq) && Re(freq) > 0 && Im(freq) >= 0)
		)){
			stop("Frequency must be a positive number or a complex number in the first quadrant")
		}
		.Call("R_S4v2_Simulation_SetFrequency", S, as.complex(freq));
		invisible()
	}

	get_epsilon <- function(p = c(0,0,0)){
		if(!(
			is.numeric(p) && 3 == length(p)
		)){
			stop("position must be a real 3-vector")
		}
		.Call("R_S4v2_Simulation_GetEpsilon", S, as.double(p));
	}

	excitation_planewave <- function(k = c(0,0,1), u = c(1,0,0), cu = 0, cv = 0){
		if(!(
			is.numeric(k) && 3 == length(k)
		)){
			stop("k must be a real 3-vector")
		}
		if(!(
			is.numeric(u) && 3 == length(u)
		)){
			stop("k must be a real 3-vector")
		}
		if(!(
			(is.numeric(cu) || is.complex(cu)) && 1 == length(cu)
		)){
			stop("cu must be a scalar real or complex number")
		}
		if(!(
			(is.numeric(cv) || is.complex(cv)) && 1 == length(cv)
		)){
			stop("cu must be a scalar real or complex number")
		}
		.Call("R_S4v2_Simulation_ExcitationPlanewave", S, as.double(k), as.double(u), as.complex(cu), as.complex(cv));
		invisible()
	}
	structure(class = "S4_Simulation", environment())
}

# Primitive functions

S4_Simulation_New <- function(lattice, bases){
	storage.mode(lattice) <- "double"
	storage.mode(bases) <- "integer"
	.Call("R_S4v2_Simulation_New", lattice, bases);
}
S4_Simulation_AddMaterial <- function(S, name = NULL, epsilon = 1){
	.Call("R_S4v2_Simulation_AddMaterial", S, name, epsilon);
}
S4_Simulation_AddLayer <- function(S, name = NULL, thickness = 0, copy = NULL, material = NULL){
	.Call("R_S4v2_Simulation_AddLayer", S, name, thickness, copy, material);
}
S4_Layer_SetRegionHalfwidths <- function(S, layer = NULL, material = NULL, shape = "rectangle", center = c(0,0), halfwidths = c(0,0), angle_frac = 0){
	.Call("R_S4v2_Layer_SetRegionHalfwidths", S, layer, material, shape, center, halfwidths, angle_frac);
}
S4_Layer_SetRegionVertices <- function(S, layer = NULL, material = NULL, shape = "rectangle", center = c(0,0), vertices = matrix(c(0,0,1,0,0,1),nrow=2), angle_frac = 0){
	.Call("S4_Layer_SetRegionVertices", S, layer, material, shape, center, vertices, angle_frac);
}
S4_Simulation_SetFrequency <- function(S, freq){
	.Call("R_S4v2_Simulation_SetFrequency", S, as.complex(freq));
}
S4_Simulation_ExcitationPlanewave <- function(S, k = c(0,0,1), u = c(1,0,0), cu = 0, cv = 0){
	.Call("R_S4v2_Simulation_ExcitationPlanewave", S, k, u, as.complex(cu), as.complex(cv));
}
S4_Simulation_GetEpsilon <- function(S, p = c(0,0,0)){
	.Call("R_S4v2_Simulation_GetEpsilon", S, p);
}

S4_Layer_GetPowerFlux <- function(S, layer = NULL, offset = 0){
	.Call("R_S4v2_Layer_GetPowerFlux", S, layer, offset);
}
S4_Layer_GetWaves <- function(S, layer = NULL){
	.Call("R_S4v2_Layer_GetWaves", S, layer);
}
