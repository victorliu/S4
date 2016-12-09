#include "mex.h"
#include "S4.h"

//#define OUTPUT_MEX_TRACE

#ifdef OUTPUT_MEX_TRACE
#define MEX_TRACE(...) mexPrintf(__VA_ARGS__)
#else
#define MEX_TRACE(...)
#endif

struct S4ptr{
	S4_Simulation *S;
	int is1d;
};

static int nS_alloc = 0;
static struct S4ptr *SS = NULL;

static int S_add(S4_Simulation *S){
	int i;
	if(0 == nS_alloc){
		nS_alloc = 8;
		SS = (struct S4ptr*)malloc(sizeof(struct S4ptr) * nS_alloc);
		SS[0].S = S;
		for(i = 1; i < nS_alloc; ++i){
			SS[i].S = NULL;
		}
		return 0;
	}
	for(i = 0; i < nS_alloc; ++i){
		if(NULL == SS[i].S){
			SS[i].S = S;
			return i;
		}
	}
	i = nS_alloc;
	nS_alloc *= 2;
	SS = (struct S4ptr*)realloc(S, sizeof(struct S4ptr) * nS_alloc);
	SS[i].S = S;
	return i;
}
static int S_del(int i){
	if(i < 0 || i >= nS_alloc){ return -1; }
	SS[i].S = NULL;
	return i;
}
static int S_check(const mxArray *prhs[]){
	int i;
	if(mxINT32_CLASS != mxGetClassID(prhs[1]) || !mxIsScalar(prhs[1])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected simulation handle argument after command.");
	}
	i = *(int*)(mxGetData(prhs[1]));
	if(i < 0 || i >= nS_alloc){
		mexErrMsgIdAndTxt("S4:invalidInput", "Invalid simulation handle argument after command.");
	}
	return i;
}
static S4_Simulation* Sptr_check(const mxArray *prhs[]){
	int i;
	if(mxINT32_CLASS != mxGetClassID(prhs[1]) || !mxIsScalar(prhs[1])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected simulation handle argument after command.");
	}
	i = *(int*)(mxGetData(prhs[1]));
	if(i < 0 || i >= nS_alloc){
		mexErrMsgIdAndTxt("S4:invalidInput", "Invalid simulation handle argument after command.");
	}
	if(NULL != SS[i].S){
		return SS[i].S;
	}else{
		mexErrMsgIdAndTxt("S4:invalidInput", "Invalid simulation handle argument after command.");
	}
}
static int ID_check(int arg, const mxArray *prhs[]){
	int i;
	if(mxINT32_CLASS != mxGetClassID(prhs[arg]) || !mxIsScalar(prhs[arg])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected handle argument after command.");
	}
	i = *(int*)(mxGetData(prhs[arg]));
	return i;
}
static double real_check(int arg, const mxArray *prhs[], const char *argname){
	if(!mxIsNumeric(prhs[arg]) || !mxIsScalar(prhs[arg]) || mxIsComplex(prhs[arg])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected number.");
	}
	return *(mxGetPr(prhs[arg]));
}

static void complex_check(int arg, const mxArray *prhs[], const char *argname, double *ret){
	if(!mxIsNumeric(prhs[arg]) || !mxIsScalar(prhs[arg])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected complex number.");
	}
	ret[0] = *(mxGetPr(prhs[arg]));
	if(mxIsComplex(prhs[arg])){
		ret[1] = *(mxGetPi(prhs[arg]));
	}
}

static void vec3_check(int arg, const mxArray *prhs[], const char *argname, double *ret){
	int i;
	const double *p;
	if(!mxIsNumeric(prhs[arg]) || 3 != mxGetNumberOfElements(prhs[arg])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected 3-vector.");
	}
	p = mxGetPr(prhs[arg]);
	for(i = 0; i < 3; ++i){
		ret[i] = p[i];
	}
}

static void vec2_check(int arg, const mxArray *prhs[], const char *argname, double *ret){
	int i;
	const double *p;
	if(!mxIsNumeric(prhs[arg]) || 2 != mxGetNumberOfElements(prhs[arg])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected 2-vector.");
	}
	p = mxGetPr(prhs[arg]);
	for(i = 0; i < 2; ++i){
		ret[i] = p[i];
	}
}

static void ExitFcn(void){
	int i;
	if(NULL == SS || nS_alloc <= 0){ return; }
	for(i = 0; i < nS_alloc; ++i){
		if(NULL != SS[i].S){
			S4_Simulation_Destroy(SS[i].S);
		}
	}
	free(SS);
}

static void S4mex_init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
}

static void S4mex_Simulation_New(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_real Lr[4] = { 0, 0, 0, 0 };
	int nbases = 1;
	int is1d = 0;
	
	if(3 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected two additional arguments to 'new'.");
	}
	
	if(!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || !((is1d = mxIsScalar(prhs[1])) || (2 == mxGetM(prhs[1]) && 2 == mxGetN(prhs[1])))){
		mexErrMsgIdAndTxt("S4:invalidInput", "First additional argument (lattice) must be scalar or a 2x2 matrix.");
	}
	if(is1d){
		Lr[0] = mxGetScalar(prhs[1]);
	}else{
		double *p = mxGetPr(prhs[1]);
		Lr[0] = p[0];
		Lr[1] = p[1];
		Lr[2] = p[2];
		Lr[3] = p[3];
		if(0 == Lr[1] && 0 == Lr[2] && 0 == Lr[3]){
			is1d = 1;
		}
	}
	
	if(!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || !mxIsScalar(prhs[2])){
		mexErrMsgIdAndTxt("S4:invalidInput", "Second additional argument (#bases) must be a positive real integer.");
	}
	if(mxIsDouble(prhs[2])){
		double *p = mxGetPr(prhs[2]);
		nbases = (int)p[0];
	}else{
		void *p = mxGetData(prhs[2]);
		switch(mxGetClassID(prhs[2])){
		case mxSINGLE_CLASS:
			{ float *pp = (float*)p; nbases = (int)pp[0]; break; }
        case mxINT8_CLASS:
			{ char *pp = (char*)p; nbases = (int)pp[0]; break; }
        case mxUINT8_CLASS:
			{ unsigned char *pp = (unsigned char*)p; nbases = (int)pp[0]; break; }
        case mxINT16_CLASS:
			{ short *pp = (short*)p; nbases = (int)pp[0]; break; }
        case mxUINT16_CLASS:
			{ unsigned short *pp = (unsigned short*)p; nbases = (int)pp[0]; break; }
        case mxINT32_CLASS:
			{ int *pp = (int*)p; nbases = (int)pp[0]; break; }
        case mxUINT32_CLASS:
			{ unsigned int *pp = (unsigned int*)p; nbases = (int)pp[0]; break; }
        case mxINT64_CLASS:
			{ long long *pp = (long long*)p; nbases = (int)pp[0]; break; }
        case mxUINT64_CLASS:
			{ unsigned long long *pp = (unsigned long long*)p; nbases = (int)pp[0]; break; }
		default:
			mexErrMsgIdAndTxt("S4:invalidInput", "Second additional argument (#bases) is an unknown class.");
			break;
		}
	}
	if(nbases < 1){
		mexErrMsgIdAndTxt("S4:invalidInput", "Second additional argument (#bases) must be a positive integer.");
	}
	
	MEX_TRACE("> S4mex_Simulation_New(lattice = [%f, %f; %f, %f], bases = %d)\n", Lr[0], Lr[1], Lr[2], Lr[3], nbases);
	
	{
		int *p;
		S4_Simulation *S = S4_Simulation_New(Lr, nbases, NULL);
		plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
		p = mxGetData(plhs[0]);
		p[0] = S_add(S);
		SS[p[0]].is1d = is1d;
		
		MEX_TRACE("< S4mex_Simulation_New %d\n", p[0]);
	}
}
static void S4mex_Simulation_Destroy(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	const int id = S_check(prhs);
	MEX_TRACE("> S4mex_Simulation_Destroy(S = %p)\n", SS[id].S);
	S4_Simulation_Destroy(SS[id].S);
	S_del(id);
	MEX_TRACE("< S4mex_Simulation_Destroy %d\n", id);
}

static void S4mex_Simulation_SetMaterial(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_real eps[18] = {0};
	S4_MaterialID M;
	
	if(3 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected two additional arguments to 'set_material'.");
	}
	
	if(mxIsScalar(prhs[2])){
		eps[0] = *(mxGetPr(prhs[2]));
		if(mxIsComplex(prhs[2])){
			eps[1] = *(mxGetPi(prhs[2]));
			M = S4_Simulation_SetMaterial(S, -1, NULL, S4_MATERIAL_TYPE_SCALAR_COMPLEX, eps);
		}else{
			M = S4_Simulation_SetMaterial(S, -1, NULL, S4_MATERIAL_TYPE_SCALAR_REAL, eps);
		}
	}else if(3 == mxGetM(prhs[2]) && 3 == mxGetN(prhs[2])){
		int i, j;
		const double *p = mxGetPr(prhs[3]);
		for(j = 0; j < 3; ++j){
			for(i = 0; i < 3; ++i){
				eps[2*(i+j*3)+0] = p[i+j*3];
			}
		}
		if(mxIsComplex(prhs[2])){
			p = mxGetPi(prhs[2]);
			for(j = 0; j < 3; ++j){
				for(i = 0; i < 3; ++i){
					eps[2*(i+j*3)+1] = p[i+j*3];
				}
			}
		}
		{
			S4_real abcde[10] = {
				eps[0], eps[1], eps[6], eps[7],
				eps[2], eps[3], eps[8], eps[9],
				eps[16], eps[17]
			};
			for(i = 0; i < 10; ++i){
				eps[i] = abcde[i];
			}
		}
		M = S4_Simulation_SetMaterial(S, -1, NULL, S4_MATERIAL_TYPE_XYTENSOR_COMPLEX, eps);
	}else{
		mexErrMsgIdAndTxt("S4:invalidInput", "Third additional argument (epsilon) must be a scalar or 3x3 matrix.");
	}
	
	MEX_TRACE("> S4mex_Simulation_SetMaterial(S = %p, epsilon = %f ...)\n", S, eps[0]);
	
	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	*(int*)(mxGetData(plhs[0])) = M;
	
	MEX_TRACE("< S4mex_Simulation_SetMaterial %d\n", M);
}

static void S4mex_Simulation_AddLayer(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_real thickness;
	S4_MaterialID M;
	S4_LayerID L;
	
	if(4 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected three additional arguments to 'add_layer'.");
	}
	
	thickness = real_check(2, prhs, "thickness");
	
	M = ID_check(3, prhs);
	
	MEX_TRACE("> S4mex_Simulation_SetMaterial(S = %p, M = %d, thickness = %f)\n", S, M, thickness);
	
	L = S4_Simulation_SetLayer(S, -1, NULL, &thickness, -1, M);

	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	*(int*)(mxGetData(plhs[0])) = L;
	
	MEX_TRACE("< S4mex_Simulation_SetMaterial %d\n", L);
}

static void S4mex_Simulation_SetLayerThickness(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_LayerID L;
	S4_real thickness;
	
	if(5 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected four additional arguments to 'add_layer'.");
	}
	
	L = ID_check(3, prhs);
	thickness = real_check(4, prhs, "thickness");

	MEX_TRACE("> S4mex_Simulation_SetLayerThickness(S = %p, L = %d, thickness = %f)\n", S, L, thickness);
	
	S4_Simulation_SetLayer(S, L, NULL, &thickness, -1, -1);
	
	MEX_TRACE("< S4mex_Simulation_SetLayerThickness\n");
}

static void S4mex_Simulation_AddLayerCopy(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_real thickness;
	S4_LayerID L;
	S4_LayerID copy;
	
	if(4 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected three additional arguments to 'add_layer'.");
	}
	
	thickness = real_check(2, prhs, "thickness");
	
	copy = ID_check(3, prhs);

	MEX_TRACE("> S4mex_Simulation_AddLayerCopy(S = %p, thickness = %f, copy = %d)\n", S, thickness, copy);
	
	L = S4_Simulation_SetLayer(S, -1, NULL, &thickness, copy, -1);

	plhs[0] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	*(int*)(mxGetData(plhs[0])) = L;
	
	MEX_TRACE("< S4mex_Simulation_AddLayerCopy %d\n", L);
}

static void S4mex_Simulation_SetLayerRegion(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	const int id = S_check(prhs);
	S4_Simulation *S = SS[id].S;
	const int is1d = SS[id].is1d;
	const S4_LayerID L = ID_check(2, prhs);
	const S4_MaterialID M = ID_check(3, prhs);
	
	if(nrhs < 6){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected at least five additional arguments to 'set_layer_region'.");
	}
	
	MEX_TRACE("> S4mex_Simulation_SetLayerRegion(S = %p, L = %d, M = %d)\n", S, L, M);
	
	if(is1d){
		const double c = real_check(4, prhs, "center");
		const double h = real_check(5, prhs, "halfwidth");
		S4_real hw[2] = { h, 0 };
		S4_real center[2] = { c, 0 };
		S4_real angle_frac = 0;
		S4_Layer_SetRegionHalfwidths(S, L, M, S4_REGION_TYPE_INTERVAL, hw, center, &angle_frac);
	}else{
		const char *shape;
		if(!(mxIsChar(prhs[4]))){
			mexErrMsgIdAndTxt("S4:invalidInput", "First argument must be of type string.");
		}
		shape = mxArrayToString(prhs[4]);
		if(0 == strcmp("circle", shape)){
			double r = real_check(5, prhs, "radius");
			S4_real hw[2] = {r,r};
			S4_real center[2] = { 0, 0 };
			S4_real angle_frac = 0;
			if(nrhs > 6){
				vec2_check(6, prhs, "center", center);
			}
			MEX_TRACE("  circle(r = %f, center = %f, %f)\n", r, center[0], center[1]);
			S4_Layer_SetRegionHalfwidths(S, L, M, S4_REGION_TYPE_CIRCLE, hw, center, &angle_frac);
		}else if(0 == strcmp("ellipse", shape)){
			double hw[2];
			S4_real center[2] = { 0, 0 };
			S4_real angle_frac = 0;
			vec2_check(5, prhs, "halfwidths", hw);
			if(nrhs > 6){
				vec2_check(6, prhs, "center", center);
			}
			if(nrhs > 7){
				angle_frac = real_check(7, prhs, "rotation");
			}
			S4_Layer_SetRegionHalfwidths(S, L, M, S4_REGION_TYPE_ELLIPSE, hw, center, &angle_frac);
		}else if(0 == strcmp("rectangle", shape)){
			double hw[2];
			S4_real center[2] = { 0, 0 };
			S4_real angle_frac = 0;
			vec2_check(5, prhs, "halfwidths", hw);
			if(nrhs > 6){
				vec2_check(6, prhs, "center", center);
			}
			if(nrhs > 7){
				angle_frac = real_check(7, prhs, "rotation");
			}
			S4_Layer_SetRegionHalfwidths(S, L, M, S4_REGION_TYPE_RECTANGLE, hw, center, &angle_frac);
		}else if(0 == strcmp("polygon", shape)){
			int nv;
			const double *v;
			S4_real center[2] = { 0, 0 };
			S4_real angle_frac = 0;
			
			if(!mxIsNumeric(prhs[5]) || 2 != mxGetM(prhs[5]) || mxIsComplex(prhs[5])){
				mexErrMsgIdAndTxt("S4:invalidInput", "Expected polygon vertices.");
			}
			nv = (int)mxGetN(prhs[5]);
			v = mxGetPr(prhs[5]);
			
			if(nrhs > 6){
				vec2_check(6, prhs, "center", center);
			}
			if(nrhs > 7){
				angle_frac = real_check(7, prhs, "rotation");
			}
			S4_Layer_SetRegionVertices(S, L, M, S4_REGION_TYPE_POLYGON, nv, v, center, &angle_frac);
		}else{
			mexErrMsgIdAndTxt("S4:invalidInput", "Unrecognized shape.");
		}
	}
	
	MEX_TRACE("< S4mex_Simulation_SetLayerRegion\n");
}

static void S4mex_Simulation_SetPlanewave(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_real k[3], u[3], cu[2] = { 0, 0 }, cv[2] = { 0, 0 };
	if(6 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected five additional arguments to 'set_planewave'.");
	}
	
	vec3_check(2, prhs, "k", k);
	vec3_check(3, prhs, "u", u);
	complex_check(4, prhs, "cu", cu);
	complex_check(5, prhs, "cv", cv);
	
	MEX_TRACE("> S4mex_Simulation_SetPlanewave(S = %p, k = [%f, %f, %f], u = [%f, %f, %f], cu = %f + %f i, cv = %f + %f i)\n", S, k[0], k[1], k[2], u[0], u[1], u[2], cu[0], cu[1], cv[0], cv[1]);
	
	S4_Simulation_ExcitationPlanewave(S, k, u, cu, cv);
	
	MEX_TRACE("< S4mex_Simulation_SetPlanewave\n");
}

static void S4mex_Simulation_SetFrequency(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_real freq[2] = {0, 0};
	
	if(3 != nrhs){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected two additional arguments to 'set_frequency'.");
	}
	if(!mxIsScalar(prhs[2])){
		mexErrMsgIdAndTxt("S4:invalidInput", "First additional argument (frequency) must be a scalar number.");
	}
	freq[0] = *(mxGetPr(prhs[2]));
	
	MEX_TRACE("> S4mex_Simulation_SetFrequency(S = %p, freq = %f)\n", S, freq[0]);
	
	S4_Simulation_SetFrequency(S, freq);
	
	MEX_TRACE("< S4mex_Simulation_SetFrequency\n");
}

static void S4mex_Simulation_GetBases(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	const int id = S_check(prhs);
	S4_Simulation *S = SS[id].S;
	const int is1d = SS[id].is1d;
	
	double *p;
	int i, n;
	int *G = NULL;
	n = S4_Simulation_GetBases(S, G);
	G = (int*)malloc(sizeof(int)*2*n);
	S4_Simulation_GetBases(S, G);
	if(is1d){
		plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
		p = mxGetPr(plhs[0]);
		for(i = 0; i < n; ++i){
			p[i] = G[2*i+0];
		}
	}else{
		plhs[0] = mxCreateDoubleMatrix(2, n, mxREAL);
		p = mxGetPr(plhs[0]);
		for(i = 0; i < n; ++i){
			p[0+i*2] = G[2*i+0];
			p[1+i*2] = G[2*i+1];
		}
	}
	free(G);
}

static void S4mex_Simulation_GetPowerFlux(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_LayerID L;
	S4_real offset = 0;
	double *p;
	S4_real power[4];
	
	if(nrhs < 3){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected at least one additional argument to 'get_power_flux'.");
	}
	
	L = ID_check(2, prhs);
	
	if(nrhs > 3){
		offset = real_check(3, prhs, "offset");
	}
	
	MEX_TRACE("> S4mex_Simulation_GetPowerFlux(S = %p, L = %d, offset = %f)\n", S, L, offset);
	
	S4_Simulation_GetPowerFlux(S, L, &offset, power);

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
	p = mxGetPr(plhs[0]);
	*p = power[0];
	p = mxGetPi(plhs[0]);
	*p = power[2];
	
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
	p = mxGetPr(plhs[1]);
	*p = power[1];
	p = mxGetPi(plhs[1]);
	*p = power[3];
	
	MEX_TRACE("< S4mex_Simulation_GetPowerFlux %f + %f i, %f + %f i\n", power[0], power[2], power[1], power[3]);
}

static void S4mex_Simulation_GetWaves(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	S4_Simulation *S = Sptr_check(prhs);
	S4_LayerID L;
	S4_real *waves = NULL;
	int n, i, j;
	
	if(nrhs < 3){
		mexErrMsgIdAndTxt("S4:invalidInput", "Expected one additional argument to 'get_waves'.");
	}
	
	L = ID_check(2, prhs);
	
	n = S4_Simulation_GetBases(S, NULL);

	waves = (S4_real*)malloc(sizeof(S4_real)*11*2*n);
	S4_Simulation_GetWaves(S, L, waves);

	plhs[0] = mxCreateDoubleMatrix(8, n, mxCOMPLEX);
	plhs[1] = mxCreateDoubleMatrix(8, n, mxCOMPLEX);
	for(j = 0; j < 2; ++j){
		double *pr = mxGetPr(plhs[j]);
		double *pi = mxGetPi(plhs[j]);
		for(i = 0; i < n; ++i){
			const double *wave = &waves[11*(i+n*j)];
			pr[0+i*8] = wave[0];
			pr[1+i*8] = wave[1];
			pr[2+i*8] = wave[2];
			pi[2+i*8] = wave[3];
			pr[3+i*8] = wave[4];
			pr[4+i*8] = wave[5];
			pr[5+i*8] = wave[6];
			pr[6+i*8] = wave[7];
			pi[6+i*8] = wave[8];
			pr[7+i*8] = wave[9];
			pi[7+i*8] = wave[10];
		}
	}
	free(waves);
}

struct dispatch_entry{
	const char *cmd;
	void (*handler)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	static const struct dispatch_entry dispatch[] = {
		{ "init", &S4mex_init },
		{ "new", &S4mex_Simulation_New },
		{ "destroy", &S4mex_Simulation_Destroy },
		{ "set_material", &S4mex_Simulation_SetMaterial },
		{ "add_layer", &S4mex_Simulation_AddLayer },
		{ "set_layer_thickness", &S4mex_Simulation_SetLayerThickness },
		{ "add_layer_copy", &S4mex_Simulation_AddLayerCopy },
		{ "set_layer_region", &S4mex_Simulation_SetLayerRegion },
		{ "set_planewave", &S4mex_Simulation_SetPlanewave },
		{ "set_frequency", &S4mex_Simulation_SetFrequency },
		{ "get_bases", &S4mex_Simulation_GetBases },
		{ "get_power_flux", &S4mex_Simulation_GetPowerFlux },
		{ "get_waves", &S4mex_Simulation_GetWaves },
		{ NULL, NULL }
	};
	const struct dispatch_entry *entry = &dispatch[0];
	const char *cmd = NULL;
	
	mexAtExit(&ExitFcn);
	if(!(mxIsChar(prhs[0]))){
		mexErrMsgIdAndTxt("S4:invalidInput", "First argument must be of type string.");
	}
	cmd = mxArrayToString(prhs[0]);
	while(NULL != entry->cmd){
		if(0 == strcmp(cmd, entry->cmd)){
			entry->handler(nlhs, plhs, nrhs, prhs);
			return;
		}
		entry++;
	}
	if(NULL == entry->cmd){
		mexErrMsgIdAndTxt("S4:invalidInput", "Unrecognized command.");
	}
}
