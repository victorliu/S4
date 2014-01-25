#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "periodic_off2.h"
#include <Cgeom/geom_circum.h>
#include <Cgeom/geom_predicates.h>

static void skip(FILE *fp);

static int read_double(FILE *fp, double *val, const char *msg){
	skip(fp);
	if(1 != fscanf(fp, "%lf", val)){
		if(NULL != msg){
			fprintf(stderr, "Expected real number: %s\n", msg);
		}else{
			fprintf(stderr, "Expected real number\n");
		}
		return 1;
	}
	return 0;
}
static int read_int(FILE *fp, int *val, const char *msg){
	skip(fp);
	if(1 != fscanf(fp, "%d", val)){
		if(NULL != msg){
			fprintf(stderr, "Expected integer: %s\n", msg);
		}else{
			fprintf(stderr, "Expected integer\n");
		}
		return 1;
	}
	return 0;
}

int POFF2_LoadFromFile(POFF2 *off, const char *filename){
	FILE *fp;
	size_t i;
	int ret = 0;
	char line[13];
	if(NULL == off){ return 1; }
	if(NULL == filename){ return 2; }
	
	fp = fopen(filename, "rt");
	if(NULL == fp){ return -1; }
	
	fscanf(fp, "%12s", line); line[12] = '\0';
	if(0 != strcmp("PeriodicOFF2", line)){ return -2; }
	
	if(!read_double(fp, &(off->L[0]), "lattice ux")){ return -10; }
	if(!read_double(fp, &(off->L[1]), "lattice uy")){ return -10; }
	if(!read_double(fp, &(off->L[2]), "lattice vx")){ return -10; }
	if(!read_double(fp, &(off->L[3]), "lattice vy")){ return -10; }
	if(!read_int(fp, &(off->n_verts), "num vertices")){ return -10; }
	if(!read_int(fp, &(off->n_faces), "num faces")){ return -10; }
	if(!read_int(fp, &(off->n_matchings), "num matchings")){ return -10; }
	POFF2_Initialize(off, off->n_verts, off->n_faces, off->n_matchings);
	
	for(i = 0; i < off->n_verts; ++i){
		if(!read_double(fp, &(off->vert[i].x), "vertex x")){ return -11; }
		if(!read_double(fp, &(off->vert[i].y), "vertex x")){ return -11; }
		if(!read_int(fp, &(off->vert[i].flags), "vertex flags")){ return -11; }
	}
	
	for(i = 0; i < off->n_faces; ++i){
		size_t j;
		if(!read_int(fp, &(off->face[i].n_verts), "face num vertices")){ return -12; }
		POFF2_InitializeFace(&(off->face[i]), off->face[i].n_verts);
		for(j = 0; j < off->face[i].n_verts; ++j){
			if(!read_int(fp, &(off->face[i].vert[j]), "face vertex index")){ return -12; }
		}
		if(!read_int(fp, &(off->face[i].flags), "face flags")){ return -12; }
	}
	
	for(i = 0; i < off->n_matchings; ++i){
		if(!read_int(fp, &(off->matching[i].face[0]), "matching face0 index")){ return -13; }
		if(!read_int(fp, &(off->matching[i].vert1[0]), "matching face0 vert0 offset")){ return -13; }
		if(!read_int(fp, &(off->matching[i].vert2[0]), "matching face0 vert1 offset")){ return -13; }
		if(!read_int(fp, &(off->matching[i].face[1]), "matching face1 index")){ return -13; }
		if(!read_int(fp, &(off->matching[i].vert1[1]), "matching face1 vert0 offset")){ return -13; }
		if(!read_int(fp, &(off->matching[i].vert2[1]), "matching face1 vert1 offset")){ return -13; }
		if(!read_int(fp, &(off->matching[i].flags), "matching flags")){ return -13; }
	}
	
	if(ferror(fp)){
		ret = -3;
	}
	fclose(fp);
	return ret;
}
int POFF2_SaveToFile(const POFF2 *off, const char *filename){
	FILE *fp;
	int i;
	int ret = 0;
	if(NULL == off){ return 1; }
	if(NULL == filename){ return 2; }
	
	fp = fopen(filename, "wt");
	if(NULL == fp){ return -1; }
	
	fprintf(fp, "PeriodicOFF2\n\n# Lattice\n%g %g %g %g\n\n# number of verts, faces, matchings\n%d %d %d\n# List of vertices: x, y, flags\n",
		off->L[0], off->L[1], off->L[2], off->L[3],
		off->n_verts, off->n_faces, off->n_matchings);

	for(i = 0; i < off->n_verts; ++i){
		fprintf(fp, "%.16f\t%.16f\t%d\n", off->vert[i].x, off->vert[i].y, off->vert[i].flags);
	}
	
	fprintf(fp, "\n\n# List of faces: num-vertices, vertex1, vertex2, ..., vertexN, flags\n");
	
	for(i = 0; i < off->n_faces; ++i){
		int j;
		fprintf(fp, "%d", off->face[i].n_verts);
		for(j = 0; j < off->face[i].n_verts; ++j){
			fprintf(fp, "\t%d", off->face[i].vert[j]);
		}
		fprintf(fp, "\t%d\n", off->face[i].flags);
	}
	
	fprintf(fp, "\n\n# List of matchings: face1, vtx11, vtx12, face2, vtx21, vtx22, flags\n");
	
	for(i = 0; i < off->n_matchings; ++i){
		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			off->matching[i].face[0],
			off->matching[i].vert1[0],
			off->matching[i].vert2[0],
			off->matching[i].face[1],
			off->matching[i].vert1[1],
			off->matching[i].vert2[1],
			off->matching[i].flags);
	}
	
	if(ferror(fp)){
		ret = -3;
	}
	fclose(fp);
	return ret;
}
void POFF2_InitializeFace(POFF2_Face *face, int nverts){
	if(NULL == face){ return; }
	face->n_verts = nverts;
	face->vert = (int*)malloc(sizeof(int) * nverts);
}
void POFF2_Initialize(POFF2 *off, int nverts, int nfaces, int nmatchings){
	if(NULL == off){ return; }
	off->n_verts = nverts;
	off->n_faces = nfaces;
	off->n_matchings = nmatchings;
	off->vert = (POFF2_Vertex*)malloc(sizeof(POFF2_Vertex) * nverts);
	off->face = (POFF2_Face*)malloc(sizeof(POFF2_Face) * nfaces);
	off->matching = (POFF2_Matching*)malloc(sizeof(POFF2_Matching) * nmatchings);
}
void POFF2_Destroy(POFF2 *off){
	int i;
	if(NULL == off){ return; }
	free(off->vert);
	for(i = 0; i < off->n_faces; ++i){
		free(off->face[i].vert);
	}
	free(off->face);
	free(off->matching);
}
int POFF2_Clone(const POFF2 *off, POFF2 *cpy){
	int i;
	if(NULL == off){ return -1; }
	if(NULL == cpy){ return -2; }
	cpy->L[0] = off->L[0];
	cpy->L[1] = off->L[1];
	cpy->L[2] = off->L[2];
	cpy->L[3] = off->L[3];
	cpy->n_verts = off->n_verts;
	cpy->n_faces = off->n_faces;
	cpy->n_matchings = off->n_matchings;
	cpy->vert = (POFF2_Vertex*)malloc(sizeof(POFF2_Vertex) * cpy->n_verts);
	memcpy(cpy->vert, off->vert, sizeof(POFF2_Vertex) * cpy->n_verts);
	cpy->face = (POFF2_Face*)malloc(sizeof(POFF2_Face) * cpy->n_faces);
	memcpy(cpy->face, off->face, sizeof(POFF2_Face) * cpy->n_faces);
	cpy->matching = (POFF2_Matching*)malloc(sizeof(POFF2_Matching) * cpy->n_matchings);
	memcpy(cpy->matching, off->matching, sizeof(POFF2_Matching) * cpy->n_matchings);
	for(i = 0; i < cpy->n_faces; ++i){
		cpy->face[i].vert = (int*)malloc(sizeof(int) * cpy->face[i].n_verts);
		memcpy(cpy->face[i].vert, off->face[i].vert, sizeof(int) * cpy->face[i].n_verts);
	}
	return 0;
}

// Skip whitespace
// Skip shell script style comments starting with #
// Skip C style single line comments, and multi-line comments
static void skip(FILE *fp){
	//while(isspace(is.peek())){ is.ignore(1); }
	//is >> std::ws;
	// parser inspired by JSON_parser.c
	enum states{
		GO,//S_DEFAULT,           // default state, not in comment
		SL,//S_GOT_SLASH,         // got the first slash, possibly starting a comment
		LN,//S_IN_LINE_COMMENT,   // in a single line comment
		BK,//S_IN_BLOCK_COMMENT,  // in a block (multi-line) comment
		ST,//S_GOT_STAR           // got (possibly) ending star of block comment,
		// these two are not real states; reaching them exits the function
		//S_ERROR_SLASH,       // non space, need to put back a slash
		//S_ERROR,             // non space
		NR_STATES
	};
	enum actions{ // these must be negative (to be different from states)
		RET = -1, // return
		RPS = -2  // return and put back a slash character
	};
	enum classes{
		C_SLASH,
		C_HASH,
		C_STAR,
		C_EOL,
		C_SPACE,
		C_ETC,
		NR_CLASSES
	};
		
	static int ascii_class[64] = {
	// This array maps the lower 64 ASCII characters into character classes.
	// The remaining characters should be mapped to C_ETC.
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_SPACE, C_EOL,  C_ETC,   C_ETC,   C_EOL,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,

		C_SPACE, C_ETC,   C_ETC,  C_HASH,  C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_STAR, C_ETC,   C_ETC,   C_ETC,   C_ETC, C_SLASH,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC,
		C_ETC,   C_ETC,   C_ETC,  C_ETC,   C_ETC,   C_ETC,   C_ETC,   C_ETC
	};

	
	/*
	State transition table:
	State\Char |  /   |  #   |  *   |  \n  | space (non \n) | anything else
	-----------------------------------------------------------------------
	       DEF | GOT/ | LINE | ERR  | DEF  | DEF            | ERR
	      GOT/ | LINE | ERR/ | BLOK | ERR/ | ERR/           | ERR/
	      LINE | LINE | LINE | LINE | DEF  | LINE           | LINE
	      BLOK | BLOK | BLOK | GOT* | BLOK | BLOK           | BLOK
	      GOT* | DEF  | BLOK | GOT* | BLOK | BLOK           | BLOK
	      ERR/ | ERR  | ERR  | ERR  | ERR  | ERR            | ERR
	       ERR | ERR  | ERR  | ERR  | ERR  | ERR            | ERR
	When hitting the error state, we may have to put back some chars into the stream
	*/
	
		
	static int state_transition_table[NR_STATES][NR_CLASSES] = {
	// The state transition table takes the current state and the current symbol,
	// and returns either a new state or an action.

	//       /    #    *   EOL space etc
	/*GO*/ {GO , LN , RET, GO , GO , RET},
	/*SL*/ {LN , RPS, BK , RPS, RPS, RPS},
	/*LN*/ {LN , LN , LN , GO , LN , LN },
	/*BK*/ {BK , BK , ST , BK , BK , BK },
	/*ST*/ {GO , BK , ST , BK , BK , BK }
	};

	int state = GO;
	while(!feof(fp)){
		int cclass;
		int c = fgetc(fp); ungetc(c, fp);
		cclass = (c < 64) ? ascii_class[c] : C_ETC;
		state = state_transition_table[state][cclass];
		if(state < 0){
			if(RET == state){
				return;
			}else if(RPS == state){
				ungetc('/', fp);
				return;
			}
		}else{
			fgetc(fp);
			if(C_EOL == cclass){
				if('\r' == c){ // possibly need to eat a \n
					if((c = fgetc(fp)) == '\n'){
						// ignore
					}else{
						ungetc(c, fp);
					}
				}
			}
		}
	}
}










typedef struct{
	int from; /* vertorg index */
	int next; /* halfedge index of next in loop */
	int face; /* face index */
	int flag;
	/* prop:
	 *   low bit specifies if this halfedge is translated by first lattice vector
	 *   2nd bit specifies if this halfedge is translated by second lattice vector
	 */
} POFF2Mesh_Halfedge;

struct POFF2Mesh_{
	/* n_verts: number of unique vertices
	 * n_vertorgs: number of original vertices in the POFF2
	 * n_halfs: number of halfedges
	 * n_faces: number of faces
	 * n_halfints = number of interior halfedges (these are at the beginning of the `half' array
	 */
	int n_verts, n_vertorgs, n_halfs, n_faces, n_halfints;
	POFF2_Vertex *vertorg;
	POFF2Mesh_Halfedge *half; /* stored in pairs of opposite halfedges */
	int *vert; /* index of an outgoing halfedge */
	int *vmap; /* length 2*n_vertorgs, map from vertorg to vert, second in pair is flags for offsets */
	/* lowest two bits is the first direction offset (0 = 0, 1 = +1, 3 = -1) */
	int *face; /* index of one halfedge in face */
	int *faceflag;
	double *center; /* list of face centers */
	int max_valence, max_dual_valence;
	double L[4];
};
/* offmap is used to map from the bits in 2nd col of vmap to actual lattice offsets */
/* Note that 0x2 is an invalid bit pattern */
static const int offmap[4] = { 0, 1, 255, -1 };
/* offmap2 is used to map from actual offset (plus 2) to 2nd-col-vmap bits */
static const int offmap2[4] = { 0x2, 0x3, 0x0, 0x1 };

static void POFF2Mesh_SwapHalfedges(const POFF2Mesh M, const int k[2]){
	int t;
	int hp[2];
	
	if(k[0] == k[1]){ return; }
	
	/* First update references in face */
	M->face[M->half[k[0]].face] = M->half[k[0]].next;
	M->face[M->half[k[1]].face] = M->half[k[1]].next;
	/* Update references in vert. Some of these could be -1's so we can't just
	 * always set it to the swapped halfedge index.
	 */
	if(M->vert[M->half[k[0]].from] == k[0]){ M->vert[M->half[k[0]].from] = k[1]; }
	if(M->vert[M->half[k[1]].from] == k[1]){ M->vert[M->half[k[1]].from] = k[0]; }
	/* update the previous halfedges' next pointers */
	for(t = 0; t < 2; ++t){
		int h = M->half[k[t]].next;
		while(M->half[h].next != k[t]){
			h = M->half[h].next;
		}
		hp[t] = h;
	}
	/* We must apply the changes atomically here in case the two halfedges
	 * are consecutive (one is the other's next)
	 */
	M->half[hp[0]].next = k[1];
	M->half[hp[1]].next = k[0];
	/* Now perform the swap */
	t = M->half[k[0]].from; M->half[k[0]].from = M->half[k[1]].from; M->half[k[1]].from = t;
	t = M->half[k[0]].next; M->half[k[0]].next = M->half[k[1]].next; M->half[k[1]].next = t;
	t = M->half[k[0]].face; M->half[k[0]].face = M->half[k[1]].face; M->half[k[1]].face = t;
	t = M->half[k[0]].flag; M->half[k[0]].flag = M->half[k[1]].flag; M->half[k[1]].flag = t;
}

/* iv[0],iv[1] are from,to vertices of halfedge h[0] */
/* iv[2],iv[3] are from,to vertices of halfedge h[1] */
static int POFF2Mesh_SetPeriodicity(const POFF2Mesh M, const double Linv[4], double tol, const int h[2], const int iv[4]){
	double v[2], u[2];
	int nonzero = 0;
	
	v[0] = M->vertorg[iv[3]].x - M->vertorg[iv[0]].x;
	v[1] = M->vertorg[iv[3]].y - M->vertorg[iv[0]].y;
	u[0] = Linv[0] * v[0] + Linv[2] * v[1];
	u[1] = Linv[1] * v[0] + Linv[3] * v[1];

	if(fabs(u[0] - 1.) < tol){
		M->half[h[1]].flag |= 0x1;
		nonzero = 1;
	}else if(fabs(u[0] + 1.) < tol){
		M->half[h[0]].flag |= 0x1;
		nonzero = 1;
	}
	if(fabs(u[1] - 1.) < tol){
		M->half[h[1]].flag |= 0x2;
		nonzero = 1;
	}else if(fabs(u[1] + 1.) < tol){
		M->half[h[0]].flag |= 0x2;
		nonzero = 1;
	}
	return (nonzero == 0);
}

POFF2Mesh POFF2Mesh_Create(const POFF2 *off, double tol){
	if(NULL == off){ return NULL; }
	POFF2Mesh M = (POFF2Mesh)malloc(sizeof(struct POFF2Mesh_));
	memset(M, 0, sizeof(struct POFF2Mesh_));
	
	if(tol <= 0){ tol = DBL_EPSILON; }
	
	M->n_vertorgs = off->n_verts;
	M->vertorg = (POFF2_Vertex*)malloc(sizeof(POFF2_Vertex) * M->n_vertorgs);
	memcpy(M->vertorg, off->vert, sizeof(POFF2_Vertex) * M->n_vertorgs);
	
	M->n_faces = off->n_faces;
	M->n_verts = 0;
	M->face = (int*)malloc(sizeof(int) * off->n_faces);
	M->faceflag = (int*)malloc(sizeof(int) * off->n_faces);
	/* We over-allocate M->vert on purpose since we need the extra space */
	M->vert = (int*)malloc(sizeof(int) * off->n_verts);
	M->L[0] = off->L[0];
	M->L[1] = off->L[1];
	M->L[2] = off->L[2];
	M->L[3] = off->L[3];

	/* First count number of halfedges and allocate */
	M->n_halfs = 0;
	M->max_dual_valence = 0;
	{
		int i;
		for(i = 0; i < off->n_faces; ++i){
			M->n_halfs += off->face[i].n_verts;
			if(off->face[i].n_verts > M->max_dual_valence){
				M->max_dual_valence = off->face[i].n_verts;
			}
		}
	}
	if(M->n_halfs % 2 == 1){ goto error; } /* We had better have an even number */
	M->half = (POFF2Mesh_Halfedge*)malloc(sizeof(POFF2Mesh_Halfedge) * M->n_halfs);
	
	/* Create all the halfedge loops for each face */
	{
		int i;
		int h = 0;
		for(i = 0; i < off->n_faces; ++i){
			int j, k;
			const int nv = off->face[i].n_verts;
			M->face[i] = h;
			M->faceflag[i] = off->face[i].flags;
			for(j = nv-1, k = 0; k < nv; j = k++){
				const int iv = off->face[i].vert[j];
				M->half[h+j].from = iv;
				M->half[h+j].next = h+k;
				M->half[h+j].face = i;
				M->half[h+j].flag = 0;
				M->vert[iv] = h+j;
			}
			h += nv;
		}
	}
	
	/* print out all halfedges */
	/*{
		int i;
		for(i = 0; i < M->n_halfs; ++i){
			printf("half\t%d: from = %d, next = %d, face = %d, flag = %d\n",
				i,
				M->half[i].from,
				M->half[i].next,
				M->half[i].face,
				M->half[i].flag
			);
		}
	}*/
	/* Now attempt to match each halfedge with its opposite */
	const int nhm = M->n_halfs/2 - off->n_matchings;
	int iend = M->n_halfs;
	{
		/* nhm is the number of halfedge pairs we should be finding */
		int pairs = 0;
		int i = 0;
		int safecount = 0; /* ensure we don't end up in infinite loop */
		while(i < iend && safecount < M->n_halfs){
			const int iv0 = M->half[i].from;
			const int iv1 = M->half[M->half[i].next].from;
			/* Look for the opposite half[i] */
			int found = 0;
			int j;
			for(j = i+1; j < iend; ++j){
				const int jv0 = M->half[j].from;
				const int jv1 = M->half[M->half[j].next].from;
				if(jv0 == iv1 && jv1 == iv0){
					found = 1;
					pairs++;
					break;
				}
			}
			
			/* Possibly swap */
			{
				int doswap = 0;
				int k[2] = { i, i };
				if(found){
					if(i+1 != j){
						k[0] = i+1; k[1] = j;
						doswap = 1;
					}
					i += 2;
				}else{
					k[0] = i; k[1] = iend-1;
					doswap = 1;
					iend--;
				}
				/* Note: cannot reference variable `i' after this point */
				if(doswap){
					POFF2Mesh_SwapHalfedges(M, k);
				}
			}
			safecount++;
		}
		if(nhm != pairs){ goto error; }
	}

	/* Check some invariants */
	{
		int i;
		for(i = 0; i < M->n_faces; ++i){
			if(!(i == M->half[M->face[i]].face)){ goto error; }
		}
		for(i = 0; i < M->n_halfs; ++i){
			if(!(M->half[i].from < M->n_vertorgs)){ goto error; }
			if(!(M->half[i].face < M->n_faces)){ goto error; }
			if(!(M->half[i].next < M->n_halfs)){ goto error; }
		}
	}
	
	/* Euler characteristic says n_verts - n_halfs/2 + n_faces = 0 */
	M->n_verts = M->n_halfs/2 - M->n_faces;
	if(M->n_verts > M->n_vertorgs){ goto error; }
	
	M->vmap = (int*)malloc(sizeof(int) * 2*off->n_verts);
	{
		int i;
		for(i = 0; i < M->n_vertorgs; ++i){
			M->vmap[2*i] = i;
		}
	}

	/* Go through all the boundary edge matchings */
	{
		double Linv[4];
		int im;
		int nextind = 2*nhm;
		
		/* Compute L inverse */
		{
			double idet = 1. / (M->L[0]*M->L[3] - M->L[1]*M->L[2]);
			Linv[0] = idet *  M->L[3];
			Linv[1] = idet * -M->L[1];
			Linv[2] = idet * -M->L[2];
			Linv[3] = idet *  M->L[0];
		}
		for(im = 0; im < off->n_matchings; ++im){
			const POFF2_Matching *m = &off->matching[im];
			int hi[2] = { -1, -1 };
			int ifi;
			for(ifi = 0; ifi < 2; ++ifi){
				const int iface = m->face[ifi];
				const int nv = off->face[iface].n_verts;
				int h = M->face[iface];
				int safecount = nv;
				while(safecount > 0 &&
					!(
						(M->half[h].from == m->vert1[ifi] && M->half[M->half[h].next].from == m->vert2[ifi]) ||
						(M->half[h].from == m->vert2[ifi] && M->half[M->half[h].next].from == m->vert1[ifi])
					)
				){
					h = M->half[h].next;
				}
				if(safecount > 0){
					hi[ifi] = h;
				}
			}
			
			/*
			printf("matching %d:\n", im);
			printf("  halfedge %d from,to = %d, %d\n", hi[0], M->half[hi[0]].from, M->half[M->half[hi[0]].next].from);
			printf("  halfedge %d from,to = %d, %d\n", hi[1], M->half[hi[1]].from, M->half[M->half[hi[1]].next].from);
			*/
			if(-1 == hi[0] || -1 == hi[1] || hi[0] < nextind || hi[1] < nextind){
				goto error;
			}
			
			/* Update the vertex map */
			{
				const int iv[4] = {
					M->half[hi[0]].from, M->half[M->half[hi[0]].next].from,
					M->half[hi[1]].from, M->half[M->half[hi[1]].next].from
				};
				if(iv[0] > iv[3]){
					M->vert[iv[0]] = -1;
					if(iv[3] < M->vmap[2*iv[0]]){ M->vmap[2*iv[0]] = iv[3]; }
					if(iv[3] < M->vmap[2*iv[3]]){ M->vmap[2*iv[3]] = iv[3]; }
				}else{
					M->vert[iv[3]] = -1;
					if(iv[0] < M->vmap[2*iv[0]]){ M->vmap[2*iv[0]] = iv[0]; }
					if(iv[0] < M->vmap[2*iv[3]]){ M->vmap[2*iv[3]] = iv[0]; }
				}
				if(iv[1] > iv[2]){
					M->vert[iv[1]] = -1;
					if(iv[2] < M->vmap[2*iv[1]]){ M->vmap[2*iv[1]] = iv[2]; }
					if(iv[2] < M->vmap[2*iv[2]]){ M->vmap[2*iv[2]] = iv[2]; }
				}else{
					M->vert[iv[2]] = -1;
					if(iv[1] < M->vmap[2*iv[1]]){ M->vmap[2*iv[1]] = iv[1]; }
					if(iv[1] < M->vmap[2*iv[2]]){ M->vmap[2*iv[2]] = iv[1]; }
				}
				/*
				printf(" 1 vmap[%d] = %d\n", iv[0], M->vmap[2*iv[0]]);
				printf(" 2 vmap[%d] = %d\n", iv[1], M->vmap[2*iv[1]]);
				printf(" 3 vmap[%d] = %d\n", iv[2], M->vmap[2*iv[2]]);
				printf(" 4 vmap[%d] = %d\n", iv[3], M->vmap[2*iv[3]]);
				*/
				/* Check that the periodicity makes sense */
				if(0 != POFF2Mesh_SetPeriodicity(M, Linv, tol, hi, iv)){ goto error; }
			}
	
			/* Swap things into place */
			{
				if(hi[0] > hi[1]){ int t = hi[0]; hi[0] = hi[1]; hi[1] = t; }
				int k[2];
				k[0] = nextind;
				k[1] = hi[0];
				POFF2Mesh_SwapHalfedges(M, k);
				k[0] = nextind+1;
				k[1] = hi[1];
				POFF2Mesh_SwapHalfedges(M, k);
				nextind += 2;
			}
		}
	}
	
	/* Compact the list of vertex halfedge pointers.
	 * While we do this, also store the reverse mapping from
	 * vertorg index to compacted vert index in the second column of vmap
	 */
	{
		int i;
		int lookahead = 0;
		for(i = 0; i < M->n_verts; ++i){
			while(-1 == M->vert[i+lookahead]){ ++lookahead; }
			M->vert[i] = M->vert[i+lookahead];
			M->vmap[2*(i+lookahead)+1] = i;
		}
	}
	
	/* Chase all vmap mappings down to their end state */
	{
		int i;
		for(i = 0; i < M->n_vertorgs; ++i){
			while(M->vmap[2*M->vmap[2*i]] != M->vmap[2*i]){
				M->vmap[2*i] = M->vmap[2*M->vmap[2*i]];
			}
		}
	}
	
	/* Convert the map into vert indices instead of vertorg indices.
	 * First column of vmap maps vertorg indices to the unique vertorg index,
	 * and the second column can then map that unique vertorg index to the
	 * new compacted vert index.
	 */
	{
		int i;
		for(i = 0; i < M->n_vertorgs; ++i){
			const int uniq_vertorg_index = M->vmap[2*i];
			const int compacted_vert_index = M->vmap[2*uniq_vertorg_index+1];
			M->vmap[2*i] = compacted_vert_index;
		}
	}
	
	
	/* Determine maximum valence by spinning around each vertex.
	 * We will also determine all the relative offsets of the non-vert
	 * vertorg indices.
	 */
	M->max_valence = 0;
	{
		int i;
		for(i = 0; i < M->n_verts; ++i){
			int h = M->vert[i];
			const int vo0 = M->half[h].from;
			int vo;
			int valence = 0;
			int safety = M->n_verts; /* to guarantee termination */
			int curoff[2] = {0,0};

			M->vmap[2*vo0+1] = 0x0;
			do{
				const int hsave = h;
				const int hopp = 1 ^ h; // hopp = h.opp
				h = M->half[hopp].next; /* h <- h.opp.next */

				/* Compute offsets */
				vo = M->half[h].from;
				if(vo0 != vo){
					M->vmap[2*vo+1] = 0x0;
					/* We spun around some periodic boundary */
					const int off0 = (int)(M->half[hsave].flag & 0x1) - (int)(M->half[hopp].flag & 0x1);
					const int off1 = (int)((M->half[hsave].flag & 0x2) >> 1) - (int)((M->half[hopp].flag & 0x2) >> 1);

					curoff[0] += off0;
					curoff[1] += off1;
					if(curoff[0] < -1 || curoff[0] > 1 || curoff[1] < -1 || curoff[1] > 1){
						goto error;
					}
					M->vmap[2*vo+1] = ((offmap2[curoff[1]+2] << 2) | offmap2[curoff[0]+2]);
				}
				
				safety--;
				++valence;
			}while(h != M->vert[i] && safety > 0);
			if(valence > M->max_valence){ M->max_valence = valence; }
		}
	}
	
	/* Compute face centers */
	M->center = (double*)malloc(sizeof(double) * 2 * M->n_faces);
	{
		int i;
		double *p = (double*)malloc(sizeof(double) * 2 * M->max_dual_valence);
		for(i = 0; i < M->n_faces; ++i){
			if(off->face[i].n_verts < 3){
				goto error;
			}else if(3 == off->face[i].n_verts){
				p[0] = M->vertorg[off->face[i].vert[0]].x;
				p[1] = M->vertorg[off->face[i].vert[0]].y;
				p[2] = M->vertorg[off->face[i].vert[1]].x;
				p[3] = M->vertorg[off->face[i].vert[1]].y;
				p[4] = M->vertorg[off->face[i].vert[2]].x;
				p[5] = M->vertorg[off->face[i].vert[2]].y;
				geom_circum_tri2d(&p[2*0], &p[2*1], &p[2*2], &M->center[2*i], NULL, NULL);
			}else{
				int j;
				double ret;
				double r;
				for(j = 0; j < off->face[i].n_verts; ++j){
					p[2*j+0] = M->vertorg[off->face[i].vert[j]].x;
					p[2*j+1] = M->vertorg[off->face[i].vert[j]].y;
				}
				ret = geom_circum_fit2d(off->face[i].n_verts, p, &M->center[2*i], &r);
				if(ret / r > tol){
					goto error;
				}
			}
		}
		free(p);
	}
	return M;
error:
	POFF2Mesh_Destroy(M);
	return NULL;
}
void POFF2Mesh_Destroy(const POFF2Mesh M){
	if(NULL == M){ return; }
	free(M->vertorg);
	free(M->half);
	free(M->vert);
	free(M->face);
	free(M->faceflag);
	free(M->center);
	free(M);
}
int POFF2Mesh_SavePostscript(const POFF2Mesh M, const char *filename){
	FILE *fp;
	
	if(NULL == M){ return -1; }
	if(NULL == filename){ return -2; }
	
	fp = fopen(filename, "wb");
	if(NULL == fp){ return 1; }
	
	/* Output header */
	{
		fprintf(fp, 
			"72 72 scale\n"
			"3 3 scale\n"
			"1 2 translate\n"
			"0.001 setlinewidth\n"
			"/Courier findfont 0.03 scalefont setfont\n"
		);
	}
	
	/* Draw vertices */
	{
		int i;
		fprintf(fp, "0.5 setgray\n");
		for(i = 0; i < M->n_vertorgs; ++i){
			fprintf(fp, "newpath %g %g 0.002 0 360 arc fill\n",
				M->vertorg[i].x, M->vertorg[i].y
			);
			fprintf(fp, "newpath %g %g moveto (vo%d) show\n",
				M->vertorg[i].x, M->vertorg[i].y, i
			);
		}
	}
	/* Draw unique vertices */
	{
		int i;
		fprintf(fp, "0 setgray\n");
		for(i = 0; i < M->n_verts; ++i){
			const int vi = M->half[M->vert[i]].from;
			fprintf(fp, "newpath %g %g 0.003 0 360 arc stroke\n",
				M->vertorg[vi].x, M->vertorg[vi].y
			);
			fprintf(fp, "newpath %g %g 0.03 add moveto (v%d) show\n",
				M->vertorg[vi].x, M->vertorg[vi].y, i
			);
		}
	}
	/* Draw face halfedges */
	{
		int i;
		fprintf(fp, "%% face halfedges\n");
		fprintf(fp, "0 0 1 setrgbcolor\n");
		for(i = 0; i < M->n_faces; ++i){
			int h = M->face[i];
			int hn = M->half[h].next;
			do{
				double v0[2] = {
					M->vertorg[M->half[h].from].x,
					M->vertorg[M->half[h].from].y
				};
				double v1[2] = {
					M->vertorg[M->half[hn].from].x,
					M->vertorg[M->half[hn].from].y
				};
				double v2[2] = {
					0.2*v0[0]+0.8*v1[0],
					0.2*v0[1]+0.8*v1[1]
				};
				v0[0] = 0.9*v0[0] + 0.1*M->center[2*i+0];
				v0[1] = 0.9*v0[1] + 0.1*M->center[2*i+1];
				v1[0] = 0.9*v1[0] + 0.1*M->center[2*i+0];
				v1[1] = 0.9*v1[1] + 0.1*M->center[2*i+1];
				v2[0] = 0.8*v2[0] + 0.2*M->center[2*i+0];
				v2[1] = 0.8*v2[1] + 0.2*M->center[2*i+1];
				fprintf(fp, "newpath %g %g moveto %g %g lineto %g %g lineto stroke\n",
					v0[0], v0[1], v1[0], v1[1], v2[0], v2[1]
				);
				h = hn;
				hn = M->half[hn].next;
			}while(h != M->face[i]);
		}
	}
	/* Draw face centers */
	{
		int i;
		fprintf(fp, "%% face centers\n");
		fprintf(fp, "1 0 0 setrgbcolor\n");
		for(i = 0; i < M->n_faces; ++i){
			fprintf(fp, "newpath %g %g 0.004 0 360 arc fill\n",
				M->center[2*i+0], M->center[2*i+1]
			);
		}
	}
	/* Draw dual edge diamonds */
	{
		int i;
		const int ne = POFF2Mesh_NumEdges(M);
		const double *L = &(M->L[0]);
		fprintf(fp, "%% dual edge diamonds\n");
		fprintf(fp, "0 1 0 setrgbcolor\n");
		for(i = 0; i < ne; ++i){
			int j;
			POFF2Mesh_Index iface[2];
			POFF2Mesh_Index ivert[2];
			double p[8], c[2];
			POFF2Mesh_GetEdgeFaces(M, i, iface);
			POFF2Mesh_GetEdgeVertices(M, i, ivert);
	
			POFF2Mesh_GetVertex(M, ivert[0].idx, &p[0]);
			p[0] += (L[0] * ivert[0].off[0] + L[2] * ivert[0].off[1]);
			p[1] += (L[1] * ivert[0].off[0] + L[3] * ivert[0].off[1]);
			POFF2Mesh_GetFaceCenter(M, iface[1].idx, &p[2]);
			p[2] += (L[0] * iface[1].off[0] + L[2] * iface[1].off[1]);
			p[3] += (L[1] * iface[1].off[0] + L[3] * iface[1].off[1]);
			POFF2Mesh_GetVertex(M, ivert[1].idx, &p[4]);
			p[4] += (L[0] * ivert[1].off[0] + L[2] * ivert[1].off[1]);
			p[5] += (L[1] * ivert[1].off[0] + L[3] * ivert[1].off[1]);
			POFF2Mesh_GetFaceCenter(M, iface[0].idx, &p[6]);
			p[6] += (L[0] * iface[0].off[0] + L[2] * iface[0].off[1]);
			p[7] += (L[1] * iface[0].off[0] + L[3] * iface[0].off[1]);
			
			c[0] = 0.5*p[0] + 0.5*p[4];
			c[1] = 0.5*p[1] + 0.5*p[5];
			
			for(j = 0; j < 4; ++j){
				p[2*j+0] = 0.9 * p[2*j+0] + 0.1 * c[0];
				p[2*j+1] = 0.9 * p[2*j+1] + 0.1 * c[1];
			}
			fprintf(fp, "newpath");
			for(j = 0; j < 4; ++j){
				fprintf(fp, " %g %g moveto", p[0], p[1]);
				fprintf(fp, " %g %g lineto", p[2], p[3]);
				fprintf(fp, " %g %g lineto", p[4], p[5]);
				fprintf(fp, " %g %g lineto closepath stroke\n", p[6], p[7]);
			}
		}
	}
	/* Draw dual vertex polygons */
	{
		int i;
		const int nv = POFF2Mesh_NumVertices(M);
		const double *L = &(M->L[0]);
		double *p = (double*)malloc(sizeof(double) * 2*POFF2Mesh_MaxValence(M));
		POFF2Mesh_Index *iface = (POFF2Mesh_Index*)malloc(sizeof(POFF2Mesh_Index) * POFF2Mesh_MaxValence(M));
		fprintf(fp, "%% dual vertex polygons\n");
		fprintf(fp, "1 0 1 setrgbcolor\n");
		for(i = 0; i < nv; ++i){
			int j;
			double c[2];
			const int nf = POFF2Mesh_GetVertexFaces(M, i, iface);
			POFF2Mesh_GetVertex(M, i, c);
			for(j = 0; j < nf; ++j){
				POFF2Mesh_GetFaceCenter(M, iface[j].idx, &p[2*j]);
				p[2*j+0] += (L[0] * iface[j].off[0] + L[2] * iface[j].off[1]);
				p[2*j+1] += (L[1] * iface[j].off[0] + L[3] * iface[j].off[1]);
				p[2*j+0] = 0.9*p[2*j+0] + 0.1*c[0];
				p[2*j+1] = 0.9*p[2*j+1] + 0.1*c[1];
			}
			
			fprintf(fp, "newpath");
			for(j = 0; j < nf; ++j){
				if(0 == j){
					fprintf(fp, " %g %g moveto", p[2*j+0], p[2*j+1]);
				}else{
					fprintf(fp, " %g %g lineto", p[2*j+0], p[2*j+1]);
				}
			}
			fprintf(fp, " closepath stroke\n");
		}
		free(iface);
		free(p);
	}
	
	fprintf(fp, "showpage\n");
	
	fclose(fp);
	return 0;
}
const double *POFF2Mesh_Lattice(const POFF2Mesh mesh){
	if(NULL == mesh){ return NULL; }
	return &(mesh->L[0]);
}
int POFF2Mesh_NumVertices(const POFF2Mesh mesh){
	if(NULL == mesh){ return 0; }
	return mesh->n_verts;
}
int POFF2Mesh_NumEdges(const POFF2Mesh mesh){
	if(NULL == mesh){ return 0; }
	return mesh->n_halfs/2;
}
int POFF2Mesh_NumFaces(const POFF2Mesh mesh){
	if(NULL == mesh){ return 0; }
	return mesh->n_faces;
}
int POFF2Mesh_MaxValence(const POFF2Mesh mesh){
	if(NULL == mesh){ return 0; }
	return mesh->max_valence;
}
int POFF2Mesh_MaxDualValence(const POFF2Mesh mesh){
	if(NULL == mesh){ return 0; }
	return mesh->max_dual_valence;
}
int POFF2Mesh_GetVertex(const POFF2Mesh mesh, int ivert, double *p){
	if(NULL == mesh){ return -1; }
	if(ivert < 0 || ivert >= mesh->n_verts){ return -2; }
	if(NULL == p){ return -3; }
	{
		const int h = mesh->vert[ivert];
		const int vo = mesh->half[h].from;
		p[0] = mesh->vertorg[vo].x;
		p[1] = mesh->vertorg[vo].y;
	}
	return 0;
}
int POFF2Mesh_GetEdgeVertices(
	const POFF2Mesh mesh, int iedge,
	POFF2Mesh_Index ivert[2]
){
	int i;
	int ivo[2];
	const int h0 = 2*iedge;
	
	if(NULL == mesh){ return -1; }
	if(iedge < 0 || iedge >= mesh->n_halfs/2){ return -2; }
	if(NULL == ivert){ return -3; }
	
	ivo[0] = mesh->half[h0].from;
	ivo[1] = mesh->half[mesh->half[h0].next].from;
	for(i = 0; i < 2; ++i){
		ivert[i].idx = mesh->vmap[2*ivo[i]];
		ivert[i].off[0] = -offmap[ mesh->vmap[2*ivo[i]+1] & 0x3      ];
		ivert[i].off[1] = -offmap[(mesh->vmap[2*ivo[i]+1] & 0xC) >> 2];
	}
	ivert[0].sgn =  1;
	ivert[1].sgn = -1;
	
	return 0;
}

int POFF2Mesh_GetFaceEdges(
	const POFF2Mesh mesh, int iface,
	POFF2Mesh_Index *iedge
){
	int i, h, h0;
	
	if(NULL == mesh){ return -1; }
	if(iface < 0 || iface >= mesh->n_faces){ return -2; }
	if(NULL == iedge){ return -3; }
	
	h0 = mesh->face[iface];
	h = h0;
	i = 0;
	do{
		const int hopp = 1 ^ h;
		iedge[i].idx = (h >> 1);
		if(h & 0x1){
			iedge[i].off[0] = (int)(mesh->half[h].flag & 0x1) - (int)(mesh->half[hopp].flag & 0x1);
			iedge[i].off[1] = (int)((mesh->half[h].flag & 0x2) >> 1) - (int)((mesh->half[hopp].flag & 0x2) >> 1);
		}else{
			iedge[i].off[0] = 0;
			iedge[i].off[1] = 0;
		}
		iedge[i].sgn = ((h & 0x1) ? -1 : 1);
		h = mesh->half[h].next;
		++i;
	}while(h != mesh->face[iface]);
	
	return i;
}

int POFF2Mesh_GetFaceVertexEdges(
	const POFF2Mesh mesh, int iface, int ivert,
	POFF2Mesh_Index iedge[2]
){
	int found = 0;
	if(NULL == mesh){ return -1; }
	if(iface < 0 || iface >= mesh->n_faces){ return -2; }
	if(ivert < 0 || ivert >= mesh->n_verts){ return -3; }
	if(NULL == iedge){ return -4; }
	
	{
		const int h0 = mesh->face[iface];
		int h[2] = {h0, -1};
		do{
			h[1] = mesh->half[h[0]].next;
			
			if(mesh->vmap[2*mesh->half[h[1]].from] == ivert){
				int i = 0;
				for(i = 0; i < 2; ++i){
					iedge[i].idx = (h[i] >> 1);
					if(h[i] & 0x1){
						iedge[i].off[0] = (int)(mesh->half[h[i]].flag & 0x1) - (int)(mesh->half[1^h[i]].flag & 0x1);
						iedge[i].off[1] = (int)((mesh->half[h[i]].flag & 0x2) >> 1) - (int)((mesh->half[1^h[i]].flag & 0x2) >> 1);
					}else{
						iedge[i].off[0] = 0;
						iedge[i].off[1] = 0;
					}
					iedge[i].sgn = ((h[i] & 0x1) ? -1 : 1);
				}
				found = 1;
				break;
			}
			h[0] = h[1];
		}while(h[0] != h0);
	}
	return !found;
}

POFF2Mesh_Index POFF2Mesh_LocatePoint(
	const POFF2Mesh mesh, const double p[2],
	POFF2Mesh_Index *iface0
){
	double q[2];
	POFF2Mesh_Index start;
	POFF2Mesh_Index ret;
	ret.idx = -1;
	ret.off[0] = 0;
	ret.off[1] = 0;
	ret.sgn = 0;

	if(NULL == mesh){ return ret; }
	if(NULL == p){ return ret; }

	/* start at iface0 if provided */
	if(NULL != iface0){
		start.idx = iface0->idx;
		start.off[0] = iface0->off[0];
		start.off[1] = iface0->off[1];
		q[0] = p[0] - (mesh->L[0] * (double)start.off[0] + mesh->L[2] * (double)start.off[1]);
		q[1] = p[1] - (mesh->L[1] * (double)start.off[0] + mesh->L[3] * (double)start.off[1]);
	}else{
		start.idx = 0;
		start.off[0] = 0;
		start.off[1] = 0;
		q[0] = p[0];
		q[1] = p[1];
	}
	
	/* Now perform hopping until we hit the face containing it */
	{
		int safety = 2*mesh->n_faces;
		int h = mesh->face[start.idx];
		do{
			int is_inside = 1;
			int fh = h;
			/* Loop around current face */
			do{
				const int fhn = mesh->half[fh].next;
				const int ivo0 = mesh->half[fh].from;
				const int ivo1 = mesh->half[fhn].from;
				const double v0[2] = {
					mesh->vertorg[ivo0].x,
					mesh->vertorg[ivo0].y
				};
				const double v1[2] = {
					mesh->vertorg[ivo1].x,
					mesh->vertorg[ivo1].y
				};
				if(geom_orient2d(v0, v1, q) < 0){
					/* Flip over this halfedge now */
					const int fhopp = 1 ^ fh;
					if(mesh->half[fh].flag != mesh->half[fhopp].flag){
						const int off0 = (int)(mesh->half[fh].flag & 0x1) - (int)(mesh->half[fhopp].flag & 0x1);
						const int off1 = (int)((mesh->half[fh].flag & 0x2) >> 1) - (int)((mesh->half[fhopp].flag & 0x2) >> 1);
						ret.off[0] += off0;
						ret.off[1] += off1;
						q[0] -= (mesh->L[0] * (double)off0 + mesh->L[2] * (double)off1);
						q[1] -= (mesh->L[1] * (double)off0 + mesh->L[3] * (double)off1);
					}
					
					is_inside = 0;
					break;
				}
				fh = fhn;
			}while(fh != h);
			if(is_inside){
				ret.idx = mesh->half[h].face;
				break;
			}
		}while(safety --> 0);
	}
	return ret;
}
int POFF2Mesh_GetFaceCenter(
	const POFF2Mesh mesh, int iface,
	double c[2]
){
	if(NULL == mesh){ return -1; }
	if(iface < 0 || iface >= mesh->n_faces){ return -2; }
	if(NULL == c){ return -3; }
	c[0] = mesh->center[2*iface+0];
	c[1] = mesh->center[2*iface+1];
	return 0;
}
int POFF2Mesh_GetVertexFaces(
	const POFF2Mesh mesh, int ivert,
	POFF2Mesh_Index *iface
){
	int i = 0;
	if(NULL == mesh){ return -1; }
	if(ivert < 0 || ivert >= mesh->n_verts){ return -2; }
	if(NULL == iface){ return -3; }
	
	{
		int curoff[2] = { 0, 0 };
		const int h0 = mesh->vert[ivert];
		int h = h0;
		do{
			iface[i].idx = mesh->half[h].face;
			iface[i].off[0] = curoff[0];
			iface[i].off[1] = curoff[1];
			
			/* Go to next face in CCW order */
			/* spin round to prev */
			{
				const int hsave = h;
				do{
					h = mesh->half[h].next;
				}while(mesh->half[h].next != hsave);
			}
			/* Flip over to opposite */
			{
				const int hopp = 1 ^ h;
				if(mesh->half[h].flag != mesh->half[hopp].flag){
					const int off0 = (int)(mesh->half[h].flag & 0x1) - (int)(mesh->half[hopp].flag & 0x1);
					const int off1 = (int)((mesh->half[h].flag & 0x2) >> 1) - (int)((mesh->half[hopp].flag & 0x2) >> 1);
					curoff[0] += off0;
					curoff[1] += off1;
				}
				h = hopp;
			}
			
			i++;
		}while(h != h0 && i < mesh->max_valence);
	}
	
	return i;
}

int POFF2Mesh_GetEdgeFaces(
	const POFF2Mesh mesh, int iedge,
	POFF2Mesh_Index iface[2]
){
	if(NULL == mesh){ return -1; }
	if(iedge < 0 || iedge >= mesh->n_halfs/2){ return -2; }
	if(NULL == iface){ return -3; }
	
	{
		const int h = 2*iedge;
		const int hopp = h+1;
		iface[0].idx = mesh->half[h].face;
		iface[0].off[0] = 0;
		iface[0].off[1] = 0;
		iface[0].sgn = -1;
		iface[1].idx = mesh->half[hopp].face;
		iface[1].off[0] = (int)(mesh->half[h].flag & 0x1) - (int)(mesh->half[hopp].flag & 0x1);
		iface[1].off[1] = (int)((mesh->half[h].flag & 0x2) >> 1) - (int)((mesh->half[hopp].flag & 0x2) >> 1);
		iface[1].sgn = -1;
	}
	return 0;
}
int POFF2Mesh_GetFaceVertices(
	const POFF2Mesh mesh, int iface,
	POFF2Mesh_Index *ivert
){
	if(NULL == mesh){ return -1; }
	if(iface < 0 || iface >= mesh->n_faces){ return -2; }
	if(NULL == ivert){ return -3; }
	
	{
		int i;
		const int h0 = mesh->face[iface];
		int h = h0;
		do{
			const int ivo = mesh->half[h].from;
			const int flags = mesh->vmap[2*ivo+1];
			ivert[i].idx = mesh->vmap[2*ivo];
			ivert[i].off[0] = -offmap[ flags & 0x3      ];
			ivert[i].off[1] = -offmap[(flags & 0xC) >> 2];

			h = mesh->half[h].next;
			++i;
		}while(h != h0);
		return i;
	}
}
