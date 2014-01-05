#ifndef PERIODIC_OFF_2_H_INCLUDED
#define PERIODIC_OFF_2_H_INCLUDED

#ifdef __cplusplus
extern "C" {
 #endif

typedef struct POFF2_Vertex_{
	double x, y;
	int flags;
} POFF2_Vertex;

typedef struct POFF2_Face_{
	int n_verts;
	int *vert;
	int flags;
} POFF2_Face;

typedef struct POFF2_Matching_{
	// First element of each array refers to the first edge of the matching
	int face[2];
	int vert1[2];
	int vert2[2];
	int flags;
	// This is output as
	//   face[0] vert1[0] vert2[0] face[1] vert1[1] vert2[1]
	// where vert2[0]-vert1[0] is the edge vector of face[0], and should
	// be aligned with vert2[1]-vert1[1], the edge vector of face[1];
	// these edge vectors should point in the same direction.
	// Therefore, the order of the vertex indexes within the arrays matters.
} POFF2_Matching;

typedef struct POFF2_{
	double L[4];
	int n_verts, n_faces, n_matchings;
	POFF2_Vertex *vert;
	POFF2_Face *face;
	POFF2_Matching *matching;
} POFF2;

// Loads a POFF2 file. Assumes *off is not initialized at all.
// Returns:
//    0 on success.
//    1 if off == NULL.
//    2 if filename == NULL.
//   -1 if the file could not be opened.
//   -2 if the header signature "PeriodicOFF2" was not found.
//   -3 if the file is missing data (more fields were expected before EOF)
int POFF2_LoadFromFile(POFF2 *off, const char *filename);

// Deallocate a POFF2 object.
void POFF2_Destroy(POFF2 *off);
int POFF2_Clone(const POFF2 *off, POFF2 *cpy);

// Saves a POFF2 file.
// Returns:
//    0 on success.
//    1 if off == NULL.
//    2 if filename == NULL.
//   -1 if the file could not be opened.
//   -3 if an error occurred during writing of the file.
int POFF2_SaveToFile(const POFF2 *off, const char *filename);

// Partially allocates a POFF2 structure.
// (The vertex arrays in each face are not initialized)
// nverts is the number of vertexes in the mesh.
// nfaces is the number of faces in the mesh.
// nmatchings is the number of boundary edge matchings in the mesh.
void POFF2_Initialize(POFF2 *off, int nverts, int nfaces, int nmatchings);

// Allocates the vertex array within a POFF2 face.
// nverts is the number of vertexes of the face.
void POFF2_InitializeFace(POFF2_Face *face, int nverts);

/* Opaque pointer to mesh object */
typedef struct POFF2Mesh_ *POFF2Mesh;
/* All indexes to elements can potentially have periodic
 * offsets associated with it.
 */
typedef struct{
	int idx;
	int off[2];
	int sgn;
} POFF2Mesh_Index;

/* Create a Mesh object from the raw POFF2 input.
 * The lattice on which the POFF2 is defined should be specified,
 * and the actual matchings are checked against this lattice to
 * ensure consistency.
 */
POFF2Mesh POFF2Mesh_Create(const POFF2 *off, double tol);
void POFF2Mesh_Destroy(const POFF2Mesh mesh);
int POFF2Mesh_SavePostscript(const POFF2Mesh mesh, const char *filename);

const double *POFF2Mesh_Lattice(const POFF2Mesh mesh);
/* Return vertex, edge, and face counts. */
int POFF2Mesh_NumVertices(const POFF2Mesh mesh);
int POFF2Mesh_NumEdges(const POFF2Mesh mesh);
int POFF2Mesh_NumFaces(const POFF2Mesh mesh);
/* Returns max vertex valence in the mesh */
int POFF2Mesh_MaxValence(const POFF2Mesh mesh);
/* Returns max vertex valence in the dual mesh.
 * This is equivalent to the maximum number of edges
 * a face can have.
 */
int POFF2Mesh_MaxDualValence(const POFF2Mesh mesh);
int POFF2Mesh_GetVertex(const POFF2Mesh mesh, int ivert, double *p);
/* Returns the endpoint vertices of an edge */
int POFF2Mesh_GetEdgeVertices(
	const POFF2Mesh mesh, int iedge,
	POFF2Mesh_Index ivert[2]
);
/* Returns the edges around a face */
int POFF2Mesh_GetFaceEdges(
	const POFF2Mesh mesh, int iface,
	POFF2Mesh_Index *iedge
);
/* Returns the vertices around a face */
int POFF2Mesh_GetFaceVertices(
	const POFF2Mesh mesh, int iface,
	POFF2Mesh_Index *ivert
);
/* Returns the circumcenter of a face */
int POFF2Mesh_GetFaceCenter(
	const POFF2Mesh mesh, int iface,
	double c[2]
);
/* Returns the face containing a point */
POFF2Mesh_Index POFF2Mesh_LocatePoint(
	const POFF2Mesh mesh, const double p[2],
	POFF2Mesh_Index *iface0
);
/* Returns the faces around a vertex */
int POFF2Mesh_GetVertexFaces(
	const POFF2Mesh mesh, int ivert,
	POFF2Mesh_Index *iface
);
/* Returns the two edges adjacent to a given vertex and face. */
int POFF2Mesh_GetFaceVertexEdges(
	const POFF2Mesh mesh, int iface, int ivert,
	POFF2Mesh_Index iedge[2]
);
/* Returns the edges around a vertex */
int POFF2Mesh_GetVertexEdges(
	const POFF2Mesh mesh, int ivert,
	POFF2Mesh_Index *iedge
);
/* Returns the faces around an edge. The face on the left is returned first. */
int POFF2Mesh_GetEdgeFaces(
	const POFF2Mesh mesh, int iedge,
	POFF2Mesh_Index iface[2]
);
/* Returns the faces neighboring a face */
int POFF2Mesh_GetFaceFaces(
	const POFF2Mesh mesh, int iface,
	POFF2Mesh_Index *ineigh
);

#ifdef __cplusplus
}
#endif

#endif // PERIODIC_OFF_2_H_INCLUDED
