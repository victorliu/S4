#include <S4r/PeriodicMesh.hpp>
using namespace S4r;

S4r::POFF2Mesh::POFF2Mesh(
	const Vec2 &u, const Vec2 &v, const char *filename
):PeriodicMesh(u,v){
	int ret;
	POFF2 off;
	ret = POFF2_LoadFromFile(&off, filename);
	if(0 != ret){ throw std::string("Could not load POFF2 mesh file"); }
	
	M = POFF2Mesh_Create(&off, 1e-6);
	if(NULL == M){ throw std::string("Error building mesh"); }
	
	POFF2_Destroy(&off);
}

S4r::POFF2Mesh::~POFF2Mesh(){
	POFF2Mesh_Destroy(M);
}

size_t S4r::POFF2Mesh::NumVertices() const{ return POFF2Mesh_NumVertices(M); }
size_t S4r::POFF2Mesh::NumEdges() const{ return POFF2Mesh_NumEdges(M); }
size_t S4r::POFF2Mesh::NumFaces() const{ return POFF2Mesh_NumFaces(M); }

Vec2 S4r::POFF2Mesh::Vertex(size_t i) const{
	double p[2];
	POFF2Mesh_GetVertex(M, i, p);
	return Vec2(p[0], p[1]);
}
Vec2 S4r::POFF2Mesh::Edge(size_t i) const{
	POFF2Mesh_Index ivert[2];
	POFF2Mesh_GetEdgeVertices(M, i, ivert);
	
	const double *L = POFF2Mesh_Lattice(M);
	double p[4];
	POFF2Mesh_GetVertex(M, ivert[0].idx, &p[0]);
	p[0] += (L[0] * ivert[0].off[0] + L[2] * ivert[0].off[1]);
	p[1] += (L[1] * ivert[0].off[0] + L[3] * ivert[0].off[1]);
	POFF2Mesh_GetVertex(M, ivert[1].idx, &p[2]);
	p[2] += (L[0] * ivert[1].off[0] + L[2] * ivert[1].off[1]);
	p[3] += (L[1] * ivert[1].off[0] + L[3] * ivert[1].off[1]);
	return Vec2(p[2]-p[0], p[3]-p[1]);
}

// Control polygons
void S4r::POFF2Mesh::VertexControlPolygon(size_t i, ConvexPolygon *poly) const{
	POFF2Mesh_Index *iface = (POFF2Mesh_Index*)malloc(sizeof(POFF2Mesh_Index) * POFF2Mesh_MaxValence(M));
	
	const int nf = POFF2Mesh_GetVertexFaces(M, i, iface);
	double *p = (double*)malloc(sizeof(double) * 2*nf);
	const double *L = POFF2Mesh_Lattice(M);
	poly->v.clear();
	poly->v.resize(nf);
	for(int i = 0; i < nf; ++i){
		double p[2];
		POFF2Mesh_GetFaceCenter(M, iface[i].idx, p);
		p[2*i+0] += (L[0] * iface[i].off[0] + L[2] * iface[i].off[1]);
		p[2*i+1] += (L[1] * iface[i].off[0] + L[3] * iface[i].off[1]);
		poly->v[i] = Vec2(p[0], p[1]);
	}
	poly->offset = Vec2(0,0);
	free(p);
	free(iface);
}
void S4r::POFF2Mesh::EdgeControlPolygon(
	size_t iedge, ConvexPolygon *self,
	size_t   *e1, ConvexPolygon *wedge1,
	size_t   *e2, ConvexPolygon *wedge2
) const{
	const int i = iedge;
	POFF2Mesh_Index iface[2];
	POFF2Mesh_Index ivert[2];
	double p[8];
	const double *L = POFF2Mesh_Lattice(M);
	
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
	
	self->v.clear(); self->v.resize(4);
	self->v[0] = Vec2(p[0],p[1]);
	self->v[1] = Vec2(p[2],p[3]);
	self->v[2] = Vec2(p[4],p[5]);
	self->v[3] = Vec2(p[6],p[7]);
	self->offset = Vec2(0,0);
	
	{
		double v[2];
		POFF2Mesh_Index ie[2];
		POFF2Mesh_GetFaceVertexEdges(M, iface[0].idx, ivert[0].idx, ie);
		wedge1->v.clear(); wedge1->v.resize(4);
		*e1 = ie[0].idx;
		wedge1->v[0] = Vec2(p[0], p[1]);
		wedge1->v[1] = Vec2(0.5*p[0]+0.5*p[4], 0.5*p[1]+0.5*p[5]);
		wedge1->v[2] = Vec2(p[6], p[7]);
		
		POFF2Mesh_GetEdgeVertices(M, ie[0].idx, ivert);
		POFF2Mesh_GetVertex(M, ivert[0].idx, v);
		v[0] += (L[0] * (ivert[0].off[0]+iface[0].off[0]) + L[2] * (ivert[0].off[1]+iface[0].off[1]));
		v[1] += (L[1] * (ivert[0].off[0]+iface[0].off[0]) + L[3] * (ivert[0].off[1]+iface[0].off[1]));
		
		wedge1->v[3] = Vec2(0.5*p[0]+0.5*v[0], 0.5*p[1]+0.5*v[1]);
		wedge1->offset = Vec2(0,0);
	}
	{
		double v[2];
		POFF2Mesh_Index ie[2];
		POFF2Mesh_GetFaceVertexEdges(M, iface[1].idx, ivert[1].idx, ie);
		wedge1->v.clear(); wedge1->v.resize(4);
		*e2 = ie[0].idx;
		wedge2->v[0] = Vec2(p[4], p[5]);
		wedge2->v[1] = Vec2(0.5*p[0]+0.5*p[4], 0.5*p[1]+0.5*p[5]);
		wedge2->v[2] = Vec2(p[2], p[3]);
		
		POFF2Mesh_GetEdgeVertices(M, ie[0].idx, ivert);
		POFF2Mesh_GetVertex(M, ivert[0].idx, v);
		v[0] += (L[0] * (ivert[0].off[0]+iface[0].off[0]) + L[2] * (ivert[0].off[1]+iface[0].off[1]));
		v[1] += (L[1] * (ivert[0].off[0]+iface[0].off[0]) + L[3] * (ivert[0].off[1]+iface[0].off[1]));
		
		wedge1->v[3] = Vec2(0.5*p[4]+0.5*v[0], 0.5*p[5]+0.5*v[1]);
		wedge2->offset = Vec2(0,0);
	}
}
void S4r::POFF2Mesh::FaceControlPolygon(size_t i, ConvexPolygon *poly) const{
	const double *L = POFF2Mesh_Lattice(M);
	POFF2Mesh_Index *ivert = (POFF2Mesh_Index*)malloc(sizeof(POFF2Mesh_Index) * POFF2Mesh_MaxDualValence(M));
	poly->v.resize(POFF2Mesh_GetFaceVertices(M, i, ivert));
	for(size_t i = 0; i < poly->v.size(); ++i){
		double p[2];
		POFF2Mesh_GetVertex(M, ivert[i].idx, p);
		p[0] += (L[0] * ivert[i].off[0] + L[2] * ivert[i].off[1]);
		p[1] += (L[1] * ivert[i].off[0] + L[3] * ivert[i].off[1]);
		poly->v[i] = Vec2(p[0], p[1]);
	}
	poly->offset = Vec2(0,0);
	free(ivert);
}

double S4r::POFF2Mesh::EdgeLength(size_t i) const{
	POFF2Mesh_Index ivert[2];
	double p[4];
	const double *L = POFF2Mesh_Lattice(M);
	
	POFF2Mesh_GetEdgeVertices(M, i, ivert);
	POFF2Mesh_GetVertex(M, ivert[0].idx, &p[0]);
	p[0] += (L[0] * ivert[0].off[0] + L[2] * ivert[0].off[1]);
	p[1] += (L[1] * ivert[0].off[0] + L[3] * ivert[0].off[1]);
	POFF2Mesh_GetVertex(M, ivert[1].idx, &p[2]);
	p[2] += (L[0] * ivert[1].off[0] + L[2] * ivert[1].off[1]);
	p[3] += (L[1] * ivert[1].off[0] + L[3] * ivert[1].off[1]);
	return hypot(p[2]-p[0], p[3]-p[1]);
}
double S4r::POFF2Mesh::FaceArea(size_t i) const{
	ConvexPolygon poly;
	FaceControlPolygon(i, &poly);
	return poly.Area();
}

double S4r::POFF2Mesh::VertexDualArea(size_t i) const{
	const double *L = POFF2Mesh_Lattice(M);
	POFF2Mesh_Index *iface = (POFF2Mesh_Index*)malloc(sizeof(POFF2Mesh_Index) * POFF2Mesh_MaxValence(M));
	const int nf = POFF2Mesh_GetVertexFaces(M, i, iface);
	double *p = (double*)malloc(sizeof(double) * 2*nf);
	for(int j = 0; j < nf; ++j){
		POFF2Mesh_GetFaceCenter(M, iface[j].idx, &p[2*j]);
		p[2*j+0] += (L[0] * iface[j].off[0] + L[2] * iface[j].off[1]);
		p[2*j+1] += (L[1] * iface[j].off[0] + L[3] * iface[j].off[1]);
	}
	ConvexPolygon poly(nf, p);
	free(p);
	free(iface);
	return poly.Area();
}
double S4r::POFF2Mesh::EdgeDualLength(size_t i) const{
	POFF2Mesh_Index iface[2];
	POFF2Mesh_Index ivert[2];
	double p[8];
	const double *L = POFF2Mesh_Lattice(M);
	
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
	
	// Compute signed dual length
	// We want to compute Cross[(p[2,3] - p[6,7]), p[4,5] - p[0,1]) / Len[p[4,5] - p[0,1]]
	p[4] -= p[0]; p[5] -= p[1];
	p[2] -= p[6]; p[3] -= p[7];
	return (p[2]*p[5]-p[3]*p[4]) / hypot(p[4],p[5]);
}

//// Utility functions

PeriodicMesh::PeriodicIndex S4r::POFF2Mesh::ContainingFace(const Vec2 &r, PeriodicMesh::PeriodicIndex *faceguess) const{
	const double p[2] = { r[0], r[1] };
	POFF2Mesh_Index guess;
	if(NULL == faceguess){
		guess.idx = 0;
		guess.off[0] = 0;
		guess.off[1] = 0;
	}else{
		guess.idx = faceguess->idx;
		guess.off[0] = faceguess->off[0];
		guess.off[1] = faceguess->off[1];
	}
	POFF2Mesh_Index iface = POFF2Mesh_LocatePoint(M, p, &guess);
	PeriodicIndex ret;
	ret.idx = iface.idx;
	ret.off[0] = iface.off[0];
	ret.off[1] = iface.off[1];
	return ret;
}
PeriodicMesh::PeriodicIndex S4r::POFF2Mesh::VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	// For now assume that we have only triangles and rectangles
	PeriodicMesh::PeriodicIndex iface = ContainingFace(r, faceguess);
	int maxnvert = POFF2Mesh_MaxDualValence(M);
	double *v = new double[2*maxnvert];
	double *c = new double[maxnvert];
	const double *L = POFF2Mesh_Lattice(M);
	
	POFF2Mesh_Index *ivert = new POFF2Mesh_Index[maxnvert];
	int nvert = POFF2Mesh_GetFaceVertices(M, iface.idx, ivert);
	
	PeriodicIndex pi;
	
	coeffs.resize(nvert);
	for(int i = 0; i < nvert; ++i){
		POFF2Mesh_GetVertex(M, ivert[i].idx, &v[2*i+0]);
	}
	const double p[2] = {
		r[0] - (L[0] * iface.off[0] + L[2] * iface.off[1]),
		r[1] - (L[1] * iface.off[0] + L[3] * iface.off[1])
	};
//	geom_interpolate_polygon(nvert, v, p, c, NULL);
	for(int i = 0; i < nvert; ++i){
		pi.idx = ivert[i].idx;
		pi.off[0] = ivert[i].off[0] + iface.off[0];
		pi.off[1] = ivert[i].off[1] + iface.off[1];
		coeffs[i] = VertexCoeff(pi, c[i]);
	}
	
	delete [] ivert;
	delete [] c;
	delete [] v;
	return iface;
}
PeriodicMesh::PeriodicIndex S4r::POFF2Mesh::EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicMesh::PeriodicIndex iface = ContainingFace(r, faceguess);
	
	return iface;
}
PeriodicMesh::PeriodicIndex S4r::POFF2Mesh::FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	// Simple piecewise constant interpolation
	PeriodicIndex iface = ContainingFace(r, faceguess);
	coeffs.resize(1);
	coeffs[0] = FaceCoeff(iface, 1.);
	return iface;
}

void S4r::POFF2Mesh::BuildMatrices(const doublecomplex phi[2]){
/*
	const size_t N = n[0]*n[1];
	d0.resize(2*N, N);
	std::vector<Triplet> triplst;
	// Each row has two nonzeros
	triplst.reserve(4*N);
	
	for(size_t j = 0; j < n[1]; ++j){
		for(size_t i = 0; i < n[0]; ++i){
			// x directed edge
			if(i+1 < n[0]){
				triplst.push_back(Triplet(
					2*(i+j*n[0])+0,
					i+j*n[0],
					-1.0
				));
				triplst.push_back(Triplet(
					2*(i+j*n[0])+0,
					i+1+j*n[0],
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					2*(i+j*n[0])+0,
					i+j*n[0],
					-std::conj(phi[0])
				));
				triplst.push_back(Triplet(
					2*(i+j*n[0])+0,
					0+j*n[0],
					phi[0]
				));
			}
			// y directed edge
			if(j+1 < n[1]){
				triplst.push_back(Triplet(
					2*(i+j*n[0])+1,
					i+j*n[0],
					-1.0
				));
				triplst.push_back(Triplet(
					2*(i+j*n[0])+1,
					i+(j+1)*n[0],
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					2*(i+j*n[0])+1,
					i+j*n[0],
					-std::conj(phi[1])
				));
				triplst.push_back(Triplet(
					2*(i+j*n[0])+1,
					i+(0)*n[0],
					phi[1]
				));
			}
		}
	}
	d0.setFromTriplets(triplst.begin(), triplst.end());
	
	d1.resize(N, 2*N);
	// Each row has four nonzeros
	triplst.clear();
	triplst.reserve(4*N);
	for(size_t j = 0; j < n[1]; ++j){
		for(size_t i = 0; i < n[0]; ++i){
			if(i+1 < n[0]){
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+j*n[0])+1,
					-1.0
				));
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+1+j*n[0])+1,
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+j*n[0])+1,
					-std::conj(phi[0])
				));
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(0+j*n[0])+1,
					phi[0]
				));
			}
			if(j+1 < n[1]){
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+j*n[0])+0,
					1.0
				));
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+(j+1)*n[0])+0,
					-1.0
				));
			}else{
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+j*n[0])+0,
					std::conj(phi[1])
				));
				triplst.push_back(Triplet(
					i+j*n[0],
					2*(i+(0)*n[0])+0,
					-phi[1]
				));
			}
		}
	}
	d1.setFromTriplets(triplst.begin(), triplst.end());
	*/
}
