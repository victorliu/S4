#include <S4r/PeriodicMesh.hpp>
using namespace S4r;

static double modf_pos(double x, int *ipart){
	double id;
	double ret = modf(x, &id);
	*ipart = (int)id;
	if(ret < 0){ ret += 1; }
	(*ipart)--;
	return ret;
}

S4r::LatticeGridRect::LatticeGridRect(
	const Vec2 &u, const Vec2 &v, const size_t ncells[2]
):PeriodicMesh(u,v){
	// Assert that u.v =~ 0
	
	n[0] = ncells[0];
	n[1] = ncells[1];
	du = u / (double)n[0];
	dv = v / (double)n[1];
	hu = 0.5*du;
	hv = 0.5*dv;
}

S4r::LatticeGridRect::~LatticeGridRect(){
}

size_t S4r::LatticeGridRect::idx(size_t ix, size_t iy, size_t iedge) const{
	ix = (ix+n[0]) % n[0];
	iy = (iy+n[1]) % n[1];
	//iedge &= 1;
	return 2*(ix + iy*n[0]) + iedge;
}

size_t S4r::LatticeGridRect::NumVertices() const{ return n[0]*n[1]; }
size_t S4r::LatticeGridRect::NumEdges() const{ return 2*n[0]*n[1]; }
size_t S4r::LatticeGridRect::NumFaces() const{ return n[0]*n[1]; }

Vec2 S4r::LatticeGridRect::Vertex(size_t i) const{
	const int iy = i / n[0];
	const int ix = i % n[0];
	const int nx = n[0];
	const int ny = n[1];
	//const int offx = n[0]%2;
	//const int offy = n[1]%2;
	double tx = (double)(-nx + 2*ix) / (double)(2*nx);
	double ty = (double)(-ny + 2*iy) / (double)(2*ny);
	// This guarantees one vertex is at the origin.
	//double tx = (double)(offx - nx + 2*ix) / (double)(2*nx);
	//double ty = (double)(offy - ny + 2*iy) / (double)(2*ny);
	return tx * Lr.col(0) + ty * Lr.col(1);
}
Vec2 S4r::LatticeGridRect::Edge(size_t i) const{
	if(0 == i%2){ // u-directed
		return Lr.col(0) / (double)n[0];
	}else{
		return Lr.col(1) / (double)n[1];
	}
}

// Control polygons
void S4r::LatticeGridRect::VertexControlPolygon(size_t i, ConvexPolygon *poly) const{
	poly->v.clear();
	poly->v.resize(4);
	poly->v[0] = -hu-hv;
	poly->v[1] =  hu-hv;
	poly->v[2] =  hu+hv;
	poly->v[3] = -hu+hv;
	poly->offset = Vertex(i);
}
void S4r::LatticeGridRect::EdgeControlPolygon(size_t iedge, ConvexPolygon *self,
                        size_t *e1, ConvexPolygon *wedge1,
                        size_t *e2, ConvexPolygon *wedge2
) const{
	size_t iy = iedge / (2*n[0]);
	size_t ix = iedge % (2*n[0]);
	self->v.clear();   self->v.resize(4);
	wedge1->v.clear(); wedge1->v.resize(4);
	wedge2->v.clear(); wedge2->v.resize(4);
	if(0 == iedge%2){ // u-directed
		/*
		self->v[0] = Vec2(0,0);
		self->v[1] = hu-hv;
		self->v[2] = du;
		self->v[3] = hu+hv;
		*/
		self->v[0] = -hv;
		self->v[1] = du-hv;
		self->v[2] = du+hv;
		self->v[3] = hv;
		*e1 = idx(ix,iy,1);
		wedge1->v[0] = Vec2(0,0);
		wedge1->v[1] = hu;
		wedge1->v[2] = hu+hv;
		wedge1->v[3] = hv;
		*e2 = idx(ix+1,iy-1,1);
		wedge2->v[0] = du;
		wedge2->v[1] = hu;
		wedge2->v[2] = hu-hv;
		wedge2->v[3] = du-hv;
	}else{
		/*
		self->v[0] = Vec2(0,0);
		self->v[1] = hu+hv;
		self->v[2] = dv;
		self->v[3] = hv-hu;
		*/
		self->v[0] = hu;
		self->v[1] = dv+hu;
		self->v[2] = dv-hu;
		self->v[3] = -hu;
		*e1 = idx(ix-1,iy,0);
		wedge1->v[0] = Vec2(0,0);
		wedge1->v[1] = hv;
		wedge1->v[2] = hv-hu;
		wedge1->v[3] = -hu;
		*e2 = idx(ix,iy+1,0);
		wedge2->v[0] = dv;
		wedge2->v[1] = hv;
		wedge2->v[2] = hv+hu;
		wedge2->v[3] = du+hu;
	}
	self->offset = Vertex(iedge/2);
	wedge1->offset = self->offset;
	wedge2->offset = self->offset;
}
void S4r::LatticeGridRect::FaceControlPolygon(size_t i, ConvexPolygon *poly) const{
	poly->v.clear();
	poly->v.resize(4);
	poly->v[0] = Vec2(0,0);
	poly->v[1] = du;
	poly->v[2] = du+dv;
	poly->v[3] = dv;
	poly->offset = Vertex(i);
}

double S4r::LatticeGridRect::EdgeLength(size_t i) const{
	if(0 == i%2){ // u-directed
		return du.norm();
	}else{
		return dv.norm();
	}
}
double S4r::LatticeGridRect::FaceArea(size_t i) const{
	return 1. / (double)(n[0]*n[1]);
}

double S4r::LatticeGridRect::VertexDualArea(size_t i) const{
	return 1. / (double)(n[0]*n[1]);
}
double S4r::LatticeGridRect::EdgeDualLength(size_t i) const{
	if(0 == i%2){ // u-directed
		return dv.norm();
	}else{
		return du.norm();
	}
}

//// Utility functions

PeriodicMesh::PeriodicIndex S4r::LatticeGridRect::ContainingFace(const Vec2 &r, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex ret;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &ret.off[0]);
	double ffy = modf_pos(q[1], &ret.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	ret.idx = ix+iy*n[0];
	return ret;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridRect::VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex pi, ci;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &pi.off[0]);
	double ffy = modf_pos(q[1], &pi.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	pi.idx = ix+iy*n[0];
	
	int off[2] = {0,0};
	size_t cx = ix+1, cy = iy+1;
	if(cx >= n[0]){ cx = 0; off[0] = 1; }
	if(cy >= n[1]){ cy = 0; off[1] = 1; }
	
	coeffs.resize(4);
	
	ci.idx = ix+iy*n[0];
	ci.off[0] = pi.off[0];
	ci.off[1] = pi.off[1];
	coeffs[0] = VertexCoeff(ci, (1.-ffx)*(1.-ffy));
	
	ci.idx = cx+iy*n[0];
	ci.off[0] = pi.off[0]+off[0];
	ci.off[1] = pi.off[1];
	coeffs[1] = VertexCoeff(ci, (   ffx)*(1.-ffy));
	
	ci.idx = ix+cy*n[0];
	ci.off[0] = pi.off[0];
	ci.off[1] = pi.off[1]+off[1];
	coeffs[2] = VertexCoeff(ci, (1.-ffx)*(   ffy));
	
	ci.idx = cx+cy*n[0];
	ci.off[0] = pi.off[0]+off[0];
	ci.off[1] = pi.off[1]+off[1];
	coeffs[3] = VertexCoeff(ci, (   ffx)*(   ffy));
	return pi;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridRect::EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex pi, ci;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &pi.off[0]);
	double ffy = modf_pos(q[1], &pi.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	pi.idx = ix+iy*n[0];
	
	int off[2] = {0,0};
	size_t cx = ix+1, cy = iy+1;
	if(cx >= n[0]){ cx = 0; off[0] = 1; }
	if(cy >= n[1]){ cy = 0; off[1] = 1; }
	
	coeffs.resize(4);
	
	ci.idx = idx(ix,iy,0);
	ci.off[0] = pi.off[0];
	ci.off[1] = pi.off[1];
	coeffs[0] = EdgeCoeff(ci, Vec2(1.-ffy, 0));
	
	ci.idx = idx(ix,iy,1);
	ci.off[0] = pi.off[0]+off[0];
	ci.off[1] = pi.off[1];
	coeffs[1] = EdgeCoeff(ci, Vec2(0, 1.-ffx));
	
	ci.idx = idx(ix,cy,0);
	ci.off[0] = pi.off[0];
	ci.off[1] = pi.off[1]+off[1];
	coeffs[2] = EdgeCoeff(ci, Vec2(ffy, 0));
	
	ci.idx = idx(cx,iy,1);
	ci.off[0] = pi.off[0]+off[0];
	ci.off[1] = pi.off[1]+off[1];
	coeffs[3] = EdgeCoeff(ci, Vec2(0, ffx));
	return pi;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridRect::FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	// Simple piecewise constant interpolation
	PeriodicIndex iface = ContainingFace(r, faceguess);
	coeffs.resize(1);
	coeffs[0] = FaceCoeff(iface, 1.);
	return iface;
}

void S4r::LatticeGridRect::BuildMatrices(const doublecomplex phi[2]){
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
}
