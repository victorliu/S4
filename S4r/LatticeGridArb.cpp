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

Vec2 Circumcenter(
	const Vec2 &a,
	const Vec2 &b,
	const Vec2 &c
){
	const Vec2 ba(b-a);
	const Vec2 ca(c-a);

	// Squares of lengths of the edges incident to `a'.
	const double balength = ba.squaredNorm();
	const double calength = ca.squaredNorm();

	// Calculate the denominator of the formulae.
	const double denominator = 0.5 / Orient2D(b, c, a);

	// Calculate offset (from `a') of circumcenter.
	return Vec2(
		(ca[1] * balength - ba[1] * calength) * denominator,
		(ba[0] * calength - ca[0] * balength) * denominator
	);
}

S4r::LatticeGridArb::LatticeGridArb(
	const Vec2 &u, const Vec2 &v, const size_t ncells[2]
):PeriodicMesh(u,v){
	// Assert that u.v < 0, u x v > 0
	
	n[0] = ncells[0];
	n[1] = ncells[1];
	du = u / (double)n[0];
	dv = v / (double)n[1];
	lu = du.norm();
	lv = dv.norm();
	luv = (du+dv).norm();
	
	c1 = Circumcenter(Vec2(0,0), du, du+dv);
	c2 = Circumcenter(Vec2(0,0), du+dv, dv);
	dlu = (c1-(c2-dv)).norm();
	dlv = (c2-(c1-du)).norm();
	dluv = (c2-c1).norm();
}

S4r::LatticeGridArb::~LatticeGridArb(){
}

size_t S4r::LatticeGridArb::idx(size_t ix, size_t iy, size_t iedge) const{
	ix = (ix+n[0]) % n[0];
	iy = (iy+n[1]) % n[1];
	//iedge = iedge % 3;
	return 3*(ix + iy*n[0]) + iedge;
}

size_t S4r::LatticeGridArb::NumVertices() const{ return n[0]*n[1]; }
size_t S4r::LatticeGridArb::NumEdges() const{ return 3*n[0]*n[1]; }
size_t S4r::LatticeGridArb::NumFaces() const{ return 2*n[0]*n[1]; }

Vec2 S4r::LatticeGridArb::Vertex(size_t i) const{
	const int iy = i / n[0];
	const int ix = i % n[0];
	const int nx = n[0];
	const int ny = n[1];
	//const int offx = n[0]%2;
	//const int offy = n[1]%2;
	// This guarantees one vertex is at the origin.
	//double tx = (double)(offx - nx + 2*ix) / (double)(2*nx);
	//double ty = (double)(offy - ny + 2*iy) / (double)(2*ny);
	double tx = (double)(-nx + 2*ix) / (double)(2*nx);
	double ty = (double)(-ny + 2*iy) / (double)(2*ny);
	return tx * Lr.col(0) + ty * Lr.col(1);
}
Vec2 S4r::LatticeGridArb::Edge(size_t i) const{
	size_t iedge = i%3;
	if(0 == iedge){
		return du;
	}else if(1 == iedge){
		return dv;
	}else{
		return du+dv;
	}
}

// Control polygons
void S4r::LatticeGridArb::VertexControlPolygon(size_t i, ConvexPolygon *poly) const{
	poly->v.clear();
	poly->v.resize(6);
	poly->v[0] = c1;
	poly->v[1] = c2;
	poly->v[2] = c1-du;
	poly->v[3] = c2-(du+dv);
	poly->v[4] = c1-(du+dv);
	poly->v[5] = c2-dv;
	poly->offset = Vertex(i);
}
void S4r::LatticeGridArb::EdgeControlPolygon(size_t i, ConvexPolygon *self,
                        size_t *e1, ConvexPolygon *wedge1,
                        size_t *e2, ConvexPolygon *wedge2
) const{
	size_t iedge = i%3;
	size_t icell = i/3;
	size_t iy = icell / n[0];
	size_t ix = icell % n[0];
	self->v.clear();   self->v.resize(4);
	wedge1->v.clear(); wedge1->v.resize(4);
	wedge2->v.clear(); wedge2->v.resize(4);
	if(0 == iedge){ // u-directed
		self->v[0] = Vec2(0,0);
		self->v[1] = c2-dv;
		self->v[2] = du;
		self->v[3] = c1;
		*e1 = idx(ix,iy,2);
		wedge1->v[0] = Vec2(0,0);
		wedge1->v[1] = 0.5*du;
		wedge1->v[2] = c1;
		wedge1->v[3] = 0.5*du+0.5*dv;
		*e2 = idx(ix,iy-1,2);
		wedge2->v[0] = du;
		wedge2->v[1] = 0.5*du;
		wedge2->v[2] = c2-dv;
		wedge2->v[3] = 0.5*du-0.5*dv;
	}else if(1 == iedge){
		self->v[0] = Vec2(0,0);
		self->v[1] = c2;
		self->v[2] = dv;
		self->v[3] = c1-du;
		*e1 = idx(ix-1,iy,0);
		wedge1->v[0] = Vec2(0,0);
		wedge1->v[1] = 0.5*dv;
		wedge1->v[2] = c1-du;
		wedge1->v[3] = -0.5*du;
		*e2 = idx(ix,iy+1,0);
		wedge2->v[0] = dv;
		wedge2->v[1] = 0.5*dv;
		wedge2->v[2] = c2;
		wedge2->v[3] = 0.5*du+0.5*dv;
	}else{
		self->v[0] = Vec2(0,0);
		self->v[1] = c1;
		self->v[2] = du+dv;
		self->v[3] = c2;
		*e1 = idx(ix,iy,1);
		wedge1->v[0] = Vec2(0,0);
		wedge1->v[1] = 0.5*du+0.5*dv;
		wedge1->v[2] = c2;
		wedge1->v[3] = 0.5*dv;
		*e2 = idx(ix+1,iy,1);
		wedge2->v[0] = du+dv;
		wedge2->v[1] = 0.5*du+0.5*dv;
		wedge2->v[2] = c1;
		wedge2->v[3] = 0.5*dv+du;
	}
	self->offset = Vertex(icell);
	wedge1->offset = self->offset;
	wedge2->offset = self->offset;
}
void S4r::LatticeGridArb::FaceControlPolygon(size_t i, ConvexPolygon *poly) const{
	poly->v.clear();
	poly->v.resize(3);
	poly->offset = Vertex(i/2);
	if(0 == i%2){
		poly->v[0] = Vec2(0,0);
		poly->v[1] = du;
		poly->v[2] = du+dv;
	}else{
		poly->v[0] = Vec2(0,0);
		poly->v[1] = du+dv;
		poly->v[2] = dv;
	}
}

double S4r::LatticeGridArb::EdgeLength(size_t i) const{
	size_t iedge = i%3;
	if(0 == iedge){
		return lu;
	}else if(1 == iedge){
		return lv;
	}else{
		return luv;
	}
}
double S4r::LatticeGridArb::FaceArea(size_t i) const{
	return 0.5 * (du[0]*dv[1]-du[1]*dv[0]);
}

double S4r::LatticeGridArb::VertexDualArea(size_t i) const{
	return (du[0]*dv[1]-du[1]*dv[0]);
}
double S4r::LatticeGridArb::EdgeDualLength(size_t i) const{
	size_t iedge = i%3;
	if(0 == iedge){ // u-directed
		return dlu;
	}else if(1 == iedge){
		return dlv;
	}else{
		return dluv;
	}
}

//// Utility functions

PeriodicMesh::PeriodicIndex S4r::LatticeGridArb::ContainingFace(const Vec2 &r, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex ret;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &ret.off[0]);
	double ffy = modf_pos(q[1], &ret.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	if(ffy > ffx){
		ret.idx = 2*(ix+iy*n[0]) + 1;
	}else{
		ret.idx = 2*(ix+iy*n[0]) + 0;
	}
	return ret;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridArb::VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex pi, ci;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &pi.off[0]);
	double ffy = modf_pos(q[1], &pi.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	
	int off[2] = {0,0};
	size_t cx = ix+1, cy = iy+1;
	if(cx >= n[0]){ cx = 0; off[0] = 1; }
	if(cy >= n[1]){ cy = 0; off[1] = 1; }
	
	coeffs.resize(3);
	if(ffy > ffx){
		// Barycentric interpolation of (ffx,ffy) in triangle (0,0), (1,1), (0,1)
		ci.idx = ix+iy*n[0];
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[0] = VertexCoeff(ci, 1.-ffy);
		
		ci.idx = cx+cy*n[0];
		ci.off[0] = pi.off[0]+off[0];
		ci.off[1] = pi.off[1]+off[1];
		coeffs[1] = VertexCoeff(ci, ffx);
		
		ci.idx = ix+cy*n[0];
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1]+off[1];
		coeffs[2] = VertexCoeff(ci, ffy-ffx);
		pi.idx = 2*(ix+iy*n[0]) + 1;
	}else{
		// Barycentric interpolation of (ffx,ffy) in triangle (0,0), (1,0), (1,1)
		ci.idx = ix+iy*n[0];
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[0] = VertexCoeff(ci, 1-ffx);
		
		ci.idx = cx+cy*n[0];
		ci.off[0] = pi.off[0]+off[0];
		ci.off[1] = pi.off[1]+off[1];
		coeffs[1] = VertexCoeff(ci, ffy);
		
		ci.idx = cx+iy*n[0];
		ci.off[0] = pi.off[0]+off[0];
		ci.off[1] = pi.off[1];
		coeffs[2] = VertexCoeff(ci, ffx-ffy);
		pi.idx = 2*(ix+iy*n[0]) + 0;
	}
	return pi;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridArb::EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	PeriodicIndex pi, ci;
	Vec2 q(Lr.partialPivLu().solve(r)); // q is in lattice coordinates
	double ffx = modf_pos(q[0], &pi.off[0]);
	double ffy = modf_pos(q[1], &pi.off[1]);
	size_t ix = (int)(ffx*(double)n[0]);
	size_t iy = (int)(ffy*(double)n[1]);
	
	int off[2] = {0,0};
	size_t cx = ix+1, cy = iy+1;
	if(cx >= n[0]){ cx = 0; off[0] = 1; }
	if(cy >= n[1]){ cy = 0; off[1] = 1; }
	
	coeffs.resize(3);
	if(ffy > ffx){
		// 1-form interpolation of (ffx,ffy) in triangle (0,0), (1,1), (0,1)
		//coeffs[0] = EdgeCoeff(idx(ix,iy,1), Vec2(0,ffy-ffx) + (1.-ffy)*Vec2(-1,1), false, false);
		ci.idx = idx(ix,iy,1);
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[0] = EdgeCoeff(ci, Vec2(ffy-1.,1.-ffx));
		
		ci.idx = idx(ix,iy,2);
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[1] = EdgeCoeff(ci, Vec2(1.-ffy,ffx));
		
		ci.idx = idx(ix,cy,0);
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1]+off[1];
		coeffs[2] = EdgeCoeff(ci, Vec2(ffy, -ffx));
		pi.idx = 2*(ix+iy*n[0]) + 1;
	}else{
		// 1-form interpolation of (ffx,ffy) in triangle (0,0), (1,0), (1,1)
		ci.idx = idx(ix,iy,0);
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[0] = EdgeCoeff(ci, Vec2(1.-ffy, ffx-1.));
		
		ci.idx = idx(ix,iy,2);
		ci.off[0] = pi.off[0];
		ci.off[1] = pi.off[1];
		coeffs[1] = EdgeCoeff(ci, Vec2(ffy, 1.-ffx));
		
		ci.idx = idx(cx,iy,1);
		ci.off[0] = pi.off[0]+off[0];
		ci.off[1] = pi.off[1];
		coeffs[2] = EdgeCoeff(ci, Vec2(-ffy,ffx));
		pi.idx = 2*(ix+iy*n[0]) + 0;
	}
	return pi;
}
PeriodicMesh::PeriodicIndex S4r::LatticeGridArb::FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicMesh::PeriodicIndex *faceguess) const{
	// Simple piecewise constant interpolation
	PeriodicIndex iface = ContainingFace(r, faceguess);
	coeffs.resize(1);
	coeffs[0] = FaceCoeff(iface, 1.);
	return iface;
}

void S4r::LatticeGridArb::BuildMatrices(const doublecomplex phi[2]){
	const size_t N = n[0]*n[1];
	d0.resize(3*N, N);
	std::vector<Triplet> triplst;
	// Each row has two nonzeros
	triplst.reserve(6*N);
	
	for(size_t j = 0; j < n[1]; ++j){
		const size_t jp1 = ((j+1 < n[1]) ? j+1 : 0);
		for(size_t i = 0; i < n[0]; ++i){
			const size_t ip1 = ((i+1 < n[0]) ? i+1 : 0);
			size_t row;
			// u directed edge
			row = 3*(i+j*n[0])+0;
			if(i+1 < n[0]){
				triplst.push_back(Triplet(
					row,
					i+j*n[0],
					-1.0
				));
				triplst.push_back(Triplet(
					row,
					ip1+j*n[0],
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					row,
					i+j*n[0],
					-std::conj(phi[0])
				));
				triplst.push_back(Triplet(
					row,
					ip1+j*n[0],
					phi[0]
				));
			}
			// v directed edge
			row = 3*(i+j*n[0])+1;
			if(j+1 < n[1]){
				triplst.push_back(Triplet(
					row,
					i+j*n[0],
					-1.0
				));
				triplst.push_back(Triplet(
					row,
					i+jp1*n[0],
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					row,
					i+j*n[0],
					-std::conj(phi[1])
				));
				triplst.push_back(Triplet(
					row,
					i+jp1*n[0],
					phi[1]
				));
			}
			// u+v directed edge
			row = 3*(i+j*n[0])+2;
			doublecomplex coef(1.);
			if(i+1 >= n[0]){
				coef *= phi[0];
			}
			if(j+1 >= n[1]){
				coef *= phi[1];
			}
			triplst.push_back(Triplet(
				row,
				i+j*n[0],
				-std::conj(coef)
			));
			triplst.push_back(Triplet(
				row,
				ip1+jp1*n[0],
				coef
			));
		}
	}
	d0.setFromTriplets(triplst.begin(), triplst.end());
	
	d1.resize(2*N, 3*N);
	// Each row has three nonzeros
	triplst.clear();
	triplst.reserve(6*N);
	for(size_t j = 0; j < n[1]; ++j){
		const size_t jp1 = ((j+1 < n[1]) ? j+1 : 0);
		for(size_t i = 0; i < n[0]; ++i){
			const size_t ip1 = ((i+1 < n[0]) ? i+1 : 0);
			size_t row;
			// face 1
			row = 2*(i+j*n[0])+0;
			if(i+1 < n[0]){
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+0,
					1.0
				));
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+2,
					-1.0
				));
				triplst.push_back(Triplet(
					row,
					3*(ip1+j*n[0])+1,
					1.0
				));
			}else{
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+0,
					std::conj(phi[0])
				));
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+2,
					-std::conj(phi[0])
				));
				triplst.push_back(Triplet(
					row,
					3*(ip1+j*n[0])+1,
					phi[0]
				));
			}
			// face 2
			row = 2*(i+j*n[0])+1;
			if(j+1 < n[1]){
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+2,
					1.0
				));
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+1,
					-1.0
				));
				triplst.push_back(Triplet(
					row,
					3*(i+jp1*n[0])+0,
					-1.0
				));
			}else{
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+2,
					std::conj(phi[1])
				));
				triplst.push_back(Triplet(
					row,
					3*(i+j*n[0])+1,
					-std::conj(phi[1])
				));
				triplst.push_back(Triplet(
					row,
					3*(i+jp1*n[0])+0,
					-phi[1]
				));
			}
		}
	}
	d1.setFromTriplets(triplst.begin(), triplst.end());
}
