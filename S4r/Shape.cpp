#include <float.h>
#include <cstdio>
#include <S4r/Shape.hpp>
#include <Eigen/Geometry>
extern "C" {
#include <S4r/intersection.h>
}
typedef Eigen::Rotation2D<double> Rot2;
using namespace S4r;

// returns true if inside or on boundary
bool TriangleContainsPoint(
	const Vec2 &org, // triangle vertices are {org,org+u,org+v}, in CCW orientation
	const Vec2 &u,
	const Vec2 &v,
	const Vec2 &p // query point
){
	if(Orient2D(org,u,p) >= 0){
		if(Orient2D(org,v,p) <= 0){
			// treat org as origin, we have points u and v, just need p-org
			Vec2 x = p-org;
			if(Orient2D(u,v,x) >= 0){
				return true;
			}
		}
	}
	return false;
}

S4r::ConvexPolygon::ConvexPolygon(size_t n, double *p):
	v(n), offset(0,0)
{
	for(size_t i = 0; i < n; ++i){
		v[i] = Vec2(p[2*i+0], p[2*i+1]);
	}
}
double S4r::ConvexPolygon::Area() const{
	size_t p, q;
	const size_t n = v.size();
	double area = 0.;
	for(p=n-1, q=0; q < n; p = q++){
		area += v[p][0]*v[q][1] - v[q][0]*v[p][1];
	}
	return 0.5*area;
}

Vec2 S4r::ConvexPolygon::ApproxCenter() const{
	Vec2 c(v[0]);
	for(size_t i = 1; i < v.size(); ++i){
		c += v[i];
	}
	return offset + (c / (double)v.size());
}


//// Circle

S4r::ShapeCircle::ShapeCircle(const Vec2 &center, const double &radius, int tag):
	Shape(center, 0, tag),
	radius(radius)
{
}

Vec2 S4r::ShapeCircle::Normal(const Vec2 &r) const{
	Vec2 d(r - this->center);
	double l = d.norm();
	if(0 == l){ return Vec2(1,0); }
	return d/l;
}

bool S4r::ShapeCircle::Inside(const Vec2 &r) const{
	return (r - this->center).squaredNorm() <= radius*radius;
}

double S4r::ShapeCircle::Area() const{
	return M_PI*radius*radius;
}

double S4r::ShapeCircle::OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const{
	double org[2] = {p0[0] - this->center[0], p0[1] - this->center[1]};
	double u[2] = {p1p0[0], p1p0[1]};
	double v[2] = {p2p0[0], p2p0[1]};
	return intersection_area_circle_triangle(radius, org, u, v);
}

void S4r::ShapeCircle::OutputPostscript(FILE *fp) const{
	fprintf(fp, "newpath 0 0 %f 0 360 arc closepath stroke\n", radius);
}

//// Ellipse

S4r::ShapeEllipse::ShapeEllipse(const Vec2 &center, const double &angle, const Vec2 &halfwidth, int tag):
	Shape(center, angle, tag), halfwidth(halfwidth), Q(Rot2(angle).matrix())
{
}

// Taking the gradient of |B(r-c)|^2 = 1, where B = inv(diag(halfwidth)) * Q,
// we get 2B(r-c)
Vec2 S4r::ShapeEllipse::Normal(const Vec2 &r) const{
	Mat2 B(halfwidth.asDiagonal() * Q);
	Vec2 d(B.transpose() * B * (r - this->center));
	double l = d.norm();
	if(0 == l){ return Vec2(1,0); }
	return d/l;
}

bool S4r::ShapeEllipse::Inside(const Vec2 &r) const{
	return (halfwidth.cwiseInverse().asDiagonal() * Q.transpose() * (r - this->center)).squaredNorm() <= 1.;
}

double S4r::ShapeEllipse::Area() const{
	return M_PI*halfwidth[0]*halfwidth[1];
}

double S4r::ShapeEllipse::OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const{
	const Vec2 o(halfwidth.cwiseInverse().asDiagonal() * Q.transpose() * (p0 - this->center));
	const double org[2] = {o[0], o[1]};
	const Vec2 a(halfwidth.cwiseInverse().asDiagonal() * Q.transpose() * p1p0);
	const double u[2] = {a[0], a[1]};
	const Vec2 b(halfwidth.cwiseInverse().asDiagonal() * Q.transpose() * p2p0);
	const double v[2] = {b[0], b[1]};
	//std::cerr << intersection_area_circle_triangle(1., org, u, v) << "\n";
	//std::cerr << "org: " << org[0] << ' ' << org[1] << "\n";
	//std::cerr << "  u: " << u[0] << ' ' << u[1] << "\n";
	//std::cerr << "  v: " << v[0] << ' ' << v[1] << "\n";
	return halfwidth[0] * halfwidth[1] * intersection_area_circle_triangle(1., org, u, v);
}

void S4r::ShapeEllipse::OutputPostscript(FILE *fp) const{
	fprintf(fp, "gsave 1 %f scale ", halfwidth[1]/halfwidth[0]);
	fprintf(fp, "newpath 0 0 %f 0 360 arc closepath stroke ", halfwidth[0]);
	fprintf(fp, "grestore\n");
}


//// Rectangle

S4r::ShapeRectangle::ShapeRectangle(const Vec2 &center, const double &angle, const Vec2 &halfwidth, int tag):
	Shape(center, angle, tag), halfwidth(halfwidth), Q(Rot2(angle).matrix())
{
}

Vec2 S4r::ShapeRectangle::Normal(const Vec2 &r) const{
	Vec2 p(Q.transpose() * (r - this->center));
	if(fabs(p[0]) - halfwidth[0] > fabs(p[1]) - halfwidth[1]){
		if(p[0] > 0){
			return Q.col(0);
		}else{
			return -Q.col(0);
		}
	}else{
		if(p[1] > 0){
			return Q.col(1);
		}else{
			return -Q.col(1);
		}
	}
}

bool S4r::ShapeRectangle::Inside(const Vec2 &r) const{
	Vec2 p(Q.transpose() * (r - this->center));
	return fabs(p[0]) <= halfwidth[0] && fabs(p[1]) <= halfwidth[1];
}

double S4r::ShapeRectangle::Area() const{
	return 4.*halfwidth[0]*halfwidth[1];
}

double S4r::ShapeRectangle::OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const{
	const double ca = cos(angle);
	const double sa = sin(angle);
	double P[6], Q[8], Pi[14], u[2], v[2];
	int ret, nPi = 7;

	P[2*0+0] = p0[0] - this->center[0];
	P[2*0+1] = p0[1] - this->center[1];
	P[2*1+0] = P[2*0+0] + p1p0[0];
	P[2*1+1] = P[2*0+1] + p1p0[1];
	P[2*2+0] = P[2*1+0] + p2p0[0];
	P[2*2+1] = P[2*1+1] + p2p0[0];
	u[0] = halfwidth[0] * ca;
	u[1] = halfwidth[0] * sa;
	v[0] = halfwidth[1] *-sa;
	v[1] = halfwidth[1] * ca;
	Q[2*0+0] = -u[0]-v[0];
	Q[2*0+1] = -u[1]-v[1];
	Q[2*1+0] = u[0]-v[0];
	Q[2*1+1] = u[1]-v[1];
	Q[2*2+0] = u[0]+v[0];
	Q[2*2+1] = u[1]+v[1];
	Q[2*3+0] = v[0]-u[0];
	Q[2*3+1] = v[1]-u[1];
	/*{int i;for(i = 0;i<4;++i){fprintf(stderr, " {%f,%f},\n", P[2*i+0], P[2*i+1]);}}*/
	/*{int i;for(i = 0;i<4;++i){fprintf(stderr, " {%f,%f},\n", Q[2*i+0], Q[2*i+1]);}}*/
	
	ret = convex_polygon_intersection(3,P,4,Q,&nPi,Pi);
	if(1 == ret){ /* pixel completely in rectangle */
		return p1p0[0]*p2p0[1]-p1p0[1]*p2p0[0];
	}
	/*{int i;for(i = 0;i<nPi;++i){fprintf(stderr, " %f %f\n", Pi[2*i+0], Pi[2*i+1]);}}*/
	return polygon_area(nPi,Pi);
}

void S4r::ShapeRectangle::OutputPostscript(FILE *fp) const{
	fprintf(fp, "gsave 1 %f scale ", halfwidth[1]/halfwidth[0]);
	fprintf(fp, "%g %g moveto %g %g lineto %g %g lineto %g %g lineto closepath stroke ",
		-halfwidth[0], -halfwidth[1],
		 halfwidth[0], -halfwidth[1],
		 halfwidth[0],  halfwidth[1],
		-halfwidth[0],  halfwidth[1]
	);
	fprintf(fp, "grestore\n");
}


//// Polygon

S4r::ShapePolygon::ShapePolygon(const Vec2 &center, const double &angle, const std::vector<Vec2> &vertex, int tag):
	Shape(center, angle, tag), vertex(vertex)
{
	const size_t n = vertex.size();
	t.resize(3*(n-2));
	// t is used to store a working copy of the currently clipped polygon vertex index, size 3 <= m <= n.
	// t must also store the resulting triangles, size 3*(n-m).
	// The total size of these two is 3*n-2*m, which must be <= 3*(n-2). This is true for m >= 3.
	// Therefore, we keep the working copy at the very end of t, and add triangles to the beginning.
	std::vector<size_t>::iterator V; // pointer to start of working copy
	size_t nv; // size of V
	size_t count;
	size_t tc = 0; // number of triangles currently in t
	
	// Make a copy of all the vertices
	V = t.begin() + 2*n-6;
	nv = n;
	for(size_t i = 0; i < n; ++i){ V[i] = i; }
	
	count = 2*nv;
	for(size_t i = nv-1; nv > 2; ){
		/*
		fprintf(stderr, "V:");
		for(j = 0; j < nv; ++j){
			fprintf(stderr, " %d", V[j]);
		}fprintf(stderr, "\n");
		fprintf(stderr, "count = %d\n", count);
		*/
		size_t u, w;
		if(0 >= (count--)){
			// bad polygon
		}
		// get 3 consecutive vertices
		u = i; //if(nv <= u){ u = 0; } // prev
		i = u+1; if(nv <= i){ i = 0; } // mid
		w = i+1; if(nv <= w){ w = 0; } // next

		// Can clip the ear?
		int can_clip = 1;
		{
			Vec2 tri_a = vertex[V[i]] - vertex[V[u]];
			Vec2 tri_b = vertex[V[w]] - vertex[V[u]];
			if(Orient2D(vertex[V[u]], vertex[V[i]], vertex[V[w]]) < 0){
				can_clip = 0;
			}else{
				// if the u-i-w triangle contains any other vertex, can't clip.
				for(size_t p = 0; p < nv; ++p){
					if((p == u) || (p == i) || (p == w)){ continue; }
					if(TriangleContainsPoint(vertex[V[u]],tri_a,tri_b, vertex[V[p]])){ can_clip = 0; break; }
				}
			}
		}

		// Clip off the ear
		if(can_clip){
			size_t tri[3] = {V[u], V[i], V[w]};
			// erase vertex i
			while(i > 0){
				V[i] = V[i-1];
				--i;
			}
			++V; --nv; count = 2*nv;
			// Add the new triangle
			t[3*tc+0] = tri[0];
			t[3*tc+1] = tri[1];
			t[3*tc+2] = tri[2];
//std::cerr << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
			++tc;
		}
	}
}

Vec2 S4r::ShapePolygon::Normal(const Vec2 &r) const{
	double mindist = DBL_MAX;
	Vec2 n;
	Mat2 Q(Rot2(angle).matrix());
	Vec2 x(Q.transpose() * (r - this->center));
	for(size_t j = 0, i = vertex.size()-1; j < vertex.size(); i = j++){
		/* compute distance from r to segment */
		double v[2], pr[2];
		v[0] = vertex[j][0] - vertex[i][0];
		v[1] = vertex[j][1] - vertex[i][1];
		
		pr[0] = x[0] - vertex[i][0];
		pr[1] = x[1] - vertex[i][1];
		
		double dist;
		double voff[2];
		{ // compute distance from point to line segment
			const double v2 = v[0]*v[0] + v[1]*v[1];
			double prj = (pr[0]*v[0] + pr[1]*v[1])/v2;
			if(prj > 1){
				voff[0] = x[0] - vertex[j][0];
				voff[1] = x[1] - vertex[j][1];
				dist = hypot(voff[0], voff[1]);
			}else if(prj < 0){
				voff[0] = pr[0];
				voff[1] = pr[1];
				dist = hypot(pr[0], pr[1]);
			}else{
				voff[0] = pr[0] - prj*v[0];
				voff[1] = pr[1] - prj*v[1];
				dist = hypot(voff[0], voff[1]);
			}
		}
		if(dist < mindist){
			mindist = dist;
			n[0] = voff[0];
			n[1] = voff[1];
			if(n[0] * v[1] - n[1] * v[0] < 0){
				// Align normals to be outward-ish
				n = -n;
			}
		}
	}
	return (Q*n).normalized();
}

bool S4r::ShapePolygon::Inside(const Vec2 &r) const{
	Mat2 Q(Rot2(angle).matrix());
	Vec2 x(Q.transpose() * (r - this->center));
	
	// From http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
	bool c = false;
	for(size_t i = 0, j = vertex.size()-1; i < vertex.size(); j = i++){
		double vix = vertex[i][0];
		double viy = vertex[i][1];
		double vjx = vertex[j][0];
		double vjy = vertex[j][1];
		if ( ((viy>x[1]) != (vjy>x[1]))
		&& (x[0] < (vjx-vix) * (x[1]-viy) / (vjy-viy) + vix) ){ c = !c; }
	}
	return c;
}

double S4r::ShapePolygon::Area() const{
	size_t p, q;
	const size_t n = vertex.size();
	double area = 0.;
	for(p=n-1, q=0; q < n; p = q++){
		area += vertex[p][0]*vertex[q][1] - vertex[q][0]*vertex[p][1];
	}
	return 0.5*area;
}

Vec2 S4r::ShapePolygon::GetSomeInteriorPoint() const{
	return this->center + (1./3.) * (vertex[t[0]] + vertex[t[1]] + vertex[t[2]]);
}

double S4r::ShapePolygon::OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const{
	const double ca = cos(angle);
	const double sa = sin(angle);
	double a = 0;
	double P[6], Q[6], Pi[12];
	
	double p0p[2] = {p0[0] - this->center[0], p0[1] - this->center[1]};
	P[2*0+0] = p0p[0]*ca + p0p[1]*sa;
	P[2*0+1] = p0p[0]*-sa + p0p[1]*ca;
	P[2*1+0] = (p1p0[0]+p0p[0])*ca + (p1p0[1]+p0p[1])*sa;
	P[2*1+1] = (p1p0[0]+p0p[0])*-sa+ (p1p0[1]+p0p[1])*ca;
	P[2*2+0] = (p2p0[0]+p0p[0])*ca + (p2p0[1]+p0p[1])*sa;
	P[2*2+1] = (p2p0[0]+p0p[0])*-sa+ (p2p0[1]+p0p[1])*ca;
	
	/*
	for(i = 0; i < s->vtab.polygon.n_vertices-2; ++i){
		fprintf(stderr, " (%d %d %d)", t[3*i+0], t[3*i+1], t[3*i+2]);
	}fprintf(stderr, "\n");
	*/
	for(size_t i = 0; i < vertex.size()-2; ++i){
		for(unsigned j = 0; j < 3; ++j){
			Q[2*j+0] = vertex[t[3*i+j]][0];
			Q[2*j+1] = vertex[t[3*i+j]][1];
		}
		int nPi = 6;
		convex_polygon_intersection(3,P,3,Q,&nPi,Pi);
		if(nPi >= 3){
			a += polygon_area(nPi,Pi);
		}
		/*
		std::cerr << "a = " << a << "\n";
		std::cerr << "  ret = " << ret << ", nPi = " << nPi << "\n";
		std::cerr << "  P =";
		for(size_t k = 0; k < 6; ++k){
			std::cerr << ' ' << P[k];
		} std::cerr << "\n";
		std::cerr << "  Q =";
		for(size_t k = 0; k < 6; ++k){
			std::cerr << ' ' << Q[k];
		} std::cerr << "\n";
		std::cerr << "  |P| = " << polygon_area(3, P) << "\n";
		std::cerr << "  |Q| = " << polygon_area(3, Q) << "\n";
		*/
	}
	//std::cerr << "\n";
	return a;
}

void S4r::ShapePolygon::OutputPostscript(FILE *fp) const{
	if(vertex.size() >= 3){
		fprintf(fp, "newpath %f %f moveto ", vertex[0][0], vertex[0][1]);
		for(size_t j = 1; j < vertex.size(); ++j){
			fprintf(fp, "%f %f lineto ", vertex[j][0], vertex[j][1]);
		}
		fprintf(fp, "closepath stroke\n");
	}
}

//// Pattern

void S4r::Pattern::Finalize(){
	if(finalized){ return; }
	
	// We will assume for now that the non-self-intersection criterion is met.
	parent.resize(shape.size());
	std::vector<double> area(shape.size());
	for(size_t i = 0; i < shape.size(); ++i){
		parent[i] = -1;
		area[i] = shape[i]->Area();
	}
	
	// Sort by area
	for(size_t k = 1; k < shape.size(); ++k){
		for(size_t j = 0; j < k; ++j){
			if(area[j] < area[k]){
				std::swap(shape[j], shape[k]);
				std::swap(area[j], area[k]);
			}
		}
	}
	
	// The following is O(n^2), we might be able to do better with sorting
	// by area but n is usually small.
	for(size_t i = 1; i < shape.size(); ++i){
		const Vec2 p(shape[i]->GetSomeInteriorPoint());

		for(size_t j = i; j--; ){
			if(shape[j]->Inside(p)){
				parent[i] = j;
				break;
			}
		}
	}
}

size_t S4r::Pattern::Overlap(
	const ConvexPolygon &poly,
	std::vector<double> &value
) const{
	value.resize(shape.size()+1);
	value[0] = 1.;
	for(size_t k = 1; k <= shape.size(); ++k){
		value[k] = 0;
	}
	const double area = poly.Area();
	const double inv_area = 1. / area;
	const Vec2 org(poly.offset + poly.v[0]);
	for(size_t j = 2; j < poly.v.size(); ++j){
		const Vec2 u(poly.v[j-1]-poly.v[0]);
		const Vec2 v(poly.v[j  ]-poly.v[0]);
		for(size_t k = 0; k < shape.size(); ++k){
			double a = shape[k]->OverlapTriangle(org, u, v) * inv_area;
			if(a > 0){
				value[k+1] += a;
				value[parent[k]+1] -= a;
			}
		}
	}
	size_t ret = 0;
	const double tol = 4.*DBL_EPSILON * area;
	for(size_t k = 0; k <= shape.size(); ++k){
		if(value[k] < tol){ value[k] = 0; }
		else if(value[k] > 1){ value[k] = 1; }
		if(value[k] > 0){ ret++; }
	}
	return ret;
}

S4r::Pattern::~Pattern(){
	for(size_t i = 0; i < shape.size(); ++i){
		delete shape[i];
	}
}

void S4r::Pattern::OutputPostscript(FILE *fp) const{
	for(size_t i = 0; i < shape.size(); ++i){
		const Shape &s = *shape[i];
		fprintf(fp, "gsave %f %f translate %f rotate\n", s.center[0], s.center[1], s.angle * 180./M_PI);
		s.OutputPostscript(fp);
		fprintf(fp, "grestore\n");
	}
}

