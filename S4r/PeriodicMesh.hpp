#ifndef S4R_PERIODIC_MESH_HPP_INCLUDED
#define S4R_PERIODIC_MESH_HPP_INCLUDED

#include <S4r/Types.hpp>
#include <S4r/Shape.hpp>

#include <S4r/periodic_off2.h>

namespace S4r{

class PeriodicMesh{ // abstract base
protected:
	Mat2 Lr;
public:
	struct PeriodicIndex{
		size_t idx;
		int off[2];
	};
	template <typename T>
	struct PeriodicCoeff{
		PeriodicIndex pi;
		T c;
		PeriodicCoeff(){}
		PeriodicCoeff(const PeriodicIndex &i, const T &co):pi(i),c(co){}
	};
	typedef PeriodicCoeff<double> VertexCoeff;
	typedef PeriodicCoeff<double> FaceCoeff;
	typedef PeriodicCoeff<Vec2> EdgeCoeff;
/*
	struct VertexCoeff{
		size_t i;
		double c;
		bool wx, wy;
		VertexCoeff(){}
		VertexCoeff(size_t i, double coeff):pi(i),c(coeff){}
	};
	struct EdgeCoeff{
		size_t i;
		Vec2 c;
		bool wx, wy;
		EdgeCoeff(){}
		EdgeCoeff(size_t idx, Vec2 coeff, bool wrapx, bool wrapy):i(idx),c(coeff),wx(wrapx),wy(wrapy){}
	};
	struct FaceCoeff{
		size_t i;
		double c;
		bool wx, wy;
		FaceCoeff(){}
		FaceCoeff(size_t idx, double coeff, bool wrapx, bool wrapy):i(idx),c(coeff),wx(wrapx),wy(wrapy){}
	};
	*/
	typedef std::vector<VertexCoeff> VertexCoeffs;
	typedef std::vector<EdgeCoeff> EdgeCoeffs;
	typedef std::vector<FaceCoeff> FaceCoeffs;
	
	// Members
	SpMat d0, d1; // Discrete exterior derivative matrices
	
	// Constructor/Destructor
	PeriodicMesh(const Vec2 &u, const Vec2 &v){
		Lr(0,0) = u[0]; Lr(0,1) = v[0];
		Lr(1,0) = u[1]; Lr(1,1) = v[1];
	}
	virtual ~PeriodicMesh(){}
	
	// Basic geometric properties
	
	virtual bool IsOrthogonal() const = 0;
	
	virtual size_t NumVertices() const = 0;
	virtual size_t NumEdges() const = 0;
	virtual size_t NumFaces() const = 0;
	virtual Vec2 Vertex(size_t i) const = 0;
	virtual Vec2 Edge(size_t i) const = 0;
	
	// Control polygons
	virtual void VertexControlPolygon(size_t i, ConvexPolygon *poly) const = 0;
	virtual void EdgeControlPolygon(size_t i, ConvexPolygon *self,
	                                size_t *e1, ConvexPolygon *wedge1,
	                                size_t *e2, ConvexPolygon *wedge2) const = 0;
	virtual void FaceControlPolygon(size_t i, ConvexPolygon *poly) const = 0;
	
	virtual double EdgeLength(size_t i) const = 0;
	virtual double FaceArea(size_t i) const = 0;
	
	virtual double VertexDualArea(size_t i) const = 0;
	virtual double EdgeDualLength(size_t i) const = 0;

	virtual PeriodicIndex ContainingFace(const Vec2 &r, PeriodicIndex *faceguess = NULL) const = 0;
	virtual PeriodicIndex VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const = 0;
	virtual PeriodicIndex EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const = 0;
	virtual PeriodicIndex FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const = 0;
	
	//// Utility functions
	
	const Mat2 &GetLattice() const{ return Lr; }
	
	// Brings an arbitrary 2D point into fundamental parallelogram.
	void Remap(Vec2 &r) const;
	// Computes d0 and d1
	virtual void BuildMatrices(const doublecomplex phi[2]) = 0;
};

class LatticeGridRect : public PeriodicMesh{
	size_t n[2];
	Vec2 du, dv, hu, hv;
	
	size_t idx(size_t ix, size_t iy, size_t vert) const;
public:
	LatticeGridRect(const Vec2 &u, const Vec2 &v, const size_t ncells[2]);
	~LatticeGridRect();
	
	// Basic geometric properties
	
	bool IsOrthogonal() const{ return true; }
	
	size_t NumVertices() const;
	size_t NumEdges() const;
	size_t NumFaces() const;
	Vec2 Vertex(size_t i) const;
	Vec2 Edge(size_t i) const;
	
	// Control polygons
	void VertexControlPolygon(size_t i, ConvexPolygon *poly) const;
	void EdgeControlPolygon(size_t i, ConvexPolygon *self,
	                        size_t *e1, ConvexPolygon *wedge1,
	                        size_t *e2, ConvexPolygon *wedge2) const;
	void FaceControlPolygon(size_t i, ConvexPolygon *poly) const;
	
	double EdgeLength(size_t i) const;
	double FaceArea(size_t i) const;
	
	double VertexDualArea(size_t i) const;
	double EdgeDualLength(size_t i) const;

	PeriodicIndex ContainingFace(const Vec2 &r, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	
	void BuildMatrices(const doublecomplex phi[2]);
};

// We assume the lattice vectors are Delaunay reduced, so that they
// form an obtuse angle:
//
//    ^-----+        .
//     \ 2 / \       .
//      \ / 1 \      .
//       +----->     .
//
class LatticeGridArb : public PeriodicMesh{
	size_t n[2];
	Vec2 c1, c2; // circumcenters
	Vec2 du, dv;
	double lu, lv, luv; // edge lengths
	double dlu, dlv, dluv; // dual edge lengths
	
	size_t idx(size_t ix, size_t iy, size_t vert) const;
public:
	LatticeGridArb(const Vec2 &u, const Vec2 &v, const size_t ncells[2]);
	~LatticeGridArb();
	
	// Basic geometric properties
	
	bool IsOrthogonal() const{ return false; }
	
	size_t NumVertices() const;
	size_t NumEdges() const;
	size_t NumFaces() const;
	Vec2 Vertex(size_t i) const;
	Vec2 Edge(size_t i) const;
	
	// Control polygons
	void VertexControlPolygon(size_t i, ConvexPolygon *poly) const;
	void EdgeControlPolygon(size_t i, ConvexPolygon *self,
	                        size_t *e1, ConvexPolygon *wedge1,
	                        size_t *e2, ConvexPolygon *wedge2) const;
	void FaceControlPolygon(size_t i, ConvexPolygon *poly) const;
	
	double EdgeLength(size_t i) const;
	double FaceArea(size_t i) const;
	
	double VertexDualArea(size_t i) const;
	double EdgeDualLength(size_t i) const;

	PeriodicIndex ContainingFace(const Vec2 &r, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;

	void BuildMatrices(const doublecomplex phi[2]);
};

class POFF2Mesh : public PeriodicMesh{
	::POFF2Mesh M;
public:
	POFF2Mesh(const Vec2 &u, const Vec2 &v, const char *filename);
	~POFF2Mesh();
	
	// Basic geometric properties
	
	bool IsOrthogonal() const{ return false; }
	
	size_t NumVertices() const;
	size_t NumEdges() const;
	size_t NumFaces() const;
	Vec2 Vertex(size_t i) const;
	Vec2 Edge(size_t i) const;
	
	// Control polygons
	void VertexControlPolygon(size_t i, ConvexPolygon *poly) const;
	void EdgeControlPolygon(size_t i, ConvexPolygon *self,
	                        size_t *e1, ConvexPolygon *wedge1,
	                        size_t *e2, ConvexPolygon *wedge2) const;
	void FaceControlPolygon(size_t i, ConvexPolygon *poly) const;
	
	double EdgeLength(size_t i) const;
	double FaceArea(size_t i) const;
	
	double VertexDualArea(size_t i) const;
	double EdgeDualLength(size_t i) const;

	PeriodicIndex ContainingFace(const Vec2 &r, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex VertexInterpolation(const Vec2 &r, VertexCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex EdgeInterpolation(const Vec2 &r, EdgeCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	PeriodicIndex FaceInterpolation(const Vec2 &r, FaceCoeffs &coeffs, PeriodicIndex *faceguess = NULL) const;
	
	void BuildMatrices(const doublecomplex phi[2]);
};

} // namespace S4r

#endif // S4R_PERIODIC_MESH_HPP_INCLUDED
