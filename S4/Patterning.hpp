#ifndef PATTERNING_HPP_INCLUDED
#define PATTERNING_HPP_INCLUDED

#include "S4.h"
#include <complex>
#include <vector>
//L:MaterialMap{ image = 'foo.png', map = { ["000000"] = mSi } }

class Patterning{
	// Represents a generic layer patterning
public:
	typedef S4_real real_type;
	typedef std::complex<real_type> complex_type;
public:
	Patterning();
	virtual ~Patterning(){}
	
	// For types of patterning in which discrete material regions are used,
	// this function allows setting the mapping from tags to complex epsilon.
	// (Each region containing a material is tagged with the material id.)
	void SetTagToValueMap(const complex_type *m, int incm);
	
	const complex_type &TagToValue(int tag) const;
	
	// Computes the Fourier series coefficient of the patterning at a particular
	// spatial frequency vector given by:
	//    [ Lk[0]  Lk[2] ] * [ ik[0] ]
	//    [ Lk[1]  Lk[3] ]   [ ik[1] ]
	//   = ( Lk[0]*ik[0]+Lk[2]*ik[1], Lk[1]*ik[0]+Lk[3]*ik[1] )
	virtual void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc) const = 0;
protected:
	const complex_type *t2v; // tag to value map
	int inct2v;
};

class PatterningByBreakpoints : public Patterning{
	// For a 1D patterning, this defines a patterning by regions as follows:
	// reg[i] is defined as containing material reg[i].tag from reg[i].x to reg[i+1].x,
	// where the x coordinate periodically wraps around by the lattice constant.
protected:
	struct Region{
		real_type xi; // in range of [0,1) for one period
		int tag;
		Region(const real_type &xi, int t):xi(xi),tag(t){}
	};
	std::vector<Region> reg;
	real_type period, iperiod;
public:
	PatterningByBreakpoints(const real_type &L);
	void SetRegion(const real_type &center, const real_type &halfwidth, int tag);
};

class PatterningByShapes : public Patterning{
public:
	// Abstract base classes for shapes
	class Shape{
		real_type center[2];
		real_type angle_frac;
	protected:
		int parent;
		int tag;
		Shape(const real_type *center = NULL, const real_type &angle_frac = 0);
		friend class PatterningByShapes;
	public:
		virtual ~Shape(){}
		bool Contains(const real_type *p) const;
		void Center(real_type *center) const;
		void Normal(const real_type *p, real_type *n) const;
		void FourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type IntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
		
		// These functions assume center = 0, angle = 0
		virtual real_type Area() const = 0;
		virtual bool BaseContains(const real_type *p) const = 0;
		virtual void BaseCenter(real_type *center) const = 0;
		virtual void BaseNormal(const real_type *p, real_type *n) const = 0;
		virtual void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const = 0;
		virtual real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const = 0;
	};
	class ShapeHalfwidths : public Shape{
	protected:
		real_type halfwidth[2];
		ShapeHalfwidths(const real_type &hx, const real_type &hy, const real_type *center = NULL, const real_type &angle_frac = 0);
	public:
		virtual ~ShapeHalfwidths(){}
		void BaseCenter(real_type *center) const;
	};
	class ShapeVertices : public Shape{
	protected:
		int n;
		real_type *v;
		int incv;
		ShapeVertices(int n, const real_type *v, int incv, const real_type *center = NULL, const real_type &angle_frac = 0);
	public:
		virtual ~ShapeVertices();
	};
	
	class Circle : public ShapeHalfwidths{
	public:
		Circle(const real_type &radius, const real_type *center = NULL);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
	
	class Ellipse : public ShapeHalfwidths{
	public:
		Ellipse(const real_type *halfwidths, const real_type *center = NULL, const real_type &ang = 0);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
	
	class Rectangle : public ShapeHalfwidths{
	public:
		Rectangle(const real_type *halfwidths, const real_type *center = NULL, const real_type &ang = 0);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
	
	class Polygon : public ShapeVertices{
	public:
		Polygon(int nv, const real_type *v, const real_type *center = NULL, const real_type &ang = 0);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseCenter(real_type *center) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
	
	class Arcpoly : public ShapeVertices{
	public:
		Arcpoly(int narcs, const real_type *v, const real_type *center = NULL, const real_type &ang = 0);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseCenter(real_type *center) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
	
	class QuadraticBezier : public ShapeVertices{
	public:
		QuadraticBezier(int nbez, const real_type *v, const real_type *center = NULL, const real_type &ang = 0);
		
		real_type Area() const;
		bool BaseContains(const real_type *p) const;
		void BaseCenter(real_type *center) const;
		void BaseNormal(const real_type *p, real_type *n) const;
		void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
		real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	};
public:
	PatterningByShapes();
	~PatterningByShapes();
	void AddShape(Shape *s, int tag = 0);
	virtual void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc) const;
protected:
	std::vector<Shape*> shape;
};

class PatterningByMap : public Patterning{
protected:
	std::vector<int> tag;
	int nu, nv;
};

class PatterningByArray : public Patterning{
protected:
	std::vector<complex_type> val;
	int nu, nv;
};

#endif // PATTERNING_HPP_INCLUDED
