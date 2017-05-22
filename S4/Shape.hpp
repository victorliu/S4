#ifndef SHAPE_HPP_INCLUDED
#define SHAPE_HPP_INCLUDED

#include "S4.h"

class PatterningByShapes;

// Abstract base classes for shapes
class Shape{
public:
	typedef S4_real real_type;
	typedef std::complex<real_type> complex_type;
protected:
	real_type center[2];
	real_type angle_frac;
protected:
	int parent;
	int tag;
	Shape(const real_type *center = NULL, const real_type &angle_frac = 0);
	friend class PatterningByShapes;
public:
	virtual ~Shape(){}
	virtual Shape* Clone() const = 0;
	bool Contains(const real_type *p) const;
	void Center(real_type *center) const;
	void Normal(const real_type *p, real_type *n) const;
	void FourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type IntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
	
	// These functions assume center = 0, angle = 0
	virtual real_type Area() const = 0;
	virtual bool BaseInversionSymmetric() const = 0;
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
	bool BaseInversionSymmetric() const{ return true; }
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
	Circle *Clone() const;
	
	real_type Area() const;
	bool BaseContains(const real_type *p) const;
	void BaseNormal(const real_type *p, real_type *n) const;
	void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
};

class Ellipse : public ShapeHalfwidths{
public:
	Ellipse(const real_type *halfwidths, const real_type *center = NULL, const real_type &ang = 0);
	Ellipse *Clone() const;
	
	real_type Area() const;
	bool BaseContains(const real_type *p) const;
	void BaseNormal(const real_type *p, real_type *n) const;
	void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
};

class Rectangle : public ShapeHalfwidths{
public:
	Rectangle(const real_type *halfwidths, const real_type *center = NULL, const real_type &ang = 0);
	Rectangle *Clone() const;
	
	real_type Area() const;
	bool BaseContains(const real_type *p) const;
	void BaseNormal(const real_type *p, real_type *n) const;
	void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
};

class Polygon : public ShapeVertices{
public:
	Polygon(int nv, const real_type *v, const real_type *center = NULL, const real_type &ang = 0);
	Polygon *Clone() const;
	
	real_type Area() const;
	bool BaseInversionSymmetric() const;
	bool BaseContains(const real_type *p) const;
	void BaseCenter(real_type *center) const;
	void BaseNormal(const real_type *p, real_type *n) const;
	void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
};

class Arcpoly : public ShapeVertices{
public:
	Arcpoly(int narcs, const real_type *v, const real_type *center = NULL, const real_type &ang = 0);
	Arcpoly *Clone() const;
	
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
	QuadraticBezier *Clone() const;
	
	real_type Area() const;
	bool BaseContains(const real_type *p) const;
	void BaseCenter(real_type *center) const;
	void BaseNormal(const real_type *p, real_type *n) const;
	void BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const;
	real_type BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv = NULL) const;
};

#endif // SHAPE_HPP_INCLUDED
