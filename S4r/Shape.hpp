#ifndef S4R_SHAPE_HPP_INCLUDED
#define S4R_SHAPE_HPP_INCLUDED

#include <S4r/Types.hpp>
#include <vector>

extern "C" double orient2d(const double*, const double*, const double*);


namespace S4r{

inline double Orient2D(const Vec2 &a, const Vec2 &b, const Vec2 &c){
	const double pa[2] = { a[0], a[1] };
	const double pb[2] = { b[0], b[1] };
	const double pc[2] = { c[0], c[1] };
	return orient2d(pa, pb, pc);
}
struct ConvexPolygon{
	std::vector<Vec2> v;
	Vec2 offset;
	ConvexPolygon():offset(0,0){}
	ConvexPolygon(size_t n, double *p);
	double Area() const;
	Vec2 ApproxCenter() const;
};

struct Shape{ // abstract base
	Vec2 center;
	double angle;
	int tag;
	
	Shape(const Vec2 &center, const double &angle, int tag):
		center(center), angle(angle), tag(tag)
	{}
	virtual ~Shape(){}
	
	virtual Vec2 Normal(const Vec2 &r) const = 0;
	virtual bool Inside(const Vec2 &r) const = 0;
	virtual double Area() const = 0;
	virtual Vec2 GetSomeInteriorPoint() const{ return center; }
	virtual double OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const = 0;
	virtual void OutputPostscript(FILE *fp) const = 0;
};

// A circle is defined by |r-c|^2 = 1
struct ShapeCircle : public Shape{
	double radius;
	
	ShapeCircle(const Vec2 &center, const double &radius, int tag);
	Vec2 Normal(const Vec2 &r) const;
	bool Inside(const Vec2 &r) const;
	double Area() const;
	double OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const;
	void OutputPostscript(FILE *fp) const;
};

// An ellipse is defined by |B(r-c)|^2 = 1, where B = inv(diag(halfwidth)) * Q,
// where Q is the rotation matrix corresponding to angle.
struct ShapeEllipse : public Shape{
	Vec2 halfwidth;
	
	ShapeEllipse(const Vec2 &center, const double &angle, const Vec2 &halfwidth, int tag);
	Vec2 Normal(const Vec2 &r) const;
	bool Inside(const Vec2 &r) const;
	double Area() const;
	double OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const;
	void OutputPostscript(FILE *fp) const;
private:
	Mat2 Q;
};

// A rectangle is defined by norminf[B(r-c)] = 1, where B = inv(diag(halfwidth)) * Q,
// where Q is the rotation matrix corresponding to angle.
struct ShapeRectangle : public Shape{
	Vec2 halfwidth;
	
	ShapeRectangle(const Vec2 &center, const double &angle, const Vec2 &halfwidth, int tag);
	Vec2 Normal(const Vec2 &r) const;
	bool Inside(const Vec2 &r) const;
	double Area() const;
	double OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const;
	void OutputPostscript(FILE *fp) const;
private:
	Mat2 Q;
};

struct ShapePolygon : public Shape{
	std::vector<Vec2> vertex;
	
	ShapePolygon(const Vec2 &center, const double &angle, const std::vector<Vec2> &vertex, int tag);
	Vec2 Normal(const Vec2 &r) const;
	bool Inside(const Vec2 &r) const;
	double Area() const;
	Vec2 GetSomeInteriorPoint() const;
	double OverlapTriangle(const Vec2 &p0, const Vec2 &p1p0, const Vec2 &p2p0) const;
	void OutputPostscript(FILE *fp) const;
private:
	std::vector<size_t> t; // list of triangles
};

class Pattern{
	std::vector<Shape*> shape;
	std::vector<size_t> parent;
	bool finalized;
public:
	Pattern():finalized(false){
	}
	~Pattern();
	void AddShape(Shape *s){ shape.push_back(s); finalized = false; }
	size_t NumShapes() const{ return shape.size(); }
	void RemoveShapes(){ shape.clear(); finalized = false; }
	const Shape& GetShape(size_t i) const{ return *(shape[i]); }
	void Finalize(); // recomputes the containment tree
	
	size_t Overlap(
		const ConvexPolygon &poly,
		std::vector<double> &value
	) const;
	
	void OutputPostscript(FILE *fp) const;
};

} // namespace S4r

#endif // S4R_SHAPE_HPP_INCLUDED
