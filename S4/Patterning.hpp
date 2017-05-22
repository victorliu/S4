#ifndef PATTERNING_HPP_INCLUDED
#define PATTERNING_HPP_INCLUDED

#include "S4.h"
#include "Shape.hpp"
#include <complex>
#include <vector>
//L:MaterialMap{ image = 'foo.png', map = { ["000000"] = mSi } }

struct Patterning{
	// Represents a generic layer patterning
public:
	typedef S4_real real_type;
	typedef std::complex<real_type> complex_type;
	static const int UNSUPPORTED_PATTERN_TYPE = -999;
public:
	Patterning();
	virtual ~Patterning(){}
	virtual Patterning *Clone() const = 0;
	
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
	virtual void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc, const real_type *shift = NULL) const = 0;
	
	// Returns true if the pattern is inversion symmetric for some shift, and optionally returns that shift.
	// The shift is the location of the point of inversion symmetry. If the pattern is a single circle located at (x,y),
	// the shift[0] = x, shift[1] = y.
	// If shift is NULL, then checks for inversion symmetry under zero shift.
	virtual bool InversionSymmetric(real_type *shift) const{ return false; }
	
	// These routines return -
	virtual int SetRegion(const real_type &center, const real_type &halfwidth, int tag = 0){ return UNSUPPORTED_PATTERN_TYPE; }
	virtual int AddShape(Shape *s, int tag = 0){ return UNSUPPORTED_PATTERN_TYPE; }
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
	PatterningByBreakpoints *Clone() const;
	int SetRegion(const real_type &center, const real_type &halfwidth, int tag = 0);
	void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc, const real_type *shift = NULL) const;
};

class PatterningByShapes : public Patterning{
public:
	PatterningByShapes();
	PatterningByShapes *Clone() const;
	~PatterningByShapes();
	int AddShape(Shape *s, int tag = 0);
	bool InversionSymmetric(real_type *shift) const;
	void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc, const real_type *shift = NULL) const;
protected:
	std::vector<Shape*> shape;
};
class PatterningByIntervals : public Patterning{
public:
	struct Interval{
		real_type center;
		real_type halfwidth;
		int tag;
		int parent;
		Interval(const real_type &center, const real_type &halfwidth, int tag = 0):
			center(center),
			halfwidth(halfwidth),
			tag(tag),
			parent(-1)
		{}
		bool Contains(const real_type &x) const{ return std::abs(x - center) < halfwidth; }
	};
	PatterningByIntervals();
	PatterningByIntervals *Clone() const;
	~PatterningByIntervals();
	int SetRegion(const real_type &center, const real_type &halfwidth, int tag = 0);
	int AddShape(Shape *s, int tag = 0){ return UNSUPPORTED_PATTERN_TYPE; }
	bool InversionSymmetric(real_type *shift) const;
	void FourierSeries(const real_type *Lk, int nik, const int *ik, complex_type *fc, const real_type *shift = NULL) const;
protected:
	std::vector<Interval> interval;
};

class PatterningByMap : public Patterning{
protected:
	std::vector<int> tag;
	int nu, nv;
public:
	PatterningByMap *Clone() const;
};

class PatterningByArray : public Patterning{
protected:
	std::vector<complex_type> val;
	int nu, nv;
public:
	PatterningByArray *Clone() const;
};

#endif // PATTERNING_HPP_INCLUDED
