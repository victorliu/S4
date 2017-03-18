#include "Patterning.hpp"
#include "SpecialFunction.hpp"
#include <limits>
#include <iostream>

typedef Patterning::real_type real_type;
typedef Patterning::complex_type complex_type;

/*********************************************************************/
/**************************** Patterning *****************************/
/*********************************************************************/

Patterning::Patterning():t2v(NULL),inct2v(0){}

const complex_type& Patterning::TagToValue(int tag) const{
	static const complex_type nan(
		std::numeric_limits<real_type>::quiet_NaN(),
		std::numeric_limits<real_type>::quiet_NaN()
	);
	if(NULL == t2v){
		return nan;
	}
	if(tag < 0){
		return t2v[0];
	}
	tag++;
	return t2v[tag*inct2v];
}

void Patterning::SetTagToValueMap(const complex_type *m, int incm){
	t2v = m;
	inct2v = incm;
}

/*********************************************************************/
/************************ PatterningByShapes *************************/
/*********************************************************************/

PatterningByShapes::PatterningByShapes(){}

PatterningByShapes::~PatterningByShapes(){
	for(
		std::vector<Shape*>::const_iterator si = shape.begin();
		si != shape.end();
		++si
	){
		delete (*si);
	}
}

void PatterningByShapes::AddShape(Shape *s, int tag){
	real_type p[2];
	s->Center(&p[0]);
	s->tag = tag;
	s->parent = -1;
	for(int i = (int)shape.size()-1; i >= 0; --i){
		if(shape[i]->Contains(p)){
			s->parent = i;
			break;
		}
	}
	shape.push_back(s);
}

void PatterningByShapes::FourierSeries(
	const real_type *Lk, int nik, const int *ik, complex_type *fc
) const{
	real_type *f = (real_type*)malloc(sizeof(real_type) * 2 * nik);
	complex_type *fc1 = (complex_type*)malloc(sizeof(complex_type) * nik);
	for(int i = 0; i < nik; ++i){
		fc[i] = TagToValue(-1);
	}
	for(
		std::vector<Shape*>::const_iterator si = shape.begin();
		si != shape.end();
		++si
	){
		const Shape *s = (*si);
		const int itag = s->tag;
		real_type cs = 1, sn = 0;
		if(0 != s->angle_frac){
			SpecialFunction::CosSin2pi(s->angle_frac, &cs, &sn);
		}
		for(int i = 0; i < nik; ++i){
			const real_type f0[2] = {
				Lk[0] * ik[2*i+0] + Lk[2] * ik[2*i+1],
				Lk[1] * ik[2*i+0] + Lk[3] * ik[2*i+1]
			};
			f[2*i+0] = f0[0] * cs + f0[1] * sn;
			f[2*i+1] = f0[1] * cs - f0[0] * sn;
		}
		complex_type deps = TagToValue(itag);
		const int iparent = s->parent;
		int iparent_tag = -1;
		if(iparent >= 0){
			iparent_tag = shape[iparent]->tag;
		}
		deps -= TagToValue(iparent_tag);
		s->BaseFourierTransform(nik, f, fc1);
		
		for(int i = 0; i < nik; ++i){
			cs = 1;
			sn = 0;
			if(real_type(0) != s->center[0] || real_type(0) != s->center[1]){
				const real_type f0[2] = {
					Lk[0] * ik[2*i+0] + Lk[2] * ik[2*i+1],
					Lk[1] * ik[2*i+0] + Lk[3] * ik[2*i+1]
				};
				SpecialFunction::CosSin2pi(-(f0[0]*s->center[0] + f0[1]*s->center[1]), &cs, &sn);
			}
			complex_type shift_phase(cs, sn);
			fc[i] += deps * shift_phase * fc1[i];
		}
	}
	free(fc1);
	free(f);
}

/*********************************************************************/
/********************** PatterningByBreakpoints **********************/
/*********************************************************************/

real_type fracpart(const real_type &x){
	double dum;
	return modf(x, &dum);
}

PatterningByBreakpoints::PatterningByBreakpoints(const real_type &L):
	period(L),
	iperiod(real_type(1) / L)
{
}
void PatterningByBreakpoints::SetRegion(
	const real_type &center, const real_type &halfwidth, int tag
){
	const real_type x0 = fracpart((center - halfwidth) * iperiod);
	real_type xw = 2*halfwidth*iperiod;
	if(xw > 1){ xw = 1; }
	if(0 == reg.size()){
		real_type x1 = fracpart(x0 + xw);
		if(x1 > x0){
			reg.push_back(PatterningByBreakpoints::Region(x0, tag));
			reg.push_back(PatterningByBreakpoints::Region(x1, -1));
		}else{
			reg.push_back(PatterningByBreakpoints::Region(x1, -1));
			reg.push_back(PatterningByBreakpoints::Region(x0, tag));
		}
		return;
	}
	// TODO: finish this
	// Find where to insert x0
}

/*********************************************************************/
/****************************** Shape ********************************/
/*********************************************************************/

PatterningByShapes::Shape::Shape(
	const real_type *c, const real_type &ang
):angle_frac(ang){
	if(NULL != c){
		center[0] = c[0];
		center[1] = c[1];
	}else{
		center[0] = 0;
		center[1] = 0;
	}
}

bool PatterningByShapes::Shape::Contains(const real_type *p) const{
	real_type cs, sn;
	SpecialFunction::CosSin2pi(angle_frac, &cs, &sn);
	const real_type p0[2] = {
		p[0] - center[0],
		p[1] - center[1]
	};
	const real_type p1[2] = {
		p0[0] * cs + p0[1] * sn,
		p0[1] * cs - p0[0] * sn
	};
	return BaseContains(p1);
}
void PatterningByShapes::Shape::Center(real_type *c) const{
	real_type c0[2];
	BaseCenter(&c0[0]);
	real_type cs, sn;
	SpecialFunction::CosSin2pi(angle_frac, &cs, &sn);
	c[0] = center[0] + c0[0] * cs - c0[1] * sn;
	c[1] = center[1] + c0[1] * cs + c0[0] * sn;
}
void PatterningByShapes::Shape::Normal(const real_type *p, real_type *n) const{
	real_type cs, sn;
	SpecialFunction::CosSin2pi(angle_frac, &cs, &sn);
	const real_type p0[2] = {
		p[0] - center[0],
		p[1] - center[1]
	};
	const real_type p1[2] = {
		p0[0] * cs + p0[1] * sn,
		p0[1] * cs - p0[0] * sn
	};
	real_type n0[2];
	BaseNormal(p1, &n0[0]);
	n[0] = n0[0] * cs - n0[1] * sn;
	n[1] = n0[1] * cs + n0[0] * sn;
}

/*********************************************************************/
/*************************** ShapeHalfwidths *************************/
/*********************************************************************/

PatterningByShapes::ShapeHalfwidths::ShapeHalfwidths(
	const real_type &hx, const real_type &hy,
	const real_type *c, const real_type &ang
):Shape(c, ang){
	halfwidth[0] = hx;
	halfwidth[1] = hy;
}

void PatterningByShapes::ShapeHalfwidths::BaseCenter(real_type *center) const{
	center[0] = 0;
	center[1] = 0;
}

/*********************************************************************/
/****************************** Circle *******************************/
/*********************************************************************/

PatterningByShapes::Circle::Circle(
	const real_type &radius, const real_type *c
):ShapeHalfwidths(radius, 0, c, 0){
}

real_type PatterningByShapes::Circle::Area() const{
	return M_PI*halfwidth[0]*halfwidth[0];
}
bool PatterningByShapes::Circle::BaseContains(const real_type *p) const{
	return hypot(p[0], p[1]) < halfwidth[0];
}
void PatterningByShapes::Circle::BaseNormal(const real_type *p, real_type *n) const{
	real_type a = hypot(p[0], p[1]);
	if(0 == a){
		n[0] = 1;
		n[1] = 0;
	}else{
		a = real_type(1) / a;
		n[0] = p[0] * a;
		n[1] = p[1] * a;
	}
}
void PatterningByShapes::Circle::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	const real_type area = Area();
	for(int i = 0; i < nf; ++i){
		fc[i] = area * SpecialFunction::FourierJinc(halfwidth[0] * hypot(f[2*i+0], f[2*i+1]));
	}
}
real_type PatterningByShapes::Circle::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/***************************** Ellipse *******************************/
/*********************************************************************/

PatterningByShapes::Ellipse::Ellipse(
	const real_type *radii, const real_type *c, const real_type &ang
):ShapeHalfwidths(radii[0], radii[1], c, ang){
}

real_type PatterningByShapes::Ellipse::Area() const{
	return M_PI*halfwidth[0]*halfwidth[1];
}
bool PatterningByShapes::Ellipse::BaseContains(const real_type *x) const{
	const real_type r = halfwidth[0]*halfwidth[0] - halfwidth[1]*halfwidth[1];
	real_type f; // focus length
	real_type L; // twice the semi-major axis length
	real_type d;
	if(r >= 0){
		f = sqrt(r);
		L = 2 * halfwidth[0];
		d = hypot(x[0]-f, x[1]) + hypot(x[0]+f, x[1]);
	}else{
		f = sqrt(-r);
		L = 2 * halfwidth[1];
		d = hypot(x[1]-f, x[0]) + hypot(x[1]+f, x[0]);
	}
	return d < L;
}
void PatterningByShapes::Ellipse::BaseNormal(const real_type *p, real_type *n) const{
	/* The ellipse is defined by the equation
	 *   Norm_2(B.{x,y})^2 == 1
	 * where
	 *   B = diag{1/halfwidth[0],1/halfwidth[1]}.[ cos(angle) sin(angle) ]
	 *                                           [-sin(angle) cos(angle) ]
	 * let ca = cos(angle), sa = sin(angle), ilx2 = 1/lx^2, ily2 = 1/ly^2
	 * Let A = B^T B = [ ca*ca*ilx2 + sa*sa*ily2  ca*sa*(ilx2-ily2) ]
	 *                 [ same as other      ca*ca*ily2 + sa*sa*ilx2 ]
	 * The gradient is then A.{x,y}
	 */
	n[0] = p[0] / (halfwidth[0]*halfwidth[0]);
	n[1] = p[1] / (halfwidth[1]*halfwidth[1]);
}
void PatterningByShapes::Ellipse::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	const real_type area = Area();
	if(halfwidth[0] >= halfwidth[1]){
		const real_type r = halfwidth[1] / halfwidth[0];
		for(int i = 0; i < nf; ++i){
			fc[i] = area * SpecialFunction::FourierJinc(halfwidth[0]*hypot(f[2*i+0], r*f[2*i+1]));
		}
	}else{
		const real_type r = halfwidth[0] / halfwidth[1];
		for(int i = 0; i < nf; ++i){
			fc[i] = area * SpecialFunction::FourierJinc(halfwidth[1]*hypot(f[2*i+1], r*f[2*i+0]));
		}
	}
}
real_type PatterningByShapes::Ellipse::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/**************************** Rectangle ******************************/
/*********************************************************************/

PatterningByShapes::Rectangle::Rectangle(
	const real_type *halfwidths, const real_type *c, const real_type &ang
):ShapeHalfwidths(halfwidths[0], halfwidths[1], c, ang){
}

real_type PatterningByShapes::Rectangle::Area() const{
	return 4*halfwidth[0]*halfwidth[1];
}
bool PatterningByShapes::Rectangle::BaseContains(const real_type *x) const{
	return (fabs(x[0]) < halfwidth[0]) && (fabs(x[1]) < halfwidth[1]);
}
void PatterningByShapes::Rectangle::BaseNormal(const real_type *p, real_type *n) const{
	/* The rectangle is defined by the equation
	 *   max(abs(px)-halfwidth[0], abs(py)-halfwidth[1]) == 0
	 */
	if((fabs(p[0])-halfwidth[0]) > (fabs(p[1])-halfwidth[1])){
		n[0] = (p[0] > 0 ? 1 : -1);
		n[1] = 0;
	}else{
		n[0] = 0;
		n[1] = (p[1] > 0 ? 1 : -1);
	}
}
void PatterningByShapes::Rectangle::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	for(int i = 0; i < nf; ++i){
		fc[i] = SpecialFunction::FourierSinc(2*halfwidth[0]*f[2*i+0]) * SpecialFunction::FourierSinc(2*halfwidth[1]*f[2*i+1]);
	}
}
real_type PatterningByShapes::Rectangle::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/************************** ShapeVertices ****************************/
/*********************************************************************/

PatterningByShapes::ShapeVertices::ShapeVertices(
	int n, const real_type *vxy, int incv,
	const real_type *c, const real_type &ang
):Shape(c, ang), n(n), incv(incv){
	v = (real_type*)malloc(sizeof(real_type)*incv*n);
	memcpy(v, vxy, sizeof(real_type)*incv*n);
}
PatterningByShapes::ShapeVertices::~ShapeVertices(){
	free(v);
}

/*********************************************************************/
/***************************** Polygon *******************************/
/*********************************************************************/

PatterningByShapes::Polygon::Polygon(
	int nv, const real_type *vxy, const real_type *center, const real_type &ang
):ShapeVertices(n, vxy, 2, center, ang){
}

real_type PatterningByShapes::Polygon::Area() const{
	real_type a = 0;
	for(int i = n-1, j = 0; j < n; i=j++){
		a += (v[2*i+0]*v[2*j+1] - v[2*i+1]*v[2*j+0]);
	}
	return (real_type(1) / real_type(2))*a;
}
bool PatterningByShapes::Polygon::BaseContains(const real_type *x) const{
	bool c = false;
	for(int i = n-1, j = 0; j < n; i=j++){
		real_type vix = v[2*i+0];
		real_type viy = v[2*i+1];
		real_type vjx = v[2*j+0];
		real_type vjy = v[2*j+1];
		if(
			((viy>x[1]) != (vjy>x[1]))
			&& (x[0] < (vjx-vix) * (x[1]-viy) / (vjy-viy) + vix)
		){ c = !c; }
	}
	return c;
}
void PatterningByShapes::Polygon::BaseCenter(real_type *center) const{
	center[0] = 0;
	center[1] = 0;
}
void PatterningByShapes::Polygon::BaseNormal(const real_type *p, real_type *nv) const{
	int i, j;
	real_type maxdist = -1;
	for(j = 0, i = n-1; j < n; i = j++){
		/* compute distance from r to segment */
		real_type vv[2], pr[2];
		vv[0] = v[2*j+0] - v[2*i+0];
		vv[1] = v[2*j+1] - v[2*i+1];

		pr[0] = p[0] - v[2*i+0];
		pr[1] = p[1] - v[2*i+1];
		{
			const real_type vv2 = vv[0]*vv[0] + vv[1]*vv[1];
			const real_type prj = (pr[0]*vv[0] + pr[1]*vv[1])/vv2;
			const real_type voff[2] = {pr[0] - prj*vv[0], pr[1] - prj*vv[1]};
			const real_type dist = hypot(voff[0], voff[1]);
			if(dist > maxdist){
				maxdist = dist;
				nv[0] = vv[1];
				nv[1] = -vv[0];
			}
		}
	}
}
void PatterningByShapes::Polygon::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	for(int i = 0; i < nf; ++i){
		if(0 == f[2*i+0] && 0 == f[2*i+1]){
			fc[i] = Area();
			continue;
		}
		const real_type ff = f[2*i+0]*f[2*i+0] + f[2*i+1]*f[2*i+1];
		/* For k != 0,
		 * S(k) = i/|k|^2 * Sum_{i=0,n-1} z.((v_{i+1}-v_{i}) x k) j0(k.(v_{i+1}-v_{i})/2) e^{ik.(v_{i+1}+v_{i})/2}
		 */
		int p,q;
		real_type num, cs, sn;
		real_type rc[2], u[2];
		for(int p = n-1, q = 0; q < n; p = q++){
			u[0] = v[2*q+0]-v[2*p+0];
			u[1] = v[2*q+1]-v[2*p+1];
			rc[0] = 0.5*(v[2*q+0]+v[2*p+0]);
			rc[1] = 0.5*(v[2*q+1]+v[2*p+1]);

			num = (u[0]*f[2*i+1]-u[1]*f[2*i+0]) * SpecialFunction::FourierSinc(f[2*i+0]*u[0]+f[2*i+1]*u[1]);
			SpecialFunction::CosSin2pi(-(f[2*i+0]*rc[0]+f[2*i+1]*rc[1]), &cs, &sn);

			// Multiplication by i means we mess up the order here
			fc[i] += num * complex_type(sn, -cs);
		}
		// Our f lacks a 2pi factor
		fc[i] /= 2*M_PI*ff;
	}
}
real_type PatterningByShapes::Polygon::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}
