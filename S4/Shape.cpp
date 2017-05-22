#include "Patterning.hpp"
#include "Shape.hpp"
#include <limits>
#include <iostream>

typedef Patterning::real_type real_type;
typedef Patterning::complex_type complex_type;

#include "SpecialFunction.hpp"

/*********************************************************************/
/****************************** Shape ********************************/
/*********************************************************************/

Shape::Shape(
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

bool Shape::Contains(const real_type *p) const{
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
void Shape::Center(real_type *c) const{
	real_type c0[2];
	BaseCenter(&c0[0]);
	real_type cs, sn;
	SpecialFunction::CosSin2pi(angle_frac, &cs, &sn);
	c[0] = center[0] + c0[0] * cs - c0[1] * sn;
	c[1] = center[1] + c0[1] * cs + c0[0] * sn;
}
void Shape::Normal(const real_type *p, real_type *n) const{
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

ShapeHalfwidths::ShapeHalfwidths(
	const real_type &hx, const real_type &hy,
	const real_type *c, const real_type &ang
):Shape(c, ang){
	halfwidth[0] = hx;
	halfwidth[1] = hy;
}

void ShapeHalfwidths::BaseCenter(real_type *center) const{
	center[0] = 0;
	center[1] = 0;
}

/*********************************************************************/
/****************************** Circle *******************************/
/*********************************************************************/

Circle::Circle(
	const real_type &radius, const real_type *c
):ShapeHalfwidths(radius, 0, c, 0){
}
Circle* Circle::Clone() const{
	return new Circle(halfwidth[0], center);
}

real_type Circle::Area() const{
	return M_PI*halfwidth[0]*halfwidth[0];
}
bool Circle::BaseContains(const real_type *p) const{
	return hypot(p[0], p[1]) < halfwidth[0];
}
void Circle::BaseNormal(const real_type *p, real_type *n) const{
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
void Circle::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	const real_type area = Area();
	for(int i = 0; i < nf; ++i){
		fc[i] = area * SpecialFunction::FourierJinc(halfwidth[0] * hypot(f[2*i+0], f[2*i+1]));
	}
}
real_type Circle::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/***************************** Ellipse *******************************/
/*********************************************************************/

Ellipse::Ellipse(
	const real_type *radii, const real_type *c, const real_type &ang
):ShapeHalfwidths(radii[0], radii[1], c, ang){
}
Ellipse* Ellipse::Clone() const{
	return new Ellipse(halfwidth, center, angle_frac);
}

real_type Ellipse::Area() const{
	return M_PI*halfwidth[0]*halfwidth[1];
}
bool Ellipse::BaseContains(const real_type *x) const{
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
void Ellipse::BaseNormal(const real_type *p, real_type *n) const{
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
void Ellipse::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
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
real_type Ellipse::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/**************************** Rectangle ******************************/
/*********************************************************************/

Rectangle::Rectangle(
	const real_type *halfwidths, const real_type *c, const real_type &ang
):ShapeHalfwidths(halfwidths[0], halfwidths[1], c, ang){
}
Rectangle* Rectangle::Clone() const{
	return new Rectangle(halfwidth, center, angle_frac);
}

real_type Rectangle::Area() const{
	return 4*halfwidth[0]*halfwidth[1];
}
bool Rectangle::BaseContains(const real_type *x) const{
	return (fabs(x[0]) < halfwidth[0]) && (fabs(x[1]) < halfwidth[1]);
}
void Rectangle::BaseNormal(const real_type *p, real_type *n) const{
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
void Rectangle::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	const real_type area = Area();
	for(int i = 0; i < nf; ++i){
		fc[i] = area * SpecialFunction::FourierSinc(2*halfwidth[0]*f[2*i+0]) * SpecialFunction::FourierSinc(2*halfwidth[1]*f[2*i+1]);
	}
}
real_type Rectangle::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}

/*********************************************************************/
/************************** ShapeVertices ****************************/
/*********************************************************************/

ShapeVertices::ShapeVertices(
	int n, const real_type *vxy, int incv,
	const real_type *c, const real_type &ang
):Shape(c, ang), n(n), incv(incv){
	v = (real_type*)malloc(sizeof(real_type)*incv*n);
	memcpy(v, vxy, sizeof(real_type)*incv*n);
}
ShapeVertices::~ShapeVertices(){
	free(v);
}

/*********************************************************************/
/***************************** Polygon *******************************/
/*********************************************************************/

Polygon::Polygon(
	int nv, const real_type *vxy, const real_type *center, const real_type &ang
):ShapeVertices(nv, vxy, 2, center, ang){
}
Polygon* Polygon::Clone() const{
	return new Polygon(n, v, center, angle_frac);
}

real_type Polygon::Area() const{
	real_type a = 0;
	for(int i = n-1, j = 0; j < n; i=j++){
		a += (v[2*i+0]*v[2*j+1] - v[2*i+1]*v[2*j+0]);
	}
	return (real_type(1) / real_type(2))*a;
}

bool Polygon::BaseInversionSymmetric() const{
	if(0 != n % 2){ return false; }
	const int h = n/2;
	for(int i = 0; i < h; ++i){
		const int j = n-i-1;
		if(v[2*i+0] != -v[2*j+0] || v[2*i+1] != -v[2*j+1]){ return false; }
	}
	return true;
}
bool Polygon::BaseContains(const real_type *x) const{
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
void Polygon::BaseCenter(real_type *center) const{
	center[0] = 0;
	center[1] = 0;
}
void Polygon::BaseNormal(const real_type *p, real_type *nv) const{
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
void Polygon::BaseFourierTransform(int nf, const real_type *f, complex_type *fc) const{
	for(int i = 0; i < nf; ++i){
		if(0 == f[2*i+0] && 0 == f[2*i+1]){
			fc[i] = Area();
			continue;
		}
		const real_type ff = f[2*i+0]*f[2*i+0] + f[2*i+1]*f[2*i+1];
		/* For k != 0,
		 * S(k) = i/|k|^2 * Sum_{i=0,n-1} z.((v_{i+1}-v_{i}) x k) j0(k.(v_{i+1}-v_{i})/2) e^{ik.(v_{i+1}+v_{i})/2}
		 */
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
real_type Polygon::BaseIntersectTriangle(const real_type *p0, const real_type *u, const real_type *v, real_type *crossuv) const{
	return 0;
}
