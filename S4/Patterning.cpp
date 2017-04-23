#include "Patterning.hpp"
#include "Shape.hpp"
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

int PatterningByShapes::AddShape(Shape *s, int tag){
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
	return 0;
}

void PatterningByShapes::FourierSeries(
	const real_type *Lk, int nik, const int *ik, complex_type *fc, real_type *shift
) const{
	real_type *f = (real_type*)malloc(sizeof(real_type) * 2 * nik);
	complex_type *fc1 = (complex_type*)malloc(sizeof(complex_type) * nik);
	for(int i = 0; i < nik; ++i){
		if(0 == ik[2*i+0] && 0 == ik[2*i+1]){
			fc[i] = TagToValue(-1);
		}else{
			fc[i] = S4_real(0);
		}
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
			real_type cen[2] = { s->center[0], s-> center[1] };
			if(NULL != shift){ cen[0] -= shift[0]; cen[1] -= shift[1]; }
			if(real_type(0) != cen[0] || real_type(0) != cen[1]){
				const real_type f0[2] = {
					Lk[0] * ik[2*i+0] + Lk[2] * ik[2*i+1],
					Lk[1] * ik[2*i+0] + Lk[3] * ik[2*i+1]
				};
				SpecialFunction::CosSin2pi(-(f0[0]*cen[0] + f0[1]*cen[1]), &cs, &sn);
			}
			complex_type shift_phase(cs, sn);
			fc[i] += deps * shift_phase * fc1[i];
		}
	}
	free(fc1);
	free(f);
}

PatterningByShapes* PatterningByShapes::Clone() const{
	PatterningByShapes *p = new PatterningByShapes();
	for(std::vector<Shape*>::const_iterator i = shape.begin(); shape.end() != i; ++i){
		p->AddShape((*i)->Clone(), (*i)->tag);
	}
	return p;
}

/*********************************************************************/
/*********************** PatterningByIntervals ***********************/
/*********************************************************************/

PatterningByIntervals::PatterningByIntervals(){}

PatterningByIntervals::~PatterningByIntervals(){
}

int PatterningByIntervals::SetRegion(const real_type &center, const real_type &halfwidth, int tag){
	interval.push_back(Interval(center, halfwidth, tag));
	for(int i = (int)interval.size()-1; i >= 0; --i){
		if(interval[i].Contains(center)){
			interval.back().parent = i;
			break;
		}
	}
	return 0;
}

void PatterningByIntervals::FourierSeries(
	const real_type *Lk, int nik, const int *ik, complex_type *fc, real_type *shift
) const{
	real_type *f = (real_type*)malloc(sizeof(real_type) * nik);
	complex_type *fc1 = (complex_type*)malloc(sizeof(complex_type) * nik);
	for(int i = 0; i < nik; ++i){
		fc[i] = TagToValue(-1);
	}
	for(
		std::vector<Interval>::const_iterator si = interval.begin();
		si != interval.end();
		++si
	){
		const Interval &s = (*si);
		const int itag = s.tag;
		for(int i = 0; i < nik; ++i){
			f[i] = Lk[0] * ik[2*i+0];
		}
		complex_type deps = TagToValue(itag);
		const int iparent = s.parent;
		int iparent_tag = -1;
		if(iparent >= 0){
			iparent_tag = interval[iparent].tag;
		}
		deps -= TagToValue(iparent_tag);
		for(int i = 0; i < nik; ++i){
			fc1[i] = SpecialFunction::FourierSinc(s.halfwidth * f[i]);
		}
		
		for(int i = 0; i < nik; ++i){
			real_type cs = 1, sn = 0;
			real_type cen = s.center;
			if(NULL != shift){ cen -= shift[0]; }
			if(real_type(0) != cen){
				const real_type f0 = Lk[0] * ik[2*i+0];
				SpecialFunction::CosSin2pi(-f0*cen, &cs, &sn);
			}
			complex_type shift_phase(cs, sn);
			fc[i] += deps * shift_phase * fc1[i];
		}
	}
	free(fc1);
	free(f);
}

PatterningByIntervals* PatterningByIntervals::Clone() const{
	PatterningByIntervals *p = new PatterningByIntervals();
	for(std::vector<Interval>::const_iterator i = interval.begin(); interval.end() != i; ++i){
		p->SetRegion(i->center, i->halfwidth, i->tag);
	}
	return p;
}

/*********************************************************************/
/********************** PatterningByBreakpoints **********************/
/*********************************************************************/

static real_type fracpart(const real_type &x){
	double dum;
	return modf(x, &dum);
}

PatterningByBreakpoints::PatterningByBreakpoints(const real_type &L):
	period(L),
	iperiod(real_type(1) / L)
{
}

PatterningByBreakpoints* PatterningByBreakpoints::Clone() const{
	PatterningByBreakpoints *p = new PatterningByBreakpoints(period);
	p->reg = this->reg;
	return p;
}

int PatterningByBreakpoints::SetRegion(
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
		return 0;
	}
	// TODO: finish this
	// Find where to insert x0
	return 0;
}

void PatterningByBreakpoints::FourierSeries(
	const real_type *Lk, int nik, const int *ik, complex_type *fc, real_type *shift
) const{
	return;
}
