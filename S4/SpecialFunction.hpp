#ifndef SPECIAL_FUNCTION_HPP_INCLUDED
#define SPECIAL_FUNCTION_HPP_INCLUDED

#include <cmath>

namespace SpecialFunction{

template <typename T>
void CosSin2pi(T x, T *cs, T *sn){
	if(0 == x){
		*cs = 1;
		*sn = 0;
	}
	T intpart, fracpart = std::modf(x, &intpart);
	if(fracpart < 0){ fracpart += 1; }
	fracpart *= 2*M_PI;
	*cs = std::cos(fracpart);
	*sn = std::sin(fracpart);
}

template <typename T>
T FourierSinc(T x){
	static const T zero(0);
	static const T one(1);
	static const T sixth(one/T(6));
	static const T half(one/T(2));
	T intpart;
	// Defined as Sin[pi x] / (pi x)
	// Handle integer case
	if(zero == x){ return one; }
	if(x < zero){ x = -x; }
	const T fracpart = modf(x, &intpart);
	if(zero == fracpart){
		return zero;
	}
	// Handle the half-integer case
	if(half == fracpart){
		T a = M_2_PI / (1+2*intpart);
		if(0 != ((int)(intpart)) % 2){
			a = -a;
		}
		return a;
	}
	// Handle the small argument case
	x *= M_PI;
	if(x < 2.44140625e-4){
		return one - sixth*x*x;
	}
	// Handle the general case
	return std::sin(x)/x;
}

template <typename T>
T FourierJinc(T x){
	static const T zero(0);
	static const T one(1);
	static const T eighth(one/T(8));
	// Defined as 2*BesselJ[1, 2 pi x] / (2 pi x)
	if(zero == x){ return one; }
	if(x < zero){ x = -x; }
	x *= (2*M_PI);
	if(x < 5e-5){
		return one - eighth*x*x;
	}
	// The following is from crbond.com
	static const double a1[] = {
		 0.1171875,
		-0.1441955566406250,
		 0.6765925884246826,
		-6.883914268109947,
		 1.215978918765359e2,
		-3.302272294480852e3,
		 1.276412726461746e5,
		-6.656367718817688e6,
		 4.502786003050393e8,
		-3.833857520742790e10,
		 4.011838599133198e12,
		-5.060568503314727e14,
		 7.572616461117958e16,
		-1.326257285320556e19};
	static const double b1[] = {
		-0.1025390625,
		 0.2775764465332031,
		-1.993531733751297,
		 2.724882731126854e1,
		-6.038440767050702e2,
		 1.971837591223663e4,
		-8.902978767070678e5,
		 5.310411010968522e7,
		-4.043620325107754e9,
		 3.827011346598605e11,
		-4.406481417852278e13,
		 6.065091351222699e15,
		-9.833883876590679e17,
		 1.855045211579828e20};

	if(x <= 12.0){
		double x2 = x*x;
		double j1 = 1.0;
		double r = 1.0;
		int k;
		for(k=1;k<=30;k++){
			r *= -0.25*x2/(k*(k+1));
			j1 += r;
			if (fabs(r) < fabs(j1)*1e-15) break;
		}
		return j1;
	}else{
		int kz;
		if (x >= 50.0) kz = 8;          /* Can be changed to 10 */
		else if (x >= 35.0) kz = 10;    /*  "       "        12 */
		else kz = 12;                   /*  "       "        14 */
		double cu = sqrt(M_2_PI/x);
		double t2 = x-0.75*M_PI;
		double p1 = 1.0;
		double q1 = 0.375/x;
		int k;
		for(k=0;k<kz;k++){
			p1 += a1[k]*pow(x,-2*k-2);
			q1 += b1[k]*pow(x,-2*k-3);
		}
		return 2.0*cu*(p1*cos(t2)-q1*sin(t2))/x;
	}
}

} // namespace SpecialFunction

#endif // SPECIAL_FUNCTION_HPP_INCLUDED
