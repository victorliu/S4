#ifndef _RNP_IO_H_
#define _RNP_IO_H_

#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>

/* Functions: (templated on T)
 *   std::ostream& Print(const T &value, std::ostream &os = std::cout);
 *   std::ostream& PrintVector(size_t n, const T *x, size_t incx = 1, std::ostream &os = std::cout);
 *   std::ostream& PrintMatrix(size_t m, size_t n, const T *a, size_t lda, std::ostream &os = std::cout);
 *
 * Preprocessor flags:
 *   RNP_OUTPUT_MATHEMATICA or RNP_OUTPUT_MATLAB or <none> - choose the output format.
 */

namespace RNP{
namespace IO{

template <class T>
class _NumberOutputter{
	std::ostream &os;
public:
	_NumberOutputter(std::ostream &o):os(o){}
	std::ostream& operator()(const T &val){
#if defined(RNP_OUTPUT_MATHEMATICA)
		std::stringstream s;
		std::string str;
		s << std::setprecision(17);
		s << val;
		str = s.str();
		size_t pos;
		if((pos = str.find_first_of('e')) != std::string::npos){
			str.replace(pos, 1, "*^");
		}
		return (os << str);
#else
		return (os << val);
#endif
	}
};

template <class T>
class _NumberOutputter<std::complex<T> >{
	std::ostream &os;
public:
	_NumberOutputter(std::ostream &o):os(o){}
	std::ostream& operator()(const std::complex<T> &val){
		_NumberOutputter<T> outputter(os);
#if defined(RNP_OUTPUT_MATHEMATICA)
		if(val.real() != 0 || 0 == val.imag()){
			outputter(val.real());
		}
		if(val.imag() != 0){
			if(val.imag() >= 0){
				if(val.real() != 0){
					os << "+";
				}
				os << "I*";
				outputter(val.imag());
			}else{
				os << "-I*";
				outputter(-val.imag());
			}
		}
#elif defined(RNP_OUTPUT_MATLAB)
		if(val.real() != 0 || 0 == val.imag()){
			outputter(val.real());
		}
		if(val.imag() != 0){
			if(val.imag() >= 0){
				if(val.real() != 0){
					os << "+";
				}
				outputter(val.imag());
			}else{
				outputter(val.imag());
			}
			os << "i";
		}
#else
		os << val;
#endif
		return os;
	}
};

template <class T>
std::ostream& Print(const T &value, std::ostream &os = std::cout){
	_NumberOutputter<T> outputter(os);
	outputter(value);
	return os;
}

template <class T>
std::ostream& PrintVector(size_t n, const T *x, size_t incx = 1, std::ostream &os = std::cout){
	_NumberOutputter<T> outputter(os);

	// Opening format
#if defined(RNP_OUTPUT_MATHEMATICA)
	os << "{";
#elif defined(RNP_OUTPUT_MATLAB)
	os << "[";
#else
#endif

	while(n --> 0){
		outputter(*x);
		x += incx;

		if(n > 0){
			// Separator format
#if defined(RNP_OUTPUT_MATHEMATICA)
			os << ", ";
#elif defined(RNP_OUTPUT_MATLAB)
			os << "; ";
#else
			os << std::endl;
#endif
		}
	}

	// Closing format
#if defined(RNP_OUTPUT_MATHEMATICA)
	os << "}";
#elif defined(RNP_OUTPUT_MATLAB)
	os << "]";
#else
#endif

	return os;
}

template <class T>
std::ostream& PrintMatrix(size_t m, size_t n, const T *a, size_t lda, std::ostream &os = std::cout){
	_NumberOutputter<T> outputter(os);

	// Opening format
#if defined(RNP_OUTPUT_MATHEMATICA)
	os << "{";
#elif defined(RNP_OUTPUT_MATLAB)
	os << "[";
#else
#endif

	for(size_t i = 0; i < m; ++i){
		// Row start format
#if defined(RNP_OUTPUT_MATHEMATICA)
		os << "{";
#elif defined(RNP_OUTPUT_MATLAB)
		os << " ";
#else
#endif
		for(size_t j = 0; j < n; ++j){
			outputter(a[i+j*lda]);

			if(j < n-1){
				// Col separator format
#if defined(RNP_OUTPUT_MATHEMATICA)
				os << ", ";
#elif defined(RNP_OUTPUT_MATLAB)
				os << ", ";
#else
				os << "\t";
#endif
			}
		}
		
		// Row end format
#if defined(RNP_OUTPUT_MATHEMATICA)
		os << "}";
#elif defined(RNP_OUTPUT_MATLAB)
#else
#endif

		if(i < m-1){
			// Row separator format
#if defined(RNP_OUTPUT_MATHEMATICA)
			os << "," << std::endl;
#elif defined(RNP_OUTPUT_MATLAB)
			os << ";" << std::endl;
#else
			os << "\n";
#endif
		}
	}

	// Closing format
#if defined(RNP_OUTPUT_MATHEMATICA)
	os << "}";
#elif defined(RNP_OUTPUT_MATLAB)
	os << "]";
#else
#endif

	return os;
}

}; // namespace IO
}; // namespace RNP

#endif // _RNP_IO_H_
