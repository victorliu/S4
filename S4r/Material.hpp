#ifndef S4R_MATERIAL_HPP_INCLUDED
#define S4R_MATERIAL_HPP_INCLUDED

#include <S4r/Types.hpp>

namespace S4r{

// This object must be simple enough to be copy constructed without much overhead
struct Material{
	static const bitfield DIAGONAL  = 0x1;
	static const bitfield HERMITIAN = 0x2;
	static const bitfield SCALAR    = 0x8;
	
	CTensor3 eps, mu;
	bitfield flags;
	
	Material();
	Material(const double eps[18], const double mu[18]);
	~Material();
};

} // namespace S4r

#endif // S4R_MATERIAL_HPP_INCLUDED
