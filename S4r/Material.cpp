#include <S4r/Material.hpp>

S4r::Material::Material(){
	eps(0,0) = 1.;
	eps(1,0) = 0.;
	eps(2,0) = 0.;
	eps(0,1) = 0.;
	eps(1,1) = 1.;
	eps(2,1) = 0.;
	eps(0,2) = 0.;
	eps(1,2) = 0.;
	eps(2,2) = 1.;
	
	mu(0,0) = 1.;
	mu(1,0) = 0.;
	mu(2,0) = 0.;
	mu(0,1) = 0.;
	mu(1,1) = 1.;
	mu(2,1) = 0.;
	mu(0,2) = 0.;
	mu(1,2) = 0.;
	mu(2,2) = 1.;
	flags |= Material::DIAGONAL;
	flags |= Material::HERMITIAN;
}

S4r::Material::Material(const double meps[18], const double mmu[18]){
	static const doublecomplex zero(0.);
	flags = Material::DIAGONAL | Material::HERMITIAN;
	for(unsigned j = 0; j < 3; ++j){
		for(unsigned i = 0; i < 3; ++i){
			eps(i,j) = doublecomplex(meps[0+(i+j*3)*2], meps[1+(i+j*3)*2]);
			if(i != j && zero != eps(i,j)){
				flags &= ~Material::DIAGONAL;
			}
			if(i < j && eps(i,j) != std::conj(eps(j,i))){
				flags &= ~Material::HERMITIAN;
			}else if(i == j && 0. != eps(i,j).imag()){
				flags &= ~Material::HERMITIAN;
			}
		}
	}
	for(unsigned j = 0; j < 3; ++j){
		for(unsigned i = 0; i < 3; ++i){
			mu(i,j) = doublecomplex(mmu[0+(i+j*3)*2], mmu[1+(i+j*3)*2]);
			if(i != j && zero != mu(i,j)){
				flags &= ~Material::DIAGONAL;
			}
			if(i < j && mu(i,j) != std::conj(mu(j,i))){
				flags &= ~Material::HERMITIAN;
			}else if(i == j && 0. != mu(i,j).imag()){
				flags &= ~Material::HERMITIAN;
			}
		}
	}
	if(flags & Material::DIAGONAL){
		if(eps(0,0) == eps(1,1) && eps(0,0) == eps(2,2) &&
		    mu(0,0) ==  mu(1,1) &&  mu(0,0) ==  mu(2,2)
		){
			flags |= Material::SCALAR;
		}
	}
}

S4r::Material::~Material(){
}
