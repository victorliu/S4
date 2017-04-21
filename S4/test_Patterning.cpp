#include "Patterning.hpp"
#include <iostream>

int main(int argc, char *argv[]){
	PatterningByShapes *pat = new PatterningByShapes();
	PatterningByShapes::Shape *s = new PatterningByShapes::Circle(0.2);
	pat->AddShape(s, 0);
	
	const Patterning::complex_type eps[2] = { 12., 1. };
	pat->SetTagToValueMap(eps, 1);
	
	const double Lk[4] = {1, 0, 0, 1};
	const int nik = 1;
	const int ik[2] = {1, 0};
	Patterning::complex_type fc[1];
	
	pat->FourierSeries(Lk, nik, ik, fc);
	std::cout << "f: " << fc[0] << std::endl;
	delete pat;
	return 0;
}
