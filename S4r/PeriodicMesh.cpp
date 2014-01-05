#include <S4r/PeriodicMesh.hpp>

void S4r::PeriodicMesh::Remap(Vec2 &r) const{
	double dummy;
	Vec2 p = Lr.inverse() * r;
	p[0] = modf(p[0], &dummy);
	if(p[0] < -0.5){ p[0] += 1.; }
	if(p[0] >= 0.5){ p[0] -= 1.; }
	p[1] = modf(p[1], &dummy);
	if(p[1] < -0.5){ p[1] += 1.; }
	if(p[1] >= 0.5){ p[1] -= 1.; }
	r = Lr * p;
}
