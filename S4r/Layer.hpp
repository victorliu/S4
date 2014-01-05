#ifndef S4R_LAYER_HPP_INCLUDED
#define S4R_LAYER_HPP_INCLUDED

#include <S4r/Types.hpp>
#include <map>
extern "C" {
#include <S4r/Shape.hpp>
}

namespace S4r{

class PeriodicMesh;
struct Material;

struct Layer{
	struct{
		static const bitfield PEC = 0x1000;
		static const bitfield PMC = 0x2000;
		
		std::string name;
		double thickness;
		std::string material;
		Pattern pattern;
		std::string copy;
		bitfield flags;
		size_t num_modes; // 0 for default (maximal)
	} description;
	
	struct{
		bool geom_force_inexact; // avoid exact geometry overlap computations
		size_t geom_num_samples; // number of samples to use when averaging material constants (actual number may go as square of this number)
		
		bool use_tvf; // use tangent vector field for polarization decomposition
	} options;
	
	struct{
		RMat tvf; // tangent vector field on the mesh
	} coord;
	
	struct Matrices{
		static const bitfield DIAGONAL  = 0x1;
		static const bitfield HERMITIAN = 0x2;
		static const bitfield SCALAR    = 0x8;
		SpMat eps, mu, ieps, imu;
		double eps_unif, mu_unif;
		bitfield flags;
		bool valid;
	} matrices;

	struct{
		size_t n; // number of modes
		CVec  q;
		SpMat kp;
		CMat  phi;
		bool valid;
	} modes;
	
	struct Solution{
		static const bitfield A_SET = 0x1;
		static const bitfield B_SET = 0x2;
		CVec a, b;
		bitfield flags;
	} solution;
	
	// Constructor/Destructor
	Layer(const std::string &name);
	~Layer();
	
	// Methods
	size_t GetNumModes() const;
	
	void GetEpsilon(const Vec2 &r, CTensor3 &eps) const;
	void GetMu(const Vec2 &r, CTensor3 &mu) const;
	int BuildMatrices(
		const PeriodicMesh *mesh,
		const std::map<std::string, size_t> &matmap,
		const std::vector<Material*> &mat
	);
	int ComputeModes(
		const doublecomplex &omega,
		const PeriodicMesh *mesh
	);
	void AmplitudesToField(
		const doublecomplex &omega,
		const PeriodicMesh *mesh,
		const CVec &a, const CVec &b,
		CVec &e, CVec &h
	) const;
	void AmplitudesAtZ( // assumes solution a and b are set
		const double &z,
		CVec &a, CVec &b
	) const;
	void GetPowerFlux(
		const doublecomplex &omega,
		const PeriodicMesh *mesh,
		const double &z,
		doublecomplex &forward, doublecomplex &backward
	) const;
	void GetFields(
		const doublecomplex &omega,
		const PeriodicMesh *mesh,
		const CVec2 &phi,
		const Vec3 &r,
		CVec3 &E, CVec3 &H
	) const;
	// Force a BuildMatrices() call before calling this
	// to also overlay the discretization
	void DumpDescription(
		const PeriodicMesh *mesh,
		FILE *fp
	) const;
private:
	bool IsUniformHermitian() const;
	void GetMaterialAverage(
		const ConvexPolygon &poly,
		CTensor2 &eps, CTensor2 &mu,
		const std::vector<Material*> &mat,
		const Material *bkmat,
		std::vector<double> &value // workspace
	);
	void GetTranslatedSolution(const double &z, CVec &a, CVec &b) const;
};

} // namespace S4r

#endif // S4R_LAYER_HPP_INCLUDED