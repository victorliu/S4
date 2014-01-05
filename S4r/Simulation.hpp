#ifndef S4R_SIMULATION_HPP_INCLUDED
#define S4R_SIMULATION_HPP_INCLUDED

#include <S4r/Types.hpp>
#include <vector>
#include <map>

#include <S4r/PeriodicMesh.hpp>
#include <S4r/Layer.hpp>
#include <S4r/Material.hpp>

namespace S4r{

struct Simulation{
	PeriodicMesh *mesh;
	
	typedef std::map<std::string, size_t> LayerMap;
	typedef std::map<std::string, size_t> MaterialMap;
	typedef std::vector<Layer*> LayerList;
	typedef std::vector<Material*> MaterialList;
	
	MaterialList materials;
	MaterialMap matmap;
	LayerList layers;
	LayerMap layer_name_to_index;
	
	doublecomplex omega;
	struct{
		double k[2]; // in the reciprocal vector basis, not in cartesian
		double pol[2];
	} excitation;
	
	struct{
		int verbosity;
		size_t num_modes; // 0 for maximal
	} options;
	
	// Constructor/Destructor
	Simulation(const double Lr[4], const size_t n[2], size_t nmodes = 0);
	~Simulation();
	
	// Methods
	int SetMaterial(const std::string &name, Material* matobj);
	Layer* AddLayer(const std::string &name, const double &thickness, const std::string &material);
	void SetK(const double k[2]);
	void SetFrequency(const doublecomplex &freq);
	
	int SetRegionCircle   (const std::string &layer, const std::string &matname, const Vec2 &center, double radius);
	int SetRegionEllipse  (const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const Vec2 &halfwidths);
	int SetRegionRectangle(const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const Vec2 &halfwidths);
	int SetRegionPolygon  (const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const std::vector<Vec2> &vert);
	int RemoveLayerRegions(const std::string &layer);
	int SetLayerThickness(const std::string &layer, const double &thickness);
	
	void GetPowerFlux(const std::string &layer, const double &z, doublecomplex &forw, doublecomplex &back);
	void GetFields(const Vec3 &r, CVec3 &E, CVec3 &H);
	
	void SolveAllLayers();
	void OutputLayerDescription(const std::string &layer, FILE *fp) const;
	
	CVec2 GetBlochPhaseFactors() const;
private:
	bool layer_modes_computed;

	int MaterialNameToIndex(const std::string &matname) const;
	int LayerNameToIndex(const std::string &name) const;
	
	void BuildExcitation();
	void SolveLayer(size_t i);
	int GetSMatrix(
		size_t layer_start, size_t layer_end,
		CMat &S
	);
	void ComputeLayerModes();
	void ResetLayerModes();
	void SolveAllLayers(size_t i0, size_t iN);
	void ResetLayerSolutions();
};

} // namespace S4r

#endif // S4R_SIMULATION_HPP_INCLUDED
