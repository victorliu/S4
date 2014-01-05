#include <S4r/Simulation.hpp>
#include <S4r/Pseudoinverse.hpp>
#include <S4r/StarProduct.hpp>
#include <S4r/Debug.hpp>
#include <iostream>
#include <float.h>

S4r::Simulation::Simulation(const double Lr[4], const size_t n[2], size_t nmodes){
	S4R_TRACE("> Simulation::Simulation(Lr = {%g,%g,%g,%g}, n = {%lu,%lu})\n",
		Lr[0], Lr[1], Lr[2], Lr[3], n[0], n[1]
	);
	
	{ // determine lattice type
		double dot = Lr[0]*Lr[2] + Lr[1]*Lr[3];
		double area = Lr[0]*Lr[3] - Lr[1]*Lr[2];
		if(-dot < 1e-10 * area){
			mesh = new LatticeGridRect(Vec2(Lr[0], Lr[1]), Vec2(Lr[2], Lr[3]), n);
		}else{
			mesh = new LatticeGridArb(Vec2(Lr[0], Lr[1]), Vec2(Lr[2], Lr[3]), n);
		}
	}
	
	omega = 1;
	excitation.k[0] = 0;
	excitation.k[1] = 0;
	excitation.pol[0] = 1.;
	excitation.pol[1] = 0.;
	
	options.verbosity = 0;
	options.num_modes = nmodes;
	
	layer_modes_computed = false;
	S4R_TRACE("< Simulation::Simulation\n");
}

S4r::Simulation::~Simulation(){
	S4R_TRACE("> Simulation::~Simulation()\n");
	for(size_t i = 0; i < layers.size(); ++i){
		delete layers[i];
	}
	for(MaterialList::iterator i = materials.begin(); i != materials.end(); ++i){
		delete *i;
	}
	delete mesh;
	S4R_TRACE("< Simulation::~Simulation\n");
}

int S4r::Simulation::SetMaterial(const std::string &name, Material* matobj){
	S4R_TRACE("> Simulation::SetMaterial(name = `%s', matobj = %p)\n",
		name.c_str(), matobj
	);
	MaterialMap::iterator i = matmap.find(name);
	if(matmap.end() == i){
		matmap[name] = materials.size();
		materials.push_back(matobj);
	}else{
		size_t idx = i->second;
		delete materials[idx];
		materials[idx] = matobj;
	}
	S4R_TRACE("< Simulation::SetMaterial\n");
	return 0;
}

S4r::Layer* S4r::Simulation::AddLayer(const std::string &name, const double &thickness, const std::string &material){
	S4R_TRACE("> Simulation::AddLayer(name = `%s', thickness = %g, material = `%s')\n",
		name.c_str(), thickness, material.c_str()
	);
	layer_name_to_index[name] = layers.size();
	Layer *ptr = new Layer(name);
	layers.push_back(ptr);
	ptr->description.thickness = thickness;
	ptr->description.material = material;
	ptr->description.num_modes = options.num_modes;
	S4R_TRACE("< Simulation::AddLayer\n");
	return ptr;
}

void S4r::Simulation::SetK(const double k[2]){
	S4R_TRACE("> Simulation::SetK(k = {%g,%g})\n",
		k[0], k[1]
	);
	// phi should be phase of half angle
	doublecomplex phi[2] = { // should resolve k into Lk components
		doublecomplex(cos(M_PI*k[0]), sin(M_PI*k[0])),
		doublecomplex(cos(M_PI*k[1]), sin(M_PI*k[1]))
	};
	mesh->BuildMatrices(phi);
	ResetLayerModes();
	ResetLayerSolutions();
	S4R_TRACE("< Simulation::SetK\n");
}
void S4r::Simulation::SetFrequency(const doublecomplex &freq){
	S4R_TRACE("> Simulation::SetFrequency(freq = %g + I %g)\n",
		freq.real(), freq.imag()
	);
	omega = 2*M_PI*freq;
	ResetLayerModes();
	ResetLayerSolutions();
	S4R_TRACE("< Simulation::SetFrequency\n");
}


int S4r::Simulation::SetRegionCircle(const std::string &layer, const std::string &matname, const Vec2 &center, double radius){
	S4R_TRACE("> Simulation::SetRegionCircle(layer = `%s', matname = `%s', center = {%g,%g}, radius = %g)\n",
		layer.c_str(), matname.c_str(), center[0], center[1], radius
	);
	int imat = MaterialNameToIndex(matname);
	int ilayer = LayerNameToIndex(layer);
	if(imat < 0 || ilayer < 0){ return 0; }
	
	Layer &cl = *(layers[ilayer]);
	Pattern &pat = cl.description.pattern;
	
	pat.AddShape(new ShapeCircle(center, radius, imat));
	
	ResetLayerSolutions();
	cl.modes.valid = false;
	cl.matrices.valid = false;
	
	S4R_TRACE("< Simulation::SetRegionCircle\n");
	return 0;
}

int S4r::Simulation::SetRegionEllipse(const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const Vec2 &halfwidths){
	S4R_TRACE("> Simulation::SetRegionEllipse(layer = `%s', matname = `%s', center = {%g,%g}, angle = %g, halfwidths = {%g,%g})\n",
		layer.c_str(), matname.c_str(), center[0], center[1], angle, halfwidths[0], halfwidths[1]
	);
	int imat = MaterialNameToIndex(matname);
	int ilayer = LayerNameToIndex(layer);
	if(imat < 0 || ilayer < 0){ return 0; }
	
	Layer &cl = *(layers[ilayer]);
	Pattern &pat = cl.description.pattern;
	
	pat.AddShape(new ShapeEllipse(center, angle, halfwidths, imat));
	
	ResetLayerSolutions();
	cl.modes.valid = false;
	cl.matrices.valid = false;
	
	S4R_TRACE("< Simulation::SetRegionEllipse\n");
	return 0;
}
int S4r::Simulation::SetRegionRectangle(const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const Vec2 &halfwidths){
	S4R_TRACE("> Simulation::SetRegionRectangle(layer = `%s', matname = `%s', center = {%g,%g}, angle = %g, halfwidths = {%g,%g})\n",
		layer.c_str(), matname.c_str(), center[0], center[1], angle, halfwidths[0], halfwidths[1]
	);
	int imat = MaterialNameToIndex(matname);
	int ilayer = LayerNameToIndex(layer);
	if(imat < 0 || ilayer < 0){ return 0; }
	
	Layer &cl = *(layers[ilayer]);
	Pattern &pat = cl.description.pattern;
	
	pat.AddShape(new ShapeRectangle(center, angle, halfwidths, imat));
	
	ResetLayerSolutions();
	cl.modes.valid = false;
	cl.matrices.valid = false;
	
	S4R_TRACE("< Simulation::SetRegionRectangle\n");
	return 0;
}
int S4r::Simulation::SetRegionPolygon  (const std::string &layer, const std::string &matname, const Vec2 &center, double angle, const std::vector<Vec2> &vert){
	S4R_TRACE("> Simulation::SetRegionPolygon(layer = `%s', matname = `%s', center = {%g,%g}, angle = %g, vert=...)\n",
		layer.c_str(), matname.c_str(), center[0], center[1], angle
	);
	int imat = MaterialNameToIndex(matname);
	int ilayer = LayerNameToIndex(layer);
	if(imat < 0 || ilayer < 0){ return 0; }
	
	Layer &cl = *(layers[ilayer]);
	Pattern &pat = cl.description.pattern;
	
	pat.AddShape(new ShapePolygon(center, angle, vert, imat));
	
	ResetLayerSolutions();
	cl.modes.valid = false;
	cl.matrices.valid = false;
	
	S4R_TRACE("< Simulation::SetRegionPolygon\n");
	return 0;
}
int S4r::Simulation::RemoveLayerRegions(const std::string &layer){
	S4R_TRACE("> Simulation::RemoveLayerRegions(layer = `%s')\n",
		layer.c_str()
	);
	int ilayer = LayerNameToIndex(layer);
	if(ilayer < 0){ return 0; }
	
	Layer &cl = *(layers[ilayer]);
	Pattern &pat = cl.description.pattern;
	
	pat.RemoveShapes();
	
	ResetLayerSolutions();
	cl.modes.valid = false;
	cl.matrices.valid = false;
	
	S4R_TRACE("< Simulation::RemoveLayerRegions\n");
	return 0;
}

int S4r::Simulation::SetLayerThickness(const std::string &layer, const double &thickness){
	S4R_TRACE("> Simulation::RemoveLayerRegions(layer = `%s', thickness = %g)\n",
		layer.c_str(), thickness
	);
	LayerMap::iterator i = layer_name_to_index.find(layer);
	if(layer_name_to_index.end() == i){
		// error
		return -1;
	}
	Layer &cl = *layers[i->second];
	cl.description.thickness = thickness;
	
	ResetLayerSolutions();
	S4R_TRACE("< Simulation::SetLayerThickness\n");
	return 0;
}

int S4r::Simulation::GetSMatrix(
	size_t layer_start, size_t layer_end,
	CMat &S
){
	S4R_TRACE("> Simulation::GetSMatrix(start = %lu, end = %lu)\n",
		layer_start, layer_end
	);

	const size_t m = mesh->NumEdges();
	
	const size_t n0 = layers[layer_start]->GetNumModes();
	
	// Ultimately, S will have dimensions:
	//
	//           n0   n1
	//     n1 [ S11  S12 ]
	// S =    [          ]
	//     n0 [ S21  S22 ]
	//
	
	S.resize(n0+n0,n0+n0); // S starts off at layer_start
	S.setIdentity();
	for(size_t l = layer_start; l < layer_end; ++l){
		const size_t lp1 = l+1;
		const size_t nl = layers[l]->GetNumModes();
		const size_t nlp1 = layers[lp1]->GetNumModes();
		
		const Layer &Ll = (*layers[l]);
		const Layer &Llp1 = (*layers[lp1]);
		
		CMat in1(nl,m);
		CMat in2(nlp1,m);
		
		// Make the interface matrices. They will have dimensions:
		//          nlp1 nlp1
		//     nl [ I11  I12 ]
		// I =    [          ]
		//     nl [ I21  I22 ]
		//
		if(false){
			// This is a trivial interface, set to identity
			//in1.setIdentity();
			//in2.setIdentity();
		}else{
			// The interface matrix is the inverse of the mode-to-field matrix of layer l
			// times the mode-to-field matrix of layer l+1 (lp1).
			// The mode-to-field matrix is of the form
			//   [ B -B ] where A = phi
			//   [ A  A ] where B = kp*phi*inv(diag(q)) = G*A/q
			// So we want
			//   Interface = 0.5 * [  iBl  iAl ] [ Blp1 -Blp1 ]
			//                     [ -iBl  iAl ] [ Alp1  Alp1 ]
			// where iBl = inv(Bl), etc.
			// Multiplying out gives
			//   0.5 * [ P+Q P-Q ] // where P = iAl*Alp1, and i in front means inverse
			//         [ P-Q P+Q ] // where Q = iBl*Blp1
			// Making P is easy, since A is a single matrix.
			// Making Q is as follows:
			//   Q = iBl*Blp1
			//     = ql*iAl*iGl * Gl*Alp1*iqlp1
			// We will only store I11 and I21

			{
				// Make Blp1 in in2
				CMat t1(Llp1.modes.kp * Llp1.modes.phi);
				
				// Make Q in in1
				// Take care when inverting t1, since it may be singular at a
				// diffraction threshold.
				Eigen::ColPivHouseholderQR<CMat> iBl(Ll.modes.kp * Ll.modes.phi);
				in1 = Ll.modes.q.asDiagonal() * iBl.solve(t1);
				{ // Take special care since an element of q may be zero
					double maxel = 0;
					for(size_t i = 0; i < nlp1; ++i){
						double el = std::abs(Llp1.modes.q[i]);
						if(el > maxel){ 
							maxel = el;
						}
					}
					for(size_t i = 0; i < nlp1; ++i){
						double el = std::abs(Llp1.modes.q[i]);
						if(el > maxel * DBL_EPSILON){
							in1.col(i) *= 1./Llp1.modes.q[i];
						}else{
							in1.col(i).setZero();
						}
					}
				}
			}
			
			// Make P in in2
			{
				Eigen::ColPivHouseholderQR<CMat> iAl(Ll.modes.phi);
				// t1 may become singular for planewaves at diffraction thresholds
				// with a z-directed polarization.
				in2 = iAl.solve(Llp1.modes.phi);
			}
			
			CMat t1(in2); // in2 = P, t1 = P, in1 = Q
			in2 -= in1; // in2 = P-Q, t1 = P, in1 = Q
			in1 += t1; // in2 = P-Q, t1 = P, in1 = P+Q
			in1 *= 0.5;
			in2 *= 0.5;
			
			// The inverse of the Interface matrix is
			//   inv(Interface) = 0.5 * [ iP+iQ    iP-iQ ]
			//                          [ iP-iQ    iP+iQ ]
			// where iP = inv(P), etc.
		}
		
		CVec d1(nl);
		for(size_t i = 0; i < nl; ++i){
			d1[i] = std::exp(Ll.modes.q[i] * std::complex<double>(0,Ll.description.thickness));
		}
		CVec d2(nlp1);
		for(size_t i = 0; i < nlp1; ++i){
			d2[i] = std::exp(Llp1.modes.q[i] * std::complex<double>(0,Llp1.description.thickness));
		}
		//std::cout << "Layer " << l   << " d1:\n" << d1 << "\n";
		//std::cout << "Layer " << lp1 << " d2:\n" << d2 << "\n";

		// At this point, S is
		//           n0   nl
		//     nl [ S11  S12 ]
		// S =    [          ]
		//     n0 [ S21  S22 ]
		
		CMat *Snew = &S;
		if(nl != nlp1){
			Snew = new CMat(n0+nlp1, n0+nlp1);
		}

		{
			// We need to do two simultaneous solves, so we will form the
			// adjoined matrix of both right hand sides:
			//   J = [ f_l S11    |    (f_l S12 I22 - I12) f_{l+1} ]
			// They have sizes:
			//          n0     nlp1
			//    nl [      |       ]
			// After the solve, we should have the matrices
			//       [ S11  |  S12 ]
			// which has sizes
			//            n0     nlp1
			//    nlp1 [      |       ]
			CMat A(in1 - d1.asDiagonal() * S.block(0,n0, nl,nl) * in2);
			size_t mx = (nl > nlp1 ? nl : nlp1);
		
			CMat J(mx, n0+nlp1);
			J.block(0, 0, nl,n0  ) = d1.asDiagonal() * S.block(0,0, nl,n0);
			J.block(0,n0, nl,nlp1) = (d1.asDiagonal() * S.block(0,n0, nl,nl) * in1 - in2) * d2.asDiagonal();
			
			LeastNormSolve(
				A.rows(), A.cols(), n0+nlp1,
				A.data(), A.outerStride(),
				J.data(), J.outerStride()
			);
			Snew->block(0, 0, nlp1,n0  ) = J.block(0, 0, nlp1,n0  );
			Snew->block(0,n0, nlp1,nlp1) = J.block(0,n0, nlp1,nlp1);
		}
		
		{
			// Make S'21 = S21 + S22 I22 S'11
			CMat t1(S.block(nl,n0, n0,nl) * in2); // t1 = S22 I21
			Snew->block(nlp1,0, n0,n0) += t1 * Snew->block(0,0, nlp1,n0);
			
			// Make S'22 = S22 I22 f_{l+1} + t1 S'12
			Snew->block(nlp1,n0, n0,nlp1) = S.block(nl,n0, n0,nl) * in1 * d2.asDiagonal() + t1 * Snew->block(0,n0, nlp1,nlp1);
		}
		if(nl != nlp1){
			S = *Snew;
			delete Snew;
		}
	}
	
	S4R_TRACE("< Simulation::GetSMatrix\n");
	return 0;
}
void S4r::Simulation::SolveLayer(size_t layer){
	S4R_TRACE("> Simulation::SolveLayer(layer = %lu)\n",
		layer
	);
	Layer &Ll = (*layers[layer]);
	if((Ll.solution.flags & Layer::Solution::A_SET) &&
	   (Ll.solution.flags & Layer::Solution::B_SET)
	){
		return;
	}
	ComputeLayerModes();
	
	size_t i0 = layer;
	size_t iN = layer;
	if(!(layers.front()->solution.flags & Layer::Solution::A_SET) ||
	   !(layers.back( )->solution.flags & Layer::Solution::B_SET)
	){
		BuildExcitation();
		i0 = 0;
		iN = layers.size()-1;
	}else{
		// search outwards for a layer with the required solutions
		while(!(layers[i0]->solution.flags & Layer::Solution::A_SET) && i0 > 0){
			--i0;
		}
		while(!(layers[iN]->solution.flags & Layer::Solution::A_SET) && iN < layers.size()-1){
			++iN;
		}
	}
	
	const size_t n0 = layers[i0]->GetNumModes();
	const size_t nN = layers[iN]->GetNumModes();
	const size_t nl = layers[layer]->GetNumModes();
	
	const CVec &a0 = layers[i0]->solution.a;
	const CVec &bN = layers[iN]->solution.b;
	
	CMat S0l(n0+nl,n0+nl); GetSMatrix(i0, layer, S0l);
	CMat SlN(nl+nN,nl+nN); GetSMatrix(layer, iN, SlN);
	
	//             n0   nl                   nl   nN
	//       nl [ S11  S12 ]           nN [ S11  S12 ]
	// S0l =    [          ]     SlN =    [          ]
	//       n0 [ S21  S22 ]           nl [ S21  S22 ]

	CVec S11a0 = S0l.block(0,0, nl,n0) * a0;
	CVec S22bN = S0l.block(nl,n0, n0,nl) * bN;
	
	// We overwrite the upper left submatrix S11(0,l) since it's not needed anymore
	// temp is set to S0l for this reason.
	
	{
		CMat temp(nl,nl);
		temp.setIdentity();
		temp -= S0l.block(0,n0, nl,nl) * SlN.block(nN,0, nl,nl); // temp = (1 - S_12(0,l)S_21(l,N))
	
//std::cout << "itemp = " << S0l.block(0,0, n,n).inverse() << "\n";
	
		Ll.solution.a = temp.partialPivLu().solve(
			S11a0 + S0l.block(0,n0, nl,nl) * S22bN
		); // al = inv(temp) * [ S_11(0,l)*a0 + S_12(0,l)S_22(l,N)bN ]
		Ll.solution.flags |= Layer::Solution::A_SET;
	}
	
	// Make the other matrix
	{
		// Compute S_21(l,N)S_12(0,l)
		CMat temp(nl, nl);
		temp.setIdentity();
		temp -= SlN.block(nN,0, nl,nl) * S0l.block(0,n0, nl,nl); // temp = (1 - S_21(l,N)S_12(0,l))
		
//std::cout << "itemp = " << S0l.block(0,0, n,n).inverse() << "\n";
	
		Ll.solution.b = temp.partialPivLu().solve(
			S22bN + SlN.block(nN,0, nl,nl) * S11a0
		); // bl = inv(temp) * [ S_21(l,N)S_11(0,l)a0 + S_22(l,N)bN ]
		Ll.solution.flags |= Layer::Solution::B_SET;
	}
	
	// Now check the outer layers to see if they require solutions,
	// since we can easily obtain those.
	Layer &L0 = (*layers[i0]);
	Layer &LN = (*layers[iN]);
	if(!(L0.solution.flags & Layer::Solution::B_SET)){
		// b0 = S_21(0,l)*a0 + S_22(0,l)*bl
		L0.solution.b = S0l.block(nl,0, n0,n0) * a0 + S0l.block(nl,n0, n0,nl) * Ll.solution.b;
		L0.solution.flags |= Layer::Solution::B_SET;
	}
	if(!(LN.solution.flags & Layer::Solution::A_SET)){
		LN.solution.a = SlN.block(0,0, nN,nl) * Ll.solution.a + SlN.block(0,nl, nN,nN) * bN;
		LN.solution.flags |= Layer::Solution::A_SET;
	}

	S4R_TRACE("< Simulation::SolveLayer\n");
}

int S4r::Simulation::MaterialNameToIndex(const std::string &matname) const{
	MaterialMap::const_iterator i = matmap.find(matname);
	if(i == matmap.end()){
		return -1;
	}
	return (int)(i->second);
}
int S4r::Simulation::LayerNameToIndex(const std::string &name) const{
	LayerMap::const_iterator i = layer_name_to_index.find(name);
	if(i == layer_name_to_index.end()){
		return -1;
	}
	return (int)(i->second);
}

void S4r::Simulation::BuildExcitation(){
	S4R_TRACE("> Simulation::BuildExcitation()\n");
	layers.front()->ComputeModes(omega, mesh);
	layers.back()->ComputeModes(omega, mesh);
	
	CVec &a = layers.front()->solution.a;
	CVec &b = layers.back()->solution.b;
	
	const size_t m = mesh->NumEdges();
	a.resize(m);
	// build some kind of planewave
	
	for(size_t i = 0; i < m; ++i){
		const Vec2 edge = mesh->Edge(i);
		//std::cout << "i = " << i << ", edge = " << edge << std::endl;
		a[i] = excitation.pol[0] * edge[0] + excitation.pol[1] * edge[1];
	}
//std::cout << "BuildEx: a = \n" << a << std::endl;
//std::cout << "phi = \n" << layers.front()->modes.phi << std::endl;
	{
		CMat phicopy(layers.front()->modes.phi);
		LeastNormSolve(
			phicopy.rows(), phicopy.cols(), 1,
			phicopy.data(), phicopy.outerStride(),
			a.data(), a.size()
		);
		a.conservativeResize(phicopy.cols());
	}
//std::cout << "BuildEx: a = \n" << a << std::endl;
//std::cout << "phi = \n" << layers.front()->modes.phi << std::endl;
//std::cout << "pol = \n" << excitation.pol[0] << ", " << excitation.pol[1] << std::endl;
	
	//a.setZero(); a[1] = 1.;
	b.resize(layers.back()->modes.n);
	b.setZero();
	layers.front()->solution.flags |= Layer::Solution::A_SET;
	layers.back()->solution.flags |= Layer::Solution::B_SET;
	
	S4R_TRACE("< Simulation::BuildExcitation\n");
}
void S4r::Simulation::SolveAllLayers(size_t i0, size_t iN){
	S4R_TRACE("> Simulation::SolveAllLayers(i0 = %lu, iN = %lu)\n", i0, iN);
	if(i0+1 <= iN){ return; }
	size_t imid = (i0+iN)/2;
	SolveLayer(imid);
	SolveAllLayers(i0, imid);
	SolveAllLayers(imid, iN);
	S4R_TRACE("< Simulation::SolveAllLayers\n");
}
void S4r::Simulation::SolveAllLayers(){
	// recursively subdivide
	SolveAllLayers(0, layers.size()-1);
}

void S4r::Simulation::ComputeLayerModes(){
	S4R_TRACE("> Simulation::ComputeLayerModes()\n");
	//if(layer_modes_computed){ return; }
	for(size_t i = 0; i < layers.size(); ++i){
		layers[i]->BuildMatrices(mesh, matmap, materials);
		layers[i]->ComputeModes(omega, mesh);
	}
	//layer_modes_computed = true;
	S4R_TRACE("< Simulation::ComputeLayerModes\n");
}
void S4r::Simulation::ResetLayerModes(){
	S4R_TRACE("> Simulation::ResetLayerModes()\n");
	for(size_t i = 0; i < layers.size(); ++i){
		layers[i]->modes.valid = false;
	}
	S4R_TRACE("< Simulation::ResetLayerModes\n");
}
void S4r::Simulation::ResetLayerSolutions(){
	S4R_TRACE("> Simulation::ResetLayerSolutions()\n");
	for(size_t i = 0; i < layers.size(); ++i){
		layers[i]->solution.flags &= ~(Layer::Solution::A_SET | Layer::Solution::B_SET);
	}
	S4R_TRACE("< Simulation::ResetLayerSolutions\n");
}
void S4r::Simulation::GetPowerFlux(const std::string &layer, const double &z, doublecomplex &forw, doublecomplex &back){
	S4R_TRACE("> Simulation::GetPowerFlux(layer = `%s', z = %g)\n",
		layer.c_str(), z
	);
	int ilayer = LayerNameToIndex(layer);
	SolveLayer(ilayer);
	layers[ilayer]->GetPowerFlux(omega, mesh, z, forw, back);
	S4R_TRACE("< Simulation::GetPowerFlux\n");
}


void S4r::Simulation::OutputLayerDescription(const std::string &layer, FILE *fp) const{
	int ilayer = LayerNameToIndex(layer);
	layers[ilayer]->BuildMatrices(mesh, matmap, materials);
	layers[ilayer]->DumpDescription(mesh, fp);
}

S4r::CVec2 S4r::Simulation::GetBlochPhaseFactors() const{
	const double angle[2] = {
		excitation.k[0] * 2.*M_PI,
		excitation.k[1] * 2.*M_PI,
	};
	return CVec2(
		doublecomplex(cos(angle[0]), sin(angle[0])),
		doublecomplex(cos(angle[1]), sin(angle[1]))
	);
}


void S4r::Simulation::GetFields(const Vec3 &r, CVec3 &E, CVec3 &H){
	if(0 == layers.size()){
		// no layers
		E.setZero();
		H.setZero();
	}
	
	size_t ilayer = 0;
	double dz = r[2];
	{
		double z = 0;
		double t = layers[ilayer]->description.thickness;
		while(r[2] > z + t){
			z += t;
			dz -= t;
			if(ilayer+1 >= layers.size()){ break; }
			ilayer++;
			t = layers[ilayer]->description.thickness;
		}
	}
//fprintf(stderr, "(%f,%f,%f) in %s: dz = %f\n", r[0], r[1], r[2], L->name, dz);
	
	SolveLayer(ilayer);
	CVec2 phi = GetBlochPhaseFactors();
	layers[ilayer]->GetFields(omega, mesh, phi, Vec3(r[0], r[1], dz), E, H);
}
