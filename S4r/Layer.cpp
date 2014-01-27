#include <S4r/Layer.hpp>
#include <S4r/Eigensystems.hpp>
#include <S4r/PeriodicMesh.hpp>
#include <S4r/Material.hpp>
#include <S4r/Debug.hpp>

#include <iostream>

using namespace S4r;

S4r::Layer::Layer(const std::string &name){
	description.name = name;
	description.thickness = 0;
	description.num_modes = 0;
	matrices.valid = false;
	modes.valid = false;
	solution.flags = 0;
}
S4r::Layer::~Layer(){
}

static void RotateTensor(const Vec2 &n, CTensor2 &M){
	// Computes R^T M R where R = [ n[0] -n[1] ]
	//                            [ n[1]  n[0] ]
	doublecomplex MR[4] = { // row major
		M(0,0)*n[0] + M(0,1)*n[1], M(0,1)*n[0] - M(0,0)*n[1],
		M(1,0)*n[0] + M(1,1)*n[1], M(1,1)*n[0] - M(1,0)*n[1]
	};
	M(0,0) = n[0]*MR[0] + n[1]*MR[2];
	M(1,0) = n[0]*MR[2] - n[1]*MR[0];
	M(0,1) = n[0]*MR[1] + n[1]*MR[3];
	M(1,1) = n[0]*MR[3] - n[1]*MR[1];
}

static void TauTransform(CTensor2 &A){
	// tau(e__) = [ -1/e11         e12/e11      ]
	//            [ e21/e11   e22 - e21 e12/e11 ]
	doublecomplex ia11 = 1. / A(0,0);
	A(0,0) = -ia11;
	A(1,0) *= ia11;
	A(1,1) -= A(1,0)*A(0,1);
	A(0,1) *= ia11;
}
static void InverseTauTransform(CTensor2 &A){
	// invtau(t__) = [ -1/t11         -t12/t11      ]
	//               [ -t21/t11   t22 - t21 t12/t11 ]
	doublecomplex ia11 = -1. / A(0,0);
	A(0,0) = ia11;
	A(1,0) *= ia11;
	A(1,1) += A(1,0)*A(0,1);
	A(0,1) *= ia11;
}
static void AnisotropicAverage(
	Vec2 &n,
	const double &f1, CTensor2 &e1,
	const double &f2, CTensor2 &e2,
	CTensor2 &e
){
	RotateTensor(n, e1);
	RotateTensor(n, e2);
//std::cerr << "rot eps1 = " << e1 << "\n";			
//std::cerr << "rot eps2 = " << e2 << "\n";			
	TauTransform(e1);
	TauTransform(e2);
//std::cerr << "tau1 = " << e1 << "\n";
//std::cerr << "tau2 = " << e2 << "\n";
	e = f1 * e1 + f2 * e2;
//std::cerr << "tauavg = " << e << "\n";
	InverseTauTransform(e);
//std::cerr << "epsavg = " << eps << "\n";
	// Undo rotation
	n[1] = -n[1];
	RotateTensor(n, e);
}
void S4r::Layer::GetMaterialAverage(
	const ConvexPolygon &poly,
	CTensor2 &eps, CTensor2 &mu,
	const std::vector<Material*> &mat,
	const Material *bkmat,
	std::vector<double> &value // workspace
){
	size_t nmat = description.pattern.Overlap(poly, value);
	if(2 != nmat){
		eps = value[0] * bkmat->eps.block(0,0,2,2);
		for(size_t s = 0; s < description.pattern.NumShapes(); ++s){
			eps += value[s+1] * mat[description.pattern.GetShape(s).tag]->eps.block(0,0,2,2);
		}
		mu = value[0] * bkmat->mu.block(0,0,2,2);
		for(size_t s = 0; s < description.pattern.NumShapes(); ++s){
			mu += value[s+1] * mat[description.pattern.GetShape(s).tag]->mu.block(0,0,2,2);
		}
	}else{
		// perform fancy averaging
		int ip1 = -1, ip2 = -1;
		for(size_t i = 0; i <= description.pattern.NumShapes(); ++i){
			if(value[i] > 0){
				if(ip1 >= 0){ ip2 = i; break; }
				else{ ip1 = i; }
			}
		}
		// Assert that ip1 and ip2 >= 0
		CTensor2 eps1, mu1;
		CTensor2 eps2, mu2;
		const double fill1 = value[ip1], fill2 = value[ip2];
		if(0 == ip1){
			eps1 = bkmat->eps.block(0,0,2,2);
			mu1  = bkmat->mu.block(0,0,2,2);
		}else{
			int imat = description.pattern.GetShape(ip1-1).tag;
			eps1 = mat[imat]->eps.block(0,0,2,2);
			mu1  = mat[imat]->mu.block(0,0,2,2);
		}
		{
			int imat = description.pattern.GetShape(ip2-1).tag;
			eps2 = mat[imat]->eps.block(0,0,2,2);
			mu2  = mat[imat]->mu.block(0,0,2,2);
		}
		Vec2 n(description.pattern.GetShape(ip2-1).Normal(poly.ApproxCenter()));

		AnisotropicAverage(n, fill1, eps1, fill2, eps2, eps);
		AnisotropicAverage(n, fill1,  mu1, fill2,  mu2, mu );
	}
}

int S4r::Layer::BuildMatrices(
	const PeriodicMesh *mesh,
	const std::map<std::string, size_t> &matmap,
	const std::vector<Material*> &mat
){
	S4R_TRACE("> Layer::BuildMatrices(mesh = %p, ...) name = %s\n", mesh, description.name.c_str());
	if(matrices.valid){
		S4R_TRACE("< Layer::BuildMatrices (early exit)\n");
		return 0;
	}
	const size_t Nv = mesh->NumVertices();
	const size_t Nf = mesh->NumFaces();
	const size_t Ne = mesh->NumEdges();
	
	matrices.ieps.resize(Nf, Nf);
	matrices.ieps.reserve(Nf);
	matrices.imu.resize(Nv, Nv);
	matrices.imu.reserve(Nv);
	matrices.eps.resize(Ne, Ne);
	matrices.mu.resize(Ne, Ne);
	
	std::vector<Triplet> triplst;
	
	// Determine background material
	const Material * bkmat;
	{
		std::map<std::string, size_t>::const_iterator i = matmap.find(description.material);
		if(matmap.end() == i){
			// error
			return -1;
		}
		bkmat = mat[i->second];
	}
	
	matrices.flags = bkmat->flags;
	bool all_scalar = true;
	for(size_t s = 0; s < description.pattern.NumShapes(); ++s){
		int imat = description.pattern.GetShape(s).tag;
		matrices.flags &= mat[imat]->flags;
		if(!(mat[imat]->flags & Material::SCALAR)){ all_scalar = false; }
	}
	
	bool apply_scaling = !IsUniformHermitian();
	if(!apply_scaling){
		matrices.eps_unif = bkmat->eps(0,0).real();
		matrices.mu_unif = bkmat->mu(0,0).real();
	}else{
		// fill in with something?
	}
	
	if(apply_scaling){
		description.pattern.Finalize();
		std::vector<double> value;
		
		// Make ieps
		triplst.clear();
		triplst.reserve(Ne);
		for(size_t iface = 0; iface < Nf; ++iface){
			ConvexPolygon poly;
			mesh->FaceControlPolygon(iface, &poly);
			description.pattern.Overlap(poly, value);
			doublecomplex eps(value[0] * bkmat->eps(2,2));
			for(size_t s = 0; s < description.pattern.NumShapes(); ++s){
				eps += value[s+1] * mat[description.pattern.GetShape(s).tag]->eps(2,2);
			}
	//std::cerr << poly.offset[0] << "\t" << poly.offset[1] << "\t" << ieps.real() << "\t" << ieps.imag() << "\n";
			triplst.push_back(Triplet(
				iface, iface, 1. / (eps * mesh->FaceArea(iface))
			));
		}
		matrices.ieps.setFromTriplets(triplst.begin(), triplst.end());
		
		// Make imu
		triplst.clear();
		triplst.reserve(Nv);
		for(size_t ivert = 0; ivert < Nv; ++ivert){
			ConvexPolygon poly;
			mesh->VertexControlPolygon(ivert, &poly);
			description.pattern.Overlap(poly, value);
			doublecomplex mu(value[0] * bkmat->mu(2,2));
			for(size_t s = 0; s < description.pattern.NumShapes(); ++s){
				mu += value[s+1] * mat[description.pattern.GetShape(s).tag]->mu(2,2);
			}
	//std::cerr << poly.offset[0] << "\t" << poly.offset[1] << "\t" << imu.real() << "\t" << imu.imag() << "\n";
			triplst.push_back(Triplet(
				ivert, ivert, 1. / (mu * mesh->VertexDualArea(ivert))
			));
		}
		matrices.imu.setFromTriplets(triplst.begin(), triplst.end());
		
		// Make eps and mu
		triplst.clear();
		triplst.reserve(Ne); // modify for offdiagonals
		std::vector<Triplet> triplstmu;
		triplstmu.reserve(Ne);
		for(size_t iedge = 0; iedge < Ne; ++iedge){
			ConvexPolygon poly[3];
			size_t eother[3];
			mesh->EdgeControlPolygon(iedge, &poly[0], &eother[1], &poly[1], &eother[2], &poly[2]);
			eother[0] = iedge;
			
			// Grab the right component
			const Vec2 edgedir(mesh->Edge(iedge).normalized());

			unsigned nother = all_scalar ? 1 : 3;
			for(unsigned iother = 0; iother < nother; ++iother){
				const size_t iedge2 = eother[iother];
				const Vec2 otherdir(mesh->Edge(iedge2).normalized());
				CTensor2 eps, mu;
				GetMaterialAverage(poly[iother], eps, mu, mat, bkmat, value);
				doublecomplex epsc(otherdir.transpose() * eps * edgedir);
				doublecomplex muc (otherdir.transpose() * mu  * edgedir);
				triplst.push_back(Triplet(
					iedge2, iedge, epsc * mesh->EdgeLength(iedge2) / mesh->EdgeDualLength(iedge)
				));
				triplstmu.push_back(Triplet(
					iedge2, iedge, muc * mesh->EdgeDualLength(iedge2) / mesh->EdgeLength(iedge)
				));
				if(iedge2 != iedge){ // 2nd offdiag contribution
					epsc = (edgedir.transpose() * eps * otherdir);
					muc  = (edgedir.transpose() * mu  * otherdir);
					triplst.push_back(Triplet(
						iedge, iedge2, epsc * mesh->EdgeLength(iedge) / mesh->EdgeDualLength(iedge2)
					));
					triplstmu.push_back(Triplet(
						iedge, iedge2, muc * mesh->EdgeDualLength(iedge) / mesh->EdgeLength(iedge2)
					));
				}
			}
		}
		matrices.eps.setFromTriplets(triplst.begin(), triplst.end());
		matrices.mu.setFromTriplets(triplstmu.begin(), triplstmu.end());
		//CMat temp(matrices.eps); temp -= matrices.eps.adjoint();
		//std::cout << "diff = " << temp.norm() << std::endl;
	}else{
		// Make ieps
		triplst.clear();
		triplst.reserve(Ne);
		for(size_t iface = 0; iface < Nf; ++iface){
			triplst.push_back(Triplet(
				iface, iface, 1. / mesh->FaceArea(iface)
			));
		}
		matrices.ieps.setFromTriplets(triplst.begin(), triplst.end());
		
		// Make imu
		triplst.clear();
		triplst.reserve(Nv);
		for(size_t ivert = 0; ivert < Nv; ++ivert){
			triplst.push_back(Triplet(
				ivert, ivert, 1. / mesh->VertexDualArea(ivert)
			));
		}
		matrices.imu.setFromTriplets(triplst.begin(), triplst.end());
		
		// Make eps and mu
		triplst.clear();
		triplst.reserve(Ne); // modify for offdiagonals
		std::vector<Triplet> triplstmu;
		triplstmu.reserve(Ne);
		for(size_t iedge = 0; iedge < Ne; ++iedge){
			double ratio = mesh->EdgeLength(iedge) / mesh->EdgeDualLength(iedge);
			triplst.push_back(Triplet(
				iedge, iedge, ratio
			));
			triplstmu.push_back(Triplet(
				iedge, iedge, 1. / ratio
			));
		}
		matrices.eps.setFromTriplets(triplst.begin(), triplst.end());
		matrices.mu.setFromTriplets(triplstmu.begin(), triplstmu.end());
	}
	matrices.valid = true;
	S4R_TRACE("< Layer::BuildMatrices\n");
	return 0;
}

bool S4r::Layer::IsUniformHermitian() const{
	return (0 == description.pattern.NumShapes() &&
	   (matrices.flags & Matrices::SCALAR) &&
	   (matrices.flags & Matrices::HERMITIAN));
}

size_t S4r::Layer::GetNumModes() const{
	return modes.n;
}

int S4r::Layer::ComputeModes(
	const doublecomplex &omega,
	const PeriodicMesh *mesh
){
	S4R_TRACE("> Layer::ComputeModes(omega = %g + I %g, mesh = %p) name = %s\n", omega.real(), omega.imag(), mesh, description.name.c_str());
	if(modes.valid){
		S4R_TRACE("< Layer::ComputeModes (early exit)\n");
		return 0;
	}
	const doublecomplex w2 = omega*omega;
	
	const size_t m = matrices.eps.rows();
	
	size_t n = description.num_modes;
	//if(0 == n){ n = this->options.num_modes; }
	if(0 == n || n > mesh->NumEdges()){ n = mesh->NumEdges(); }
	modes.n = n;
	modes.q.resize(n);
	
	//std::cout << "Layer `" << description.name << "' has " << n << " modes" << std::endl;
	
	if(IsUniformHermitian() && omega.imag() == 0.){
		SpMat A(
			(w2 * matrices.eps_unif * matrices.mu_unif) * matrices.mu
			- matrices.mu * mesh->d0 * matrices.imu * mesh->d0.adjoint() * matrices.mu
			- mesh->d1.adjoint() * matrices.ieps * mesh->d1
		);
		modes.kp = (w2*matrices.mu_unif) * matrices.mu - (1. / matrices.eps_unif) * mesh->d1.adjoint() * matrices.ieps * mesh->d1;
		
		//CMat test = modes.phi - modes.phi.adjoint();
		//std::cout << "norm = " << test.norm() << std::endl;
		//std::cout << "eps_unif = " << matrices.eps_unif << std::endl;
		//std::cout << "mu_unif = " << matrices.mu_unif << std::endl;
		//std::cout << "ieps = " << matrices.ieps << std::endl;
		//std::cout << "imu = " << matrices.imu << std::endl;
		//std::cout << "mu = " << matrices.mu << std::endl;
		//std::cout << "A = " << modes.phi << std::endl;
		
		int ret = HermitianEigensystem(n, A, modes.q, modes.phi);
		if(0 != ret){
			std::cerr << "Layer `" << description.name << "': HermitianEigensystem returned error code: " << ret << std::endl;
		}
		modes.phi = matrices.eps * modes.phi;
		
		/*
		modes.kp = w2*matrices.mu - mesh->d1.adjoint() * matrices.ieps * mesh->d1;
		modes.phi.resize(n,n);
		CMat A(
			w2 * matrices.eps * matrices.mu
			- mesh->d0 * matrices.imu * mesh->d0.adjoint() * matrices.mu
			- matrices.eps * mesh->d1.adjoint() * matrices.ieps * mesh->d1
		);
		//CMat Asave(A);
		Eigensystem(
			m, n, A.data(), A.outerStride(),
			modes.q.data(),
			NULL, 1, modes.phi.data(), modes.phi.outerStride()
		);*/
	}else{
		modes.kp = w2*matrices.mu - mesh->d1.adjoint() * matrices.ieps * mesh->d1;
		modes.phi.resize(m,n);
		SpMat A(
			w2 * matrices.eps * matrices.mu
			- mesh->d0 * matrices.imu * mesh->d0.adjoint() * matrices.mu
			- matrices.eps * mesh->d1.adjoint() * matrices.ieps * mesh->d1
		);
		//CMat Asave(A);
		int ret = Eigensystem(n, A, modes.q, modes.phi);
		if(0 != ret){
			std::cerr << "Layer `" << description.name << "': Eigensystem returned error code: " << ret << std::endl;
		}
		//CMat temp = Asave * modes.phi - modes.phi * modes.q.asDiagonal();
		//std::cout << "# should be zero: " << temp.norm() << "\n";
	}
	
	// Sorting now happens in the Eigensystem functions
	/*
	// Sort eigenvalues with selection sort to minimize number of swaps
	for(size_t j = 0; j+1 < n; ++j){
		size_t imax = j;
		for(size_t i = j+1; i < n; ++i){
			if(modes.q[i].real() > modes.q[imax].real()){
				imax = i;
			}
		}
		if(imax != j){
			std::swap(modes.q[j], modes.q[imax]);
			modes.phi.col(j).swap(modes.phi.col(imax));
		}
	}
	*/
	
	//std::cout << "Layer: " << description.name << "\n";
	//std::cout << "q:\n" << modes.q << "\n";
	//std::cout << "phi:\n" << modes.phi << "\n";
	
	//std::cerr << "q[0]: " << modes.q[0] << "\n";
	
	// Take the right branch cut
	for(size_t i = 0; i < n; ++i){
		// Set the \hat{q} vector (diagonal matrix) while we're at it
		if(0 == omega.imag()){ // Not bandsolving
			modes.q[i] = std::sqrt(modes.q[i]);
			if(modes.q[i].imag() < -1e-10*fabs(modes.q[i].real())){
				modes.q[i] = -modes.q[i];
			}
		}else{ // performing some kind of bandsolving, need to choose the appropriate branch
			if(modes.q[i].real() < 0){
				// branch cut should be just below positive real axis
				modes.q[i] = doublecomplex(0,1) * std::sqrt(-modes.q[i]);
			}else{
				// branch cut should be just below negative real axis
				// This is the default behavior for sqrt(std::complex)
				modes.q[i] = std::sqrt(modes.q[i]);
			}
		}
	}
	
	// The columns come out normalized in the max-element norm
	/*
	for(size_t i = 0; i < n; ++i){
		modes.phi.col(i).normalize();
	}
	*/
	modes.valid = true;
	S4R_TRACE("< Layer::ComputeModes\n");
	return 0;
}

void S4r::Layer::GetTranslatedSolution(
	const double &z, CVec &a, CVec &b
) const{
	const size_t n = solution.a.size();
	a = solution.a;
	b = solution.b;
	const double mz = description.thickness - z;
	for(size_t i = 0; i < n; ++i){
		doublecomplex iq = doublecomplex(0.,1.) * modes.q[i];
		a[i] *= std::exp(iq * z);
		b[i] *= std::exp(iq * mz);
	}
}

void S4r::Layer::GetPowerFlux(
	const doublecomplex &omega,
	const PeriodicMesh *mesh,
	const double &z,
	doublecomplex &forward, doublecomplex &backward
) const{
	S4R_TRACE("> Layer::GetPowerFlux(omega = %g + I %g, mesh = %p, z = %g) name = %s\n", omega.real(), omega.imag(), mesh, z, description.name.c_str());
	//std::cout << "  flags = " << solution.flags << "\n";
	//std::cout << "  a = " << solution.a << "\n";
	//std::cout << "  b = " << solution.b << "\n";
	CVec ea, eb, ha, hb;
	if(0 == description.thickness && 0 == z){
		ea = modes.kp * modes.phi * modes.q.cwiseInverse().asDiagonal() * solution.a;
		ha = modes.phi * solution.a;
		eb = modes.kp * modes.phi * modes.q.cwiseInverse().asDiagonal() * -solution.b;
		hb = modes.phi * solution.b;
	}else{
		CVec a, b;
		GetTranslatedSolution(z, a, b);
		ea = modes.kp * modes.phi * modes.q.cwiseInverse().asDiagonal() * a;
		ha = modes.phi * a;
		eb = modes.kp * modes.phi * modes.q.cwiseInverse().asDiagonal() * -b;
		hb = modes.phi * b;
	}
	forward = ea.adjoint() * ha;
	backward = eb.adjoint() * hb;
	forward /= omega;
	backward /= omega;
	S4R_TRACE("< Layer::GetPowerFlux\n");
}

void S4r::Layer::DumpDescription(
	const PeriodicMesh *mesh,
	FILE *fp
) const{
	FILE *f = fp;
	if(NULL == fp){ f = stdout; }
	
	double scale[2] = {4*72, 4*72};
	fprintf(f, "%f %f scale\n", scale[0], scale[1]);
	fprintf(f, "%f setlinewidth\n", 1./scale[0]);
	fprintf(f, "%f %f translate\n", 8.5*0.5*72/scale[0], 11*0.5*72/scale[1]);
	// Draw the origin cross
	fprintf(f, "%% origin cross\n0.5 setgray\n");
	fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", -0.05, 0.00, 0.05, 0.00);
	fprintf(f, "newpath %f %f moveto %f %f lineto stroke\n", 0.00, -0.05, 0.00, 0.05);
	
	fprintf(f, "%% mesh polygons\n0.5 setgray\n");
	for(size_t i = 0; i < mesh->NumFaces(); ++i){
		ConvexPolygon poly;
		mesh->FaceControlPolygon(i, &poly);
		
		// output the polygon
		fprintf(f, "newpath ");
		for(size_t j = 0; j < poly.v.size(); ++j){
			Vec2 v(poly.offset + poly.v[j]);
			fprintf(f, "%g %g", v[0], v[1]);
			if(j > 0){
				fprintf(f, " lineto ");
			}else{
				fprintf(f, " moveto ");
			}
		}
		fprintf(f, "closepath stroke\n");
	}
	
	/*
	Mat2 L = mesh->GetLattice();
	// Draw the clip marks
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		0.55*L(0,0)+0.50*L(0,1), 0.55*L(1,0)+0.50*L(1,1),
		0.50*L(0,0)+0.50*L(0,1), 0.50*L(1,0)+0.50*L(1,1),
		0.50*L(0,0)+0.55*L(0,1), 0.50*L(1,0)+0.55*L(1,1));
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		-0.50*L(0,0)+0.55*L(0,1), -0.50*L(1,0)+0.55*L(1,1),
		-0.50*L(0,0)+0.50*L(0,1), -0.50*L(1,0)+0.50*L(1,1),
		-0.55*L(0,0)+0.50*L(0,1), -0.55*L(1,0)+0.50*L(1,1));
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		-0.55*L(0,0)-0.50*L(0,1), -0.55*L(1,0)-0.50*L(1,1),
		-0.50*L(0,0)-0.50*L(0,1), -0.50*L(1,0)-0.50*L(1,1),
		-0.50*L(0,0)-0.55*L(0,1), -0.50*L(1,0)-0.55*L(1,1));
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto\nstroke\n",
		0.50*L(0,0)-0.55*L(0,1), 0.50*L(1,0)-0.55*L(1,1),
		0.50*L(0,0)-0.50*L(0,1), 0.50*L(1,0)-0.50*L(1,1),
		0.55*L(0,0)-0.50*L(0,1), 0.55*L(1,0)-0.50*L(1,1));
	*/
	
	fprintf(f, "%% pattern\n0 setgray\n");
	description.pattern.OutputPostscript(f);
	
	fprintf(f, "showpage\n");
}


void S4r::Layer::GetFields(
	const doublecomplex &omega,
	const PeriodicMesh *mesh,
	const CVec2 &phi,
	const Vec3 &r,
	CVec3 &E, CVec3 &H
) const{
/*
	CVec e, h;
	{
		CVec a, b;
		GetTranslatedSolution(r[2], a, b);
		e = modes.kp * modes.phi * modes.q.cwiseInverse().asDiagonal() * (a-b);
		h = modes.phi * (a+b);
	}
	
	E.setZero();
	H.setZero();
	
	PeriodicMesh::PeriodicIndex iface; // cache the face index for warm-start
	
	// Perform 1-form interpolation on the transverse fields
	{
		S4r::PeriodicMesh::EdgeCoeffs coeffs;
		iface = mesh->EdgeInterpolation(Vec2(r[0], r[1]), coeffs);
		
		for(S4r::PeriodicMesh::EdgeCoeffs::const_iterator i = coeffs.begin();
			i != coeffs.end(); ++i
		){
			doublecomplex t;
			for(unsigned j = 0; j < 2; ++j){
				t = e[i->idx] * i->c[j]; if(i->wx){ t *= phi[0]; } if(i->wy){ t *= phi[1]; }
				E[j] += t;
				t = h[i->idx] * i->c[j]; if(i->wx){ t *= phi[0]; } if(i->wy){ t *= phi[1]; }
				H[j] += t;
			}
		}
	}

	// The normal components are defined by
	//   curl(E_transverse) =  i omega mu_z  Hz
	//   curl(H_transverse) = -i omega eps_z Ez
	// There is one value of Ez per face, and one value of Hz per vertex
	
	{ // piecewise linear interpolation
		CVec hz = (doublecomplex(0,  1)/omega) * matrices.imu  * mesh->d0.adjoint() * e;
	
		S4r::PeriodicMesh::VertexCoeffs coeffs;
		iface = mesh->VertexInterpolation(Vec2(r[0], r[1]), coeffs, iface);
		
		for(S4r::PeriodicMesh::VertexCoeffs::const_iterator i = coeffs.begin();
			i != coeffs.end(); ++i
		){
			doublecomplex t = hz[i->i] * i->c;
			if(i->wx){ t *= phi[0]; } if(i->wy){ t *= phi[1]; }
			H[2] += t;
		}
	}
	{ // piecewise constant interpolation
		CVec ez = (doublecomplex(0, -1)/omega) * matrices.ieps * mesh->d1 * h;
		
		S4r::PeriodicMesh::FaceCoeffs coeffs;
		iface = mesh->FaceInterpolation(Vec2(r[0], r[1]), coeffs, iface);
		
		for(S4r::PeriodicMesh::FaceCoeffs::const_iterator i = coeffs.begin();
			i != coeffs.end(); ++i
		){
			doublecomplex t = ez[i->i] * i->c;
			if(i->wx){ t *= phi[0]; } if(i->wy){ t *= phi[1]; }
			E[2] += t;
		}
	}*/
}
