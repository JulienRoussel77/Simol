#pragma once
//#include "SparseMatrix.hpp"
//#include <Eigen/Dense>

#include "tools.hpp"
#include "input.hpp"
#include "potential.hpp"
#include "basis.hpp"


//#include "eigen.hpp"

namespace simol
{
	
	SMat kron(const SMat& A, const SMat& B);
	DMat kron(const DMat& A, const DMat& B);
	
	class Galerkin;
  
  Galerkin* createLangevinGalerkin(Input const& input);
	
	class Galerkin
	{
	friend Galerkin* createLangevinGalerkin(Input const& input);
	
	protected:
		//DenseMatrix<double> A;
		int nbOfParticles_;
		size_t nbOfFourier_, nbOfHermite_, maxOfFourier_;
		size_t sizeOfBasis_;
		SMat SIdQ_, SIdP_;
		DMat DIdQ_, DIdP_;
		SMat Q_, P_;
		SMat tQ_, tP_;
		SMat Lthm0_;
		SMat Lthm_, Lham_;
		SMat L1_;
		SMat Leq_, Leta_;
		double beta_, gamma_;
		double amplitude_;
		double externalForce_;
		size_t nbOfIntegrationNodes_;
		//vector<double> expFourierCoeffs_;
		DMat trigToExpMat_, expToTrigMat_;
		DMat trigToExpTens_, expToTrigTens_;
		Potential* potential_;
		ExpFourierHermiteBasis basis_;
	public:
		Galerkin(Input const& input);
		
		virtual int nbOfVariables() const;
		virtual int nbOfParticles() const;
		
		const double& expFourierCoeffs(int iOfElt) const;
		size_t iTens(size_t iOfFourier2, size_t iOfHermite) const;
		DMat shapeSaddle(const DMat& A) const;
		DMat unshapeSaddle(const DMat& Asad) const;
		DVec shapeSaddle(const DVec& X) const;
		DVec unshapeSaddle(const DVec& Xsad) const;
		DVec solveWithSaddle(const DMat& A, const DVec& X) const;
		DVec solveWithSaddle(const SMat& A, const DVec& X) const;
		DMat invWithSaddle(const DMat& A) const;
		
		//void computeFourierCoeffsExp();
		void computeExpToTrigMat();
		virtual void computeExpToTrigTens() = 0;
		DMat convertToTrigBasis(const DMat& X);
		void createQ();
		void createP();
		void createLthm0();
		virtual void createLthm() = 0;
		virtual void compute();
		
		DVec gettGiHj(int i, int j) const;
		DVec gettGiHjTrig(int i, int j) const;
		DVec getLtGiHj(int i, int j) const;
		DVec getLtGiHjTrig(int i, int j) const;
		DVec getLinvtGiHj(int i, int j) const;
		DVec getLinvtGiHjTrig(int i, int j) const;
		SMat CVcoeffs() const;
	};
	
	class LangevinGalerkin : public Galerkin
	{
	public:
		LangevinGalerkin(Input const& input);
		
		virtual void computeExpToTrigTens();
		virtual void createLthm();
	};
	
	class BoundaryLangevinGalerkin : public Galerkin
	{
		SMat SId_;
	public:
		BoundaryLangevinGalerkin(Input const& input);
		size_t iTens(size_t iOfFourier2, size_t iOfHermite, int iOfParticle) const;
		SMat PMatToTens(SMat const& PMat, int iOfParticleP);
		SMat doubleMatToTens(SMat const& QMat, SMat const& PMat, int iOfParticleQ, int iOfParticleP);
		virtual void computeExpToTrigTens();
		void createLham();
		virtual void createLthm();
		
		virtual void compute();
	};
	

}