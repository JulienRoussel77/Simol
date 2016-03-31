#ifndef SIMOL_GALERKIN_HPP
#define SIMOL_GALERKIN_HPP

#include "Tools.hpp"
#include "Input.hpp"
#include "Potential.hpp"
#include "Basis.hpp"
#include "core/linalg/DenseMatrix.hpp"

namespace simol
{

	SMat kron(const SMat& A, const SMat& B);
	DenseMatrix<double> kron(const DenseMatrix<double>& A, const DenseMatrix<double>& B);

	class Galerkin;

  Galerkin* createLangevinGalerkin(Input const& input);

	class Galerkin
	{
	friend Galerkin* createLangevinGalerkin(Input const& input);

	protected:
		int nbOfParticles_;
		size_t nbOfFourier_, nbOfHermite_, maxOfFourier_;
		size_t sizeOfBasis_;
		SMat SIdQ_, SIdP_;
		DenseMatrix<double> DIdQ_, DIdP_;
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
		DenseMatrix<double> trigToExpMat_, expToTrigMat_;
		DenseMatrix<double> trigToExpTens_, expToTrigTens_;
		Potential* potential_;
		ExpFourierHermiteBasis basis_;
	public:
		Galerkin(Input const& input);

		virtual int nbOfVariables() const;
		virtual int nbOfParticles() const;

		const double& expFourierCoeffs(int iOfElt) const;
		size_t iTens(size_t iOfFourier2, size_t iOfHermite) const;
		DenseMatrix<double> shapeSaddle(const DenseMatrix<double>& A) const;
		DenseMatrix<double> unshapeSaddle(const DenseMatrix<double>& Asad) const;
		DVec shapeSaddle(const DVec& X) const;
		DVec unshapeSaddle(const DVec& Xsad) const;
		DVec solveWithSaddle(const DenseMatrix<double>& A, const DVec& X) const;
		DVec solveWithSaddle(const SMat& A, const DVec& X) const;
    DenseMatrix<double> invWithSaddle(const SparseMatrix<double>& A) const;
		DenseMatrix<double> invWithSaddle(const DenseMatrix<double>& A) const;

		void computeExpToTrigMat();
		virtual void computeExpToTrigTens() = 0;
		DenseMatrix<double> convertToTrigBasis(const DenseMatrix<double>& X);
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
		//SMat SId_;
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

#endif
