#ifndef SIMOL_BASIS_HPP
#define SIMOL_BASIS_HPP

#include "tools.hpp"
#include "potential.hpp"
#include "particle.hpp"

namespace simol
{
	class TVec
	{
		vector<size_t> nbOfElts_;
	public:
		//virtual const double& operator()(vector<size_t>& vecOfElt) const = 0;
		//virtual double& operator()(vector<size_t>& vecOfElt) = 0;
		TVec(vector<size_t>& nbOfElts0);
		size_t nbOfVariables();
		size_t const& nbOfElts(size_t iOfVariable) const;
		size_t& nbOfElts(size_t iOfVariable);
		vector<size_t> const& nbOfElts() const;
		vector<size_t>& nbOfElts();
		virtual size_t size() const = 0;
		size_t iTens(vector<size_t>& vecOfElt) const;
	};

	class DTVec : public TVec
	{
		Vector<double> data_;
	public:
		DTVec(vector<size_t>& nbOfElts0);
		virtual size_t size() const;
		const double& operator()(vector<size_t>& vecIndex) const;
		double& operator()(vector<size_t>& vecIndex);
		const double& operator()(size_t iTensOfElt) const;
		double& operator()(size_t iTensOfElt);
	};

	size_t product(vector<size_t>& nbOfElts);
	//DTVec product(DMat& A, DTVec& X, size_t iOfVariable);

	class Basis
	{
	protected:
		size_t nbOfElts_;
	public:
		Basis(const size_t nbOfElts);
		virtual ~Basis(){};
		virtual size_t const& nbOfElts() const;
		virtual size_t& nbOfElts();
		virtual const double& expFourierCoeffs(int /*iOfElt*/) const {assert(false);};
		virtual double value(double variable, const int iOfElt) const = 0;
		virtual Vector<double> gradient(double variable, const int iOfElt) const = 0;
		virtual double laplacian(double variable, const int iOfElt) const = 0;
	};

	class FourierBasis : public Basis
	{
	public:
		FourierBasis(const size_t nbOfElts);
		virtual size_t nbOfFreq() const;
		virtual double value(double variable, const int iOfElt) const;
		virtual Vector<double> gradient(double variable, const int iOfElt) const;
		virtual double laplacian(double variable, const int iOfElt) const;
	};

	class ExpFourierBasis : public Basis
	{
		double beta_;
		Potential* potential_;
		size_t nbOfIntegrationNodes_;
		double qRepartitionFct_, basisCoefficient_;
		vector<double> expFourierCoeffs_;
	public:
		ExpFourierBasis(const size_t nbOfElts, double beta0, Potential* potential);
		virtual size_t nbOfFreq() const;
		const double& expFourierCoeffs(int iOfElt) const;
		virtual double potential(double variable) const;
		virtual double potDeriv(double variable) const;
		virtual double potLapla(double variable) const;
		virtual double value(double variable, const int iOfElt) const;
		virtual Vector<double> gradient(double variable, const int iOfElt) const;
		virtual double laplacian(double variable, const int iOfElt) const;
	};



	class HermiteBasis : public Basis
	{
		double beta_;
		DMat polyCoeffs_;
	public:
		HermiteBasis(const size_t nbOfElts, double beta0);
		virtual double value(double variable, const int iOfElt) const;
		virtual Vector<double> gradient(double variable, const int iOfElt) const;
		virtual double laplacian(double variable, const int iOfElt) const;
	};

	class TensorBasis
	{
	protected:
		vector<Basis*> bases_;
		//vector<size_t> nbOfElts_;
	public:
		TensorBasis(const size_t nbOfVariables);
		virtual ~TensorBasis();
		virtual size_t nbOfVariables() const;
		virtual size_t const& nbOfElts(const int iOfVariable) const;
		virtual size_t& nbOfElts(const int iOfVariable);
		virtual vector<size_t> nbOfElts() const;
		virtual double value(vector<Particle> const& configuration, const size_t iOfElt) const = 0;
		virtual double value(vector<Particle> const& configuration, vector<size_t>& vecIndex) const = 0;
		virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const = 0;
		virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const = 0;
		virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const = 0;
		virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const = 0;
		virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const = 0;
		virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const = 0;
		virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const = 0;
		virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const = 0;
	};

	class QPBasis : public TensorBasis
	{
	public:
		QPBasis();
		size_t& nbOfFourier();
		const size_t& nbOfFourier() const;
		size_t& nbOfHermite();
		const size_t& nbOfHermite() const;
		size_t iTens(size_t iOfFourier2, size_t iOfHermite) const;
		vector<size_t> vecTens(size_t iTens0) const;
		virtual double value(vector<Particle> const& configuration, const size_t iOfElt) const;
		virtual double value(vector<Particle> const& configuration, vector<size_t>& vecIndex) const;
		virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const;
		virtual Vector<double> gradientQ(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const;
		virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const;
		virtual double laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const;
		virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const;
		virtual Vector<double> gradientP(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const;
		virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const;
		virtual double laplacianP(vector<Particle> const& configuration, size_t iOfParticle, vector<size_t>& vecIndex) const;
	};

	class FourierHermiteBasis : public QPBasis
	{
	public:
		FourierHermiteBasis(Input const& input);
	};

	class ExpFourierHermiteBasis : public QPBasis
	{
	public:
		ExpFourierHermiteBasis(Input const& input, Potential* potential);
		const double& expFourierCoeffs(int iOfElt) const;
	};

	/*class TrigBasis : public TensorBasis
	{
	public:
		TensorBasis(const size_t nbOfElts1, const size_t nbOfElts2);
		virtual double valueOfElt(double variable, const int iOfElt, const int iOfVariable);
	};*/
}

#endif
