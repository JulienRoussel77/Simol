#include "Basis.hpp"

namespace simol
{
	
	TVec::TVec(vector<size_t>& nbOfElts0):
		nbOfElts_(nbOfElts0)
	{}
	
	size_t TVec::nbOfVariables()
	{
		return nbOfElts_.size();
	}
	
	size_t const& TVec::nbOfElts(size_t iOfVariable) const
	{
		return nbOfElts_[iOfVariable];
	}
	
	size_t& TVec::nbOfElts(size_t iOfVariable)
	{
		return nbOfElts_[iOfVariable];
	}
	
	vector<size_t> const& TVec::nbOfElts() const
	{
		return nbOfElts_;
	}
	
	vector<size_t>& TVec::nbOfElts()
	{
		return nbOfElts_;
	}
	
	size_t TVec::iTens(vector<size_t>& vecOfElt) const
	{
		size_t iTensOfElt = vecOfElt[vecOfElt.size() - 1];
		for (int iOfVariable = vecOfElt.size() - 2; iOfVariable >= 0; iOfVariable--)
		{
			iTensOfElt *= nbOfElts(iOfVariable+1);
			iTensOfElt += vecOfElt[iOfVariable];
		}
		return iTensOfElt;
	}
	
	
	
	DTVec::DTVec(vector<size_t>& nbOfElts0):
		TVec(nbOfElts0),
		data_(product(nbOfElts0))
	{}
	
	size_t DTVec::size() const
	{
		return data_.size();
	}
	
	const double& DTVec::operator()(vector<size_t>& vecIndex) const
	{
		return data_(iTens(vecIndex));
	}
	
	double& DTVec::operator()(vector<size_t>& vecIndex)
	{
		return data_(iTens(vecIndex));
	}
	
	const double& DTVec::operator()(size_t iTensOfElt) const
	{
		return data_(iTensOfElt);
	}
	
	double& DTVec::operator()(size_t iTensOfElt)
	{
		return data_(iTensOfElt);
	}
	
	size_t product(vector<size_t>& nbOfElts)
	{
		size_t result = 1;
		for (int iOfVariable = 0; iOfVariable < (int)nbOfElts.size(); iOfVariable++)
		{
			result *= nbOfElts[iOfVariable];
		}
		return result;
	}
	
	/*DTVec product(DMat& A, DTVec& X, size_t iOfVariable)
	{
		DTVec AX(X.nbOfElts());
		return AX;
	}*/
	
	
	
	Basis::Basis(const size_t nbOfElts):
		nbOfElts_(nbOfElts)
	{}

	size_t const& Basis::nbOfElts() const
	{
		return nbOfElts_;
	}
	
	size_t& Basis::nbOfElts()
	{
		return nbOfElts_;
	}
	
	
	FourierBasis::FourierBasis(const size_t nbOfElts):
		Basis(nbOfElts)
	{
		assert(nbOfElts % 2 == 1);
	}
	
	size_t FourierBasis::nbOfFreq() const
	{
		return (nbOfElts_-1)/2;
	}
	
	double FourierBasis::value(double variable, const int iOfElt) const
	{
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return 1;
		else if (iOfElt % 2 == 1)
			return sqrt(2.) * sin(iOfFreq * variable);
		else
			return sqrt(2.) * cos(iOfFreq * variable);
	}
	
	Vector<double> FourierBasis::gradient(double variable, const int iOfElt) const
	{
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return Vector<double>(1, 0);
		else if (iOfElt % 2 == 1)
			return Vector<double>(1,   sqrt(2.) * iOfFreq * cos(iOfFreq * variable));
		else
			return Vector<double>(1, - sqrt(2.) * iOfFreq * sin(iOfFreq * variable));
	}
	
	double FourierBasis::laplacian(double variable, const int iOfElt) const
	{
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return 1;
		else if (iOfElt % 2 == 1)
			return - sqrt(2.) * pow(iOfFreq, 2) * sin(iOfFreq * variable);
		else
			return - sqrt(2.) * pow(iOfFreq, 2) * cos(iOfFreq * variable);
	}
	
	ExpFourierBasis::ExpFourierBasis(const size_t nbOfElts0, double beta0, Potential& potential0):
		Basis(nbOfElts0),
		beta_(beta0),
		potential_(&potential0),
		nbOfIntegrationNodes_(10000),
		qRepartitionFct_(0),
		expFourierCoeffs_(2 * nbOfElts0)
	{
		//cout << "ExpFourierBasis(const size_t nbOfElts0, double beta0, Potential& potential0)" << endl;
		assert(nbOfElts0 % 2 == 1);
				
	//We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
	//where C = ExpFourierBasis::basisCoefficient_
		double step = 2 * M_PI / (double)nbOfIntegrationNodes_;
		for (size_t iOfNode = 0; iOfNode < nbOfIntegrationNodes_; iOfNode++)
		{
			double q = - M_PI + iOfNode * step;
			qRepartitionFct_ += exp(-beta_ * potential(q));
			expFourierCoeffs_[0] += sqrt(2) * exp(-beta_ * potential(q)/2);
			for (size_t iOfFourier=1; iOfFourier <= 2*nbOfFreq(); iOfFourier++)
			{
				//cout << q << " -> " << potential(q) << endl;
				expFourierCoeffs_[2 * iOfFourier] += sqrt(2) * cos(iOfFourier * q) * exp(-beta_ * potential(q)/2);
				//if (iOfFourier == 1)
					//cout << q << " " << expFourierCoeffs_[2 * iOfFourier] << endl;
			}
			for (size_t iOfFourier=1; iOfFourier <= 2*nbOfFreq(); iOfFourier++)
				expFourierCoeffs_[2 * iOfFourier-1] += sqrt(2) * sin(iOfFourier * q) * exp(-beta_ * potential(q)/2);
		}
		qRepartitionFct_ *= step / (2*M_PI);
		basisCoefficient_ = sqrt(qRepartitionFct_);
		for (size_t iOfFourier2=0; iOfFourier2 <= 4 * nbOfFreq(); iOfFourier2++)
		{
			expFourierCoeffs_[iOfFourier2] *= step / (2*M_PI*basisCoefficient_);
			//cout << "coeff n°"<< iOfFourier2 << " = " << expFourierCoeffs_[iOfFourier2] << endl;
		}
		//cout << "qRepartitionFct_ = " << qRepartitionFct_ << endl;
		//cout << "basisCoefficient_ = " << basisCoefficient_ << endl;
	}
	
	size_t ExpFourierBasis::nbOfFreq() const
	{
		return (nbOfElts_-1)/2;
	}
	
	const double& ExpFourierBasis::expFourierCoeffs(int iOfElt) const
	{
		return expFourierCoeffs_[iOfElt];
	}
	
	double ExpFourierBasis::potential(double variable) const
	{
		return potential_->value(variable);
	}
	
	double ExpFourierBasis::potDeriv(double variable) const
	{
		return (potential_->gradient(variable))(0);
	}
	
	double ExpFourierBasis::potLapla(double variable) const
	{
		return (potential_->laplacian(variable));
	}
	
	double ExpFourierBasis::value(double variable, const int iOfElt) const
	{
		double expo = basisCoefficient_ * exp(beta_*potential(variable)/2);
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return expo;
		else if (iOfElt % 2 == 1)
			return sqrt(2) * sin(iOfFreq * variable) * expo;
		else
			return sqrt(2) * cos(iOfFreq * variable) * expo;
	}
	
	Vector<double> ExpFourierBasis::gradient(double variable, const int iOfElt) const
	{
		//cout << basisCoefficient_ << " " << exp(beta_*potential(variable)/2) << " " << potDeriv(variable) << endl;
		double expo = basisCoefficient_ * exp(beta_*potential(variable)/2);
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return Vector<double>(1, beta_/2 * potDeriv(variable) * expo);
		else if (iOfElt % 2 == 1)
			return Vector<double>(1, sqrt(2.) * ( iOfFreq*cos(iOfFreq * variable) + sin(iOfFreq * variable) * beta_/2 * potDeriv(variable)) * expo);
		else
			return Vector<double>(1, sqrt(2.) * (-iOfFreq*sin(iOfFreq * variable) + cos(iOfFreq * variable) * beta_/2 * potDeriv(variable)) * expo);
	}
	
	double ExpFourierBasis::laplacian(double variable, const int iOfElt) const
	{
		double expo = basisCoefficient_ * exp(potential(variable)/2);
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return beta_/2 * (potLapla(variable) + beta_/2 * pow(potDeriv(variable) , 2)) * expo;
		else if (iOfElt % 2 == 1)
			return sqrt(2.) * (beta_ * iOfFreq * potDeriv(variable) * cos(iOfFreq * variable) 
												+ (beta_/2 * potLapla(variable) + pow(beta_, 2)/4 * pow(potDeriv(variable), 2) - pow(iOfFreq, 2)) * sin(iOfFreq * variable)) * expo;
		else
			return sqrt(2.) * (- beta_ * iOfFreq * potDeriv(variable) * sin(iOfFreq * variable) 
												+ (beta_/2 * potLapla(variable) + pow(beta_, 2)/4 * pow(potDeriv(variable), 2) - pow(iOfFreq, 2)) * cos(iOfFreq * variable)) * expo;
	}
	
	HermiteBasis::HermiteBasis(const size_t nbOfElts0, double beta0)
  : Basis(nbOfElts0),
		beta_(beta0),
		polyCoeffs_(zero<double>(nbOfElts0, nbOfElts0))
		//polyCoeffs_(DenseMatrix<double>::Zero(nbOfElts0, nbOfElts0))
	{
		//cout << "HermiteBasis(const size_t nbOfElts0)" << endl;
		polyCoeffs_(0,0) = 1;
		polyCoeffs_(1,1) = beta_;
		for (int iOfElt = 2; iOfElt < (int)nbOfElts0; iOfElt++)
			for (int iOfCoeff=0; iOfCoeff <= iOfElt; iOfCoeff++)
			{
				polyCoeffs_(iOfElt, iOfCoeff) = (iOfCoeff?(sqrt(beta_/iOfElt)*polyCoeffs_(iOfElt-1, iOfCoeff-1)):0) - sqrt((iOfElt-1) / (double)iOfElt) * polyCoeffs_(iOfElt-2, iOfCoeff);
				//polyCoeffs_(iOfElt, iOfCoeff) /= sqrt(iOfElt);				
			}
			//cout << "Hermite coeffs : " << endl;
			//cout << polyCoeffs_ << endl;
	}
	
	double HermiteBasis::value(double variable, const int iOfElt) const
	{
		double result = 0;
		for (int iOfCoeff=0; iOfCoeff <= iOfElt; iOfCoeff++)
		{
			result += polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff);
			//cout << polyCoeffs_(iOfElt, iOfCoeff) << endl;
		}
		//cout << variable << " , " << iOfElt << " -> " << result << endl;
		return result;
	}
	
	Vector<double> HermiteBasis::gradient(double variable, const int iOfElt) const
	{
		double result = 0;
		for (int iOfCoeff=1; iOfCoeff <= iOfElt; iOfCoeff++)
			result += iOfCoeff * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff-1);
		return Vector<double>(1, result);
	}
	
	double HermiteBasis::laplacian(double variable, const int iOfElt) const
	{
		double result = 0;
		for (int iOfCoeff=2; iOfCoeff <= iOfElt; iOfCoeff++)
			result += iOfCoeff * (iOfCoeff-1) * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff-2);
		return result;
	}
	
			
	
	TensorBasis::TensorBasis(const size_t nbOfVariables0):
		bases_(nbOfVariables0)
	{
	}
	
	TensorBasis::~TensorBasis()
	{
		for (int iOfVariable = 0; iOfVariable < (int)nbOfVariables(); iOfVariable++)
			delete bases_[iOfVariable];
	}
	
	size_t TensorBasis::nbOfVariables() const
	{
		return bases_.size();
	}

	size_t const& TensorBasis::nbOfElts(const int iOfVariable) const
	{
		return bases_[iOfVariable]->nbOfElts();
	}
	
	size_t& TensorBasis::nbOfElts(const int iOfVariable)
	{
		return bases_[iOfVariable]->nbOfElts();
	}
	
	vector<size_t> TensorBasis::nbOfElts() const
	{
		vector<size_t> nbOfElts0(nbOfVariables());
		for (int iOfVariable = 0; iOfVariable < (int)nbOfVariables(); iOfVariable++)
			nbOfElts0[iOfVariable] = bases_[iOfVariable]->nbOfElts();
		return nbOfElts0;
	}
	
	/*double TensorBasis::value(double variable, const int iOfElt)
	{
		int iOfFreq = (iOfElt+1)/2;
		if (iOfElt == 0)
			return 1;
		else if (iOfElt % 2 == 1)
			return sqrt(2.) * sin(iOfFreq * variable);
		else
			return sqrt(2.) * cos(iOfFreq * variable);
	}*/
	
	QPBasis::QPBasis():
		TensorBasis(2)
	{
		//cout << "QPBasis()" << endl;
	}
	
	/*QPBasis::QPBasis(const size_t nbOfFourier0, const size_t nbOfHermite0):
		TensorBasis(2)
	{
		bases_[0] = new ExpFourierBasis(nbOfFourier0);
		bases_[1] = new HermiteBasis(nbOfHermite0);
	}*/
	
	size_t& QPBasis::nbOfFourier()
	{
		return nbOfElts(0);
	}
	
	const size_t& QPBasis::nbOfFourier() const
	{
		return nbOfElts(0);
	}
	
	size_t& QPBasis::nbOfHermite()
	{
		return nbOfElts(1);
	}
	
	const size_t& QPBasis::nbOfHermite() const
	{
		return nbOfElts(1);
	}
	
	//psi = (1,1  1,2  ...  1,N_H  2,1 ... )
	//N_H blocks of size N_G (we concatene the columns of the matrix)
	size_t QPBasis::iTens(size_t iOfFourier2, size_t iOfHermite) const
	{
		assert(iOfFourier2 < nbOfFourier()	&& iOfHermite < nbOfHermite());
		return nbOfFourier() * iOfHermite + iOfFourier2;
	}
	
	vector<size_t> QPBasis::vecTens(size_t iTens0) const
	{
		assert(iTens0 < nbOfHermite() * nbOfFourier());
		vector<size_t> vecIndex(2);
		vecIndex[0] = iTens0 % nbOfFourier();
		vecIndex[1] = iTens0 / nbOfFourier();
		return vecIndex;
	}
	
	double QPBasis::value(vector<Particle> const& configuration, const size_t iOfElt) const
	{
		vector<size_t> vecIndex = vecTens(iOfElt);
		return value(configuration, vecIndex);
	}
	
	double QPBasis::value(vector<Particle> const& configuration, vector<size_t>& vecIndex) const
	{
		return bases_[0]->value(configuration[0].position(0), vecIndex[0]) * bases_[1]->value(configuration[0].momentum(0), vecIndex[1]);
	}
	
	Vector<double> QPBasis::gradientQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, vector<size_t>& vecIndex) const
	{
		//cout << "QPBasis::gradientQ" << endl;
		//cout << vecIndex[0] << " " << vecIndex[1] << endl;
		//cout << bases_[0]->gradient(configuration[0].position(0), vecIndex[0]) << " X "
		//		<< bases_[1]->value(configuration[0].momentum(0), vecIndex[1]) << endl;
		return bases_[0]->gradient(configuration[0].position(0), vecIndex[0]) 
					* bases_[1]->value(configuration[0].momentum(0), vecIndex[1]);
	}	
		
	Vector<double> QPBasis::gradientQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const
	{
		vector<size_t> vecIndex = vecTens(iOfCoeff);
		return gradientQ(configuration, iOfParticle, vecIndex);
	}
		
	double QPBasis::laplacianQ(vector<Particle> const& configuration, size_t /*iOfParticle*/, vector<size_t>& vecIndex) const
	{
		return bases_[0]->laplacian(configuration[0].position(0), vecIndex[0]) 
					* bases_[1]->value(configuration[0].momentum(0), vecIndex[1]);
	}
	
	double QPBasis::laplacianQ(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const
	{
		vector<size_t> vecIndex = vecTens(iOfCoeff);
		return laplacianQ(configuration, iOfParticle, vecIndex);
	}
	
	Vector<double> QPBasis::gradientP(vector<Particle> const& configuration, size_t /*iOfParticle*/, vector<size_t>& vecIndex) const
	{
		//cout << bases_[0]->value(configuration[0].position(0), vecIndex[0]) << " X "
		//		<< bases_[1]->gradient(configuration[0].momentum(0), vecIndex[1]) << endl;
		return bases_[0]->value(configuration[0].position(0), vecIndex[0]) 
					* bases_[1]->gradient(configuration[0].momentum(0), vecIndex[1]);
	}
	
	Vector<double> QPBasis::gradientP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const
	{
		vector<size_t> vecIndex = vecTens(iOfCoeff);
		return gradientP(configuration, iOfParticle, vecIndex);
	}
	
	double QPBasis::laplacianP(vector<Particle> const& configuration, size_t /*iOfParticle*/, vector<size_t>& vecIndex) const
	{
		return bases_[0]->value(configuration[0].position(0), vecIndex[0]) 
					* bases_[1]->laplacian(configuration[0].momentum(0), vecIndex[1]);
	}
	
	double QPBasis::laplacianP(vector<Particle> const& configuration, size_t iOfParticle, size_t iOfCoeff) const
	{
		vector<size_t> vecIndex = vecTens(iOfCoeff);
		return laplacianP(configuration, iOfParticle, vecIndex);	
	}
	
	
	FourierHermiteBasis::FourierHermiteBasis(Input const& input):
		QPBasis()
	{
		bases_[0] = new FourierBasis(input.nbOfFourier());
		bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
	}
	
	ExpFourierHermiteBasis::ExpFourierHermiteBasis(Input const& input, Potential& potential):
		QPBasis()
	{
		//cout << "ExpFourierHermiteBasis(Input const& input)" << endl;
		bases_[0] = new ExpFourierBasis(input.nbOfFourier(), input.beta(), potential);
		bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
		//cout << "end ExpFourierHermiteBasis(Input const& input)" << endl;
	}
	
	const double& ExpFourierHermiteBasis::expFourierCoeffs(int iOfElt) const
	{
		return bases_[0]->expFourierCoeffs(iOfElt);
	}
}












