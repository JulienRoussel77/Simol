#include "simol/statphys/controlVariate/Basis.hpp"

namespace simol
{

  TVec::TVec(vector<int>& nbOfElts0):
    nbOfElts_(nbOfElts0)
  {}

  int TVec::nbOfVariables()
  {
    return nbOfElts_.size();
  }

  int const& TVec::nbOfElts(int iOfVariable) const
  {
    return nbOfElts_[iOfVariable];
  }

  int& TVec::nbOfElts(int iOfVariable)
  {
    return nbOfElts_[iOfVariable];
  }

  vector<int> const& TVec::nbOfElts() const
  {
    return nbOfElts_;
  }

  vector<int>& TVec::nbOfElts()
  {
    return nbOfElts_;
  }

  int TVec::iTens(vector<int>& vecOfElt) const
  {
    int iTensOfElt = vecOfElt[vecOfElt.size() - 1];
    for (int iOfVariable = vecOfElt.size() - 2; iOfVariable >= 0; iOfVariable--)
    {
      iTensOfElt *= nbOfElts(iOfVariable + 1);
      iTensOfElt += vecOfElt[iOfVariable];
    }
    return iTensOfElt;
  }



  DTVec::DTVec(vector<int>& nbOfElts0):
    TVec(nbOfElts0),
    data_(product(nbOfElts0))
  {}

  int DTVec::size() const
  {
    return data_.size();
  }

  double const& DTVec::operator()(vector<int>& vecIndex) const
  {
    return data_(iTens(vecIndex));
  }

  double& DTVec::operator()(vector<int>& vecIndex)
  {
    return data_(iTens(vecIndex));
  }

  double const& DTVec::operator()(int iTensOfElt) const
  {
    return data_(iTensOfElt);
  }

  double& DTVec::operator()(int iTensOfElt)
  {
    return data_(iTensOfElt);
  }

  int product(vector<int>& nbOfElts)
  {
    int result = 1;
    for (int iOfVariable = 0; iOfVariable < (int)nbOfElts.size(); iOfVariable++)
    {
      result *= nbOfElts[iOfVariable];
    }
    return result;
  }

  /*DTVec product(DMat& A, DTVec& X, int iOfVariable)
  {
    DTVec AX(X.nbOfElts());
    return AX;
  }*/
  
  
  
  
  
  

  /*Basis* createBasis(Input const& input, Potential& potential)
  {
    if (input.galerkinElts() == "ExpFourierHermite")
      return new ExpFourierHermiteBasis(input, potential);
    //throw runtime_error("Invalid type of basis !");
    return nullptr;
  }*/

  Basis::Basis(const int nbOfElts):
    nbOfElts_(nbOfElts)
  {}

  int const& Basis::nbOfElts() const
  {
    return nbOfElts_;
  }

  int& Basis::nbOfElts()
  {
    return nbOfElts_;
  }

  double const& Basis::expFourierMeans(int /*iOfElt*/) const
  {
    throw std::runtime_error("Basis::expFourierMeans is not implemented");
  };
  
  DVec const& Basis::gVector() const
  {
    throw std::runtime_error("Basis::expFourierMeans is not implemented");
  };
  
  double const& Basis::norm2gVector() const
  {
    throw std::runtime_error("Basis::sumExpFourierMeans is not implemented");
  };
  
  
  // Compute <H_n, \partial_p^* \partial_p H_m>
  double Basis::xGradStarY(const int iOfEltLeft, const int iOfEltRight) const
  {
    return xGradY(iOfEltRight, iOfEltLeft);
  }
  
  void Basis::gradStarMatrix(SMat& A) const
  {
    gradMatrix(A);
    A = A.adjoint();
  }
  
  ExpFourierBasis::ExpFourierBasis(const int nbOfElts0, double beta0, Potential& potential0):
    Basis(nbOfElts0),
    beta_(beta0),
    potential_(&potential0),
    nbOfIntegrationNodes_(1000*nbOfElts()),
    //nbOfIntegrationNodes_(10),
    qRepartitionFct_(0),
    expFourierMeans_(DVec::Zero(2 * nbOfElts0)),
    gVector_(DVec::Zero(nbOfElts0)),
    norm2gVector_(0)
  {
    //cout << "ExpFourierBasis(const int nbOfElts0, double beta0, Potential& potential0)" << endl;
    //assert(nbOfElts0 % 2 == 1);
    //cout << "ExpFourierBasis : nbOfElts = " << nbOfElts0 << " and " <<  nbOfElts0 % 2 << endl << endl;
    
    computeExpFourierMeans();
  }

  /*int ExpFourierBasis::nbOfFreq() const
  {
    return (nbOfElts_ - 1) / 2;
  }*/
  
  const double& ExpFourierBasis::amplitude() const
  {
    return potential_->parameter1();
  }
  

  double const& ExpFourierBasis::expFourierMeans(int iOfElt) const
  {
    return expFourierMeans_(iOfElt);
  }
  
  DVec const& ExpFourierBasis::gVector() const
  {
    return gVector_;
  }
  
  double const& ExpFourierBasis::norm2gVector() const
  {
    return norm2gVector_;
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
  
  ///We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  ///More precisely a_k = 1 / pi int cos(k q) C^-1 exp(-beta V / 2) dq and b_k = 1 / pi int sin(k q) C^-1 exp(-beta V / 2) dq
  ///where C = ExpFourierBasis::basisCoefficient_ Note that a_0 has the same normalization as a_k so that a_{k-l} = a_0 when k=l
  void ExpFourierBasis::computeExpFourierMeans()
  {
    double step = 2 * M_PI / (double)nbOfIntegrationNodes_;
    for (int iOfNode = 0; iOfNode < nbOfIntegrationNodes_; iOfNode++)
    {
      double q = - M_PI + iOfNode * step;
      qRepartitionFct_ += exp(-beta_ * potential(q)) * step;
      //expFourierMeans_(0) += exp(-beta_ * potential(q) / 2);
      for (int iOfFourier = 0; iOfFourier < 2 * nbOfElts(); iOfFourier++)
        if (iOfFourier % 2 == 0) expFourierMeans_(iOfFourier) += cos(iOfFourier/2 * q) * exp(-beta_ * potential(q) / 2);
        else expFourierMeans_(iOfFourier) += sin((iOfFourier+1)/2 * q) * exp(-beta_ * potential(q) / 2);
      
      /*for (int iOfFourier = 1; iOfFourier <= 2 * nbOfFreq(); iOfFourier++)
        expFourierMeans_(2 * iOfFourier) += sqrt(2) * cos(iOfFourier * q) * exp(-beta_ * potential(q) / 2);
      for (int iOfFourier = 1; iOfFourier <= 2 * nbOfFreq(); iOfFourier++)
        expFourierMeans_(2 * iOfFourier - 1) += sqrt(2) * sin(iOfFourier * q) * exp(-beta_ * potential(q) / 2);*/
    }
    basisCoefficient_ = sqrt(qRepartitionFct_ / (2 * M_PI));
    
    expFourierMeans_ *= step / (M_PI * basisCoefficient_);
    
    // gVector is the projection of the constant fct 1, ie it contains the mean of the expfourier elements
    gVector_ = expFourierMeans_.head(nbOfElts()) / sqrt(2.);
    gVector_(0) /= sqrt(2.);
    norm2gVector_ = pow(gVector_.norm(), 2);
    cout << "Squared norm of vector g = " << norm2gVector_ << endl;
  }
  
  

  double ExpFourierBasis::value(double variable, const int iOfElt) const
  {
    double expo = basisCoefficient_ * exp(beta_ * potential(variable) / 2);
    int iOfFreq = (iOfElt + 1) / 2;
    if (iOfElt == 0)
      return expo;
    else if (iOfElt % 2 == 1)
      return sqrt(2) * sin(iOfFreq * variable) * expo;
    else
      return sqrt(2) * cos(iOfFreq * variable) * expo;
  }

  DVec ExpFourierBasis::gradient(double variable, const int iOfElt) const
  {
    //cout << basisCoefficient_ << " " << exp(beta_*potential(variable)/2) << " " << potDeriv(variable) << endl;
    double expo = basisCoefficient_ * exp(beta_ * potential(variable) / 2);
    int iOfFreq = (iOfElt + 1) / 2;
    if (iOfElt == 0)
      return DVec::Constant(1, beta_ / 2 * potDeriv(variable) * expo);
    else if (iOfElt % 2 == 1)
      return DVec::Constant(1, sqrt(2.) * ( iOfFreq * cos(iOfFreq * variable) + sin(iOfFreq * variable) * beta_ / 2 * potDeriv(variable)) * expo);
    else
      return DVec::Constant(1, sqrt(2.) * (-iOfFreq * sin(iOfFreq * variable) + cos(iOfFreq * variable) * beta_ / 2 * potDeriv(variable)) * expo);
  }

  double ExpFourierBasis::laplacian(double variable, const int iOfElt) const
  {
    double expo = basisCoefficient_ * exp(potential(variable) / 2);
    int iOfFreq = (iOfElt + 1) / 2;
    if (iOfElt == 0)
      return beta_ / 2 * (potLapla(variable) + beta_ / 2 * pow(potDeriv(variable) , 2)) * expo;
    else if (iOfElt % 2 == 1)
      return sqrt(2.) * (beta_ * iOfFreq * potDeriv(variable) * cos(iOfFreq * variable)
                         + (beta_ / 2 * potLapla(variable) + pow(beta_, 2) / 4 * pow(potDeriv(variable), 2) - pow(iOfFreq, 2)) * sin(iOfFreq * variable)) * expo;
    else
      return sqrt(2.) * (- beta_ * iOfFreq * potDeriv(variable) * sin(iOfFreq * variable)
                         + (beta_ / 2 * potLapla(variable) + pow(beta_, 2) / 4 * pow(potDeriv(variable), 2) - pow(iOfFreq, 2)) * cos(iOfFreq * variable)) * expo;
  }
  
     // Compute <G_k, \partial_q G_l>
  double ExpFourierBasis::xGradY(const int iOfEltLeft, const int iOfEltRight) const
  {
    if ((iOfEltLeft - iOfEltRight) % 2 == 0 || iOfEltLeft < 0 || iOfEltRight < 0) return 0;
    double result = 0;
    if (iOfEltRight % 2 == 1)
    {
      if (iOfEltLeft == iOfEltRight-1) result = amplitude() * beta_ / 4;
      else if (iOfEltLeft == iOfEltRight+1) result = (iOfEltRight+1)/2;
      else if (iOfEltLeft == iOfEltRight+3) result = -amplitude() * beta_ / 4;
      else return 0;
    }
    else if (iOfEltRight % 2 == 0)
    {
      if (iOfEltLeft == iOfEltRight-3) result = -amplitude() * beta_ / 4;
      else if (iOfEltLeft == iOfEltRight-1) result = -iOfEltRight/2;
      else if (iOfEltLeft == iOfEltRight+1) result = amplitude() * beta_ / 4;
      else return 0;
    }
    
    if (iOfEltLeft == 0 || iOfEltRight == 0)
      return result * sqrt(2.);
    else return result;
  }
  
  void ExpFourierBasis::gradMatrix(SMat& A) const
  {
    for (int iOfFourier = 0; iOfFourier < nbOfElts_; iOfFourier++)
      for (int jOfFourier = max(0, iOfFourier-3); jOfFourier < min(nbOfElts_, iOfFourier+4); jOfFourier++)
        A.insert(iOfFourier, jOfFourier) = xGradY(iOfFourier, jOfFourier);
  }
  
  // Compute <G_k, \partial_q^* \partial_q G_l>
  double ExpFourierBasis::xLaplacianY(const int iOfEltLeft, const int iOfEltRight) const
  {
    if ((iOfEltLeft - iOfEltRight) % 2 == 1 || iOfEltLeft < 0 || iOfEltRight < 0) return 0;
    double result = 0;
    for (int iOfElt = max(iOfEltLeft, iOfEltRight) - 3; iOfElt < min(iOfEltLeft, iOfEltRight) + 4; iOfElt++)
    //for (int iOfElt = 0; iOfElt < nbOfElts(); iOfElt++)
      result += xGradY(iOfElt, iOfEltLeft) * xGradY(iOfElt, iOfEltRight);
    return result;
    
    /*if      (iOfEltLeft == iOfEltRight-4) result = - pow(amplitude() * beta_, 2) / 16;
    else if (iOfEltLeft == iOfEltRight-2) result = -amplitude() * beta_/4;
    else if (iOfEltLeft == iOfEltRight)   result = pow(floor((iOfEltRight+1)/2), 2) + pow(beta_, 2) / 8;      
    else if (iOfEltLeft == iOfEltRight+2) result = -amplitude() * beta_/4;
    else if (iOfEltLeft == iOfEltRight+4) result = - pow(amplitude() * beta_, 2) / 16;
    else result = 0;
    
    return result * (iOfEltLeft == 0?sqrt(2.):1) * (iOfEltRight == 0?sqrt(2.):1);*/
  }
  
  void ExpFourierBasis::laplacianMatrix(SMat& A) const
  {
    cout << "laplacianMatrix" << endl;
    for (int iOfFourier = 0; iOfFourier < nbOfElts_; iOfFourier++)
      for (int jOfFourier = max(0, iOfFourier-4); jOfFourier < min(nbOfElts_, iOfFourier+5); jOfFourier++)
      {
        cout << "coucou" << endl;
        A.insert(iOfFourier, jOfFourier) = xLaplacianY(iOfFourier, jOfFourier);
      }
  }
    
 
    
  
  
  

  HermiteBasis::HermiteBasis(const int nbOfElts0, double beta0)
    : Basis(nbOfElts0),
      beta_(beta0),
      polyCoeffs_(DMat::Zero(nbOfElts0, nbOfElts0))
      //polyCoeffs_(DMat::Zero(nbOfElts0, nbOfElts0))
  {
    //cout << "HermiteBasis(const int nbOfElts0)" << endl;
    polyCoeffs_(0, 0) = 1;
    polyCoeffs_(1, 1) = beta_;
    for (int iOfElt = 2; iOfElt < (int)nbOfElts0; iOfElt++)
      for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
      {
        polyCoeffs_(iOfElt, iOfCoeff) = (iOfCoeff ? (sqrt(beta_ / iOfElt) * polyCoeffs_(iOfElt - 1, iOfCoeff - 1)) : 0) - sqrt((iOfElt - 1) / (double)iOfElt) * polyCoeffs_(iOfElt - 2, iOfCoeff);
        //polyCoeffs_(iOfElt, iOfCoeff) /= sqrt(iOfElt);
      }
    //cout << "Hermite coeffs : " << endl;
    //cout << polyCoeffs_ << endl;
  }

  double HermiteBasis::value(double variable, const int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
    {
      result += polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff);
      //cout << polyCoeffs_(iOfElt, iOfCoeff) << endl;
    }
    //cout << variable << " , " << iOfElt << " -> " << result << endl;
    return result;
  }

  DVec HermiteBasis::gradient(double variable, const int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 1; iOfCoeff <= iOfElt; iOfCoeff++)
      result += iOfCoeff * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 1);
    return DVec::Constant(1, result);
  }

  double HermiteBasis::laplacian(double variable, const int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 2; iOfCoeff <= iOfElt; iOfCoeff++)
      result += iOfCoeff * (iOfCoeff - 1) * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 2);
    return result;
  }
  
  // Compute <H_n, \partial_p H_m>
  double HermiteBasis::xGradY(const int iOfEltLeft, const int iOfEltRight) const
  {
    return (iOfEltLeft == iOfEltRight+1) * sqrt(beta_*iOfEltRight);
  }
  
  void HermiteBasis::gradMatrix(SMat& A) const
  {
    for (int iOfHermite = 0; iOfHermite< nbOfElts_-1; iOfHermite++)
      A.insert(iOfHermite+1, iOfHermite) = xGradY(iOfHermite+1, iOfHermite);
  }
  
  // Compute <H_n, \partial_p H_m>
  double HermiteBasis::xLaplacianY(const int iOfEltLeft, const int iOfEltRight) const
  {
    return (iOfEltLeft == iOfEltRight) * beta_*iOfEltRight;
  }
 
  void HermiteBasis::laplacianMatrix(SMat& A) const
  {
    for (int iOfHermite = 0; iOfHermite< nbOfElts_; iOfHermite++)
      A.insert(iOfHermite, iOfHermite) = xLaplacianY(iOfHermite, iOfHermite);
  }
  
  
 

  TensorBasis::TensorBasis(const int nbOfVariables0):
    bases_(nbOfVariables0)
  {
  }

  TensorBasis::~TensorBasis()
  {
    for (int iOfVariable = 0; iOfVariable < (int)nbOfVariables(); iOfVariable++)
      delete bases_[iOfVariable];
  }

  int TensorBasis::nbOfVariables() const
  {
    return bases_.size();
  }

  int const& TensorBasis::nbOfElts(const int iOfVariable) const
  {
    return bases_[iOfVariable]->nbOfElts();
  }

  int& TensorBasis::nbOfElts(const int iOfVariable)
  {
    return bases_[iOfVariable]->nbOfElts();
  }

  vector<int> TensorBasis::nbOfElts() const
  {
    vector<int> nbOfElts0(nbOfVariables());
    for (int iOfVariable = 0; iOfVariable < (int)nbOfVariables(); iOfVariable++)
      nbOfElts0[iOfVariable] = bases_[iOfVariable]->nbOfElts();
    return nbOfElts0;
  }
  
  int TensorBasis::totalNbOfElts() const
  {
    int totalNbOfElts0 = 1;
    for (int iOfVariable = 0; iOfVariable < (int)nbOfVariables(); iOfVariable++)
      totalNbOfElts0 *= bases_[iOfVariable]->nbOfElts();
    return totalNbOfElts0;
  }

  const Basis* TensorBasis::operator()(const int iOfBasis) const
  {
    return bases_[iOfBasis];
  }
  
  

  QPBasis::QPBasis():
    TensorBasis(2)
  {
    //cout << "QPBasis()" << endl;
  }

  /*QPBasis::QPBasis(const int nbOfFourier0, const int nbOfHermite0):
    TensorBasis(2)
  {
    bases_[0] = new ExpFourierBasis(nbOfFourier0);
    bases_[1] = new HermiteBasis(nbOfHermite0);
  }*/

  int& QPBasis::nbOfFourier()
  {
    return nbOfElts(0);
  }

  const int& QPBasis::nbOfFourier() const
  {
    return nbOfElts(0);
  }

  int& QPBasis::nbOfHermite()
  {
    return nbOfElts(1);
  }

  const int& QPBasis::nbOfHermite() const
  {
    return nbOfElts(1);
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  int QPBasis::iTens(int iOfFourier2, int iOfHermite) const
  {
    assert(iOfFourier2 < nbOfFourier()  && iOfHermite < nbOfHermite());
    return nbOfFourier() * iOfHermite + iOfFourier2;
  }

  vector<int> QPBasis::vecTens(int iTens0) const
  {
    assert(iTens0 < nbOfHermite() * nbOfFourier());
    vector<int> vecIndex(2);
    vecIndex[0] = iTens0 % nbOfFourier();
    vecIndex[1] = iTens0 / nbOfFourier();
    return vecIndex;
  }
  
  /*DVec QPBasis::values(System const& syst) const
  {
    DVec val(nbOfElts(), 0);
    for (int iOfElt=0; iOfElt<nbOfElts(); iOfElt++)
      val(iOfElt) = values
  }*/

  double QPBasis::value(System const& syst, const int iOfElt) const
  {
    vector<int> vecIndex = vecTens(iOfElt);
    return value(syst, vecIndex);
  }

  double QPBasis::value(System const& syst, vector<int>& vecIndex) const
  {
    //cout << syst(0).position(0) << " " << bases_[0]->value(syst(0).position(0), vecIndex[0]) << " " << bases_[1]->value(syst(0).momentum(0), vecIndex[1]) << endl;
    return bases_[0]->value(syst(0).position(0), vecIndex[0]) * bases_[1]->value(syst(0).momentum(0), vecIndex[1]);
  }

  DVec QPBasis::gradientQ(System const& syst, int /*iOfParticle*/, vector<int>& vecIndex) const
  {
    //cout << "QPBasis::gradientQ" << endl;
    //cout << vecIndex[0] << " " << vecIndex[1] << endl;
    //cout << bases_[0]->gradient(syst(0).position(0), vecIndex[0]) << " X "
    //    << bases_[1]->value(syst(0).momentum(0), vecIndex[1]) << endl;
    return bases_[0]->gradient(syst(0).position(0), vecIndex[0])
           * bases_[1]->value(syst(0).momentum(0), vecIndex[1]);
  }

  DVec QPBasis::gradientQ(System const& syst, int iOfParticle, int iOfCoeff) const
  {
    vector<int> vecIndex = vecTens(iOfCoeff);
    return gradientQ(syst, iOfParticle, vecIndex);
  }

  double QPBasis::laplacianQ(System const& syst, int /*iOfParticle*/, vector<int>& vecIndex) const
  {
    return bases_[0]->laplacian(syst(0).position(0), vecIndex[0])
           * bases_[1]->value(syst(0).momentum(0), vecIndex[1]);
  }

  double QPBasis::laplacianQ(System const& syst, int iOfParticle, int iOfCoeff) const
  {
    vector<int> vecIndex = vecTens(iOfCoeff);
    return laplacianQ(syst, iOfParticle, vecIndex);
  }

  DVec QPBasis::gradientP(System const& syst, int /*iOfParticle*/, vector<int>& vecIndex) const
  {
    //cout << bases_[0]->value(syst(0).position(0), vecIndex[0]) << " X "
    //    << bases_[1]->gradient(syst(0).momentum(0), vecIndex[1]) << endl;
    return bases_[0]->value(syst(0).position(0), vecIndex[0])
           * bases_[1]->gradient(syst(0).momentum(0), vecIndex[1]);
  }

  DVec QPBasis::gradientP(System const& syst, int iOfParticle, int iOfCoeff) const
  {
    vector<int> vecIndex = vecTens(iOfCoeff);
    return gradientP(syst, iOfParticle, vecIndex);
  }

  double QPBasis::laplacianP(System const& syst, int /*iOfParticle*/, vector<int>& vecIndex) const
  {
    return bases_[0]->value(syst(0).position(0), vecIndex[0])
           * bases_[1]->laplacian(syst(0).momentum(0), vecIndex[1]);
  }

  double QPBasis::laplacianP(System const& syst, int iOfParticle, int iOfCoeff) const
  {
    vector<int> vecIndex = vecTens(iOfCoeff);
    return laplacianP(syst, iOfParticle, vecIndex);
  }

  
  

  /*FourierHermiteBasis::FourierHermiteBasis(Input const& input):
    QPBasis()
  {
    bases_[0] = new FourierBasis(input.nbOfFourier());
    bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
  }*/

  ExpFourierHermiteBasis::ExpFourierHermiteBasis(Input const& input, Potential& potential):
    QPBasis()
  {
    //cout << "ExpFourierHermiteBasis(Input const& input)" << endl;
    bases_[0] = new ExpFourierBasis(input.nbOfFourier(), input.beta(), potential);
    bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
    //cout << "end ExpFourierHermiteBasis(Input const& input)" << endl;
  }

  double const& ExpFourierHermiteBasis::expFourierMeans(int iOfElt) const
  {
    return bases_[0]->expFourierMeans(iOfElt);
  }
  
  DVec const& ExpFourierHermiteBasis::gVector() const
  {
    return bases_[0]->gVector();
  }
  
  double const& ExpFourierHermiteBasis::norm2gVector() const
  {
    return bases_[0]->norm2gVector();
  }
    
}












