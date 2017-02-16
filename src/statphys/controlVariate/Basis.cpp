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

  double DTVec::operator()(vector<int>& vecIndex) const
  {
    return data_(iTens(vecIndex));
  }

  double& DTVec::operator()(vector<int>& vecIndex)
  {
    return data_(iTens(vecIndex));
  }

  double DTVec::operator()(int iTensOfElt) const
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


  
  

  Basis::Basis(const int nbOfElts0):
    nbOfElts_(nbOfElts0),
    basisMeans2_(DVec::Zero(nbOfElts0)),
    norm2meansVec_(0)
  {}

  int const& Basis::nbOfElts() const
  {
    return nbOfElts_;
  }

  int& Basis::nbOfElts()
  {
    return nbOfElts_;
  }

  DVec const& Basis::basisMeans2() const
  {
    return basisMeans2_;
    //throw std::runtime_error("Basis::basisMean2 is not implemented");
  };
  
  double Basis::basisMean2(int iOfElt) const
  {
    return basisMeans2_[iOfElt];
    //throw std::runtime_error("Basis::basisMean2 is not implemented");
  };
  
  /*DVec const& Basis::gVector() const
  {
    throw std::runtime_error("Basis::gVector is not implemented");
  }
  
  double Basis::gVector(int) const
  {
    throw std::runtime_error("Basis::gVector(int) is not implemented");
  }*/
  
  double Basis::norm2meansVec() const
  {
    throw std::runtime_error("Basis::norm2meansVec is not implemented");
  }
  
  
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
  
  // ##### QBasis ####
  
  QBasis::QBasis(const int nbOfElts0, double beta0, Potential& potential0):
    Basis(nbOfElts0),
    beta_(beta0),
    potential_(&potential0),
    nbOfIntegrationNodes_(1000*nbOfElts()),
    qRepartitionFct_(0)
    //basisMeans2_(DVec::Zero(nbOfElts0)),

  {}
  
  double QBasis::potential(double variable) const
  {
    return potential_->value(variable);
  }

  double QBasis::potDeriv(double variable) const
  {
    return (potential_->gradient(variable))(0);
  }

  double QBasis::potLapla(double variable) const
  {
    return (potential_->laplacian(variable));
  }
  
  const double& QBasis::amplitude() const
  {
    return potential_->parameter1();
  }
  

  /*double QBasis::basisMean2(int iOfElt) const
  {
    return basisMean2_(iOfElt);
  }
  
  DVec const& QBasis::gVector() const
  {
    return basisMeans2_;
  }
  
  double QBasis::gVector(int iOfElt) const
  {
    return basisMeans2_[iOfElt];
  }*/
  
  /*double QBasis::norm2meansVec() const
  {
    return norm2meansVec_;
  }*/
  
  // #### ExpFourierBasis ####
  
  ExpFourierBasis::ExpFourierBasis(const int nbOfElts0, double beta0, Potential& potential0):
    QBasis(nbOfElts0, beta0, potential0),
    expFourierCoeffs_(DVec::Zero(2 * nbOfElts0))
  {
    computeBasisMeans();
  }
  
  DVec const& ExpFourierBasis::expFourierCoeffs() const
  {
    return expFourierCoeffs_;
  }
  
  double ExpFourierBasis::expFourierCoeff(int iOfElt) const
  {
    return expFourierCoeffs_[iOfElt];
  }

  
  ///We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  ///More precisely a_k = 1 / pi int cos(k q) C^-1 exp(-beta V / 2) dq and b_k = 1 / pi int sin(k q) C^-1 exp(-beta V / 2) dq
  ///where C = ExpFourierBasis::basisCoefficient_ Note that a_0 has the same normalization as a_k so that a_{k-l} = a_0 when k=l
  void ExpFourierBasis::computeBasisMeans()
  {
    double step = 2 * M_PI / (double)nbOfIntegrationNodes_;
    for (int iOfNode = 0; iOfNode < nbOfIntegrationNodes_; iOfNode++)
    {
      double q = - M_PI + iOfNode * step;
      qRepartitionFct_ += exp(-beta_ * potential(q)) * step;
      for (int iOfFourier = 0; iOfFourier < 2 * nbOfElts(); iOfFourier++)
        if (iOfFourier % 2 == 0) expFourierCoeffs_(iOfFourier) += cos(iOfFourier/2 * q) * exp(-beta_ * potential(q) / 2);
        else expFourierCoeffs_(iOfFourier) += sin((iOfFourier+1)/2 * q) * exp(-beta_ * potential(q) / 2);
    }
    basisCoefficient_ = sqrt(qRepartitionFct_ / (2 * M_PI));
    
    expFourierCoeffs_ *= step / (M_PI * basisCoefficient_);
    
    // gVector is the projection of the constant fct 1, ie it contains the mean of the expfourier elements
    basisMeans2_ = expFourierCoeffs_.head(nbOfElts()) / sqrt(2.);
    basisMeans2_(0) /= sqrt(2.);
    norm2meansVec_ = pow(basisMeans2_.norm(), 2);
    cout << "Squared norm of vector g = " << norm2meansVec_ << endl;
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
  
    //We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  //where C = ExpFourierBasis::basisCoefficient_
  void ExpFourierBasis::computeExpToTrigMat()
  {
    for (int iOfFourier2 = 1; iOfFourier2 < nbOfElts(); iOfFourier2++)
    {
      trigToExpMat_(0, iOfFourier2) = expFourierCoeff(iOfFourier2) / sqrt(2.);
      trigToExpMat_(iOfFourier2, 0) = expFourierCoeff(iOfFourier2) / sqrt(2.);
    }
    trigToExpMat_(0, 0) = expFourierCoeff(0) / 2;

    int maxOfFourier = (nbOfElts() + 1) / 2;
    for (int iOfFourier = 1; iOfFourier < maxOfFourier; iOfFourier++)
      for (int jOfFourier = 1; jOfFourier < maxOfFourier; jOfFourier++)
      {
        //cosine times cosine
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier) =
          expFourierCoeff(2 * (iOfFourier + jOfFourier)) / 2
          + expFourierCoeff(2 * abs(iOfFourier - jOfFourier)) / 2;

        //sinus times sinus
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier - 1) =
          - expFourierCoeff(2 * (iOfFourier + jOfFourier)) / 2
          + expFourierCoeff(2 * abs(iOfFourier - jOfFourier)) / 2;

        int eps = (iOfFourier >= jOfFourier) - (iOfFourier <= jOfFourier); // -1 if i smaller, 0 if equal, 1 if i larger
        //cosine times sinus   /!\ doing the test in this order avoids to read expFourierCoeff(-1)
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier - 1) =
          expFourierCoeff(2 * (iOfFourier + jOfFourier) - 1) / 2
          + (!eps)?0:eps * expFourierCoeff(2 * abs(iOfFourier - jOfFourier) - 1) / 2;
          
        //sinus times cosine
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier) =
          expFourierCoeff(2 * (iOfFourier + jOfFourier) - 1) / 2
          - (!eps)?0:eps * expFourierCoeff(2 * abs(iOfFourier - jOfFourier) - 1) / 2;
      }


    ofstream out_trigToExpMat("output/Galerkin/trigToExpMat");
    //display(trigToExpMat_, out_trigToExpMat);

    expToTrigMat_  = trigToExpMat_.inverse();
  }
  
  
// #### HermiteBasis ####

  HermiteBasis::HermiteBasis(const int nbOfElts0, double beta0)
    : Basis(nbOfElts0),
      beta_(beta0),
      polyCoeffs_(DMat::Zero(nbOfElts0, nbOfElts0))
  {
    polyCoeffs_(0, 0) = 1;
    polyCoeffs_(1, 1) = beta_;
    for (int iOfElt = 2; iOfElt < (int)nbOfElts0; iOfElt++)
      for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
        polyCoeffs_(iOfElt, iOfCoeff) = (iOfCoeff ? (sqrt(beta_ / iOfElt) * polyCoeffs_(iOfElt - 1, iOfCoeff - 1)) : 0) - sqrt((iOfElt - 1) / (double)iOfElt) * polyCoeffs_(iOfElt - 2, iOfCoeff);
  
    basisMeans2_[0] = 1;
  }

  double HermiteBasis::value(double variable, const int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
      result += polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff);
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
  
  double TensorBasis::basisMean2(int iOfVariable, int iOfElt) const
  {
    return bases_[iOfVariable]->basisMean2(iOfElt);
  }
  
  DVec const& TensorBasis::basisMeans2(int iOfVariable) const
  {
    return bases_[iOfVariable]->basisMeans2();
  }
  
  double TensorBasis::norm2meansVec(int iOfVariable) const
    {
    return bases_[iOfVariable]->norm2meansVec();
  }

  const Basis* TensorBasis::operator()(int iOfBasis) const
  {
    return bases_[iOfBasis];
  }
  
  Basis* TensorBasis::operator()(int iOfBasis)
  {
    return bases_[iOfBasis];
  }
  
  
  // #### QPBasis ####

  QPBasis::QPBasis():
    TensorBasis(2)
  {
    //cout << "QPBasis()" << endl;
  }

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
  
  DVec QPBasis::getBasisElement(int iOfElt1, int iOfElt2) const
  {
    cout << "getBasisElement" << endl;
    DVec vect = DVec::Zero(totalNbOfElts());
    vect(iTens(iOfElt1, iOfElt2)) = 1;
    cout << "end of getBasisElement" << endl;
    cout << vect << endl;
    return vect;
  }
  
  DVec QPBasis::getPartialElement(int iOfVariable, int iOfElt) const
  {
    cout << "getPartialElement" << endl;
    DVec vect = DVec::Zero(totalNbOfElts());
    if (iOfVariable == 0)
      for (int iOfElt2 = 0; iOfElt2 < nbOfElts(1); iOfElt2++)
        vect(iTens(iOfElt, iOfElt2)) = basisMean2(1, iOfElt2);
    else
      for (int iOfElt2 = 0; iOfElt2 < nbOfElts(0); iOfElt2++)
        vect(iTens(iOfElt2, iOfElt)) = basisMean2(0, iOfElt2);
    cout << "end of getPartialElement" << endl;
    cout << vect << endl;
    return vect;
  }

  // #### ExpFourierHermiteBasis ####

  ExpFourierHermiteBasis::ExpFourierHermiteBasis(Input const& input, Potential& potential):
    QPBasis()
  {
    bases_[0] = new ExpFourierBasis(input.nbOfFourier(), input.beta(), potential);
    bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
  }
  


  /*double ExpFourierHermiteBasis::basisMean2(int iOfVariable, int iOfElt) const
  {
    return bases_[iOfVariable]->basisMeans(iOfElt);
  }
  
  DVec const& ExpFourierHermiteBasis::gVector() const
  {
    return bases_[0]->gVector();
  }
  
  double ExpFourierHermiteBasis::gVector(int iOfElt) const
  {
    return bases_[0]->gVector(iOfElt);
  }
  
  double ExpFourierHermiteBasis::norm2meansVec() const
  {
    return bases_[0]->norm2meansVec();
  }*/
  
  DMat ExpFourierHermiteBasis::convertToTrigBasis(const DMat& X)
  {
    return bases_[0]->expToTrigMat() * X;
  }
  
  // #### HermiteHermiteBasis ####
    
  HermiteHermiteBasis::HermiteHermiteBasis(Input const& input, Potential& /*potential*/):
    QPBasis()
  {
    bases_[0] = new HermiteBasis(input.nbOfFourier(), input.beta()/*, potential*/);
    bases_[1] = new HermiteBasis(input.nbOfHermite(), input.beta());
  }
}












