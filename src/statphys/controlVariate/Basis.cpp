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


  
  

  Basis::Basis(int nbOfElts0, double beta0):
    nbOfElts_(nbOfElts0),
    basisMeans2_(DVec::Zero(nbOfElts0)),
    beta_(beta0),
    gramMatrix_(nbOfElts_, nbOfElts0),
    gradMatrix_(nbOfElts_, nbOfElts0),
    laplacianMatrix_(nbOfElts_, nbOfElts0)
  {
    if (nbOfElts0 == 0) throw runtime_error("The Galerkin basis is empty, check the input file !");
  }

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
  
  double Basis::norm2meansVec() const
  {
    return pow(basisMeans2().norm(), 2);
    //throw std::runtime_error("Basis::norm2meansVec is not implemented");
  }
  
  void Basis::computeGramMatrix()
  {
    for (int iOfElt1 = 0; iOfElt1< nbOfElts_; iOfElt1++)
      for (int iOfElt2 = iOfElt1; iOfElt2< nbOfElts_; iOfElt2++)
      {
        double scalarProd = xY(iOfElt1, iOfElt2);
        gramMatrix_.insert(iOfElt1, iOfElt2) = scalarProd;
        if (iOfElt1 != iOfElt2) gramMatrix_.insert(iOfElt2, iOfElt1) = scalarProd;
      }
  }
  
  void Basis::computeGradMatrix()
  {
    for (int iOfElt1 = 0; iOfElt1< nbOfElts_; iOfElt1++)
      for (int iOfElt2 = 0; iOfElt2< nbOfElts_; iOfElt2++)
        gradMatrix_.insert(iOfElt1, iOfElt2) = xGradY(iOfElt1, iOfElt2);
  }
  
  void Basis::computeLaplacianMatrix()
  {
    for (int iOfElt1 = 0; iOfElt1< nbOfElts_; iOfElt1++)
      for (int iOfElt2 = iOfElt1; iOfElt2< nbOfElts_; iOfElt2++)
      {
        double scalarProd = xLaplacianY(iOfElt1, iOfElt2);
        laplacianMatrix_.insert(iOfElt1, iOfElt2) = scalarProd;
        if (iOfElt1 != iOfElt2) laplacianMatrix_.insert(iOfElt2, iOfElt1) = scalarProd;
      }
  }
  
  /*// Compute <H_n, \partial_p^* \partial_p H_m>
  double Basis::xGradStarY(int iOfEltLeft, int iOfEltRight) const
  {
    return xGradY(iOfEltRight, iOfEltLeft);
  }
  
  void Basis::gradStarMatrix(SMat& A) const
  {
    gradMatrix(A);
    A = A.adjoint();
  }*/
  
  // ##### QBasis ####
  
  QBasis::QBasis(Input const& input, Potential& potential0):
    Basis(input.nbOfQModes(), input.beta()),
    potential_(&potential0),
    integrationStep_(input.integrationStep()),
    qRepartitionFct_(0),
    basisCoefficient_(0),
    measureMomenta_(DVec::Zero(2*nbOfElts_-1))
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
  
  
  int QBasis::nbOfIntegrationSteps() const
  {
    //cout << "length = " << length() << " integrationStep_ = " << integrationStep_ << "nbOfIntegrationSteps = " << (int) (length() / integrationStep_) << endl;
    return (int) (length() / integrationStep_);
  }
  
  // #### ExpFourierBasis ####
  
  ExpFourierBasis::ExpFourierBasis(Input const& input, Potential& potential0):
    QBasis(input, potential0),
    expFourierCoeffs_(DVec::Zero(2 * nbOfElts_))
  {
    /*cout << "nb :" << nbOfElts0 << endl;
    cout << "Constructor expFourierCoeffs_ : " << endl << expFourierCoeffs_ << endl << endl;
    cout << "zero :" << DVec::Zero(2 * nbOfElts0)<< endl << endl;*/
    computeBasisMeans();
    computeGramMatrix();
    computeGradMatrix();
    computeLaplacianMatrix();  
  }
  
  double ExpFourierBasis::length() const
  {
    return 2*M_PI;
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
    double step = integrationStep_;
    //double step = 2 * M_PI / (double)nbOfIntegrationNodes_;
    //nbOfIntegrationNodes = (int)
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationSteps(); iOfNode++)
    {
      double q = - M_PI + iOfNode * step;
      qRepartitionFct_ += exp(-beta_ * potential(q)) * step;
      for (int iOfFourier = 0; iOfFourier < 2 * nbOfElts(); iOfFourier++)
        if (iOfFourier % 2 == 0) expFourierCoeffs_(iOfFourier) += cos(iOfFourier/2 * q) * exp(-beta_ * potential(q) / 2);
        else expFourierCoeffs_(iOfFourier) += sin((iOfFourier+1)/2 * q) * exp(-beta_ * potential(q) / 2);
    }
    basisCoefficient_ = sqrt(qRepartitionFct_ / (2 * M_PI));
    
    expFourierCoeffs_ *= step / (M_PI * basisCoefficient_);
    
    cout << "expFourierCoeffs_ :" << endl << expFourierCoeffs_ << endl;
    
    // gVector is the projection of the constant fct 1, ie it contains the mean of the expfourier elements
    basisMeans2_ = expFourierCoeffs_.head(nbOfElts()) / sqrt(2.);
    basisMeans2_(0) /= sqrt(2.);
    //norm2meansVec_ = pow(basisMeans2_.norm(), 2);
    //cout << "Squared norm of vector g = " << norm2meansVec_ << endl;
  }

  double ExpFourierBasis::value(double variable, int iOfElt) const
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

  DVec ExpFourierBasis::gradient(double variable, int iOfElt) const
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

  double ExpFourierBasis::laplacian(double variable, int iOfElt) const
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
  
  // Compute <G_k, G_k'>
  double ExpFourierBasis::xY(int iOfEltLeft, int iOfEltRight) const
  {
    return (iOfEltLeft == iOfEltRight);
  }
  
  /*void ExpFourierBasis::gramMatrix(SMat& A) const
  {
    A.setIdentity();
  }*/
  
  // Compute <G_k, \partial_q G_l>
  double ExpFourierBasis::xGradY(int iOfEltLeft, int iOfEltRight) const
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
  
  // Compute <G_k, \partial_q^* \partial_q G_l>
  double ExpFourierBasis::xLaplacianY(int iOfEltLeft, int iOfEltRight) const
  {
    if ((iOfEltLeft - iOfEltRight) % 2 == 1 || iOfEltLeft < 0 || iOfEltRight < 0) return 0;
    double result = 0;
    for (int iOfElt = max(iOfEltLeft, iOfEltRight) - 3; iOfElt < min(iOfEltLeft, iOfEltRight) + 4; iOfElt++)
    //for (int iOfElt = 0; iOfElt < nbOfElts(); iOfElt++)
      result += xGradY(iOfElt, iOfEltLeft) * xGradY(iOfElt, iOfEltRight);
    return result;
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
  
  // #### HermiteQBasis ####

  ///
  ///Builds the basis composed of the Hermite polynomials
  ///The measure is any nu so that the basis is not orthonormal
  HermiteQBasis::HermiteQBasis(Input const& input, Potential& potential0)
    : QBasis(input, potential0),
      length_(input.integrationLength()),
      qMin_(input.integrationQMin()),
      //qMin_(0),
      omega_(input.omegaHermite()),
      polyCoeffs_(DMat::Zero(nbOfElts_, nbOfElts_))
  {    
    if (length_ == 0)
      throw runtime_error("The integration length is set to zero !");
    cout << "HermiteQBasis::HermiteQBasis" << endl;
    polyCoeffs_(0, 0) = 1;
    polyCoeffs_(1, 1) = sqrt(omega_);
    for (int iOfElt = 2; iOfElt < (int)nbOfElts_; iOfElt++)
      for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
        polyCoeffs_(iOfElt, iOfCoeff) = (iOfCoeff ? (sqrt(omega_ / iOfElt) * polyCoeffs_(iOfElt - 1, iOfCoeff - 1)) : 0) - sqrt((iOfElt - 1) / (double)iOfElt) * polyCoeffs_(iOfElt - 2, iOfCoeff);
  
    cout << "Initialization of HermiteQBasis..."; cout.flush();
    computeMeasureMomenta();
    computeGramMatrix();
    computeGradMatrix();
    computeLaplacianMatrix();  
    computeBasisMeans(); 
    cout << "OK !" << endl;
    
    //double step = integrationStep_;
    //double step = length_ / (double)nbOfIntegrationSteps();
    ofstream qbasis("output/Galerkin/qbasis", std::ofstream::app);
    ofstream qbasisgrad("output/Galerkin/qbasisgrad", std::ofstream::app);
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationSteps(); iOfNode++)
    {
      double q = qMin_ + iOfNode * integrationStep_;
      qbasis << q << " " << value(q, 0) << " " << value(q, 1) << " " << value(q, 2) <<" " <<exp(-beta_ * potential(q)) << " " << potential(q) << endl; 
      qbasisgrad << q << " " << gradient(q, 0)(0) << " " << gradient(q, 1)(0) << " " << gradient(q, 2)(0) <<" " <<exp(-beta_ * potential(q)) << " " << potential(q) << endl; 
   }
  }
  
  double HermiteQBasis::length() const
  {
    return length_;
  }
  
  void HermiteQBasis::computeMeasureMomenta()
  {
    //double step = (length_/2) / (double)nbOfIntegrationNodes_;
    // Triangle integration
    measureMomenta_(0) += exp(-beta_ * potential(0));
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationSteps()/2; iOfNode++)
    {
      double q = (iOfNode+1) * integrationStep_;
      double betaPotPos = beta_ * potential(q);
      double betaPotNeg = beta_ * potential(-q);
      //double logq = (q!=0)?log(pow(q,2))/2:-1000;
      double logq = log(q);
      //cout << "q=" << q << " logq = " << logq << " betaPotPos=" << betaPotPos << " betaPotNeg=" << betaPotNeg<< endl;
      for (int iOfElt = 0; iOfElt < 2*nbOfElts()-1; iOfElt++)
      {
        if (iOfElt%2==0)
          measureMomenta_(iOfElt) += exp(iOfElt * logq - betaPotNeg) + exp(iOfElt * logq - betaPotPos);
        else
          measureMomenta_(iOfElt) += -exp(iOfElt * logq - betaPotNeg) + exp(iOfElt * logq - betaPotPos);
        //measureMomenta_(iOfElt) += ((q<0 && iOfElt%2==1)?-1:1) * exp(iOfElt * logq - betaPot);
        //if (iOfElt == 2*nbOfElts()-2)
        //  cout << "q=" << q << " betaPot=" << betaPot << " arg = " << ((q<0 && iOfElt%2==1)?-1:1) * exp(iOfElt * logq - betaPot) << endl;
        //if (iOfElt == 2)
          //cout << "q=" << q << " exp(iOfElt * logq - betaPotNeg) = " << exp(iOfElt * logq - betaPotPos) << endl;
          //cout << "q=" << q << " A " << exp(iOfElt * logq - betaPotNeg) << " B " << measureMomenta_(iOfElt)*step << endl;
      }
    }
    measureMomenta_ *= integrationStep_;
    qRepartitionFct_ = measureMomenta_[0];
    measureMomenta_ /= qRepartitionFct_;
    cout << "end computeMeasureMomenta :" << endl << measureMomenta_ << endl;
    cout << "qRepartitionFct_ = " << qRepartitionFct_ << endl;
  }
  
  void HermiteQBasis::computeBasisMeans()
  {   
    for (int iOfElt = 0; iOfElt < nbOfElts(); iOfElt++)
      basisMeans2_[iOfElt] = dot(polyCoeffs_.row(iOfElt), measureMomenta_.head(nbOfElts()+1));
    
    /*cout << "basisMeans2_1 :" << endl << basisMeans2_ << endl;
    
    double step = length_ / (double)nbOfIntegrationNodes_;
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationNodes_; iOfNode++)
    {
      double q = qMin_ + iOfNode * step;
      qRepartitionFct_ += exp(-beta_ * potential(q)) * step;
      for (int iOfElt = 0; iOfElt < nbOfElts(); iOfElt++)
        basisMeans2_(iOfElt) += value(q, iOfElt) * exp(-beta_ * potential(q)) * step;
    }
    basisMeans2_ /= qRepartitionFct_;
    
    cout << "basisMeans2_2 :" << endl << basisMeans2_ << endl;
    //cout << "Squared norm of vector g = " << norm2meansVec_ << endl;*/
  }

  double HermiteQBasis::value(double variable, int iOfElt) const
  {
    if (iOfElt < 0) return 0;
    double result = 0;
    for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
      result += polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff);
    return result;
  }

  DVec HermiteQBasis::gradient(double variable, int iOfElt) const
  {
    double result = 0;
    //for (int iOfCoeff = 1; iOfCoeff <= iOfElt; iOfCoeff++)
    //  result += iOfCoeff * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 1);
    if (iOfElt == 0) result = 0;
    else result = sqrt(omega_ * iOfElt) * value(variable, iOfElt-1);
    return DVec::Constant(1, result);
  }

  double HermiteQBasis::laplacian(double variable, int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 2; iOfCoeff <= iOfElt; iOfCoeff++)
      result += iOfCoeff * (iOfCoeff - 1) * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 2);
    return result;
  }
  
  // Computes <H_n, \partial_p H_m> for any the measure nu on [0, length]
  double HermiteQBasis::xY(int iOfEltLeft, int iOfEltRight) const
  {
    DVec hermiteProduct = polynomialProduct(polyCoeffs_.row(iOfEltLeft), polyCoeffs_.row(iOfEltRight));
    /*cout << iOfEltLeft << " x " << iOfEltRight  << endl;
    cout << "X =" << endl << polyCoeffs_.row(iOfEltLeft) << endl;
    cout << "Y =" << endl << polyCoeffs_.row(iOfEltRight) << endl;
    cout << "hermiteProduct = " << endl << hermiteProduct << endl;
    cout << "measureMomenta_ = " << endl << measureMomenta_ << endl;*/
    return dot(hermiteProduct, measureMomenta_);
    
    /*double scalarProd = 0;
    int nbOfIntegrationNodes = nbOfIntegrationNodes_;
    double step = length_ / (double)nbOfIntegrationNodes;
    if(exp(-beta_ * potential(-qMin_)) > 1e-12 || exp(-beta_ * potential(qMin_)) > 1e-12)
      throw runtime_error("Integration length too small !");
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationNodes; iOfNode++)
    {
      double q = qMin_ + iOfNode * step;
      double expo = exp(-beta_ * potential(q));
      if (expo > 1e-15)
        scalarProd += value(q, iOfEltLeft) * value(q, iOfEltRight) * expo / qRepartitionFct_ * step;
    }
    //cout << iOfEltLeft << " xY " << iOfEltRight << " -> " << scalarProd << endl;
    return 1e-9*round(1e9*scalarProd);*/
  }
  
  // Compute <H_n, \partial_q H_m>
  double HermiteQBasis::xGradY(int iOfEltLeft, int iOfEltRight) const
  {
    cout << "xGradY should not be called !" << endl;
    /*double scalarProd = 0;
    double step = length_ / (double)nbOfIntegrationNodes_;
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationNodes_; iOfNode++)
    {
      double q = qMin_ + iOfNode * step;
      scalarProd += value(q, iOfEltLeft) * gradient(q, iOfEltRight)(0) * exp(-beta_ * potential(q)) / qRepartitionFct_ * step;
    }
    return 1e-9*round(1e9*scalarProd);*/
    return sqrt(omega_ * iOfEltRight) * xY(iOfEltLeft, iOfEltRight-1);
  }
  
  // Compute <partial_q H_n, \partial_q H_m>
  double HermiteQBasis::xLaplacianY(int iOfEltLeft, int iOfEltRight) const
  {
    cout << "xLaplacianY should not be called !" << endl;
    /*double scalarProd = 0;
    double step = length_ / (double)nbOfIntegrationNodes_;
    for (int iOfNode = 0; iOfNode <= nbOfIntegrationNodes_; iOfNode++)
    {
      double q = qMin_ + iOfNode * step;
      scalarProd += gradient(q, iOfEltLeft)(0) * gradient(q, iOfEltRight)(0) * exp(-beta_ * potential(q)) / qRepartitionFct_ * step;
    }
    return 1e-9*round(1e9*scalarProd);*/
    return omega_ * sqrt(iOfEltLeft * iOfEltRight) * xY(iOfEltLeft-1, iOfEltRight-1);
  }  
  
  void HermiteQBasis::computeGradMatrix()
  {
    for (int iOfElt1 = 0; iOfElt1< nbOfElts_; iOfElt1++)
      for (int iOfElt2 = 1; iOfElt2< nbOfElts_; iOfElt2++)
        gradMatrix_.insert(iOfElt1, iOfElt2) = sqrt(omega_ * iOfElt2) * gramMatrix_.coeff(iOfElt1, iOfElt2-1);
  }
  
  void HermiteQBasis::computeLaplacianMatrix()
  {
    for (int iOfElt1 = 1; iOfElt1< nbOfElts_; iOfElt1++)
      for (int iOfElt2 = 1; iOfElt2< nbOfElts_; iOfElt2++)
        laplacianMatrix_.insert(iOfElt1, iOfElt2) = omega_ * sqrt(iOfElt1 * iOfElt2) * gramMatrix_.coeff(iOfElt1-1, iOfElt2-1);
  }
  
// #### HermiteBasis ####

  HermiteBasis::HermiteBasis(Input const& input)
    : Basis(input.nbOfPModes(), input.beta()),
      polyCoeffs_(DMat::Zero(nbOfElts_, nbOfElts_))
  {
    polyCoeffs_(0, 0) = 1;
    polyCoeffs_(1, 1) = sqrt(beta_);
    for (int iOfElt = 2; iOfElt < (int)nbOfElts_; iOfElt++)
      for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
        polyCoeffs_(iOfElt, iOfCoeff) = (iOfCoeff ? (sqrt(beta_ / iOfElt) * polyCoeffs_(iOfElt - 1, iOfCoeff - 1)) : 0) - sqrt((iOfElt - 1) / (double)iOfElt) * polyCoeffs_(iOfElt - 2, iOfCoeff);
  
    basisMeans2_[0] = 1;
    computeGramMatrix();
    computeGradMatrix();
    computeLaplacianMatrix();  
  }

  double HermiteBasis::value(double variable, int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 0; iOfCoeff <= iOfElt; iOfCoeff++)
      result += polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff);
    return result;
  }

  DVec HermiteBasis::gradient(double variable, int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 1; iOfCoeff <= iOfElt; iOfCoeff++)
      result += iOfCoeff * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 1);
    return DVec::Constant(1, result);
  }

  double HermiteBasis::laplacian(double variable, int iOfElt) const
  {
    double result = 0;
    for (int iOfCoeff = 2; iOfCoeff <= iOfElt; iOfCoeff++)
      result += iOfCoeff * (iOfCoeff - 1) * polyCoeffs_(iOfElt, iOfCoeff) * pow(variable, iOfCoeff - 2);
    return result;
  }
  
  // Computes <H_k, H_k'> for the measure kappa
  double HermiteBasis::xY(int iOfEltLeft, int iOfEltRight) const
  {
    return (iOfEltLeft == iOfEltRight);
  }
  
  // Computes <H_n, \partial_p H_m>
  double HermiteBasis::xGradY(int iOfEltLeft, int iOfEltRight) const
  {
    return (iOfEltLeft+1 == iOfEltRight) * sqrt(beta_*iOfEltRight);
  }
  
  // Compute <H_n, \partial_p H_m>
  double HermiteBasis::xLaplacianY(int iOfEltLeft, int iOfEltRight) const
  {
    return (iOfEltLeft == iOfEltRight) * beta_*iOfEltRight;
  }  
  
 

  TensorBasis::TensorBasis(int nbOfVariables0):
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

  int const& TensorBasis::nbOfElts(int iOfVariable) const
  {
    return bases_[iOfVariable]->nbOfElts();
  }

  int& TensorBasis::nbOfElts(int iOfVariable)
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
  
  DVec TensorBasis::basisMeans() const
  {
    DVec vect = basisMeans2(0);
    for (int iOfVariable = 1; iOfVariable < nbOfVariables(); iOfVariable++)
      vect = kron(basisMeans2(iOfVariable), vect);
    return vect;
    //throw runtime_error("TensorBasis::allBasisMeans not implemented");
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
  
  ///
  ///Computes the Gram matrix of the global basis
  ///Recall that the variable order is reversed for the operators
  SMat TensorBasis::gramMatrix() const
  {
    SMat mat = bases_[0]->gramMatrix();
    for (int iOfVariable = 1; iOfVariable < nbOfVariables(); iOfVariable++)
      mat = kron(bases_[iOfVariable]->gramMatrix(), mat);
    return mat;
    //throw runtime_error("TensorBasis::allBasisMeans not implemented");
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
  {}

  int& QPBasis::nbOfQModes()
  {
    return nbOfElts(0);
  }

  const int& QPBasis::nbOfQModes() const
  {
    return nbOfElts(0);
  }

  int& QPBasis::nbOfPModes()
  {
    return nbOfElts(1);
  }

  const int& QPBasis::nbOfPModes() const
  {
    return nbOfElts(1);
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  int QPBasis::iTens(int iOfFourier2, int iOfHermite) const
  {
    assert(iOfFourier2 < nbOfQModes()  && iOfHermite < nbOfPModes());
    return nbOfQModes() * iOfHermite + iOfFourier2;
  }

  vector<int> QPBasis::vecTens(int iTens0) const
  {
    assert(iTens0 < nbOfPModes() * nbOfQModes());
    vector<int> vecIndex(2);
    vecIndex[0] = iTens0 % nbOfQModes();
    vecIndex[1] = iTens0 / nbOfQModes();
    return vecIndex;
  }

  double QPBasis::value(System const& syst, int iOfElt) const
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
    DVec vect = DVec::Zero(totalNbOfElts());
    vect(iTens(iOfElt1, iOfElt2)) = 1;
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
    bases_[0] = new ExpFourierBasis(input, potential);
    bases_[1] = new HermiteBasis(input);
  }
  
  DMat ExpFourierHermiteBasis::convertToTrigBasis(const DMat& X)
  {
    return bases_[0]->expToTrigMat() * X;
  }
  
  // #### HermiteHermiteBasis ####
    
  HermiteHermiteBasis::HermiteHermiteBasis(Input const& input, Potential& potential):
    QPBasis()
  {
    bases_[0] = new HermiteQBasis(input, potential);
    bases_[1] = new HermiteBasis(input);
  }
}












