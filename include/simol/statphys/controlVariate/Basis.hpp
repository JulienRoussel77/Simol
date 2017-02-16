#ifndef SIMOL_BASIS_HPP
#define SIMOL_BASIS_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/system/System.hpp"

namespace simol
{
  class TVec
  {
      vector<int> nbOfElts_;
    public:
      //virtual double operator()(vector<int>& vecOfElt) const = 0;
      //virtual double& operator()(vector<int>& vecOfElt) = 0;
      TVec(vector<int>& nbOfElts0);
      int nbOfVariables();
      int const& nbOfElts(int iOfVariable) const;
      int& nbOfElts(int iOfVariable);
      vector<int> const& nbOfElts() const;
      vector<int>& nbOfElts();
      virtual int size() const = 0;
      int iTens(vector<int>& vecOfElt) const;
  };

  
  class DTVec : public TVec
  {
      DVec data_;
    public:
      DTVec(vector<int>& nbOfElts0);
      virtual int size() const;
      double operator()(vector<int>& vecIndex) const;
      double& operator()(vector<int>& vecIndex);
      double operator()(int iTensOfElt) const;
      double& operator()(int iTensOfElt);
  };

  int product(vector<int>& nbOfElts);
  //DTVec product(DMat& A, DTVec& X, int iOfVariable);

  //class Basis;
  //Basis* createBasis(Input const& input, Potential& potential);
  
  class Basis
  {
  //friend   Basis* createBasis(Input const& input, Potential& potential);
  protected:
    int nbOfElts_;
    DVec basisMeans2_;
    double norm2meansVec_;
  public:
    Basis(const int nbOfElts);
    virtual ~Basis() {};
    virtual int const& nbOfElts() const;
    virtual int& nbOfElts();
    virtual double basisMean2(int iOfElt) const;
    virtual DVec const& basisMeans2() const;
    virtual double norm2meansVec() const;
    
    virtual DMat const& expToTrigMat() const {throw runtime_error("expToTrigMat not defined !");}
    virtual DMat const& trigToExpMat() const {throw runtime_error("trigToExpMat not defined !");}
    
    virtual double value(double variable, const int iOfElt) const = 0;
    virtual DVec gradient(double variable, const int iOfElt) const = 0;
    virtual double laplacian(double variable, const int iOfElt) const = 0;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const = 0;
    virtual void gradMatrix(SMat& A) const = 0;
    virtual double xGradStarY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradStarMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const = 0;
    virtual void laplacianMatrix(SMat& A) const = 0;
  };

  class QBasis : public Basis
  {
  protected:
    double beta_;
    Potential* potential_;
    int nbOfIntegrationNodes_;
    double qRepartitionFct_, basisCoefficient_;
  public:
    QBasis(const int nbOfElts, double beta0, Potential& potential);
    virtual double potential(double variable) const;
    virtual double potDeriv(double variable) const;
    virtual double potLapla(double variable) const;
    
    const double& amplitude() const;
  };

  class ExpFourierBasis : public QBasis
  {
  public:
    DVec expFourierCoeffs_;
    DMat trigToExpMat_, expToTrigMat_;
  public:
    ExpFourierBasis(const int nbOfElts, double beta0, Potential& potential);
    
    virtual DVec const& expFourierCoeffs() const;
    virtual double expFourierCoeff(int iOfElt) const;
        
    void computeBasisMeans();
    void computeExpToTrigMat();
    
    DMat const& expToTrigMat() const {return expToTrigMat_;}
    DMat const& trigToExpMat() const {return trigToExpMat_;}
      
    virtual double value(double variable, const int iOfElt) const;
    virtual DVec gradient(double variable, const int iOfElt) const;
    virtual double laplacian(double variable, const int iOfElt) const;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void laplacianMatrix(SMat& A) const;
  };

  class HermiteBasis : public Basis
  {
    double beta_;
    DMat polyCoeffs_;
  public:
    HermiteBasis(const int nbOfElts, double beta0);
    virtual double value(double variable, const int iOfElt) const;
    virtual DVec gradient(double variable, const int iOfElt) const;
    virtual double laplacian(double variable, const int iOfElt) const;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void laplacianMatrix(SMat& A) const;
  };

  class TensorBasis
  {
  public:
  //protected:
    vector<Basis*> bases_;
    //vector<int> nbOfElts_;
  public:
    TensorBasis(const int nbOfVariables);
    virtual ~TensorBasis();
    virtual int nbOfVariables() const;
    virtual int const& nbOfElts(const int iOfVariable) const;
    virtual int& nbOfElts(const int iOfVariable);
    virtual vector<int> nbOfElts() const;
    virtual int totalNbOfElts() const;
    virtual double basisMean2(int iOfVariable, int iOfElt) const;
    virtual DVec const& basisMeans2(int iOfVariable) const;
    virtual double norm2meansVec(int iOfVariable) const;
    
    const Basis* operator()(int iOfBasis) const;
    Basis* operator()(int iOfBasis);
    virtual double value(System const& syst, const int iOfElt) const = 0;
    virtual double value(System const& syst, vector<int>& vecIndex) const = 0;
    virtual DVec gradientQ(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual DVec gradientQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual double laplacianQ(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual double laplacianQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual DVec gradientP(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual DVec gradientP(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual double laplacianP(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual double laplacianP(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    
    ///
    /// For iOfVariable = 1 and iOfElt = l returns (q,p) -> H_l(p)
    virtual DVec getPartialElement(int iOfVariable, int iOfElt) const = 0;
  };

  class QPBasis : public TensorBasis
  {
  public:
    QPBasis();
    int& nbOfFourier();
    const int& nbOfFourier() const;
    int& nbOfHermite();
    const int& nbOfHermite() const;
    int iTens(int iOfFourier2, int iOfHermite) const;
    vector<int> vecTens(int iTens0) const;
    virtual double value(System const& syst, const int iOfElt) const;
    virtual double value(System const& syst, vector<int>& vecIndex) const;
    virtual DVec gradientQ(System const& syst, int iOfParticle, int iOfCoeff) const;
    virtual DVec gradientQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const;
    virtual double laplacianQ(System const& syst, int iOfParticle, int iOfCoeff) const;
    virtual double laplacianQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const;
    virtual DVec gradientP(System const& syst, int iOfParticle, int iOfCoeff) const;
    virtual DVec gradientP(System const& syst, int iOfParticle, vector<int>& vecIndex) const;
    virtual double laplacianP(System const& syst, int iOfParticle, int iOfCoeff) const;
    virtual double laplacianP(System const& syst, int iOfParticle, vector<int>& vecIndex) const;
    
    virtual DVec getBasisElement(int iOfElt1, int iOfElt2) const;
    virtual DVec getPartialElement(int iOfVariable, int iOfElt) const;
  };

  class ExpFourierHermiteBasis : public QPBasis
  {
  public:
    ExpFourierHermiteBasis(Input const& input, Potential& potential);

    /*virtual DVec const& gVector() const;
    virtual double gVector(int iOfElt) const;
    virtual double norm2meansVec() const;*/
    virtual DMat convertToTrigBasis(const DMat& X);
  };
  
  class HermiteHermiteBasis : public QPBasis
  {
  public:
    HermiteHermiteBasis(Input const& input, Potential& potential);
  };


}

#endif
