#ifndef SIMOL_BASIS_HPP
#define SIMOL_BASIS_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/TensorTools.hpp"
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
    //double norm2meansVec_;
    double beta_;
    
    SMat gramMatrix_, gradMatrix_, laplacianMatrix_;
  public:
    Basis(int nbOfElts, double beta0);
    virtual ~Basis() {};
    virtual int const& nbOfElts() const;
    virtual int& nbOfElts();
    virtual double basisMean2(int iOfElt) const;
    virtual DVec const& basisMeans2() const;
    virtual double norm2meansVec() const;
    SMat const& gramMatrix() const {return gramMatrix_;}
    SMat const& gradMatrix() const {return gradMatrix_;}
    SMat const& laplacianMatrix() const {return laplacianMatrix_;}

    
    virtual DMat const& expToTrigMat() const {throw runtime_error("expToTrigMat not defined !");}
    virtual DMat const& trigToExpMat() const {throw runtime_error("trigToExpMat not defined !");}
    
    virtual double value(double variable, int iOfElt) const = 0;
    virtual DVec gradient(double variable, int iOfElt) const = 0;
    virtual double laplacian(double variable, int iOfElt) const = 0;
    
    virtual double xY(int iOfEltLeft, int iOfEltRight) const = 0;
    virtual void computeGramMatrix(); 
    virtual double xGradY(int iOfElementLeft, int iOfElementRight) const = 0;
    virtual void computeGradMatrix();
    //virtual double xGradStarY(int iOfElementLeft, int iOfElementRight) const;
    virtual double xLaplacianY(int iOfElementLeft, int iOfElementRight) const = 0;
    virtual void computeLaplacianMatrix();
    virtual double const& omega() const {throw runtime_error("gamma not defined !");};
  };
  
  class QBasis : public Basis
  {
  protected:
    Potential* potential_;
    double integrationStep_;
    double qRepartitionFct_, basisCoefficient_;
    DVec measureMomenta_;
  public:
    QBasis(Input const& input, Potential& potential0);
    virtual double length() const = 0;
    
    virtual double potential(double variable) const;
    virtual double potDeriv(double variable) const;
    virtual double potLapla(double variable) const;
    
    const double& amplitude() const;
    int nbOfIntegrationSteps() const;
   
    virtual void computeBasisMeans() = 0;
  };

  class ExpFourierBasis : public QBasis
  {
  public:
    DVec expFourierCoeffs_;
    DMat trigToExpMat_, expToTrigMat_;
  public:
    ExpFourierBasis(Input const& input, Potential& potential);
    
    virtual double length() const;
    virtual DVec const& expFourierCoeffs() const;
    virtual double expFourierCoeff(int iOfElt) const;
        
    virtual void computeBasisMeans();
    void computeExpToTrigMat();
    
    DMat const& expToTrigMat() const {return expToTrigMat_;}
    DMat const& trigToExpMat() const {return trigToExpMat_;}
      
    virtual double value(double variable, int iOfElt) const;
    virtual DVec gradient(double variable, int iOfElt) const;
    virtual double laplacian(double variable, int iOfElt) const;
    
    virtual double xY(int iOfEltLeft, int iOfEltRight) const;
    virtual double xGradY(int iOfElementLeft, int iOfElementRight) const;
    virtual double xLaplacianY(int iOfElementLeft, int iOfElementRight) const;
  };
  
  class HermiteQBasis : public QBasis
  {
    double length_;
    double qMin_;
    double omega_;
    DMat polyCoeffs_;
  public:
    HermiteQBasis(Input const& input, Potential& potential0);
    virtual double length() const;
    virtual void computeMeasureMomenta();
    virtual void computeBasisMeans();
    virtual double const& omega() const {return omega_;}
    
    virtual double value(double variable, int iOfElt) const;
    virtual DVec gradient(double variable, int iOfElt) const;
    virtual double laplacian(double variable, int iOfElt) const;
    
    virtual double xY(int iOfEltLeft, int iOfEltRight) const;
    virtual double xGradY(int iOfElementLeft, int iOfElementRight) const;
    virtual double xLaplacianY(int iOfElementLeft, int iOfElementRight) const;
    
    virtual void computeGradMatrix();
    virtual void computeLaplacianMatrix();
  };

  class HermiteBasis : public Basis
  {
    DMat polyCoeffs_;
  public:
    HermiteBasis(Input const& input);
    virtual double value(double variable, int iOfElt) const;
    virtual DVec gradient(double variable, int iOfElt) const;
    virtual double laplacian(double variable, int iOfElt) const;
    
    virtual double xY(int iOfEltLeft, int iOfEltRight) const;
    virtual double xGradY(int iOfElementLeft, int iOfElementRight) const;
    virtual double xLaplacianY(int iOfElementLeft, int iOfElementRight) const;
  };

  class TensorBasis
  {
  public:
  //protected:
    vector<Basis*> bases_;
    //vector<int> nbOfElts_;
    //SMat gramMatrix_;
  public:
    TensorBasis(int nbOfVariables);
    virtual ~TensorBasis();
    virtual int nbOfVariables() const;
    virtual int const& nbOfElts(int iOfVariable) const;
    virtual int& nbOfElts(int iOfVariable);
    virtual vector<int> nbOfElts() const;
    virtual int totalNbOfElts() const;
    virtual DVec basisMeans() const;
    virtual double basisMean2(int iOfVariable, int iOfElt) const;
    virtual DVec const& basisMeans2(int iOfVariable) const;
    virtual double norm2meansVec(int iOfVariable) const;
    SMat gramMatrix() const;
    
    const Basis* operator()(int iOfBasis) const;
    Basis* operator()(int iOfBasis);
    virtual double value(System const& syst, int iOfElt) const = 0;
    virtual double value(System const& syst, vector<int>& vecIndex) const = 0;
    virtual DVec gradientQ(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual DVec gradientQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual double laplacianQ(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual double laplacianQ(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual DVec gradientP(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual DVec gradientP(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    virtual double laplacianP(System const& syst, int iOfParticle, int iOfCoeff) const = 0;
    virtual double laplacianP(System const& syst, int iOfParticle, vector<int>& vecIndex) const = 0;
    
    virtual DVec getBasisElement(int iOfElt1, int iOfElt2) const = 0;
    ///
    /// For iOfVariable = 1 and iOfElt = l returns (q,p) -> H_l(p)
    virtual DVec getPartialElement(int iOfVariable, int iOfElt) const = 0;
  };

  class QPBasis : public TensorBasis
  {
  public:
    QPBasis();
    int& nbOfQModes();
    const int& nbOfQModes() const;
    int& nbOfPModes();
    const int& nbOfPModes() const;
    int iTens(int iOfFourier2, int iOfHermite) const;
    vector<int> vecTens(int iTens0) const;
    virtual double value(System const& syst, int iOfElt) const;
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
