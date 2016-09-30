#ifndef SIMOL_BASIS_HPP
#define SIMOL_BASIS_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/potential/Potential.hpp"
#include "simol/statphys/system/Particle.hpp"

namespace simol
{
  class TVec
  {
      vector<int> nbOfElts_;
    public:
      //virtual double const& operator()(vector<int>& vecOfElt) const = 0;
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
      Vector<double> data_;
    public:
      DTVec(vector<int>& nbOfElts0);
      virtual int size() const;
      double const& operator()(vector<int>& vecIndex) const;
      double& operator()(vector<int>& vecIndex);
      double const& operator()(int iTensOfElt) const;
      double& operator()(int iTensOfElt);
  };

  int product(vector<int>& nbOfElts);
  //DTVec product(DMat& A, DTVec& X, int iOfVariable);

  class Basis
  {
  protected:
    int nbOfElts_;
  public:
    Basis(const int nbOfElts);
    virtual ~Basis() {};
    virtual int const& nbOfElts() const;
    virtual int& nbOfElts();
    virtual double const& expFourierMeans(int iOfElt) const;
    virtual Vector<double> const& gVector() const;
    virtual double const& norm2gVector() const;
    virtual double value(double variable, const int iOfElt) const = 0;
    virtual Vector<double> gradient(double variable, const int iOfElt) const = 0;
    virtual double laplacian(double variable, const int iOfElt) const = 0;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const = 0;
    virtual void gradMatrix(SMat& A) const = 0;
    virtual double xGradStarY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradStarMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const = 0;
    virtual void laplacianMatrix(SMat& A) const = 0;
  };

  class FourierBasis : public Basis
  {
  public:
    FourierBasis(const int nbOfElts);
    //virtual int nbOfFreq() const;
    virtual double value(double variable, const int iOfElt) const;
    virtual Vector<double> gradient(double variable, const int iOfElt) const;
    virtual double laplacian(double variable, const int iOfElt) const;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void laplacianMatrix(SMat& A) const;
  };

  class ExpFourierBasis : public Basis
  {
  protected:
    double beta_;
    Potential* potential_;
    int nbOfIntegrationNodes_;
    double qRepartitionFct_, basisCoefficient_;
    Vector<double> expFourierMeans_;
    Vector<double> gVector_;
    double norm2gVector_;
  public:
    ExpFourierBasis(const int nbOfElts, double beta0, Potential& potential);
    //virtual int nbOfFreq() const;
    const double& amplitude() const;

    double const& expFourierMeans(int iOfElt) const;
    Vector<double> const& gVector() const;
    double const& norm2gVector() const;
    
    virtual double potential(double variable) const;
    virtual double potDeriv(double variable) const;
    virtual double potLapla(double variable) const;
    
    double computeExpFourierMeans();
      
    virtual double value(double variable, const int iOfElt) const;
    virtual Vector<double> gradient(double variable, const int iOfElt) const;
    virtual double laplacian(double variable, const int iOfElt) const;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void laplacianMatrix(SMat& A) const;
  };
  
  /*class SinExpFourierBasis : public ExpFourierBasis
  {
    public:
      SinExpFourierBasis(const int nbOfElts, double beta0, Potential& potential);
      const double& amplitude() const;
      
      virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
      virtual double xGradStarY(const int iOfElementLeft, const int iOfElementRight) const;
      virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
  };*/

  class HermiteBasis : public Basis
  {
    double beta_;
    DMat polyCoeffs_;
  public:
    HermiteBasis(const int nbOfElts, double beta0);
    virtual double value(double variable, const int iOfElt) const;
    virtual Vector<double> gradient(double variable, const int iOfElt) const;
    virtual double laplacian(double variable, const int iOfElt) const;
    
    virtual double xGradY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void gradMatrix(SMat& A) const;
    virtual double xLaplacianY(const int iOfElementLeft, const int iOfElementRight) const;
    virtual void laplacianMatrix(SMat& A) const;
  };

  class TensorBasis
  {
    protected:
      vector<Basis*> bases_;
      //vector<int> nbOfElts_;
    public:
      TensorBasis(const int nbOfVariables);
      virtual ~TensorBasis();
      virtual int nbOfVariables() const;
      virtual int const& nbOfElts(const int iOfVariable) const;
      virtual int& nbOfElts(const int iOfVariable);
      virtual vector<int> nbOfElts() const;
      const Basis* operator()(const int iOfBasis) const;
      virtual double value(vector<Particle*> const& configuration, const int iOfElt) const = 0;
      virtual double value(vector<Particle*> const& configuration, vector<int>& vecIndex) const = 0;
      virtual Vector<double> gradientQ(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const = 0;
      virtual Vector<double> gradientQ(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const = 0;
      virtual double laplacianQ(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const = 0;
      virtual double laplacianQ(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const = 0;
      virtual Vector<double> gradientP(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const = 0;
      virtual Vector<double> gradientP(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const = 0;
      virtual double laplacianP(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const = 0;
      virtual double laplacianP(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const = 0;
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
      virtual double value(vector<Particle*> const& configuration, const int iOfElt) const;
      virtual double value(vector<Particle*> const& configuration, vector<int>& vecIndex) const;
      virtual Vector<double> gradientQ(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const;
      virtual Vector<double> gradientQ(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const;
      virtual double laplacianQ(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const;
      virtual double laplacianQ(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const;
      virtual Vector<double> gradientP(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const;
      virtual Vector<double> gradientP(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const;
      virtual double laplacianP(vector<Particle*> const& configuration, int iOfParticle, int iOfCoeff) const;
      virtual double laplacianP(vector<Particle*> const& configuration, int iOfParticle, vector<int>& vecIndex) const;
  };

  class FourierHermiteBasis : public QPBasis
  {
    public:
      FourierHermiteBasis(Input const& input);
  };

  class ExpFourierHermiteBasis : public QPBasis
  {
  public:
    ExpFourierHermiteBasis(Input const& input, Potential& potential);
    virtual double const& expFourierMeans(int iOfElt) const;
    virtual Vector<double> const& gVector() const;
    virtual double const& norm2gVector() const;
  };

  /*class TrigBasis : public TensorBasis
  {
  public:
    TensorBasis(const int nbOfElts1, const int nbOfElts2);
    virtual double valueOfElt(double variable, const int iOfElt, const int iOfVariable);
  };*/
}

#endif
