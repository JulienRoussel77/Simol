#ifndef SIMOL_GALERKIN_HPP
#define SIMOL_GALERKIN_HPP

#include "simol/statphys/Tools.hpp"
#include "simol/statphys/input/Input.hpp"
#include "simol/statphys/potential/AllPotentials.hpp"
#include "simol/statphys/controlVariate/Basis.hpp"
#include "simol/statphys/controlVariate/CVBasis.hpp"
//#include "simol/statphys/controlVariate/MatrixFree.hpp"
#include "simol/statphys/TensorTools.hpp"


#include <unsupported/Eigen/KroneckerProduct>

namespace simol
{

  /*SMat kron(const SMat& A, const SMat& B);
  DMat kron(const DMat& A, const DMat& B);
  SMat kron(const DMat& A, const SMat& B);
  double computeSpectralGap(SMat const& A);*/

  //class Galerkin;
  //Galerkin* createGalerkin(Input const& input);

  class Galerkin
  {
    //friend Galerkin* createGalerkin(Input const& input);

    protected:
      int nbOfParticles_;
      int nbOfFourier_, nbOfHermite_, maxOfFourier_;
      int sizeOfBasis_;
      SMat SIdQ_, SIdP_;
      DMat DIdQ_, DIdP_;
      SMat Q_, P_;
      SMat tQ_, tP_;
      SMat Lthm0_;
      SMat Lthm_, Lham_;
      SMat Lrep_;
      SMat Leq_, Leta_;
      double beta_, gamma_;
      double amplitude_;
      double externalForce_;
      bool doNonequilibrium_;
      int nbOfIntegrationNodes_;
      DMat trigToExpMat_, expToTrigMat_;
      SMat trigToExpTens_, expToTrigTens_;
      Potential* potential_;
      ExpFourierHermiteBasis basis_;
      //Basis* basis_;
    public:
      Galerkin(Input const& input);
      virtual ~Galerkin();

      virtual int nbOfVariables() const;
      virtual int nbOfParticles() const;
      
      virtual int nbOfFourier() const;
      virtual int nbOfHermite() const;
      virtual int sizeOfBasis() const;
      
      const double& gamma() const;
      const bool& doNonequilibrium() const {return doNonequilibrium_;}
      ExpFourierHermiteBasis const& basis() const {return basis_;}
      ExpFourierHermiteBasis& basis() {return basis_;}
      
      SMat const& Leq() const {return Leq_;}
      SMat const& Leta() const {return Leta_;}

      const double& expFourierMeans(int iOfElt) const;
      DVec const& gVector() const;
      double gVector(int iOfElt) const;
      const double& norm2gVector() const;
      int iTens(int iOfFourier2, int iOfHermite) const;
      DMat shapeSaddle(const DMat& A) const;
      SMat shapeSaddle(const SMat& A) const;
      //DMat shapePrec(const DMat& A) const;
      DMat unshapeSaddle(const DMat& Asad) const;
      SMat unshapeSaddle(const SMat& Asad) const;
      DVec shapeSaddle(const DVec& X) const;
      DVec unshapeSaddle(const DVec& Xsad) const;
      DVec solve(const DMat& A, const DVec& X) const;
      DVec solve(const SMat& A, const DVec& X) const;
      DVec solveWithGuess(const SMat& A, const DVec& X, DVec const& X0) const;
      DVec solveWithSaddle(const DMat& A, const DVec& X) const;
      DVec solveWithSaddle(const SMat& A, const DVec& X) const;
      DVec solveWithSaddleAndGuess(const SMat& A, const DVec& X, const DVec& X0) const;
      DMat invWithSaddle(const SMat& A) const;
      DMat invWithSaddle(const DMat& A) const;
      
      Eigen::EigenSolver<Eigen::MatrixXd> getEigenSolver() const;

      void computeExpToTrigMat();
      virtual void computeExpToTrigTens() = 0;
      DMat convertToTrigBasis(const DMat& X);
      DVec projectionOrthoG(DVec const& X) const;
      void createQ();
      void createP();
      virtual void compute() = 0;

      DVec gettGiHj(int i, int j) const;
      DVec gettGiHjTrig(int i, int j) const;
      DVec getLtGiHj(int i, int j) const;
      DVec getLtGiHjTrig(int i, int j) const;
      DVec getLinvtGiHj(int i, int j) const;
      DVec getLinvtGiHjTrig(int i, int j) const;
      //SMat CVcoeffs() const;
      virtual DVec CVcoeffsVec() const = 0;
      CVBasis makeCvBasis();
        
      void computeEigen() const;
  };
  
  class OverdampedGalerkin : public Galerkin
  {
    public:
      OverdampedGalerkin(Input const& input);
      void createLeta();
      virtual void computeExpToTrigTens();
      virtual void compute();
      DVec getGradV() const;
      DVec CVcoeffsVec() const;
  };
}

#endif
