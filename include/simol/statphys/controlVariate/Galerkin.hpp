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

  class Galerkin
  {

    protected:
      int nbOfParticles_;
      int nbOfQModes_, nbOfPModes_;
      int sizeOfBasis_;
      SMat SIdQ_, SIdP_;
      DMat DIdQ_, DIdP_;
      SMat Q_, P_;
      SMat Qt_, Pt_;
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
      TensorBasis* tensorBasis_;
      DVec SU_;
      
      string outputFolderName_;
      //Basis* basis_;
    public:
      Galerkin(Input const& input);
      virtual ~Galerkin();

      virtual int nbOfVariables() const;
      virtual int nbOfParticles() const;
      
      virtual int nbOfQModes() const;
      virtual int nbOfPModes() const;
      virtual int sizeOfBasis() const;
      
      double gamma() const;
      bool doNonequilibrium() const {return doNonequilibrium_;}
      TensorBasis* const& tensorBasis() const {return tensorBasis_;}
      TensorBasis* tensorBasis() {return tensorBasis_;}
      Basis const* tensorBasis(int iOfVariable) const {return tensorBasis_->bases_[iOfVariable];}
      Basis* tensorBasis(int iOfVariable) {return tensorBasis_->bases_[iOfVariable];}
      
      SMat const& Leq() const {return Leq_;}
      SMat const& Leta() const {return Leta_;}

      double basisMean(int iOfVariable, int iOfElt) const;
      const DVec& basisMeans(int iOfVariable) const;
      /*DVec const& gVector() const;
      double gVector(int iOfElt) const;*/
      double norm2meansVec(int iOfVariable) const;
      int iTens(int iOfFourier2, int iOfHermite) const;
      string outputFolderName() const {return outputFolderName_;}
      
      DMat shapeSaddle(const DMat& A, const DMat& C) const;
      SMat shapeSaddle(const SMat& A, const DMat& C) const;
      DMat unshapeSaddle(const DMat& Asad, const DMat& C) const;
      SMat unshapeSaddle(const SMat& Asad, const DMat& C) const;
      DVec shapeSaddle(const DVec& X, const DMat& C) const;
      DVec unshapeSaddle(const DVec& Xsad, const DMat& C) const;
      DVec solve(const DMat& A, const DVec& X) const;
      DVec solve(const SMat& A, const DVec& X) const;
      DVec solveWithGuess(const SMat& A, const DVec& X, DVec const& X0) const;
      DVec solveWithSaddle(const DMat& A, const DVec& X, const DMat& C) const;
      DVec solveWithSaddle(const SMat& A, const DVec& X, const DMat& C) const;
      DVec solveWithSaddleAndGuess(const SMat& A, const DVec& X, const DVec& X0, const DMat& C) const;
      DMat invWithSaddle(const SMat& A, const DMat& C) const;
      DMat invWithSaddle(const DMat& A, const DMat& C) const;
      
      Eigen::EigenSolver<Eigen::MatrixXd> getEigenSolver() const;

      DVec projectionOrthoG(DVec const& X) const;
      virtual void compute() = 0;
      virtual DMat computeConstraints(double tol) const;
      void studyLangevinErrors(bool doComputeRef);
      void studyColloidLangevinErrors(bool doComputeRef);
      
      virtual DVec CVObservable() const = 0;
      virtual DVec CVcoeffsVec() const = 0;
      CVBasis makeCvBasis();
        
      void computeEigen() const;
  };
  
}

#endif
