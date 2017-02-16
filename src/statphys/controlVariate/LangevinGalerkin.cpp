#include "simol/statphys/controlVariate/LangevinGalerkin.hpp"

//#### LangevinGalerkin ####

namespace simol
{
  LangevinGalerkin::LangevinGalerkin(Input const& input):
    Galerkin(input)
  {
    tensorBasis_ = new ExpFourierHermiteBasis(input, *potential_);
    double totalTime = clock();
    
    tensorBasis(0)->gradMatrix(Q_);
    tQ_ = Q_.adjoint();
    
    for (int iOfHermite = 1; iOfHermite < nbOfHermite_; iOfHermite++)
      P_.insert(iOfHermite - 1, iOfHermite) = sqrt(beta_ * iOfHermite);
    tP_ = P_.adjoint();
    
    
    //computeExpToTrigTens();
    //cout << "- computeExp : " << clock() - totalTime << endl; totalTime = clock();
    SMat PQt = kron(tP_, Q_);
    SMat PtQ = PQt.transpose();
    Lham_ = -PQt + PtQ;   // The matricies and the operators have opposite signs
    //Lham_ = kron(Q_, tP_) - kron(tQ_, P_);
    //display(Lham_, outputFolderName()+"Lham");
    createLthm();
    //cout << "- Lham and Lthm : " << clock() - totalTime << endl; totalTime = clock();
    //display(Lthm_, outputFolderName()+"Lthm");
    Leq_  = Lham_ + gamma_ * Lthm_;
    //cout << "- Leq : " << clock() - totalTime << endl; totalTime = clock();
    //Lrep_ = kron(SMat::Identity(nbOfFourier_) , P_);
    Lrep_ = -kron(P_, SIdQ_);
    Leta_ = Leq_ + externalForce_ * Lrep_;
    cout << "- LangevinGalerkin : " << clock() - totalTime << endl; totalTime = clock();
  }
  
  ///
  /// The matricies Leq, Leta are positive !
  void LangevinGalerkin::createLthm0()
  {
    tensorBasis(1)->laplacianMatrix(Lthm0_);
    Lthm0_ /= beta_;
    /*cout << "createLthm0" << endl;
    for (int iOfHermite = 1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      Lthm0_(iOfHermite, iOfHermite) = - iOfHermite;
    cout << "end createLthm0" << endl;*/
  }

  void LangevinGalerkin::createLthm()
  {
    createLthm0();
    Lthm_ = kron(Lthm0_, SIdQ_);
  }

  /*void LangevinGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = kron(SIdP_, trigToExpMat_);
    expToTrigTens_ = kron(SIdP_, expToTrigMat_);
  }*/
  
  //### PeriodicLangevinGalerkin ###
  
  PeriodicLangevinGalerkin::PeriodicLangevinGalerkin(Input const& input):
    LangevinGalerkin(input)
  {}
  
  void PeriodicLangevinGalerkin::compute()
  {

    //DMat DLeq(Leq_);
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    //DMat LeqInv = invWithSaddle(DLeq);
    DMat DLeqInv = invWithSaddle(DMat(Leq_));
    

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");


    DVec H1Trig = DVec::Zero(sizeOfBasis_);
    H1Trig(iTens(0, 1)) = 1;
    
    

    //Map<DMat, Eigen::Aligned> H1TrigMat(H1Trig, nbOfFourier_, nbOfHermite_);
    DMat H1TrigMat = reshape(H1Trig, nbOfFourier_, nbOfHermite_);

    //DMat H1Mat = trigToExpMat_ * H1TrigMat;
    //DVec H1 = reshape(H1Mat, nbOfFourier_ * nbOfHermite_, 1);

    //DVec H1 = trigToExpTens_ * H1Trig;
    DVec H1 = tensorBasis()->getPartialElement(1,1);

    /*DVec H1(sizeOfBasis_);
    for (int iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
      H1(iTens(iOfFourier2, 1)) = expFourierMeans_[iOfFourier2];*/

    cout << "############ H1Mat ############" << endl;
    display(H1, outputFolderName()+"H1");
    //DMat H1Mat = reshape(gettGiHj(0,1), nbOfFourier_, nbOfHermite_);
    DMat H1Mat = reshape(H1, nbOfFourier_, nbOfHermite_);
    display(H1Mat, outputFolderName()+"H1Mat");

    display(H1Trig, outputFolderName()+"H1Trig");
    display(H1TrigMat, outputFolderName()+"H1TrigMat");
    
    DVec LinvH1 = solveWithSaddle(Leq(), H1);
    //display(LinvH1, outputFolderName()+"LinvH1");
    DMat LinvH1Mat = reshape(LinvH1, nbOfFourier_, nbOfHermite_);
    display(LinvH1Mat, outputFolderName()+"LinvH1Mat");
    
    DVec LLinvH1 = Leq() * H1;
    //display(LLinvH1, outputFolderName()+"LLinvH1");
    DMat LLinvH1Mat = reshape(LLinvH1, nbOfFourier_, nbOfHermite_);
    display(LLinvH1Mat, outputFolderName()+"LLinvH1Mat");

    double varOfH1 = 2 * dot(H1, LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, Lrep_ * LinvH1);
    //double varCoeff = -.5 * dot( Lrep_ * LinvH1, LeqInv * Lrep_ * LinvH1);
    double varCoeff = 2 * dot( Lrep_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }
  
  ///
  ///Returns the control variate in the Galerkin basis
  DVec PeriodicLangevinGalerkin::CVcoeffsVec() const
  {
    if (doNonequilibrium())
    {
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      SMat Keta = Leta_.transpose() * Leta_;
      //DVec CV = -solve(Keta, Leta_.transpose() * gettGiHj(0,1));
      DVec CV = -solve(Keta, Leta_.transpose() * tensorBasis()->getPartialElement(1,1));
      return CV;
    }
    else
    {
      cout << "Computing the control variate by an equilibrium LLT" << endl;
      //SMat Keq = Leq_.transpose() * Leq_;
      //return -solve(Keq, Leq_.transpose() * gettGiHj(0,1));
      
      Eigen::BiCGSTAB<SMat> solver(Leq_);
      //return -solver.solve(-projectionOrthoG(gettGiHj(0,1)));
      return -solver.solve(-projectionOrthoG(tensorBasis()->getPartialElement(1,1)));
    }
  }
  
  //### ColloidLangevinGalerkin
  
  ColloidLangevinGalerkin::ColloidLangevinGalerkin(Input const& input):
    LangevinGalerkin(input)
  {}
  
  ///
  ///
  void ColloidLangevinGalerkin::compute()
  {
    throw runtime_error("ColloidLangevinGalerkin::compute() not defined !");
  }

  ///
  ///Returns the control variate in the Galerkin basis
  DVec ColloidLangevinGalerkin::CVcoeffsVec() const
  {
    throw runtime_error("ColloidLangevinGalerkin::CVcoeffsVec() not defined !");
  }
  
  
  //#### BoundaryLangevinGalerkin ####
  
  
  

  SMat tensorPower(SMat const& A, int power)
  {
    SMat temp = A;
    for (int i = 1; i < power; i++)
      temp = kron(temp, A);
    return temp;
  }

  DMat tensorPower(DMat const& A, int power)
  {
    DMat temp = A;
    for (int i = 1; i < power; i++)
      temp = kron(temp, A);
    return temp;
  }

  SMat BoundaryLangevinGalerkin::PMatToTens(SMat const& PMat, int iOfParticleP)
  {
    assert(iOfParticleP < nbOfParticles_);
    //SMat tempId = speye(PMat.rows(), PMat.cols());
    SMat res = DMat::Identity(1, 1).sparseView();
    for (int i = 0; i < nbOfParticles_; i++)
    {
      res = kron(res, SIdQ_);
      if (i == iOfParticleP) res = kron(res, PMat);
      else res = kron(res, SIdP_);
    }
    return res;
  }

  ///
  ///Creates the tensor of size sizeOfBasis X sizeOfBasis that is identity
  ///except for indices iOfParticleQ and iOfParticleP
  SMat BoundaryLangevinGalerkin::doubleMatToTens(SMat const& QMat, SMat const& PMat, int iOfParticleQ, int iOfParticleP)
  {
    cout << "doubleMatToTens(SMat const& QMat, SMat const& PMat, int iOfParticleQ, int iOfParticleP)" << endl;
    /*assert(iOfVariableA < iOfVariableB);
    assert(iOfVariableB < nbOfVariables);
    assert(A.size() == B.size());*/
    //SMat tempId = speye(A.rows(), A.cols());
    SMat res = DMat::Identity(1, 1).sparseView();
    for (int i = 0; i < nbOfParticles_; i++)
    {
      if (i == iOfParticleQ) res = kron(res, QMat);
      else res = kron(res, SIdQ_);
      if (i == iOfParticleP) res = kron(res, PMat);
      else res = kron(res, SIdP_);
    }
    cout << "res : " << res.rows() << " x " << res.cols() << endl;
    return res;
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  //Allows to access elements involving a single particle !
  int BoundaryLangevinGalerkin::iTens(int iOfFourier2, int iOfHermite, int iOfParticle) const
  {
    assert(iOfFourier2 < nbOfFourier_ && iOfHermite < nbOfHermite_ && iOfParticle < nbOfParticles_);
    return pow(nbOfFourier_ * nbOfHermite_ , iOfParticle) * ( nbOfFourier_ * iOfHermite + iOfFourier2);
  }


  BoundaryLangevinGalerkin::BoundaryLangevinGalerkin(Input const& input):
    Galerkin(input)
  {
    //computeExpToTrigTens();
    
    createLham();
    createLthm();

    Leq_ = Lham_ + gamma_ * Lthm_;
  }
  
  void BoundaryLangevinGalerkin::createLthm0()
  {
    cout << "createLthm0" << endl;
    for (int iOfHermite = 1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      Lthm0_.insert(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
    cout << "end createLthm0" << endl;
  }

  void BoundaryLangevinGalerkin::createLham()
  {
    cout << "BoundaryLangevinGalerkin::createLham()" << endl;
    //Lham_ = arma::zeros(sizeOfBasis_, sizeOfBasis_);
    for (int i = 0; i < nbOfParticles_; i++)
    {
      Lham_ += doubleMatToTens(Q_, tP_, i, i);
      if (i > 0) Lham_ -= doubleMatToTens(Q_, tP_, i, i - 1);
      Lham_ -= doubleMatToTens(tQ_, P_, i, i);
      if (i < nbOfParticles_ - 1) Lham_ += doubleMatToTens(Q_, tP_, i + 1, i);
    }
    Lham_ /= beta_;
    display(Lham_, "output/Galerkin/Lham");
    cout << "end BoundaryLangevinGalerkin::createLham()" << endl;
  }

  void BoundaryLangevinGalerkin::createLthm()
  {
    createLthm0();
    //SMat Lthm1_ = kron(SIdQ_, Lthm0_);
    //Lthm_ = arma::zeros(sizeOfBasis_, sizeOfBasis_);
    for (int i = 0; i < nbOfParticles_; i++)
    {
      //cout << i << " < " << nbOfParticles_ << endl;
      Lthm_ += PMatToTens(Lthm0_, i);
    }
    Lthm_ /= beta_;
    display(Lthm_, "output/Galerkin/Lthm");
  }


  void BoundaryLangevinGalerkin::compute()
  {
  }
}