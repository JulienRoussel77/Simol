#include "simol/statphys/controlVariate/LangevinGalerkin.hpp"

//#### LangevinGalerkin ####

namespace simol
{
  LangevinGalerkin::LangevinGalerkin(Input const& input):
    Galerkin(input)
  {}
  
  ///
  ///Should be called in all the constructors of the daughters of LangevinGalerkin
  /// /!\ By convention for the operators the variables are in the order (p,q) and for the states the variables are (q,p)
  void LangevinGalerkin::createOperators()
  {
    Q_ = tensorBasis(0)->gradMatrix();
    Qt_ = Q_.adjoint();
    display(Q_, outputFolderName()+"Q");
    
    P_ = tensorBasis(1)->gradMatrix();
    /*for (int iOfHermite = 1; iOfHermite < nbOfPModes_; iOfHermite++)
      P_.insert(iOfHermite - 1, iOfHermite) = sqrt(beta_ * iOfHermite);*/
    Pt_ = P_.adjoint();
    display(P_, outputFolderName()+"P");
    
    SMat PQt = kron(Pt_, Q_);
    SMat PtQ = PQt.transpose();
    Lham_ = -PQt + PtQ;   // The matrices and the operators have opposite signs
    createLthm();
    display(Lthm0_, outputFolderName()+"Lthm0");
    Leq_  = Lham_ + gamma_ * Lthm_;
    Lrep_ = -kron(P_, SIdQ_);
    Leta_ = Leq_ + externalForce_ * Lrep_;
    display(Leq_, outputFolderName()+"Leq");
    cout << "Leq : " << Leq_.rows() << " " << Leq_.cols() << endl;
    SMat Lsad = shapeSaddle(Leq_);
    display(Lsad, outputFolderName()+"Lsad");
    cout << "Lsad : " << Lsad.rows() << " " << Lsad.cols() << endl;
    display(tensorBasis_->gramMatrix(), outputFolderName()+"gramMat");
    display(tensorBasis(0)->gramMatrix(), outputFolderName()+"qGramMat");
    display(tensorBasis(1)->gramMatrix(), outputFolderName()+"pGramMat");
    
    //cout << "Leq spectral gap : " << computeSpectralGap(Leq_) << endl;
    cout << "Lsad spectral gap : " << computeSpectralGap(Lsad) << endl;
  }
  
  ///
  /// The matricies Leq, Leta are positive !
  void LangevinGalerkin::createLthm0()
  {
    Lthm0_ = tensorBasis(1)->laplacianMatrix();
    Lthm0_ /= beta_;
  }

  void LangevinGalerkin::createLthm()
  {
    createLthm0();
    Lthm_ = kron(Lthm0_, SIdQ_);
  }
  
  ///
  ///Returns the control variate in the Galerkin basis
  DVec LangevinGalerkin::CVcoeffsVec() const
  {
    if (doNonequilibrium())
    {
      throw runtime_error("Nonequilibrium Galerkin not implemented !");
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      SMat Keta = Leta_.transpose() * Leta_;
      //DVec CV = -solve(Keta, Leta_.transpose() * gettGiHj(0,1));
      DVec CV = -solve(Keta, Leta_.transpose() * CVObservable());
      return CV;
    }
    else
    {
      cout << "Computing the control variate by an equilibrium LLT" << endl;
      
      /*SMat qGramMat(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
      tensorBasis(0)->gramMatrix(qGramMat);
      SMat gramMat = kron(qGramMat, SIdP_);*/
      SMat gramMat = tensorBasis()->gramMatrix();
      //cout << "CVObs :" << endl << CVObservable() << endl;
      cout << "basisMeansQ :" << endl << basisMeans2(0) << endl;
      //cout << "basisMeansP :" << endl << basisMeans2(1) << endl;
      DVec LinvRObs = solveWithSaddle(Leq(), gramMat * CVObservable());
      return LinvRObs;
    }
  }
  
  //### PeriodicLangevinGalerkin ###
  
  PeriodicLangevinGalerkin::PeriodicLangevinGalerkin(Input const& input):
    LangevinGalerkin(input)
  {
    tensorBasis_ = new ExpFourierHermiteBasis(input, *potential_);
    createOperators();
  }
  
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
    
    //Map<DMat, Eigen::Aligned> H1TrigMat(H1Trig, nbOfQModes_, nbOfPModes_);
    DMat H1TrigMat = reshape(H1Trig, nbOfQModes_, nbOfPModes_);

    DVec H1 = tensorBasis()->getPartialElement(1,1);

    /*DVec H1(sizeOfBasis_);
    for (int iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
      H1(iTens(iOfFourier2, 1)) = expFourierMeans_[iOfFourier2];*/

    cout << "############ H1Mat ############" << endl;
    display(H1, outputFolderName()+"H1");
    //DMat H1Mat = reshape(gettGiHj(0,1), nbOfQModes_, nbOfPModes_);
    DMat H1Mat = reshape(H1, nbOfQModes_, nbOfPModes_);
    display(H1Mat, outputFolderName()+"H1Mat");

    display(H1Trig, outputFolderName()+"H1Trig");
    display(H1TrigMat, outputFolderName()+"H1TrigMat");
    
    DVec LinvH1 = solveWithSaddle(Leq(), H1);
    //display(LinvH1, outputFolderName()+"LinvH1");
    DMat LinvH1Mat = reshape(LinvH1, nbOfQModes_, nbOfPModes_);
    display(LinvH1Mat, outputFolderName()+"LinvH1Mat");
    
    DVec LLinvH1 = Leq() * H1;
    //display(LLinvH1, outputFolderName()+"LLinvH1");
    DMat LLinvH1Mat = reshape(LLinvH1, nbOfQModes_, nbOfPModes_);
    display(LLinvH1Mat, outputFolderName()+"LLinvH1Mat");

    double varOfH1 = 2 * dot(H1, LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, Lrep_ * LinvH1);
    //double varCoeff = -.5 * dot( Lrep_ * LinvH1, LeqInv * Lrep_ * LinvH1);
    double varCoeff = 2 * dot( Lrep_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }
  
  DVec PeriodicLangevinGalerkin::CVObservable() const
  {
    return tensorBasis()->getPartialElement(1,1);
  }
  

  
  //### ColloidLangevinGalerkin
  
  ColloidLangevinGalerkin::ColloidLangevinGalerkin(Input const& input):
    LangevinGalerkin(input)
  {
    tensorBasis_ = new HermiteHermiteBasis(input, *potential_);
    createOperators();
  }
  
  ///
  ///
  void ColloidLangevinGalerkin::compute()
  {
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    DMat DLeqInv = invWithSaddle(DMat(Leq_));
    
    DMat Lsad = shapeSaddle(Leq_);
    display(Lsad, outputFolderName()+"Lsad");

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");


    DVec rObs = tensorBasis()->getBasisElement(1,0);
    //rObs(0) = -2;
    //DVec rObs = DVec::Zero(sizeOfBasis());
    DMat rObsMat = reshape(rObs, nbOfQModes_, nbOfPModes_);
    display(rObsMat, outputFolderName()+"rObsMat");
    
    /*SMat qGramMat(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
    tensorBasis(0)->gramMatrix(qGramMat);
    display(qGramMat, outputFolderName()+"qGramMat");
    SMat gramMat = kron(qGramMat, SIdP_);*/
    
    SMat gramMat = tensorBasis_->gramMatrix();
    display(tensorBasis(0)->gramMatrix(), outputFolderName()+"qGramMat");
    DVec qMeans = tensorBasis(0)->basisMeans2();
    display(qMeans, outputFolderName()+"qMeans");
    DVec pMeans = tensorBasis(1)->basisMeans2();
    display(pMeans, outputFolderName()+"pMeans");
    DVec SR = gramMat * rObs;
    display(SR, outputFolderName()+"SR");
    
    DVec LinvRObs = solveWithSaddle(Leq(), gramMat * rObs);
    //cout << "solution vec" << endl << LinvRObs << endl;
    
    cout << "coucou" << endl;
    DMat LinvRObsMat = reshape(LinvRObs, nbOfQModes_, nbOfPModes_);
    cout << "solution mat" << endl << LinvRObsMat << endl;
    display(LinvRObsMat, outputFolderName()+"LinvRObsMat");
    
    DVec LLinvRObs = Leq() * LinvRObs;
    DMat LLinvRObsMat = reshape(LLinvRObs, nbOfQModes_, nbOfPModes_);
    cout << "L of solution mat" << endl << LLinvRObsMat << endl;
  }
  
  DVec ColloidLangevinGalerkin::CVObservable() const
  {
    return tensorBasis()->getBasisElement(1,0) / sqrt(tensorBasis(0)->omega());
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
    assert(iOfFourier2 < nbOfQModes_ && iOfHermite < nbOfPModes_ && iOfParticle < nbOfParticles_);
    return pow(nbOfQModes_ * nbOfPModes_ , iOfParticle) * ( nbOfQModes_ * iOfHermite + iOfFourier2);
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
    for (int iOfHermite = 1; iOfHermite < (int)nbOfPModes_; iOfHermite++)
      Lthm0_.insert(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
    cout << "end createLthm0" << endl;
  }

  void BoundaryLangevinGalerkin::createLham()
  {
    cout << "BoundaryLangevinGalerkin::createLham()" << endl;
    //Lham_ = arma::zeros(sizeOfBasis_, sizeOfBasis_);
    for (int i = 0; i < nbOfParticles_; i++)
    {
      Lham_ += doubleMatToTens(Q_, Pt_, i, i);
      if (i > 0) Lham_ -= doubleMatToTens(Q_, Pt_, i, i - 1);
      Lham_ -= doubleMatToTens(Qt_, P_, i, i);
      if (i < nbOfParticles_ - 1) Lham_ += doubleMatToTens(Q_, Pt_, i + 1, i);
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