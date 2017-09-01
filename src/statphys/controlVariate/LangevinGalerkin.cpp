#include "simol/statphys/controlVariate/LangevinGalerkin.hpp"

//#### LangevinGalerkin ####

namespace simol
{
  LangevinGalerkin::LangevinGalerkin(Input const& input, shared_ptr<CVBasis> cvBasis0):
    Galerkin(input, cvBasis0)
  {}
  
  ///
  ///Should be called in all the constructors of the daughters of LangevinGalerkin
  /// /!\ By convention for the operators the variables are in the order (p,q) and for the states the variables are (q,p)
  /// /!\ L_ij = < e_i, -calL e_j > is a positive matrix !
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
    
    SMat PQt = kron(P_, Qt_);
    SMat PtQ = PQt.transpose();
    Lham_ = PQt - PtQ;   // The matrices and the operators have opposite signs
    createLthm();
    display(Lthm0_, outputFolderName()+"Lthm0");
    Leq_  = Lham_ + gamma_ * Lthm_;
    Lrep_ = -kron(P_, SIdQ_);
    Leta_ = Leq_ + nonEqAmplitude_ * Lrep_;
    display(Leq_, outputFolderName()+"Leq");
    
    SU_ = tensorBasis()->basisMeans();
    cout << "Leq : " << Leq_.rows() << " " << Leq_.cols() << endl;
    SMat Lsad = shapeSaddle(Leq_, SU_);
    display(Lsad, outputFolderName()+"Lsad");
    cout << "Lsad : " << Lsad.rows() << " " << Lsad.cols() << endl;
    display(tensorBasis_->gramMatrix(), outputFolderName()+"gramMat");
    display(tensorBasis(0)->gramMatrix(), outputFolderName()+"qGramMat");
    display(tensorBasis(1)->gramMatrix(), outputFolderName()+"pGramMat");
    
    //cout << "Leq spectral gap : " << computeSpectralGap(Leq_) << endl;
    //cout << "Lsad spectral gap : " << computeSpectralGap(Lsad) << endl;
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
      /*throw runtime_error("Nonequilibrium Galerkin not implemented !");
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      SMat Keta = Leta_.transpose() * Leta_;
      //DVec CV = -solve(Keta, Leta_.transpose() * gettGiHj(0,1));
      DVec CV = -solve(Keta, Leta_.transpose() * CVObservable());
      return CV;*/
    }
    else
    {
      /*DVec rObs = CVObservable();
      DVec LinvRObs = solveResilient(Leq(), rObs);
      display(LinvRObs, outputFolderName()+"LinvRObs.txt");
      return LinvRObs;*/
      
      SMat gramMat = tensorBasis()->gramMatrix();
      cout << "CVObs : " << CVObservable().rows() << " " << CVObservable().cols() << endl << CVObservable() << endl;
      cout << "--> norm = " << CVObservable().norm() << endl;
      cout << "--> obs iid variance = " << pow(CVObservable().norm(), 2) << endl;
      cout << "Starting solveWithSaddle" << endl;
      DVec LinvObs = solveWithSaddle(Leq(), gramMat * CVObservable(), SU_);
      cout << "LinvObs : " << endl << LinvObs << endl;
      cout << "--> LinvObs norm = " << LinvObs.norm() << endl;
      cout << "--> obs asy variance = " << 2*dot(CVObservable(), LinvObs) << endl;
      
      return LinvObs;
    }
  }
  
  
  
  ///
  ///Finds the eigenvectors of S corresponding to eigenvectors > tol and solves AX=Y on the spanned subspace
  DVec LangevinGalerkin::solveResilient(const SMat& A, const DVec& Y) const
  {    
    double tol = 1e-6;
    DMat C = computeConstraints(tol);
    SMat S = tensorBasis()->gramMatrix();
    DVec X = solveWithSaddle(A, S*Y, C);
    cout << "X = " << endl << reshape(X, nbOfQModes(), nbOfPModes()) << endl;
    return X;
  }
  
  //### PeriodicLangevinGalerkin ###
  
  PeriodicLangevinGalerkin::PeriodicLangevinGalerkin(Input const& input, shared_ptr<CVBasis> cvBasis0):
    LangevinGalerkin(input, cvBasis0)
  {
    //tensorBasis_ = new ExpFourierHermiteBasis(input, *potential_);
    createOperators();
    computeCVBasisCoeffs();
  }
  
  void PeriodicLangevinGalerkin::compute()
  {

    //DMat DLeq(Leq_);
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    //DMat LeqInv = invWithSaddle(DLeq);
    //DVec U = tensorBasis()->basisMeans();
    DMat DLeqInv = invWithSaddle(DMat(Leq_), SU_);
    

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");
    
    DVec H1Trig = DVec::Zero(sizeOfBasis_);
    H1Trig(iTens(0, 1)) = 1;
    
    //Map<DMat, Eigen::Aligned> H1TrigMat(H1Trig, nbOfQModes_, nbOfPModes_);
    DMat H1TrigMat = reshape(H1Trig, nbOfQModes_, nbOfPModes_);

    DVec H1 = tensorBasis()->getPartialElement(1,1);

    cout << "############ H1Mat ############" << endl;
    display(H1, outputFolderName()+"H1");
    //DMat H1Mat = reshape(gettGiHj(0,1), nbOfQModes_, nbOfPModes_);
    DMat H1Mat = reshape(H1, nbOfQModes_, nbOfPModes_);
    display(H1Mat, outputFolderName()+"H1Mat");

    display(H1Trig, outputFolderName()+"H1Trig");
    display(H1TrigMat, outputFolderName()+"H1TrigMat");
    
    DVec LinvH1 = solveWithSaddle(Leq(), H1, SU_);
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

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, Lrep_ * LinvH1, SU_);
    //double varCoeff = -.5 * dot( Lrep_ * LinvH1, LeqInv * Lrep_ * LinvH1);
    double varCoeff = 2 * dot( Lrep_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }
  
  DVec PeriodicLangevinGalerkin::CVObservable() const
  {
    return tensorBasis()->getPartialElement(1,1);
  }
  

  
  //### ColloidLangevinGalerkin
  
  ColloidLangevinGalerkin::ColloidLangevinGalerkin(Input const& input, shared_ptr<CVBasis> cvBasis0):
    LangevinGalerkin(input, cvBasis0)
  {
    //tensorBasis_ = new HermiteHermiteBasis(input, *potential_);
    createOperators();
    computeCVBasisCoeffs();
  }
  
  ///
  ///
  void ColloidLangevinGalerkin::compute()
  {
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    DVec SU = tensorBasis()->gramMatrix() * tensorBasis()->basisMeans();
    DMat DLeqInv = invWithSaddle(DMat(Leq_), SU_);
    
    DMat Lsad = shapeSaddle(Leq_, SU_);
    display(Lsad, outputFolderName()+"Lsad");

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");


    DVec rObs = CVObservable();
    //rObs(0) = -2;
    //DVec rObs = DVec::Zero(sizeOfBasis());
    DMat rObsMat = reshape(rObs, nbOfQModes_, nbOfPModes_);
    display(rObsMat, outputFolderName()+"rObsMat");
    
    /*SMat qGramMat(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
    tensorBasis(0)->gramMatrix(qGramMat);
    display(qGramMat, outputFolderName()+"qGramMat");
    SMat gramMat = kron(qGramMat, SIdP_);*/
    
    cout << "rObs : " << endl << reshape(rObs, nbOfQModes(), nbOfPModes()) << endl;
    DVec LinvRObs = solveResilient(Leq(), rObs);
    
    //DVec Xexact = kron(tensorBasis(1)->getMonome0(), tensorBasis(0)->getMonome1());
    
    /*SMat gramMat = tensorBasis_->gramMatrix();
    display(tensorBasis(0)->gramMatrix(), outputFolderName()+"qGramMat");
    DVec qMeans = tensorBasis(0)->basisMeans();
    display(qMeans, outputFolderName()+"qMeans");
    DVec pMeans = tensorBasis(1)->basisMeans();
    display(pMeans, outputFolderName()+"pMeans");
    DVec SR = gramMat * rObs;
    display(SR, outputFolderName()+"SR");
    
    DVec LinvRObs = solveWithSaddle(Leq(), gramMat * rObs, SU);
    //cout << "solution vec" << endl << LinvRObs << endl;*/
    
    DMat LinvRObsMat = reshape(LinvRObs, nbOfQModes_, nbOfPModes_);
    display(LinvRObsMat, outputFolderName()+"LinvRObsMat");
    
    DVec LLinvRObs = Leq() * LinvRObs;
    DMat LLinvRObsMat = reshape(LLinvRObs, nbOfQModes_, nbOfPModes_);
    cout << "L of solution mat" << endl << LLinvRObsMat << endl;
  }
  
  DVec ColloidLangevinGalerkin::CVObservable() const
  {
    cout << "ColloidLangevinGalerkin::CVObservable" << endl;
    //return tensorBasis()->getBasisElement(1,0) / sqrt(tensorBasis(0)->omega());
    return kron(tensorBasis(1)->getMonome0(), tensorBasis(0)->getMonome1());
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
}