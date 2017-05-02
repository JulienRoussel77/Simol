#include "simol/statphys/controlVariate/OverdampedGalerkin.hpp"


namespace simol
{
 
  //########## OverdampedGalerkin #########
  
  OverdampedGalerkin::OverdampedGalerkin(Input const& input):
    Galerkin(input)
  {}
  
  ///
  /// The matrix Leq, Leta are positive !
  void OverdampedGalerkin::createOperators()
  {
    Leq_ = tensorBasis(0)->laplacianMatrix();
    Leq_ /= beta_;
    
    Lrep_ = tensorBasis(0)->gradMatrix();
    Leta_ = Leq_ - externalForce_ * Lrep_;
    
    // U is the componants of 1 in the basis and SU is the scalar products of 1 with these elements
    SU_ = tensorBasis()->basisMeans();
    
    display(tensorBasis_->gramMatrix(), outputFolderName()+"gramMat");
    display(Lrep_, outputFolderName()+"Q");
    display(Leq_, outputFolderName()+"Leq");
    
    SMat Lsad = shapeSaddle(Leq_, SU_);
    display(Lsad, outputFolderName()+"Lsad");
    //cout << "Leq spectral gap : " << computeSpectralGap(Leq_) << endl;
    //cout << "Lsad spectral gap : " << computeSpectralGap(Lsad) << endl;
  }
  
  DVec OverdampedGalerkin::CVcoeffsVec() const
  {
    //return basisMeans();
    //return getGradV();
    
    /*if (doNonequilibrium())
      return -projectionOrthoG(solve(Leta_, -projectionOrthoG(getGradV())));
    else
      return -projectionOrthoG(solve(Leq_, -projectionOrthoG(getGradV())));*/
        
    if (doNonequilibrium())
    {
      throw runtime_error("Nonequilibrium Galerkin not implemented !");
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      
      /*///TODO : use a matrix free approach to avoid dense matrices
      SMat LtL = Leta_.transpose() * Leta_;
      DMat DKeta = DMat(LtL) + basisMeans(0) * basisMeans(0).transpose();
      SMat Keta = DKeta.sparseView();
      //Eigen::LeastSquaresConjugateGradient<SMat> solver(Leta_);
      Eigen::ConjugateGradient<SMat> solver0(Keta);
      //solver0.setTolerance(1e-12);
      DVec CV0 = -solver0.solve(-Leta_.transpose() * getGradV());
      cout << "CG : " << solver0.info() << " " << solver0.error() << " " << solver0.iterations()<< endl << Keta * CV0 - Leta_.transpose() * getGradV() << endl;
      Eigen::BiCGSTAB<SMat> solver(Keta);
      //solver.setTolerance(1e-12);
      DVec CV = -solver.solve(-Leta_.transpose() * getGradV());
      cout << "BC : " << solver.info() << " " << solver.error() << " " << solver0.iterations() << endl <<Keta * CV - Leta_.transpose() * getGradV() << endl;
      return CV;
      //return -projectionOrthoG(solve(Keta, -Leta_.transpose() * getGradV()));*/
    }
    else
    {
      SMat gramMat = tensorBasis()->gramMatrix();
      
      //cout << "gramMat : " << gramMat.rows() << " " << gramMat.cols() << endl << gramMat << endl;
      //cout << "Leq : " << Leq().rows() << " " << Leq().cols() << endl << Leq() << endl;
      cout << "CVObs : " << CVObservable().rows() << " " << CVObservable().cols() << endl << CVObservable() << endl;
      cout << "--> norm = " << CVObservable().norm() << endl;
      cout << "Starting solveWithSaddle" << endl;
      DVec LinvObs = solveWithSaddle(Leq(), gramMat * CVObservable(), SU_);
      cout << "LinvObs : " << endl << LinvObs << endl;
      //return CVObservable();
      return LinvObs;
    }
  }
  
  ///
  ///Finds the eigenvectors of S corresponding to eigenvectors > tol and solves AX=Y on the spanned subspace
  DVec OverdampedGalerkin::solveResilient(const SMat& A, const DVec& Y) const
  {    
    DMat S = tensorBasis()->gramMatrix();
    ComplexEigenSolver<DMat> ces(S);
    DVec vapS = ces.eigenvalues().real();
    DMat vepS = ces.eigenvectors().real();
    cout << "Eigenvalues of S :" << vapS.adjoint() << endl;
    //cout << "Eigenvectors of S :" << endl << vepS << endl;
    DMat C(DMat::Zero(sizeOfBasis(), sizeOfBasis()));
    C.col(0) = tensorBasis_->basisMeans();
    double tol = 1e-6;
    int nbOfSmall = 1;    

    ofstream outSpGramMat(outputFolderName()+"spGramMat", ofstream::app);
    for (int iOfElt = 0; iOfElt < sizeOfBasis(); iOfElt++)
    {
      outSpGramMat << sizeOfBasis() << " " << vapS(iOfElt) << endl;
      if (vapS(iOfElt) < tol)
      {
        //cout << "Addeed col " << iOfElt << " as constraint n" << nbOfSmall << endl;
        C.col(nbOfSmall) = vepS.col(iOfElt);
        nbOfSmall++;
      }
    }
    // Becareful the block function does not work "in place"
    DMat Ctemp = C.leftCols(nbOfSmall); 
    C = Ctemp;
    //cout << "C = " << endl << C << endl;
    //cout << "S = " << endl << S << endl;
    DVec X = solveWithSaddle(A, S*Y, C);
    cout << "X = " << endl << X << endl;
    return X;
  }
  
  // #### PeriodicOverdampedGalerkin ####
  
  PeriodicOverdampedGalerkin::PeriodicOverdampedGalerkin(Input const& input):
    OverdampedGalerkin(input)
  {    
    tensorBasis_ = new ExpFourierHermiteBasis(input, *potential_);
    
    createOperators();
  }
  
  void PeriodicOverdampedGalerkin::compute()
  {
    cout << "start OverdampedGalerkin::compute()" << endl;
    cout << "############ Leq ############" << endl;
    //cout << Leq_ << endl << endl;
    
    cout << "Leq_ size : " << Leq_.rows() << " x " << Leq_.cols() << endl;
    cout << "Leq_ nnz : " << Leq_.nonZeros() << endl;

    DMat DLeq(Leq_);
    display(Leq_, "output/Galerkin/Overdamped/Leq");


    //DMat DLeqSad = shapeSaddle(DLeq);
    //display(DLeqSad, "output/Galerkin/Overdamped/LeqSad");

    cout << "############ DLeqSad ############" << endl;
    //cout << DLeqSad << endl << endl;

    //DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);

    cout << "Computing DLeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat DLeqInv = invWithSaddle(DLeq, SU_);
    cout << "OK" << endl;
    
    cout << "Norm2basisMeans : " << norm2meansVec(0) << endl;
    cout << "Norm of LeqinvSad : " << DLeqInv.norm() << endl;
    cout << "Norm of Leqinv : " << DLeq.inverse().norm() << endl;
    //cout << "Norm of Sinv Leqinv : " << shapePrec(DLeq).inverse().norm() << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat LeqInv = invWithSaddle(Leq_, SU_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "g   = " << basisMeans(0) << endl << endl;
    cout << "L g = " << Leq_ * basisMeans(0) << endl << endl;
    
    DVec gradV = CVObservable();
    display(gradV, "output/Galerkin/Overdamped/gradV");
    
    DVec LinvGradV = solve(Leq_, gradV);
    display(LinvGradV, "output/Galerkin/Overdamped/LinvGradV");
    
    DVec LinvGradVsad = solveWithSaddle(Leq_, gradV, SU_);
    display(LinvGradVsad, "output/Galerkin/Overdamped/LinvGradVsad");
  }
  
  ///
  /// Returns the coefficients of the sinus functions in the G_i * H_j basis
  DVec PeriodicOverdampedGalerkin::CVObservable() const
  {
    throw runtime_error("PeriodicOverdampedGalerkin::getGradV not implemented !");
    //return amplitude_ * gettGiHj(1,0) / sqrt(2.);    
  }
  
  
  // #### ColloidOverdampedGalerkin ####
  
  ColloidOverdampedGalerkin::ColloidOverdampedGalerkin(Input const& input):
    OverdampedGalerkin(input)
  {    
    //tensorBasis_ = new HermiteHermiteBasis(input, *potential_);
    tensorBasis_ = new ExpHermiteHermiteBasis(input, *potential_);

    createOperators();
  }
  
  void ColloidOverdampedGalerkin::compute()
  {
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    DMat DLeqInv = invWithSaddle(DMat(Leq_), SU_);
    

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");


    DVec rObs = CVObservable();
    DMat rObsMat = reshape(rObs, nbOfQModes_, nbOfPModes_);
    display(rObsMat, outputFolderName()+"rObsMat");
    
    //SMat gramMat(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
    //tensorBasis(0)->gramMatrix(gramMat);
    SMat gramMat = tensorBasis()->gramMatrix();
    display(gramMat, outputFolderName()+"gramMat");
    
    DVec qMeans = tensorBasis(0)->basisMeans();
    display(qMeans, outputFolderName()+"qMeans");
    
    cout << "gramMat : " << gramMat.rows() << " " << gramMat.cols() << endl << gramMat << endl;
    cout << "Leq : " << Leq().rows() << " " << Leq().cols() << endl << Leq() << endl;
    cout << "rObs : " << rObs.rows() << " " << rObs.cols() << endl << rObs << endl;
    DVec tU = basisMeans(0);
    cout << "tU : " << endl << tU << endl;
    //DMat gramMatZero = (DMat)gramMat - U * U.adjoint();
    
    DVec LinvRObs = solveResilient(Leq(), rObs);
    
    DMat LinvRObsMat = reshape(LinvRObs, nbOfQModes_, nbOfPModes_);
    display(LinvRObsMat, outputFolderName()+"LinvRObsMat");
  }
  
  ///
  /// Returns the coefficients of the function q,p -> q in the G_i * H_j basis
  DVec ColloidOverdampedGalerkin::CVObservable() const
  {
    //throw runtime_error("PeriodicOverdampedGalerkin::getGradV not implemented !");
    //return tensorBasis()->getBasisElement(1,0);  
    return tensorBasis(0)->getMonome1();
  }
}
