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
    
    display(tensorBasis_->gramMatrix(), outputFolderName()+"gramMat");
    display(Leq_, outputFolderName()+"Leq");
    SMat Lsad = shapeSaddle(Leq_);
    display(Lsad, outputFolderName()+"Lsad");
    //cout << "Leq spectral gap : " << computeSpectralGap(Leq_) << endl;
    cout << "Lsad spectral gap : " << computeSpectralGap(Lsad) << endl;
  }
  
  DVec OverdampedGalerkin::CVcoeffsVec() const
  {
    //return basisMeans2();
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
      DMat DKeta = DMat(LtL) + basisMeans2(0) * basisMeans2(0).transpose();
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
      //(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
      //tensorBasis(0)->gramMatrix(gramMat);
      cout << "gramMat : " << gramMat.rows() << " " << gramMat.cols() << endl << gramMat << endl;
      cout << "Leq : " << Leq().rows() << " " << Leq().cols() << endl << Leq() << endl;
      cout << "CVObs : " << CVObservable().rows() << " " << CVObservable().cols() << endl << CVObservable() << endl;
      DVec LinvObs = solveWithSaddle(Leq(), gramMat * CVObservable());
      return LinvObs;
    }
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
    DMat DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;
    
    cout << "Norm2basisMeans2 : " << norm2meansVec(0) << endl;
    cout << "Norm of LeqinvSad : " << DLeqInv.norm() << endl;
    cout << "Norm of Leqinv : " << DLeq.inverse().norm() << endl;
    //cout << "Norm of Sinv Leqinv : " << shapePrec(DLeq).inverse().norm() << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat LeqInv = invWithSaddle(Leq_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "g   = " << basisMeans2(0) << endl << endl;
    cout << "L g = " << Leq_ * basisMeans2(0) << endl << endl;
    
    DVec gradV = CVObservable();
    display(gradV, "output/Galerkin/Overdamped/gradV");
    
    DVec LinvGradV = solve(Leq_, gradV);
    display(LinvGradV, "output/Galerkin/Overdamped/LinvGradV");
    
    DVec LinvGradVsad = solveWithSaddle(Leq_, gradV);
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
    tensorBasis_ = new HermiteHermiteBasis(input, *potential_);
    
    createOperators();
  }
  
  void ColloidOverdampedGalerkin::compute()
  {
    display(Leq_, outputFolderName()+"Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    DMat DLeqInv = invWithSaddle(DMat(Leq_));
    

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, outputFolderName()+"DLeqInv");


    DVec rObs = CVObservable();
    DMat rObsMat = reshape(rObs, nbOfQModes_, nbOfPModes_);
    display(rObsMat, outputFolderName()+"rObsMat");
    
    //SMat gramMat(nbOfQModes(), nbOfQModes());  // /!\ qGramMat should be dense !
    //tensorBasis(0)->gramMatrix(gramMat);
    SMat gramMat = tensorBasis()->gramMatrix();
    display(gramMat, outputFolderName()+"gramMat");
    
    DVec qMeans = tensorBasis(0)->basisMeans2();
    display(qMeans, outputFolderName()+"qMeans");
    
    cout << "gramMat : " << gramMat.rows() << " " << gramMat.cols() << endl << gramMat << endl;
    cout << "Leq : " << Leq().rows() << " " << Leq().cols() << endl << Leq() << endl;
    cout << "rObs : " << rObs.rows() << " " << rObs.cols() << endl << rObs << endl;
    
    /*DMat DGramMat = gramMat;
    JacobiSVD<DMat> svd(DGramMat, ComputeThinU | ComputeThinV);
    display(, outputFolderName()+"LinvRObsMat");*/
      
    DVec LinvRObs = solveWithSaddle(Leq(), gramMat * rObs);
    cout << "solution vec" << endl << LinvRObs << endl;
    
    cout << "coucou" << endl;
    DMat LinvRObsMat = reshape(LinvRObs, nbOfQModes_, nbOfPModes_);
    cout << "solution mat" << endl << LinvRObsMat << endl;
    display(LinvRObsMat, outputFolderName()+"LinvRObsMat");
  }
  
  ///
  /// Returns the coefficients of the sinus functions in the G_i * H_j basis
  DVec ColloidOverdampedGalerkin::CVObservable() const
  {
    //throw runtime_error("PeriodicOverdampedGalerkin::getGradV not implemented !");
    return tensorBasis()->getBasisElement(1,0);  
  }
}
