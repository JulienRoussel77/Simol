#include "simol/statphys/controlVariate/Galerkin.hpp"

/*using std::cout;
using std::endl;
using std::ostream;

#include <iomanip>
using std::setprecision;
using std::setw;

#include "simol/statphys/Tools.hpp"*/

namespace simol
{

  
  

  DMat Galerkin::shapeSaddle(const DMat& A) const
  {
    DMat Asad = DMat::Zero(A.rows() + 1, A.cols() + 1);    
    Asad.block(0, 0, A.rows(), A.cols()) = A.block(0, 0, A.rows(), A.cols());
    for (int iOfFourier = 0; iOfFourier < nbOfFourier(); iOfFourier++)
    {
      Asad(A.rows(), iTens(iOfFourier, 0)) = basisMean2(0, iOfFourier); //expFourierMeans(iOfFourier2);
      Asad(iTens(iOfFourier, 0), A.cols()) = basisMean2(0, iOfFourier); //expFourierMeans(iOfFourier2);
    }
    return Asad;
  }
  
  SMat Galerkin::shapeSaddle(const SMat& A) const
  {
    //SMat Asad = SMat::Zero(A.rows() + 1, A.cols() + 1);
    SMat Asad(A.rows() + 1, A.cols() + 1);
    //Asad.block(0, 0, A.rows(), A.cols()) = A.block(0, 0, A.rows(), A.cols());
    
    vector<Trid> coeffs;
    coeffs.reserve(A.nonZeros());
    for (int jOfA = 0; jOfA < (int)A.cols(); ++jOfA)
      for (SMat::InnerIterator it(A, jOfA); it; ++it)
      {
        int iOfA = it.row();
        double valOfA = it.value();
        coeffs.push_back(Trid(iOfA, jOfA, valOfA));
        //Asad.insert(iOfA, jOfA) = valOfA;
      }
    
    for (int iOfFourier = 0; iOfFourier < nbOfFourier(); iOfFourier++)
    {
      coeffs.push_back(Trid(A.rows(), iTens(iOfFourier, 0), basisMean2(0, iOfFourier)));
      coeffs.push_back(Trid(iTens(iOfFourier, 0), A.cols(), basisMean2(0, iOfFourier)));
    }
    Asad.setFromTriplets(coeffs.begin(), coeffs.end());
    return Asad;
  }

  DMat Galerkin::unshapeSaddle(const DMat& Asad) const
  { return Asad.block(0, 0, Asad.rows() - 1, Asad.cols() - 1); }
  
  SMat Galerkin::unshapeSaddle(const SMat& Asad) const
  { return Asad.block(0, 0, Asad.rows() - 1, Asad.cols() - 1); }

  DVec Galerkin::shapeSaddle(const DVec& X) const
  {
    DVec Xsad = DVec::Zero(X.size() + 1);
    Xsad.head(X.size()) = X;
    return Xsad;
  }

  DVec Galerkin::unshapeSaddle(const DVec& Xsad) const
  {
    return Xsad.head(Xsad.size() - 1);
  }
  
  DVec Galerkin::solve(SMat const& A, DVec const& Y) const
  {
    SMat Acomp = A;
    Acomp.makeCompressed();
    Eigen::SparseLU<SMat > solver(Acomp);
    /*solver.compute(Acomp);
    if(solver.info()!=Eigen::Success) {
      // decomposition failed
      cout << "failed !" << endl;
    }*/
    DVec X = solver.solve(Y);
    return X;
    
    //Eigen::BiCGSTAB<SMat> solver(A);
    //solver.compute(A);
    //DVec X = solver.solve(Y);
    //if(solver.info() != Eigen::Success) throw runtime_error("System resolution failed !");
    //return X;
  }
  
  DVec Galerkin::solveWithGuess(SMat const& A, DVec const& Y, DVec const& X0) const
  {
    //SMat Acomp = A;
    //Acomp.makeCompressed();
    Eigen::GMRES<SMat, Eigen::IncompleteLUT<double>> solver(A);
    DVec X = solver.solveWithGuess(Y, X0);
    if (solver.info() != Eigen::Success)
      throw runtime_error("solveWithGuess didnt succeed !");
    DVec residual = A*X-Y;
    cout << "Residual = " << residual.norm() << endl << endl;
    return X;
  }
  
  
  DVec Galerkin::solve(const DMat& A, const DVec& X) const
  {    
    return A.fullPivLu().solve(X);
  }
  
  DVec Galerkin::solveWithSaddle(SMat const& A, DVec const& Y) const
  {
    //return solveWithSaddle(DMat(A), Y);
    cout << "Sparse solveWithSaddle..."; cout.flush();
    
    DVec Ysad = shapeSaddle(Y);
    SMat Asad = shapeSaddle(A);
    DVec Xsad = solve(Asad,Ysad);  
    cout << "Done !" << endl;
    return unshapeSaddle(Xsad);
  }

  DVec Galerkin::solveWithSaddleAndGuess(SMat const& A, DVec const& Y, const DVec& X0) const
  {
    //return solveWithSaddle(DMat(A), Y);
    cout << "Sparse solveWithSaddle..."; cout.flush();
    
    DVec Ysad = shapeSaddle(Y);
    SMat Asad = shapeSaddle(A);
    //cout << Asad << endl;
    DVec X0sad = shapeSaddle(X0);
    DVec Xsad = solveWithGuess(Asad,Ysad, X0sad);
    
    return unshapeSaddle(Xsad);
  }

  DVec Galerkin::solveWithSaddle(const DMat& A, const DVec& Y) const
  {
    cout << "Dense solveWithSaddle..."; cout.flush();
    cout << "Y -> " << Y.size() << endl;
    DVec Ysad = shapeSaddle(Y);
    cout << "Xsad -> " << Ysad.size() << endl;
    cout << "A -> " << A.rows() << "x" << A.cols() << endl;
    DMat Asad = shapeSaddle(A);
    cout << "Asad -> " << Asad.rows() << "x" << Asad.cols() << endl;
      
    DVec Xsad = Asad.fullPivLu().solve(Ysad);
     
    cout << "Xsad -> " << Xsad.size() << endl;
    cout << "OK ! (lambda = " << Xsad(Xsad.size() - 1) << ")" << endl;
    return unshapeSaddle(Xsad);
  }

  DMat Galerkin::invWithSaddle(const SMat& A) const
  {
    DMat DA(A);
    return invWithSaddle(DA);
  }

  DMat Galerkin::invWithSaddle(const DMat& A) const
  {
    cout << "invWithSaddle..."; cout.flush();
    DMat Asad = shapeSaddle(A);
    DMat Bsad = Asad.inverse();
    cout << "fin invWithSaddle..."; cout.flush();
    return unshapeSaddle(Bsad);
  }

  
  
  


  
  DVec Galerkin::projectionOrthoG(DVec const& X) const
  {
    return X - dot(X, basisMeans2(0)) / norm2meansVec(0) * basisMeans2(0);
  }

  using namespace Eigen;

  Galerkin::Galerkin(Input const& input):       //ex : [0:4]
    nbOfParticles_(input.nbOfParticles()),
    nbOfFourier_(input.nbOfFourier()),  //ex : 5
    nbOfHermite_(input.nbOfHermite()),
    maxOfFourier_((nbOfFourier_ + 1) / 2),  //ex : 3
    sizeOfBasis_(pow(nbOfFourier_ * nbOfHermite_, nbOfParticles_)),
    SIdQ_(nbOfFourier_, nbOfFourier_),
    SIdP_(nbOfHermite_, nbOfHermite_),
    DIdQ_(DMat::Identity(nbOfFourier_, nbOfFourier_)),
    DIdP_(DMat::Identity(nbOfHermite_, nbOfHermite_)),
    Q_(nbOfFourier_, nbOfFourier_),
    P_(nbOfHermite_, nbOfHermite_),
    tQ_(nbOfFourier_, nbOfFourier_),
    tP_(nbOfHermite_, nbOfHermite_),
    Lthm0_(nbOfHermite_, nbOfHermite_),
    Lthm_(sizeOfBasis_, sizeOfBasis_),
    Lham_(sizeOfBasis_, sizeOfBasis_),
    Lrep_(sizeOfBasis_, sizeOfBasis_),
    Leq_(sizeOfBasis_, sizeOfBasis_),
    Leta_(sizeOfBasis_, sizeOfBasis_),
    beta_(input.beta()),
    gamma_(input.gamma()),
    amplitude_(input.amplitude()),
    externalForce_(input.externalForce()),
    doNonequilibrium_(input.doGalerkinNonequilibrium()),
    nbOfIntegrationNodes_(1000),
    //expFourierMeans_(2 * nbOfFourier_, 0),
    trigToExpMat_(DMat::Zero(nbOfFourier_, nbOfFourier_)),
    expToTrigMat_(DMat::Zero(nbOfFourier_, nbOfFourier_)),
    trigToExpTens_(sizeOfBasis_, sizeOfBasis_),
    expToTrigTens_(sizeOfBasis_, sizeOfBasis_),
    potential_(createPotential(input)),
    tensorBasis_(nullptr),
    outputFolderName_(input.outputFolderName())
  {
    assert(nbOfFourier_ % 2 == 1);
    cout << endl << "Number of modes : " << nbOfFourier_ << " x " << nbOfHermite_ << endl;
    
    SIdQ_.setIdentity();
    SIdP_.setIdentity();

    // Computation of the passage matrix
    /*computeExpToTrigMat();
    createQ();
    createP();*/
  }

  Galerkin::~Galerkin()
  {
    delete potential_;
  }
  
  int Galerkin::nbOfVariables() const
  {
    return 2 * nbOfParticles_;
  }

  int Galerkin::nbOfParticles() const
  {
    return nbOfParticles_;
  }
  
  int Galerkin::sizeOfBasis() const
  {
    return sizeOfBasis_;
  }
  
  int Galerkin::nbOfFourier() const
  {
    return nbOfFourier_;
  }
  
  int Galerkin::nbOfHermite() const
  {
    return nbOfHermite_;
  }

  double Galerkin::gamma() const
  {
    return gamma_;
  }
  
  double Galerkin::basisMean2(int iOfVariable, int iOfElt) const
  {
    return tensorBasis(iOfVariable)->basisMean2(iOfElt);
  }
  
  const DVec& Galerkin::basisMeans2(int iOfVariable) const
  {
    return tensorBasis(iOfVariable)->basisMeans2();
  }
  
  double Galerkin::norm2meansVec(int iOfVariable) const
  {
    return tensorBasis(iOfVariable)->norm2meansVec();
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  int Galerkin::iTens(int iOfFourier2, int iOfHermite) const
  {
    assert(iOfFourier2 < nbOfFourier_ && iOfHermite < nbOfHermite_);
    return nbOfFourier_ * iOfHermite + iOfFourier2;
  }


 
  
  EigenSolver<MatrixXd> Galerkin::getEigenSolver() const
  {
    cout << "Getting eigen elements by a dense matrix method !" << endl;
    //DMat DLeq(Leq_);
    DMat DLeta(Leta_);
    return EigenSolver<MatrixXd>(DLeta);
  }
  
  /*///Returns te coefficients of tG_i * Hj in the tG_i * H_j basis
  ///tG designs the Fourier functions (V=0)
  DVec Galerkin::gettGiHjTrig(int i, int j) const
  {
    DVec GiHjTrig = DVec::Zero(sizeOfBasis_);
    GiHjTrig(iTens(i, j)) = 1;
    return GiHjTrig;
  }

  ///Returns te coefficients of tG_i * Hj in the G_i * H_j basis
  ///tG designs the Fourier functions (V=0)
  DVec Galerkin::gettGiHj(int i, int j) const
  {
    
    return trigToExpTens_ * gettGiHjTrig(i,j);
  }



  DVec Galerkin::getLtGiHj(int i, int j) const
  { return Leta_ * gettGiHj(i, j); }

  DVec Galerkin::getLtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * getLtGiHj(i, j); }

  DVec Galerkin::getLinvtGiHj(int i, int j) const
  { return solveWithSaddle(Leta_, gettGiHj(i, j)); }

  DVec Galerkin::getLinvtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * getLinvtGiHj(i,j); }*/

  
  CVBasis Galerkin::makeCvBasis()
  {
    return CVBasis(dynamic_cast<TensorBasis*>(tensorBasis_), make_shared<DVec>(CVcoeffsVec()));
  }
  
  void Galerkin::computeEigen() const
  {    
    Eigen::EigenSolver<DMat> eigSol = getEigenSolver();
    //const cplx* eigVal = eigSol.eigenvalues().data();
    vector<cplx> eigVal;
    for (int i=0; i<sizeOfBasis(); i++)
      eigVal.push_back(eigSol.eigenvalues().data()[i]);
    
    std::sort(eigVal.begin(), eigVal.end(), hasSmallerNorm);
    
    ofstream out_eigvalLeq("output/Galerkin/eigvalLeq");
    for (int i = 0; i < (int) eigVal.size(); i++)
      out_eigvalLeq << real(eigVal[i]) << " " << imag(eigVal[i]) << endl;
    
    cout << "computeEigen ok" << endl;
  }
  
  void Galerkin::studyLangevinErrors(bool doComputeRef)
  {
    string outPath = "output/Galerkin/Langevin/H2/";
    //string outPath = "output/Galerkin/Langevin/velocity/";
    
    
    cout << "gamma = " << gamma_ << endl;
    
    SMat Lsad = shapeSaddle(Leq());
    //cout << Lsad << endl;
    
    int regOfq = 2;
    int regOfp = 2;
    DVec z(nbOfFourier()* nbOfHermite());       // y is weakly H^1/2+regOfq in q and H^1/2+regOfp in p
    for (int iOfFourier=0; iOfFourier < nbOfFourier(); iOfFourier++)
      for (int iOfHermite=0; iOfHermite < nbOfHermite(); iOfHermite++)
        z(iOfHermite*nbOfFourier() + iOfFourier) = (iOfFourier?pow(iOfFourier, -.5-regOfq):1) * (iOfHermite?pow(iOfHermite, -.5-regOfp/2.):1);
    
    //DVec z = gettGiHj(0,1); 
    
    string strRefLinvZ = outPath+"ref/refLinvZ"+doubleToString(gamma_)+".txt";
    string strRefZ = outPath+"ref/refZ"+doubleToString(gamma_)+".txt";
    string strRefSpGap = outPath+"ref/refSpGap.txt";
    
    if (doComputeRef)
    {
      DVec refLinvZ = solveWithSaddle(Leq(), z);
      double refSpGap = computeSpectralGap(Lsad);
      //DMat refLinvZMat = reshape(refLinvZ, nbOfFourier(), nbOfHermite());
      
      ofstream outRefZ(strRefZ);
      if (!outRefZ.is_open())
        throw runtime_error("Reference output file in strRefLinvZ cannot be open");
      ofstream outRefLinvZ(strRefLinvZ);      
      ofstream outRefSpGap(strRefSpGap, std::ofstream::app);
      // Closing the file is very important if it is read right after !
      outRefZ << setprecision(15); outRefZ << "# " << gamma_ << " " << nbOfFourier() << " " << nbOfHermite() << endl << z; outRefZ.close();
      outRefLinvZ << setprecision(15); outRefLinvZ << "# " << gamma_ << " " << nbOfFourier() << " " << nbOfHermite() << endl << refLinvZ; outRefLinvZ.close();
      outRefSpGap << setprecision(15); outRefSpGap << gamma_<< " " << refSpGap << " " << dot(refLinvZ, z) << endl; outRefSpGap.close();
      
      cout << "Reference written in " << strRefLinvZ << endl;
    }
    else
    {    

      
      vector<int> dimensions;
      
      cout << "Reading reference in " << strRefZ << endl;
      DVec refZ = scanTensor(strRefZ, dimensions);
      cout << "Reading reference in " << strRefLinvZ << endl;
      DVec refLinvZ = scanTensor(strRefLinvZ, dimensions);
      /*ifstream inSpGap(strRefSpGap);
      if (!inSpGap.is_open())
        throw runtime_error(strRefSpGap+" is not a valid path !");*/
      cout << "Reading reference in " << strRefSpGap << endl;
      map<double, double> gammaToSpGap = scanMap(strRefSpGap);
      //double refSpGap = readItem(inSpGap);      
      if (gammaToSpGap.find(gamma_) == gammaToSpGap.end())
        throw runtime_error("gamma = "+doubleToString(gamma_)+" not found in refSpGap !");
      double refSpGap = gammaToSpGap[gamma_];
      int refNbOfFourier = dimensions[0];
      int refNbOfHermite = dimensions[1];
      
      cout << "Reference solution : " << refNbOfFourier << " x " << refNbOfHermite << endl;
      if (refNbOfFourier < nbOfFourier() || refNbOfHermite < nbOfHermite())
        throw runtime_error("The reference solution must be bigger !");
      
      DMat refLinvZMat = reshape(refLinvZ, refNbOfFourier, refNbOfHermite);
      DMat truncRefLinvZMat = refLinvZMat.block(0,0,nbOfFourier(), nbOfHermite());
      DVec truncRefLinvZ = reshape(truncRefLinvZMat, 1, nbOfFourier() * nbOfHermite());
      
      cout << "Starting to solve" << endl;
      DVec LinvZ = solveWithSaddleAndGuess(Leq(), z, truncRefLinvZ);
      
      DMat LinvZMat = reshape(LinvZ, nbOfFourier(), nbOfHermite());
      cout << "Solved !" << endl;
      
      DMat approxErrorMat = refLinvZMat; approxErrorMat.block(0,0,nbOfFourier(), nbOfHermite()) = DMat::Zero(nbOfFourier(), nbOfHermite()); // Xref - Pi Xref
      DMat consistErrorMat = LinvZMat - refLinvZMat.block(0,0,nbOfFourier(), nbOfHermite());   // X - Pi Xref
      
      double approxError = approxErrorMat.norm();
      double consistError = consistErrorMat.norm();
      double totalError = sqrt(pow(approxError, 2) + pow(consistError, 2));
      double mobilityError = dot(refLinvZ, refZ) - dot(LinvZ, z);
      cout << "Reference mobility : " << dot(refLinvZ, refZ) << endl;
      
      double spGap = refSpGap;
      //if (false)
        spGap = computeSpectralGap(Lsad);
      double spGapError = refSpGap - spGap;
      cout << setprecision(15);
      cout << "spGapError = " << refSpGap << " - " << spGap << " = " << spGapError << endl;

      
      ofstream outError(outPath+"errors.txt", std::ofstream::app);
      outError << gamma_ << " " << nbOfFourier() << " " << nbOfHermite() << " " << approxError << " " << consistError << " " << totalError << " " << mobilityError << " " << spGap << " " << spGapError << endl;;
    } 
  }

  
  
  
  //########## OverdampedGalerkin #########
  
  OverdampedGalerkin::OverdampedGalerkin(Input const& input):
    Galerkin(input)
  {    
    tensorBasis_ = new ExpFourierHermiteBasis(input, *potential_);
    
    //computeExpToTrigMat();
    //computeExpToTrigTens();
    
    tensorBasis(0)->gradMatrix(Q_);
    tQ_ = Q_.adjoint();
    
    createLeta();
  }
  
  ///
  /// The matrix Leq, Leta are positive !
  void OverdampedGalerkin::createLeta()
  {
    tensorBasis(0)->laplacianMatrix(Leq_);
    Leq_ /= beta_;
    
    tensorBasis(0)->gradMatrix(Lrep_);
    Leta_ = Leq_ - externalForce_ * Lrep_;
    /*for (int iOfHermite = 0; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      for (int jOfHermite = 0; jOfHermite < (int)nbOfHermite_; jOfHermite++)
        Leq_(iOfHermite, iOfHermite) = - 1 / beta_ * basis(0)->xLaplacianY(iOfHermite+1, jOfHermite+1);*/
  }
  
  /*void OverdampedGalerkin::computeExpToTrigTens()
  {
    cout << "Suboptimal implementation of OverdampedGalerkin::computeExpToTrigTens !!!" << endl;
    trigToExpTens_ = trigToExpMat_.sparseView();
    expToTrigTens_ = expToTrigMat_.sparseView();
  }*/
  
  void OverdampedGalerkin::compute()
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
    
    DVec gradV = getGradV();
    display(gradV, "output/Galerkin/Overdamped/gradV");
    
    DVec LinvGradV = solve(Leq_, gradV);
    display(LinvGradV, "output/Galerkin/Overdamped/LinvGradV");
    
    DVec LinvGradVsad = solveWithSaddle(Leq_, gradV);
    display(LinvGradVsad, "output/Galerkin/Overdamped/LinvGradVsad");
  }
  
  ///
  /// Returns the coefficients of the sinus functions in the G_i * H_j basis
  DVec OverdampedGalerkin::getGradV() const
  { 
    throw runtime_error("OverdampedGalerkin::getGradV not implemented !");
    //return amplitude_ * gettGiHj(1,0) / sqrt(2.);    // /!\ the basis elements are normalized in L2, not Linfty, so   sin = gettGiHj(1,0) / sqrt(2.)
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
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      
      ///TODO : use a matrix free approach to avoid dense matrices
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
      //return -projectionOrthoG(solve(Keta, -Leta_.transpose() * getGradV()));
    }
    else
    {
      cout << "Computing the control variate by an equilibrium CG" << endl;
      //SMat Keq = Leq_.transpose() * Leq_;
      //return -projectionOrthoG(solve(Keq, -Leq_.transpose() * getGradV()));
      Eigen::ConjugateGradient<SMat> solver(Leq_);
      return -solver.solve(-projectionOrthoG(getGradV()));
    }
  }


  



}



