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

  
  

  DMat Galerkin::shapeSaddle(const DMat& A, const DMat& C) const
  {
    DMat Asad = DMat::Zero(A.rows() + C.cols(), A.cols() + C.cols());    
    Asad.topLeftCorner(A.rows(), A.cols()) = A.topLeftCorner(A.rows(), A.cols());
    Asad.rightCols(C.cols()) = C;
    Asad.bottomRows(C.cols()) = C.adjoint();

    /*DVec SU = tensorBasis()->gramMatrix() * tensorBasis()->basisMeans();
    cout << "SU :" << SU << endl << endl;
    for (int iOfElt = 0; iOfElt < sizeOfBasis(); iOfElt++)
    {
      Asad(A.rows(), iOfElt) = SU(iOfElt);
      Asad(iOfElt, A.cols()) = SU(iOfElt);
    }*/
    return Asad;
  }
  
  SMat Galerkin::shapeSaddle(const SMat& A, const DMat& C) const
  {
    //SMat Asad = SMat::Zero(A.rows() + 1, A.cols() + 1);
    SMat Asad(A.rows() + C.cols(), A.cols() + C.cols());
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
      
    //DVec SU = tensorBasis()->gramMatrix() * tensorBasis()->basisMeans();
    //cout << "SU :" << endl << SU << endl << endl;
    for (int iOfElt = 0; iOfElt < sizeOfBasis(); iOfElt++)
      for (int iOfCol = 0; iOfCol < C.cols(); iOfCol++)
        if (C(iOfElt, iOfCol) != 0)
        {
          coeffs.push_back(Trid(A.rows() + iOfCol, iOfElt, C(iOfElt, iOfCol)));
          coeffs.push_back(Trid(iOfElt, A.cols() + iOfCol, C(iOfElt, iOfCol)));
          //coeffs.push_back(Trid(A.rows(), iTens(iOfQModes, 0), basisMean(0, iOfQModes)));
          //coeffs.push_back(Trid(iTens(iOfQModes, 0), A.cols(), basisMean(0, iOfQModes)));
        }
    Asad.setFromTriplets(coeffs.begin(), coeffs.end());
    return Asad;
  }

  DMat Galerkin::unshapeSaddle(const DMat& Asad, const DMat& C) const
  { return Asad.topLeftCorner(Asad.rows() - C.cols(), Asad.cols() - C.cols()); }
  
  SMat Galerkin::unshapeSaddle(const SMat& Asad, const DMat& C) const
  { return Asad.topLeftCorner(Asad.rows() - C.cols(), Asad.cols() - C.cols()); }

  DVec Galerkin::shapeSaddle(const DVec& X, const DMat& C) const
  {
    DVec Xsad = DVec::Zero(X.size() + C.cols());
    Xsad.head(X.size()) = X;
    return Xsad;
  }

  DVec Galerkin::unshapeSaddle(const DVec& Xsad, const DMat& C) const
  {
    return Xsad.head(Xsad.size() - C.cols());
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
  
  DVec Galerkin::solveWithSaddle(SMat const& A, DVec const& Y, const DMat& C) const
  {
    //return solveWithSaddle(DMat(A), Y);
    cout << "Sparse solveWithSaddle..."; cout.flush();
    //cout << "Constraint column : " << endl << C << endl;
    DVec Ysad = shapeSaddle(Y, C);
    //cout << "Ysad = " << endl << Ysad << endl;
    SMat Asad = shapeSaddle(A, C);
    //cout << "Asad = " << endl << Asad << endl;
    ComplexEigenSolver<DMat> ces(A);
    cout << endl << "Sp A : " << ces.eigenvalues().real().adjoint() << endl;
    ComplexEigenSolver<DMat> ces2(Asad);
    cout << endl << "Sp Asad : " << ces2.eigenvalues().real().adjoint() << endl;
    DVec Xsad = solve(Asad,Ysad);  
    //cout << "Xsad = " << endl << Xsad << endl;
    
    //cout << "Done !" << endl;
    //cout << "Asad =" << endl << Asad << endl << endl << "Xsad =" << endl << Xsad << endl << endl <<"Ysad =" << endl << Ysad << endl << endl;
    //cout << "Asad * Xsad =" << Asad * Xsad << endl << endl;
    return unshapeSaddle(Xsad, C);
  }

  DVec Galerkin::solveWithSaddleAndGuess(SMat const& A, DVec const& Y, const DVec& X0, const DMat& C) const
  {
    //return solveWithSaddle(DMat(A), Y);
    cout << "Sparse solveWithSaddle..."; cout.flush();
    
    DVec Ysad = shapeSaddle(Y, C);
    cout << "Ysad = " << endl << Ysad << endl;
    SMat Asad = shapeSaddle(A, C);
    cout << "Asad = " << endl << Asad << endl;
    
    DVec X0sad = shapeSaddle(X0, C);
    DVec Xsad = solveWithGuess(Asad,Ysad, X0sad);
    
    return unshapeSaddle(Xsad, C);
  }

  DVec Galerkin::solveWithSaddle(const DMat& A, const DVec& Y, const DMat& C) const
  {
    cout << "Dense solveWithSaddle..."; cout.flush();
    cout << "Y -> " << Y.size() << endl;
    DVec Ysad = shapeSaddle(Y, C);
    cout << "Xsad -> " << Ysad.size() << endl;
    cout << "A -> " << A.rows() << "x" << A.cols() << endl;
    DMat Asad = shapeSaddle(A, C);
    cout << "Asad -> " << Asad.rows() << "x" << Asad.cols() << endl;
      
    DVec Xsad = Asad.fullPivLu().solve(Ysad);
     
    cout << "Xsad -> " << Xsad.size() << endl;
    cout << "OK ! (lambda = " << Xsad(Xsad.size() - 1) << ")" << endl;
    return unshapeSaddle(Xsad, C);
  }
  
  /*/// 
  /// Finds by SVD the eigenvectors of S corresponding to an eigenvalue > tol and pseudo solve AX=SY on the space spanned by these eigenvectors
  DVec Galerkin::pseudoSolve(const SMat& A, const DVec& Y, const DMat& C)
  {
    //JacobiSVD<DMat> svd(S, ComputeThinU | ComputeThinV);
    //DVec M = 
    
  }*/

  DMat Galerkin::invWithSaddle(const SMat& A, const DMat& C) const
  {
    DMat DA(A);
    return invWithSaddle(DA, C);
  }

  DMat Galerkin::invWithSaddle(const DMat& A, const DMat& C) const
  {
    cout << "invWithSaddle..."; cout.flush();
    DMat Asad = shapeSaddle(A, C);
    DMat Bsad = Asad.inverse();
    cout << "fin invWithSaddle..."; cout.flush();
    return unshapeSaddle(Bsad, C);
  }
  
  
  


  
  DVec Galerkin::projectionOrthoG(DVec const& X) const
  {
    return X - dot(X, basisMeans(0)) / norm2meansVec(0) * basisMeans(0);
  }

  using namespace Eigen;

  Galerkin::Galerkin(Input const& input):       //ex : [0:4]
    nbOfParticles_(input.nbOfParticles()),
    nbOfQModes_(input.nbOfQModes()),  //ex : 5
    nbOfPModes_(input.nbOfPModes()),
    sizeOfBasis_(nbOfQModes_ * nbOfPModes_),
    SIdQ_(nbOfQModes_, nbOfQModes_),
    SIdP_(nbOfPModes_, nbOfPModes_),
    DIdQ_(DMat::Identity(nbOfQModes_, nbOfQModes_)),
    DIdP_(DMat::Identity(nbOfPModes_, nbOfPModes_)),
    Q_(nbOfQModes_, nbOfQModes_),
    P_(nbOfPModes_, nbOfPModes_),
    Qt_(nbOfQModes_, nbOfQModes_),
    Pt_(nbOfPModes_, nbOfPModes_),
    Lthm0_(nbOfPModes_, nbOfPModes_),
    Lthm_(sizeOfBasis_, sizeOfBasis_),
    Lham_(sizeOfBasis_, sizeOfBasis_),
    Lrep_(sizeOfBasis_, sizeOfBasis_),
    Leq_(sizeOfBasis_, sizeOfBasis_),
    Leta_(sizeOfBasis_, sizeOfBasis_),
    beta_(input.beta()),
    gamma_(input.gamma()),
    amplitude_(input.amplitude()),
    nonEqAmplitude_(input.nonEqAmplitude()),
    doNonequilibrium_(input.doGalerkinNonequilibrium()),
    nbOfIntegrationNodes_(1000),
    //expFourierMeans_(2 * nbOfQModes_, 0),
    trigToExpMat_(DMat::Zero(nbOfQModes_, nbOfQModes_)),
    expToTrigMat_(DMat::Zero(nbOfQModes_, nbOfQModes_)),
    trigToExpTens_(sizeOfBasis_, sizeOfBasis_),
    expToTrigTens_(sizeOfBasis_, sizeOfBasis_),
    potential_(createPotential(input, input.galerkinPotentialName())),
    //tensorBasis_(nullptr),
    outputFolderName_(input.outputFolderName())
  {
    cout << endl << "Number of modes : " << nbOfQModes_ << " x " << nbOfPModes_ << endl;
    cout << "Output written in " << outputFolderName() << endl;
    SIdQ_.setIdentity();
    SIdP_.setIdentity();
    
    //cout << "GalerkinPotential : " << input.secondPotentialName() << endl;
    //cout << potential_->value(-1) << " " << potential_->value(0) << " " << potential_->value(1) << " " << endl;

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
  
  int Galerkin::nbOfQModes() const
  {
    return nbOfQModes_;
  }
  
  int Galerkin::nbOfPModes() const
  {
    return nbOfPModes_;
  }

  double Galerkin::gamma() const
  {
    return gamma_;
  }
  
  double Galerkin::basisMean(int iOfVariable, int iOfElt) const
  {
    return tensorBasis(iOfVariable)->basisMean(iOfElt);
  }
  
  const DVec& Galerkin::basisMeans(int iOfVariable) const
  {
    return tensorBasis(iOfVariable)->basisMeans();
  }
  
  double Galerkin::norm2meansVec(int iOfVariable) const
  {
    return tensorBasis(iOfVariable)->norm2meansVec();
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  int Galerkin::iTens(int iOfFourier2, int iOfHermite) const
  {
    assert(iOfFourier2 < nbOfQModes_ && iOfHermite < nbOfPModes_);
    return nbOfQModes_ * iOfHermite + iOfFourier2;
  }


 
  
  EigenSolver<MatrixXd> Galerkin::getEigenSolver() const
  {
    cout << "Getting eigen elements by a dense matrix method !" << endl;
    //DMat DLeq(Leq_);
    DMat DLeta(Leta_);
    return EigenSolver<MatrixXd>(DLeta);
  }
  
  shared_ptr<CVBasis> Galerkin::createCvBasis(Input const& input)
  {
    TensorBasis* tensBasis = tensorBasis_;
    shared_ptr<DVec> cvcoeffs = make_shared<DVec>(CVcoeffsVec());
    //return new IsolatesCVBasis(dynamic_cast<TensorBasis*>(tensorBasis_), make_shared<DVec>(CVcoeffsVec()));
    
    if (input.systemName() == "Colloid")
      return make_shared<ColloidCVBasis>(tensBasis, cvcoeffs);
    else
      return make_shared<IsolatedCVBasis>(tensBasis, cvcoeffs);
      
    /*if (input.doGalerkinCV())
    {
        if (input.dynamicsName() == "Overdamped")
        {
          if (input.potentialName() == "Sinusoidal") return new PeriodicOverdampedGalerkin(input);
          else if (input.potentialName() == "DoubleWell" || input.potentialName() == "Harmonic" || input.potentialName() == "TwoTypes") return new ColloidOverdampedGalerkin(input);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        else if (input.dynamicsName() == "Langevin")
        {
          if (input.potentialName() == "Sinusoidal") return new PeriodicLangevinGalerkin(input);
          else if (input.potentialName() == "DoubleWell" || input.potentialName() == "Harmonic") return new ColloidLangevinGalerkin(input);
          else throw runtime_error("This potential matches no Galerkin method !");
        }
        //else if (input.dynamicsName() == "BoundaryLangevin") return new BoundaryLangevinGalerkin(input);
        else throw runtime_error("This dynamics matches no Galerkin method !");
      }
    else
      return nullptr;
  }*/
  }
  
  ///
  /// Builds the columns corresponding to degenerate directions of the Galerkin space, plus the constant function
  DMat Galerkin::computeConstraints(double tol) const
  {
    DMat Sq = tensorBasis(0)->gramMatrix();
    ComplexEigenSolver<DMat> ces(Sq);
    DVec vapSq = ces.eigenvalues().real();
    DMat vepSq = ces.eigenvectors().real();
    cout << "Eigenvalues of Sq :" << vapSq.adjoint() << endl;
    //cout << "Eigenvectors of S :" << endl << vepS << endl;
    DMat Cq(DMat::Zero(nbOfQModes(), nbOfQModes()));

    int nbOfSmallQ = 0;    

    ofstream outSpGramMat(outputFolderName()+"spGramMatQ", ofstream::app);
    for (int iOfElt = 0; iOfElt < nbOfQModes(); iOfElt++)
    {
      outSpGramMat << nbOfQModes() << " " << vapSq(iOfElt) << endl;
      if (vapSq(iOfElt) < tol)
      {
        //cout << "Addeed col " << iOfElt << " as constraint n" << nbOfSmall << endl;
        Cq.col(nbOfSmallQ) = vepSq.col(iOfElt);
        nbOfSmallQ++;
      }
    }
    int nbOfSmall = nbOfSmallQ*nbOfPModes();
    DMat C(DMat::Zero(sizeOfBasis(), nbOfSmall+1));
    C.rightCols(nbOfSmall) = kron(Cq, DIdP_);
    C.col(0) = tensorBasis()->basisMeans();
    cout << nbOfSmallQ << " / " << nbOfQModes() << " Q modes deleted" << endl;
    return C;
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
    
    SMat Lsad = shapeSaddle(Leq(), SU_);
    //cout << Lsad << endl;
    
    int regOfq = 2;
    int regOfp = 2;
    DVec z(nbOfQModes()* nbOfPModes());       // y is weakly H^1/2+regOfq in q and H^1/2+regOfp in p
    for (int iOfFourier=0; iOfFourier < nbOfQModes(); iOfFourier++)
      for (int iOfHermite=0; iOfHermite < nbOfPModes(); iOfHermite++)
        z(iOfHermite*nbOfQModes() + iOfFourier) = (iOfFourier?pow(iOfFourier, -.5-regOfq):1) * (iOfHermite?pow(iOfHermite, -.5-regOfp/2.):1);
    
    //DVec z = gettGiHj(0,1); 
    
    string strRefLinvZ = outPath+"ref/refLinvZ"+doubleToString(gamma_)+".txt";
    string strRefZ = outPath+"ref/refZ"+doubleToString(gamma_)+".txt";
    string strRefSpGap = outPath+"ref/refSpGap.txt";
    
    if (doComputeRef)
    {
      DVec refLinvZ = solveWithSaddle(Leq(), z, SU_);
      double refSpGap = computeSpectralGap(Lsad);
      //DMat refLinvZMat = reshape(refLinvZ, nbOfQModes(), nbOfPModes());
      
      ofstream outRefZ(strRefZ);
      if (!outRefZ.is_open())
        throw runtime_error("Reference output file in strRefLinvZ cannot be open");
      ofstream outRefLinvZ(strRefLinvZ);      
      ofstream outRefSpGap(strRefSpGap, std::ofstream::app);
      // Closing the file is very important if it is read right after !
      outRefZ << setprecision(15); outRefZ << "# " << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << endl << z; outRefZ.close();
      outRefLinvZ << setprecision(15); outRefLinvZ << "# " << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << endl << refLinvZ; outRefLinvZ.close();
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
      int refNbOfQModes = dimensions[0];
      int refNbOfPModes = dimensions[1];
      
      cout << "Reference solution : " << refNbOfQModes << " x " << refNbOfPModes << endl;
      if (refNbOfQModes < nbOfQModes() || refNbOfPModes < nbOfPModes())
        throw runtime_error("The reference solution must be bigger !");
      
      DMat refLinvZMat = reshape(refLinvZ, refNbOfQModes, refNbOfPModes);
      DMat truncRefLinvZMat = refLinvZMat.topLeftCorner(nbOfQModes(), nbOfPModes());
      DVec truncRefLinvZ = reshape(truncRefLinvZMat, 1, nbOfQModes() * nbOfPModes());
      
      cout << "Starting to solve" << endl;
      DVec LinvZ = solveWithSaddleAndGuess(Leq(), z, truncRefLinvZ, SU_);
      
      DMat LinvZMat = reshape(LinvZ, nbOfQModes(), nbOfPModes());
      cout << "Solved !" << endl;
      
      DMat approxErrorMat = refLinvZMat; approxErrorMat.topLeftCorner(nbOfQModes(), nbOfPModes()) = DMat::Zero(nbOfQModes(), nbOfPModes()); // Xref - Pi Xref
      DMat consistErrorMat = LinvZMat - refLinvZMat.topLeftCorner(nbOfQModes(), nbOfPModes());   // X - Pi Xref
      
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
      outError << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << " " << approxError << " " << consistError << " " << totalError << " " << mobilityError << " " << spGap << " " << spGapError << endl;;
    } 
  }
  
  void Galerkin::studyColloidLangevinErrors(bool doComputeRef)
  {
    //string outPath = "output/Galerkin/Langevin/DoubleWell/";
    //string outPath = "output/Galerkin/Langevin/velocity/";
    
    
    cout << "gamma = " << gamma_ << endl;
    
    DMat C = computeConstraints(1e-6);
    SMat Lsad = shapeSaddle(Leq(), C);
    SMat S = tensorBasis_->gramMatrix();
    //cout << Lsad << endl;
    
    DVec Z = CVObservable();
    
    string strRefLinvZ = outputFolderName()+"ref/refLinvZ"+doubleToString(gamma_)+".txt";
    string strRefZ = outputFolderName()+"ref/refZ"+doubleToString(gamma_)+".txt";
    string strRefSpGap = outputFolderName()+"ref/refSpGap.txt";
    
    if (doComputeRef)
    {
      DVec refLinvZ = solveWithSaddle(Leq(), S*Z, C);
      //double refSpGap = computeSpectralGap(Lsad);
      //DMat refLinvZMat = reshape(refLinvZ, nbOfQModes(), nbOfPModes());
      DMat DLsad = Lsad;
      ComplexEigenSolver<DMat> ces(DLsad);
      DVec eigLsad = ces.eigenvalues().real();
      //cout << endl << "Sp A : " << ces.eigenvalues().real().adjoint() << endl;
      double refSpGap = fabs(eigLsad[0]);
      /*for (int iOfVap = 1; iOfVap < (int)eigLsad.size(); iOfVap++)
        if (fabs(eigLsad[iOfVap]) < refSpGap)
          refSpGap = eigLsad[iOfVap];*/
      
      ofstream outRefZ(strRefZ);
      if (!outRefZ.is_open())
        throw runtime_error("Reference output file in strRefLinvZ cannot be open");
      ofstream outRefLinvZ(strRefLinvZ);      
      ofstream outRefSpGap(strRefSpGap, std::ofstream::app);
      // Closing the file is very important if it is read right after !
      outRefZ << setprecision(15); outRefZ << "# " << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << endl << Z; outRefZ.close();
      outRefLinvZ << setprecision(15); outRefLinvZ << "# " << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << endl << refLinvZ; outRefLinvZ.close();
      outRefSpGap << setprecision(15); outRefSpGap << gamma_<< " " << refSpGap << " " << dot(refLinvZ, Z) << endl; outRefSpGap.close();
      
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
      int refNbOfQModes = dimensions[0];
      int refNbOfPModes = dimensions[1];
      
      cout << "Reference solution : " << refNbOfQModes << " x " << refNbOfPModes << endl;
      if (refNbOfQModes < nbOfQModes() || refNbOfPModes < nbOfPModes())
        throw runtime_error("The reference solution must be bigger !");
      
      DMat refLinvZMat = reshape(refLinvZ, refNbOfQModes, refNbOfPModes);
      DMat truncRefLinvZMat = refLinvZMat.topLeftCorner(nbOfQModes(), nbOfPModes());
      DVec truncRefLinvZ = reshape(truncRefLinvZMat, 1, nbOfQModes() * nbOfPModes());
      
      cout << "Starting to solve" << endl;
      DVec LinvZ = solveWithSaddleAndGuess(Leq(), Z, truncRefLinvZ, SU_);
      
      DMat LinvZMat = reshape(LinvZ, nbOfQModes(), nbOfPModes());
      cout << "Solved !" << endl;
      
      DMat approxErrorMat = refLinvZMat; approxErrorMat.topLeftCorner(nbOfQModes(), nbOfPModes()) = DMat::Zero(nbOfQModes(), nbOfPModes()); // Xref - Pi Xref
      DMat consistErrorMat = LinvZMat - refLinvZMat.topLeftCorner(nbOfQModes(), nbOfPModes());   // X - Pi Xref
      
      /*DVec truncLinvZ= DVec::Zero(refNbOfQModes*refNbOfPModes); 
      truncLinvZTrunc.head(nbOfQModes()* nbOfPModes()) = LinvZ;
      DVec consistImErrorMat = refL*/
      
      double approxError = approxErrorMat.norm();
      double consistError = consistErrorMat.norm();
      double totalError = sqrt(pow(approxError, 2) + pow(consistError, 2));
      double mobilityError = dot(refLinvZ, refZ) - dot(LinvZ, Z);
      cout << "Reference mobility : " << dot(refLinvZ, refZ) << endl;
      
      double spGap = refSpGap;
      //if (false)
        spGap = computeSpectralGap(Lsad);
      double spGapError = refSpGap - spGap;
      cout << setprecision(15);
      cout << "spGapError = " << refSpGap << " - " << spGap << " = " << spGapError << endl;

      
      ofstream outError(outputFolderName()+"errors.txt", std::ofstream::app);
      outError << gamma_ << " " << nbOfQModes() << " " << nbOfPModes() << " " << approxError << " " << consistError << " " << totalError << " " << mobilityError << " " << spGap << " " << spGapError << endl;;
    } 
  }


}



