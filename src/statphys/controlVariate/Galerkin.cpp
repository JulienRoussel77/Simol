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
    DVec SU = tensorBasis()->gramMatrix() * tensorBasis()->basisMeans();
    cout << "SU :" << SU << endl << endl;
    for (int iOfElt = 0; iOfElt < sizeOfBasis(); iOfElt++)
    {
      Asad(A.rows(), iOfElt) = SU(iOfElt);
      Asad(iOfElt, A.cols()) = SU(iOfElt);
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
      
    DVec SU = tensorBasis()->gramMatrix() * tensorBasis()->basisMeans();
    cout << "SU :" << endl << SU << endl << endl;
    for (int iOfElt = 0; iOfElt < sizeOfBasis(); iOfElt++)
    {
      coeffs.push_back(Trid(A.rows(), iOfElt, SU(iOfElt)));
      coeffs.push_back(Trid(iOfElt, A.cols(), SU(iOfElt)));
      //coeffs.push_back(Trid(A.rows(), iTens(iOfQModes, 0), basisMean2(0, iOfQModes)));
      //coeffs.push_back(Trid(iTens(iOfQModes, 0), A.cols(), basisMean2(0, iOfQModes)));
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
    //cout << "Done !" << endl;
    //cout << "Asad =" << endl << Asad << endl << endl << "Xsad =" << endl << Xsad << endl << endl <<"Ysad =" << endl << Ysad << endl << endl;
    //cout << "Asad * Xsad =" << Asad * Xsad << endl << endl;
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
    nbOfQModes_(input.nbOfQModes()),  //ex : 5
    nbOfPModes_(input.nbOfPModes()),
    sizeOfBasis_(pow(nbOfQModes_ * nbOfPModes_, nbOfParticles_)),
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
    externalForce_(input.externalForce()),
    doNonequilibrium_(input.doGalerkinNonequilibrium()),
    nbOfIntegrationNodes_(1000),
    //expFourierMeans_(2 * nbOfQModes_, 0),
    trigToExpMat_(DMat::Zero(nbOfQModes_, nbOfQModes_)),
    expToTrigMat_(DMat::Zero(nbOfQModes_, nbOfQModes_)),
    trigToExpTens_(sizeOfBasis_, sizeOfBasis_),
    expToTrigTens_(sizeOfBasis_, sizeOfBasis_),
    potential_(createPotential(input)),
    tensorBasis_(nullptr),
    outputFolderName_(input.outputFolderName())
  {
    assert(nbOfQModes_ % 2 == 1);
    cout << endl << "Number of modes : " << nbOfQModes_ << " x " << nbOfPModes_ << endl;
    cout << "Output written in " << outputFolderName() << endl;
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
      DVec refLinvZ = solveWithSaddle(Leq(), z);
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
      int refNbOfFourier = dimensions[0];
      int refNbOfHermite = dimensions[1];
      
      cout << "Reference solution : " << refNbOfFourier << " x " << refNbOfHermite << endl;
      if (refNbOfFourier < nbOfQModes() || refNbOfHermite < nbOfPModes())
        throw runtime_error("The reference solution must be bigger !");
      
      DMat refLinvZMat = reshape(refLinvZ, refNbOfFourier, refNbOfHermite);
      DMat truncRefLinvZMat = refLinvZMat.block(0,0,nbOfQModes(), nbOfPModes());
      DVec truncRefLinvZ = reshape(truncRefLinvZMat, 1, nbOfQModes() * nbOfPModes());
      
      cout << "Starting to solve" << endl;
      DVec LinvZ = solveWithSaddleAndGuess(Leq(), z, truncRefLinvZ);
      
      DMat LinvZMat = reshape(LinvZ, nbOfQModes(), nbOfPModes());
      cout << "Solved !" << endl;
      
      DMat approxErrorMat = refLinvZMat; approxErrorMat.block(0,0,nbOfQModes(), nbOfPModes()) = DMat::Zero(nbOfQModes(), nbOfPModes()); // Xref - Pi Xref
      DMat consistErrorMat = LinvZMat - refLinvZMat.block(0,0,nbOfQModes(), nbOfPModes());   // X - Pi Xref
      
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


}



