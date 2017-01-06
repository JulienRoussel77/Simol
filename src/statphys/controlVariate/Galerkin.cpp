#include "simol/statphys/controlVariate/Galerkin.hpp"

using std::cout;
using std::endl;
using std::ostream;

#include <iomanip>
using std::setprecision;
using std::setw;

#include "simol/statphys/Tools.hpp"

namespace simol
{
  Galerkin* createGalerkin(Input const& input)
  {
    if (input.doGalerkinCV())
    {
        if (input.dynamicsName() == "Overdamped") return new OverdampedGalerkin(input);
        else if (input.dynamicsName() == "Langevin") return new LangevinGalerkin(input);
        else if (input.dynamicsName() == "BoundaryLangevin") return new BoundaryLangevinGalerkin(input);
        else throw runtime_error("This dynamics matches no Galerkin method !");
      }
    else
      return nullptr;
  }
  
  

  DMat Galerkin::shapeSaddle(const DMat& A) const
  {
    DMat Asad = DMat::Zero(A.rows() + 1, A.cols() + 1);    
    Asad.block(0, 0, A.rows(), A.cols()) = A.block(0, 0, A.rows(), A.cols());
    for (int iOfFourier = 0; iOfFourier < nbOfFourier(); iOfFourier++)
    {
      Asad(A.rows(), iTens(iOfFourier, 0)) = gVector(iOfFourier); //expFourierMeans(iOfFourier2);
      Asad(iTens(iOfFourier, 0), A.cols()) = gVector(iOfFourier); //expFourierMeans(iOfFourier2);
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
    
    /*for (int iOfFourier = 0; iOfFourier < nbOfFourier(); iOfFourier++)
      for (int iOfHermite = 0; iOfHermite < nbOfHermite(); iOfHermite++)
        Asad.insert(iOfFourier, iOfHermite) = A(iOfFourier, iOfHermite);*/
    for (int iOfFourier = 0; iOfFourier < nbOfFourier(); iOfFourier++)
    {
      coeffs.push_back(Trid(A.rows(), iTens(iOfFourier, 0), gVector(iOfFourier)));
      coeffs.push_back(Trid(iTens(iOfFourier, 0), A.cols(), gVector(iOfFourier)));
      //Asad.insert(A.rows(), iTens(iOfFourier, 0)) = gVector(iOfFourier); //expFourierMeans(iOfFourier);
      //Asad.insert(iTens(iOfFourier, 0), A.cols()) = gVector(iOfFourier); //expFourierMeans(iOfFourier);
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

  /// C = A \otimes B = ( A_11 B    A_12 B   ...
  ///                     A_21 B    A_22 B   ...
  ///                                            )
  /// Note that (A \otimes B) (x \otimes y) = (Bx) \otimes (Ay)
  SMat kron(const SMat& A, const SMat& B)
  {
    double totalTime = clock();
    
    
    vector<Trid> Ccoeffs;
    Ccoeffs.reserve(A.nonZeros() * B.nonZeros());
    for (int jOfA = 0; jOfA < (int)A.cols(); ++jOfA)
    {
      for (SMat::InnerIterator it(A, jOfA); it; ++it)
      {
        int iOfA = it.row();
        double valOfA = it.value();
        for (int jOfB = 0; jOfB < (int)B.cols(); jOfB++)
        {
          for (SMat::InnerIterator it2(B, jOfB); it2; ++it2)
          {
            int iOfB = it2.row();
            double valOfB = it2.value();
            Ccoeffs.push_back(Trid(B.rows() * iOfA + iOfB, B.cols() * jOfA + jOfB, valOfA * valOfB));
          }
        }
      }
    }
    
    SMat C(A.rows()*B.rows(), A.cols()*B.cols());
    C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

    cout << "- - Sparse kron : " << clock() - totalTime << endl;
    return C;
  }

  DMat kron(const DMat& A, const DMat& B)
  {
    double totalTime = clock();
    
    DMat C(A.rows()*B.rows(), A.cols()*B.cols());
    //C = kroneckerProduct(A,B);
    for (int iOfA = 0; iOfA < (int) A.rows(); iOfA++)
      for (int jOfA = 0; jOfA < (int) A.cols(); jOfA++)
        C.block(B.rows() * iOfA, B.cols() * jOfA, B.rows(), B.cols()) = A(iOfA, jOfA) * B;
          
    cout << "- - Dense kron : " << clock() - totalTime << endl;
    return C;
  }
  
  SMat kron(const SMat& A, const DMat& B)
  {
    double totalTime = clock(); 
    
    // Fastest and no memory issue !
    vector<Trid> Ccoeffs;
    Ccoeffs.reserve(A.nonZeros() * B.nonZeros());
    for (int jOfA = 0; jOfA < (int)A.cols(); jOfA++)
    {
      for (SMat::InnerIterator it2(A, jOfA); it2; ++it2)
      {
        int iOfA = it2.row();
        double valOfA = it2.value();
        for (int iOfB = 0; iOfB < (int) B.rows(); iOfB++)
          for (int jOfB = 0; jOfB < (int) B.cols(); jOfB++)
            Ccoeffs.push_back(Trid(B.rows() * iOfA + iOfB, B.cols() * jOfA + jOfB, valOfA * B(iOfB, jOfB)));
      }
    }
    SMat SC(A.rows()*B.rows(), A.cols()*B.cols());
    SC.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    
    cout << "- - Sparse Dense kron : " << clock() - totalTime << endl;
    return SC;
  }
  
  SMat kron(const DMat& A, const SMat& B)
  {
    double totalTime = clock();
    
    vector<Trid> Ccoeffs;
    Ccoeffs.reserve(A.nonZeros() * B.nonZeros());
    for (int iOfA = 0; iOfA < (int) A.rows(); iOfA++)
       for (int jOfA = 0; jOfA < (int) A.cols(); jOfA++)
    for (int jOfB = 0; jOfB < (int)B.cols(); jOfB++)
      for (SMat::InnerIterator it2(B, jOfB); it2; ++it2)
      {
        int iOfB = it2.row();
        double valOfB = it2.value();
        Ccoeffs.push_back(Trid(B.rows() * iOfA + iOfB, B.cols() * jOfA + jOfB, A(iOfA, jOfA) * valOfB));
      }
    SMat SC(A.rows()*B.rows(), A.cols()*B.cols());
    SC.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    
    cout << "- - Dense Sparse kron : " << clock() - totalTime << endl;
    return SC;
  }
  
  ///
  /// Becareful by default the preconditioner is diagonal, and not working for our matrices !!
  /// Becareful X = solver.solve(X) is not valid !!
  double computeSpectralGap(SMat const& A)
  {
    cout << "computeSpectralGap" << endl;
     
    double tol = 1e-15;
    int nbOfIter = 0;
    //Eigen::SparseLU<SMat> solver(A);
    Eigen::GMRES<SMat, Eigen::IncompleteLUT<double>> solver(A); //solver.setTolerance(1e-15);
    //Eigen::GMRES<SMat, Eigen::IdentityPreconditioner> solver(A);
    //Eigen::BiCGSTAB<SMat, Eigen::IncompleteLUT<double>> solver(A);

    DVec X = DVec::Random(A.rows());
    DVec Y;
    double eigVal = X.norm();
    double eigValDiff = 1;
    double prevEigValDiff = 2;
    
    while (fabs(eigValDiff) > tol && fabs(eigValDiff) < fabs(prevEigValDiff))
    {
      Y = X / eigVal;
      //X = solver.solve(Y);
      X = solver.solveWithGuess(Y, Y * eigVal);
      //std::cout << "#iterations:     " << solver.iterations() << std::endl;
      //std::cout << "estimated error: " << solver.error()      << std::endl;
      DVec Xdiff = X - Y*eigVal;
      prevEigValDiff = eigValDiff;
      eigValDiff = X.norm() - eigVal;
      eigVal = X.norm();   // The sign is due to the positive sign of the matrix A !
      nbOfIter++;
      
      cout << nbOfIter << " : eigVal = " << eigVal << ", " << eigValDiff << " > " << tol << " / XdiffNorm = " << Xdiff.norm() << endl;
    }
    if (fabs(eigValDiff) > tol)
      cout << "##################### !! Simulation not converged !! #########################" << endl;
    cout << "Lanczos algo : " << nbOfIter << " iterations" << endl;
    return 1/eigVal;
  }
  
  

  void displayCplx(const DVecCplx& X, ostream& out)
  {
    for (int i = 0; i < (int) X.size(); i++)
      out << real(X(i)) << " " << imag(X(i)) << endl;
  }

  void display(const DVec& A, string path)
  {
    ofstream out(path);
    //display(A, out);
    out << A;
  }

  void display(const DMat& A, ostream& out)
  {
    for (int i = 0; i < (int) A.rows(); i++)
    {
      for (int j = 0; j < (int) A.cols(); j++)
      {
        if (true)//fabs(A(i,j)) > 1e-15)
        {
          //cout << i << " " << j << " " << A(i,j) << endl;
          out << setw(12) << A(i, j) << " ";
        }
        else
          out << setw(12) << "nan ";
      }
      out << endl;
    }
  }

  void display(const DMat& A, string path)
  {
    ofstream out(path);
    assert(out.is_open());
    display(A, out);
  }

  void display(const SMat& A, ostream& out)
  {
    DMat DA(A);;
    display(DA, out);
  }

  void display(const SMat& A, string path)
  {
    DMat DA(A);;
    display(DA, path);
  }
  
  

  //We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  //where C = ExpFourierBasis::basisCoefficient_
  void Galerkin::computeExpToTrigMat()
  {
    for (int iOfFourier2 = 1; iOfFourier2 < nbOfFourier(); iOfFourier2++)
    {
      trigToExpMat_(0, iOfFourier2) = expFourierMeans(iOfFourier2) / sqrt(2.);
      trigToExpMat_(iOfFourier2, 0) = expFourierMeans(iOfFourier2) / sqrt(2.);
    }
    trigToExpMat_(0, 0) = expFourierMeans(0) / 2;

    for (int iOfFourier = 1; iOfFourier < maxOfFourier_; iOfFourier++)
      for (int jOfFourier = 1; jOfFourier < maxOfFourier_; jOfFourier++)
      {
        //cosine times cosine
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier) =
          expFourierMeans(2 * (iOfFourier + jOfFourier)) / 2
          + expFourierMeans(2 * abs(iOfFourier - jOfFourier)) / 2;

        //sinus times sinus
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier - 1) =
          - expFourierMeans(2 * (iOfFourier + jOfFourier)) / 2
          + expFourierMeans(2 * abs(iOfFourier - jOfFourier)) / 2;

        int eps = (iOfFourier >= jOfFourier) - (iOfFourier <= jOfFourier); // -1 if i smaller, 0 if equal, 1 if i larger
        //cosine times sinus   /!\ doing the test in this order avoids to read expFourierMeans(-1)
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier - 1) =
          expFourierMeans(2 * (iOfFourier + jOfFourier) - 1) / 2
          + (!eps)?0:eps * expFourierMeans(2 * abs(iOfFourier - jOfFourier) - 1) / 2;
          
        //sinus times cosine
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier) =
          expFourierMeans(2 * (iOfFourier + jOfFourier) - 1) / 2
          - (!eps)?0:eps * expFourierMeans(2 * abs(iOfFourier - jOfFourier) - 1) / 2;
      }


    ofstream out_trigToExpMat("output/Galerkin/trigToExpMat");
    display(trigToExpMat_, out_trigToExpMat);

    expToTrigMat_  = trigToExpMat_.inverse();
  }
  
  DVec Galerkin::projectionOrthoG(DVec const& X) const
  {
    return X - dot(X, gVector()) / norm2gVector() * gVector();
  }

  DMat Galerkin::convertToTrigBasis(const DMat& X)
  {
    return expToTrigMat_ * X;
  }

  void Galerkin::createQ()
  {
    basis_(0)->gradMatrix(Q_);

    tQ_ = Q_.adjoint();
  }

  void Galerkin::createP()
  {
    for (int iOfHermite = 1; iOfHermite < nbOfHermite_; iOfHermite++)
      P_.insert(iOfHermite - 1, iOfHermite) = sqrt(beta_ * iOfHermite);

    tP_ = P_.adjoint();
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
    basis_(input, *potential_)
  {
    assert(nbOfFourier_ % 2 == 1);
    cout << endl << "Number of modes : " << nbOfFourier_ << " x " << nbOfHermite_ << endl;
    
    SIdQ_.setIdentity();
    SIdP_.setIdentity();

    //computeFourierCoeffsExp();
    // Computation of the passage matrix
    double totalTime = clock();
    computeExpToTrigMat();
    createQ();
    createP();
    cout << "- Galerkin::Galerkin : " << clock() - totalTime << endl;
    //display(P_, "output/Galerkin/P");
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

  const double& Galerkin::gamma() const
  {
    return gamma_;
  }
  
  double const& Galerkin::expFourierMeans(int iOfElt) const
  {
    return basis_.expFourierMeans(iOfElt);
  }
  
  const DVec& Galerkin::gVector() const
  {
    return basis_.gVector();
  }
  
  double Galerkin::gVector(int iOfElt) const
  {
    return basis_.gVector(iOfElt);
  }
  
  const double& Galerkin::norm2gVector() const
  {
    return basis_.norm2gVector();
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
  
  ///Returns te coefficients of tG_i * Hj in the tG_i * H_j basis
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
    /*DVec GiHjTrig = DVec::Zero(sizeOfBasis_);
    GiHjTrig(iTens(i, j)) = 1;
    DVec GiHj = trigToExpTens_ * GiHjTrig;
    return GiHj;*/
    
    return trigToExpTens_ * gettGiHjTrig(i,j);
  }



  DVec Galerkin::getLtGiHj(int i, int j) const
  { return Leta_ * gettGiHj(i, j); }

  DVec Galerkin::getLtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * getLtGiHj(i, j); }

  DVec Galerkin::getLinvtGiHj(int i, int j) const
  { return solveWithSaddle(Leta_, gettGiHj(i, j)); }

  DVec Galerkin::getLinvtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * getLinvtGiHj(i,j); }

  /*SMat Galerkin::CVcoeffs() const
  {
    DVec vecCoeffs = getLinvtGiHj(0, 1);

    return SMat(vecCoeffs, vecCoeffs.size(), 1);
  }*/
  

  
  CVBasis Galerkin::makeCvBasis()
  {
    return CVBasis(dynamic_cast<TensorBasis*>(&basis_), make_shared<DVec>(CVcoeffsVec()));
  }
  
  void Galerkin::computeEigen() const
  {    
    Eigen::EigenSolver<DMat> eigSol = getEigenSolver();
    //const cplx* eigVal = eigSol.eigenvalues().data();
    vector<cplx> eigVal;
    for (int i=0; i<sizeOfBasis(); i++)
      eigVal.push_back(eigSol.eigenvalues().data()[i]);
    
    std::sort(eigVal.begin(), eigVal.end(), hasSmallerNorm);
    
    //for (int i=0; i<sizeOfBasis(); i++)
    //  cout << eigVal[i] << endl;
    
    
    /*ofstream outFinalGap("output/Galerkin/finalGap.txt",  std::ofstream::app);
    //outFinalGap << externalForce_ << " " <<  amplitude_ << " " << gamma() << " " << nbOfFourier() << " " << nbOfHermite() << " " << abs(eigVal[0]) << " " << abs(eigVal[1]) << endl;
    DVec gradV = amplitude_ * gettGiHj(2,0) / sqrt(2.);
    cout << "gradV : " << gradV << endl << endl;
    DVec gradVz = gradV - dot(gradV, gVector()) * gVector();
    cout << "gradVz : " << gradVz << endl;
   
    DVec CV = solve(Leq_, gradV);
    cout << "CV before proj : " << CV << endl;
    DVec CVzero = CV - dot(CV, gVector()) * gVector();
    cout << "CV after proj : " << CVzero << endl;
    DVec LCV = Leq_*CV;
    DVec LCVz = Leq_*CVzero;
    cout << "LCV : " << LCV << endl;
    cout << "LCVz : " << LCVz << endl;
    DVec CVneq = CV - solve(Leq_, gradV);
    outFinalGap << externalForce_ << " " <<  amplitude_ << " " << gamma() << " " << nbOfFourier() << " " << nbOfHermite() << " " << abs(eigVal[0]) << " " << abs(eigVal[1]) << " " << CVneq.norm() << endl;*/
    
    ofstream out_eigvalLeq("output/Galerkin/eigvalLeq");
    for (int i = 0; i < (int) eigVal.size(); i++)
      out_eigvalLeq << real(eigVal[i]) << " " << imag(eigVal[i]) << endl;
    
    cout << "computeEigen ok" << endl;

    //displayCplx(eigVal, out_eigvalLeq);
  }
  


  
  
  
  //########## OverdampedGalerkin #########
  
  OverdampedGalerkin::OverdampedGalerkin(Input const& input):
    Galerkin(input)
  {    
    computeExpToTrigTens();
    createLeta();
  }
  
  ///
  /// The matrix Leq, Leta are positive !
  void OverdampedGalerkin::createLeta()
  {
    basis_(0)->laplacianMatrix(Leq_);
    Leq_ /= beta_;
    
    basis_(0)->gradMatrix(Lrep_);
    Leta_ = Leq_ - externalForce_ * Lrep_;
    /*for (int iOfHermite = 0; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      for (int jOfHermite = 0; jOfHermite < (int)nbOfHermite_; jOfHermite++)
        Leq_(iOfHermite, iOfHermite) = - 1 / beta_ * basis_(0)->xLaplacianY(iOfHermite+1, jOfHermite+1);*/
  }
  
  void OverdampedGalerkin::computeExpToTrigTens()
  {
    cout << "Suboptimal implementation of OverdampedGalerkin::computeExpToTrigTens !!!" << endl;
    trigToExpTens_ = trigToExpMat_.sparseView();
    expToTrigTens_ = expToTrigMat_.sparseView();
  }
  
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
    
    cout << "Norm2gVector : " << norm2gVector() << endl;
    cout << "Norm of LeqinvSad : " << DLeqInv.norm() << endl;
    cout << "Norm of Leqinv : " << DLeq.inverse().norm() << endl;
    //cout << "Norm of Sinv Leqinv : " << shapePrec(DLeq).inverse().norm() << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat LeqInv = invWithSaddle(Leq_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/Overdamped/DLeqInv");
    
    cout << "g   = " << gVector() << endl << endl;
    cout << "L g = " << Leq_ * gVector() << endl << endl;
    
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
    return amplitude_ * gettGiHj(1,0) / sqrt(2.);    // /!\ the basis elements are normalized in L2, not Linfty, so   sin = gettGiHj(1,0) / sqrt(2.)
  }
  
  DVec OverdampedGalerkin::CVcoeffsVec() const
  {
    //return gVector();
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
      DMat DKeta = DMat(LtL) + gVector() * gVector().transpose();
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


  //#### LangevinGalerkin ####


  LangevinGalerkin::LangevinGalerkin(Input const& input):
    Galerkin(input)
  {
    double totalTime = clock();
    computeExpToTrigTens();
    //cout << "- computeExp : " << clock() - totalTime << endl; totalTime = clock();
    SMat PQt = kron(tP_, Q_);
    SMat PtQ = PQt.transpose();
    Lham_ = -PQt + PtQ;   // The matricies and the operators have opposite signs
    //Lham_ = kron(Q_, tP_) - kron(tQ_, P_);
    //display(Lham_, "output/Galerkin/Langevin/Lham");
    createLthm();
    //cout << "- Lham and Lthm : " << clock() - totalTime << endl; totalTime = clock();
    //display(Lthm_, "output/Galerkin/Langevin/Lthm");
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
    basis_(1)->laplacianMatrix(Lthm0_);
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

  void LangevinGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = kron(SIdP_, trigToExpMat_);
    expToTrigTens_ = kron(SIdP_, expToTrigMat_);
  }
  
  void LangevinGalerkin::compute()
  {

    //DMat DLeq(Leq_);
    display(Leq_, "output/Galerkin/Langevin/Leq");

    cout << "Computing DLeqInv..."; cout.flush();
    //DMat LeqInv = invWithSaddle(DLeq);
    DMat DLeqInv = invWithSaddle(DMat(Leq_));
    

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, "output/Galerkin/Langevin/DLeqInv");


    DVec H1Trig = DVec::Zero(sizeOfBasis_);
    H1Trig(iTens(0, 1)) = 1;

    //Map<DMat, Eigen::Aligned> H1TrigMat(H1Trig, nbOfFourier_, nbOfHermite_);
    DMat H1TrigMat = reshape(H1Trig, nbOfFourier_, nbOfHermite_);

    //DMat H1Mat = trigToExpMat_ * H1TrigMat;
    //DVec H1 = reshape(H1Mat, nbOfFourier_ * nbOfHermite_, 1);

    DVec H1 = trigToExpTens_ * H1Trig;


    /*DVec H1(sizeOfBasis_);
    for (int iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
      H1(iTens(iOfFourier2, 1)) = expFourierMeans_[iOfFourier2];*/

    cout << "############ H1Mat ############" << endl;
    display(gettGiHj(0, 1), "output/Galerkin/H1");
    //DMat H1Mat = reshape(gettGiHj(0,1), nbOfFourier_, nbOfHermite_);
    DMat H1Mat = reshape(gettGiHj(0, 1), nbOfFourier_, nbOfHermite_);
    display(H1Mat, "output/Galerkin/Langevin/H1Mat");

    display(H1Trig, "output/Galerkin/Langevin/H1Trig");
    display(H1TrigMat, "output/Galerkin/Langevin/H1TrigMat");


    display(gettGiHj(0, 0), "output/Galerkin/Langevin/G0Trig");

    DMat G0Mat = reshape(gettGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(G0Mat, "output/Galerkin/Langevin/G0Mat");
    DMat LG0Mat = reshape(getLtGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(LG0Mat, "output/Galerkin/Langevin/LG0Mat");

    display(gettGiHj(1, 2), "output/Galerkin/Langevin/G1H2Trig");
    display(getLtGiHj(1, 2), "output/Galerkin/Langevin/LG1H2Trig");

    display(getLtGiHj(0, 1), "output/Galerkin/Langevin/LH1");
    DVec LinvH1 = getLinvtGiHj(0, 1);
    display(LinvH1, "output/Galerkin/Langevin/LinvH1");

    //double lambda = LinvH1Sad(sizeOfBasis_);

    //cout << "############ lambda ############" << endl;
    //cout << lambda << endl << endl;

    DMat LinvH1Mat = reshape(LinvH1, nbOfFourier_, nbOfHermite_);

    cout << "############ LinvH1Mat ############" << endl;
    display(LinvH1Mat, "output/Galerkin/Langevin/LinvH1Mat");

    DMat LinvH1MatTrig = convertToTrigBasis(LinvH1Mat);
    display(LinvH1MatTrig, "output/Galerkin/Langevin/LinvH1MatTrig");

    DVec H1back = Leq_ * LinvH1;
    DMat H1backMat = reshape(H1back, nbOfFourier_, nbOfHermite_);
    display(H1backMat, "output/Galerkin/Langevin/H1backMat");

    /*cx_vec eigvalLeq;
    cx_mat eigvecLeq;
    eig_gen(eigvalLeq, eigvecLeq, -DLeq);*/

    //DMat DLeqW = Leq_.wrapped_;
    EigenSolver<MatrixXd> eigSol(DMat(Leq_));

    /*/eigs_gen(eigvalLeq, eigvecLeq, Leq_, 20);

    ofstream out_eigvalLeq("output/Galerkin/eigvalLeq");
    //out_eigvalLeq << eigvalLeq << endl;
    displayCplx(eigvalLeq, out_eigvalLeq);*/

    double varOfH1 = -2 * dot(gettGiHj(0, 1), LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, Lrep_ * LinvH1);
    //double varCoeff = -.5 * dot( Lrep_ * LinvH1, LeqInv * Lrep_ * LinvH1);
    double varCoeff = -2 * dot( Lrep_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }
  
  ///
  ///Returns 
  DVec LangevinGalerkin::CVcoeffsVec() const
  {
    /*cout << "p    : " << gettGiHj(0,1) << endl << endl;
    cout << "Lp   : " << getLinvtGiHj(0, 1) << endl << endl;
    cout << "Linvp: " << getLtGiHj(0, 1) << endl << endl;*/
    
    //return -getLinvtGiHj(0, 1);
    //return gettGiHj(0,1);     // returns "p"
    
    
    if (doNonequilibrium())
    {
      cout << "Computing the control variate by a nonequilibrium LLT" << endl;
      SMat Keta = Leta_.transpose() * Leta_;
      DVec CV = -solve(Keta, Leta_.transpose() * gettGiHj(0,1));
      return CV;
    }
    else
    {
      cout << "Computing the control variate by an equilibrium LLT" << endl;
      //SMat Keq = Leq_.transpose() * Leq_;
      //return -solve(Keq, Leq_.transpose() * gettGiHj(0,1));
      
      Eigen::BiCGSTAB<SMat> solver(Leq_);
      return -solver.solve(-projectionOrthoG(gettGiHj(0,1)));
    }
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
    computeExpToTrigTens();
    
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

  void BoundaryLangevinGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = kron(trigToExpMat_, SIdP_);
    trigToExpTens_ = tensorPower(trigToExpTens_, nbOfParticles_);
    expToTrigTens_ = kron(expToTrigMat_, SIdP_);
    expToTrigTens_ = tensorPower(expToTrigTens_, nbOfParticles_);
    
    display(trigToExpTens_, "output/Galerkin/trigToExpTens");
  }


  void BoundaryLangevinGalerkin::compute()
  {
    cout << "start Galerkin::compute()" << endl;
    cout << "############ Leq ############" << endl;
    //cout << Leq_ << endl << endl;
    
    cout << "Leq_ size : " << Leq_.rows() << " x " << Leq_.cols() << endl;
    cout << "Leq_ nnz : " << Leq_.nonZeros() << endl;

    DMat DLeq(Leq_);
    display(Leq_, "output/Galerkin/Leq");

    //DMat DLeqSad = shapeSaddle(DLeq);
    //display(DLeqSad, "output/Galerkin/LeqSad");

    cout << "############ DLeqSad ############" << endl;
    //cout << DLeqSad << endl << endl;

    //DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);

    cout << "Computing DLeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DMat LeqInvSad = inv(DLeqSad);
    DMat LeqInv = invWithSaddle(Leq_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");

    DVec N0H2Trig = DVec::Zero(sizeOfBasis_);
    N0H2Trig(iTens(0, 2, 0)) = 1;
    display(N0H2Trig, "output/Galerkin/N0H2Trig");
    DVec N0H2 = trigToExpTens_ * N0H2Trig;
    display(N0H2, "output/Galerkin/N0H2");
    
    //DMat N0H2TrigMat(N0H2Trig, nbOfFourier_, nbOfHermite_); 
    //DMat N0H2Mat = trigToExpMat_ * N0H2TrigMat;
    
    DMat N0H2Mat = reshape(N0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));
    
    display(N0H2Mat, "output/Galerkin/N0H2Mat");
    
    DVec LinvN0H2 = solveWithSaddle(Leq_, N0H2);
    DMat LinvN0H2Mat = reshape(LinvN0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));
    display(LinvN0H2Mat, "output/Galerkin/LinvN0H2Mat");
    
    /*cx_vec eigvalLeq;
    cx_mat eigvecLeq;
    eig_gen(eigvalLeq, eigvecLeq, -DLeq);
    //eigs_gen(eigvalLeq, eigvecLeq, Leq_, 3);

    ofstream out_eigvalLeq("output/Galerkin/eigvalLeq");
    //out_eigvalLeq << eigvalLeq << endl;
    displayCplx(eigvalLeq, out_eigvalLeq);*/

    /*double varOfH1 = -2 * dot(gettGiHj(0,1), LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;*/

    /*DVec LinvL1LinvH1 = solveWithSaddle(Leq_, Lrep_ * LinvH1);
    //double varCoeff = -.5 * dot( Lrep_ * LinvH1, LeqInv * Lrep_ * LinvH1);
    double varCoeff = -2 * dot( Lrep_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;*/
  }



}



