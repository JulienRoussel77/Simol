#include "simol/statphys/controlVariate/Galerkin.hpp"

using std::cout;
using std::endl;
using std::ostream;

namespace simol
{

  Galerkin* createLangevinGalerkin(Input const& input)
  {
    if (input.doGalerkinCV())
      return new LangevinGalerkin(input);
    else
      return nullptr;
  }

  Galerkin* createBoundaryLangevinGalerkin(Input const& input)
  {
    if (input.doGalerkinCV())
      return new BoundaryLangevinGalerkin(input);
    else
      return nullptr;
  }

  DenseMatrix<double> Galerkin::shapeSaddle(const DenseMatrix<double>& A) const
  {
    DenseMatrix<double> Asad = DenseMatrix<double>::Zero(A.number_of_rows() + 1, A.number_of_columns() + 1);
    Asad.block(0, 0, A.number_of_rows(), A.number_of_columns()) = A.block(0, 0, A.number_of_rows(), A.number_of_columns());
    for (int iOfFourier2 = 0; iOfFourier2 <= 2 * maxOfFourier_; iOfFourier2++)
    {
      Asad(A.number_of_rows(), iTens(iOfFourier2, 0)) = expFourierCoeffs(iOfFourier2);
      Asad(iTens(iOfFourier2, 0), A.number_of_columns()) = expFourierCoeffs(iOfFourier2);
    }
    return Asad;
  }

  DenseMatrix<double> Galerkin::unshapeSaddle(const DenseMatrix<double>& Asad) const
  { return Asad.block(0, 0, Asad.number_of_rows() - 1, Asad.number_of_columns() - 1); }

  DVec Galerkin::shapeSaddle(const DVec& X) const
  {
    DVec Xsad = Vector<double>::Zero(X.size() + 1);
    Xsad.subvector(0, X.size()) = X;
    return Xsad;
  }

  DVec Galerkin::unshapeSaddle(const DVec& Xsad) const
  {
    return Xsad.subvector(0, Xsad.size() - 2);
  }

  DVec Galerkin::solveWithSaddle(const SMat& A, const DVec& X) const
  {
    DenseMatrix<double> DA = A;
    return solveWithSaddle(DA, X);
  }

  DVec Galerkin::solveWithSaddle(const DenseMatrix<double>& A, const DVec& X) const
  {
    cout << "solveWithSaddle...";
    DVec Xsad = shapeSaddle(X);
    DenseMatrix<double> Asad = shapeSaddle(A);
    DVec Bsad = Asad.solve(Xsad);
    cout << "OK ! (lambda = " << Bsad(Bsad.size() - 1) << ")" << endl;
    return unshapeSaddle(Bsad);
  }

  DenseMatrix<double> Galerkin::invWithSaddle(const SparseMatrix<double>& A) const
  {
    DenseMatrix<double> DA(A);
    return invWithSaddle(DA);
  }

  DenseMatrix<double> Galerkin::invWithSaddle(const DenseMatrix<double>& A) const
  {
    cout << "invWithSaddle...";
    DenseMatrix<double> Asad = shapeSaddle(A);
    DenseMatrix<double> Bsad = Asad.inverse();
    return unshapeSaddle(Bsad);
  }

  // C = ( A B_11    A B_12   ...
  //       A B_21    A B_22   ...
  //                              )
  SMat kron(const SMat& A, const SMat& B)
  {
    cout << "kron(const SMat& A, const SMat& B)" << endl;

    SMat C(A.number_of_rows()*B.number_of_rows(), A.number_of_columns()*B.number_of_columns());
    cout << "A : " << A.number_of_rows() << "x" << A.number_of_columns() << endl;
    cout << A(0, 0) << endl;
    cout << "B : " << B.number_of_rows() << "x" << B.number_of_columns() << endl;




    for (int jOfA = 0; jOfA < (int)A.number_of_columns(); ++jOfA)
    {
      for (SMat::iterator it(A, jOfA); it; ++it)
      {
        int iOfA = it.row();
        double valOfA = it.value();
        cout << "truc : " << it.row() << " " << it.col() << " " << it.value() << endl;
        for (int jOfB = 0; jOfB < (int)B.number_of_columns(); jOfB++)
        {
          for (SMat::iterator it2(B, jOfB); it2; ++it2)
          {
            int iOfB = it2.row();
            double valOfB = it2.value();
            cout << iOfA << "+" << A.number_of_rows() << "*" << iOfB << endl;
            cout << iOfA + A.number_of_rows() * iOfB << " , " << jOfA + A.number_of_columns() * jOfB << " ->" << valOfA*valOfB << endl;
            C(iOfA + A.number_of_rows() * iOfB, jOfA + A.number_of_columns() * jOfB) = valOfA * valOfB;
          }
        }
      }
    }
    cout << "end  kron(const SMat& A, const SMat& B)" << endl;
    return C;
  }

  DenseMatrix<double> kron(const DenseMatrix<double>& A, const DenseMatrix<double>& B)
  {
    cout << "kron(DenseMatrix<double>& A, DenseMatrix<double>& B)" << endl;
    DenseMatrix<double> C(A.number_of_rows()*B.number_of_rows(), A.number_of_columns()*B.number_of_columns());
    for (int iOfA = 0; iOfA < (int) A.number_of_rows(); iOfA++)
      for (int jOfA = 0; jOfA < (int) A.number_of_columns(); jOfA++)
        for (int iOfB = 0; iOfB < (int) B.number_of_rows(); iOfB++)
          for (int jOfB = 0; jOfB < (int) B.number_of_columns(); jOfB++)
            C(iOfA + A.number_of_rows() * iOfB, jOfA + A.number_of_columns() * jOfB) = A(iOfA, jOfA) * B(iOfB, jOfB);
    return C;
  }

  void displayCplx(const Vector<cplx>& X, ostream& out)
  {
    for (int i = 0; i < (int) X.size(); i++)
      out << real(X(i)) << " " << imag(X(i)) << endl;
  }

  void display(const Vector<double>& A, string path)
  {
    ofstream out(path);
    out << A;
  }

  void display(const DenseMatrix<double>& A, ostream& out)
  {
    for (int i = 0; i < (int) A.number_of_rows(); i++)
    {
      for (int j = 0; j < (int) A.number_of_columns(); j++)
      {
        if (true)//fabs(A(i,j)) > 1e-15)
        {
          out << A(i, j) << " ";
        }
        else
          out << "nan ";
      }
      out << endl;
    }
  }

  void display(const DenseMatrix<double>& A, string path)
  {
    ofstream out(path);
    display(A, out);
  }

  void display(const SMat& A, ostream& out)
  {
    DenseMatrix<double> DA = A;
    display(DA, out);
  }

  void display(const SMat& A, string path)
  {
    DenseMatrix<double> DA = A;
    display(DA, path);
  }

  //We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  //where C = ExpFourierBasis::basisCoefficient_
  void Galerkin::computeExpToTrigMat()
  {
    for (int iOfFourier2 = 1; iOfFourier2 <=  2 * (int)maxOfFourier_; iOfFourier2++)
    {
      trigToExpMat_(0, iOfFourier2) = expFourierCoeffs(iOfFourier2);
      trigToExpMat_(iOfFourier2, 0) = expFourierCoeffs(iOfFourier2);
    }
    trigToExpMat_(0, 0) = expFourierCoeffs(0) / sqrt(2.);

    for (int iOfFourier = 1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
      for (int jOfFourier = 1; jOfFourier <= (int) maxOfFourier_; jOfFourier++)
      {
        //cosine times cosine
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier) =
          expFourierCoeffs(2 * (iOfFourier + jOfFourier)) / sqrt(2.)
          + expFourierCoeffs(2 * abs(iOfFourier - jOfFourier)) / sqrt(2.);

        //sinus times sinus
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier - 1) =
          - expFourierCoeffs(2 * (iOfFourier + jOfFourier)) / sqrt(2.)
          + expFourierCoeffs(2 * abs(iOfFourier - jOfFourier)) / sqrt(2.);

        int eps = (iOfFourier >= jOfFourier) - (iOfFourier <= jOfFourier); // -1 if i smaller, 0 if equal, 1 if i larger
        //cosine times sinus
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier - 1) =
          expFourierCoeffs(2 * (iOfFourier + jOfFourier) - 1) / sqrt(2.)
          + eps * expFourierCoeffs(2 * abs(iOfFourier - jOfFourier) - 1) / sqrt(2.);

        //sinus times cosine
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier) =
          expFourierCoeffs(2 * (iOfFourier + jOfFourier) - 1) / sqrt(2.)
          - eps * expFourierCoeffs(2 * abs(iOfFourier - jOfFourier) - 1) / sqrt(2.);
      }


    ofstream out_trigToExpMat("../output/Galerkin/trigToExpMat");
    display(trigToExpMat_, out_trigToExpMat);

    expToTrigMat_  = trigToExpMat_.inverse();
  }

  DenseMatrix<double> Galerkin::convertToTrigBasis(const DenseMatrix<double>& X)
  {
    return expToTrigMat_ * X;
  }

  void Galerkin::createQ()
  {
    Q_(1, 0) = amplitude_ * beta_ / (2 * sqrt(2.));

    for (int iOfFourier = 1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
    {
      Q_(2 * iOfFourier - 2, 2 * iOfFourier - 1) = amplitude_ * beta_ / 4;   //overwritten if i=1
      Q_(2 * iOfFourier    , 2 * iOfFourier - 1) = iOfFourier;
      if (iOfFourier != (int) maxOfFourier_)
        Q_(2 * iOfFourier + 2, 2 * iOfFourier - 1) = - amplitude_ * beta_ / 4;
    }
    Q_(0, 1) = amplitude_ * beta_ / (2 * sqrt(2.));

    for (int iOfFourier = 1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
    {
      if (iOfFourier != 1)
        Q_(2 * iOfFourier - 3, 2 * iOfFourier) = - amplitude_ * beta_ / 4;
      Q_(2 * iOfFourier - 1, 2 * iOfFourier) = - iOfFourier;
      if (iOfFourier != (int) maxOfFourier_)
        Q_(2 * iOfFourier + 1, 2 * iOfFourier) =   amplitude_ * beta_ / 4;
    }

    tQ_ = Q_.adjoint();
  }

  void Galerkin::createP()
  {
    for (int iOfHermite = 1; iOfHermite < nbOfHermite_; iOfHermite++)
      P_(iOfHermite - 1, iOfHermite) = sqrt(beta_ * iOfHermite);

    tP_ = P_.adjoint();
  }

  void Galerkin::createLthm0()
  {
    cout << "createLthm0" << endl;
    for (int iOfHermite = 1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
    cout << "end createLthm0" << endl;
  }

  using namespace Eigen;

  Galerkin::Galerkin(Input const& input):       //ex : [0:4]
    nbOfParticles_(input.nbOfParticles()),
    nbOfFourier_(input.nbOfFourier()),  //ex : 5
    nbOfHermite_(input.nbOfHermite()),
    maxOfFourier_((nbOfFourier_ - 1) / 2),  //ex : 2
    sizeOfBasis_(pow(nbOfFourier_ * nbOfHermite_, nbOfParticles_)),
    SIdQ_(speye<double>(nbOfFourier_, nbOfFourier_)),
    SIdP_(speye<double>(nbOfHermite_, nbOfHermite_)),
    DIdQ_(DenseMatrix<double>::Identity(nbOfFourier_)),
    DIdP_(DenseMatrix<double>::Identity(nbOfHermite_)),
    Q_(nbOfFourier_, nbOfFourier_),
    P_(nbOfHermite_, nbOfHermite_),
    tQ_(nbOfFourier_, nbOfFourier_),
    tP_(nbOfHermite_, nbOfHermite_),
    Lthm0_(nbOfHermite_, nbOfHermite_),
    Lthm_(sizeOfBasis_, sizeOfBasis_),
    Lham_(sizeOfBasis_, sizeOfBasis_),
    L1_(sizeOfBasis_, sizeOfBasis_),
    Leq_(sizeOfBasis_, sizeOfBasis_),
    Leta_(sizeOfBasis_, sizeOfBasis_),
    beta_(input.beta()),
    gamma_(input.gamma()),
    amplitude_(input.amplitude()),
    externalForce_(input.externalForce()),
    nbOfIntegrationNodes_(1000),
    trigToExpMat_(DenseMatrix<double, eigen>::Zero(nbOfFourier_, nbOfFourier_)),
    expToTrigMat_(DenseMatrix<double, eigen>::Zero(nbOfFourier_, nbOfFourier_)),
    trigToExpTens_(sizeOfBasis_, sizeOfBasis_),
    expToTrigTens_(sizeOfBasis_, sizeOfBasis_),
    potential_(createPotential(input)),
    basis_(input, *potential_)
  {

    Eigen::SparseMatrix<double> A(2, 3);
    A.insert(0, 1) = 34;
    A.insert(1, 2) = 56;
    SMat::iterator it(A, 1);
    cout << it.row() << "\t";
    cout << it.col() << "\t";
    cout << it.value() << endl;
    cout << "---------------------" << endl;

    SMat B(2, 3);
    B.insert(0, 1) = 24;
    B.insert(1, 2) = 36;
    SMat::iterator it2(B, 1);
    cout << it2.row() << "\t";
    cout << it2.col() << "\t";
    cout << it2.value() << endl;

    SMat C = speye<double>(3, 3);
    for (int jOfC = 0; jOfC < (int) C.number_of_columns(); ++jOfC)
      for (int iOfC = 0; iOfC < (int) C.number_of_rows(); ++iOfC)
        cout << C(iOfC, jOfC) << endl;



    SMat::iterator itTest(C, 0);
    cout << "bool : " << (bool)itTest << " (" << itTest.row() << " , " << itTest.col() << ") -> " << itTest.value() << endl;
    ++itTest;
    cout << "bool : " << (bool)itTest << " (" << itTest.row() << " , " << itTest.col() << ") -> " << itTest.value() << endl;
    ++itTest;
    cout << "bool : " << (bool)itTest << " (" << itTest.row() << " , " << itTest.col() << ") -> " << itTest.value() << endl;
    ++itTest;
    cout << "bool : " << (bool)itTest << " (" << itTest.row() << " , " << itTest.col() << ") -> " << itTest.value() << endl;

    for (SMat::iterator it(C, 0); it; ++it)
      cout << "test : (" << it.row() << " , " << it.col() << ") -> " << it.value() << endl;

    for (SMat::iterator it(C, 1); it; ++it)
      cout << "test2 : (" << it.row() << " , " << it.col() << ") -> " << it.value() << endl;


    assert(nbOfFourier_ % 2 == 1);
    cout << endl << "Number of modes : " << nbOfFourier_ << " x " << nbOfHermite_ << endl;

    //computeFourierCoeffsExp();
    // Computation of the passage matrix
    computeExpToTrigMat();

    cout << "Computing Q...";
    createQ();
    display(Q_, "../output/Galerkin/Q");
    cout << "OK" << endl;

    cout << "Computing Q...";
    createP();
    display(P_, "../output/Galerkin/P");
    cout << "OK" << endl;


  }

  int Galerkin::nbOfVariables() const
  {
    return 2 * nbOfParticles_;
  }

  int Galerkin::nbOfParticles() const
  {
    return nbOfParticles_;
  }

  const double& Galerkin::expFourierCoeffs(int iOfElt) const
  {
    return basis_.expFourierCoeffs(iOfElt);
  }

  //psi = (1,1  1,2  ...  1,N_H  2,1 ... )
  //N_H blocks of size N_G (we concatene the columns of the matrix)
  int Galerkin::iTens(int iOfFourier2, int iOfHermite) const
  {
    assert(iOfFourier2 < nbOfFourier_ && iOfHermite < nbOfHermite_);
    return nbOfFourier_ * iOfHermite + iOfFourier2;
  }


  void Galerkin::compute()
  {
    cout << "start Galerkin::compute()" << endl;

    cout << "############ Leq ############" << endl;

    DenseMatrix<double> DLeq(Leq_);
    display(Leq_, "../output/Galerkin/Leq");

    cout << "############ DLeqSad ############" << endl;

    cout << "Computing DLeqInv...";
    DenseMatrix<double> DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, "../output/Galerkin/DLeqInv");


    DVec H1Trig = Vector<double>::Zero(sizeOfBasis_);
    H1Trig(iTens(0, 1)) = 1;

    DenseMatrix<double> H1TrigMat(H1Trig, nbOfFourier_, nbOfHermite_);

    DVec H1 = trigToExpTens_ * H1Trig;

    cout << "############ H1Mat ############" << endl;
    display(gettGiHj(0, 1), "../output/Galerkin/H1");
    DenseMatrix<double> H1Mat(gettGiHj(0, 1), nbOfFourier_, nbOfHermite_);
    display(H1Mat, "../output/Galerkin/H1Mat");

    display(H1Trig, "../output/Galerkin/H1Trig");
    display(H1TrigMat, "../output/Galerkin/H1TrigMat");


    display(gettGiHj(0, 0), "../output/Galerkin/G0Trig");

    DenseMatrix<double> G0Mat(gettGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(G0Mat, "../output/Galerkin/G0Mat");
    DenseMatrix<double> LG0Mat(getLtGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(LG0Mat, "../output/Galerkin/LG0Mat");

    display(gettGiHj(1, 2), "../output/Galerkin/G1H2Trig");
    display(getLtGiHj(1, 2), "../output/Galerkin/LG1H2Trig");

    display(getLtGiHj(0, 1), "../output/Galerkin/LH1");
    DVec LinvH1 = getLinvtGiHj(0, 1);
    display(LinvH1, "../output/Galerkin/LinvH1");

    DenseMatrix<double> LinvH1Mat(LinvH1, nbOfFourier_, nbOfHermite_);

    cout << "############ LinvH1Mat ############" << endl;
    display(LinvH1Mat, "../output/Galerkin/LinvH1Mat");

    DenseMatrix<double> LinvH1MatTrig = convertToTrigBasis(LinvH1Mat);
    display(LinvH1MatTrig, "../output/Galerkin/LinvH1MatTrig");

    DVec H1back = Leq_ * LinvH1;
    DenseMatrix<double> H1backMat(H1back, nbOfFourier_, nbOfHermite_);
    display(H1backMat, "../output/Galerkin/H1backMat");

    double varOfH1 = -2 * inner_product(gettGiHj(0, 1), LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, L1_ * LinvH1);
    double varCoeff = -2 * inner_product( L1_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }

  DVec Galerkin::gettGiHj(int i, int j) const
  {
    DVec GiHjTrig = Vector<double>::Zero(sizeOfBasis_);
    GiHjTrig(iTens(i, j)) = 1;
    DVec GiHj = trigToExpTens_ * GiHjTrig;
    return GiHj;
  }

  DVec Galerkin::gettGiHjTrig(int i, int j) const
  {
    DVec GiHjTrig = Vector<double>::Zero(sizeOfBasis_);
    GiHjTrig(iTens(i, j)) = 1;
    return GiHjTrig;
  }

  DVec Galerkin::getLtGiHj(int i, int j) const
  { return Leq_ * gettGiHj(i, j); }

  DVec Galerkin::getLtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * getLtGiHj(i, j); }

  DVec Galerkin::getLinvtGiHj(int i, int j) const
  { return solveWithSaddle(Leq_, gettGiHj(i, j)); }

  DVec Galerkin::getLinvtGiHjTrig(int i, int j) const
  { return expToTrigTens_ * solveWithSaddle(Leq_, gettGiHj(i, j)); }

  SMat Galerkin::CVcoeffs() const
  {
    DVec vecCoeffs = getLinvtGiHj(0, 1);

    return SparseMatrix<double>(vecCoeffs, vecCoeffs.size(), 1);
  }



  //#### LangevinGalerkin ####


  LangevinGalerkin::LangevinGalerkin(Input const& input):
    Galerkin(input)
  {
    cout << "############ Lham ############" << endl;
    Lham_ = kron(Q_, tP_) - kron(tQ_, P_);
    display(Lham_, "../output/Galerkin/Lham");

    cout << "############ Lthm ############" << endl;
    createLthm();
    display(Lthm_, "../output/Galerkin/Lthm");

    cout << "############ Leq ############" << endl;
    Leq_  = Lham_ + gamma_ * Lthm_;

    cout << "############ L1 ############" << endl;
    //L1_ = kron(SMat::Identity(nbOfFourier_) , P_);
    L1_ = kron(SIdQ_, P_);

    cout << "############ Leta ############" << endl;
    Leta_ = Leq_ + externalForce_ * L1_;
  }

  void LangevinGalerkin::createLthm()
  {
    createLthm0();
    Lthm_ = kron(SIdQ_, Lthm0_);
  }

  void LangevinGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = simol::kron(trigToExpMat_, DIdP_);
    expToTrigTens_ = simol::kron(expToTrigMat_, DIdP_);
  }


  //#### BoundaryLangevinGalerkin ####

  SMat tensorPower(SMat const& A, int power)
  {
    SMat temp = A;
    for (int i = 1; i < power; i++)
      temp = kron(temp, A);
    return temp;
  }

  DenseMatrix<double> tensorPower(DenseMatrix<double> const& A, int power)
  {
    DenseMatrix<double> temp = A;
    for (int i = 1; i < power; i++)
      temp = kron(temp, A);
    return temp;
  }

  SMat BoundaryLangevinGalerkin::PMatToTens(SMat const& PMat, int iOfParticleP)
  {
    assert(iOfParticleP < nbOfParticles_);
    SMat res = speye<double>(1, 1);
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
    SMat res = speye<double>(1, 1);
    for (int i = 0; i < nbOfParticles_; i++)
    {
      cout << "a" << endl;
      if (i == iOfParticleQ) res = kron(res, QMat);
      else res = kron(res, SIdQ_);
      cout << "b" << endl;
      if (i == iOfParticleP) res = kron(res, PMat);
      else res = kron(res, SIdP_);
    }
    cout << "res : " << res.number_of_rows() << " x " << res.number_of_columns() << endl;
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
    createLham();
    createLthm();

    Leq_ = Lham_ + gamma_ * Lthm_;
  }

  void BoundaryLangevinGalerkin::createLham()
  {
    cout << "BoundaryLangevinGalerkin::createLham()" << endl;
    for (int i = 0; i < nbOfParticles_; i++)
    {
      Lham_ += doubleMatToTens(Q_, tP_, i, i);
      if (i > 0) Lham_ -= doubleMatToTens(Q_, tP_, i, i - 1);
      Lham_ -= doubleMatToTens(tQ_, P_, i, i);
      if (i < nbOfParticles_ - 1) Lham_ += doubleMatToTens(Q_, tP_, i + 1, i);
    }
    Lham_ /= beta_;
    display(Lham_, "../output/Galerkin/Lham");
    cout << "end BoundaryLangevinGalerkin::createLham()" << endl;
  }

  void BoundaryLangevinGalerkin::createLthm()
  {
    createLthm0();
    for (int i = 0; i < nbOfParticles_; i++)
    {
      Lthm_ += PMatToTens(Lthm0_, i);
    }
    Lthm_ /= beta_;
    display(Lthm_, "../output/Galerkin/Lthm");
  }

  void BoundaryLangevinGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = simol::kron(trigToExpMat_, DIdP_);
    trigToExpTens_ = tensorPower(trigToExpTens_, nbOfParticles_);
    expToTrigTens_ = simol::kron(expToTrigMat_, DIdP_);
    expToTrigTens_ = tensorPower(expToTrigTens_, nbOfParticles_);
  }


  void BoundaryLangevinGalerkin::compute()
  {
    cout << "start Galerkin::compute()" << endl;
    cout << "############ Leq ############" << endl;

    DenseMatrix<double> DLeq(Leq_);
    display(Leq_, "../output/Galerkin/Leq");

    cout << "Computing DLeqInv...";
    DenseMatrix<double> DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;
    display(DLeqInv, "../output/Galerkin/DLeqInv");

    DVec N0H2Trig = Vector<double>::Zero(sizeOfBasis_);
    N0H2Trig(iTens(0, 2, 0)) = 1;

    DVec N0H2 = trigToExpTens_ * N0H2Trig;

    DenseMatrix<double> N0H2Mat(N0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));

  }



}



