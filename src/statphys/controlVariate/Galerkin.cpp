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
  
  //
  //Create a unitary matrix rectangular matrix representing a bon of the orthogonal of u
  /*DMat createOrthoBasis(Vector<double> const& u)
  {
    Vector<double> u0 = u / u.norm();
    DMat A = DMat::Random(u0.size(), u0.size() - 1);
    for (int jOfCol = 0; jOfCol < (int) A.numberOfColumns(); jOfCol++)
    {
      A.column(jOfCol) -= dot(A.column(jOfCol), u0) * u0;
      for (int jOfCol2 = 0; jOfCol2 < jOfCol; jOfCol2++)
        A.column(jOfCol) -= dot(A.column(jOfCol), A.column(jOfCol2)) * A.column(jOfCol2);
      assert(A.column(jOfCol).norm() != 0);           // can be zero if the initial random family is not free
      A.column(jOfCol) /= A.column(jOfCol).norm();
    }
    return A;
  }*/
  
  

  DenseMatrix<double> Galerkin::shapeSaddle(const DenseMatrix<double>& A) const
  {
    DenseMatrix<double> Asad = DenseMatrix<double>::Zero(A.numberOfRows() + 1, A.numberOfColumns() + 1);
    Asad.block(0, 0, A.numberOfRows(), A.numberOfColumns()) = A.block(0, 0, A.numberOfRows(), A.numberOfColumns());
    for (int iOfFourier2 = 0; iOfFourier2 <= 2 * maxOfFourier_; iOfFourier2++)
    {
      Asad(A.numberOfRows(), iTens(iOfFourier2, 0)) = expFourierMeans(iOfFourier2);
      Asad(iTens(iOfFourier2, 0), A.numberOfColumns()) = expFourierMeans(iOfFourier2);
    }
    return Asad;
  }
  
  /*DenseMatrix<double> Galerkin::shapePrec(const DenseMatrix<double>& A) const
  {
    return A + norm2gVector() / (1 - norm2gVector()) * extProduct(gVector(), A.adjoint() * gVector());
  }*/

  DenseMatrix<double> Galerkin::unshapeSaddle(const DenseMatrix<double>& Asad) const
  { return Asad.block(0, 0, Asad.numberOfRows() - 1, Asad.numberOfColumns() - 1); }

  DVec Galerkin::shapeSaddle(const DVec& X) const
  {
    DVec Xsad = Vector<double>::Zero(X.size() + 1);
    Xsad.subvec(0, X.size()) = X;
    return Xsad;
  }

  DVec Galerkin::unshapeSaddle(const DVec& Xsad) const
  {
    return Xsad.subvec(0, Xsad.size() - 1);
  }

  // TODO: utiliser des solveurs lineaires plutôt qu'inverser
  // de plus, on n'inverse jamais une matrice creuse car
  // son inverse est generalement dense
  // par consequent, pas de fonction membre inverse() dans SparseMatrix

  /*DenseMatrix<double> inverse(const SMat& A)
  {
    DenseMatrix<double> Id = DenseMatrix<double>::Identity(A.numberOfRows());
    DenseMatrix<double> C = spsolve(A, Id);
    return C;
  }

  DenseMatrix<double> inverse(const DenseMatrix<double>& A)
  {
    DenseMatrix<double> Id = DenseMatrix<double>::Identity(A.numberOfRows());
    DenseMatrix<double> C = solve(A, Id);
    return C;
  }*/
  
  //Denote P the projector on the orthogonal of u
  //Compute a vector x such that (1-P) A x = (1-P) b and dot(u, x) = 0 
  /*DVec Galerkin::pseudoSolve(SMat const& A, DVec const& b, DVec const& u) const
  {
    DMat B = computeOrthoBasis(u);
    assert(false);
    return b;
  }*/
  
  DVec Galerkin::solve(SMat const& A, DVec const& Y) const
  {
    //cout << "wiiiiiiiiiiiiiiiiiiii" << endl;
    //DenseMatrix<double> DA = A;
    //return solveWithSaddle(DA, Y);
    
    Eigen::SparseMatrix<double> eigA = A.wrapped_;
    Eigen::VectorXd eigY = Y.wrapped_;
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > solver;
    solver.compute(eigA);
    Eigen::VectorXd eigX = solver.solve(eigY);
    cout << "Solver info : " << solver.info() << endl;
    assert(solver.info() == Eigen::Success);
    return DVec(eigX);
  }
  
  DVec Galerkin::solve(const DenseMatrix<double>& A, const DVec& X) const
  {    
    return A.solve(X);
  }

  DVec Galerkin::solveWithSaddle(SMat const& A, DVec const& Y) const
  {
    
    Eigen::SparseMatrix<double> eigA = A.wrapped_;
    Eigen::VectorXd eigY = Y.wrapped_;
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > solver;
    solver.compute(eigA);
    Eigen::VectorXd eigX = solver.solve(eigY);
    cout << "Solver info : " << solver.info() << endl;
    assert(solver.info() == Eigen::Success);
    return DVec(eigX);
  }

  DVec Galerkin::solveWithSaddle(const DenseMatrix<double>& A, const DVec& X) const
  {
    cout << "solveWithSaddle..."; cout.flush();
    cout << "X -> " << X.size() << endl;
    DVec Xsad = shapeSaddle(X);
    cout << "Xsad -> " << Xsad.size() << endl;
    cout << "A -> " << A.numberOfRows() << "x" << A.numberOfColumns() << endl;
    DenseMatrix<double> Asad = shapeSaddle(A);
    cout << "Asad -> " << Asad.numberOfRows() << "x" << Asad.numberOfColumns() << endl;
    DVec Bsad = Asad.solve(Xsad);
    cout << "Bsad -> " << Bsad.size() << endl;
    cout << "OK ! (lambda = " << Bsad(Bsad.size() - 1) << ")" << endl;
    return unshapeSaddle(Bsad);
  }

  DenseMatrix<double> Galerkin::invWithSaddle(const SparseMatrix<double>& A) const
  {
    DenseMatrix<double> DA = A.dense();
    return invWithSaddle(DA);
  }

  DenseMatrix<double> Galerkin::invWithSaddle(const DenseMatrix<double>& A) const
  {
    cout << "invWithSaddle..."; cout.flush();
    DenseMatrix<double> Asad = shapeSaddle(A);
    DenseMatrix<double> Bsad = Asad.inverse();
    cout << "fin invWithSaddle..."; cout.flush();
    return unshapeSaddle(Bsad);
  }

  // C = ( A B_11    A B_12   ...
  //       A B_21    A B_22   ...
  //                              )
  SMat kron(const SMat& A, const SMat& B)
  {
    /*double a_1 = A.numberOfColumns() * (B.numberOfColumns() - 1)/(A.numberOfColumns() - 1);
    double b_1 = a_1 - B.numberOfColumns();
    double a_2 = A.numberOfRows() * (B.numberOfRows() - 1)/(A.numberOfRows() - 1);
    double b_2 = a_2 - B.numberOfRows();*/

    /*cout << "A : O <= i < " << A.numberOfRows() << ", O <= jOfA < " << A.numberOfColumns() << endl;
    cout << "B : O <= i2 < " << B.numberOfRows() << ", O <= jOfB < " << B.numberOfColumns() << endl;
    cout << "We keep the coefficients such that j * (jOfB + " << b_1 << ") <= " << a_1
      << " and such that i * (i2 + " << b_2 << ") <= " << a_2 << endl; */
    SMat C(A.numberOfRows()*B.numberOfRows(), A.numberOfColumns()*B.numberOfColumns());

    //SMat C(A.size() % B.size());          //element-wise product of the dimensions
    //

    //C(0,0) = 1;

    for (int jOfA = 0; jOfA < (int)A.numberOfColumns(); ++jOfA)
    {
      for (SMat::iterator it(A, jOfA); it; ++it)
      {
        int iOfA = it.row();
        double valOfA = it.value();
        for (int jOfB = 0; jOfB < (int)B.numberOfColumns(); jOfB++)
        {
          for (SMat::iterator it2(B, jOfB); it2; ++it2)
          {
            int iOfB = it2.row();
            double valOfB = it2.value();
            C(iOfA + A.numberOfRows() * iOfB, jOfA + A.numberOfColumns() * jOfB) = valOfA * valOfB;
          }
        }
      }
    }
    return C;
  }

  DenseMatrix<double> kron(const DenseMatrix<double>& A, const DenseMatrix<double>& B)
  {
    DenseMatrix<double> C(A.numberOfRows()*B.numberOfRows(), A.numberOfColumns()*B.numberOfColumns());
    for (int iOfA = 0; iOfA < (int) A.numberOfRows(); iOfA++)
      for (int jOfA = 0; jOfA < (int) A.numberOfColumns(); jOfA++)
        for (int iOfB = 0; iOfB < (int) B.numberOfRows(); iOfB++)
          for (int jOfB = 0; jOfB < (int) B.numberOfColumns(); jOfB++)
            C(iOfA + A.numberOfRows() * iOfB, jOfA + A.numberOfColumns() * jOfB) = A(iOfA, jOfA) * B(iOfB, jOfB);
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
    //display(A, out);
    out << A;
  }

  void display(const DenseMatrix<double>& A, ostream& out)
  {
    for (int i = 0; i < (int) A.numberOfRows(); i++)
    {
      for (int j = 0; j < (int) A.numberOfColumns(); j++)
      {
        if (true)//fabs(A(i,j)) > 1e-15)
        {
          //cout << i << " " << j << " " << A(i,j) << endl;
          out << setw(10) << A(i, j) << " ";
        }
        else
          out << setw(10) << "nan ";
      }
      out << endl;
    }
  }

  void display(const DenseMatrix<double>& A, string path)
  {
    ofstream out(path);
    assert(out.is_open());
    display(A, out);
  }

  void display(const SMat& A, ostream& out)
  {
    DenseMatrix<double> DA = A.dense();
    display(DA, out);
  }

  void display(const SMat& A, string path)
  {
    DenseMatrix<double> DA = A.dense();
    display(DA, path);
  }
  
  

  //We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
  //where C = ExpFourierBasis::basisCoefficient_
  void Galerkin::computeExpToTrigMat()
  {
    for (int iOfFourier2 = 1; iOfFourier2 <=  2 * (int)maxOfFourier_; iOfFourier2++)
    {
      trigToExpMat_(0, iOfFourier2) = expFourierMeans(iOfFourier2);
      trigToExpMat_(iOfFourier2, 0) = expFourierMeans(iOfFourier2);
    }
    trigToExpMat_(0, 0) = expFourierMeans(0) / sqrt(2.);

    for (int iOfFourier = 1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
      for (int jOfFourier = 1; jOfFourier <= (int) maxOfFourier_; jOfFourier++)
      {
        //cosine times cosine
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier) =
          expFourierMeans(2 * (iOfFourier + jOfFourier)) / sqrt(2.)
          + expFourierMeans(2 * abs(iOfFourier - jOfFourier)) / sqrt(2.);

        //sinus times sinus
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier - 1) =
          - expFourierMeans(2 * (iOfFourier + jOfFourier)) / sqrt(2.)
          + expFourierMeans(2 * abs(iOfFourier - jOfFourier)) / sqrt(2.);

        int eps = (iOfFourier >= jOfFourier) - (iOfFourier <= jOfFourier); // -1 if i smaller, 0 if equal, 1 if i larger
        //cosine times sinus
        trigToExpMat_(2 * iOfFourier, 2 * jOfFourier - 1) =
          expFourierMeans(2 * (iOfFourier + jOfFourier) - 1) / sqrt(2.)
          + eps * expFourierMeans(2 * abs(iOfFourier - jOfFourier) - 1) / sqrt(2.);

        //sinus times cosine
        trigToExpMat_(2 * iOfFourier - 1, 2 * jOfFourier) =
          expFourierMeans(2 * (iOfFourier + jOfFourier) - 1) / sqrt(2.)
          - eps * expFourierMeans(2 * abs(iOfFourier - jOfFourier) - 1) / sqrt(2.);
      }


    ofstream out_trigToExpMat("output/Galerkin/trigToExpMat");
    display(trigToExpMat_, out_trigToExpMat);

    expToTrigMat_  = trigToExpMat_.inverse();
  }

  DenseMatrix<double> Galerkin::convertToTrigBasis(const DenseMatrix<double>& X)
  {
    return expToTrigMat_ * X;
  }

  void Galerkin::createQ()
  {
    basis_(0)->gradMatrix(Q_);
    
    /*for (int iOfFourier2 = 0; iOfFourier2 < (int) nbOfFourier_; iOfFourier2++)
      for (int jOfFourier2 = 0; jOfFourier2 < (int) nbOfFourier_; jOfFourier2++)
        Q_(iOfFourier2, jOfFourier2) = basis_(0)->xGradY(iOfFourier2, jOfFourier2);*/
    
      
    /*Q_(1, 0) = amplitude_ * beta_ / (2 * sqrt(2.));

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
    }*/

    tQ_ = Q_.adjoint();
  }

  void Galerkin::createP()
  {
    for (int iOfHermite = 1; iOfHermite < nbOfHermite_; iOfHermite++)
      //P_(iOfHermite, iOfHermite+1) = cplx(sqrt(beta_) * sqrt(iOfHermite+1.), 0.);
      P_(iOfHermite - 1, iOfHermite) = sqrt(beta_ * iOfHermite);

    tP_ = P_.adjoint();
  }

  /*void Galerkin::createLthm0()
  {
    cout << "createLthm0" << endl;
    for (int iOfHermite = 1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
    cout << "end createLthm0" << endl;
  }*/

  using namespace Eigen;

  Galerkin::Galerkin(Input const& input):       //ex : [0:4]
    nbOfParticles_(input.nbOfParticles()),
    nbOfFourier_(input.nbOfFourier()),  //ex : 5
    nbOfHermite_(input.nbOfHermite()),
    maxOfFourier_((nbOfFourier_ - 1) / 2),  //ex : 2
    sizeOfBasis_(pow(nbOfFourier_ * nbOfHermite_, nbOfParticles_)),
    SIdQ_(speye<double>(nbOfFourier_, nbOfFourier_)),
    SIdP_(speye<double>(nbOfHermite_, nbOfHermite_)),
    DIdQ_(eye<double>(nbOfFourier_, nbOfFourier_)),
    DIdP_(eye<double>(nbOfHermite_, nbOfHermite_)),
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
    //expFourierMeans_(2 * nbOfFourier_, 0),
    trigToExpMat_(DenseMatrix<double, eigen>::Zero(nbOfFourier_, nbOfFourier_)),
    expToTrigMat_(DenseMatrix<double, eigen>::Zero(nbOfFourier_, nbOfFourier_)),
    trigToExpTens_(sizeOfBasis_, sizeOfBasis_),
    expToTrigTens_(sizeOfBasis_, sizeOfBasis_),
    potential_(createPotential(input)),
    basis_(input, *potential_)
  {
    assert(nbOfFourier_ % 2 == 1);
    cout << endl << "Number of modes : " << nbOfFourier_ << " x " << nbOfHermite_ << endl;

    //computeFourierCoeffsExp();
    // Computation of the passage matrix
    computeExpToTrigMat();
    
    cout << "Computing Q..."; cout.flush();
    createQ();
    display(Q_, "output/Galerkin/Q");
    cout << "OK" << endl;

    cout << "Computing P..."; cout.flush();
    createP();
    display(P_, "output/Galerkin/P");
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
  
  const double& Galerkin::expFourierMeans(int iOfElt) const
  {
    return basis_.expFourierMeans(iOfElt);
  }
  
  const Vector<double>& Galerkin::gVector() const
  {
    return basis_.gVector();
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


  void Galerkin::compute()
  {
    cout << "start Galerkin::compute()" << endl;

    cout << "############ Leq ############" << endl;
    //cout << Leq_ << endl << endl;

    DenseMatrix<double> DLeq = Leq_.dense();
    display(Leq_, "output/Galerkin/Leq");

    //DenseMatrix<double> DLeqSad = shapeSaddle(DLeq);
    //display(DLeqSad, "output/Galerkin/LeqSad");

    cout << "############ DLeqSad ############" << endl;
    //cout << DLeqSad << endl << endl;

    //DenseMatrix<double> IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);

    cout << "Computing DLeqInv..."; cout.flush();
    //DenseMatrix<double> LeqInv = invWithSaddle(DLeq);
    DenseMatrix<double> DLeqInv = invWithSaddle(DLeq);

    cout << "############ LeqInv ############" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");


    DVec H1Trig = Vector<double>::Zero(sizeOfBasis_);
    H1Trig(iTens(0, 1)) = 1;

    //DenseMatrix<double> H1TrigMat = reshape(H1Trig, nbOfFourier_, nbOfHermite_);
    DenseMatrix<double> H1TrigMat(H1Trig, nbOfFourier_, nbOfHermite_);

    //DenseMatrix<double> H1Mat = trigToExpMat_ * H1TrigMat;
    //DVec H1 = reshape(H1Mat, nbOfFourier_ * nbOfHermite_, 1);

    DVec H1 = trigToExpTens_ * H1Trig;


    /*DVec H1(sizeOfBasis_);
    for (int iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
      H1(iTens(iOfFourier2, 1)) = expFourierMeans_[iOfFourier2];*/

    cout << "############ H1Mat ############" << endl;
    display(gettGiHj(0, 1), "output/Galerkin/H1");
    //DenseMatrix<double> H1Mat = reshape(gettGiHj(0,1), nbOfFourier_, nbOfHermite_);
    DenseMatrix<double> H1Mat(gettGiHj(0, 1), nbOfFourier_, nbOfHermite_);
    display(H1Mat, "output/Galerkin/H1Mat");

    display(H1Trig, "output/Galerkin/H1Trig");
    display(H1TrigMat, "output/Galerkin/H1TrigMat");


    display(gettGiHj(0, 0), "output/Galerkin/G0Trig");

    DenseMatrix<double> G0Mat(gettGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(G0Mat, "output/Galerkin/G0Mat");
    DenseMatrix<double> LG0Mat(getLtGiHj(0, 0), nbOfFourier_, nbOfHermite_);
    display(LG0Mat, "output/Galerkin/LG0Mat");

    display(gettGiHj(1, 2), "output/Galerkin/G1H2Trig");
    display(getLtGiHj(1, 2), "output/Galerkin/LG1H2Trig");

    display(getLtGiHj(0, 1), "output/Galerkin/LH1");
    DVec LinvH1 = getLinvtGiHj(0, 1);
    display(LinvH1, "output/Galerkin/LinvH1");

    //double lambda = LinvH1Sad(sizeOfBasis_);

    //cout << "############ lambda ############" << endl;
    //cout << lambda << endl << endl;

    DenseMatrix<double> LinvH1Mat(LinvH1, nbOfFourier_, nbOfHermite_);

    cout << "############ LinvH1Mat ############" << endl;
    display(LinvH1Mat, "output/Galerkin/LinvH1Mat");

    DenseMatrix<double> LinvH1MatTrig = convertToTrigBasis(LinvH1Mat);
    display(LinvH1MatTrig, "output/Galerkin/LinvH1MatTrig");

    DVec H1back = Leq_ * LinvH1;
    DenseMatrix<double> H1backMat(H1back, nbOfFourier_, nbOfHermite_);
    display(H1backMat, "output/Galerkin/H1backMat");

    /*cx_vec eigvalLeq;
    cx_mat eigvecLeq;
    eig_gen(eigvalLeq, eigvecLeq, -DLeq);*/

    //DenseMatrix<double> DLeqW = Leq_.wrapped_;
    EigenSolver<MatrixXd> eigSol(DLeq.wrapped_);

    /*/eigs_gen(eigvalLeq, eigvecLeq, Leq_, 20);

    ofstream out_eigvalLeq("output/Galerkin/eigvalLeq");
    //out_eigvalLeq << eigvalLeq << endl;
    displayCplx(eigvalLeq, out_eigvalLeq);*/

    double varOfH1 = -2 * dot(gettGiHj(0, 1), LinvH1);
    cout << "varOfH1 = " << varOfH1 << endl;
    cout << "conductivity = " << .5 * varOfH1 << endl;

    DVec LinvL1LinvH1 = solveWithSaddle(Leq_, L1_ * LinvH1);
    //double varCoeff = -.5 * dot( L1_ * LinvH1, LeqInv * L1_ * LinvH1);
    double varCoeff = -2 * dot( L1_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;
  }
  
  EigenSolver<MatrixXd> Galerkin::getEigenSolver() const
  {
    cout << "Getting eigen elements by a dense matrix method !" << endl;
    DenseMatrix<double> DLeq = Leq_.dense();
    return EigenSolver<MatrixXd>(DLeq.wrapped_);
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
  
  void Galerkin::computeEigen() const
  {
    Eigen::EigenSolver<Eigen::MatrixXd> eigSol = getEigenSolver();
    //const cplx* eigVal = eigSol.eigenvalues().data();
    vector<cplx> eigVal;
    for (int i=0; i<sizeOfBasis(); i++)
      eigVal.push_back(eigSol.eigenvalues().data()[i]);
    
    std::sort(eigVal.begin(), eigVal.end(), hasSmallerNorm);
    
    for (int i=0; i<sizeOfBasis(); i++)
      cout << eigVal[i] << endl;
    
    ofstream outFinalGap("output/Galerkin/finalGap.txt",  std::ofstream::app);
    outFinalGap << gamma() << " " << nbOfFourier() << " " << nbOfHermite() << " " << abs(eigVal[0]) << " " << abs(eigVal[1]) << endl;
    
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
    createLeq();
  }
  
  void OverdampedGalerkin::createLeq()
  {
    basis_(0)->laplacianMatrix(Leq_);
    Leq_ /= -beta_;
    
    /*for (int iOfHermite = 0; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      for (int jOfHermite = 0; jOfHermite < (int)nbOfHermite_; jOfHermite++)
        Leq_(iOfHermite, iOfHermite) = - 1 / beta_ * basis_(0)->xLaplacianY(iOfHermite+1, jOfHermite+1);*/
  }
  
  void OverdampedGalerkin::computeExpToTrigTens()
  {
    trigToExpTens_ = trigToExpMat_;
    expToTrigTens_ = expToTrigMat_;
  }
  
  void OverdampedGalerkin::compute()
  {
    cout << "start OverdampedGalerkin::compute()" << endl;
    cout << "############ Leq ############" << endl;
    //cout << Leq_ << endl << endl;
    
    cout << "Leq_ size : " << Leq_.numberOfRows() << " x " << Leq_.numberOfColumns() << endl;
    cout << "Leq_ nnz : " << Leq_.nonZeros() << endl;

    DenseMatrix<double> DLeq = Leq_.dense();
    display(Leq_, "output/Galerkin/Leq");


    //DenseMatrix<double> DLeqSad = shapeSaddle(DLeq);
    //display(DLeqSad, "output/Galerkin/LeqSad");

    cout << "############ DLeqSad ############" << endl;
    //cout << DLeqSad << endl << endl;

    //DenseMatrix<double> IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);

    cout << "Computing DLeqInv..."; cout.flush();
    //DenseMatrix<double> LeqInvSad = inv(DLeqSad);
    DenseMatrix<double> DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;
    
    cout << "Norm2gVector : " << norm2gVector() << endl;
    cout << "Norm of LeqinvSad : " << DLeqInv.norm() << endl;
    cout << "Norm of Leqinv : " << DLeq.inverse().norm() << endl;
    //cout << "Norm of Sinv Leqinv : " << shapePrec(DLeq).inverse().norm() << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DenseMatrix<double> LeqInvSad = inv(DLeqSad);
    DenseMatrix<double> LeqInv = invWithSaddle(Leq_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");
    
    cout << "L g = " << Leq_ * gVector() << endl << endl;
    
    DVec gradV = gVector(); //getGradV();
    display(gradV, "output/Galerkin/gradV");
    
    DVec LinvGradV = solve(Leq_, gradV);
    display(LinvGradV, "output/Galerkin/LinvGradV");
    
    DVec LinvGradVsad = solveWithSaddle(Leq_, gradV);
    display(LinvGradVsad, "output/Galerkin/LinvGradVsad");
  }
  
  DVec OverdampedGalerkin::getGradV() const
  {
    DVec gradV(sizeOfBasis());
    gradV(0) = 0;
    for (int iOfFourier = 1; iOfFourier < nbOfFourier(); iOfFourier++)
      if (iOfFourier % 2 == 0)  gradV(iOfFourier) = iOfFourier / beta_ * expFourierMeans(iOfFourier-1);
      else gradV(iOfFourier) = -(iOfFourier+1) / beta_ * expFourierMeans(iOfFourier+1);
    
    assert(dot(gradV, gVector()) == 0);
    return gradV;
  }


  //#### LangevinGalerkin ####


  LangevinGalerkin::LangevinGalerkin(Input const& input):
    Galerkin(input)
  {
    computeExpToTrigTens();
    
    cout << "############ Lham ############" << endl;
    Lham_ = kron(Q_, tP_) - kron(tQ_, P_);
    display(Lham_, "output/Galerkin/Lham");

    cout << "############ Lthm ############" << endl;
    createLthm();
    display(Lthm_, "output/Galerkin/Lthm");

    cout << "############ Leq ############" << endl;
    Leq_  = Lham_ + gamma_ * Lthm_;

    cout << "############ L1 ############" << endl;
    //L1_ = kron(SMat::Identity(nbOfFourier_) , P_);
    L1_ = kron(SIdQ_, P_);

    cout << "############ Leta ############" << endl;
    Leta_ = Leq_ + externalForce_ * L1_;
  }
  
  void LangevinGalerkin::createLthm0()
  {
    basis_(1)->laplacianMatrix(Lthm0_);
    Lthm0_ /= -beta_;
    /*cout << "createLthm0" << endl;
    for (int iOfHermite = 1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
      Lthm0_(iOfHermite, iOfHermite) = - iOfHermite;
    cout << "end createLthm0" << endl;*/
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
    //SMat tempId = speye(PMat.numberOfRows(), PMat.numberOfColumns());
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
    /*assert(iOfVariableA < iOfVariableB);
    assert(iOfVariableB < nbOfVariables);
    assert(A.size() == B.size());*/
    //SMat tempId = speye(A.numberOfRows(), A.numberOfColumns());
    SMat res = speye<double>(1, 1);
    for (int i = 0; i < nbOfParticles_; i++)
    {
      if (i == iOfParticleQ) res = kron(res, QMat);
      else res = kron(res, SIdQ_);
      if (i == iOfParticleP) res = kron(res, PMat);
      else res = kron(res, SIdP_);
    }
    cout << "res : " << res.numberOfRows() << " x " << res.numberOfColumns() << endl;
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
      Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
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
    trigToExpTens_ = simol::kron(trigToExpMat_, DIdP_);
    trigToExpTens_ = tensorPower(trigToExpTens_, nbOfParticles_);
    expToTrigTens_ = simol::kron(expToTrigMat_, DIdP_);
    expToTrigTens_ = tensorPower(expToTrigTens_, nbOfParticles_);
    
    display(trigToExpTens_, "output/Galerkin/trigToExpTens");
  }


  void BoundaryLangevinGalerkin::compute()
  {
    cout << "start Galerkin::compute()" << endl;
    cout << "############ Leq ############" << endl;
    //cout << Leq_ << endl << endl;
    
    cout << "Leq_ size : " << Leq_.numberOfRows() << " x " << Leq_.numberOfColumns() << endl;
    cout << "Leq_ nnz : " << Leq_.nonZeros() << endl;

    DenseMatrix<double> DLeq = Leq_.dense();
    display(Leq_, "output/Galerkin/Leq");

    //DenseMatrix<double> DLeqSad = shapeSaddle(DLeq);
    //display(DLeqSad, "output/Galerkin/LeqSad");

    cout << "############ DLeqSad ############" << endl;
    //cout << DLeqSad << endl << endl;

    //DenseMatrix<double> IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);

    cout << "Computing DLeqInv..."; cout.flush();
    //DenseMatrix<double> LeqInvSad = inv(DLeqSad);
    DenseMatrix<double> DLeqInv = invWithSaddle(DLeq);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");
    
    cout << "Computing LeqInv..."; cout.flush();
    //DenseMatrix<double> LeqInvSad = inv(DLeqSad);
    DenseMatrix<double> LeqInv = invWithSaddle(Leq_);
    cout << "OK" << endl;
    display(DLeqInv, "output/Galerkin/DLeqInv");

    DVec N0H2Trig = Vector<double>::Zero(sizeOfBasis_);
    N0H2Trig(iTens(0, 2, 0)) = 1;
    display(N0H2Trig, "output/Galerkin/N0H2Trig");
    DVec N0H2 = trigToExpTens_ * N0H2Trig;
    display(N0H2, "output/Galerkin/N0H2");
    
    //DenseMatrix<double> N0H2TrigMat(N0H2Trig, nbOfFourier_, nbOfHermite_); 
    //DenseMatrix<double> N0H2Mat = trigToExpMat_ * N0H2TrigMat;
    
    DenseMatrix<double> N0H2Mat(N0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));
    
    display(N0H2Mat, "output/Galerkin/N0H2Mat");
    
    DVec LinvN0H2 = solveWithSaddle(Leq_, N0H2);
    DenseMatrix<double> LinvN0H2Mat(LinvN0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));
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

    /*DVec LinvL1LinvH1 = solveWithSaddle(Leq_, L1_ * LinvH1);
    //double varCoeff = -.5 * dot( L1_ * LinvH1, LeqInv * L1_ * LinvH1);
    double varCoeff = -2 * dot( L1_ * LinvH1, LinvL1LinvH1);
    cout << "varCoeff = " << varCoeff << endl;*/
  }



}



