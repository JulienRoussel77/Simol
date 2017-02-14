#include "simol/statphys/TensorTools.hpp"



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
  Eigen::SparseLU<SMat> solver(A);
  //Eigen::GMRES<SMat, Eigen::IncompleteLUT<double>> solver(A); //solver.setTolerance(1e-15);
  //Eigen::GMRES<SMat, Eigen::IdentityPreconditioner> solver(A);
  //Eigen::BiCGSTAB<SMat, Eigen::IncompleteLUT<double>> solver(A);

  DVec X = DVec::Random(A.rows());
  DVec Y;
  double eigVal = X.norm();
  double eigValDiff = 1;
  double prevEigValDiff = 2;
  
  while (fabs(eigValDiff) > tol && (fabs(eigValDiff) < fabs(prevEigValDiff) || nbOfIter < 100))
  {
    Y = X / eigVal;
    X = solver.solve(Y);
    //X = solver.solveWithGuess(Y, Y * eigVal);
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