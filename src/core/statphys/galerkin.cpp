#include "galerkin.hpp"

using std::cout; 
using std::endl; 

namespace simol
{
	
	SpMat tensor(SpMat A, SpMat B)
	{
		SpMat C(A.rows()*B.rows(), A.cols()*B.cols());
		for (int j=0; j < A.cols(); j++)
			for (SpMat::InnerIterator it(A,j); it; ++it)
				for (int j2=0; j2 < B.cols(); j2++)
					for (SpMat::InnerIterator it2(B,j2); it2; ++it2)
						C.insert(B.cols() * it.row() + it2.row(), B.rows() * j + j2) = it.value()*it2.value();
		return C;
	}
	
	SpMat identity(size_t nbRow)
	{
		SpMat A(nbRow, nbRow);
		for (size_t i=0; i<nbRow; i++)
			A.insert(i,i) = cplx(1.,0.);
		return A;
	}
	
	using namespace Eigen;
	
	SpMat inverse(SpMat& A)
	{
		cout << "start inverse" << endl;
		assert(A.rows() == A.cols());
		A.makeCompressed();
		Eigen::Matrix<cplx, Eigen::Dynamic, 1> x(A.rows()), b(A.rows());
		//Eigen::SparseMatrix<cplx, A.rows()> A;
		Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> >   solver;
		// fill A and b;
		// Compute the ordering permutation vector from the structural pattern of A
		solver.analyzePattern(A); 
		// Compute the numerical factorization 
		solver.factorize(A); 
		cout << "start solve" << endl;
		//Use the factors to solve the linear system 
		x = solver.solve(b);
		cout << x << endl;
		
		/*SpMat I(A.rows(),A.cols());
		I.setIdentity();
		return solver.solve(I);*/
		
		/*SparseLU<SpMat > solver;
		<SparseMatrix<scalar, ColMajor>, COLAMDOrdering<Index> >
		solver.compute(A);
		SpMat I(A.rows(),A.cols());
		I.setIdentity();
		return solver.solve(I);*/
		//auto A_inv = solver.solve(I);
	}
	
	Galerkin::Galerkin(Input const& input):				//ex : [-2:2]
		numberOfFourier_(input.numberOfFourier()),	//ex : 5
		numberOfHermite_(input.numberOfHermite()),
		maxOfFourier_((numberOfFourier_-1)/2),			//ex : 2
		Q_(numberOfFourier_, numberOfFourier_),
		P_(numberOfHermite_, numberOfHermite_),
		Lthm0_(numberOfHermite_, numberOfHermite_),
		Lthm_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		Lham_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		Leq_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		beta_(input.beta())
	{
		assert(numberOfFourier_ % 2 == 1);
		cout << maxOfFourier_ << " " << numberOfFourier_ << endl;
		for (int iOfFourier=-maxOfFourier_; iOfFourier <= maxOfFourier_; iOfFourier++)
		{
			cout << iPosiF(iOfFourier) << " " << cplx(0., iOfFourier) << endl;
			Q_.insert(iPosiF(iOfFourier), iPosiF(iOfFourier)) = cplx(0., iOfFourier);
		}
		
		for (int iOfHermite=0; iOfHermite < (int)numberOfHermite_-1; iOfHermite++)
			P_.insert(iOfHermite, iOfHermite+1) = cplx(sqrt(beta_) * sqrt(iOfHermite+1.), 0.);
		
		Lham_ = tensor(Q_, P_) - tensor(Q_.transpose(), P_.transpose());
		
		cout << "############ Lham ############" << endl;
		cout << Lham_ << endl << endl;
		
		for (int iOfHermite=1; iOfHermite < (int)numberOfHermite_; iOfHermite++)
			Lthm0_.insert(iOfHermite, iOfHermite) = cplx(-iOfHermite, 0.);
		
		SpMat IdQ = identity(numberOfFourier_);
		Lthm_ = tensor(IdQ, Lthm0_);
		
		cout << "############ IdQ ############" << endl;
		cout << IdQ << endl;
		
		cout << "############ Lthm0 ############" << endl;
		cout << Lthm0_ << endl << endl;
		
		cout << "############ Lthm ############" << endl;
		cout << Lthm_ << endl << endl;
		
		Leq_  = Lham_ + Lthm_;
	}
	
	size_t Galerkin::iPosi(int iOfFourier, size_t iOfHermite) const
	{
		return numberOfHermite_ * (iOfFourier + maxOfFourier_) + iOfHermite;
	}
	
	size_t Galerkin::iSign(int iOfFourier, size_t iOfHermite) const
	{
		return numberOfHermite_ * (iOfFourier + maxOfFourier_) + iOfHermite;
	}
	
	size_t Galerkin::iPosiF(int iOfFourier) const
	{
		return iOfFourier + maxOfFourier_;
	}
	
	size_t Galerkin::iSignF(int iOfFourier) const
	{
		return iOfFourier - maxOfFourier_;
	}
	
	void Galerkin::solve()
	{
		DsMat Linv = inverse(Leq_);
		cout << "############ Linv ############" << endl;
		cout << Linv << endl << endl;
	}
	
}