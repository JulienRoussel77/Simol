#include "galerkin.hpp"

using std::cout; 
using std::endl; 
using std::ostream; 

using namespace arma;

namespace simol
{
	double Galerkin::potential(double q)
	{
		return amplitude_ * (1 - cos(q));
	}
	
	DMat Galerkin::shapeSaddle(DMat& A)
	{
		//DMat Asad(A.n_rows+1, A.n_cols+1, fill::ones);
		DMat Asad(A.n_rows+1, A.n_cols+1, fill::zeros);
		Asad.submat(0, 0, A.n_rows-1, A.n_cols-1) = A;
		//Asad(A.n_rows, A.n_cols) = 0;
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
		{
			Asad(A.n_rows, iTens(iOfFourier2, 0)) = fourierCoeffsExp_[iOfFourier2];
			Asad(iTens(iOfFourier2, 0), A.n_cols) = fourierCoeffsExp_[iOfFourier2];
		}
		return Asad;
	}
	
	DMat Galerkin::unshapeSaddle(DMat& Asad)
	{
		return Asad.submat(0, 0, Asad.n_rows-2, Asad.n_cols-2);
	}
	
	DVec Galerkin::shapeSaddle(DVec& X)
	{
		DVec Xsad(X.n_rows+1, fill::zeros);
		Xsad.subvec(0, X.n_rows-1) = X;
		return Xsad;
	}
	
	DVec Galerkin::unshapeSaddle(DVec& Xsad)
	{
		return Xsad.subvec(0, Xsad.n_rows-2);
	}
	
	void Galerkin::computeFourierCoeffsExp()
	{		
		double step = 2 * M_PI / (double)nbIntegrationNodes_;
		for (size_t iOfNode = 0; iOfNode < nbIntegrationNodes_; iOfNode++)
		{
			double q = - M_PI + iOfNode * step;
			fourierCoeffsExp_[0] += exp(-beta_ * potential(q)/2);
			for (size_t iOfFourier=1; iOfFourier <= maxOfFourier_; iOfFourier++)
				fourierCoeffsExp_[2 * iOfFourier] += sqrt(2) * cos(iOfFourier * q) * exp(-beta_ * potential(q)/2);
			for (size_t iOfFourier=1; iOfFourier <= maxOfFourier_; iOfFourier++)
				fourierCoeffsExp_[2 * iOfFourier-1] += sqrt(2) * sin(iOfFourier * q) * exp(-beta_ * potential(q)/2);
		}
		for (size_t iOfFourier2=0; iOfFourier2 <= 2 * maxOfFourier_; iOfFourier2++)
		{
			fourierCoeffsExp_[iOfFourier2] *= step / (2*M_PI);
			//cout << "coeff n°"<< iOfFourier2 << " = " << fourierCoeffsExp_[iOfFourier2] << endl;
		}
	}
	
	DMat inverse(SMat& A)
	{
		//DMat Id = eye<DMat>(A.size());
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = spsolve(A, Id);
		return C;
	}
	
		DMat inverse(DMat& A)
	{
		//DMat Id = eye<DMat>(A.size());
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = solve(A, Id);
		return C;
	}
	
	// C = ( A B_11    A B_12   ... 
	//			 A B_21		 A B_22   ...
	//															)
	SMat kron(SMat& A, SMat& B)
	{
		double a_1 = A.n_cols * (B.n_cols - 1)/(A.n_cols - 1);
		double b_1 = a_1 - B.n_cols;
		double a_2 = A.n_rows * (B.n_rows - 1)/(A.n_rows - 1);
		double b_2 = a_2 - B.n_rows;
		cout << "A : O <= i < " << A.n_rows << ", O <= jOfA < " << A.n_cols << endl;
		cout << "B : O <= i2 < " << B.n_rows << ", O <= jOfB < " << B.n_cols << endl;
		cout << "We keep the coefficients such that j * (jOfB + " << b_1 << ") <= " << a_1 
			<< " and such that i * (i2 + " << b_2 << ") <= " << a_2 << endl; 
		SMat C(A.n_rows*B.n_rows, A.n_cols*B.n_cols);
		//cout << A.size() << endl << B.size() << endl;
		//SMat C(A.size() % B.size());					//element-wise product of the dimensions
		for (int jOfA=0; (size_t) jOfA < A.n_cols; jOfA++)
			for (SMat::iterator it = A.begin_col(jOfA); it != A.end_col(jOfA); ++it)
			{
				int iOfA = it.row();
				double valOfA = *it;
				for (int jOfB=0; (size_t) jOfB < B.n_cols; jOfB++)
					for (SMat::iterator it2 = B.begin_col(jOfB); it2 != B.end_col(jOfB); ++it2)
					{
						int iOfB = it2.row();
						double valOfB = *it2;
						//cout << iOfA << " " << j << " " << iOfB << " " << jOfB << endl;
						//if (jOfA * (b_1 + jOfB) <= a_1 && iOfA * (b_2 + iOfB) <= a_2)
						//{						
							C(iOfA + A.n_rows * iOfB, jOfA + A.n_cols * jOfB) = valOfA*valOfB;
							//cout << "ok"<<endl;
						//}
					}
			}
		return C;
	}

	void display(cx_vec& X, ostream& out)
	{
		for (int i=0; (size_t) i < X.n_rows; i++)
			out << real(X(i)) << " " << imag(X(i)) << endl;
	}
	
	void displayWithNan(DMat& A, ostream& out)
	{
		for (int i=0; (size_t) i < A.n_rows; i++)
		{
			for (int j=0; (size_t) j < A.n_cols; j++)			
			{
				if (true)//fabs(A(i,j)) > 1e-15)
				{
					//cout << i << " " << j << " " << A(i,j) << endl;
					out << A(i,j) << " ";
				}
				else
					out << "nan ";
			}
			out << endl;
		}
	}
	
	void displayWithNan(SMat& A, ostream& out)
	{
		DMat DA = conv_to<DMat>::from(A);
		displayWithNan(DA, out);
	}
	
	/*void displaySparse(SMat& A, ofstream& out)
	{
		for (int j=0; (size_t) j < A.n_cols; j++)
			for (SMat::iterator it = A.begin_col(j); it != A.end_col(j); ++it)
				out << it.row() << " " << j << " " << *it << endl;
	}*/
	
	using namespace Eigen;
	
	Galerkin::Galerkin(Input const& input):				//ex : [0:4]
		numberOfFourier_(input.numberOfFourier()),	//ex : 5
		numberOfHermite_(input.numberOfHermite()),
		maxOfFourier_((numberOfFourier_-1)/2),			//ex : 2
		sizeOfBasis_(numberOfFourier_ * numberOfHermite_),
		Q_(numberOfFourier_, numberOfFourier_),
		P_(numberOfHermite_, numberOfHermite_),
		Lthm0_(numberOfHermite_, numberOfHermite_),
		Lthm_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		Lham_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		Leq_(numberOfFourier_*numberOfHermite_, numberOfFourier_*numberOfHermite_),
		beta_(input.beta()),
		amplitude_(input.amplitude()),
		nbIntegrationNodes_(1000),
		fourierCoeffsExp_(1000, 0)
	{		
		assert(numberOfFourier_ % 2 == 1);
		cout << endl << "Number of modes : " << numberOfFourier_ << " x " << numberOfHermite_ << endl;
		
		computeFourierCoeffsExp();
		
		Q_(1,0) = amplitude_ * beta_ / 2;
		for (int iOfFourier=1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
		{
			Q_(2 * iOfFourier - 2, 2 * iOfFourier - 1) = amplitude_ * beta_ / 4;
			Q_(2 * iOfFourier    , 2 * iOfFourier - 1) = iOfFourier;
			if (iOfFourier != (int) maxOfFourier_)
				Q_(2 * iOfFourier + 2, 2 * iOfFourier - 1) = - amplitude_ * beta_ / 4;
		}
		for (int iOfFourier=1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
		{
			if (iOfFourier != 1)
				Q_(2 * iOfFourier - 3, 2 * iOfFourier) = - amplitude_ * beta_ / 4;
			Q_(2 * iOfFourier - 1, 2 * iOfFourier) = - iOfFourier;
			if (iOfFourier != (int) maxOfFourier_)
				Q_(2 * iOfFourier + 1, 2 * iOfFourier) =   amplitude_ * beta_ / 4;
		}
		tQ_ = Q_.t();
		
		cout << "############ Q ############" << endl;
		//cout << Q_ << endl << endl;
		ofstream out_Q("../output/Galerkin/Q");
		//out_Q << Q_ << endl;
		DMat DQ = conv_to<DMat>::from(Q_);
		displayWithNan(DQ, out_Q);
		
		for (size_t iOfHermite=1; iOfHermite < numberOfHermite_; iOfHermite++)
			//P_(iOfHermite, iOfHermite+1) = cplx(sqrt(beta_) * sqrt(iOfHermite+1.), 0.);
			P_(iOfHermite-1, iOfHermite) = sqrt(beta_*iOfHermite);
		
		tP_ = P_.t();
		
		cout << "############ P ############" << endl;
		cout << P_ << endl << endl;
		ofstream out_P("../output/Galerkin/P");
		//out_P << P_ << endl;
		DMat DP = conv_to<DMat>::from(P_);
		displayWithNan(DP, out_P);
		
		Lham_ = kron(Q_, tP_) - kron(tQ_, P_);

		
		cout << "############ Lham ############" << endl;
		//cout << Lham_ << endl << endl;
		ofstream out_Lham("../output/Galerkin/Lham");
		//out_Leq << DLeq << endl;
		displayWithNan(Lham_, out_Lham);
		
		for (int iOfHermite=1; iOfHermite < (int)numberOfHermite_; iOfHermite++)
			//Lthm0_(iOfHermite, iOfHermite) = cplx(-iOfHermite, 0.);
			Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
		
		SMat IdQ = speye<SMat>(numberOfFourier_, numberOfFourier_);
		Lthm_ = kron(IdQ, Lthm0_);
		
		cout << "############ Lthm0 ############" << endl;
		//cout << Lthm0_ << endl << endl;
		
		cout << "############ Lthm ############" << endl;
		//cout << Lthm_ << endl << endl;
		ofstream out_Lthm("../output/Galerkin/Lthm");
		//out_Leq << DLeq << endl;
		displayWithNan(Lthm_, out_Lthm);
		
		Leq_  = Lham_ + Lthm_;
	}
	
	//psi = (1,1  1,2  ...  1,N_H  2,1 ... )
	//N_H blocks of size N_G (we concatene the columns of the matrix)
	size_t Galerkin::iTens(size_t iOfFourier2, size_t iOfHermite) const
	{
		assert(iOfFourier2 < numberOfFourier_	&& iOfHermite < numberOfHermite_);
		return numberOfFourier_ * iOfHermite + iOfFourier2;
	}
	
	
	void Galerkin::compute()
	{
		cout << "start Galerkin::compute()" << endl;
		cout << "############ Leq ############" << endl;
		//cout << Leq_ << endl << endl;
		
		DMat DLeq = conv_to<DMat>::from(Leq_);
		ofstream out_Leq("../output/Galerkin/Leq");
		//out_Leq << DLeq << endl;
		displayWithNan(DLeq, out_Leq);
		
		DMat DLeqSad = shapeSaddle(DLeq);
		ofstream out_LeqSad("../output/Galerkin/LeqSad");
		displayWithNan(DLeqSad, out_LeqSad);
		
		cout << "############ DLeqSad ############" << endl;
		//cout << DLeqSad << endl << endl;
		
		DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);
		//DMat LeqInvSad = solve(DLeq, IdSad);
		DMat LeqInvSad = inv(DLeqSad);
		DMat LeqInv = LeqInvSad.submat(0, 0, sizeOfBasis_-1, sizeOfBasis_-1);
		
		cout << "############ LeqInv ############" << endl;
		//cout << LeqInv << endl << endl;
		ofstream out_LeqInv("../output/Galerkin/LeqInv");
		//out_LeqInv << LeqInv << endl;
		displayWithNan(LeqInv, out_LeqInv);

		
		DVec H1(sizeOfBasis_);
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
			H1(iTens(iOfFourier2, 1)) = fourierCoeffsExp_[iOfFourier2];
		
		displayWithNan(H1, cout);
		
		DVec H1sad = shapeSaddle(H1);
		
		DMat H1Mat = reshape(H1, numberOfFourier_, numberOfHermite_);
		cout << "############ H1Mat ############" << endl;
		//cout << H1Mat << endl << endl;
		ofstream out_H1("../output/Galerkin/H1");
		//out_H1 << H1Mat << endl;
		displayWithNan(H1Mat, out_H1);
		
		//SMat IdFourier = speye(numberOfFourier_, numberOfFourier_);
		//SMat test = kron(IdFourier, P_);
		
		DVec LH1 = Lthm_ * H1;
		//DVec LH1 = test * H1;
		
		DMat LH1Mat = reshape(LH1, numberOfFourier_, numberOfHermite_);
		ofstream out_LH1("../output/Galerkin/LH1");
		displayWithNan(LH1Mat, out_LH1);
		
		
		DVec H1invSad = solve(DLeqSad, H1sad);
		cout << "solve ok !" << endl;
		DVec H1inv = unshapeSaddle(H1invSad);
		
		double lambda = H1invSad(sizeOfBasis_);
		
		cout << "############ lambda ############" << endl;
		cout << lambda << endl << endl;
		
		DMat H1invMat = reshape(H1inv, numberOfFourier_, numberOfHermite_);
		
		cout << "############ H1invMat ############" << endl;
		ofstream out_H1inv("../output/Galerkin/H1inv");
		//out_H1inv << H1invMat << endl;
		displayWithNan(H1invMat, out_H1inv);
		
		DVec H1back = Leq_ * H1inv;
		ofstream out_H1back("../output/Galerkin/H1back");
		DMat H1backMat = reshape(H1back, numberOfFourier_, numberOfHermite_);
		displayWithNan(H1backMat, out_H1back);
		
		
		cx_vec eigvalLeq;
		cx_mat eigvecLeq;
		eig_gen(eigvalLeq, eigvecLeq, -DLeq);
		//eigs_gen(eigvalLeq, eigvecLeq, Leq_, 20);
		
		ofstream out_eigvalLeq("../output/Galerkin/eigvalLeq");
		//out_eigvalLeq << eigvalLeq << endl;
		display(eigvalLeq, out_eigvalLeq);
		
		ofstream out_eigvecLeq0("../output/Galerkin/eigvecLeq0");
		DMat eigvecLeq0 = reshape(real(eigvecLeq.col(eigvecLeq.n_cols-1)), numberOfFourier_, numberOfHermite_);
		//display(eigvecLeq0, out_eigvecLeq0);
		out_eigvecLeq0 << eigvecLeq0 << endl;
		
		ofstream out_eigvecLeq1("../output/Galerkin/eigvecLeq1");
		DMat eigvecLeq1 = reshape(real(eigvecLeq.col(eigvecLeq.n_cols-2)), numberOfFourier_, numberOfHermite_);
		//cx_vec eigvecLeq1 = eigvecLeq.col(1);
		//display(eigvecLeq1, out_eigvecLeq1);
		out_eigvecLeq1 << eigvecLeq1 << endl;
		

		//DMat Linv = inverse(DsLeq);
		
		//cout << "############ Linv ############" << endl;
		//cout << Linv << endl << endl;
		
		//ofstream out_re("../output/Galerkin/Linv_re");
		//out_re << real(Linv) << endl;		
		//ofstream out_im("../output/Galerkin/Linv_im");
		//out_im << imag(Linv) << endl;
		
		//ofstream out("../output/Galerkin/Linv");
		//out << Linv << endl;	
	}
	
}