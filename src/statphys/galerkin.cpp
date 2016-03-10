#include "galerkin.hpp"

using namespace arma;

namespace simol
{
	double Galerkin::potential(double q)
	{ return amplitude_ * (1 - cos(q)); }

	DMat Galerkin::shapeSaddle(DMat& A)
	{
		DMat Asad(A.n_rows+1, A.n_cols+1, fill::zeros);
		Asad.submat(0, 0, A.n_rows-1, A.n_cols-1) = A;
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
		{
			Asad(A.n_rows, iTens(iOfFourier2, 0)) = fourierCoeffsExp_[iOfFourier2];
			Asad(iTens(iOfFourier2, 0), A.n_cols) = fourierCoeffsExp_[iOfFourier2];
		}
		return Asad;
	}

	DMat Galerkin::unshapeSaddle(DMat& Asad)
	{ return Asad.submat(0, 0, Asad.n_rows-2, Asad.n_cols-2); }

	DVec Galerkin::shapeSaddle(DVec& X)
	{
		DVec Xsad(X.n_rows+1, fill::zeros);
		Xsad.subvec(0, X.n_rows-1) = X;
		return Xsad;
	}

	DVec Galerkin::unshapeSaddle(DVec& Xsad)
	{ return Xsad.subvec(0, Xsad.n_rows-2); }

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
		}
	}

	DMat inverse(SMat& A)
	{
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = spsolve(A, Id);
		return C;
	}

		DMat inverse(DMat& A)
	{
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = solve(A, Id);
		return C;
	}

	SMat kron(SMat& A, SMat& B)
	{
		double a_1 = A.n_cols * (B.n_cols - 1)/(A.n_cols - 1);
		double b_1 = a_1 - B.n_cols;
		double a_2 = A.n_rows * (B.n_rows - 1)/(A.n_rows - 1);
		double b_2 = a_2 - B.n_rows;
        std::cout << "A : O <= i < " << A.n_rows << ", O <= jOfA < " << A.n_cols << std::endl;
        std::cout << "B : O <= i2 < " << B.n_rows << ", O <= jOfB < " << B.n_cols << std::endl;
        std::cout << "We keep the coefficients such that j * (jOfB + " << b_1 << ") <= " << a_1
			<< " and such that i * (i2 + " << b_2 << ") <= " << a_2 << std::endl;
		SMat C(A.n_rows*B.n_rows, A.n_cols*B.n_cols);
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
							C(iOfA + A.n_rows * iOfB, jOfA + A.n_cols * jOfB) = valOfA*valOfB;
					}
			}
		return C;
	}

	void display(cx_vec& X, std::ostream& out)
	{
		for (int i=0; (size_t) i < X.n_rows; i++)
			out << real(X(i)) << " " << imag(X(i)) << std::endl;
	}

	void displayWithNan(DMat& A, std::ostream& out)
	{
		for (int i=0; (size_t) i < A.n_rows; i++)
		{
			for (int j=0; (size_t) j < A.n_cols; j++)
			{
				if (true)//fabs(A(i,j)) > 1e-15)
				{
					out << A(i,j) << " ";
				}
				else
					out << "nan ";
			}
			out << std::endl;
		}
	}

	void displayWithNan(SMat& A, std::ostream& out)
	{
		DMat DA = conv_to<DMat>::from(A);
		displayWithNan(DA, out);
	}

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
        std::cout << std::endl << "Number of modes : " << numberOfFourier_ << " x " << numberOfHermite_ << std::endl;

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

        std::cout << "############ Q ############" << std::endl;
		//std::cout << Q_ << std::endl << std::endl;
        std::ofstream out_Q("../output/Galerkin/Q");
		//out_Q << Q_ << std::endl;
		DMat DQ = conv_to<DMat>::from(Q_);
		displayWithNan(DQ, out_Q);

		for (size_t iOfHermite=1; iOfHermite < numberOfHermite_; iOfHermite++)
			//P_(iOfHermite, iOfHermite+1) = cplx(sqrt(beta_) * sqrt(iOfHermite+1.), 0.);
			P_(iOfHermite-1, iOfHermite) = sqrt(beta_*iOfHermite);

		tP_ = P_.t();

        std::cout << "############ P ############" << std::endl;
        std::cout << P_ << std::endl << std::endl;
        std::ofstream out_P("../output/Galerkin/P");
		//out_P << P_ << std::endl;
		DMat DP = conv_to<DMat>::from(P_);
		displayWithNan(DP, out_P);

		Lham_ = kron(Q_, tP_) - kron(tQ_, P_);


        std::cout << "############ Lham ############" << std::endl;
		//std::cout << Lham_ << std::endl << std::endl;
        std::ofstream out_Lham("../output/Galerkin/Lham");
		//out_Leq << DLeq << std::endl;
		displayWithNan(Lham_, out_Lham);

		for (int iOfHermite=1; iOfHermite < (int)numberOfHermite_; iOfHermite++)
			//Lthm0_(iOfHermite, iOfHermite) = cplx(-iOfHermite, 0.);
			Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;

		SMat IdQ = speye<SMat>(numberOfFourier_, numberOfFourier_);
		Lthm_ = kron(IdQ, Lthm0_);

        std::cout << "############ Lthm0 ############" << std::endl;
		//std::cout << Lthm0_ << std::endl << std::endl;

        std::cout << "############ Lthm ############" << std::endl;
		//std::cout << Lthm_ << std::endl << std::endl;
        std::ofstream out_Lthm("../output/Galerkin/Lthm");
		//out_Leq << DLeq << std::endl;
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
        std::cout << "start Galerkin::compute()" << std::endl;
        std::cout << "############ Leq ############" << std::endl;

		DMat DLeq = conv_to<DMat>::from(Leq_);
        std::ofstream out_Leq("../output/Galerkin/Leq");
		displayWithNan(DLeq, out_Leq);

		DMat DLeqSad = shapeSaddle(DLeq);
        std::ofstream out_LeqSad("../output/Galerkin/LeqSad");
		displayWithNan(DLeqSad, out_LeqSad);

        std::cout << "############ DLeqSad ############" << std::endl;

		DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);
		DMat LeqInvSad = inv(DLeqSad);
		DMat LeqInv = LeqInvSad.submat(0, 0, sizeOfBasis_-1, sizeOfBasis_-1);

        std::cout << "############ LeqInv ############" << std::endl;
        std::ofstream out_LeqInv("../output/Galerkin/LeqInv");
		displayWithNan(LeqInv, out_LeqInv);


		DVec H1(sizeOfBasis_);
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
			H1(iTens(iOfFourier2, 1)) = fourierCoeffsExp_[iOfFourier2];

		displayWithNan(H1, std::cout);

		DVec H1sad = shapeSaddle(H1);

		DMat H1Mat = reshape(H1, numberOfFourier_, numberOfHermite_);
        std::cout << "############ H1Mat ############" << std::endl;
        std::ofstream out_H1("../output/Galerkin/H1");
		displayWithNan(H1Mat, out_H1);

		DVec LH1 = Lthm_ * H1;

		DMat LH1Mat = reshape(LH1, numberOfFourier_, numberOfHermite_);
        std::ofstream out_LH1("../output/Galerkin/LH1");
		displayWithNan(LH1Mat, out_LH1);


		DVec H1invSad = solve(DLeqSad, H1sad);
        std::cout << "solve ok !" << std::endl;
		DVec H1inv = unshapeSaddle(H1invSad);

		double lambda = H1invSad(sizeOfBasis_);

        std::cout << "############ lambda ############" << std::endl;
        std::cout << lambda << std::endl << std::endl;

		DMat H1invMat = reshape(H1inv, numberOfFourier_, numberOfHermite_);

        std::cout << "############ H1invMat ############" << std::endl;
        std::ofstream out_H1inv("../output/Galerkin/H1inv");
		displayWithNan(H1invMat, out_H1inv);

		DVec H1back = Leq_ * H1inv;
        std::ofstream out_H1back("../output/Galerkin/H1back");
		DMat H1backMat = reshape(H1back, numberOfFourier_, numberOfHermite_);
		displayWithNan(H1backMat, out_H1back);


		cx_vec eigvalLeq;
		cx_mat eigvecLeq;
		eig_gen(eigvalLeq, eigvecLeq, -DLeq);

        std::ofstream out_eigvalLeq("../output/Galerkin/eigvalLeq");
		display(eigvalLeq, out_eigvalLeq);

        std::ofstream out_eigvecLeq0("../output/Galerkin/eigvecLeq0");
		DMat eigvecLeq0 = reshape(real(eigvecLeq.col(eigvecLeq.n_cols-1)), numberOfFourier_, numberOfHermite_);
		out_eigvecLeq0 << eigvecLeq0 << std::endl;

        std::ofstream out_eigvecLeq1("../output/Galerkin/eigvecLeq1");
		DMat eigvecLeq1 = reshape(real(eigvecLeq.col(eigvecLeq.n_cols-2)), numberOfFourier_, numberOfHermite_);
		out_eigvecLeq1 << eigvecLeq1 << std::endl;

	}

}
