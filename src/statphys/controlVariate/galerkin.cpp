#include "galerkin.hpp"

using std::cout;
using std::endl;
using std::ostream;

using namespace arma;


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
	
	DMat Galerkin::shapeSaddle(const DMat& A) const
	{
		//DMat Asad(A.n_rows+1, A.n_cols+1, fill::ones);
		DMat Asad(A.n_rows+1, A.n_cols+1, fill::zeros);
		Asad.submat(0, 0, A.n_rows-1, A.n_cols-1) = A;
		//Asad(A.n_rows, A.n_cols) = 0;
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
		{
			Asad(A.n_rows, iTens(iOfFourier2, 0)) = expFourierCoeffs(iOfFourier2);
			Asad(iTens(iOfFourier2, 0), A.n_cols) = expFourierCoeffs(iOfFourier2);
		}
		return Asad;
	}
	
	DMat Galerkin::unshapeSaddle(const DMat& Asad) const
	{
		return Asad.submat(0, 0, Asad.n_rows-2, Asad.n_cols-2);
	}
	
	DVec Galerkin::shapeSaddle(const DVec& X) const
	{
		DVec Xsad(X.n_rows+1, fill::zeros);
		Xsad.subvec(0, X.n_rows-1) = X;
		return Xsad;
	}
	
	DVec Galerkin::unshapeSaddle(const DVec& Xsad) const
	{
		return Xsad.subvec(0, Xsad.n_rows-2);
	}
	
	DMat inverse(const SMat& A)
	{
		//DMat Id = eye<DMat>(A.size());
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = spsolve(A, Id);
		return C;
	}
	
	DMat inverse(const DMat& A)
	{
		//DMat Id = eye<DMat>(A.size());
		DMat Id(A.n_rows, A.n_cols, fill::eye);
		DMat C = solve(A, Id);
		return C;
	}
	
	DVec Galerkin::solveWithSaddle(const SMat& A, const DVec& X) const
	{
		DMat DA = conv_to<DMat>::from(A);
		//DMat DA = conv<DMat>(A);
		return solveWithSaddle(DA, X);
	}
	
	DVec Galerkin::solveWithSaddle(const DMat& A, const DVec& X) const
	{
		cout << "solveWithSaddle...";
		DVec Xsad = shapeSaddle(X);
		DMat Asad = shapeSaddle(A);
		DVec Bsad = solve(Asad, Xsad);
		cout << "OK ! (lambda = " << Bsad(Bsad.size()-1) << ")" << endl;
		return unshapeSaddle(Bsad);
	}
	
	DMat Galerkin::invWithSaddle(const DMat& A) const
	{
		cout << "invWithSaddle...";
		DMat Asad = shapeSaddle(A);
		DMat Bsad = inv(Asad);
		//cout << "OK ! (lambda = " << Bsad(Bsad.size()-1) << ")" << endl;
		return unshapeSaddle(Bsad);
	}
	
	// C = ( A B_11    A B_12   ... 
	//			 A B_21		 A B_22   ...
	//															)
	SMat kron(const SMat& A, const SMat& B)
	{
		/*double a_1 = A.n_cols * (B.n_cols - 1)/(A.n_cols - 1);
		double b_1 = a_1 - B.n_cols;
		double a_2 = A.n_rows * (B.n_rows - 1)/(A.n_rows - 1);
		double b_2 = a_2 - B.n_rows;*/
		
		/*cout << "A : O <= i < " << A.n_rows << ", O <= jOfA < " << A.n_cols << endl;
		cout << "B : O <= i2 < " << B.n_rows << ", O <= jOfB < " << B.n_cols << endl;
		cout << "We keep the coefficients such that j * (jOfB + " << b_1 << ") <= " << a_1 
			<< " and such that i * (i2 + " << b_2 << ") <= " << a_2 << endl; */
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
// 							//cout << "ok"<<endl;
						//}
					}
			}
		return C;
	}
	
	DMat kron(const DMat& A, const DMat& B)
	{
		cout << "kron(DMat& A, DMat& B)" << endl;
		DMat C(A.n_rows*B.n_rows, A.n_cols*B.n_cols);
		for (int iOfA = 0; iOfA < (int) A.n_rows; iOfA++)
			for (int jOfA = 0; jOfA < (int) A.n_cols; jOfA++)
				for (int iOfB = 0; iOfB < (int) B.n_rows; iOfB++)
					for (int jOfB = 0; jOfB < (int) B.n_cols; jOfB++)
						C(iOfA + A.n_rows * iOfB, jOfA + A.n_cols * jOfB) = A(iOfA, jOfA) * B(iOfB, jOfB);
		return C;
	}

	void displayCplx(const cx_vec& X, ostream& out)
	{
		for (int i=0; (size_t) i < X.n_rows; i++)
			out << real(X(i)) << " " << imag(X(i)) << endl;
	}
	
	void displayMat(const DMat& A, ostream& out)
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
	
	void displayMat(const DMat& A, string path)
	{
		ofstream out(path);
		displayMat(A, out);
	}
	
	void displayMat(const SMat& A, ostream& out)
	{
		DMat DA = conv_to<DMat>::from(A);
		displayMat(DA, out);
	}
	
	void displayMat(const SMat& A, string path)
	{
		DMat DA = conv_to<DMat>::from(A);
		displayMat(DA, path);
	}
	
	//We compute the real Fourier coefficients of the function C^-1 exp(-\beta V(q)/2)
	//where C = ExpFourierBasis::basisCoefficient_
	void Galerkin::computeExpToTrigMat()
	{					
		for (int iOfFourier2=1; iOfFourier2 <=  2*(int)maxOfFourier_; iOfFourier2++)
		{
			trigToExpMat_(0, iOfFourier2) = expFourierCoeffs(iOfFourier2);
			trigToExpMat_(iOfFourier2, 0) = expFourierCoeffs(iOfFourier2);
		}
		trigToExpMat_(0, 0) = expFourierCoeffs(0) / sqrt(2.);
		
		for (int iOfFourier=1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
			for (int jOfFourier=1; jOfFourier <= (int) maxOfFourier_; jOfFourier++)
			{
				//cosine times cosine
				trigToExpMat_(2*iOfFourier, 2*jOfFourier) = 
					expFourierCoeffs(2*(iOfFourier+jOfFourier)) / sqrt(2.)
					+ expFourierCoeffs(2*abs(iOfFourier-jOfFourier)) / sqrt(2.);
					
				//sinus times sinus
				trigToExpMat_(2*iOfFourier-1, 2*jOfFourier-1) = 
					- expFourierCoeffs(2*(iOfFourier+jOfFourier)) / sqrt(2.)
					+ expFourierCoeffs(2*abs(iOfFourier-jOfFourier)) / sqrt(2.);
					
				int eps = (iOfFourier >= jOfFourier) - (iOfFourier <= jOfFourier); // -1 if i smaller, 0 if equal, 1 if i larger
				//cosine times sinus
				trigToExpMat_(2*iOfFourier, 2*jOfFourier-1) = 
					expFourierCoeffs(2*(iOfFourier+jOfFourier)-1) / sqrt(2.)				
					+ eps * expFourierCoeffs(2*abs(iOfFourier-jOfFourier)-1) / sqrt(2.);
					
				//sinus times cosine
				trigToExpMat_(2*iOfFourier-1, 2*jOfFourier) = 
					expFourierCoeffs(2*(iOfFourier+jOfFourier)-1) / sqrt(2.)
					- eps * expFourierCoeffs(2*abs(iOfFourier-jOfFourier)-1) / sqrt(2.);
			}				
			
					
		ofstream out_trigToExpMat("../output/Galerkin/trigToExpMat");
		displayMat(trigToExpMat_, out_trigToExpMat);
		
		expToTrigMat_  = inv(trigToExpMat_);
	}
	
	DMat Galerkin::convertToTrigBasis(const DMat& X)
	{
		return expToTrigMat_ * X;
	}
	
	void Galerkin::createQ()
	{
		Q_(1,0) = amplitude_ * beta_ / (2 * sqrt(2.));
		
		for (int iOfFourier=1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
		{
			Q_(2 * iOfFourier - 2, 2 * iOfFourier - 1) = amplitude_ * beta_ / 4;   //overwritten if i=1
			Q_(2 * iOfFourier    , 2 * iOfFourier - 1) = iOfFourier;
			if (iOfFourier != (int) maxOfFourier_)
				Q_(2 * iOfFourier + 2, 2 * iOfFourier - 1) = - amplitude_ * beta_ / 4;
		}
		Q_(0,1) = amplitude_ * beta_ / (2 * sqrt(2.));
		
		for (int iOfFourier=1; iOfFourier <= (int) maxOfFourier_; iOfFourier++)
		{
			if (iOfFourier != 1)
				Q_(2 * iOfFourier - 3, 2 * iOfFourier) = - amplitude_ * beta_ / 4;
			Q_(2 * iOfFourier - 1, 2 * iOfFourier) = - iOfFourier;
			if (iOfFourier != (int) maxOfFourier_)
				Q_(2 * iOfFourier + 1, 2 * iOfFourier) =   amplitude_ * beta_ / 4;
		}
		tQ_ = Q_.t();
	}
	
	void Galerkin::createP()
	{
		for (size_t iOfHermite=1; iOfHermite < nbOfHermite_; iOfHermite++)
			//P_(iOfHermite, iOfHermite+1) = cplx(sqrt(beta_) * sqrt(iOfHermite+1.), 0.);
			P_(iOfHermite-1, iOfHermite) = sqrt(beta_*iOfHermite);
		
		tP_ = P_.t();
	}
	
	void Galerkin::createLthm0()
	{
		cout << "createLthm0" << endl;
		for (int iOfHermite=1; iOfHermite < (int)nbOfHermite_; iOfHermite++)
			Lthm0_(iOfHermite, iOfHermite) = -beta_ * iOfHermite;
		cout << "end createLthm0" << endl;
	}
	
	using namespace Eigen;
	
	Galerkin::Galerkin(Input const& input):				//ex : [0:4]
		nbOfParticles_(input.nbOfParticles()),
		nbOfFourier_(input.nbOfFourier()),	//ex : 5
		nbOfHermite_(input.nbOfHermite()),
		maxOfFourier_((nbOfFourier_-1)/2),			//ex : 2
		sizeOfBasis_(pow(nbOfFourier_ * nbOfHermite_, nbOfParticles_)),
		SIdQ_(speye<SMat>(nbOfFourier_, nbOfFourier_)),
		SIdP_(speye<SMat>(nbOfHermite_, nbOfHermite_)),
		DIdQ_(eye<SMat>(nbOfFourier_, nbOfFourier_)),
		DIdP_(eye<SMat>(nbOfHermite_, nbOfHermite_)),
		Q_(nbOfFourier_, nbOfFourier_),
		P_(nbOfHermite_, nbOfHermite_),
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
		//expFourierCoeffs_(2 * nbOfFourier_, 0),
		trigToExpMat_(nbOfFourier_, nbOfFourier_, fill::zeros),
		expToTrigMat_(nbOfFourier_, nbOfFourier_, fill::zeros),
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
		
		cout << "Computing Q...";
		createQ();
		displayMat(Q_, "../output/Galerkin/Q");
		cout << "OK" << endl;
		
		cout << "############ P ############" << endl;
		createP();
		displayMat(P_, "../output/Galerkin/P");
		

		
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
	size_t Galerkin::iTens(size_t iOfFourier2, size_t iOfHermite) const
	{
		assert(iOfFourier2 < nbOfFourier_	&& iOfHermite < nbOfHermite_);
		return nbOfFourier_ * iOfHermite + iOfFourier2;
	}
	
	
	void Galerkin::compute()
	{
		cout << "start Galerkin::compute()" << endl;
		cout << "############ Leq ############" << endl;
		//cout << Leq_ << endl << endl;
		
		DMat DLeq = conv_to<DMat>::from(Leq_);
		displayMat(Leq_, "../output/Galerkin/Leq");
		
		//DMat DLeqSad = shapeSaddle(DLeq);
		//displayMat(DLeqSad, "../output/Galerkin/LeqSad");
		
		cout << "############ DLeqSad ############" << endl;
		//cout << DLeqSad << endl << endl;
		
		//DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);
		
		cout << "Computing LeqInv...";
		//DMat LeqInvSad = inv(DLeqSad);
		DMat LeqInv = invWithSaddle(DLeq);
		cout << "OK" << endl;
		
		//DMat LeqInv = LeqInvSad.submat(0, 0, sizeOfBasis_-1, sizeOfBasis_-1);
		
		cout << "############ LeqInv ############" << endl;
		displayMat(LeqInv, "../output/Galerkin/LeqInv");
		

		DVec H1Trig(sizeOfBasis_, fill::zeros);
		H1Trig(iTens(0, 1)) = 1;
		
		DMat H1TrigMat = reshape(H1Trig, nbOfFourier_, nbOfHermite_);
		
		//DMat H1Mat = trigToExpMat_ * H1TrigMat;
		//DVec H1 = reshape(H1Mat, nbOfFourier_ * nbOfHermite_, 1);
		
		DVec H1 = trigToExpTens_ * H1Trig;
		
		
		/*DVec H1(sizeOfBasis_);
		for (size_t iOfFourier2=0; iOfFourier2 <= 2*maxOfFourier_; iOfFourier2++)
			H1(iTens(iOfFourier2, 1)) = expFourierCoeffs_[iOfFourier2];*/
		
		cout << "############ H1Mat ############" << endl;
		displayMat(gettGiHj(0,1), "../output/Galerkin/H1");
		DMat H1Mat = reshape(gettGiHj(0,1), nbOfFourier_, nbOfHermite_);
		displayMat(H1Mat, "../output/Galerkin/H1Mat");
		
		displayMat(H1Trig, "../output/Galerkin/H1Trig");
		displayMat(H1TrigMat, "../output/Galerkin/H1TrigMat");
		
		
		displayMat(gettGiHj(0,0), "../output/Galerkin/G0Trig");
		
		DMat G0Mat = reshape(gettGiHj(0,0), nbOfFourier_, nbOfHermite_);
		displayMat(G0Mat, "../output/Galerkin/G0Mat");
		DMat LG0Mat = reshape(getLtGiHj(0,0), nbOfFourier_, nbOfHermite_);
		displayMat(LG0Mat, "../output/Galerkin/LG0Mat");
		
		displayMat(gettGiHj(1,2), "../output/Galerkin/G1H2Trig");
		displayMat(getLtGiHj(1,2), "../output/Galerkin/LG1H2Trig");
		
		displayMat(getLtGiHj(0,1), "../output/Galerkin/LH1");
		DVec LinvH1 = getLinvtGiHj(0,1);
		displayMat(LinvH1, "../output/Galerkin/LinvH1");
		
		//double lambda = LinvH1Sad(sizeOfBasis_);
		
		//cout << "############ lambda ############" << endl;
		//cout << lambda << endl << endl;
		
		DMat LinvH1Mat = reshape(LinvH1, nbOfFourier_, nbOfHermite_);
		
		cout << "############ LinvH1Mat ############" << endl;
		displayMat(LinvH1Mat, "../output/Galerkin/LinvH1Mat");
		
		DMat LinvH1MatTrig = convertToTrigBasis(LinvH1Mat);
		displayMat(LinvH1MatTrig, "../output/Galerkin/LinvH1MatTrig");
		
		DVec H1back = Leq_ * LinvH1;
		DMat H1backMat = reshape(H1back, nbOfFourier_, nbOfHermite_);
		displayMat(H1backMat, "../output/Galerkin/H1backMat");
		
		cx_vec eigvalLeq;
		cx_mat eigvecLeq;
		eig_gen(eigvalLeq, eigvecLeq, -DLeq);
		//eigs_gen(eigvalLeq, eigvecLeq, Leq_, 20);
		
		ofstream out_eigvalLeq("../output/Galerkin/eigvalLeq");
		//out_eigvalLeq << eigvalLeq << endl;
		displayCplx(eigvalLeq, out_eigvalLeq);
		
		double varOfH1 = -2 * dot(gettGiHj(0,1), LinvH1);
		cout << "varOfH1 = " << varOfH1 << endl;
		cout << "conductivity = " << .5 * varOfH1 << endl;
		
		DVec LinvL1LinvH1 = solveWithSaddle(Leq_, L1_ * LinvH1);
		//double varCoeff = -.5 * dot( L1_ * LinvH1, LeqInv * L1_ * LinvH1);
		double varCoeff = -2 * dot( L1_ * LinvH1, LinvL1LinvH1);
		cout << "varCoeff = " << varCoeff << endl;
	}
	
	DVec Galerkin::gettGiHj(int i, int j) const
	{
		DVec GiHjTrig(sizeOfBasis_, fill::zeros);
		GiHjTrig(iTens(i, j)) = 1;
		DVec GiHj = trigToExpTens_ * GiHjTrig;		
		return GiHj;
	}
	
	DVec Galerkin::gettGiHjTrig(int i, int j) const
	{
		DVec GiHjTrig(sizeOfBasis_, fill::zeros);
		GiHjTrig(iTens(i, j)) = 1;
		return GiHjTrig;
	}
	
	DVec Galerkin::getLtGiHj(int i, int j) const
	{
		return Leq_ * gettGiHj(i,j);
	}
	
	DVec Galerkin::getLtGiHjTrig(int i, int j) const
	{
		return expToTrigTens_ * getLtGiHj(i,j);
	}
	
	DVec Galerkin::getLinvtGiHj(int i, int j) const
	{
		return solveWithSaddle(Leq_, gettGiHj(i,j));
	}
	
	DVec Galerkin::getLinvtGiHjTrig(int i, int j) const
	{
		return expToTrigTens_ * solveWithSaddle(Leq_, gettGiHj(i,j));
	}
	
	SMat Galerkin::CVcoeffs() const
	{
		return conv_to<SMat>::from(getLinvtGiHj(0,1));
	}
	
	
	
	//#### LangevinGalerkin ####
	
	
	LangevinGalerkin::LangevinGalerkin(Input const& input):
		Galerkin(input)
	{
		cout << "############ Lham ############" << endl;
		Lham_ = kron(Q_, tP_) - kron(tQ_, P_);
		displayMat(Lham_, "../output/Galerkin/Lham");
		
		cout << "############ Lthm ############" << endl;
		createLthm();
		displayMat(Lthm_, "../output/Galerkin/Lthm");
		
		cout << "############ Leq ############" << endl;
		Leq_  = Lham_ + gamma_ * Lthm_;
		
		cout << "############ L1 ############" << endl;
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
		//SMat tempId = speye(PMat.n_rows, PMat.n_cols);
		SMat res = speye(1,1);
		for (int i = 0; i < nbOfParticles_; i++)
		{
			res = kron(res, SIdQ_);
			if (i==iOfParticleP) res = kron(res, PMat);
			else res = kron(res, SIdP_);
		}
		return res;
	}
	
	SMat BoundaryLangevinGalerkin::doubleMatToTens(SMat const& QMat, SMat const& PMat, int iOfParticleQ, int iOfParticleP)
	{
		/*assert(iOfVariableA < iOfVariableB);
		assert(iOfVariableB < nbOfVariables);
		assert(A.size() == B.size());*/
		//SMat tempId = speye(A.n_rows, A.n_cols);
		SMat res = speye(1,1);
		for (int i = 0; i < nbOfParticles_; i++)
		{
			if (i==iOfParticleQ) res = kron(res, QMat);
			else res = kron(res, SIdQ_);
			if (i==iOfParticleP) res = kron(res, PMat);
			else res = kron(res, SIdP_);
		}
		cout << "res : " << res.n_rows << " x " << res.n_cols << endl;
		return res;
	}
	
	/*SMat doubleMatToTens(SMat const& A, SMat const& B, int iOfVariableA, int iOfVariableB, int nbOfVariables)
	{
		assert(iOfVariableA < iOfVariableB);
		assert(iOfVariableB < nbOfVariables);
		assert(A.size() == B.size());
		SMat tempId = speye(A.n_rows, A.n_cols);
		SMat res = speye(1,1);
		for (int i = 0; i < iOfVariableA; i++)
			res = kron(res, tempId);
		res = kron(res, A);
		for (int i = iOfVariableA+1; i < iOfVariableB; i++)
			res = kron(res, tempId);
		res = kron(res, B);
		for (int i = iOfVariableB+1; i < nbOfVariables; i++)
			res = kron(res, tempId);
		cout << "res : " << res.n_rows << " x " << res.n_cols << endl;
		return res;
	}*/
	
		//psi = (1,1  1,2  ...  1,N_H  2,1 ... )
	//N_H blocks of size N_G (we concatene the columns of the matrix)
	//Allows to access elements involving a single particle !
	size_t BoundaryLangevinGalerkin::iTens(size_t iOfFourier2, size_t iOfHermite, int iOfParticle) const
	{
		assert(iOfFourier2 < nbOfFourier_	&& iOfHermite < nbOfHermite_ && iOfParticle < nbOfParticles_);
		return pow(nbOfFourier_*nbOfHermite_ ,iOfParticle) * ( nbOfFourier_ * iOfHermite + iOfFourier2);
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
		Lham_ = arma::zeros(sizeOfBasis_, sizeOfBasis_);
		for (int i=0; i<nbOfParticles_; i++)
		{
			Lham_ += doubleMatToTens(Q_, tP_, i, i);
			if (i>0) Lham_ -= doubleMatToTens(Q_, tP_, i, i-1);
			Lham_ -= doubleMatToTens(tQ_, P_, i, i);
			if (i<nbOfParticles_-1) Lham_ += doubleMatToTens(Q_, tP_, i+1, i);
		}
		Lham_ /= beta_;
		displayMat(Lham_, "../output/Galerkin/Lham");
	}
	
	void BoundaryLangevinGalerkin::createLthm()
	{
		createLthm0();
		//SMat Lthm1_ = kron(SIdQ_, Lthm0_);
		Lthm_ = arma::zeros(sizeOfBasis_, sizeOfBasis_);
		for (int i=0; i<nbOfParticles_; i++)
		{
			//cout << i << " < " << nbOfParticles_ << endl;
			Lthm_ += PMatToTens(Lthm0_, i);
		}
		Lthm_ /= beta_;
		displayMat(Lthm_, "../output/Galerkin/Lthm");
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
		//cout << Leq_ << endl << endl;
		
		DMat DLeq = conv_to<DMat>::from(Leq_);
		displayMat(Leq_, "../output/Galerkin/Leq");
		
		//DMat DLeqSad = shapeSaddle(DLeq);
		//displayMat(DLeqSad, "../output/Galerkin/LeqSad");
		
		cout << "############ DLeqSad ############" << endl;
		//cout << DLeqSad << endl << endl;
		
		//DMat IdSad(sizeOfBasis_+1, sizeOfBasis_, fill::eye);
		
		cout << "Computing LeqInv...";
		//DMat LeqInvSad = inv(DLeqSad);
		DMat LeqInv = invWithSaddle(DLeq);
		cout << "OK" << endl;
		displayMat(LeqInv, "../output/Galerkin/LeqInv");
		
		DVec N0H2Trig(sizeOfBasis_, fill::zeros);
		N0H2Trig(iTens(0, 2, 0)) = 1;
		
		DVec N0H2 = trigToExpTens_ * N0H2Trig;
		
		DMat N0H2Mat = reshape(N0H2, pow(nbOfFourier_, nbOfParticles_), pow(nbOfHermite_, nbOfParticles_));
		
		//SparseMatrix<double> Atest;
		cx_vec eigvalLeq;
		cx_mat eigvecLeq;
		eig_gen(eigvalLeq, eigvecLeq, -DLeq);
		//eigs_gen(eigvalLeq, eigvecLeq, Leq_, 3);
		
		ofstream out_eigvalLeq("../output/Galerkin/eigvalLeq");
		//out_eigvalLeq << eigvalLeq << endl;
		displayCplx(eigvalLeq, out_eigvalLeq);
		
		/*double varOfH1 = -2 * dot(gettGiHj(0,1), LinvH1);
		cout << "varOfH1 = " << varOfH1 << endl;
		cout << "conductivity = " << .5 * varOfH1 << endl;
		
		DVec LinvL1LinvH1 = solveWithSaddle(Leq_, L1_ * LinvH1);
		//double varCoeff = -.5 * dot( L1_ * LinvH1, LeqInv * L1_ * LinvH1);
		double varCoeff = -2 * dot( L1_ * LinvH1, LinvL1LinvH1);
		cout << "varCoeff = " << varCoeff << endl;*/
	}
	

	
}



