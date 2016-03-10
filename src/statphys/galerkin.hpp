#pragma once
//#include "SparseMatrix.hpp"
//#include <Eigen/Dense>

#include "tools.hpp"
#include "input.hpp"

//#include "eigen.hpp"

namespace simol
{
	SMat tensor(SMat A, SMat B);

	class Galerkin
	{
		//DenseMatrix<double> A;
		size_t numberOfFourier_, numberOfHermite_, maxOfFourier_;
		size_t sizeOfBasis_;
		SMat Q_, P_;
		SMat tQ_, tP_;
		SMat Lthm0_, Lthm_, Lham_;
		SMat Leq_;
		double beta_;
		double amplitude_;
		size_t nbIntegrationNodes_;
        std::vector<double> fourierCoeffsExp_;
	public:
		Galerkin(Input const& input);
		double potential(double q);
		size_t iTens(size_t iOfFourier2, size_t iOfHermite) const;
		DMat shapeSaddle(DMat& A);
		DMat unshapeSaddle(DMat& Asad);
		DVec shapeSaddle(DVec& X);
		DVec unshapeSaddle(DVec& Xsad);
		void computeFourierCoeffsExp();
		void compute();
	};
}
