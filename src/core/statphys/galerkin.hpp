#pragma once
//#include "SparseMatrix.hpp"
//#include <Eigen/Dense>

#include "tools.hpp"
#include "input.hpp"

//#include "eigen.hpp"

namespace simol
{
	SpMat tensor(SpMat A, SpMat B);
	
	class Galerkin
	{
		//DenseMatrix<double> A;
		size_t numberOfFourier_, numberOfHermite_;
		int maxOfFourier_;
		SpMat Q_, P_;
		SpMat Lthm0_, Lthm_, Lham_;
		SpMat Leq_;
		double beta_;
	public:
		Galerkin(Input const& input);
		size_t iPosi(int iOfFourier, size_t iOfHermite) const;
		size_t iSign(int iOfFourier, size_t iOfHermite) const;
		size_t iPosiF(int iOfFourier) const;
		size_t iSignF(int iOfFourier) const;
		void solve();
	};
}