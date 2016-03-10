#ifndef SIMOL_TOOLS_HPP
#define SIMOL_TOOL_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Sparse>
#include<armadillo>

typedef std::complex<double> cplx;

typedef arma::sp_mat SMat;
typedef arma::mat DMat;
typedef arma::vec DVec;
typedef arma::rowvec DRow;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Matrix<cplx, Dynamic, Dynamic> DsMat;
typedef Eigen::Triplet<double> T;

#endif
