#ifndef SIMOL_TOOLS_HPP
#define SIMOL_TOOLS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <Eigen/Sparse>

#include <stdexcept>

using std::ofstream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::to_string;
using std::vector;
using std::ifstream;

using std::vector;

typedef std::complex<double> cplx;

typedef arma::sp_mat SMat;
typedef arma::mat DMat;
typedef arma::vec DVec;
typedef arma::rowvec DRow;

using arma::fill::zeros;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Matrix<cplx, Dynamic, Dynamic> DsMat;
typedef Eigen::Triplet<double> T;

double modulo(double variable, double mini, double maxi);

#endif
