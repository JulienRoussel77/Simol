#pragma once

//#define UNUSED(expr) (void)(expr)

#include <iostream>
#include <fstream>
#include <vector>

using std::ofstream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::vector;
using std::ifstream;

using std::vector;

#include <Eigen/Sparse>
typedef std::complex<double> cplx;

#include<armadillo>
//typedef Eigen::SparseMatrix<cplx> SpMat;
typedef arma::sp_mat SMat;
typedef arma::mat DMat;
typedef arma::vec DVec;
typedef arma::rowvec DRow;

//typedef arma::cx_mat DCxMat;
//typedef arma::cx_vec DCxVec;

using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Matrix<cplx, Dynamic, Dynamic> DsMat;
typedef Eigen::Triplet<double> T;