#pragma once

//#define UNUSED(expr) (void)(expr)

#include <iostream>
#include <fstream>
#include <vector>

using std::cout;
using std::endl;
using std::min;
using std::max;

using std::vector;

#include <Eigen/Sparse>
typedef std::complex<double> cplx;
typedef Eigen::SparseMatrix<cplx> SpMat;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Matrix<cplx, Dynamic, Dynamic> DsMat;
typedef Eigen::Triplet<double> T;