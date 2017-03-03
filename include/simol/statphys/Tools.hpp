#ifndef SIMOL_TOOLS_HPP
#define SIMOL_TOOLS_HPP

#include <iostream>
#include <fstream>
#include <locale>   //contains std::tolower
#include <vector>
//#include "simol/core/linalg/SparseMatrix.hpp"
//#include "simol/core/linalg/DenseMatrix.hpp"
//#include "simol/core/linalg/Vector.hpp"
#include <stdexcept>
#include <list>
#include <iomanip>
#include <fstream>
#include <vector>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/IterativeSolvers>

using std::setw;
using std::setprecision;
using std::ofstream;
using std::ostream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::to_string;
using std::ifstream;
using std::tolower;
using std::shared_ptr;
using std::make_shared;
using std::runtime_error;

using std::vector;
using std::list;
using std::map;

typedef std::complex<double> cplx;

//--- if using the Simol wrappers ---
//typedef simol::SparseMatrix<double> SMat;
//typedef simol::DenseMatrix<double> DMat;
//typedef simol::Vector<double> DVec;

//--- standard eigen functions ---
typedef Eigen::SparseMatrix<double> SMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DVec;

typedef Eigen::Matrix<cplx, Eigen::Dynamic, 1> DVecCplx;

typedef Eigen::SparseMatrix<long int> SMatInt;
typedef Eigen::Matrix<long int, Eigen::Dynamic, Eigen::Dynamic> DMatInt;
typedef Eigen::Matrix<long int, Eigen::Dynamic, 1> DVecInt;

typedef Eigen::Triplet<double> Trid;

double modulo(double variable, double mini, double maxi);
int intModulo(int variable, int maxi);
void displayTime(double time);
int getNbOfLines(ifstream& file);
//bool hasSmallerNorm(const cplx& a, const cplx& b);
bool hasSmallerNorm(cplx a, cplx b);

double dot(DVec const& u, DVec const& v);
DMat reshape(DVec const& u, int nbOfRows, int nbOfCols);
DVec rint(DVec const& u);
DVec polynomialProduct(DVec const& P, DVec const& Q);

const int idKineticEnergy = 0;
const int idPotentialEnergy = 1;
const int idPressure = 2;
const int idInternalEnergy = 3;
const int idInternalTemperature = 4;
const int idVelocity = 5;
const int idForce = 6;
const int idLength = 7;
const int idMidFlow = 8;
const int idSumFlow = 9;
const int idModiFlow = 10;
const int nbOfIdObs = 11;

#endif
