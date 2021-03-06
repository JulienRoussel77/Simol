#ifndef SIMOL_TOOLS_HPP
#define SIMOL_TOOLS_HPP

// Disables certain warnings when including the following files
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include <iostream>
#include <fstream>
#include <locale>   //contains std::tolower
#include <vector>
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

#pragma GCC diagnostic pop

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


//--- standard eigen functions ---
typedef Eigen::SparseMatrix<double> SMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> DVec;

typedef Eigen::Matrix<cplx, Eigen::Dynamic, 1> DVecCplx;

typedef Eigen::SparseMatrix<long int> SMatInt;
typedef Eigen::Matrix<long int, Eigen::Dynamic, Eigen::Dynamic> DMatInt;
typedef Eigen::Matrix<long int, Eigen::Dynamic, 1> DVecInt;

typedef Eigen::Triplet<double> Trid;

using Eigen::JacobiSVD;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;

using Eigen::ComplexEigenSolver;
using Eigen::real;
using Eigen::imag;

using Eigen::FullPivLU;
using Eigen::HouseholderQR;

double modulo(double variable, double mini, double maxi);
int intModulo(int variable, int maxi);
void displayTime(double time);

double dot(DVec const& u, DVec const& v);
DVec rint(DVec const& u);
DVec polynomialDerivative(DVec const& P);

// integer identifiers of the different observables, 
const int idKineticEnergy = 0;
const int idPotentialEnergy = 1;
const int idTotalEnergy = 2;
const int idPressure =3;
const int idInternalEnergy = 4;
const int idInternalTemperature = 5;
const int idVelocity = 6;
const int idForce = 7;
const int idLength = 8;
const int idMidFlux = 9;
const int idSumFlux = 10;
const int idModiFlux = 11;
const int idLagrangeMultiplier = 12;
const int nbOfIdObs = 13;

#endif
