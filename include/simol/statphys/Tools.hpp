#ifndef SIMOL_TOOLS_HPP
#define SIMOL_TOOLS_HPP

#include <iostream>
#include <fstream>
#include <locale>   //contains std::tolower
#include <vector>
#include "simol/core/linalg/SparseMatrix.hpp"
#include "simol/core/linalg/DenseMatrix.hpp"
#include "simol/core/linalg/Vector.hpp"
#include <stdexcept>
#include <list>

using std::ofstream;
using std::cout;
using std::endl;
using std::min;
using std::max;
using std::string;
using std::to_string;
using std::ifstream;
using std::tolower;

using std::vector;
using std::list;

typedef std::complex<double> cplx;

typedef simol::SparseMatrix<double> SMat;
typedef simol::DenseMatrix<double> DMat;
typedef simol::Vector<double> DVec;

double modulo(double variable, double mini, double maxi);
int intModulo(int variable, int maxi);
void displayTime(double time);
int getNbOfLines(ifstream& file);

#endif
