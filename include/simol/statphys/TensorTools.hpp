#ifndef SIMOL_TENSORTOOLS_HPP
#define SIMOL_TENSORTOOLS_HPP

#include "Tools.hpp"


SMat kron(const SMat& A, const SMat& B);
DMat kron(const DMat& A, const DMat& B);
SMat kron(const SMat& A, const DMat& B);
SMat kron(const DMat& A, const SMat& B);
double computeSpectralGap(SMat const& A);
void displayCplx(const DVecCplx& X, ostream& out);
void display(const DVec& A, string path);
void display(const DMat& A, ostream& out);
void display(const DMat& A, string path);
void display(const SMat& A, ostream& out);
void display(const SMat& A, string path);

#endif
