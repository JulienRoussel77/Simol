#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "core/linalg/Vector.hpp"
#include "core/io/MatrixMarketFile.hpp"
#include "eigen.hpp"

#include "dsaupd.hpp"
#include "dseupd.hpp"

#include <fstream>
#include <string>

namespace simol
{

  template<class ScalarType, template<class> class WrappedLibrary = eigen>
  class SparseMatrix;

  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead,
                             SparseMatrix<ScalarType, Library> & matrixToWrite);


}

#include "SparseMatrix_eigen.hpp"


#endif
