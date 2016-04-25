#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "simol/core/linalg/Vector.hpp"
#include "simol/core/io/MatrixMarketFile.hpp"
#include "simol/core/linalg/eigen.hpp"

#include "simol/core/linalg/dsaupd.hpp"
#include "simol/core/linalg/dseupd.hpp"

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

#include "simol/core/linalg/SparseMatrix_eigen.hpp"


#endif
