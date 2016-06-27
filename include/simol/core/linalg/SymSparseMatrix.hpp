#ifndef SIMOL_SYMSPARSEMATRIX_HPP
#define SIMOL_SYMSPARSEMATRIX_HPP

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
  class SymSparseMatrix;
}

#include "simol/core/linalg/SymSparseMatrix_eigen.hpp"


#endif
