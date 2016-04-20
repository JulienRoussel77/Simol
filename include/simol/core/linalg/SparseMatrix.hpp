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

  template<typename Scalar, typename Wrapper = eigen>
  class SparseMatrix;

  template<typename Scalar, typename Wrapper>
  std::ifstream & operator>>(std::ifstream & fileToRead,
                             SparseMatrix<Scalar, Wrapper> & matrixToWrite);


}

#include "simol/core/linalg/SparseMatrix_eigen.hpp"


#endif
