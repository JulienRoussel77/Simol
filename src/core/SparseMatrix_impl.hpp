#ifndef SIMOL_SPARSEMATRIX_IMPL_HPP
#define SIMOL_SPARSEMATRIX_IMPL_HPP

#include <cstdlib>


namespace simol
{
  template<class ScalarType> inline
  SparseMatrix<ScalarType,eigen>::SparseMatrix(size_t const numberOfRows,
                                               size_t const numberOfColumns)
  :wrapped_(numberOfRows,numberOfColumns)
  {}
}

#endif
