#ifndef SPARSEMATRIX_IMPL_HPP
#define SPARSEMATRIX_IMPL_HPP

template<class ScalarType, template<class> class WrappingPolicy>
inline SparseMatrix<ScalarType,WrappingPolicy>::SparseMatrix(size_t numberOfRows, size_t numberOfColumns)
:wrapped_(numberOfRows,numberOfColumns)
{}


#endif
