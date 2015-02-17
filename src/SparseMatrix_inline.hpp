#ifndef SPARSEMATRIX_INLINE_HPP
#define SPARSEMATRIX_INLINE_HPP

template<class ScalarType>
inline SparseMatrix<ScalarType>::SparseMatrix(size_t numberOfRows, size_t numberOfColumns)
:wrapped_(numberOfRows,numberOfColumns)
{}


#endif
