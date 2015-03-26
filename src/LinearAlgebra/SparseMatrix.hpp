#ifndef SIMOL_SPARSEMATRIX_HPP
#define SIMOL_SPARSEMATRIX_HPP

#include "SparseMatrix_fwd.hpp"

namespace simol
{
  //=====================
  // FORWARD DECLARATIONS
  //=====================

  template<class ScalarType, template<class> class Library>
  std::ifstream & operator>>(std::ifstream & fileToRead, SparseMatrix<ScalarType,Library> & matrixToWrite);

  //====================
  // SPARSE MATRIX CLASS
  //====================

  template<class ScalarType>
  class SparseMatrix<ScalarType,eigen>
  {

    //=================
    // FRIEND FUNCTIONS
    //=================

    friend std::ifstream & operator>> <>(std::ifstream & fileToRead, SparseMatrix<ScalarType,eigen> & matrixToWrite);

    //=============
    // CONSTRUCTORS
    //=============

    public:

SparseMatrix(size_t const numberOfRows, size_t const numberOfColumns);
    
    //=============
    // DATA MEMBERS
    //=============

    private:

      typename eigen<ScalarType>::SparseMatrixType wrapped_;
  };
}

#include "SparseMatrix_impl.hpp"

#endif
