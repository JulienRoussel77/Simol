#ifndef SIMOL_MATRIXMARKETFILE_HPP
#define SIMOL_MATRIXMARKETFILE_HPP

#include <cstddef>
#include <cstdio>

extern "C"
{
  #include "mmio.h"
}


namespace simol
{
  class MatrixMarketFile
  {
    public:
      MatrixMarketFile(char * filename);
      ~MatrixMarketFile();
    public:
      std::size_t numberOfRows() const;
      std::size_t numberOfColumns() const;
      std::size_t numberOfNonzeros() const;
    private:
      FILE * file_;
      int numberOfRows_;
      int numberOfColumns_;
      int numberOfNonzeros_;
  };

  inline
  MatrixMarketFile::MatrixMarketFile(char * filename)
  : file_(fopen(filename, "r"))
  {
    MM_typecode matcode;
    mm_read_banner(file_, &matcode);
    mm_read_mtx_crd_size(file_, 
                         &numberOfRows_,
                         &numberOfColumns_,
                         &numberOfNonzeros_);
  }

  inline
  MatrixMarketFile::~MatrixMarketFile()
  { fclose(file_); }

  inline std::size_t 
  MatrixMarketFile::numberOfRows() const
  { return numberOfRows_; }

  inline std::size_t 
  MatrixMarketFile::numberOfColumns() const
  { return numberOfColumns_; }
  
  inline std::size_t 
  MatrixMarketFile::numberOfNonzeros() const
  { return numberOfNonzeros_; }
}


#endif
