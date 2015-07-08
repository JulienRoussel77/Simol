#ifndef SIMOL_MATRIXMARKETFILE_HPP
#define SIMOL_MATRIXMARKETFILE_HPP

#include <cstddef>
#include <cstdio>
#include <string>

extern "C"
{
  #include "mmio.h"
}


namespace simol
{
  class MatrixMarketFile
  {
    public:
      MatrixMarketFile(std::string const & filename);
      ~MatrixMarketFile();
    public:
      FILE * content() const;
      std::size_t numberOfRows() const;
      std::size_t numberOfColumns() const;
      std::size_t numberOfNonzeros() const;
    private:
      FILE * content_;
      int numberOfRows_;
      int numberOfColumns_;
      int numberOfNonzeros_;
  };

  inline
  MatrixMarketFile::MatrixMarketFile(std::string const & filename)
  : content_(fopen(filename.c_str(), "r"))
  {
    MM_typecode matcode;
    mm_read_banner(content_, &matcode);
    mm_read_mtx_crd_size(content_, 
                         &numberOfRows_,
                         &numberOfColumns_,
                         &numberOfNonzeros_);
  }

  inline FILE *
  MatrixMarketFile::content() const
  { return content_; }

  inline
  MatrixMarketFile::~MatrixMarketFile()
  { fclose(content_); }

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
