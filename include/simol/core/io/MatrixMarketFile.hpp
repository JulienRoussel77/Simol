#ifndef SIMOL_MATRIXMARKETFILE_HPP
#define SIMOL_MATRIXMARKETFILE_HPP

#include <cstddef>
#include <cstdio>
#include <ios>
#include <iostream>
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
    if(!content_)
      throw std::ios_base::failure("Unable to open " + filename);
    MM_typecode matcode;
    if (mm_read_banner(content_, &matcode) != 0)
      throw new std::ios_base::failure("Unable to read the banner");
    int ret_code = mm_read_mtx_crd_size(content_,
                                        &numberOfRows_,
                                        &numberOfColumns_,
                                        &numberOfNonzeros_);
    if( ret_code != 0 )
      throw new std::ios_base::failure("Unable to read the data");

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
