/*
 * File:   TensorFile.hpp
 * Author: cdoucet
 *
 * Created on 28 octobre 2015, 15:10
 */

#ifndef TENSORFILE_HPP
#define TENSORFILE_HPP

#include <cassert>

namespace simol
{
  class TensorFile
  {
    public:

      explicit
      TensorFile(std::string const & filename)
        : content_(filename, std::ios::in)
      {
        read_size();
        read_nonzeros();
      }

      ~TensorFile()
      { content_.close(); }


    private:

      void
      read_size()
      { content_ >> size_; }

      void
      read_nonzeros()
      {
        std::size_t i, j, k, l;
        double nonzero;
        while(content_ >> i >> j >> k >> l >> nonzero)
          nonzeros_.push_back(NonZero(vector_index(i, k), vector_index(j, l), nonzero));
      }

      std::size_t
      vector_index(std::size_t const rowIndex,
                   std::size_t const columnIndex)
      {
        assert(rowIndex < size_ && columnIndex < size_);
        return rowIndex * size_ + columnIndex;
      }


    private:
      typedef Eigen::Triplet<double, size_t> NonZero;
      std::size_t size_;
      std::vector<NonZero> nonzeros_;
      std::ifstream content_;
  };
}

#endif  /* TENSORFILE_HPP */

