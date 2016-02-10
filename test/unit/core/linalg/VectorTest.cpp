#include "gtest/gtest.h"
#include "Vector.hpp"

namespace simol
{

  TEST(VectorTest, size_is_correct)
  {
    size_t expected_size = 10;
    Vector<double> vector(expected_size);
    size_t actual_size = vector.size();
    EXPECT_EQ(expected_size, actual_size);
  }

  TEST(VectorTest, min_index_is_correct)
  {
      Vector<double> vec(3);
      vec(0) = 0;
      vec(1) = -4;
      vec(2) = 6;

      size_t expected = 1;
      size_t actual = vec.min_index();

      EXPECT_EQ(expected, actual);

  }

  TEST(VectorTest, sort_is_correct)
  {
      size_t size = 5;
      Vector<double> vec(size);
      vec(0) = 4;
      vec(1) = 1;
      vec(2) = 2;
      vec(3) = 0;
      vec(4) = 3;

      Vector<double> actual = vec.sort();

      Vector<double> expected(size);
      for(size_t index=0; index<size; ++index)
          expected(index) = index;

      // penser à utiliser google mock a la place
      for(size_t index=0; index<size; ++index)
          EXPECT_EQ( actual(index), expected(index));

  }

}

