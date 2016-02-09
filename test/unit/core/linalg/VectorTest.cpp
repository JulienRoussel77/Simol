#include "gtest/gtest.h"
#include "Vector.hpp"

namespace simol
{

  TEST(VectorTest, SizeIsCorrect)
  {
    size_t expected_size = 10;
    Vector<double> vector(expected_size);
    size_t actual_size = vector.size();
    EXPECT_EQ(expected_size, actual_size);
  }

}

