#include "gtest/gtest.h"

#include "simol/core/linalg/DenseMatrix.hpp"

namespace simol
{

  TEST(MatrixTest, NumberOfColumnsIsCorrect)
  {
    DenseMatrix<double, eigen> mat(3, 5);
    size_t expected = 5;
    EXPECT_EQ(expected, mat.numberOfColumns());
  }

  TEST(MatrixTest, NumberOfRowsIsCorrect)
  {
    DenseMatrix<double, eigen> mat(3, 5);
    size_t expected = 3;
    EXPECT_EQ(expected, mat.numberOfRows());
  }
}
