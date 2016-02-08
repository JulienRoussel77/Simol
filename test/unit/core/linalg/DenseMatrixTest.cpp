#include "gtest/gtest.h"

#include "DenseMatrix.hpp"

namespace simol
{

    TEST(MatrixTest, NumberOfRowsIsCorrect)
    {
        DenseMatrix<double, eigen> mat(3,5);
        size_t expected = 3;
        EXPECT_EQ(expected, mat.numberOfRows());
    }
}
