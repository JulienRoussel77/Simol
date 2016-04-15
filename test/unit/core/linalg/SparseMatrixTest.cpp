#include "gtest/gtest.h"

#include "simol/core/linalg/SparseMatrix.hpp"

namespace simol
{

    TEST(SparseMatrix, IteratorIsCorrect)
    {
        SparseMatrix<double> B(2, 3);
        B.insert(0,1) = 24;
        B.insert(1,2) = 36;
        SparseMatrix<double>::iterator it(B, 1);
        std::cout << it.row() <<"\t";
        std::cout << it.col() << "\t";
        std::cout << it.value() << std::endl;
        EXPECT_EQ(0, it.row());
    }

}
