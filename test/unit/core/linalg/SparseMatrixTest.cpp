#include "gtest/gtest.h"

#include "simol/core/linalg/SparseMatrix.hpp"

namespace simol
{

  TEST(SparseMatrix, IteratorIsCorrect)
  {
    SparseMatrix<double> B(2, 3);
    B.insert(0, 1) = 24;
    B.insert(1, 2) = 36;
    SparseMatrix<double>::iterator it(B, 1);
    std::cout << it.row() << "\t";
    std::cout << it.col() << "\t";
    std::cout << it.value() << std::endl;
    SparseMatrix<double> C = B*(2.0); 
    SparseMatrix<double> D = 2.0*B; 
    std::cout << B(0,1) << std::endl;
    std::cout << C(0,1) << std::endl;
    B += C; 
    std::cout << B(0,1) << std::endl; 
    EXPECT_EQ(0, it.row());
    EXPECT_EQ(72, B(0,1));
    EXPECT_EQ(48, D(0,1));
  }

}
