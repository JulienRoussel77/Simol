#include "gtest/gtest.h"

#include "simol/core/linalg/SymSparseMatrix.hpp"
#include "simol/core/linalg/DenseMatrix.hpp"

namespace simol
{

  TEST(SymSparseMatrix, DenseConversionIsCorrect)
  {


    SymSparseMatrix<double> A(3, 3);
    A.insert(0, 1) = 24;
    A.insert(1, 2) = 36;
    A.insert(1, 1) = 12;

    DenseMatrix<double> B = A.dense();
    
    std::size_t si  = A.numberOfRows(); 
    std::size_t si2 = 3; 

    EXPECT_EQ(24, B(1,0));
    EXPECT_EQ(36, B(2,1));
    EXPECT_EQ(12, B(1,1));
    EXPECT_EQ(0, B(0,0));
    EXPECT_EQ(si2,si);

  }

}
