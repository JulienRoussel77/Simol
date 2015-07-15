#include "SparseMatrixTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION(SparseMatrixTest);


  void SparseMatrixTest::setUp()
  {
    identity_ = new SparseMatrix<double>(MatrixMarketFile("../src/test/data/identity3.mtx"));
  }

  void SparseMatrixTest::tearDown()
  {
    delete identity_;
  }

  void SparseMatrixTest::testEigenvaluesOfIdentityAreAllOne()
  {
    CPPUNIT_ASSERT_EQUAL(10, 1);
  }

}



