#include "SparseMatrixTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

#include "EigenSolver.hpp"

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
    EigenSolver<double> arpack;
    double * eigenvalues = arpack.solve(*identity_, 2, 1e-6);
    CPPUNIT_ASSERT_EQUAL(1.0, eigenvalues[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, eigenvalues[1]);
  }

}



