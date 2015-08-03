#include "EigenSolverTest.hpp"

#include <cppunit/extensions/HelperMacros.h>


namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION(EigenSolverTest);


  void EigenSolverTest::setUp()
  {
    identity_ = new SparseMatrix<double>(MatrixMarketFile("../src/test/data/identity5.mtx"));
  }

  void EigenSolverTest::tearDown()
  {
    delete identity_;
  }

  void EigenSolverTest::testEigenvaluesOfIdentityAreAllOne()
  {
    EigenSolver<double> arpack;
    double * eigenvalues = arpack.solve(*identity_, 2, 1e-6);
    CPPUNIT_ASSERT_EQUAL(1.0, eigenvalues[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, eigenvalues[1]);
  }

}



