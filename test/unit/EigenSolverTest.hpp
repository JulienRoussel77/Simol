#ifndef SIMOL_EIGENSOLVERTEST_HPP
#define SIMOL_EIGENSOLVERTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "../../core/linalg/EigenSolver.hpp"

namespace simol
{
  class EigenSolverTest : public CppUnit::TestFixture
  {
      CPPUNIT_TEST_SUITE(EigenSolverTest);
      CPPUNIT_TEST(testEigenvaluesOfIdentityAreAllOne);
      CPPUNIT_TEST_SUITE_END();

    public:
      void setUp();
      void tearDown();
    public:
      void testEigenvaluesOfIdentityAreAllOne();
    private:
      SparseMatrix<double> * identity_;
  };
}
#endif

