#ifndef SIMOL_SPARSEMATRIXTEST_HPP
#define SIMOL_SPARSEMATRIXTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "../../core/linalg/SparseMatrix.hpp"

namespace simol
{
  class SparseMatrixTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(SparseMatrixTest);
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp();
    void tearDown();
  };
}
#endif 

