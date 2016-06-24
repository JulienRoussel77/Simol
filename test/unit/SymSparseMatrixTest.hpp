#ifndef SIMOL_SYMSPARSEMATRIXTEST_HPP
#define SIMOL_SYMSPARSEMATRIXTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "../../core/linalg/SymSparseMatrix.hpp"

namespace simol
{
  class SymSparseMatrixTest : public CppUnit::TestFixture
  {
      CPPUNIT_TEST_SUITE(SymSparseMatrixTest);
      CPPUNIT_TEST_SUITE_END();

    public:
      void setUp();
      void tearDown();
  };
}
#endif

