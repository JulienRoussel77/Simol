#ifndef SIMOL_DENSEMATRIXTEST_HPP
#define SIMOL_DENSEMATRIXTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "../../core/linalg/Matrix.hpp"

namespace simol
{
  class DenseMatrixTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(DenseMatrixTest);
      CPPUNIT_TEST(test_number_of_rows_is_correct);
      CPPUNIT_TEST(test_number_of_columns_is_correct);
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp();
    void tearDown();

    void test_number_of_rows_is_correct();
    void test_number_of_columns_is_correct();

  private:
    size_t randomNumberOfRows;
    size_t randomNumberOfColumns;
    DenseMatrix<double>* randomMatrix;
  };
}
#endif 

