#include "DenseMatrixTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION(DenseMatrixTest);

  void DenseMatrixTest::setUp()
  {
    randomNumberOfRows = 12;
    randomNumberOfColumns = 7;
    randomMatrix = new DenseMatrix<double>(randomNumberOfRows, randomNumberOfColumns);
  }

  void DenseMatrixTest::tearDown()
  {
    delete randomMatrix;
  }

  void DenseMatrixTest::test_number_of_rows_is_correct()
  {
    CPPUNIT_ASSERT_EQUAL(randomMatrix->number_of_rows(), randomNumberOfRows);
  }

  void DenseMatrixTest::test_number_of_columns_is_correct()
  {
    CPPUNIT_ASSERT_EQUAL(randomMatrix->number_of_columns(), randomNumberOfColumns);
  }


}



