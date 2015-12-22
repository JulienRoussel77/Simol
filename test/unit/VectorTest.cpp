
#include "VectorTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION( VectorTest );

  void VectorTest::test_size()
  {
    size_t expected_size = 10;
    Vector<double> vector(expected_size);
    size_t actual_size = vector.size();
    CPPUNIT_ASSERT_EQUAL(expected_size, actual_size);
  }

  void VectorTest::test_sum()
  {
    Vector<double> left(3);
    Vector<double> right(3);
    Vector<double> sum(3);

    for (size_t i=0;i<3;++i)
    {
      left(i) = 1;
      right(i) = 2;
      sum(i) = 3;
    }
    //CPPUNIT_ASSERT_EQUAL(sum, left+right);
  }

}

