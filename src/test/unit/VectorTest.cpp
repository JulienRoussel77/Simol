
#include "VectorTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION( VectorTest );

  void VectorTest::test_size()
  {
    randomSize = 10;
    Vector<double> vector(randomSize);
    CPPUNIT_ASSERT_EQUAL(randomSize, vector.size());
  }

}

