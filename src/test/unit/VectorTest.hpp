#ifndef SIMOL_VECTORTEST_HPP
#define SIMOL_VECTORTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <ctime>

#include "Vector.hpp"

namespace simol
{

  class VectorTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(VectorTest);
      CPPUNIT_TEST(test_size);
      CPPUNIT_TEST(test_sum);
    CPPUNIT_TEST_SUITE_END();

    public:
      void test_size();
      void test_sum();
 
  };
}


#endif
