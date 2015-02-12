#ifndef SIMOL_VECTORTEST_HPP
#define SIMOL_VECTORTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <ctime>

#include "vector.hpp"

class VectorTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(VectorTest);
  CPPUNIT_TEST_SUITE_END();

public:

  void setUp();
  void tearDown();

private:
 
};


#endif 
