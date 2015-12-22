#ifndef SIMOL_MATRIXMARKETFILETEST_HPP
#define SIMOL_MATRIXMARKETFILETEST_HPP

#include <cppunit/extensions/HelperMacros.h>

#include "../../core/io/MatrixMarketFile.hpp"

namespace simol
{
  class MatrixMarketFileTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(MatrixMarketFileTest);
      CPPUNIT_TEST(testNumberOfRows);
      CPPUNIT_TEST(testNonExistentFilenameThrowsException);
    CPPUNIT_TEST_SUITE_END();

  public:
    void setUp();
    void tearDown();
  public:
    void testNonExistentFilenameThrowsException();
    void testNumberOfRows();
  private:
    MatrixMarketFile * kinetic_;
  };
}
#endif 

