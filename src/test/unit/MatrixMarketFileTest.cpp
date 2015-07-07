#include "MatrixMarketFileTest.hpp"

#include <cppunit/extensions/HelperMacros.h>

CPPUNIT_TEST_SUITE_REGISTRATION(MatrixMarketFileTest);


void MatrixMarketFileTest::setUp()
{
  kinetic_ = new MatrixMarketFile("../../../data/quantum_chemistry/kinetic_matrix.mtx");
}

void MatrixMarketFileTest::tearDown()
{
  delete kinetic_;
}

void MatrixMarketFileTest::testNumberOfRows()
{
  CPPUNIT_ASSERT_EQUAL(10, kinetic_.numberOfRows());
}



