#include "MatrixMarketFileTest.hpp"

#include <cppunit/extensions/HelperMacros.h>


namespace simol
{

  CPPUNIT_TEST_SUITE_REGISTRATION(MatrixMarketFileTest);


  void MatrixMarketFileTest::setUp()
  {
    kinetic_ = new simol::MatrixMarketFile("../src/test/data/identity5.mtx");
  }

  void MatrixMarketFileTest::tearDown()
  {
  //  delete kinetic_;
  }

  void MatrixMarketFileTest::testNonExistentFilenameThrowsException()
  {
    bool exceptionCatched  = false;
    try
    {
      MatrixMarketFile nonexistent("Does not exist");
    }
    catch(std::ios_base::failure & error)
    {
      exceptionCatched = true;
    }
    CPPUNIT_ASSERT(exceptionCatched);
  }

  void MatrixMarketFileTest::testNumberOfRows()
  {
    std::size_t actual_size = kinetic_->numberOfRows();
    std::size_t expected_size = 5;
    CPPUNIT_ASSERT_EQUAL(expected_size, actual_size);
  }

}
