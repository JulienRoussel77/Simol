#include "gtest/gtest.h"
#include "hartree_fock.hpp"

namespace simol
{

  TEST(HartreeFockTest, LiH)
  {
    std::string const dataPath = "/local/cdoucet/simol/test/functional/quantchem/LiH/";
    simol::MatrixMarketFile kinetic_file(dataPath + "kinetic_matrix.mtx");
    simol::MatrixMarketFile overlap_file(dataPath + "overlap_matrix.mtx");
    simol::MatrixMarketFile potential_file(dataPath + "potential_matrix.mtx");

    simol::DenseMatrix<double> kinetic_matrix(kinetic_file);
    simol::DenseMatrix<double> overlap_matrix(overlap_file);
    simol::DenseMatrix<double> potential_matrix(potential_file);
    //simol::SparseTensor<double> twoelectron_matrix(dataPath + "twoelectron_matrix.txt", M_disc);*/


    EXPECT_EQ(1,1);

  }

}

