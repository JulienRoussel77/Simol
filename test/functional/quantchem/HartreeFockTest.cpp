#include "gtest/gtest.h"
#include "hartree_fock.hpp"

namespace simol
{

  TEST(HartreeFockTest, LiH)
  {
    std::string const dataPath = "/local/cdoucet/simol/test/functional/quantchem/LiH/";
    MatrixMarketFile kinetic_file(dataPath + "kinetic_matrix.mtx");
    MatrixMarketFile overlap_file(dataPath + "overlap_matrix.mtx");
    MatrixMarketFile potential_file(dataPath + "potential_matrix.mtx");

    DenseMatrix<double> kinetic_matrix(kinetic_file);
    std::cout << "K apres lecture:" << std::endl;
    std::cout << kinetic_matrix << std::endl;
    DenseMatrix<double> overlap_matrix(overlap_file);
    DenseMatrix<double> potential_matrix(potential_file);

    std::string filename = dataPath + "n.txt";
    FILE* size_file = fopen(filename.c_str(), "r");
    size_t M_disc;
    fscanf(size_file, "%lu", &M_disc);

    std::cout << "M_disc = " << M_disc << std::endl;

    SparseTensor<double> twoelectron_matrix(dataPath + "twoelectron_matrix.txt", M_disc);
   std::size_t numberOfElectrons = 2;

   DenseMatrix<double> matrix(M_disc,M_disc);
   matrix.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc, numberOfElectrons);
   SlaterDeterminant initial_solution(matrix);

   std::size_t numberOfIterations = 20;

   SlaterDeterminant solution = hartree_fock(M_disc, numberOfElectrons,
                                             twoelectron_matrix,
                                             kinetic_matrix,
                                             overlap_matrix,
                                             potential_matrix,
                                             initial_solution,
                                             numberOfIterations);



    EXPECT_EQ(1,1);

  }

}

