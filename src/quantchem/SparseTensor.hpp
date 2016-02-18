

#ifndef SIMOL_SPARSETENSOR_HPP
#define	SIMOL_SPARSETENSOR_HPP

#include <cassert>
#include <fstream>
#include <string>

#include <Eigen/SparseCore>

#include "core/linalg/SparseMatrix.hpp"

namespace simol
{
    std::size_t
    getInd(std::size_t const M_disc,
           std::size_t const rowIndex,
           std::size_t const columnIndex);


    template<typename ScalarType = double>
    class SparseTensor
    {
        public:

            SparseTensor(std::string const & filename,
                         size_t const M_disc);

            SparseMatrix<ScalarType> const &
            nonzeros() const;

            ScalarType const operator()(std::size_t const rowIndex,
                                          std::size_t const columnIndex) const
            { return nonzeros_(rowIndex, columnIndex); }

            std::size_t numberOfRows() const
            { return nonzeros_.numberOfRows(); }

        public:
            SparseMatrix<ScalarType> nonzeros_;
    };

    template<typename ScalarType>
    SparseTensor<ScalarType>::SparseTensor(std::string const & filename,
                                           size_t const M_disc)
    : nonzeros_(M_disc * M_disc, M_disc * M_disc)
    {
        // Attention: le fichier contient les valeurs de \int chi_i(x) chi_j(x) \chi_k(y) \chi_l(y)/|x-y|
        // uniquement pour les valeurs de i>=j et k>=l
        // Donc, comme on veut une matrice qui en contiennent que la partie supérieure de la matrice, il faut
        // que la ligne soit donnée par (j,l) et la colonne par (i,k).

        FILE* fichier = fopen(filename.c_str(),"r" ); //ON ouvre le fichier en lecture seule

        std::ifstream in(filename); //Ouverture en mode lecture de "bdd.txt"
        std::string ligne; //Création d'une chaine de caractere
        int nbLignes = 0;
        int test = 0;

        typedef Eigen::Triplet<ScalarType, size_t> NonZero;
        std::vector<NonZero> nonzeros;
        while(std::getline(in, ligne))
        {
            int i;
            int j;
            int k;
            int l;
            long double t;
            test = 0;
            test = fscanf(fichier, "%d %d %d %d %Lf", &i , &j, &k, &l, &t);
            assert(test>0);
            double t0 =t;

            // A checker
            nonzeros.push_back(NonZero(getInd(M_disc, i,k), getInd(M_disc, j,l),t0));
            nbLignes++;
        }
        nonzeros_.wrapped_.setFromTriplets(nonzeros.begin(), nonzeros.end());
        in.close(); //On ferme le fichier
    }

    template<typename ScalarType>
    SparseMatrix<ScalarType> const &
    SparseTensor<ScalarType>::nonzeros() const
    { return nonzeros_; }
}




#endif	/* SIMOL_SPARSETENSOR_HPP */

