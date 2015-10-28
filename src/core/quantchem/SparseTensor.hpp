

#ifndef SIMOL_SPARSETENSOR_HPP
#define	SIMOL_SPARSETENSOR_HPP

#include <cassert>
#include <fstream>
#include <string>

#include <Eigen/SparseCore>

namespace simol
{
    std::size_t 
    getInd(std::size_t const M_disc,
           std::size_t const rowIndex,
           std::size_t const columnIndex)
    {
        assert(rowIndex < M_disc && columnIndex < M_disc);
        return rowIndex * M_disc + columnIndex;
    }
    
    template<typename ScalarType = double>
    class SparseTensor
    {
        public:
            SparseTensor(std::string const & filename, size_t const M_disc)
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
                    nonzeros_.push_back(CoefficientType(getInd(M_disc, i,k), getInd(M_disc, j,l),t0));
                    nbLignes++;
                }
                in.close(); //On ferme le fichier
            }
            
        private:
            typedef Eigen::Triplet<ScalarType, size_t> CoefficientType;
            std::vector<CoefficientType> nonzeros_;
    };
}




#endif	/* SIMOL_SPARSETENSOR_HPP */

