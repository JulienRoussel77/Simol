
#ifndef SIMOL_HARTREE_FOCK_HPP
#define	SIMOL_HARTREE_FOCK_HPP

#include "SlaterDeterminant.hpp"

#include <vector>
#include "../linalg/Vector.hpp"

namespace simol
{
    
    size_t getIndMin(Vector<double> const & W)
    {
        size_t ind = 0;
        double minw = W(0);

        size_t w = W.size();

        for (size_t i = 1; i< w; i++)
        {
            if (W(i)<minw)
            {
                ind = i;
                minw = W(i);
            }
        }
        return ind;
    }

    std::vector<size_t> 
    getIndMin(Vector<double> const & W0, size_t const N)
    {
        std::vector<size_t> vec(N);

        Vector<double> W = W0;
        size_t num = 0;
        while (num < N)
        {
            size_t ind = 0;
            double minw = W(0);

            size_t w = W.size();

            for (size_t i = 1; i< w; i++)
            {
                if (W(i)<minw)
                {
                    ind = i;
                    minw = W(i);
                }
            }
            vec[num] = ind;
            W(ind) = 1e32;

            num = num +1;

        }
        return vec;
    }
    
    std::size_t 
    getInd(std::size_t const M_disc,
           std::size_t const rowIndex,
           std::size_t const columnIndex)
    {
        assert(rowIndex < M_disc && columnIndex < M_disc);
        return rowIndex * M_disc + columnIndex;
    }
  
    
       DenseMatrix<double>
    FockMat(std::size_t const M_disc,
            DenseMatrix<double> const & Phi, 
            DenseMatrix<double> const & E)
    {

        //A améliorer pour faire que la matrice de Fock soit une matrice sparse symmétrique
        DenseMatrix<double> G0(M_disc, M_disc);
        G0.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc, M_disc);
        DenseMatrix<double> D(Phi.number_of_rows(), Phi.number_of_rows());
        D.wrapped_ = Phi.wrapped_ * Phi.wrapped_.adjoint();

        for (size_t k=0; k< M_disc; k++)
        {
            for (size_t l=k; l< M_disc; l++)
            {
                for (size_t k1=0; k1 < M_disc; k1++)
                {
                    for (size_t k2=0; k2 < M_disc; k2++)
                        G0(k,l) += D(k1,k2)*(E(getInd(M_disc,k,l),getInd(M_disc,k1,k2)) - E(getInd(M_disc,k,k2),getInd(M_disc,k1,l)));
                }
                G0(l,k) = G0(k,l);
            }
        }

        return G0;

    }
    
       
    double 
    elint(std::size_t const M_disc,
          Vector<double> const & psi1, 
          Vector<double> const & psi2, 
          Vector<double> const & psi3,
          Vector<double> const & psi4, 
          SparseMatrix<double> const & E)
    {

        Vector<double> U12v(M_disc*M_disc);
        U12v.wrapped_ = eigen<double>::VectorType::Zero(M_disc*M_disc);
        Vector<double> U34v(M_disc*M_disc);
        U34v.wrapped_ = eigen<double>::VectorType::Zero(M_disc*M_disc);

        for (size_t i=0; i< M_disc; i++)
        {
            for (size_t j=0; j< M_disc; j++)
            {
                U12v(getInd(M_disc, i,j)) = psi1(i)*psi3(j);
                U34v(getInd(M_disc, i,j)) = psi2(i)*psi4(j);
            }
        }

        Vector<double> temp(E.numberOfRows());
        temp.wrapped_ = E.wrapped_.selfadjointView<Eigen::Upper>() * U12v.wrapped_;
        return U34v.wrapped_.adjoint() * temp.wrapped_;
    }
    
    
    double
    H1_slat(SlaterDeterminant const & Phi, 
            SlaterDeterminant const & Psi, 
            SparseMatrix<double> const & O, 
            SparseMatrix<double> const & H,
            size_t const M_disc_,
            size_t const numberOfElectrons,
            double ratio_ = 1e-12)
    {
        DenseMatrix<double> U = Phi.matrix();
        DenseMatrix<double> V = Psi.matrix();

        DenseMatrix<double> S = Smat(Phi,Psi,O);
        DenseMatrix<double> A = Smat(Phi,Psi,H);

        //Si le ratio n'est pas trop petit, c'est la formule habituelle
        if (S.rcond() > ratio_)
        {
            DenseMatrix<double> temp(S.number_of_rows(), A.number_of_columns());
            temp.wrapped_ = S.wrapped_.inverse() * A.wrapped_;
            return (S.wrapped_.determinant())*(temp.wrapped_.trace());
        }
        else
        {
            Eigen::JacobiSVD<eigen<double>::DenseMatrixType> svd(S.wrapped_, Eigen::ComputeThinU | Eigen::ComputeThinV);

            DenseMatrix<double> Uvec = svd.matrixU();
            DenseMatrix<double> Vvec = svd.matrixV();

            Vector<double> D = svd.singularValues();
            double lmax = D.max();

            D.wrapped_ = (1.0/lmax) * D.wrapped_;

            //On compte la multiplicité de la valeur propre nulle
            int mult = 0;
            for (size_t i= 0; i< D.size(); i++)
            {
                if (D(i) < ratio_) 
                    mult = mult+1;
            }

            //si la multiplicité de 0 est plus grande que 1, la valeur est 0
            if (mult >1.5) 
                return 0;
            //Sinon, on utilise une formule patriculière
            else
            {
                int indmin = getIndMin(D);

                Vector<double> xV(Vvec.number_of_rows());
                xV.wrapped_ = Vvec.wrapped_.col(indmin);
                Vector<double> xU(Uvec.number_of_rows());
                xU.wrapped_ = Uvec.wrapped_.col(indmin);

                Vector<double> orthU(U.number_of_rows());
                orthU.wrapped_ = U.wrapped_ * xU.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                Vector<double> orthV(V.number_of_rows());
                orthV.wrapped_ = V.wrapped_ * xV.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par U

                Vector<double> temp(O.numberOfRows());
                temp.wrapped_ = O.wrapped_.selfadjointView<Eigen::Upper>() * orthU.wrapped_;
                double n2orthU = orthU.wrapped_.adjoint() * temp.wrapped_;
                orthU.wrapped_ = (1.0/sqrt(n2orthU))* orthU.wrapped_;

                temp.wrapped_ = O.wrapped_.selfadjointView<Eigen::Upper>() * orthV.wrapped_;
                double n2orthV = orthV.wrapped_.adjoint() * temp.wrapped_;
                orthV.wrapped_ = (1.0/sqrt(n2orthV)) * orthV.wrapped_;

                DenseMatrix<double> PsiU(M_disc_, numberOfElectrons);
                PsiU.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: le deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                DenseMatrix<double> PsiV(M_disc_, numberOfElectrons);
                PsiV.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_ , numberOfElectrons);

                Vector<double> muU(numberOfElectrons);
                muU.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);
                Vector<double> muV(numberOfElectrons);
                muV.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);

                for(size_t k=0; k<numberOfElectrons; k++)
                {
                    temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (U.wrapped_.col(k));
                    double PS = orthU.wrapped_.adjoint()*temp.wrapped_;
                    PsiU.wrapped_.col(k) = U.wrapped_.col(k) - PS*orthU.wrapped_;
                    muU(k) = PS;

                    temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (V.wrapped_.col(k));
                    PS = orthV.wrapped_.adjoint()*temp.wrapped_;
                    PsiV.wrapped_.col(k) = V.wrapped_.col(k) - PS*orthV.wrapped_;
                    muV(k) = PS;
                }

                    //On crée une base orthonormale de l'intersection des deux
                //sous-espaces: par exemple à partir de xU
                DenseMatrix<double> xUbas(numberOfElectrons, numberOfElectrons-1);
                xUbas.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons, numberOfElectrons-1);
                xUbas.wrapped_.block(0,0, numberOfElectrons, indmin) = Uvec.wrapped_.block(0,0, numberOfElectrons, indmin);
                xUbas.wrapped_.block(0,indmin, numberOfElectrons, numberOfElectrons-indmin-1) = Uvec.wrapped_.block(0, indmin+1, numberOfElectrons, numberOfElectrons-indmin-1);

                DenseMatrix<double> coeffs_bas(U.number_of_rows(), xUbas.number_of_columns());
                coeffs_bas.wrapped_ = U.wrapped_*xUbas.wrapped_;
                //Puis on orthonormalise les vecteurs de coeffs_bas
                for (size_t i=0; i<(numberOfElectrons-1); i++)
                {
                    Vector<double> temp2(O.numberOfRows());
                    temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                    double PS2 = orthU.wrapped_.adjoint()*temp2.wrapped_;
                    coeffs_bas.wrapped_.col(i) = coeffs_bas.wrapped_.col(i) - PS2*orthU.wrapped_;

                    for (size_t j=0; j<i; j++)
                    {
                        temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(j));
                        PS2 = ((coeffs_bas.wrapped_.col(i)).adjoint())*temp2.wrapped_;
                        coeffs_bas.wrapped_.col(i) = coeffs_bas.wrapped_.col(i) - PS2*coeffs_bas.wrapped_.col(j);
                    }

                    temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                    PS2 = ((coeffs_bas.wrapped_.col(i)).adjoint())*temp2.wrapped_;
                    coeffs_bas.wrapped_.col(i) = (1.0/sqrt(PS2))*coeffs_bas.wrapped_.col(i);

                }

                DenseMatrix<double> CU(numberOfElectrons-1,numberOfElectrons);
                CU.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1,numberOfElectrons);
                DenseMatrix<double> CV(numberOfElectrons-1,numberOfElectrons);
                CV.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1,numberOfElectrons);

                for (size_t k=0; k< numberOfElectrons; k++)
                {
                    for (size_t l=0; l< (numberOfElectrons-1); l++)
                    {
                        Vector<double> temp(O.numberOfRows());
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(l));
                        double PS = ((PsiU.wrapped_.col(k)).adjoint())*temp.wrapped_;
                        CU.wrapped_(l,k) = PS;

                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(l));
                        PS = ((PsiV.wrapped_.col(k)).adjoint())*temp.wrapped_;
                        CV.wrapped_(l,k) = PS;
                    }
                }

                Vector<double> temp3(H.numberOfRows());
                temp3.wrapped_ = (H.wrapped_.selfadjointView<Eigen::Upper>()) *orthV.wrapped_;
                double sum2 = (orthU.wrapped_.adjoint())*temp3.wrapped_;

                //On a alors pour tout k U(:,k) = muU(k)*orthU + PsiU(:,k)

                double sum = 0;

                for (size_t k1= 0; k1 < numberOfElectrons; k1++)
                {
                    for (size_t k2 =0; k2 < numberOfElectrons; k2++)
                    {
                        DenseMatrix<double> aU(numberOfElectrons-1, numberOfElectrons-1);
                        aU.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1, numberOfElectrons-1);
                        aU.wrapped_.block(0,0, numberOfElectrons-1, k1) = CU.wrapped_.block(0,0, numberOfElectrons-1, k1);
                        aU.wrapped_.block(0,k1, numberOfElectrons-1, numberOfElectrons-1-k1) = CU.wrapped_.block(0,k1+1, numberOfElectrons-1, numberOfElectrons -1 -k1);

                        DenseMatrix<double> aV(numberOfElectrons-1, numberOfElectrons-1);
                        aV.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1, numberOfElectrons-1);
                        aV.wrapped_.block(0,0, numberOfElectrons-1,k2) = CV.wrapped_.block(0,0, numberOfElectrons-1, k2);
                        aV.wrapped_.block(0,k2, numberOfElectrons-1, numberOfElectrons-1-k2) = CV.wrapped_.block(0,k2+1, numberOfElectrons-1, numberOfElectrons-1 -k2);

                         sum = sum + muU(k1)*muV(k2)*(pow(-1, k1))*(pow(-1, k2))*sum2*(aU.wrapped_.determinant())*(aV.wrapped_.determinant());
                    }
                }
                return sum;

            // *********************** Fin de la formule particulière
            }

        }
    }
   

    
    
 
    

/*
    SlaterDeterminant
    hartree_fock(std::size_t const M_disc,
                 std::size_t const numberOfElectrons,
                 DenseMatrix<double> const & E,
                 DenseMatrix<double> const & K,
                 DenseMatrix<double> const & O,
                 DenseMatrix<double> const & Nu,
                 SlaterDeterminant<double> const & initial_solution,
                 std::size_t const numberOfIterations)
    {
        //A améliorer: pour l'instant j'utilise des matrices pleines car il n'y a pas de solveur eigenvalue pour les matrices sparses

        DenseMatrix<double> Phi0 = initial_solution.Phi_;
        DenseMatrix<double> F0 = K + Nu;
        
        //Roothan parce que c'est le plus simple: ToDo coder ODA
        for( std::size_t iteration = 0; iteration < numberOfIterations; ++iteration )
        {
            DenseMatrix<double> F = F0 + FockMat(Phi0,E);


 //           GeneralizedSelfAdjointEigenSolver<MatrixXd> es(F, O, ComputeEigenvectors|Ax_lBx);

 //           Vector<double> D = es.eigenvalues();
//            DenseMatrix<double> V = es.eigenvectors();
            
            Vector<double> D;
            DenseMatrix<double> V;

            std::vector<int> Itab = getIndMin(D, numberOfElectrons);

            DenseMatrix<double> Phinew = DenseMatrix<double>::Zero(M_disc, numberOfElectrons);
            for (int i=0; i< numberOfElectrons; i++)
                Phinew.col(i) = V.col(Itab[i]);
            Phi0 = Phinew;


            double lambda = 0;
            for (std::size_t i = 0; i < numberOfElectrons; ++i)
            {
                lambda += D(Itab[i]);
                for (std::size_t j = 0; j< numberOfElectrons; ++j)
                    lambda += 0.5*(elint(M_disc, Phi0.col(i), Phi0.col(i), Phi0.col(j), Phi0.col(j), E) - elint(M_disc, Phi0.col(i), Phi0.col(j), Phi0.col(j), Phi0.col(i), E));
            }

            SlaterDeterminant<double> sol(Phi0);
            double lambda2 = calc_.H1_slat(Phi0, Phi0, O, K + Nu) + calc_.H2_slat(Phi0, Phi0, O, E);
            lambda2 /= calc_.over_slat(Phi0, Phi0, O);


        }

        SlaterDeterminant sol(Phi0);

        //sol_ = sol;
    }
*/
}


#endif	

