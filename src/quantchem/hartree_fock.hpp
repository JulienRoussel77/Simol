
#ifndef SIMOL_HARTREE_FOCK_HPP
#define	SIMOL_HARTREE_FOCK_HPP

#include "SlaterDeterminant.hpp"
#include "SparseTensor.hpp"
#include <vector>
#include "core/linalg/Vector.hpp"
#include "DiscreteHamiltonian.hpp"
#include "core/linalg/EigenDecomposition.hpp"
namespace simol
{

    double
    electric_integral(std::size_t const M_disc,
                      Vector<double> const & psi1,
                      Vector<double> const & psi2,
                      Vector<double> const & psi3,
                      Vector<double> const & psi4,
                      SparseTensor<double> const & E)
    {

        Vector<double> U12v = Vector<double>::Zero(M_disc*M_disc);
        Vector<double> U34v = Vector<double>::Zero(M_disc*M_disc);

        for (size_t i=0; i< M_disc; i++)
        {
            for (size_t j=0; j< M_disc; j++)
            {
                U12v(getInd(M_disc, i,j)) = psi1(i)*psi3(j);
                U34v(getInd(M_disc, i,j)) = psi2(i)*psi4(j);
            }
        }

        Vector<double> temp = E * U12v;
        return (U34v, temp);
    }

    double
    electric_integral(std::size_t const M_disc,
                      size_t const index1,
                      size_t const index2,
                      size_t const index3,
                      size_t const index4,
                      SparseTensor<double> const & E)
    {
        return E(getInd(M_disc, index2, index4), getInd(M_disc, index1, index3));
    }

    DenseMatrix<double>
    FockMat(std::size_t const M_disc,
            DenseMatrix<double> const & Phi,
            SparseTensor<double> const & E)
    {

        //A améliorer pour faire que la matrice de Fock soit une matrice sparse symmétrique
        DenseMatrix<double> G0 = DenseMatrix<double>::Zero(M_disc, M_disc);
        DenseMatrix<double> D = Phi * Phi.adjoint();

        for (size_t k=0; k< M_disc; k++)
        {
            for (size_t l=k; l< M_disc; l++)
            {
                for (size_t k1=0; k1 < M_disc; k1++)
                {
                    for (size_t k2=0; k2 < M_disc; k2++)
                        G0(k,l) += D(k1,k2) * ( electric_integral(M_disc, k, l, k1, k2, E) - electric_integral(M_disc, k, k2, k1, l, E) );
                }
                G0(k,l) *= 0.5;
                G0(l,k) = G0(k,l);
            }
        }

        return G0;

    }



    double
    H1_slat(SlaterDeterminant const & Phi,
            SlaterDeterminant const & Psi,
            DenseMatrix<double> const & O,
            DenseMatrix<double> const & H,
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
            DenseMatrix<double> temp = S.inverse() * A;
            return S.determinant() * temp.trace();
        }
        else
        {
            Eigen::JacobiSVD<eigen<double>::DenseMatrixType> svd(S.wrapped_, Eigen::ComputeThinU | Eigen::ComputeThinV);

            DenseMatrix<double> Uvec = svd.matrixU();
            DenseMatrix<double> Vvec = svd.matrixV();

            Vector<double> D = svd.singularValues();
            double lmax = D.max();

            D = (1.0/lmax) * D;

            //On compte la multiplicité de la valeur propre nulle
            int mult = 0;
            for (size_t i= 0; i< D.size(); i++)
            {
                if (D(i) < ratio_)
                    ++mult;
            }

            //si la multiplicité de 0 est plus grande que 1, la valeur est 0
            if (mult >1.5)
                return 0;
            //Sinon, on utilise une formule patriculière
            else
            {
                int indmin = D.index_of_minimum();

                Vector<double> xV = Vvec.column(indmin);
                Vector<double> xU = Uvec.column(indmin);

                Vector<double> orthU = U * xU;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                Vector<double> orthV = V * xV;  //Vecteur dans l'orthogonal du sous-espace engendré par U

                Vector<double> temp = O * orthU;
                double n2orthU = (orthU, temp);
                orthU = ( 1.0 / sqrt(n2orthU) ) * orthU;

                temp = O * orthV;
                double n2orthV = (orthV, temp);
                orthV = ( 1.0 / sqrt(n2orthV) ) * orthV;

                DenseMatrix<double> PsiU = DenseMatrix<double>::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: le deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                DenseMatrix<double> PsiV = DenseMatrix<double>::Zero(M_disc_ , numberOfElectrons);

                Vector<double> muU = Vector<double>::Zero(numberOfElectrons);
                Vector<double> muV = Vector<double>::Zero(numberOfElectrons);

                for(size_t k=0; k<numberOfElectrons; k++)
                {
                    temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (U.wrapped_.col(k));
                    double PS = (orthU, temp);
                    PsiU.column(k) = U.column(k) - PS * orthU;
                    muU(k) = PS;

                    temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (V.wrapped_.col(k));
                    PS = (orthV, temp);
                    PsiV.column(k) = V.column(k) - PS * orthV;
                    muV(k) = PS;
                }

                    //On crée une base orthonormale de l'intersection des deux
                //sous-espaces: par exemple à partir de xU
                DenseMatrix<double> xUbas = DenseMatrix<double>::Zero(numberOfElectrons, numberOfElectrons-1);
                xUbas.wrapped_.block(0,0, numberOfElectrons, indmin) = Uvec.wrapped_.block(0,0, numberOfElectrons, indmin);
                xUbas.wrapped_.block(0,indmin, numberOfElectrons, numberOfElectrons-indmin-1) = Uvec.wrapped_.block(0, indmin+1, numberOfElectrons, numberOfElectrons-indmin-1);

                DenseMatrix<double> coeffs_bas = U * xUbas;
                //Puis on orthonormalise les vecteurs de coeffs_bas
                for (size_t i=0; i<(numberOfElectrons-1); i++)
                {
                    Vector<double> temp2(O.numberOfRows());
                    temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                    double PS2 = (orthU, temp2);
                    coeffs_bas.column(i) = coeffs_bas.column(i) - PS2 * orthU;

                    for (size_t j=0; j<i; j++)
                    {
                        temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(j));
                        PS2 = (coeffs_bas.column(i), temp2);
                        coeffs_bas.column(i) = coeffs_bas.column(i) - PS2 * coeffs_bas.column(j);
                    }

                    temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                    PS2 = ((coeffs_bas.wrapped_.col(i)).adjoint())*temp2.wrapped_;
                    coeffs_bas.wrapped_.col(i) = (1.0/sqrt(PS2))*coeffs_bas.wrapped_.col(i);

                }

                DenseMatrix<double> CU = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons);
                DenseMatrix<double> CV = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons);

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

                Vector<double> temp3 = H * orthV;
                double sum2 = (orthU.wrapped_.adjoint())*temp3.wrapped_;

                //On a alors pour tout k U(:,k) = muU(k)*orthU + PsiU(:,k)

                double sum = 0;

                for (size_t k1= 0; k1 < numberOfElectrons; k1++)
                {
                    for (size_t k2 =0; k2 < numberOfElectrons; k2++)
                    {
                        DenseMatrix<double> aU = DenseMatrix<double>::Zero(numberOfElectrons-1, numberOfElectrons-1);
                        aU.wrapped_.block(0,0, numberOfElectrons-1, k1) = CU.wrapped_.block(0,0, numberOfElectrons-1, k1);
                        aU.wrapped_.block(0,k1, numberOfElectrons-1, numberOfElectrons-1-k1) = CU.wrapped_.block(0,k1+1, numberOfElectrons-1, numberOfElectrons -1 -k1);

                        DenseMatrix<double> aV = DenseMatrix<double>::Zero(numberOfElectrons-1, numberOfElectrons-1);
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

    double H2_slat_N2(const SlaterDeterminant& Phi,
                      const SlaterDeterminant& Psi,
                      const SparseTensor<double>& E,
                      size_t const M_disc)
    {
        assert(Phi.number_of_electrons()==2);
        assert(Psi.number_of_electrons()==2);

        DenseMatrix<double> U = Phi.matrix();
        DenseMatrix<double> V = Psi.matrix();

        // On calcule directement

        Vector<double> Vcol0 = V.column(0);
        Vector<double> Vcol1 = V.column(1);
        Vector<double> Ucol0 = U.column(0);
        Vector<double> Ucol1 = U.column(1);

        return 0.5 * ( electric_integral(M_disc, Ucol0, Vcol0, Ucol1, Vcol1, E)
                     - electric_integral(M_disc, Ucol1, Vcol0, Ucol0, Vcol1, E)
                     - electric_integral(M_disc, Ucol0, Vcol1, Ucol1, Vcol0, E)
                     + electric_integral(M_disc, Ucol1, Vcol1, Ucol0, Vcol0, E) );

    }

    double H2_slat(const SlaterDeterminant& Phi,
                   const SlaterDeterminant& Psi,
                   const DenseMatrix<double>& O,
                   const SparseTensor<double>& E,
                   size_t const numberOfElectrons,
                   size_t const M_disc_,
                   double ratio_ = 1e-12)
    {
        if (numberOfElectrons==2)
            return H2_slat_N2(Phi,Psi,E,M_disc_ );

        else
        {
            DenseMatrix<double> S = Smat(Phi, Psi, O);

            DenseMatrix<double> U = Phi.matrix();
            DenseMatrix<double> V = Psi.matrix();

            Eigen::JacobiSVD<eigen<double>::DenseMatrixType> svd(S.wrapped_);

            Vector<double> D = svd.singularValues();
            double lmin = D.min();
            double lmax = D.max();

            double ratio = lmin/lmax;

            //Si le ratio n'est pas trop petit, c'est la formule habituelle
            if (ratio > ratio_)
            {
                DenseMatrix<double> Sinv(S.numberOfRows(), S.numberOfColumns());
                Sinv = S.inverse();
                double scal0 = 0;

                for (size_t i=0; i < numberOfElectrons; i++)
                {
                    Vector<double> Ucol_i = U.column(i);

                    for (size_t j=0; j < numberOfElectrons; j++)
                    {
                        Vector<double> Ucol_j = U.column(j);

                        for (size_t k = 0; k < numberOfElectrons; k++)
                        {
                            Vector<double> Vcol_k = V.column(k);

                            for (size_t l=0; l < numberOfElectrons; l++)
                            {
                                Vector<double> Vcol_l = V.column(l);

                                 //Ici je prends les notations de Friedrichs
                                scal0 += 0.5 * electric_integral(M_disc_, Ucol_i, Vcol_k, Ucol_j, Vcol_l, E)
                                             * (Sinv(k,i)*Sinv(l,j) - Sinv(k,j)*Sinv(l,i));

                            }
                        }
                    }
                }

                return scal0 * (S.determinant());
            }
            else
            {
                Eigen::JacobiSVD<eigen<double>::DenseMatrixType> svd(S.wrapped_, Eigen::ComputeThinU | Eigen::ComputeThinV);

                DenseMatrix<double> Uvec = svd.matrixU();
                DenseMatrix<double> Vvec = svd.matrixV();

                Vector<double> D = svd.singularValues();
                double lmax = D.max();

                D *= (1.0/lmax);

                //On compte la multiplicité de la valeur propre nulle
                int mult = 0;
                for (size_t i= 0; i< numberOfElectrons; i++)
                {
                    if (D(i) < ratio_)
                        ++mult;
                }

                //si la multiplicité de 0 est plus grande que 2, la valeur est 0
                if (mult >2.5)
                    return 0;

                else
                {
                    if ((1.5 > mult) && (mult > 0.5)) //mult = 1
                    {

                        int indmin = D.index_of_minimum();
                        Vector<double> xV = Vvec.column(indmin);
                        Vector<double> xU = Uvec.column(indmin);

                        Vector<double> orthU = U * xU;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                        Vector<double> orthV = V * xV;  //Vecteur dans l'orthogonal du sous-espace engendré par U

                        Vector<double> temp = O * orthU;
                        orthU = 1.0 / sqrt( (orthU, temp) ) * orthU;
                        temp = O * orthV;
                        orthV = 1.0 / sqrt( (orthV, temp) ) * orthV;


                        DenseMatrix<double> PsiU = DenseMatrix<double>::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: le deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                        DenseMatrix<double> PsiV = DenseMatrix<double>::Zero(M_disc_,numberOfElectrons);

                        Vector<double> muU = Vector<double>::Zero(numberOfElectrons);
                        Vector<double> muV = Vector<double>::Zero(numberOfElectrons);

                        for(size_t k=0; k<numberOfElectrons; k++)
                        {
                            temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (U.wrapped_.col(k));
                            double PS = (orthU, temp);
                            PsiU.column(k) = U.column(k) - PS * orthU;
                            muU(k) = PS;

                            temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (V.wrapped_.col(k));
                            PS = (orthV, temp);
                            PsiV.column(k) = V.column(k) - PS * orthV;
                            muV(k) = PS;
                        }

                        //On crée une base orthonormale de l'intersection des deux
                        //sous-espaces: par exemple à partir de xU
                        DenseMatrix<double> xUbas = DenseMatrix<double>::Zero(numberOfElectrons, numberOfElectrons-1);
                        xUbas.wrapped_.block(0,0, numberOfElectrons, indmin) = Uvec.wrapped_.block(0,0, numberOfElectrons, indmin);
                        xUbas.wrapped_.block(0,indmin, numberOfElectrons, numberOfElectrons-indmin-1) = Uvec.wrapped_.block(0, indmin+1, numberOfElectrons, numberOfElectrons-indmin-1);

                        DenseMatrix<double> coeffs_bas = U * xUbas;
                        //Puis on orthonormalise les vecteurs de coeffs_bas
                        for (size_t i=0; i<(numberOfElectrons-1); i++)
                        {
                            Vector<double> temp2(O.numberOfColumns());
                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                            double PS2 = (orthU, temp2);
                            coeffs_bas.column(i) = coeffs_bas.column(i) - PS2 * orthU;

                            for (size_t j=0; j<i; j++)
                            {
                                temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(j));
                                PS2 = (coeffs_bas.column(i), temp2);
                                coeffs_bas.column(i) = coeffs_bas.column(i) - PS2 * coeffs_bas.column(j);
                            }

                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(i));
                            PS2 = (coeffs_bas.column(i), temp2);
                            coeffs_bas.column(i) = (1.0/sqrt(PS2)) * coeffs_bas.column(i);

                        }

                        DenseMatrix<double> CU = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons);
                        DenseMatrix<double> CV = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons);

                        for (size_t k=0; k< numberOfElectrons; k++)
                        {
                            for (size_t l=0; l< (numberOfElectrons-1); l++)
                            {
                                Vector<double> temp(O.numberOfRows());
                                temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(l));
                                double PS = (PsiU.column(k), temp);
                                CU(l,k) = PS;

                                temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(l));
                                PS = (PsiV.column(k), temp);
                                CV(l,k) = PS;
                            }
                        }

                        //On construit la somme
                        double sum2 = 0;
                        for (size_t l=0; l< numberOfElectrons-1; l++)
                            sum2 += electric_integral(M_disc_,coeffs_bas.column(l),coeffs_bas.column(l),orthU,orthV,E)
                                  - electric_integral(M_disc_,coeffs_bas.column(l),orthV,coeffs_bas.column(l),orthU,E);


                        //On a alors pour tout k U(:,k) = muU(k)*orthU + PsiU(:,k)
                        double sum = 0;

                        for (size_t k1= 0; k1 < numberOfElectrons; k1++)
                        {
                            for (size_t k2 = 0; k2 < numberOfElectrons; k2++)
                            {
                                DenseMatrix<double> aU = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons-1);
                                aU.wrapped_.block(0,0,numberOfElectrons-1,k1) = CU.wrapped_.block(0,0,numberOfElectrons-1,k1);
                                aU.wrapped_.block(0,k1,numberOfElectrons-1,numberOfElectrons-1-k1) = CU.wrapped_.block(0,(k1+1),numberOfElectrons-1, numberOfElectrons-1-k1);

                                DenseMatrix<double> aV = DenseMatrix<double>::Zero(numberOfElectrons-1,numberOfElectrons-1);
                                aV.wrapped_.block(0,0,numberOfElectrons-1,k2) = CV.wrapped_.block(0,0,numberOfElectrons-1,k2);
                                aV.wrapped_.block(0,k2,numberOfElectrons-1, numberOfElectrons-1-k2) = CV.wrapped_.block(0,k2+1,numberOfElectrons-1, numberOfElectrons-1-k2);

                                sum += muU(k1)*muV(k2)*pow(-1,k1)*pow(-1,k2)*sum2*(aU.determinant())*(aV.determinant());
                            }

                        }

                        return sum;
                    }

                    else
                    {	//mult = 2
                        //La dimension de l'intersection des deux espaces est N-2: il
                        //faut récupérer pour U et pour V les deux vecteurs qui sont
                        //orthonormaux à l'espace

                        size_t indmin = D.index_of_minimum();
                        Vector<double> xV1 = Vvec.column(indmin);
                        Vector<double> xU1 = Uvec.column(indmin);

                        D(indmin) = 1e20;

                        size_t indmin2 = D.index_of_minimum();
                        Vector<double> xV2 = Vvec.column(indmin2);
                        Vector<double> xU2 = Uvec.column(indmin2);

                        Vector<double> orthU1 = U * xU1;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                        Vector<double> orthU2 = U * xU2;
                        Vector<double> orthV1 = V * xV1;
                        Vector<double> orthV2 = V * xV2;

                        Vector<double> temp = O * orthU1;
                        orthU1 = ( 1.0 / sqrt( (orthU1, temp) ) ) * orthU1;

                        temp = O * orthU2;
                        double PS = (orthU1, temp);
                        orthU2 -= PS * orthU1;
                        temp = O * orthU2;
                        orthU2 = ( 1.0 / sqrt( (orthU2, temp) ) ) * orthU2;

                        temp = O * orthV1;
                        orthV1 = ( 1.0 / sqrt( (orthV1, temp) ) ) * orthV1;

                        temp = O * orthV2;
                        PS = (orthV1, temp);
                        orthV2 -= PS * orthV1;
                        temp = O * orthV2;
                        orthV2 = ( 1.0 / sqrt( (orthV2, temp) ) ) * orthV2;

                        DenseMatrix<double> PsiU = DenseMatrix<double>::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: les deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                        DenseMatrix<double> PsiV = DenseMatrix<double>::Zero(M_disc_,numberOfElectrons);

                        Vector<double> muU1 = Vector<double>::Zero(numberOfElectrons);
                        Vector<double> muV1 = Vector<double>::Zero(numberOfElectrons);
                        Vector<double> muU2 = Vector<double>::Zero(numberOfElectrons);
                        Vector<double> muV2 = Vector<double>::Zero(numberOfElectrons);

                        for (size_t k=0; k< numberOfElectrons; k++)
                        {

                            Vector<double> temp2(O.numberOfRows());
                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(U.wrapped_.col(k));
                            muU1(k) = (orthU1, temp2);
                            muU2(k) = (orthU2, temp2);
                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(V.wrapped_.col(k));
                            muV1(k) = (orthV1, temp2);
                            muV2(k) = (orthV2, temp2);


                            PsiU.column(k) = U.column(k) - muU1(k) * orthU1 - muU2(k) * orthU2;
                            PsiV.column(k) = V.column(k) - muV1(k) * orthV1 - muV2(k) * orthV2;
                        }

                        //On crée une base orthonormale de l'intersection des deux
                        //sous-espaces: par exemple à partir de xU
                        DenseMatrix<double> xUbas = DenseMatrix<double>::Zero(numberOfElectrons, numberOfElectrons-2);
                        int ind = 0;
                        for (size_t k=0; k< numberOfElectrons; k++)
                        {
                            if ( (k!=indmin) && (k!=indmin2))
                            {
                                xUbas.column(ind) = Uvec.column(k);
                                ++ind;
                            }
                        }

                        DenseMatrix<double> coeffs_bas = U * xUbas;

                        //Puis on orthonormalise les vecteurs de coeffs_bas
                        for (size_t i= 0; i < numberOfElectrons-2; i++)
                        {
                            //On commence par orthonormaliser par rapport à orthU1 et orthU2
                            Vector<double> temp3(O.numberOfRows());
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            double PS3 = (orthU1, temp3);
                            coeffs_bas.column(i) = coeffs_bas.column(i) - PS3 * orthU1;
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            PS3 = (orthU2, temp3);
                            coeffs_bas.column(i) = coeffs_bas.column(i) - PS3 * orthU2;
                            for (size_t j=0; j<i; j++)
                            {
                                temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                                PS3 = (coeffs_bas.column(j), temp3);
                                coeffs_bas.column(i) = coeffs_bas.column(i) - PS3 * coeffs_bas.column(j);
                            }
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            coeffs_bas.column(i) = ( 1.0 / sqrt( (coeffs_bas.column(i), temp3) ) ) * coeffs_bas.column(i);
                        }

                        DenseMatrix<double> CU = DenseMatrix<double>::Zero(numberOfElectrons-2,numberOfElectrons);
                        DenseMatrix<double> CV = DenseMatrix<double>::Zero(numberOfElectrons-2,numberOfElectrons);

                        for (size_t k=0; k< numberOfElectrons; k++)
                        {
                            for (size_t l=0; l< numberOfElectrons-2; l++)
                            {
                                Vector<double> temp4(O.numberOfRows());
                                temp4.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*coeffs_bas.wrapped_.col(l);
                                CU(l,k) = (PsiU.column(k), temp4);
                                CV(l,k) = (PsiV.column(k), temp4);
                            }
                        }

                        //On construit la somme
                        double sum2 = electric_integral(M_disc_, orthU1,orthV1,orthU2,orthV2,E) - electric_integral(M_disc_, orthU1,orthV2,orthU2,orthV1, E);

                        //On a alors pour tout k U(:,k) = muU1(k)*orthU1 + muU2(k)*orthU2 + PsiU(:,k)
                        double sum = 0;

                        for (size_t iu = 0; iu < numberOfElectrons-1; iu++)
                        {
                            for (size_t ju = iu+1 ; ju <numberOfElectrons; ju++)
                            {
                                DenseMatrix<double> aU = DenseMatrix<double>::Zero(numberOfElectrons-2,numberOfElectrons-2);
                                aU.wrapped_.block(0,0,numberOfElectrons-2, iu) = CU.wrapped_.block(0,0,numberOfElectrons-2,iu);
                                aU.wrapped_.block(0,iu,numberOfElectrons-2, ju-iu-1) = CU.wrapped_.block(0,iu+1, numberOfElectrons-2, ju-iu-1);
                                aU.wrapped_.block(0,ju-1,numberOfElectrons-2, numberOfElectrons-ju-1) = CU.wrapped_.block(0,ju+1,numberOfElectrons-2, numberOfElectrons-ju-1);


                                for (size_t iv = 0; iv < numberOfElectrons-1; iv++)
                                {
                                    for (size_t jv = iv+1; jv < numberOfElectrons; jv++)
                                    {

                                        DenseMatrix<double> aV = DenseMatrix<double>::Zero(numberOfElectrons-2,numberOfElectrons-2);
                                        aV.wrapped_.block(0,0,numberOfElectrons-2,iv) = CV.wrapped_.block(0,0, numberOfElectrons-2, iv);
                                        aV.wrapped_.block(0,iv,numberOfElectrons-2,jv-iv-1) = CV.wrapped_.block(0,iv+1,numberOfElectrons-2,jv-iv-1);
                                        aV.wrapped_.block(0,jv-1,numberOfElectrons-2,numberOfElectrons-jv-1) = CV.wrapped_.block(0,jv+1,numberOfElectrons-2,numberOfElectrons-jv-1);


                                        sum += pow(-1,iu) * pow(-1,ju+1) *pow(-1,iv) *pow(-1,jv+1)
                                             * ((muU1(iu)*muU2(ju) -muU1(ju)*muU2(iu)) * (muV1(iv)*muV2(jv) -muV1(jv)*muV2(iv)))
                                             * sum2*aU.determinant()*aV.determinant();

                                    }
                                }

                            }
                        }

                        return sum;
                    }
                }
            }
        }

    }


    double
    over_slat(SlaterDeterminant const & Phi, SlaterDeterminant const & Psi, DenseMatrix<double> const & O)
    { return Smat(Phi, Psi, O).determinant(); }

    SlaterDeterminant
    hartree_fock(DiscreteHamiltonian H,
                 std::size_t const numberOfElectrons,
                 SlaterDeterminant const & initial_solution,
                 std::size_t const numberOfIterations)
    {
        //A améliorer: pour l'instant j'utilise des matrices pleines car il n'y a pas de solveur eigenvalue pour les matrices sparses

        DenseMatrix<double> Phi0 = initial_solution.matrix();

        DenseMatrix<double> K = H.kinetic();
        DenseMatrix<double> O = H.overlap();
        DenseMatrix<double> Nu = H.potential();

        DenseMatrix<double> F0 = K + Nu;

        //Roothan parce que c'est le plus simple: ToDo coder ODA
        for( std::size_t iteration = 0; iteration < numberOfIterations; ++iteration )
        {
            DenseMatrix<double> F = F0 + FockMat(H.basisDimension(), Phi0, H.two_electrons());

            DenseMatrix<double> densmat = Phi0 * Phi0.adjoint();
	        DenseMatrix<double> prod = F * densmat;
	        double lambda = prod.trace();
            std::cout << "lambda = " << lambda << std::endl;

            EigenDecomposition<double> es(F, O);

            Vector<double> D = es.eigenvalues();
            DenseMatrix<double> V = es.eigenvectors();

            std::vector<size_t> Itab = D.indices_of_smallest(numberOfElectrons);

            Phi0 = V.permute_columns(Itab);

            double lambda2 = H1_slat(Phi0, Phi0, O, F0, numberOfElectrons, H.basisDimension()) + H2_slat(Phi0, Phi0, O, H.two_electrons(), numberOfElectrons, H.basisDimension());
            lambda2 /= over_slat(Phi0, Phi0, O);

            std::cout << "lambda2 = " << lambda2 << std::endl;

        }

        return SlaterDeterminant(Phi0);
    }

}


#endif

