
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
    
    
  
    
       DenseMatrix<double>
    FockMat(std::size_t const M_disc,
            DenseMatrix<double> const & Phi, 
            DenseMatrix<double> const & E)
    {

        //A améliorer pour faire que la matrice de Fock soit une matrice sparse symmétrique
        DenseMatrix<double> G0(M_disc, M_disc);
        G0.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc, M_disc);
        DenseMatrix<double> D(Phi.numberOfRows(), Phi.numberOfRows());
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
          DenseMatrix<double> const & E)
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
            DenseMatrix<double> temp(S.numberOfRows(), A.numberOfColumns());
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

                Vector<double> xV(Vvec.numberOfRows());
                xV.wrapped_ = Vvec.wrapped_.col(indmin);
                Vector<double> xU(Uvec.numberOfRows());
                xU.wrapped_ = Uvec.wrapped_.col(indmin);

                Vector<double> orthU(U.numberOfRows());
                orthU.wrapped_ = U.wrapped_ * xU.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                Vector<double> orthV(V.numberOfRows());
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

                DenseMatrix<double> coeffs_bas(U.numberOfRows(), xUbas.numberOfColumns());
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
    
    double H2_slat_N2(const SlaterDeterminant& Phi, 
                      const SlaterDeterminant& Psi, 
                      const DenseMatrix<double>& E,
                      size_t const numberOfElectrons,
                      size_t const M_disc)
    {
        assert(numberOfElectrons==2);

        DenseMatrix<double> U = Phi.matrix();
        DenseMatrix<double> V = Psi.matrix();

        // On calcule directement
        
        Vector<double> Vcol0(V.numberOfRows());
        Vector<double> Vcol1(V.numberOfRows());
        Vector<double> Ucol0(V.numberOfRows());
        Vector<double> Ucol1(V.numberOfRows());
        
        Ucol0.wrapped_ = U.wrapped_.col(0);
        Ucol1.wrapped_ = U.wrapped_.col(1);
        Vcol0.wrapped_ = V.wrapped_.col(0);
        Vcol1.wrapped_ = V.wrapped_.col(1);

        return 0.5 * ( elint(M_disc, Ucol0, Vcol0, Ucol1, Vcol1, E) 
                     - elint(M_disc, Ucol1, Vcol0, Ucol0, Vcol1, E)
                     - elint(M_disc, Ucol0, Vcol1, Ucol1, Vcol0, E) 
                     + elint(M_disc, Ucol1, Vcol1, Ucol0, Vcol0, E) );

    }
   
    double H2_slat(const SlaterDeterminant& Phi, 
                   const SlaterDeterminant& Psi, 
                   const DenseMatrix<double>& O, 
                   const DenseMatrix<double>& E,
                   size_t const numberOfElectrons,
                   size_t const M_disc_,
                   double ratio_ = 1e-12)
    {
        if (numberOfElectrons==2) 
            return H2_slat_N2(Phi,Psi,E,numberOfElectrons,M_disc_ );

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
                Sinv.wrapped_ = S.wrapped_.inverse();
                double scal0 = 0;

                for (size_t i=0; i < numberOfElectrons; i++)
                {
                    Vector<double> Ucol_i(U.numberOfRows());
                    Ucol_i.wrapped_ = U.wrapped_.col(i);

                    for (size_t j=0; j < numberOfElectrons; j++)
                    {
                        Vector<double> Ucol_j(U.numberOfRows());
                        Ucol_j.wrapped_ = U.wrapped_.col(j);

                        for (size_t k = 0; k < numberOfElectrons; k++)
                        {
                            Vector<double> Vcol_k(U.numberOfRows());
                            Vcol_k.wrapped_ = V.wrapped_.col(k);

                            for (size_t l=0; l < numberOfElectrons; l++)
                            {
                                Vector<double> Vcol_l(U.numberOfRows());
                                Vcol_l.wrapped_ = V.wrapped_.col(l);

                                 //Ici je prends les notations de Friedrichs
                                scal0 += 0.5 * elint(M_disc_, Ucol_i, Vcol_k, Ucol_j, Vcol_l, E)
                                             * (Sinv(k,i)*Sinv(l,j) - Sinv(k,j)*Sinv(l,i));

                            }
                        }
                    }
                }

                return scal0 * (S.wrapped_.determinant());
            }
            else
            {
                Eigen::JacobiSVD<eigen<double>::DenseMatrixType> svd(S.wrapped_, Eigen::ComputeThinU | Eigen::ComputeThinV);

                DenseMatrix<double> Uvec = svd.matrixU();
                DenseMatrix<double> Vvec = svd.matrixV();

                Vector<double> D = svd.singularValues();
                double lmax = D.max();

                D.wrapped_ *= (1.0/lmax);

                //On compte la multiplicité de la valeur propre nulle
                int mult = 0;
                for (size_t i= 0; i< numberOfElectrons; i++)
                {
                    if (D(i) < ratio_) 
                        ++mult;
                }


                //cout << "mult = " << mult << endl;

                //si la multiplicité de 0 est plus grande que 2, la valeur est 0
                if (mult >2.5) 
                    return 0;

                else
                {
                    if ((1.5 > mult) && (mult > 0.5)) //mult = 1
                    {

                        int indmin = getIndMin(D);
                        Vector<double> xV(Vvec.numberOfRows());
                        xV.wrapped_ = Vvec.wrapped_.col(indmin);
                        Vector<double> xU(Uvec.numberOfRows());
                        xU.wrapped_ = Uvec.wrapped_.col(indmin);

                        Vector<double> orthU(U.numberOfRows());
                        orthU.wrapped_ = U.wrapped_ * xU.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                        Vector<double> orthV(V.numberOfRows());
                        orthV.wrapped_ = V.wrapped_ * xV.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par U

                        Vector<double> temp(O.numberOfColumns());
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * orthU.wrapped_;
                        orthU.wrapped_ = 1.0 / sqrt( (orthU.wrapped_.adjoint() * temp.wrapped_) ) * orthU.wrapped_;
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * orthV.wrapped_;
                        orthV.wrapped_ = 1.0 / sqrt( (orthV.wrapped_.adjoint()) *temp.wrapped_) *orthV.wrapped_;


                        DenseMatrix<double> PsiU(M_disc_,numberOfElectrons);
                        PsiU.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: le deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                        DenseMatrix<double> PsiV(M_disc_,numberOfElectrons);
                        PsiV.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_,numberOfElectrons);

                        Vector<double> muU(numberOfElectrons);
                        muU.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);
                        Vector<double> muV(numberOfElectrons);
                        muV.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);

                        for(size_t k=0; k<numberOfElectrons; k++)
                        {
                            temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (U.wrapped_.col(k));
                            double PS = orthU.wrapped_.adjoint() * temp.wrapped_;
                            PsiU.wrapped_.col(k) = U.wrapped_.col(k) - PS * orthU.wrapped_;
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

                        DenseMatrix<double> coeffs_bas(U.numberOfRows(), xUbas.numberOfColumns());
                        coeffs_bas.wrapped_ = U.wrapped_ * xUbas.wrapped_;
                        //Puis on orthonormalise les vecteurs de coeffs_bas
                        for (size_t i=0; i<(numberOfElectrons-1); i++)
                        {
                            Vector<double> temp2(O.numberOfColumns());
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
                            PS2 = ((coeffs_bas.wrapped_.col(i)).adjoint()) * temp2.wrapped_;
                            coeffs_bas.wrapped_.col(i) = (1.0/sqrt(PS2)) * coeffs_bas.wrapped_.col(i);

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
                                double PS = ((PsiU.wrapped_.col(k)).adjoint()) * temp.wrapped_;
                                CU(l,k) = PS;

                                temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * (coeffs_bas.wrapped_.col(l));
                                PS = ((PsiV.wrapped_.col(k)).adjoint())*temp.wrapped_;
                                CV(l,k) = PS;
                            }
                        }

                        //On construit la somme
                        double sum2 = 0;
                        for (size_t l=0; l< numberOfElectrons-1; l++)
                            sum2 += elint(M_disc_,coeffs_bas.column(l),coeffs_bas.column(l),orthU,orthV,E) 
                                  - elint(M_disc_,coeffs_bas.column(l),orthV,coeffs_bas.column(l),orthU,E);


                        //On a alors pour tout k U(:,k) = muU(k)*orthU + PsiU(:,k)
                        double sum = 0;

                        for (size_t k1= 0; k1 < numberOfElectrons; k1++)
                        {
                            for (size_t k2 = 0; k2 < numberOfElectrons; k2++)
                            {
                                DenseMatrix<double> aU(numberOfElectrons-1,numberOfElectrons-1);
                                aU.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1,numberOfElectrons-1);
                                aU.wrapped_.block(0,0,numberOfElectrons-1,k1) = CU.wrapped_.block(0,0,numberOfElectrons-1,k1);
                                aU.wrapped_.block(0,k1,numberOfElectrons-1,numberOfElectrons-1-k1) = CU.wrapped_.block(0,(k1+1),numberOfElectrons-1, numberOfElectrons-1-k1);

                                DenseMatrix<double> aV(numberOfElectrons-1,numberOfElectrons-1);
                                aV.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-1,numberOfElectrons-1);
                                aV.wrapped_.block(0,0,numberOfElectrons-1,k2) = CV.wrapped_.block(0,0,numberOfElectrons-1,k2);
                                aV.wrapped_.block(0,k2,numberOfElectrons-1, numberOfElectrons-1-k2) = CV.wrapped_.block(0,k2+1,numberOfElectrons-1, numberOfElectrons-1-k2);

                                sum += muU(k1)*muV(k2)*pow(-1,k1)*pow(-1,k2)*sum2*(aU.wrapped_.determinant())*(aV.wrapped_.determinant());
                            }

                        }

                        return sum;
                    }

                    else
                    {	//mult = 2
                        //La dimension de l'intersection des deux espaces est N-2: il
                        //faut récupérer pour U et pour V les deux vecteurs qui sont
                        //orthonormaux à l'espace

                        size_t indmin = getIndMin(D);
                        Vector<double> xV1(Vvec.numberOfRows());
                        xV1.wrapped_ = Vvec.wrapped_.col(indmin);
                        Vector<double> xU1(Uvec.numberOfRows());
                        xU1.wrapped_ = Uvec.wrapped_.col(indmin);

                        D(indmin) = 1e20;

                        size_t indmin2 = getIndMin(D);
                        Vector<double> xV2 = Vvec.column(indmin2);
                        Vector<double> xU2 = Uvec.column(indmin2);

                        Vector<double> orthU1(U.numberOfRows());
                        orthU1.wrapped_ = U.wrapped_ * xU1.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par V
                        Vector<double> orthU2(U.numberOfRows());
                        orthU1.wrapped_ = U.wrapped_ * xU2.wrapped_;
                        Vector<double> orthV1(V.numberOfRows());
                        orthV1.wrapped_ = V.wrapped_ * xV1.wrapped_;  //Vecteur dans l'orthogonal du sous-espace engendré par U
                        Vector<double> orthV2(V.numberOfRows());
                        orthV2.wrapped_ = V.wrapped_ * xV2.wrapped_;

                        Vector<double> temp(O.numberOfRows());
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * orthU1.wrapped_;
                        orthU1.wrapped_ = (1.0/sqrt((orthU1.wrapped_.adjoint()) * temp.wrapped_)) * orthU1.wrapped_;

                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * orthU2.wrapped_;
                        double PS = (orthU1.wrapped_.adjoint()) * temp.wrapped_;
                        orthU2.wrapped_ -= PS * orthU1.wrapped_;
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>()) * orthU2.wrapped_;
                        orthU2.wrapped_ = (1.0/sqrt((orthU2.wrapped_.adjoint())*temp.wrapped_))*orthU2.wrapped_;

                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*orthV1.wrapped_;
                        orthV1.wrapped_ = (1.0/sqrt((orthV1.wrapped_.adjoint())*temp.wrapped_))*orthV1.wrapped_;

                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*orthV2.wrapped_;
                        PS = (orthV1.wrapped_.adjoint())*temp.wrapped_;
                        orthV2.wrapped_ -= PS*orthV1.wrapped_;
                        temp.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*orthV2.wrapped_;
                        orthV2.wrapped_ = (1.0/sqrt((orthV2.wrapped_.adjoint())*temp.wrapped_))*orthV2.wrapped_;

                        DenseMatrix<double> PsiU(M_disc_,numberOfElectrons);
                        PsiU.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_,numberOfElectrons); //Le reste des fonctions: les deux sous-espaces engendrés sont les mêmes, égaux à l'intersection de deux sous-espaces de départ
                        DenseMatrix<double> PsiV(M_disc_,numberOfElectrons);
                        PsiV.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc_,numberOfElectrons);

                        Vector<double> muU1(numberOfElectrons); 
                        muU1.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);
                        Vector<double> muV1(numberOfElectrons); 
                        muV1.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);
                        Vector<double> muU2(numberOfElectrons); 
                        muU2.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);
                        Vector<double> muV2(numberOfElectrons); 
                        muV2.wrapped_ = eigen<double>::VectorType::Zero(numberOfElectrons);

                        for (size_t k=0; k< numberOfElectrons; k++)
                        {

                            Vector<double> temp2(O.numberOfRows());
                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(U.wrapped_.col(k));
                            muU1(k) = (orthU1.wrapped_.adjoint())*temp2.wrapped_;
                            muU2(k) = (orthU2.wrapped_.adjoint())*temp2.wrapped_;
                            temp2.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(V.wrapped_.col(k));
                            muV1(k) = orthV1.wrapped_.adjoint()*temp2.wrapped_;
                            muV2(k) = orthV2.wrapped_.adjoint()*temp2.wrapped_;


                            PsiU.wrapped_.col(k) = U.wrapped_.col(k) - muU1(k)*orthU1.wrapped_ - muU2(k)*orthU2.wrapped_;
                            PsiV.wrapped_.col(k) = V.wrapped_.col(k) - muV1(k)*orthV1.wrapped_ - muV2(k)*orthV2.wrapped_;
                        }

                        //On crée une base orthonormale de l'intersection des deux
                        //sous-espaces: par exemple à partir de xU
                        DenseMatrix<double> xUbas(numberOfElectrons, numberOfElectrons-2);
                        xUbas.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons, numberOfElectrons-2);
                        int ind = 0;
                        for (size_t k=0; k< numberOfElectrons; k++)
                        {
                            if ( (k!=indmin) && (k!=indmin2))
                            {
                                xUbas.wrapped_.col(ind) = Uvec.wrapped_.col(k);
                                ind = ind+1;
                            }
                        }

                        DenseMatrix<double> coeffs_bas(U.numberOfRows(), xUbas.numberOfColumns());
                        coeffs_bas.wrapped_ = U.wrapped_ * xUbas.wrapped_;

                        //Puis on orthonormalise les vecteurs de coeffs_bas
                        for (size_t i= 0; i < numberOfElectrons-2; i++)
                        {
                            //On commence par orthonormaliser par rapport à orthU1 et orthU2
                            Vector<double> temp3(O.numberOfRows());
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            double PS3 = (orthU1.wrapped_.adjoint())*temp3.wrapped_;
                            coeffs_bas.wrapped_.col(i) = coeffs_bas.wrapped_.col(i) - PS3*orthU1.wrapped_;
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            PS3 = (orthU2.wrapped_.adjoint())*temp3.wrapped_;
                            coeffs_bas.wrapped_.col(i) = coeffs_bas.wrapped_.col(i) - PS3*orthU2.wrapped_;
                            for (size_t j=0; j<i; j++)
                            {
                                temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                                PS3 = ((coeffs_bas.wrapped_.col(j)).adjoint())*temp3.wrapped_;
                                coeffs_bas.wrapped_.col(i) = coeffs_bas.wrapped_.col(i) - PS3*coeffs_bas.wrapped_.col(j);
                            }
                            temp3.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*(coeffs_bas.wrapped_.col(i));
                            coeffs_bas.wrapped_.col(i) = 1.0/sqrt((coeffs_bas.wrapped_.col(i).adjoint())*temp3.wrapped_)*coeffs_bas.wrapped_.col(i);
                        }


                        DenseMatrix<double> CU(numberOfElectrons-2,numberOfElectrons); 
                        CU.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-2,numberOfElectrons);
                        DenseMatrix<double> CV(numberOfElectrons-2,numberOfElectrons); 
                        CV.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-2,numberOfElectrons);

                        for (size_t k=0; k< numberOfElectrons; k++)
                        {
                            for (size_t l=0; l< numberOfElectrons-2; l++)
                            {
                                Vector<double> temp4(O.numberOfRows());
                                temp4.wrapped_ = (O.wrapped_.selfadjointView<Eigen::Upper>())*coeffs_bas.wrapped_.col(l);
                                CU(l,k) = (PsiU.wrapped_.col(k)).adjoint()*temp4.wrapped_;
                                CV(l,k) = (PsiV.wrapped_.col(k)).adjoint()*temp4.wrapped_;
                            }
                        }

                        //On construit la somme
                        double sum2 = elint(M_disc_, orthU1,orthV1,orthU2,orthV2,E) - elint(M_disc_, orthU1,orthV2,orthU2,orthV1, E);

                        //On a alors pour tout k U(:,k) = muU1(k)*orthU1 + muU2(k)*orthU2 + PsiU(:,k)
                        double sum = 0;

                        for (size_t iu = 0; iu < numberOfElectrons-1; iu++)
                        {
                            for (size_t ju = iu+1 ; ju <numberOfElectrons; ju++)
                            {
                                DenseMatrix<double> aU(numberOfElectrons-2,numberOfElectrons-2);
                                aU.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-2,numberOfElectrons-2);
                                aU.wrapped_.block(0,0,numberOfElectrons-2, iu) = CU.wrapped_.block(0,0,numberOfElectrons-2,iu);
                                aU.wrapped_.block(0,iu,numberOfElectrons-2, ju-iu-1) = CU.wrapped_.block(0,iu+1, numberOfElectrons-2, ju-iu-1);
                                aU.wrapped_.block(0,ju-1,numberOfElectrons-2, numberOfElectrons-ju-1) = CU.wrapped_.block(0,ju+1,numberOfElectrons-2, numberOfElectrons-ju-1);


                                for (size_t iv = 0; iv < numberOfElectrons-1; iv++)
                                {
                                    for (size_t jv = iv+1; jv < numberOfElectrons; jv++)
                                    {

                                        DenseMatrix<double> aV(numberOfElectrons-2,numberOfElectrons-2);
                                        aV.wrapped_ = eigen<double>::DenseMatrixType::Zero(numberOfElectrons-2,numberOfElectrons-2);
                                        aV.wrapped_.block(0,0,numberOfElectrons-2,iv) = CV.wrapped_.block(0,0, numberOfElectrons-2, iv);
                                        aV.wrapped_.block(0,iv,numberOfElectrons-2,jv-iv-1) = CV.wrapped_.block(0,iv+1,numberOfElectrons-2,jv-iv-1);
                                        aV.wrapped_.block(0,jv-1,numberOfElectrons-2,numberOfElectrons-jv-1) = CV.wrapped_.block(0,jv+1,numberOfElectrons-2,numberOfElectrons-jv-1);


                                        sum += pow(-1,iu) * pow(-1,ju+1) *pow(-1,iv) *pow(-1,jv+1) 
                                             * ((muU1(iu)*muU2(ju) -muU1(ju)*muU2(iu)) * (muV1(iv)*muV2(jv) -muV1(jv)*muV2(iv))) 
                                             * sum2*(aU.wrapped_.determinant())*(aV.wrapped_.determinant());

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
    {
        DenseMatrix<double> S = Smat(Phi, Psi, O);
        return S.wrapped_.determinant();
    }
    
     


    SlaterDeterminant
    hartree_fock(std::size_t const M_disc,
                 std::size_t const numberOfElectrons,
                 DenseMatrix<double> const & E,
                 DenseMatrix<double> const & K,
                 DenseMatrix<double> const & O,
                 DenseMatrix<double> const & Nu,
                 SlaterDeterminant const & initial_solution,
                 std::size_t const numberOfIterations)
    {
        //A améliorer: pour l'instant j'utilise des matrices pleines car il n'y a pas de solveur eigenvalue pour les matrices sparses

        DenseMatrix<double> Phi0 = initial_solution.matrix();
        DenseMatrix<double> F0(K.numberOfRows(), K.numberOfColumns());
        F0.wrapped_ = K.wrapped_ + Nu.wrapped_;
        
        //Roothan parce que c'est le plus simple: ToDo coder ODA
        for( std::size_t iteration = 0; iteration < numberOfIterations; ++iteration )
        {
            DenseMatrix<double> F(F0.numberOfRows(), F0.numberOfColumns());
            DenseMatrix<double> Fock(F0.numberOfRows(), F0.numberOfColumns());
            Fock = FockMat(M_disc, Phi0, E);
            F.wrapped_ = F0.wrapped_ + Fock.wrapped_;


            Eigen::GeneralizedSelfAdjointEigenSolver<eigen<double>::DenseMatrixType> es(F.wrapped_, O.wrapped_, Eigen::ComputeEigenvectors|Eigen::Ax_lBx);

            Vector<double> D = es.eigenvalues();
            DenseMatrix<double> V = es.eigenvectors();

            std::vector<size_t> Itab = getIndMin(D, numberOfElectrons);

            DenseMatrix<double> Phinew(M_disc, numberOfElectrons);
            Phinew.wrapped_ = eigen<double>::DenseMatrixType::Zero(M_disc, numberOfElectrons);
            for (size_t i=0; i< numberOfElectrons; i++)
                Phinew.wrapped_.col(i) = V.wrapped_.col(Itab[i]);
            Phi0 = Phinew;


            double lambda = 0;
            for (std::size_t i = 0; i < numberOfElectrons; ++i)
            {
                lambda += D(Itab[i]);
                for (std::size_t j = 0; j< numberOfElectrons; ++j)
                    lambda += 0.5 * ( elint(M_disc, Phi0.column(i), Phi0.column(i), Phi0.column(j), Phi0.column(j), E) 
                                    - elint(M_disc, Phi0.column(i), Phi0.column(j), Phi0.column(j), Phi0.column(i), E) );
            }

            SlaterDeterminant sol(Phi0);
            DenseMatrix<double> sum(K.numberOfRows(), K.numberOfColumns());
            sum.wrapped_ = K.wrapped_ + Nu.wrapped_;
            double lambda2 = H1_slat(Phi0, Phi0, O, sum, numberOfElectrons, M_disc) + H2_slat(Phi0, Phi0, O, E, numberOfElectrons, M_disc);
            lambda2 /= over_slat(Phi0, Phi0, O);


        }

        SlaterDeterminant sol(Phi0);

        return sol;
    }

}


#endif	

