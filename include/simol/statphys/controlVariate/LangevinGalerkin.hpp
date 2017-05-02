#ifndef SIMOL_LANGEVINGALERKIN_HPP
#define SIMOL_LANGEVINGALERKIN_HPP

#include "Galerkin.hpp"

namespace simol
{
  class LangevinGalerkin : public Galerkin
  {
    public:
      LangevinGalerkin(Input const& input);
      virtual void createOperators();
  
      //virtual void computeExpToTrigTens();
      virtual void createLthm0();
      virtual void createLthm();
      virtual void compute() = 0;
      virtual DVec CVcoeffsVec() const;
      virtual DVec CVObservable() const = 0;
      virtual DVec solveResilient(const SMat& A, const DVec& Y) const;
  };
  
  class PeriodicLangevinGalerkin : public LangevinGalerkin
  {
  public:
    PeriodicLangevinGalerkin(Input const& input);

    //maxOfFourier_((nbOfQModes_ + 1) / 2),  //ex : 3
    virtual void compute();
    virtual DVec CVObservable() const;
  };
  
  class ColloidLangevinGalerkin : public LangevinGalerkin
  {
  public:
    ColloidLangevinGalerkin(Input const& input);

    virtual void compute();
    virtual DVec CVObservable() const;
  };

  /*class BoundaryLangevinGalerkin : public Galerkin
  {
      //SMat SId_;
    public:
      BoundaryLangevinGalerkin(Input const& input);
      int iTens(int iOfFourier2, int iOfHermite, int iOfParticle) const;
      SMat PMatToTens(SMat const& PMat, int iOfParticleP);
      SMat doubleMatToTens(SMat const& QMat, SMat const& PMat, int iOfParticleQ, int iOfParticleP);
      //virtual void computeExpToTrigTens();
      void createLham();
      virtual void createLthm0();
      virtual void createLthm();

      virtual void compute();
      DVec CVcoeffsVec() const {throw runtime_error("CVcoeffs not implemented for the chain !");}
  };*/


}

#endif