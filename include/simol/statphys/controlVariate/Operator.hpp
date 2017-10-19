#ifndef SIMOL_OPERATOR_HPP
#define SIMOL_OPERATOR_HPP

#include "simol/statphys/potential/Potential.hpp"
//#include "simol/statphys/system/Particle.hpp"
#include "simol/statphys/input/Input.hpp"
//#include "simol/statphys/output/Output.hpp"
//#include "simol/core/random/RNG.hpp"
//#include "simol/statphys/controlVariate/ControlVariate.hpp"
//#include "simol/statphys/controlVariate/CVBasis.hpp"
# include <iostream>
//#include "simol/statphys/controlVariate/AllGalerkins.hpp"

#include "simol/statphys/system/System.hpp"
#include "simol/statphys/dynamics/DynamicsParameters.hpp"
#include "simol/statphys/controlVariate/Basis.hpp"

namespace simol
{
  class Operator;
  Operator* createCVOperator(Input const& input);
  
  class Operator
  {
  public:
    Operator(Input const& input);
    //virtual void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const = 0;
    virtual DVec value(TensorBasis* basis, System const& syst) const = 0;
    virtual DMat qVariable(const System& syst) const;
    virtual DMat pVariable(const System& syst) const;
    virtual DMat forces(const System& syst) const;
    virtual DVec basisVariables(const System& syst) const;
    
    DynamicsParameters parameters_;
  };
  
  class HamiltonianGenerator : public Operator
  {
  public:
    HamiltonianGenerator(Input const& input);
    //void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const;
    DVec value(TensorBasis* basis, System const& syst) const;
  };
  
  class OverdampedGenerator : public Operator
  {
  public:
    OverdampedGenerator(Input const& input);
    //void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const;
    DVec value(TensorBasis* basis, System const& syst) const;
  };
  
  class OvdColloidGenerator : public Operator
  {
  public:
    OvdColloidGenerator(Input const& input);
    //void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const;
    DVec value(TensorBasis* basis, System const& syst) const;
    virtual DMat qVariable(const System& syst) const;
    virtual DMat pVariable(const System& syst) const;
    virtual DMat forces(const System& syst) const;
    virtual DVec basisVariables(const System& syst) const;
  };
  
  class LangevinGenerator : public Operator
  {
  public:
    LangevinGenerator(Input const& input);
    //void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const;
    DVec value(TensorBasis* basis, System const& syst) const;
  };
  
  class UnderdampedGenerator : public Operator
  {
  public:
    UnderdampedGenerator(Input const& input);
    //void value(shared_ptr<CVBasis> cvBasis, System const& syst, shared_ptr<DynamicsParameters> dynaPara) const;
    DVec value(TensorBasis* basis, System const& syst) const;
    virtual double uVariable(const System& syst) const;
    virtual DVec basisVariables(const System& syst) const;
  };
}

#endif



