#include "simol/statphys/controlVariate/Operator.hpp"

using std::cout;
using std::endl;

namespace simol
{
  Operator* createCVOperator(Input const& input)
  {
    cout << "input.CVOperator() = " << input.CVOperator() << endl;
    if (input.CVOperator() == "None")
      return nullptr;
    if (input.CVOperator() == "Hamiltonian")
      return new HamiltonianGenerator(input);
    if (input.CVOperator() == "Overdamped")
      return new OverdampedGenerator(input);
    if (input.CVOperator() == "Langevin")
      return new LangevinGenerator(input);
    if (input.CVOperator() == "OvdColloid")
      return new OvdColloidGenerator(input);
    if (input.CVOperator() == "Underdamped")
      return new UnderdampedGenerator(input);

    else
      throw runtime_error(input.CVOperator() + " is not a valid operator !");
    return 0; 

  }
  
  Operator::Operator(Input const& input):
    parameters_(input)
  {}
  
  DMat Operator::qVariable(System const& syst) const
  {
    DMat qVariableMat= DVec::Zero(1);
    qVariableMat(0,0) = syst(0).position(0);
    return qVariableMat;
  }
  
  DMat Operator::pVariable(System const& syst) const
  {    
    DMat qVariableMat = DVec::Zero(1);
    qVariableMat(0,0) = syst(0).momentum(0);
    return qVariableMat;
  }
  
  DMat Operator::forces(System const& syst) const
  {    
    DMat forceMat = DVec::Zero(1);
    forceMat(0,0) = syst(0).force(0);
    return forceMat;
  }
  
  DVec Operator::basisVariables(System const& syst) const
  {
    DVec variables = DVec::Zero(2);
    variables(0) = syst(0).position(0);
    variables(1) = syst(0).momentum(0);
    return variables;
  }
  
  
  HamiltonianGenerator::HamiltonianGenerator(Input const& input):
    Operator(input)
  {}
  
  DVec HamiltonianGenerator::value(TensorBasis* basis, System const& syst) const
  {
    //if (variables.rows() != 2)
    //  std::throw runtime_error("Call to HamiltonianGenerator::value takes two values!");
    //double qVariable = variables[0];
    //double pVariable = variables[1];
    DVec variables = basisVariables(syst);
    DVec generatorOnBasisValues = DVec::Zero(basis->totalNbOfElts());
    for (int iOfFunction = 0; iOfFunction < basis->totalNbOfElts(); iOfFunction++)
      //generatorOnBasisValues(iOfFunction) += dot( syst.momenta() , basis->gradientQ(syst, iOfFunction))
      //                        + dot( syst.forces() , basis->gradientP(syst, iOfFunction));
      generatorOnBasisValues(iOfFunction) += dot( pVariable(syst) , basis->gradientQ(variables, iOfFunction))
                              + dot( forces(syst) , basis->gradientP(variables, iOfFunction));
    return generatorOnBasisValues;
  }
  
  
  
  OverdampedGenerator::OverdampedGenerator(Input const& input):
    Operator(input)
  {}
  
  DVec OverdampedGenerator::value(TensorBasis* basis, System const& syst) const
  {
    DVec generatorOnBasisValues = DVec::Zero(basis->totalNbOfElts());
    DVec variables = basisVariables(syst);
    for (int iOfFunction = 0; iOfFunction < basis->totalNbOfElts(); iOfFunction++)
    {
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      generatorOnBasisValues(iOfFunction) += basis->laplacianQ(variables, iOfFunction) / parameters_.beta()
                               + dot( forces(syst), basis->gradientQ(variables, iOfFunction));
    }
    return generatorOnBasisValues;
  }
  
  OvdColloidGenerator::OvdColloidGenerator(Input const& input):
    Operator(input)
  {}
  
  DVec OvdColloidGenerator::basisVariables(System const& syst) const
  {
    DVec variables = DVec::Zero(2);
    DVec r01 = syst(1).position() - syst(0).position();
    //DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    DVec p01 = syst(1).momentum() - syst(0).momentum();
    variables(0) = d01;
    variables(1) = dot(p01, r01)/d01;
    if (!std::isfinite(d01))
      throw runtime_error("The colloid distance is not valid : d = " + std::to_string(d01) + " !");
    return variables;
  }
  
  DMat OvdColloidGenerator::qVariable(System const& syst) const
  {
    DMat distanceMat = DVec::Zero(1);
    DVec r01 = syst(1).position() - syst(0).position();
    //DVec r01 = syst.periodicDistance(r01np);
    distanceMat(0,0) = r01.norm();
    return distanceMat;
  }
  
  DMat OvdColloidGenerator::pVariable(System const& syst) const
  {
    DMat momentumMat = DVec::Zero(1);
    DVec r01 = syst(1).position() - syst(0).position();
    //DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    momentumMat(0,0) = dot(syst(1).momentum() - syst(0).momentum(), r01) / d01;
    return momentumMat;
  }
  
  DMat OvdColloidGenerator::forces(System const& syst) const
  {
    DMat forcesMat = DVec::Zero(1);
    DVec r01 = syst(1).position() - syst(0).position();
    //DVec r01 = syst.periodicDistance(r01np);
    double d01 = r01.norm();
    forcesMat(0,0) = dot(syst(1).force() - syst(0).force(), r01) / (2*d01) + (syst.dimension()-1) / d01;
    //cout << r01.adjoint() << " -> "  << (syst(1).force() - syst(0).force()).adjoint() << endl;
    //cout  << d01 << " -> " << (syst(1).force() - syst(0).force()).norm() << endl;
    
    //ofstream test("test", std::ofstream::app);
    //test << d01 << " " << dot(syst(1).force() - syst(0).force(), r01) / (2*d01) << " " << (syst.dimension()-1) / d01 << " " << dot(syst(1).force(), r01) / (2*d01) << " " << dot(-syst(0).force(), r01) / (2*d01) << endl;
    
    //cout << "CVB force : " << d01 << " " << dot(syst(1).force() - syst(0).force(), r01) / (2*d01) << endl;
    return forcesMat;
  }
    

  
  DVec OvdColloidGenerator::value(TensorBasis* basis, System const& syst) const
  {
    DVec generatorOnBasisValues = DVec::Zero(basis->totalNbOfElts());
    DVec variables = basisVariables(syst);
    for (int iOfFunction = 0; iOfFunction < basis->totalNbOfElts(); iOfFunction++)
    {
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      generatorOnBasisValues(iOfFunction) += basis->laplacianQ(variables, iOfFunction) / parameters_.beta()
                               + dot( forces(syst), basis->gradientQ(variables, iOfFunction));
    }
    return generatorOnBasisValues;
  }
  
  LangevinGenerator::LangevinGenerator(Input const& input):
    Operator(input)
  {}
  
  DVec LangevinGenerator::value(TensorBasis* basis, System const& syst) const
  {
    DVec generatorOnBasisValues = DVec::Zero(basis->totalNbOfElts());
    DVec variables = basisVariables(syst);
    for (int iOfFunction = 0; iOfFunction < basis->totalNbOfElts(); iOfFunction++)
    {
      //for (int iOfParticle = 0; iOfParticle < syst.nbOfParticles(); iOfParticle++)
      generatorOnBasisValues(iOfFunction) += dot( pVariable(syst) , basis->gradientQ(variables, iOfFunction))
                              + dot( forces(syst) , basis->gradientP(variables, iOfFunction))
                              - parameters_.gamma() / parameters_.mass() * dot( pVariable(syst) , basis->gradientP(variables, iOfFunction))
                               + parameters_.gamma() / parameters_.beta() * basis->laplacianP(variables, iOfFunction);
    }
    return generatorOnBasisValues;
  }
  
  UnderdampedGenerator::UnderdampedGenerator(Input const& input):
    Operator(input)
  {}
  
  double UnderdampedGenerator::uVariable(System const& syst) const
  {
    DVec projPosition = DVec::Zero(2);
    projPosition(0) = syst(0).position(0);
    //cout << "uVar = " << syst.externalPotential().value(projPosition) << " + " << pow(syst(0).momentum(0), 2)/2 << endl;
    //return syst(0).energy();
    return syst.externalPotential().value(projPosition) + pow(syst(0).momentum(0), 2)/2;
  }
  
  DVec UnderdampedGenerator::basisVariables(System const& syst) const
  {
    DVec variables = DVec::Zero(2);
    variables(0) = uVariable(syst);
    variables(1) = syst(0).momentum(0);
    return variables;
  }
  
  ///
  /// For simplicity we use the function gradientQ instead of a "gradientU" to avoid writing more code...
  DVec UnderdampedGenerator::value(TensorBasis* basis, System const& syst) const
  {
    DVec generatorOnBasisValues = DVec::Zero(basis->totalNbOfElts());
    DVec variables = basisVariables(syst);
    double pVar = pVariable(syst)(0,0);
    for (int iOfFunction = 0; iOfFunction < basis->totalNbOfElts(); iOfFunction++)
      generatorOnBasisValues(iOfFunction) += (pVar > 0? 1:-1)*( (1./parameters_.beta() - pow(pVar, 2) + pVar * parameters_.nonEqAmplitude()/parameters_.gamma()) * basis->gradientQ(variables, iOfFunction)(0,0)
                                              + 1./parameters_.beta() * pow(pVar, 2) * basis->laplacianQ(variables, iOfFunction));
    return generatorOnBasisValues;
  }
  
}