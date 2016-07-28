#include "simol/statphys/output/Output.hpp"

using std::cout;
using std::endl;
using std::vector;

namespace simol
{

  void Output::setControlVariates(Input& input, Potential& potential, Galerkin* galerkin)
  {
    /*//velocityCV_ = createControlVariate(input, potential, galerkin);
    obsForce_ = createControlVariate(input, potential, galerkin);
    obsLength_ = createControlVariate(input, potential, galerkin);
    midFlowCV_ = createControlVariate(input, potential, galerkin);
    sumFlowCV_ = createControlVariate(input, potential, galerkin);
    modiFlowCV_ = createControlVariate(input, potential, galerkin);*/
  }
  
  Observable* Output::addControlVariate(const Input& input, const string& outPath, Potential& potential, Galerkin* galerkin)
  {
    Observable* obsPtr = createControlVariate(input, outPath, potential, galerkin);
    observables_.push_back(obsPtr);
    return obsPtr;
  }

  /*void Output::displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time)
  {
    out << time << " " << modulo(configuration[0].position(0), -M_PI, M_PI) << " " << configuration[0].momentum(0) << " " << controlVariate.lastGeneratorOnBasis()(0) << " " << controlVariate.basisFunction(configuration) << endl;
  }*/

}
