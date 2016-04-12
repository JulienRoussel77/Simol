#include "Output.hpp"

using std::cout; 
using std::endl; 
using std::vector;

namespace simol{

  void Output::setControlVariates(Input& input, Potential& potential, Galerkin* galerkin)
  {
    velocityCV_ = createControlVariate(input, potential, galerkin);
    forceCV_ = createControlVariate(input, potential, galerkin);
    lengthCV_ = createControlVariate(input, potential, galerkin);
    midFlowCV_ = createControlVariate(input, potential, galerkin);
    sumFlowCV_ = createControlVariate(input, potential, galerkin);
  }

  void Output::displayGeneratorOnBasis(ofstream& out, vector<Particle> const& configuration, ControlVariate& controlVariate, double time)
  {
    out << time << " " << modulo(configuration[0].position(0), -M_PI, M_PI) << " " << configuration[0].momentum(0) << " " << controlVariate.lastGeneratorOnBasis()(0) << " " << controlVariate.basisFunction(configuration) << endl;
  }
  
}
