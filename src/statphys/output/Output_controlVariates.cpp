#include "simol/statphys/output/Output.hpp"

using std::cout;
using std::endl;
using std::vector;

namespace simol
{

  /*void Output::setControlVariates(Input& input, Potential& potential, Galerkin*)
  {
    if (input.controlVariateName() == "ExpFourierHermite")
      cvBasis_.basis_ = new ExpFourierHermiteBasis(input, potential);
  }*/
  
  /*Observable* Output::addControlVariate(const Input& input, const string& outPath, Galerkin* galerkin)
  {
    Observable* obsPtr = createControlVariate(input, outPath, cvBasis_, galerkin);
    observables_.push_back(obsPtr);
    return obsPtr;
  }*/
  
  /*void Output::makeCVObservables(const Input& input, Galerkin* galerkin)
  {
    if (input.CVObservable() == "kineticEnergy")
    for (iOfObservable = 0; iOfObservable < nbOfObservables(); iOfObservable++)
    {
      
    }*/

  /*void Output::displayGeneratorOnBasis(ofstream& out ControlVariate& controlVariate, double time)
  {
    out << time << " " << modulo(configuration[0].position(0), -M_PI, M_PI) << " " << configuration[0].momentum(0) << " " << controlVariate.lastGeneratorOnBasis()(0) << " " << controlVariate.basisFunction(configuration) << endl;
  }*/

}
