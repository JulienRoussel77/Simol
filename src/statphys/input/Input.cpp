#include "simol/statphys/input/Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::max;
using std::min;

namespace simol
{
  bool sameLetters(const string& x, const string& y)
  {
    if (x.size() != y.size()) return false;
    for (int i=0; i < (int)x.size(); i++)
      if (tolower(x[i]) != tolower(y[i])) return false;
    return true;
  }

  ///
  ///Reads the next element in the file read
  double readItem(ifstream& in)
  {
    double item;
    in >> item;
    return item;
  }

  ///Reads the input file using the YAML library
  ///Reads the initial conditions in the "settings" file, if indicated
  Input::Input(CommandLine cmd):
    data(YAML::LoadFile(cmd.inputFileName())),
    inputPath_(cmd.inputFileName()),
    inputFlux_(inputPath()),
    initialPositions_(nbOfParticles(), Vector<double>(dimension())),
    initialMomenta_(nbOfParticles(), Vector<double>(dimension())),
    initialInternalEnergies_(nbOfParticles(), 0)
  {
    assert(data.IsDefined());
    if (doRestart())
    {
      cout << " - reading from restart file " << restartFileName() << endl;
      ifstream initialConditions(restartFileName());
      if (!initialConditions.is_open()) throw std::runtime_error("RestartFile not open !");
      int Dim = dimension();
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      {
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> initialPositions_[iOfParticle](dim) >> std::ws; 
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> initialMomenta_[iOfParticle](dim) >> std::ws;
        initialConditions >> initialInternalEnergies_[iOfParticle] >> std::ws;
      }
      cout << " Successful reading of the settings files... " << endl;
    }
  }
    
  

  const std::string& Input::inputPath() const
  {
    return inputPath_;
  }

  const std::ifstream& Input::inputFlux() const
  {
    return inputFlux_;
  }


}
