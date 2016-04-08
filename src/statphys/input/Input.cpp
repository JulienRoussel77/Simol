#include "Input.hpp"

using std::cout;
using std::endl;
using std::string;
using std::max;
using std::min;

namespace simol {
  
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
    inputSettings_(settingsPath())
  {
    if (doFileSettings())
      {
	assert(inputFlux_.is_open());
	cout << "Reading the settings from " << settingsPath() << "...";
	initialPositions_ = vector<Vector<double>>(nbOfParticles());
	initialMomenta_ = vector<Vector<double>>(nbOfParticles());
	for (int iOfParticle=0; iOfParticle < (int)nbOfParticles(); iOfParticle++)
	  {
	    readItem(inputSettings_);
	    assert( (int) readItem(inputSettings_) == (int) iOfParticle);
	    for (int i=0; i < dimension(); i++)
	      initialPositions_[iOfParticle](i) = readItem(inputSettings_);
	    for (int i=0; i < dimension(); i++)
	      initialMomenta_[iOfParticle](i) = readItem(inputSettings_);
	    for (int i=0; i<4; i++)
	      readItem(inputSettings_);
	    //cout << iOfParticle << " " << initialPositions_[iOfParticle] << " " << initialMomenta_[iOfParticle] << endl;
	  }
	cout <<" Successful reading of the settings files... " << endl;
      }
      cout << "Successful reading of input variables" << endl;
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
