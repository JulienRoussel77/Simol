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
  
  bool isYes(const string& s)
  {
    return sameLetters(s, "yes");
  }

  ///
  ///Reads the next element in the file read
  double readItem(ifstream& in)
  {
    assert(in.is_open());
    double item;
    in >> item;
    return item;
  }

  DVec scanTensor(const string path, vector<int>& dimensions)
  {
    cout << "scanTensor" << endl;
    ifstream inTensor(path);
    if (!inTensor.is_open())
      throw runtime_error("scanTensor couldn't read from " + path + "!");
    string str;
    int dim;
      
    inTensor >> str >> str >> dim;
    cout << "First dim is " << dim << endl;
    dimensions.push_back(dim);
    inTensor >> dim;
    cout << "Second dim is " << dim << endl;
    dimensions.push_back(dim);
    
    int nbOfValues = dimensions[0]*dimensions[1];
    DVec tensor(nbOfValues);
    double  value;
    int iOfVal = 0;
    while (inTensor >> value)
    {
      //cout << "value : " << value << endl;
      tensor[iOfVal] = value;
      iOfVal++;
    }
    if (iOfVal != nbOfValues)
      throw runtime_error("In file " + path + " the header dimensions are not compatible with the numer of lines!");
    return tensor;
  }
  
  map<double, double> scanMap(const string path)
  {
    cout << "Reading the spectral gaps :" << endl;
    ifstream inMap(path);
    if (!inMap.is_open())
      throw runtime_error(path+" is not a valid path!");
    string str;
    inMap >> str >> str >> str;
    double key, value;
    map<double, double> map;
    while (inMap >> key >> value)
    {
      cout << key << " -> " << value << endl;
      map[key] = value;
    }
    
    return map;
  }
  
  ///Reads the input file using the YAML library
  ///Reads the initial conditions in the "settings" file, if indicated
  Input::Input(CommandLine cmd):
    data(YAML::LoadFile(cmd.inputFileName())),
    inputPath_(cmd.inputFileName()),
    inputFlux_(inputPath()),
    initialPositions_(nbOfParticles(), DVec(dimension())),
    initialMomenta_(nbOfParticles(), DVec(dimension())),
    initialInternalEnergies_(nbOfParticles(), 0)
  {
    assert(data.IsDefined());
    if (doRestart())
    {
      cout << " -- reading from restart file " << restartFileName() << endl;
      ifstream initialConditions(restartFileName());
      if (!initialConditions.is_open()) throw std::runtime_error("RestartFile not open!");
      int Dim = dimension();
      for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
      {
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> initialPositions_[iOfParticle](dim) >> std::ws; 
        for (int dim = 0; dim < Dim; dim++)
          initialConditions >> initialMomenta_[iOfParticle](dim) >> std::ws;
        if (doDPDE()) // reading an additional field: internal energy
	  initialConditions >> initialInternalEnergies_[iOfParticle] >> std::ws;
      }
      cout << " -- successful reading of the settings files... " << endl;
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
