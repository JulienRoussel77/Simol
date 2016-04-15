#include "simol/statphys/output/Output.hpp"

namespace simol{
 
  ofstream & Output::outParticlesXMakeMol()
  {
    return *outParticlesXMakeMol_;
  }
  
  /// 
  ///-- display the current configuration in XMakemol format ; specific for NBody systems ---
  void Output::displayParticlesXMakeMol(vector<Particle> const& configuration, int iOfStep, double domainSize)
  {
    outParticlesXMakeMol() << nbOfParticles_ << endl;
    outParticlesXMakeMol() << "Time = " << iOfStep * timeStep() << endl;
    double coordinate = 0;
    int Dim = dimension_;
    for (int i = 0; i < nbOfParticles_; i++)
    {
      outParticlesXMakeMol() << " O  ";
      if (Dim == 3)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
	  //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
	  coordinate = configuration[i].position(dim);
	  coordinate -= rint(coordinate/domainSize)*domainSize;
	  outParticlesXMakeMol() << coordinate << " "; 
        }
      }
      else if (Dim == 2)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
	      //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
	      coordinate = configuration[i].position(dim);
	      coordinate -= rint(coordinate/domainSize)*domainSize;
	      outParticlesXMakeMol() << coordinate << " "; 
        }
        outParticlesXMakeMol() << 0 << " ";  
      }
      outParticlesXMakeMol() << endl;
    }

  }
  
}
