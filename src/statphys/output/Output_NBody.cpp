#include "simol/statphys/output/Output.hpp"

namespace simol
{
  ///
  ///-- display the current configuration in XMakemol format ; specific for NBody systems ---
  void Output::displayXMakeMol(vector<Particle*> const& configuration, long int iOfStep, double domainSize)
  {
    outXMakeMol() << nbOfParticles_ << endl;
    outXMakeMol() << "Time = " << iOfStep * timeStep() << endl;
    double coordinate = 0;
    int Dim = dimension_;
    for (int i = 0; i < nbOfParticles_; i++)
    {
      outXMakeMol() << " O  ";
      if (Dim == 3)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
          //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
          coordinate = configuration[i]->position(dim);
          coordinate -= rint(coordinate / domainSize) * domainSize;
          outXMakeMol() << coordinate << " ";
        }
      }
      else if (Dim == 2)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
          //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
          coordinate = configuration[i]->position(dim);
          coordinate -= rint(coordinate / domainSize) * domainSize;
          outXMakeMol() << coordinate << " ";
        }
        outXMakeMol() << 0 << " ";
      }
      outXMakeMol() << endl;
    }

  }

}
