#include "simol/output/Output.hpp"

namespace simol
{
  ///
  ///-- display the current configuration in XMakemol format ; specific for NBody systems ---
  void Output::displayXMakeMol(System const& syst, long int iOfStep, double /*domainSize*/)
  {
    outXMakeMol() << nbOfParticles_ << endl;
    outXMakeMol() << "Time = " << iOfStep * timeStep() << endl;
    //double coordinate = 0;
    //int Dim = dimension_;
    DVec pos = DVec::Zero(3);
    for (int iOfParticle = 0; iOfParticle < nbOfParticles_; iOfParticle++)
    {
      if (syst(iOfParticle).type() == 0)
        outXMakeMol() << " O  ";
      else if (syst(iOfParticle).type() == -1)
        outXMakeMol() << " He  ";
      else
        outXMakeMol() << " Ar  ";
      
      pos.head(dimension()) = syst.periodicImage(syst(iOfParticle).position());
      outXMakeMol() << pos.adjoint() << endl;
      
      /*if (Dim == 3)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
          //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
          coordinate = syst(iOfParticle).position(dim);
          coordinate -= rint(coordinate / domainSize) * domainSize;
          outXMakeMol() << coordinate << " ";
        }
      }
      else if (Dim == 2)
      {
        for (int dim = 0; dim < Dim; dim++)
        {
          //-- recenter all the coordinates in the interval [-domainSize/2, domainSize/2] --
          coordinate = syst(iOfParticle).position(dim);
          outXMakeMol() << coordinate << " ";
          coordinate -= rint(coordinate / domainSize) * domainSize;

        }
        outXMakeMol() << 0 << " ";
      }
      outXMakeMol() << endl;
    }*/

    }
  }
}
