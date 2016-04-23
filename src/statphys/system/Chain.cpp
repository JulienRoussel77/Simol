//#include <core/io/CommandLine.hpp>

#include "simol/statphys/system/Chain.hpp"

namespace simol
{
  //### Chain ###

  Chain::Chain(Input const& input):
    System(input)
  {
    assert(configuration_.size() > 1);
    for (int i = 0; i < input.nbOfParticles(); i++)
    {
      configuration_[i] = Particle(input.mass(), input.initialPosition(i), input.initialMomentum(i));
      //std::cout << configuration_[i].force() << std::endl;
    }
  }

  void Chain::thermalize(Dynamics& /*dyna*/) {}

  //void Chain::computeAllForces(Dynamics const& /*model*/)
  //{throw std::invalid_argument("computeAllForces : Function undefined");};*/

  void Chain::thermalize(LangevinBase& dyna)
  {
    //for (auto&& particle : configuration_)
    //dyna.updateBefore(particle);
    for (int i = 0; i < nbOfParticles(); i++)
      dyna.updateBefore(getParticle(i));

    computeAllForces();

    for (int i = 0; i < nbOfParticles(); i++)
      dyna.updateAfter(getParticle(i));

    for (int i = 0; i < nbOfParticles(); i++)
    {
      double localTemperature = dyna.temperatureLeft() + i * dyna.deltaTemperature() / nbOfParticles();
      dyna.updateOrsteinUhlenbeck(getParticle(0), 1 / localTemperature);
    }
  }





  //###### BiChain ######

  BiChain::BiChain(Input const& input):
    Chain(input),
    ancorParticle_(input.dimension())
  {
    //ancorParticle_.position(0) = 2 * input.initialPosition(0) - input.initialPosition(1);
    ancorParticle_.position(0) = 0;
  }







  void BiChain::computeAllForces()
  {
    //std::cout << "BiChain::computeAllForces" << std::endl;
  for (auto && particle : configuration_)
      particle.resetForce(potential());
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    interaction(ancorParticle_, configuration_[0]);
    for (int i = 0; i < nbOfParticles() - 1; i++)
    {
      interaction(configuration_[i], configuration_[i + 1]);
    }
  }

  void BiChain::computeProfile(Output& /*output*/, Dynamics const& /*dyna*/, int /*iOfStep*/) const
  {}

  void BiChain::computeProfile(Output& output, LangevinBase const& dyna, int iOfStep) const
  {
    output.energySumFlow() = 0;
    int midNb = (nbOfParticles() - 1) / 2;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      //Particle& particle = configuration_[iOfParticle];
      double dist = 0;
      double flow = 0;
      double potTempTop = 0;
      double potTempBot = 0;

      if (iOfParticle == 0)
      {
        dist = getParticle(0).position(0) - ancorParticle_.position(0);
        flow = dyna.gamma() * (dyna.temperatureLeft() - 2.0 * getParticle(0).kineticEnergy());
      }
      else
      {
        // dist is r_iOfParticle
        dist = getParticle(iOfParticle).position(0) - getParticle(iOfParticle - 1).position(0);
        // flow is j_iOfParticle
        flow = - getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle - 1).momentum(0);
        if (iOfParticle != nbOfParticles() - 1)
        {
          output.energySumFlow() += flow;
          if (iOfParticle == midNb)
            output.energyMidFlow() = flow;
        }
      }

      if (iOfParticle == nbOfParticles() - 1)
      {
        potTempTop = pow(getParticle(iOfParticle).energyGrad(0), 2);
        potTempBot = getParticle(iOfParticle).energyLapla();
      }
      else
      {
        potTempTop = pow(getParticle(iOfParticle).energyGrad(0) - getParticle(iOfParticle + 1).energyGrad(0), 2);
        potTempBot = getParticle(iOfParticle).energyLapla() + getParticle(iOfParticle + 1).energyLapla();
      }

      output.appendBendistProfile(dist , iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * getParticle(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFlowProfile(flow, iOfStep, iOfParticle);
    }
    output.energySumFlow() /= (nbOfParticles() - 2.);
    //cout << output.energySumFlow() << endl;
  }




  //###### TriChain ######


  TriChain::TriChain(Input const& input):
    Chain(input),
    ancorParticle1_(input.dimension()),
    ancorParticle2_(input.dimension()),
    isOfFixedVolum_(input.isOfFixedVolum())
  {
    ancorParticle1_.position(0) = 0;//3 * input.initialPosition(0) - 2*input.initialPosition(1);
    ancorParticle2_.position(0) = 0;//2 * input.initialPosition(0) - input.initialPosition(1);
  }

  bool const& TriChain::isOfFixedVolum() const
  {
    return isOfFixedVolum_;
  }

  ///Computes the force and the energy associated to this triplet interaction, and updates these 2 fields
  ///The first 2 derivates of the potential are stored in "particle2"
  void TriChain::triInteraction(Particle& particle1, Particle& particle2, Particle& particle3) const
  {
    Vector<double> delta = particle3.position() - 2.0 * particle2.position() + particle1.position();
    //double d12 = r12.norm();
    double energy123 = potential(delta);
    Vector<double> force123 = potentialForce(delta);    // = - v'(r_2)
    double lapla123 = laplacian(delta);

    particle2.potentialEnergy() = energy123;
    particle1.force() += force123;
    particle2.force() -= 2.0 * force123;
    particle3.force() += force123;
    particle2.energyGrad() = -force123;    // - v'(r_2)
    particle2.energyLapla() = lapla123;    // v''(r_2)
  }

  void TriChain::computeAllForces()
  {
    //std::cout << "TriChain::computeAllForces" << std::endl;
  for (auto && particle : configuration_)
      particle.resetForce(potential());
    //for (auto&& particle : configuration_)
    //  dyna.computeForce(particle);
    triInteraction(ancorParticle1_, ancorParticle2_, configuration_[0]);
    triInteraction(ancorParticle2_, configuration_[0], configuration_[1]);
    for (int i = 0; i < nbOfParticles() - 2; i++)
      triInteraction(configuration_[i], configuration_[i + 1], configuration_[i + 2]);
    //dyna.bending(configuration_[nbOfParticles() - 2], configuration_[nbOfParticles() - 1]);
    if (isOfFixedVolum_)
      triInteraction(configuration_[nbOfParticles() - 1], ancorParticle1_, ancorParticle2_);
  }

  double TriChain::boundaryPotEnergy() const
  {return ancorParticle1_.potentialEnergy() + ancorParticle2_.potentialEnergy();}

  void TriChain::computeProfile(Output& /*output*/, Dynamics const& /*dyna*/, int /*iOfStep*/) const
  {}

  void TriChain::computeProfile(Output& output, LangevinBase const& dyna, int iOfStep) const
  {
    output.energySumFlow() = 0;
    for (int iOfParticle = 0; iOfParticle < nbOfParticles(); iOfParticle++)
    {
      //Particle& particle = configuration_[iOfParticle];
      double bending = 0;
      double flow = 0;
      double potTempTop = 0;
      double potTempBot = 0;

      if (iOfParticle == 0)
      {
        bending = ancorParticle2_.position(0) - 2.0 * getParticle(0).position(0) + getParticle(1).position(0);
        flow = dyna.gamma() * (dyna.temperatureLeft() - 2.0 * getParticle(0).kineticEnergy())
               - ancorParticle2_.energyGrad(0) * getParticle(0).momentum(0);
        potTempTop = pow(- ancorParticle2_.energyGrad(0) + 2.0 * getParticle(0).energyGrad(0) - getParticle(1).energyGrad(0), 2);
        potTempBot = ancorParticle2_.energyLapla() + 4.0 * getParticle(0).energyLapla() + getParticle(1).energyLapla();

      }
      else if (iOfParticle < (int)configuration_.size() - 1)
      {
        // bending is k_iOfParticle
        bending = getParticle(iOfParticle - 1).position(0) - 2.0 * getParticle(iOfParticle).position(0) + configuration_[iOfParticle + 1].position(0);
        // flow is j_iOfParticle
        flow = - getParticle(iOfParticle - 1).energyGrad(0) * getParticle(iOfParticle).momentum(0)
               + getParticle(iOfParticle).energyGrad(0) * getParticle(iOfParticle - 1).momentum(0);

        potTempTop = pow(- getParticle(iOfParticle - 1).energyGrad(0) + 2.0 * getParticle(iOfParticle).energyGrad(0) - getParticle(iOfParticle + 1).energyGrad(0), 2);
        potTempBot = getParticle(iOfParticle - 1).energyLapla() + 4.0 * getParticle(iOfParticle).energyLapla() + getParticle(iOfParticle + 1).energyLapla();
        output.energySumFlow() += flow;
      }
      else
      {
        bending = nan("");
        flow = - getParticle(iOfParticle - 1).energyGrad(0) * getParticle(iOfParticle).momentum(0);
        potTempTop = pow(- getParticle(iOfParticle - 1).energyGrad(0), 2);
        potTempBot = getParticle(iOfParticle - 1).energyLapla();
      }

      int midNb = (nbOfParticles() - 1) / 2;
      if (iOfParticle == midNb)
        output.energyMidFlow() = flow;

      output.appendBendistProfile(bending , iOfStep, iOfParticle);
      output.appendKinTempProfile(2 * getParticle(iOfParticle).kineticEnergy(), iOfStep, iOfParticle);
      output.appendPotTempTopProfile(potTempTop, iOfStep, iOfParticle);
      output.appendPotTempBotProfile(potTempBot, iOfStep, iOfParticle);
      output.appendFlowProfile(flow, iOfStep, iOfParticle);
    }
    output.energySumFlow() /= (nbOfParticles() - 2.);
    //cout << output.energySumFlow() << endl;
  }




}
