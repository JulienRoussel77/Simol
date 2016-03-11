#ifndef SIMOL_POTENTIAL_HPP
#define SIMOL_POTENTIAL_HPP

#include "tools.hpp"
#include "Vector.hpp"	
#include "input.hpp"
#include "RNG.hpp"


namespace simol
{
  class Potential;
  
  Potential* createPotential(Input const& input, int const& indexOfReplica=1);

  class Potential
  {
	friend Potential* createPotential(Input const& input, int const& indexOfReplica);
	
	public:
		Potential();
		virtual ~Potential(){};
		Potential(Input const& input, int const& indexOfReplica=1);
		virtual double operator()(dvec const & position) const;
		virtual double operator()(double position) const;
		double value(dvec const& position) const;
		double value(double position) const;
		virtual dvec derivative(dvec const & position) const;
		virtual dvec derivative(double position) const;
		virtual dvec force(dvec const & position) const;
		virtual dvec force(double position) const;
		virtual double laplacian(dvec const & position) const;
		virtual double laplacian(double position) const;
		virtual double ratioToHarmonic() const {assert(false); return 0;};
		virtual double drawLaw(double /*localBeta*/, RNG* /*rng*/){assert(false); return 0;};

  };
  
  class Sinusoidal : public Potential{
    public:
      Sinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
      double laplacian(double position) const;

    private:
      double amplitude_;
      double pulsation_;
  };
  
    class SumSinusoidal : public Potential{
    public:
      SumSinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
      virtual double laplacian(double position) const;


    private:
      double amplitude_;
      double pulsation_;
  };
  
  class FracSinusoidal : public Potential{
    public:
      FracSinusoidal(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
      virtual double laplacian(double position) const;


    private:
      double amplitude_;
      double pulsation_;
  };
  
  class DoubleWell : public Potential{
    public:
      DoubleWell(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;

    private:
      double height_;
      double interWell_;
  };
  
    class HarmonicWell : public Potential{
    public:
      HarmonicWell(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;

    private:
      double stiffness_;
  };

  
  
  
  class Harmonic : public Potential{
    public:
      Harmonic(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
			double laplacian(double position) const;
			double drawLaw(double localBeta, RNG* rng);
    private:
      double stiffness_;
  };
  
  class Rotor : public Potential{
    public:
      Rotor(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
      double laplacian(double position) const;
			double drawLaw(double localBeta, RNG* rng);
    private:
      double stiffness_;
  };

  
    class Quadratic : public Potential{
    public:
      Quadratic(Input const& input, int const& indexOfReplica=1);
      double operator()(double position) const;
      dvec derivative(double position) const;
			double laplacian(double position) const;
			virtual double ratioToHarmonic() const;
			double drawLaw(double localBeta, RNG* rng);
    private:
      double stiffness_, alpha_, beta_;
  };
	
	class SpaceSinus : public Potential{
    public:
      SpaceSinus(Input const& input, int const& indexOfReplica=1);
      double operator()(dvec const& position) const;
      dvec derivative(dvec const& position) const;
      double laplacian(dvec const& position) const;

    private:
      double amplitude_;
      double pulsation_;
  };
}

//#include "potential.ipp"

#endif
