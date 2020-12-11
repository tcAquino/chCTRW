//
//  Reaction.h
//  Stochastic
//
//  Created by Tomas Aquino on 5/22/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

//  Reaction handlers to use with different algorithms
//  Note that different algorithms may require implementing different methods

#ifndef Reaction_h
#define Reaction_h

#include <cmath>
#include "Stochastic/Stoichiometry.h"
#include "general/Operations.h"

namespace stochastic
{
  // Generic mass-action reactions
  class Reaction_MassAction
  {
  public:
    using Stoichiometry = stochastic::Stoichiometry;
    const Stoichiometry stoichiometry;  // Reaction rate and stoichiometric coefficients
    
    Reaction_MassAction(Stoichiometry stoichiometry)
    : stoichiometry(stoichiometry)
    {
      double factor = 1.;
      for (auto const& react : stoichiometry.reactants)
        factor *= operation::factorial(react.second);
      reaction_rate_scaled = stoichiometry.reaction_rate/factor;
    }

    // For continuous concentration
    double rate(std::vector<double> const& concentration) const
    {
      double combinations = 1.;
      for (auto const& sto : stoichiometry.reactants)
        combinations *= std::pow(concentration[sto.first], sto.second);

      return reaction_rate_scaled*combinations;
    }
    
    // For discrete particle numbers
    double rate(std::vector<std::size_t> const& numbers) const
    {
      std::size_t combinations = 1;
      for (auto const& sto : stoichiometry.reactants)
        combinations *= operation::factorial_incomplete(numbers[sto.first], sto.second);

      return reaction_rate_scaled*combinations;
    }

    // For continuous concentration
    void react(std::vector<double>& concentration, double time_step) const
    {
      double rate_val = rate(concentration);

      for (auto const& sto : stoichiometry.reactants)
        concentration[sto.first] -= sto.second*rate_val*time_step;
      for (auto const& sto : stoichiometry.products)
        concentration[sto.first] += sto.second*rate_val*time_step;
    }

    // For discrete particle numbers
    void react(std::vector<std::size_t>& numbers) const
    {
      for (auto const& sto : stoichiometry.reactants)
        numbers[sto.first] -= sto.second;
      for (auto const& sto : stoichiometry.products)
        numbers[sto.first] += sto.second;
    }

    // Generic interface for both discrete and continuous
    template <typename concentration_type>
    void operator()
    (std::vector<concentration_type>& concentration, double time_step, double time = 0)
    { React(concentration, time_step); }
    
    private:
      double reaction_rate_scaled;
  };

  class Reaction_concentration_bimolecular_analytical
  {
  public:
    const std::size_t nr_types{ 2 };
    const double reaction_rate;
    
    Reaction_concentration_bimolecular_analytical
    (double reaction_rate, double tol = 1.e-10)
    : reaction_rate{ reaction_rate }
    , tol{ tol }
    {}
    
    Reaction_concentration_bimolecular_analytical
    (double reaction_rate, double mass0, double mass1, double tol = 1.e-10)
    : reaction_rate{ reaction_rate }
    , masses{ mass0, mass1 }
    , tol{ tol }
    {}
    
    Reaction_concentration_bimolecular_analytical
    (double reaction_rate, std::vector<double> const& concentration, double tol = 1.e-10)
    : reaction_rate{ reaction_rate }
    , masses{ concentration }
    , tol{ tol }
    {}

    void set(std::size_t type, double val)
    { masses[type]=val; }
    
    void set(std::vector<double> const& concentration)
    {
      for (std::size_t ii = 0; ii < nr_types; ++ii)
        set(ii, concentration[ii]);
    }

    void evolve(double time_max)
    {
      double time_step = time_max - time_current;
      time_current = time_max;

      std::size_t max_idx = masses[0] > masses[1] ? 0 : 1;
      double mass_max = masses[max_idx];
      double mass_min = masses[!max_idx];
      double diff = mass_max - mass_min;

      if (diff > tol)
      {
        double exp_val = std::exp(-reaction_rate*time_step*diff);
        double sol_base = diff/(mass_max - exp_val*mass_min);
        masses[max_idx] = mass_max*sol_base;
        masses[!max_idx] = mass_min*sol_base*std::exp(-reaction_rate*time_step*diff);
      }
      else
      {
        double solution_equal = mass_max/(1. + reaction_rate*mass_max*time_step);
        masses[0] = solution_equal;
        masses[1] = solution_equal;
      }
    }

    void time(double val)
    { time_current = val; }
    
    double mass(std::size_t type)
    { return masses[type]; }

    double particles(std::size_t type)
    { return mass(type); }

  private:
    std::vector<double> masses{ 0., 0. };
    double time_current{ 0. };
    const double tol;
  };
  
  class Reaction_concentration_decay_analytical
  {
  public:
    const std::size_t nr_types{ 1 };
    const double reaction_rate;
    
    Reaction_concentration_decay_analytical
    (double reaction_rate)
    : reaction_rate{ reaction_rate }
    {}
    
    Reaction_concentration_decay_analytical
    (double reaction_rate, double mass)
    : reaction_rate{ reaction_rate }
    , masses{ mass }
    {}
    
    Reaction_concentration_decay_analytical
    (double reaction_rate, std::vector<double> const& concentration)
    : reaction_rate{ reaction_rate }
    , masses{ concentration[0] }
    {}

    void set(double val)
    { masses=val; }
    
    void set(std::size_t type, double val)
    { masses=val; }
    
    void set(std::vector<double> const& concentration)
    { set(0, concentration[0]); }

    void evolve(double time_max)
    {
      double time_step = time_max - time_current;
      time_current = time_max;
      masses *= std::exp(-reaction_rate*time_step);
    }

    void time(double val)
    { time_current = val; }

    double mass(std::size_t type)
    { return masses; }

    double particles(std::size_t type)
    { return mass(type); }

  private:
    double masses{ 0. };
    double time_current{ 0. };
  };
}

#endif /* Reaction_h */
