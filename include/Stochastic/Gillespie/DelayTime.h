//
//  DelayTime.h
//  Stochastic
//
//  Created by Tomas Aquino on 1/31/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef DelayTime_h
#define DelayTime_h

//  Delay times for e.g. Gillespie
//  Delay classes must implement a double operator() (double time)
//  which returns the delay given a time window
//  Compound delays use NumberProcess classes
//  NumberProcess classes must implement a std::size_t operator() (double time)
//  which returns the number of i.i.d. delay events given a time window

#include <random>
#include <vector>
#include "Stochastic/Random.h"

namespace stochastic
{
  class DelayTime_NoDelay
  {
  public:
    double operator() (double) const
    { return 0.; }
  };

  class DelayTime_Exponential
  {
  public:
    const double mean;

    DelayTime_Exponential(double mean = 1.)
    : mean{ mean }
    , exp_distribution{ 1./mean }
    {}

    double operator() (double time = 0.)
    { return exp_distribution(rng); }

  private:
    std::mt19937 rng{ std::random_device{}() };
    std::exponential_distribution< double > exp_distribution;
  };

  class DelayTime_SkewedLevyStable
  {
  public:
    const double alpha;
    const double sigma;
    const double mu;

    DelayTime_SkewedLevyStable(double alpha, double sigma = 1., double mu = 0.)
    : alpha(alpha)
    , sigma(sigma)
    , mu(mu)
    {}

    double operator() (double time = 0.)
    { return stochastic::skewedlevystable_distribution<double>{ alpha, sigma, mu }(rng); }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
  };

  class DelayTime_Gamma
  {
  public:
    const double gamma;
    const double mu;

    DelayTime_Gamma(double gamma, double mu = 1.)
    : gamma(gamma)
    , mu(mu)
    {}

    double operator() (double time = 0.)
    { return std::gamma_distribution< double >{ gamma, mu }(rng); }

  private:
    std::mt19937 rng{ std::random_device{}() };
  };

  // Poisson process
  class NumberProcess_Poisson
  {
  public:
    const double rate;

    NumberProcess_Poisson(double rate)
    : rate(rate)
    {}

    std::size_t operator() (double time)
    {
      return std::poisson_distribution<std::size_t>{ rate*time }(rng);
    }

  private:
    std::mt19937 rng{ std::random_device{}() };
  };

  //  Generic compound waiting time
  template <typename Number_process, typename Waiting_process>
  class DelayTime_Compound
  {
  public:
    DelayTime_Compound(Number_process number_process, Waiting_process waiting_process)
    : number_process(number_process)
    , waiting_process(waiting_process)
    {}

    double operator() (double time)
    {
      double delay = 0.;
      std::size_t number = number_process(time);
      for (std::size_t ii = 0; ii < number; ++ii)
      {
        delay += waiting_process();
      }
      return delay;
    }

  private:
    Number_process number_process;
    Waiting_process waiting_process;
  };

  //  Compound (Number-Process)-Exponential
  template <typename Number_process>
  class DelayTime_CompoundExponential
  {
  public:
    const double gamma;
    const double mu;

    DelayTime_CompoundExponential
    (Number_process number_process, double gamma, double mu = 1.)
    : mu(mu)
    , gamma(gamma)
    , number_process(number_process)
    {}

    double operator() (double time) const
    {
      std::size_t number = number_process(time);
      return (number != 0
          ? std::gamma_distribution< double >{ number, mu }(rng)
          : 0.);
    }

  private:
    Number_process number_process;
    mutable std::mt19937 rng{ std::random_device{}() };
  };

  // Compound (Number-Process)-SkewedLevyStable
  template <typename Number_process>
  class DelayTime_CompoundSkewedLevyStable
  {
  public:
    const double alpha;
    const double sigma;
    const double mu;

    DelayTime_CompoundSkewedLevyStable
    (Number_process number_process, double alpha, double sigma = 1., double mu = 0.)
    : alpha(alpha)
    , sigma(sigma)
    , mu(mu)
    , number_process(number_process)
    {}

    double operator() (double time)
    {
      std::size_t number = number_process(time);
      return
      (number != 0
       ? stochastic::skewedlevystable_distribution<double>{
        alpha, std::pow(number, 1./alpha)*sigma, number*mu }(rng)
       : 0.);
    }

  private:
    Number_process number_process;
    std::mt19937 rng{ std::random_device{}() };
  };

  // Subordinator formulation of skewed-levy-stable delay
  class DelayTime_Subordinator_SkewedLevyStable
  {
  public:
    const double alpha;
    const double gamma;
    const double sigma;
    const double mu;

    DelayTime_Subordinator_SkewedLevyStable
    (double alpha, double gamma = 1., double sigma = 1., double mu = 0.)
    : alpha(alpha)
    , gamma(gamma)
    , sigma(sigma)
    , mu(mu)
    {}

    double operator() (double delta_time)
    {
      return std::pow(gamma * delta_time, 1./alpha)*stable_dist(rng) + mu;
    }

  private:
    std::mt19937 rng{ std::random_device{}() };
    stochastic::skewedlevystable_distribution<double> stable_dist{ alpha, sigma, 0. };
  };

  // Subordinator formulation of skewed-levy-stable delay
  // Remove the contribution of regular reaction time and keep just delay
  class DelayTime_Subordinator_SkewedLevyStable_JustDelay
  {
  public:
    const double alpha;
    const double gamma;
    const double sigma;
    const double mu;

    DelayTime_Subordinator_SkewedLevyStable_JustDelay
    (double alpha, double gamma = 1., double sigma = 1., double mu = 0.)
    : alpha(alpha)
    , gamma(gamma)
    , sigma(sigma)
    , mu(mu)
    {}

    double operator() (double delta_time)
    {
      return -delta_time + std::pow(gamma*delta_time, 1./alpha)*stable_dist(rng) + mu;
    }

  private:
    std::mt19937 rng{ std::random_device{}() };
    stochastic::skewedlevystable_distribution<double> stable_dist{ alpha, sigma, 0. };
  };
}

#endif /* DelayTime_h */


