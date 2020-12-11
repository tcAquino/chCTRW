//
//  batch_delay.cpp
//
//  Created by Tomas Aquino on 3/6/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include "general/Constants.h"
#include "general/Operations.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/Reaction.h"
#include "Stochastic/Stoichiometry.h"
#include "Stochastic/Gillespie/Gillespie_Stoichiometric.h"

int main(int argc, const char * argv[])
{
  //  Initial particle numbers of each species type
  std::vector<std::size_t> particles_initial{ 100000, 100000 };
  
  //  Reaction stoichiometries
  double reaction_rate = 1./std::pow(particles_initial[0], particles_initial.size()-1);
    stochastic::Stoichiometry stoichiometry_1{ reaction_rate, { { 0, 1 }, { 1, 1 } }, {} };

  //  Delay properties
  double delay_exponent = 0.75;
  double delay_characteristic_time = 0.1;
  double delay_rate = 10.*std::pow(delay_characteristic_time, -delay_exponent);
  double delay_characteristic_time_scaled = std::pow(
                                                     std::cos(constants::pi*delay_exponent/2.)*
                                                      delay_characteristic_time,
                                                     1./delay_exponent);
  using NumberProcess = stochastic::NumberProcess_Poisson;
  using Delay = stochastic::DelayTime_CompoundSkewedLevyStable<NumberProcess>;

  //  Max simulation time
  double time_max = 1e5;

  //  Nr of ensembles to average over
  std::size_t nr_ensembles = 100;

  //  Prepare stuff for data
  double time_min = 1e-2;
  std::size_t nr_measures = 30;
  std::vector< double > measure_times = range::logspace<std::vector<double>>
    (time_min, time_max, nr_measures);
  std::vector< double > concentration(nr_measures);

  //  Make Gillespie simulator
  auto gillespie = gillespie::make_Gillespie_MassAction_Delay(
                           particles_initial,
                           Delay{
                            delay_rate,
                            delay_exponent, delay_characteristic_time_scaled },
                           stoichiometry_1);
  
  //  Run each ensemble of particles
  //  Measure number concentration over time of species 0
  for (std::size_t ensemble = 0; ensemble < nr_ensembles; ++ensemble)
  {
    std::cout << "ensemble = " << ensemble << "\n";
    gillespie.set(particles_initial);
    double particles = 0.;
    for (std::size_t measure = 0; measure < nr_measures; ++measure)
    {
      std::cout << "\ttime = " << gillespie.time() << "\n";
      while (gillespie.time() < measure_times[ measure ])
      {
        particles = gillespie.particles(0);
        gillespie.evolve();
      }
      concentration[measure] += particles;
    }
  }
  operation::div_scalar_InPlace(concentration, nr_ensembles);

  //  Output
  std::string output_dir = "../output";
  std::string filename{ "Data_Gillespie_Delay_Example_CompoundStable.dat" };
  std::ofstream output{ output_dir + "/" + filename };
  if (!output.is_open())
    throw useful::open_write_error(filename);
  output << std::scientific << std::setprecision(8);
  output << 0. << "\t";
  useful::print(output, measure_times);
  output << "\n";
  output << double(particles_initial[0]) << "\t";
  useful::print(output, concentration);
  output << "\n";
  output.close();
  
  return 0;
}
