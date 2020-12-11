//
//  streamtube_gillespie.cpp
//
//  Created by Tomas Aquino on 3/4/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <typeinfo>
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/Gillespie/Gillespie_Stoichiometric.h"
#include "Stochastic/Reaction.h"
#include "Stochastic/Stoichiometry.h"
#include "Stochastic/Streamtube/Models.h"
#include "Stochastic/Streamtube/Measurer.h"

int main(int argc, const char* argv[])
{
  if (argc == 0)
  {
    std::cout << "streamtube_gillespie\n";
    std::cout << "Parameters (default value in []):\n"
              << "characteristic_length_reactive : Characteristic length of reactive patches\n"
              << "exp_length_reactive : Exponent for stable reactive patches (ignored for exponential)\n"
              << "characteristic_length_conservative : Characteristic length of conservative patches\n"
              << "exp_length_conservative : Exponent for stable conservative patches (ignored for exponential)\n"
              << "mean_advection : Mean advection across streamtubes\n"
              << "var_advection : Variance of advection across streamtubes (ignored for one-parameters dists)\n"
              << "reaction_rate : Macroscopic reaction rate\n"
              << "measure_min : Minimum time or distance for output\n"
              << "measure_max : Maximum time or distance for output\n"
              << "nr_measures : Number of outputs\n"
              << "flux_weighted : 0 - Homogeneous injection\n"
              << "                1 - Flux-weighted injection\n"
              << "particles_mobile_each : Initial particle numbers for each mobile species\n"
              << "particles_immobile_each : Inital particles numbers for each immobile species\n"
              << "nr_fixed_velocity : Number of streamtubes for each velocity value\n"
              << "nr_velocities : Number of separate velocity samples\n"
              << "run_nr : Tag to record same-parameter realizations to different files\n"
              << "output_dir : Directory to output to [../output]";
    return 0;
  }
  
  if (argc != 17 && argc != 18)
    throw useful::bad_parameters();

	using namespace streamtube::model_uniform_exp_exp;

  //  Parameters
	std::size_t arg = 1;
	double characteristic_length_reactive = atof(argv[arg++]);
	double exp_length_reactive = atof(argv[arg++]);
	double characteristic_length_conservative = atof(argv[arg++]);
	double exp_length_conservative = atof(argv[arg++]);  //
	double mean_advection = atof(argv[arg++]);
	double var_advection = atof(argv[arg++]);
	double reaction_rate = atof(argv[arg++]);
	double measure_min = atof(argv[arg++]);
	double measure_max = atof(argv[arg++]);
	std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
	bool flux_weighted = atoi(argv[arg++]);
	std::size_t particles_mobile_each = strtoul(argv[arg++], NULL, 0);
	std::size_t particles_immobile_each = strtoul(argv[arg++], NULL, 0);
	std::size_t nr_fixed_velocity = strtoul(argv[arg++], NULL, 0);
	std::size_t nr_velocities = strtoul(argv[arg++], NULL, 0);
	std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string output_dir = argc > arg ? argv[arg++] : "../output";
  
  double tortuosity = 1.;

	std::vector< std::size_t > average_initial_mobile_particles{ particles_mobile_each };
	std::vector< std::size_t > average_initial_immobile_particles{ particles_immobile_each };
	std::size_t mobile_types = average_initial_mobile_particles.size();
	std::size_t immobile_types = average_initial_immobile_particles.size();
	std::size_t types = mobile_types + immobile_types;
	double particles_characteristic = (mobile_types*particles_mobile_each
                                     + immobile_types*particles_immobile_each)
                                    /double(types);
  
  using Stoichiometry = stochastic::Stoichiometry;
  // A+B\to\varnothing
  Stoichiometry stoichiometry{
    reaction_rate/particles_characteristic,
    { { 0, 1 }, { 1, 1 } }, {} };

  //  Normalization of measure times or distances
	double alpha =
    characteristic_length_conservative/characteristic_length_reactive;
	double mu = characteristic_length_reactive / mean_advection;
	double characteristic_val = 0.;
	if (typeid(Mean_tag) == typeid(streamtube::Finite_tag))
		characteristic_val = (1.+alpha)/reaction_rate*particles_characteristic
      /particles_immobile_each;
	if (typeid(Mean_tag) == typeid(streamtube::Infinite_tag))
		characteristic_val = alpha * mu;
	if (typeid(Evolution_tag) == typeid(streamtube::Space_tag))
		characteristic_val *= mean_advection;

  //  Measure times or distances
	std::valarray<double> measure_points;
	if (typeid(Mean_tag) == typeid(streamtube::Finite_tag))
		measure_points =
      range::linspace<std::valarray<double>>(measure_min, measure_max, nr_measures);
	if (typeid(Mean_tag) == typeid(streamtube::Infinite_tag))
		measure_points =
      range::logspace<std::valarray<double>>(measure_min, measure_max, nr_measures);
	measure_points *= characteristic_val;

  // Setup dynamics
	using Mass = std::size_t;
	using ImmobileSpecies =
    useful::StoreConst<std::vector<Mass>, std::vector<Mass> const&>;
	using MobileSpecies = streamtube::Species_initial<Mass>;
	using PatchGenerator = streamtube::PatchGenerator_alternating
    <Length_reactive, Length_conservative, ImmobileSpecies, Mass>;
	using Reactor = decltype(gillespie::make_Gillespie_MassAction(
    std::vector<std::size_t>(types), 0., stoichiometry));
	using StreamTubeDynamics =
    streamtube::StreamTubeDynamics<PatchGenerator, Advection, Reactor, Mass>;
	using Evolver = streamtube::Evolver<StreamTubeDynamics, Evolution_tag>;
	AdvectionGenerator advection_generator =
    make_AdvectionGenerator(mean_advection, var_advection);

  //  Dynamics
	streamtube::Measurer<Evolution_tag> measurer{
    measure_points, nr_fixed_velocity, nr_velocities,
    particles_characteristic };
  //  Run each ensemble
	for (std::size_t streamtube = 0; streamtube < nr_velocities; ++streamtube)
	{
		std::cout << "velocity = " << streamtube << "\n";
		Advection advection{ advection_generator() };
    //  Multiple ensembles for each velocity
		for (std::size_t run = 0; run < nr_fixed_velocity; ++run)
		{
      std::cout << "\trun = " << run << "\n";
			StreamTubeDynamics streamtube_dynamics{
				{ make_LengthReactive(characteristic_length_reactive,
                              exp_length_reactive),
					make_LengthConservative(characteristic_length_conservative,
                                  exp_length_conservative),
					{ average_initial_immobile_particles } },
				advection,
				gillespie::make_Gillespie_MassAction(std::vector<std::size_t>(types),
                                             0., stoichiometry),
				MobileSpecies{ average_initial_mobile_particles,
          mean_advection }(advection(), flux_weighted) };
			for (std::size_t measure = 0; measure < measure_points.size(); ++measure)
			{
				Evolver::evolve(streamtube_dynamics,
                        measure_points[measure], tortuosity);
				measurer.collect(streamtube_dynamics,
                         measure, run + streamtube * nr_fixed_velocity);
			}
		}
	}

  //  Output
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
	stream << characteristic_length_reactive << "_"
         << exp_length_reactive << "_"
         << characteristic_length_conservative << "_"
         << exp_length_conservative << "_"
         << mean_advection << "_"
         << var_advection << "_"
         << reaction_rate << "_"
         << measure_min << "_"
         << measure_max << "_"
         << nr_measures << "_"
         << flux_weighted << "_"
         << particles_mobile_each << "_"
         << particles_immobile_each << "_"
         << nr_fixed_velocity << "_"
         << nr_velocities << "_"
         << run_nr;
  std::string filename_params = stream.str();
  
  std::string filename_mass{ output_dir + "/" +
    measurer.filename_base + "_" + filename_model + "_"
    + streamtube::Evolution_filename<Evolution_tag>{}.filename + "_"
    + filename_params + ".dat"};
  std::ofstream output{ filename_mass };
  if (!output.is_open())
    throw useful::open_write_error(filename_mass);
  output << std::scientific << std::setprecision(8);
  measurer.normalize();
	measurer(output);
  output.close();

	return 0;
}
