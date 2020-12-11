//
//  streamtube_concentration.cpp
//
//  Created by Tomas Aquino on 5/22/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <typeinfo>
#include <valarray>
#include "general/useful.h"
#include "general/Ranges.h"
#include "general/useful.h"
#include "Stochastic/Reaction.h"
#include "Stochastic/Streamtube/Models.h"
#include "Stochastic/Streamtube/Measurer.h"

int main(int argc, const char * argv[])
{
  if (argc == 0)
  {
    std::cout << "streamtube_concentration\n";
    std::cout << "Parameters (default value in []):\n"
              << "length_reactive : Characteristic length of reactive patches\n"
              << "alpha : Ratio between characterstic conservative and reactive lengths\n"
              << "beta : Exponent for stable conservative patches (ignored for exponential)\n"
              << "mean_advection : Mean advection across streamtubes\n"
              << "var_advection : Variance of advection across streamtubes (ignored for one-parameters dists)\n"
              << "reaction_rate : Macroscopic reaction rate\n"
              << "measure_min : Minimum time or distance for output\n"
              << "measure_max : Maximum time or distance for output\n"
              << "nr_measures : Number of outputs\n"
              << "c01 : Initial concentration of first species\n"
              << "c02 : Initial concentration of second species\n"
              << "flux_weighted : 0 - Homogeneous injection\n"
              << "                1 - Flux-weighted injection\n"
              << "dist : 0 - Measure average mass only\n"
              << "       1 - Measure average mass and mass distribution across particles\n"
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
  double length_reactive = atof(argv[arg++]);
  double alpha = atof(argv[arg++]);
  double beta = atof(argv[arg++]);
  double mean_advection = atof(argv[arg++]);
  double var_advection = atof(argv[arg++]);
  double reaction_rate = atof(argv[arg++]);
  double measure_min = atof(argv[arg++]);
  double measure_max = atof(argv[arg++]);
  std::size_t nr_measures = strtoul(argv[arg++], NULL, 0);
  double c01 = atof(argv[arg++]);
  double c02 = atof(argv[arg++]);
  bool flux_weighted = atoi(argv[arg++]);
  bool dist = atoi(argv[arg++]);
  std::size_t nr_fixed_velocity = strtoul(argv[arg++], NULL, 0);
  std::size_t nr_velocities = strtoul(argv[arg++], NULL, 0);
  std::size_t run_nr = strtoul(argv[arg++], NULL, 0);
  std::string output_dir = argc > arg ? argv[arg++] : "../output";
  
  double tortuosity = 1.;

  //  Nondimensionalization of output times or distances
  double mu = length_reactive/mean_advection;
  double characteristic_val = 0.;
  if (typeid(Mean_tag) == typeid(streamtube::Finite_tag))
    characteristic_val = (1. + alpha)/(reaction_rate * c02);
  if (typeid(Mean_tag) == typeid(streamtube::Infinite_tag))
    characteristic_val = (alpha * mu)/std::pow(mu * reaction_rate * c02, 1./beta);
  if (typeid(Evolution_tag) == typeid(streamtube::Space_tag))
    characteristic_val *= mean_advection;

  //  Output times or distances
  std::valarray< double > measure_points;
  if (typeid(Mean_tag) == typeid(streamtube::Finite_tag))
    measure_points = range::linspace< std::valarray< double > >(measure_min, measure_max, nr_measures);
  if (typeid(Mean_tag) == typeid(streamtube::Infinite_tag))
    measure_points = range::logspace< std::valarray< double > >(measure_min, measure_max, nr_measures);
  measure_points *= characteristic_val;
  nr_measures = measure_points.size();

  //  Setup dynamics
  using Mass = double;
  using ImmobileSpecies = useful::StoreConst<std::vector<Mass>, std::vector<Mass> const&>;
  using MobileSpecies = streamtube::Species_initial<Mass>;
  using PatchGenerator = streamtube::PatchGenerator_alternating
    <Length_reactive, Length_conservative, ImmobileSpecies, Mass>;
  using Reactor = stochastic::Reaction_concentration_bimolecular_analytical;
  using StreamTubeDynamics = streamtube::StreamTubeDynamics
    <PatchGenerator, Advection, Reactor, Mass>;
  using Evolver = streamtube::Evolver<StreamTubeDynamics, Evolution_tag>;
  AdvectionGenerator advection_generator = make_AdvectionGenerator(mean_advection, var_advection);

  //  Dynamics
  streamtube::Measurer<Evolution_tag> measurer{
    measure_points, nr_fixed_velocity, nr_velocities, 1., dist };
  for (std::size_t streamtube = 0; streamtube < nr_velocities; ++streamtube)
  {
    Advection advection{ advection_generator() };
    printf("velocity = %zu %.2e\n", streamtube, advection());
    for (std::size_t run = 0; run < nr_fixed_velocity; ++run)
    {
      printf("\trun = %zu\n", run);
      StreamTubeDynamics streamtube_dynamics{
        { make_LengthReactive(length_reactive),
          make_LengthConservative(alpha * length_reactive, beta),
          { { c02 } } },
        advection,
        { reaction_rate },
        MobileSpecies{ { c01 }, mean_advection }(advection(), flux_weighted) };
      for (std::size_t measure = 0; measure < measure_points.size(); ++measure)
      {
        Evolver::evolve(streamtube_dynamics, measure_points[measure], tortuosity);
        measurer.collect(streamtube_dynamics, measure, streamtube);
      }
    }
  }

  //  Output
  std::stringstream stream;
  stream << std::scientific << std::setprecision(2);
  stream << length_reactive << "_"
         << alpha << "_"
         << beta << "_"
         << mean_advection << "_"
         << var_advection << "_"
         << reaction_rate << "_"
         << c01 << "_"
         << c02 << "_"
         << measure_min << "_"
         << measure_max << "_"
         << nr_measures << "_"
         << flux_weighted << "_"
         << nr_fixed_velocity << "_"
         << nr_velocities << "_"
         << run_nr;
  std::string filename_params = stream.str();
  
  std::string filename_mass{ output_dir + "/" +
    measurer.filename_base + "_concentration_"
    + measurer.filename_base + "_" + filename_model + "_"
    + streamtube::Evolution_filename<Evolution_tag>{}.filename + "_"
    + filename_params + ".dat" };
  std::ofstream output_mass{ filename_mass };
  if (!output_mass.is_open())
    throw useful::open_write_error(filename_mass);
  output_mass << std::scientific << std::setprecision(8);
  
  std::string filename_dist{ output_dir + "/" +
    measurer.filename_base + "_dist_"
    + measurer.filename_base + "_" + filename_model + "_"
    + streamtube::Evolution_filename<Evolution_tag>{}.filename + "_"
    + filename_params + ".dat" };
  std::ofstream output_dist{ filename_dist };
  if (!output_dist.is_open())
    throw useful::open_write_error(filename_mass);
  output_dist << std::scientific << std::setprecision(8);
    
  measurer.normalize();
  measurer(output_mass, output_dist);
  output_mass.close();
  output_dist.close();
  
  return 0;
}

