//
//  Models.h
//  Streamtube
//
//  Created by Tomas Aquino on 4/3/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

//  Type and helper function definitions for different streamtube models
//  namespace naming convention is model_<velocity_dist>_<reactive_length_dist>_<conservative_length_dist>

#ifndef Models_Streamtube_h
#define Models_Streamtube_h

#include <cmath>
#include <cstdio>
#include <valarray>
#include "Streamtube.h"
#include "Patch.h"
#include "Stochastic/Random.h"
#include "general/useful.h"

namespace streamtube
{
  // Tags for compile-time feature identification
	struct Space_tag{};
	struct Time_tag{};
	struct Finite_tag{};
	struct Infinite_tag{};

  // Filename modifiers for output
	template <typename Evolution_tag>
	struct Evolution_filename;
	template <>
	struct Evolution_filename< Space_tag >
	{ char filename[32] = { "space" }; };
	template <>
	struct Evolution_filename< Time_tag >
	{ char filename[32] = { "time" }; };

  // Dynamics in space or time
	template <typename StreamTubeDynamics, typename Evolution_tag>
	struct Evolver;
	template <typename StreamTubeDynamics>
	struct Evolver<StreamTubeDynamics, Space_tag>
	{
		static void evolve(StreamTubeDynamics& streamtube_dynamics, double final_val, double tortuosity)
		{ streamtube_dynamics.evolve_position(tortuosity * final_val); }
	};
	template <typename StreamTubeDynamics>
	struct Evolver<StreamTubeDynamics, Time_tag>
	{
		static void evolve(StreamTubeDynamics& streamtube_dynamics, double final_val, double = 0.)
		{ streamtube_dynamics.evolve_time(final_val); }
	};

  //  Initial condition
	template <typename Mass>
	struct Species_initial
	{
		std::vector<Mass> operator()(double advection, bool flux_weighted)
		{
			if (flux_weighted)
			{
				std::vector<Mass> particles_fw;
				particles_fw.reserve(particles.size());
				for (auto const& part : particles)
					particles_fw.push_back(advection*part/mean_advection);
				return particles_fw;
			}
			else
				return particles;
		}
		const std::vector< Mass > particles;
		const double mean_advection;
	};

	struct Advection_uniform
	{
		const double advection;
		void generate(){};
		double operator()()
		{ return advection; }
	};

	namespace model_uniform_exp_exp
	{
		using Evolution_tag = Time_tag;
		using Mean_tag = Finite_tag;
		using Advection = Advection_uniform;
		using Tortuosity = useful::StoreConst<double>;
		using Length_reactive = stochastic::RNG<std::exponential_distribution<double>>;
		using Length_conservative = stochastic::RNG<std::exponential_distribution<double>>;
		using AdvectionGenerator = useful::StoreConst<double>;
		AdvectionGenerator make_AdvectionGenerator(double advection, double = 0.)
		{ return AdvectionGenerator{ advection }; }
		Length_reactive make_LengthReactive(double length, double = 0.)
		{ return Length_reactive{ 1./length }; }
		Length_conservative make_LengthConservative(double length, double = 0.)
		{ return Length_conservative{ 1./length }; }
		char filename_model[32] = { "uniform_exp_exp" };
	}

	namespace model_uniform_exp_power
	{
		using Evolution_tag = Time_tag;
		using Mean_tag = Infinite_tag;
		using Advection = Advection_uniform;
		using Tortuosity = useful::StoreConst<double>;
		using Length_reactive = stochastic::RNG<std::exponential_distribution<double>>;
		using Length_conservative = stochastic::RNG<stochastic::skewedlevystable_distribution<double>>;
		using AdvectionGenerator = useful::StoreConst<double>;
		AdvectionGenerator make_AdvectionGenerator(double advection, double = 0.)
		{ return AdvectionGenerator{ advection }; }
		Length_reactive make_LengthReactive(double length, double = 0.)
		{ return Length_reactive{ 1./length }; }
		Length_conservative make_LengthConservative(double length, double alpha)
		{
      return Length_conservative{
      Length_conservative::param_type{
        alpha,
        std::pow(std::cos(constants::pi * alpha / 2.), 1. / alpha) * length } };
    }
		char filename_model[32] = { "uniform_exp_power" };
	}

	namespace model_uniform_uniform_uniform
	{
		using Evolution_tag = Time_tag;
		using Mean_tag = Finite_tag;
		using Advection = Advection_uniform;
		using Tortuosity = useful::StoreConst<double>;
		using Length_reactive = useful::StoreConst<double>;
		using Length_conservative = useful::StoreConst<double>;
		using AdvectionGenerator = useful::StoreConst<double>;
		AdvectionGenerator make_AdvectionGenerator(double advection, double = 0.)
		{ return AdvectionGenerator{ advection }; }
		Length_reactive make_LengthReactive(double length, double = 0.)
		{ return Length_reactive{ 1. / length }; }
		Length_conservative make_LengthConservative(double length, double = 0.)
		{ return Length_conservative{ 1. / length }; }
		char filename_model[32] = { "uniform_uniform_uniform" };
	}

	namespace model_gamma_exp_exp
	{
		using Evolution_tag = Space_tag;
		using Mean_tag = Finite_tag;
		using Advection = Advection_uniform;
		using Tortuosity = useful::StoreConst<double>;
		using Length_reactive = stochastic::RNG<std::exponential_distribution<double>>;
		using Length_conservative = stochastic::RNG<std::exponential_distribution<double>>;
		using AdvectionGenerator = stochastic::RNG<std::gamma_distribution<double>>;
		AdvectionGenerator make_AdvectionGenerator(double mean, double var)
		{ return AdvectionGenerator{
        AdvectionGenerator::param_type{ mean * mean / var, var / mean } }; }
		Length_reactive make_LengthReactive(double length, double = 0.)
		{ return Length_reactive{ 1. / length }; }
		Length_conservative make_LengthConservative(double length, double = 0.)
		{ return Length_conservative{ 1. / length }; }
		char filename_model[32] = { "gamma_exp_exp" };
	}

	namespace model_gamma_exp_power
	{
		using Evolution_tag = Space_tag;
		using Mean_tag = Infinite_tag;
		using Advection = Advection_uniform;
		using Tortuosity = useful::StoreConst<double>;
		using Length_reactive = stochastic::RNG<std::exponential_distribution<double>>;
		using Length_conservative = stochastic::RNG< stochastic::skewedlevystable_distribution<double> >;
		using AdvectionGenerator = stochastic::RNG<std::gamma_distribution<double>>;
		AdvectionGenerator make_AdvectionGenerator(double mean, double var)
		{ return AdvectionGenerator{
        AdvectionGenerator::param_type{ mean * mean / var, var / mean } }; }
		Length_reactive make_LengthReactive(double length, double = 0.)
		{ return Length_reactive{ 1. / length }; }
		Length_conservative make_LengthConservative(double length, double alpha)
		{ return Length_conservative{
        Length_conservative::param_type{ alpha,
          std::pow(std::cos(constants::pi * alpha / 2.), 1. / alpha) * length } }; }
		char filename_model[32] = { "gamma_exp_power" };
	}
}


#endif /* Models_Streamtube_h */
