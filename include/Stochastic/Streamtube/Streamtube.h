//
//  Streamtube.h
//  Streamtube
//
//  Created by Tomas Aquino on 3/4/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

//  Transport and reaction dynamics in a streamtube model
//  For use with streamtube models, reactions handler classes implement:
//  void evolve(double time_max);
//  void set(std::size_t type, double val);
//  void time(double val);
//  double particles(std::size_t type);

#ifndef Streamtube_h
#define Streamtube_h

#include <vector>

namespace streamtube
{
  //  Patch generates the current patch
  //  Advection implements double()() returning advection in current streamtube
  //  Reactor handles reactions in reactive patches (e.g., Gillespie algorithm or a rate law)
  //  Mass is the mass type, typically std::size_t or double
	template <typename Patch, typename Advection, typename Reactor, typename Mass>
	class StreamTubeDynamics
	{
	public:
		StreamTubeDynamics
		(Patch patch, Advection advection, Reactor reactor,
		 std::vector< Mass > mass,
		 double position = 0., double time = 0.)
		: patch{ patch }
		, reactor{ reactor }
		, advection{ advection }
		, mass_mobile{ mass }
		, current_position{ position }
		, current_time{ time }
		{}

		void evolve_position(double final_position)
		{
      //  If the first patch is not the last one
			if (current_position + patch.length() - position_in_patch
          < final_position)
			{
        // React in current patch
				react(patch.length() - position_in_patch);
        // Prepare the next patch
				generate();
			}
      //  If the first patch is the last one
			else
			{
				double position_increment = final_position-current_position;
				react(position_increment);
				position_in_patch += position_increment;
				return;
			}
      //  While the current patch is not the last one
			while (current_position + patch.length()
             < final_position)
			{
				react(patch.length());
				generate();
			}
      //  Final patch if it was not the first one
			double position_increment = final_position-current_position;
			position_in_patch = position_increment;
			react(position_increment);
		}

		void evolve_time(double final_time)
		{ evolve_position(advection() * final_time); }

		auto const& particles() const
		{ return mass_mobile; }

		Mass mass(std::size_t type) const
		{ return mass_mobile[ type ]; }

		Mass mass_immobile(std::size_t type) const
		{ return patch.mass(type); }

		double time() const
		{ return current_time; }

		double position() const
		{ return current_position; }

	private:
		Patch patch;
		Reactor reactor;
		Advection advection;
		std::vector<Mass> mass_mobile;
		double current_position;
		double current_time;
		double position_in_patch{ 0. };
		std::size_t nr_types_mobile{ mass_mobile.size() };
		std::size_t nr_types_immobile{ patch.nr_types() };
		std::size_t nr_types{ nr_types_mobile + nr_types_immobile };

    //  Set the numbers of particles and time in a reactive patch
		void set_reactor()
		{
			for (std::size_t type = 0; type < nr_types_mobile; ++type)
				reactor.set(type, mass_mobile[type]);
			for (std::size_t type = 0; type < nr_types_immobile; ++type)
				reactor.set(nr_types_mobile + type, patch.mass(type));
			reactor.time(current_time);
		}

    //  Set the numbers of particles in the overall state
		void set_state()
		{
			for (std::size_t type = 0; type < nr_types_mobile; ++type)
				mass_mobile[type] = reactor.particles(type);
			for (std::size_t type = 0; type < nr_types_immobile; ++type)
				patch.mass(type, reactor.particles(nr_types_mobile + type));
		}

    // React within current patch
		void react(double position_increment)
		{
			double time_increment = position_increment/advection();
			if (patch.reactive())
			{
				set_reactor();
				reactor.evolve(current_time + time_increment);
				set_state();
			}
			current_position += position_increment;
			current_time += time_increment;
		}

    //  Generate the next patch
		void generate()
		{
			patch.generate();
			advection.generate();
		}
	};
}

#endif /* Streamtube_h */
