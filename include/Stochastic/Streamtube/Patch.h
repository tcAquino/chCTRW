//
//  Patch.h
//  Streamtube
//
//  Created by Tomas Aquino on 4/4/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

//  Patch generators for streamtube models

#ifndef Patch_h
#define Patch_h

namespace streamtube
{
  //  New patch alternatingly reactive and non-reactive
  //  Reactive_length must implement double operator()() generating a reactive patch length
  //  Conservative_length must implement double operator()() generating a conservative patch length
  //  Particle_generator must implement std::vector<Mass> operator()()
  //  generating masses of each type in a reactive patche
  //  Mass is the mass type, typically std::size_t or double
	template <typename Reactive_length, typename Conservative_length, typename Particle_generator, typename Mass>
	class PatchGenerator_alternating
	{
	public:
		PatchGenerator_alternating
		(Reactive_length reactive_length, Conservative_length conservative_length, Particle_generator particle_generator)
		: reactive_length{ reactive_length }
		, conservative_length{ conservative_length }
		, particle_generator{ particle_generator }
		{ generate(); }

    //  Generate a new patch
		void generate()
		{
			current_reactive = !current_reactive;
			if (reactive())
			{
				current_length = reactive_length();
				current_particles = particle_generator();
			}
			else
				current_length = conservative_length();
		}

		void mass(std::size_t type, Mass particles)
		{ current_particles[type] = particles; }

		double length() const
		{ return current_length; }

		bool reactive() const
		{ return current_reactive; }

		Mass mass(std::size_t type) const
		{ return current_particles[type]; }

		std::size_t nr_types() const
		{ return current_particles.size(); }

	private:
		Reactive_length reactive_length;
		Conservative_length conservative_length;
		Particle_generator particle_generator;
		std::vector<Mass> current_particles;
		double current_length;
		bool current_reactive{ 0 };
	};
}

#endif /* Patch_h */
