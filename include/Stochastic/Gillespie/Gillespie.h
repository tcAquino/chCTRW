//
//  Gillespie.h
//  Stochastic
//
//  Created by Tomas Aquino on 1/31/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef Gillespie_h
#define Gillespie_h

#include <tuple>
#include <limits>
#include <vector>
#include <random>
#include <utility>
#include "general/useful.h"

//  For use with Gillespie algorithm, reaction handler classes should implement:
//  double rate(std::vector<std::size_t> const& numbers) const;
//  void operator()(std::vector<std::size_t>& concentration, double time_step, double time);
//  A visible type Stoichiometry

namespace gillespie
{
  template<typename WaitingTime, typename DelayTime, typename... Reactions>
  class Gillespie
  {
  public:
    using ReactantStoichiometry = typename std::tuple_element<0, std::tuple<Reactions...>>::type::Stoichiometry::ReactantStoichiometry;
    using Part_Container = std::vector<std::size_t>;

    Gillespie(Part_Container particles, double time, WaitingTime waiting_time, DelayTime delay_time, Reactions... reactions)
    : particle_container(particles)
    , time_current(time)
    , waiting_time(waiting_time)
    , delay_time(delay_time)
    , reactions(reactions...)
    , reaction_table{ make_reaction_table() }
    , reactant_table{ make_reactant_table() }
    , product_table{ make_product_table() }
    {}

    void set(Part_Container const& particles, double time = 0.)
    {
      particle_container = particles;
      time_current = time;
    }

    //  Set particle numbers of a type
    void set(std::size_t type, std::size_t particle_nr)
    { particle_container[type] = particle_nr; }

    //  Set particle numbers of designated types
    void set(std::vector<std::size_t> const& types, Part_Container const& particles)
    {
      for (std::size_t type = 0; type <particles.size(); ++type)
        set(types[type], particles[type]);
    }

    void time(double time)
    { time_current = time; }

    // Remove all particles
    void clear()
    { std::fill(particle_container.begin(), particle_container.end(), 0); }

    void add(std::size_t type, std::size_t increment = 1)
    { particle_container[type] += increment; }

    void remove(std::size_t type, std::size_t increment = 1)
    {
      particle_container[type] > increment
      ? particle_container[type] -= increment
      : 0;
    }

    //  Update state to just after next reaction
    void evolve()
    {
      reacted = 0;
      rates();
      if (*std::max_element(rate_container.begin(), rate_container.end()) == 0.)
        time_next_reaction = std::numeric_limits<double>::infinity();
      else
      {
        pick_reaction();
        compute_time_next_reaction();
        react(next_reaction);
        reacted = 1;
      }
      time_current = time_next_reaction;
    }

    //  Update state to time_max
    //  A record of the next reaction time and reaction is kept
    void evolve(double time_max)
    {
      reacted = 0;
      while (1)
      {
        // Compute rates
        rates();
        // If all rates are zero
        if (*std::max_element(rate_container.begin(), rate_container.end()) == 0.)
          time_next_reaction = std::numeric_limits<double>::infinity();
        else
        {
          pick_reaction();
          compute_time_next_reaction();
        }
        if (time_next_reaction <time_max)
        {
          time_current = time_next_reaction;
          react(next_reaction);
          reacted = 1;
        }
        else
        {
          time_current = time_max;
          break;
        }
      }
    }

    double rate_sum() const
    {
      rates();
      return operation::sum(rate_container);
    }

    double time() const
    { return time_current; }

    double time_last() const
    { return time_last_reaction; }

    double time_next() const
    { return time_next_reaction; }

    double last()
    { return last_reaction; }

    double next()
    { return next_reaction; }

    bool reaction() const
    { return reacted; }

    Part_Container const& particles() const
    { return particle_container; }

    std::size_t particles(std::size_t nr) const
    { return particle_container[ nr ]; }

    std::size_t nr_types() const
    { return particle_container.size(); }

    ReactantStoichiometry const& reactants(std::size_t reaction)
    { return reactant_table[reaction](); }

    ReactantStoichiometry const& products(std::size_t reaction)
    { return product_table[reaction](); }

  private:

    // Auxiliary types for runtime dispatch implementations
    using function_type_reaction = void (*)(std::tuple<Reactions...> const&, std::vector<std::size_t>&);
    using function_array_reaction = std::array<function_type_reaction, sizeof...(Reactions)>;
    using function_type_stoichiometry = ReactantStoichiometry const& (*)(std::tuple<Reactions...> const& reactions);
    using function_array_stoichiometry = std::array<function_type_stoichiometry, sizeof...(Reactions)>;
    using array_type = std::array<double, sizeof...(Reactions)>;

    mutable std::mt19937 rng{ std::random_device{}() };  // RNG
    Part_Container particle_container;                  //Numbers of particles of each type
    double time_current;
    WaitingTime waiting_time;                           // Intrinsic inter-reaction time
    DelayTime delay_time;                               // Overall delay
    std::tuple<Reactions...> reactions;
    
    // Runtime dispatches
    const function_array_reaction reaction_table;
    const function_array_stoichiometry reactant_table;
    const function_array_stoichiometry product_table;
    
    double time_last_reaction;
    double time_next_reaction;
    std::size_t last_reaction;
    std::size_t next_reaction;
    bool reacted = 0;                   // True if reacted during the last evolution

    mutable array_type rate_container;  // State-dependent rates for each reaction

    // Compile-time check if there is more than one reaction
    constexpr static bool more_than_one_reaction{ bool(std::minus<std::size_t>{}(sizeof...(Reactions), 1)) };

    void react(std::size_t index)
    {
      last_reaction = next_reaction;
      time_last_reaction = time_next_reaction;
      reaction_table[index](reactions, particle_container);
    }

    //  Cumulative reaction rates based on current state
    void rates() const
    {
      rates_impl(useful::Selector<bool, more_than_one_reaction>{});
    }

    void pick_reaction()
    {
      next_reaction = pick_reaction_impl(useful::Selector<bool, more_than_one_reaction>{});
    }

    void compute_time_next_reaction()
    {
      double waiting = waiting_time(rate_container);
      time_next_reaction = time_current + waiting + delay_time(waiting);
    }

    // If there is more than one reaction
    void rates_impl(useful::Selector<bool, 1>) const
    {
      std::size_t ii = 0;
      useful::for_each(reactions, [ &ii, this ](auto& reaction)
      { rate_container[ ii++ ] = reaction.rate(particle_container); });
    }
    
    // If there is only one reaction
    void rates_impl(useful::Selector<bool, 0>) const
    { rate_container[0] = std::get<0>(reactions).rate(particle_container); }
    
    // If there is more than one reaction
    std::size_t pick_reaction_impl(useful::Selector<bool, 1>) const
    {
      return std::discrete_distribution<std::size_t>{
        rate_container.begin(), rate_container.end() }(rng);
    }

    // If there is more than one reaction
    std::size_t pick_reaction_impl(useful::Selector<bool, 0>) const
    { return 0; }

    function_array_reaction make_reaction_table()
    {
      return make_reaction_table(std::make_index_sequence<sizeof...(Reactions)>{});
    }

    function_array_stoichiometry make_reactant_table()
    {
      return make_reactant_table(std::make_index_sequence<sizeof...(Reactions)>{});
    }

    function_array_stoichiometry make_product_table()
    {
      return make_product_table(std::make_index_sequence<sizeof...(Reactions)>{});
    }

    //  Runtime dispatch to execute reactions
    template<std::size_t... Indices>
    function_array_reaction make_reaction_table(std::index_sequence<Indices...>)
    {
      return {
        { [](std::tuple<Reactions...> const& reactions, std::vector<std::size_t>& particles){
          std::get<Indices>(reactions).react(particles); }... } };
    }

    // Runtime dispatch to access reactants
    template<std::size_t... Indices>
    function_array_stoichiometry make_reactant_table(std::index_sequence<Indices...>)
    {
      return {
        { [](std::tuple<Reactions...> const& reactions) -> auto const& {
          return std::get<Indices>(reactions).stoichiometry.reactants; }... } };
    }

    //  Runtime dispatch to access products
    template<std::size_t... Indices>
    function_array_stoichiometry make_product_table(std::index_sequence<Indices...>)
    {
      return {
        { [](std::tuple<Reactions...> const& reactions) -> auto const& {
          return std::get<Indices>(reactions).stoichiometry.products; }... } };
    }
  };
}



#endif /* Gillespie_h */

