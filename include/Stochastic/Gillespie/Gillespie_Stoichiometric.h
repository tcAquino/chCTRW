//
//  Gillespie_Stoichiometric.h
//  Stochastic
//
//  Created by Tomas Aquino on 1/31/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

//  High-level helpers to built instances of generalized Gillespie
//  algorithm handlers for mass-action reactions

#ifndef Gillespie_Stoichiometric_h
#define Gillespie_Stoichiometric_h

#include <array>
#include "Stochastic/Reaction.h"
#include "WaitingTime.h"
#include "DelayTime.h"
#include "Gillespie.h"

namespace gillespie
{
  //  Make a Gillespie for mass action reactions with overall delay
  template <typename DelayTime, typename... Stoichiometry>
  auto make_Gillespie_MassAction_Delay(std::vector<std::size_t> numbers, double time, DelayTime delay_time, Stoichiometry&&... stoichiometry)
  {
    return
    Gillespie<WaitingTime_Exponential, DelayTime, decltype(stochastic::Reaction_MassAction{ stoichiometry })...>
    { numbers, time, {}, delay_time,
      stochastic::Reaction_MassAction{ std::forward<Stoichiometry>(stoichiometry) }... };
  }

  //  Make a Gillespie for mass action reactions with overall delay
  //  Start time at 0.
  template <typename DelayTime, typename... Stoichiometry>
  auto make_Gillespie_MassAction_Delay
  (std::vector<std::size_t> numbers, DelayTime delay_time, Stoichiometry&&... stoichiometry)
  {
    return
    make_Gillespie_MassAction_Delay
    (numbers, 0., delay_time, std::forward<Stoichiometry>(stoichiometry)...);
  }

  //  Make a Gillespie for regular mass action reactions
  template <typename... Stoichiometry>
  auto make_Gillespie_MassAction(std::vector<std::size_t> numbers, double time, Stoichiometry&&... stoichiometry)
  {
    return
    make_Gillespie_MassAction_Delay
    (numbers, time, stochastic::DelayTime_NoDelay{}, std::forward<Stoichiometry>(stoichiometry)...);
  }

  //  Make a Gillespie for regular mass action reactions
  //  Start time at 0.
  template <typename... Stoichiometry>
  auto make_Gillespie_MassAction(std::vector<std::size_t> numbers, Stoichiometry&&... stoichiometry)
  {
    return
    make_Gillespie_MassAction
    (numbers, 0., std::forward<Stoichiometry>(stoichiometry)...);
  }
}


#endif /* Gillespie_Stoichiometric_h */

