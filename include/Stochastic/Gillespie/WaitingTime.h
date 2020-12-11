//
//  WaitingTime.h
//  Stochastic
//
//  Created by Tomas Aquino on 1/31/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

//  Inter-reaction waiting times
//  WaitingTime classes must implement an operator with the signature
//  template <typename Container>
//  double operator() (Container const& rates, std::size_t reaction)
//  which return the waiting time given the rates and the reaction

#ifndef WaitingTime_h
#define WaitingTime_h

#include <cmath>
#include <vector>
#include <random>
#include "general/Operations.h"

namespace gillespie
{
  //  Standard Gillespie exponential waiting time
	class WaitingTime_Exponential
	{
	public:
		template <typename Container>
		double operator() (Container const& rates, std::size_t reaction = 0)
		{ return dist(rng)/operation::sum(rates); }
    
	private:
		std::exponential_distribution<double> dist{ 1. };
		std::mt19937 rng{ std::random_device{}() };
	};
}

#endif /* WaitingTime_h */
