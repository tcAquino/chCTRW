//
//  Stoichiometry.h
//  Stochastic
//
//  Created by Tomas Aquino on 1/31/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

//  Stoichiometry classes to use with reactions
//  Handle reaction properties such as rate parameters,
//  reactants, and products
//  Stoichiometry classes must have a visible ReactantStoichiometry type

#ifndef Stoichiometry_h
#define Stoichiometry_h

#include <unordered_map>
#include <utility>
#include <vector>

namespace stochastic
{
  //  Standard stoichiometry class
  //  Keeps a reaction rate, reactants, and products
  //  The reactants and products are encoded as vectors of pairs
  //  E.g. reactants = { { 0, 1 }, { 2, 3 } } and products = { { 0, 2 }, { 3, 1 } }
  //  corresponds to the reaction A + 3B -> 2A + C
	class Stoichiometry
	{
	public:
		using ReactantStoichiometry = std::vector<std::pair<std::size_t,std::size_t>>;
    const double reaction_rate;       // Pure (state-independent) reaction rate
    ReactantStoichiometry reactants;  /* Reaction stoichiometry for reactants
                                       * Each pair refers to a species and coefficient */
    ReactantStoichiometry products;   /* Reaction stoichiometry for products
                                       * Each pair refers to a species and coefficient */

		Stoichiometry
		(double reaction_rate, std::vector<std::pair<std::size_t, std::size_t>> const& reactants, std::vector<std::pair<std::size_t, std::size_t>> const& products)
		: reaction_rate(reaction_rate)
		, reactants(reactants)
		, products(products)
		{
			for (auto const& rea : reactants)
        reactants_map[rea.first] = rea.second;
			for (auto const& pro : products)
        products_map[pro.first] = pro.second;
		}

    //  Stoichiometric coefficient associated with a given reactant
		std::size_t reactant_coefficient(std::size_t reactant) const
		{ return reactants_map.at(reactant); }

    //  Stoichiometric coefficient associated with a given product
		std::size_t product_coefficient(std::size_t product) const
		{ return products_map.at(product); }

	private:
		std::unordered_map<std::size_t, std::size_t> reactants_map;
		std::unordered_map<std::size_t, std::size_t> products_map;
	};
}

#endif /* Stoichiometry_h */
