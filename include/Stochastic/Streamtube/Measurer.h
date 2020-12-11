//
//  Measurer.h
//  Streamtube
//
//  Created by Tomás Aquino on 13/03/2020.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

//  Measurer classes for streamtube models

#ifndef Measurer_Streamtube_h
#define Measurer_Streamtube_h

#include <string>
#include <valarray>
#include <vector>

namespace streamtube
{
  template <typename Evolution_tag>
  class Measurer;

  //  Measures average mass of each species and average product of masses as a function of time,
  //  and if dist = 1., measures particle positions and fixed-velocity mass averages
  template <>
  class Measurer<Time_tag>
  {
  public:
    std::string filename_base{ "Data_StreamTube_Mass" };

    Measurer(std::valarray<double> measure_times, std::size_t nr_runs,  std::size_t nr_streamtubes, double particles_characteristic, bool dist = 0)
    : measure_times{ measure_times }
    , nr_runs{ nr_runs }
    , nr_streamtubes{ nr_streamtubes }
    , particles_characteristic{ particles_characteristic }
    , average_of_mass_1(measure_times.size())
    , average_of_mass_2(measure_times.size())
    , average_of_product(measure_times.size())
    , average_of_mass_dist(measure_times.size(), std::valarray<double>(nr_streamtubes))
    , positions(measure_times.size(), std::valarray<double>(nr_streamtubes))
    , dist(dist)
    {}

    template < typename StreamTubeDynamics >
    void collect(StreamTubeDynamics const& streamtube_dynamics, std::size_t measure, std::size_t streamtube)
    {
      average_of_mass_1[measure] += streamtube_dynamics.mass(0);
      average_of_mass_2[measure] += streamtube_dynamics.mass_immobile(0);
      average_of_product[measure] += streamtube_dynamics.mass(0) * streamtube_dynamics.mass_immobile(0);
      if (dist)
      {
        average_of_mass_dist[measure][streamtube] += streamtube_dynamics.mass(0);
        positions[measure][streamtube] += streamtube_dynamics.position();
      }
    }
    
    void normalize()
    {
      average_of_mass_1 /= particles_characteristic*nr_runs*nr_streamtubes;
      average_of_mass_2 /= particles_characteristic*nr_runs*nr_streamtubes;
      average_of_product /= particles_characteristic*particles_characteristic*nr_runs*nr_streamtubes;
      
      if (dist)
        for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
        {
          average_of_mass_dist[tt] /= particles_characteristic*nr_runs;
          positions[tt] /= nr_runs;
        }
    }
    
    template <typename Stream>
    void operator()(Stream& output_mass, Stream& output_dist)
    {
      for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
      {
        output_mass << measure_times[tt] << "\t"
                    << average_of_mass_1[tt] << "\t"
                    << average_of_mass_2[tt] << "\t"
                    << average_of_product[tt] << "\n";
        if (dist)
        {
          output_dist << measure_times[tt] << "\t";
          for (std::size_t ss = 0; ss < nr_streamtubes; ++ss)
            output_dist << positions[tt][ss] << "\t"
                        << average_of_mass_dist[tt][ss];
          output_dist << "\n";
        }
      }
    }
    
    template <typename Stream>
    void operator()(Stream& output_mass)
    {
      for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
        output_mass << measure_times[tt] << "\t"
                    << average_of_mass_1[tt] << "\t"
                    << average_of_mass_2[tt] << "\t"
                    << average_of_product[tt] << "\n";
    }

  private:
    const std::valarray<double> measure_times;
    const std::size_t nr_runs;
    const std::size_t nr_streamtubes;
    const double particles_characteristic;
    std::valarray<double> average_of_mass_1;
    std::valarray<double> average_of_mass_2;
    std::valarray<double> average_of_product;
    std::vector<std::valarray<double>> average_of_mass_dist;
    std::vector<std::valarray<double>> positions;
    bool dist;
  };

  //  Measures average mass of first species as a function of space,
  //  and if dist = 1., measures crossing times and fixed-velocity mass averages
  template <>
  class Measurer<Space_tag>
  {
  public:
    const std::string filename_base{ "Data_StreamTube_BTC" };

    Measurer(std::valarray<double> measure_distances, std::size_t nr_runs, std::size_t nr_streamtubes, double particles_characteristic, bool dist)
    : measure_distances{ measure_distances }
    , nr_runs(nr_runs)
    , nr_streamtubes{ nr_streamtubes }
    , particles_characteristic{ particles_characteristic }
    , average_of_mass(measure_distances.size())
    , average_of_mass_dist(measure_distances.size(), std::valarray<double>(nr_streamtubes))
    , crossing_times(measure_distances.size(), std::valarray<double>(nr_streamtubes))
    , dist(dist)
    {}

    template <typename StreamTubeDynamics>
    void collect(StreamTubeDynamics const& streamtube_dynamics, std::size_t measure, std::size_t streamtube)
    {
      average_of_mass[measure] += streamtube_dynamics.mass(0);
      if (dist==1)
      {
        average_of_mass_dist[measure][streamtube] += streamtube_dynamics.mass(0);
        crossing_times[measure][streamtube] += streamtube_dynamics.time();
      }
    }
    
    void normalize()
    {
      average_of_mass /= particles_characteristic*nr_runs*nr_streamtubes;
      
      if (dist)
         for (std::size_t xx = 0; xx < measure_distances.size(); ++xx)
         {
           average_of_mass_dist[xx] /= particles_characteristic*nr_runs;
           crossing_times[xx] /= nr_runs;
         }
    }

    template <typename Stream>
    void operator()(Stream& output_mass, Stream& output_dist)
    {
      for (std::size_t xx = 0; xx < measure_distances.size(); ++xx)
      {
        output_mass << measure_distances[xx] << "\t"
                    << average_of_mass[xx] << "\n";
        if (dist)
        {
          output_mass << measure_distances[xx] << "\t";
          for (std::size_t ss = 0; ss < nr_streamtubes; ++ss)
            output_mass << crossing_times[xx][ss] << "\t"
                        << average_of_mass_dist[xx][ss] << "\n";
        }
      }
    }
    
    template <typename Stream>
    void operator()(Stream& output_mass)
    {
      for (std::size_t xx = 0; xx < measure_distances.size(); ++xx)
        output_mass << measure_distances[xx] << "\t"
                    << average_of_mass[xx] << "\n";
    }

  private:
    const std::valarray<double> measure_distances;
    const std::size_t nr_runs;
    const std::size_t nr_streamtubes;
    const double particles_characteristic;
    std::valarray<double> average_of_mass;
    std::vector<std::valarray<double>> average_of_mass_dist;
    std::vector<std::valarray<double>> crossing_times;
    bool dist;
  };
}

#endif /* Measurer_Streamtube_h */


