/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of single-spin simulated annealing algorithm for
Ising spin glasses with general interactions with magnetic field
and any number of neighbors.

---------------------------------------------------------------------

Copyright (C) 2013 by Ilia Zintchenko <zintchenko@itp.phys.ethz.ch>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <cmath>
#include <random>
#include <vector>
#include <string>
#include <functional>

#include "lattice.h"

#define OMP_VERSION_2

template<typename T = uint64_t>
  class Algorithm
  {
  public:
	
  typedef double value_type;
  typedef unsigned index_type;

  static const std::size_t word_size = 1;

  struct site_type{
    int spin;
    value_type hzv;
    std::vector<value_type> jzv;
    value_type de;
    index_type nneighbs;
    std::vector<index_type> neighbs;
  };

  typedef Lattice<value_type, index_type> lattice_type;

  Algorithm() {}

  template <typename SE>
  Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
  : generator(41)
  {
    lattice.init_sites(sites);
    // std::cout << "sites.size() = " << sites.size() << std::endl;
    // std::cout << "Print the spin state ... " << std::endl;
    // for(auto& site : sites) std::cout << site.spin << " ";  // print the spin state
    // std::cout << std::endl;

    rng_ = std::bind(std::uniform_real_distribution<double>(0, 1), std::ref(generator));
    // std::cout << "rng_ = " << rng_() << std::endl;

    bound_array.resize(sched0.size());

    auto ba = bound_array.begin();

    for(const auto& s : sched0){

      ba->resize(sites.size());
      for(auto& a : *ba)
        a = -std::log(rng_()) / (s.beta * 2);

      ++ba;
    }

  }

  void reset_sites(const std::size_t rep)
  {
    // generator.seed(rep+1);
    generator.seed(time(nullptr)+rep+1);  //  Rand depends on timestamp
    

    for(auto& site : sites)
      site.spin = 2 * ((generator() >> 29) & 1) - 1;  // rand assigned as {-1,1}

    for(auto& site : sites){
      value_type tmp = site.hzv;
      for(index_type k = 0; k < site.nneighbs; ++k)
	      tmp += site.jzv[k] * sites[site.neighbs[k]].spin;
        site.de = -tmp * site.spin;
    }
  }  

  void flip_spin(site_type& site)
  {
    site.spin = -site.spin;
    site.de = -site.de;

    for (index_type k = 0; k<site.nneighbs; ++k) {
      site_type& neighbor = sites[site.neighbs[k]];
      neighbor.de -= 2 * neighbor.spin * site.jzv[k] * site.spin;
    }    
  }

  void do_sweep(const std::size_t sweep)
  {
    const std::size_t l = generator() % sites.size();
    const auto& ba = bound_array[sweep];

    for(std::size_t i = 0; i<l; ++i)
      if(sites[i].de<  ba[i + sites.size() - l])
        flip_spin(sites[i]);

    for(std::size_t i = l; i<sites.size(); ++i)
      if(sites[i].de < ba[i - l])
        flip_spin(sites[i]);

  } 

  void print_J_h()
  {
    std::cout << "h vector: " << std::endl;
    for(auto& site : sites) std::cout << site.hzv << " "; 
    std::cout << std::endl;
    std::cout << "J vector: " << std::endl;
    for(auto& site : sites) std::cout << site.nneighbs << std::endl;
  }

  void print_final_spin()
  {
    // std::cout << "sites[0].spin = " << sites[0].spin << std::endl;
    // std::cout << "sites[0].de = " << sites[0].de << std::endl;
    // std::cout << "sites[0].nneighbs = " << sites[0].nneighbs << std::endl;
    // std::cout << "sites[0].hzv = " << sites[0].hzv << std::endl;
    // std::cout << "sites[0].jzv = " ;
    // for(auto& it : sites[0].jzv ) std::cout << it << " ";
    // std::cout << std::endl;

    std::cout << "sites.size() = " << sites.size() << std::endl;
    std::cout << "Print the spin state ... " << std::endl;
    for(auto& site : sites) std::cout << site.spin << " ";  // print the spin state
    std::cout << std::endl;
  }

  std::size_t get_energies(std::vector<value_type>& en, const std::size_t offs) const
  {
    value_type energy = 0;
    for(const auto& site : sites){
      // value_type tmp = site.hzv;
      for(index_type k = 0; k < site.nneighbs; ++k)
	      // tmp += sites[site.neighbs[k]].spin * site.jzv[k] / 2;
        // tmp += ((sites[site.neighbs[k]].spin+1)/2) * site.jzv[k] / 2;
        energy += ((site.spin*sites[site.neighbs[k]].spin-1)/2) * site.jzv[k] / 2;

      // energy += tmp * ((site.spin+1)/2);
    }

    en[offs] = energy;
    return offs+1;
  }
 
  std::string get_info() const {return "algorithm: single-spin generic, variable degree";}

  private:

  std::vector<site_type> sites;
  std::vector<std::vector<double> > bound_array;

  std::mt19937 generator;
  std::function<double()> rng_; 

  };

#endif
