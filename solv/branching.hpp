#ifndef BRANCHING_HPP
#define BRANCHING_HPP

#include "../util/statistics.hpp"
#include "solv_opts.hpp"
#include "defs.hpp"
#include <list>
#include <vector>




namespace cr{
  // prevent branches containing non-relevant A-bridges
  inline void add_branch(branch_op& bop, const edgelist& el, const mod_type& mt = Del){
    modlist_t ml;
    for(edge_ppc e = el.begin(); e != el.end(); ++e)
      if(!(*e)->is_Abridge()) ml.push_back(graph_mod_t(*e, mt));
    bop.branches.push_back(ml);
  }
  

  // run the complete branching recursively and return the number of operation it took
  solution_t run_branching_algo(instance& I, stats_t& stats, const solv_options& opts = default_opts, uint depth = 0);

  // solve a small instance (|V|<7) to save some branching
  inline solution_t solv_small_instance(instance& I){
    solution_t fes(edgelist_to_solution(get_a_FES(I.g)));
    I.k -= fes.size();
    I.g.vertices.clear();
    return fes;
  }



}


inline ostream& operator<<(ostream& os, const cr::branch_op& bop){
  return os << "[" << bop.type << ": "<< bop.branches << "] ("<<bop.bnum<<")";
}

inline ostream& operator<<(ostream& os, const cr::claw_leg& cl){
  return os << "{ "<< cl.head << ", "<<cl.E<<" }";
}

#endif
