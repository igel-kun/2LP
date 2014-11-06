#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include "defs.hpp"
#include "../solv/defs.hpp"
#include "../reduction/defs.hpp"
#include "b_vector.hpp"
#include "../solv/defs.hpp"

using namespace std;

#define TAKE_STATS

#ifdef TAKE_STATS
#define DO_STAT(x) x;
#else
#define DO_STAT(x) {}
#endif

namespace cr {
  struct stats_t;


  // if branching number x occurs a times and branching number y occurs b times,
  // then, in average, this is equal to branching number (ax+by)/(a+b) occuring a+b times
  inline pair<uint, float> combine(const pair<uint, float> p1, const pair<uint, float> p2){
    const uint num = p1.first + p2.first;
    const float avg = (float) (p1.second * p1.first + p2.second * p2.first ) / num;
    return make_pair(num, avg);
  }

 
  struct stats_t{
    uint input_vertices;
    uint input_edges;
    uint input_FES;

    uint searchtree_nodes;
    uint searchtree_depth;
    // the number of applications for each reduction rule
    unordered_map<reduction_type, uint, hash<int> > reduct_application;
    // the number of applications and average branching number for each branching rule
    unordered_map<branch_type, pair<uint, float>, hash<int> > bnum_avg;

    stats_t():input_vertices(0), input_edges(0), input_FES(0),searchtree_nodes(0), searchtree_depth(0){}
    
    stats_t(graph& g):searchtree_nodes(0), searchtree_depth(0){
      input_vertices = g.vertices.size();
      input_edges = g.num_edges();
      input_FES = get_FES(g);
    }

    // add the application of a branching rule of type t and branching vector b_vec to the statistic
    void add_BRule(const branch_type& t, const uint* b_vec, const uint entries){
      pair<uint, float>& entry(bnum_avg[t]);
      entry = combine(entry, make_pair(1U, branch_number(b_vec, entries)));
    }
    void add_BRule(const branch_op& bo){
      pair<uint, float>& entry(bnum_avg[bo.type]);
      entry = combine(entry, make_pair(1U, branch_number(bo)));
    }
  
    // compute the overall average branching number
    float get_avg_bnum() const {
      pair<uint, double> accu(make_pair(0U, 0.0));
      for(auto i : bnum_avg) accu = combine(accu, i.second);
      return accu.second;
    }
  };

}

ostream& operator<<(ostream& os, const cr::stats_t& stat);

void output_parser_friendly(ostream& os, cr::stats_t& stat);

#endif
