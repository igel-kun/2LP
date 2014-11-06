#include "statistics.hpp"
#include <cmath>

#define PRECISION 20
// get the branching number that would create a searchtree of this size and depth
// the correct formula is sum_{i=0}^d x^i = n, which equals  (x^{d+1} - 1)/(x-1) = n,
float get_bnum_from_ST(const uint size, const uint depth){
  // make sure inputs are sane
  if(!size || !depth) return 0.0;
  // we'll use bracketing to solve it
  float lower = 0.0;
  float upper = 4.0;
  // start with an approximate x = 2
  for(uint i = 0; i < PRECISION; ++i){
    const float x = (lower + upper)/2;
    if( (pow(x, depth + 1)-1)/(x-1) > size ) upper = x; else lower = x;
  }
  return (lower + upper)/2;
}

ostream& operator<<(ostream& os, const cr::stats_t& stat){
  os << "=== statistics: ==="<<endl;
  os << "fes: "<<stat.input_FES<< " ST nodes: "<<stat.searchtree_nodes<< " ST depth: "<<stat.searchtree_depth<<endl;
  os << "Reductions: "<< stat.reduct_application << endl;
  os << "Branchings: ";
  for(pair<cr::branch_type, pair<uint, float> > i : stat.bnum_avg)
    os << "(" << i.first << ": " << i.second << ") ";
  os << endl;
  os << "Overall average branching number: "<< stat.get_avg_bnum()<<endl;
  os << "branching number from ST-size vs depth: "<< get_bnum_from_ST(stat.searchtree_nodes, stat.searchtree_depth) << endl;
  os << "branching number from ST-size vs fes: "<< get_bnum_from_ST(stat.searchtree_nodes, stat.input_FES) << endl;
  return os;
}

void output_parser_friendly(ostream& os, cr::stats_t& stat){
  os << stat.input_vertices << '\t' << stat.input_edges << '\t' << stat.input_FES<< '\t' << stat.searchtree_nodes << '\t' << stat.searchtree_depth << '\t' << stat.reduct_application[cr::TRR1] << '\t'<< stat.reduct_application[cr::TRR2] << '\t'<< stat.reduct_application[cr::TRR3] << '\t'<< stat.reduct_application[cr::TRR4] << '\t'<< stat.reduct_application[cr::TRR5] << '\t'<< stat.reduct_application[cr::TRR6] << '\t'<<stat.reduct_application[cr::PRR1] << '\t'<< stat.reduct_application[cr::PRR2] << '\t'<< stat.reduct_application[cr::PRR3] << '\t'<< stat.reduct_application[cr::PRR4] << '\t'<< stat.reduct_application[cr::PRR5] << '\t'<< stat.reduct_application[cr::PRR6] << '\t'<< stat.reduct_application[cr::PRR7] << '\t' << stat.reduct_application[cr::Fin] <<  '\t' << stat.reduct_application[cr::YL] << '\t'<< stat.bnum_avg[cr::Triangle].first << '\t' << stat.bnum_avg[cr::Triangle].second << '\t'<< stat.bnum_avg[cr::Claw0].first << '\t' << stat.bnum_avg[cr::Claw0].second << '\t'<< stat.bnum_avg[cr::Claw1].first << '\t' << stat.bnum_avg[cr::Claw1].second << '\t'<< stat.bnum_avg[cr::Claw2].first << '\t' << stat.bnum_avg[cr::Claw2].second << '\t'<< stat.bnum_avg[cr::Claw3].first << '\t' << stat.bnum_avg[cr::Claw3].second << '\t'<< stat.bnum_avg[cr::Deg2Path].first << '\t' << stat.bnum_avg[cr::Deg2Path].second << '\t'<<stat.bnum_avg[cr::Token].first << '\t' << stat.bnum_avg[cr::Token].second << '\t'<< stat.bnum_avg[cr::Bbridge].first << '\t' << stat.bnum_avg[cr::Bbridge].second << '\t' << stat.get_avg_bnum() << endl;
}
