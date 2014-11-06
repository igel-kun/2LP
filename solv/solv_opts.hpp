#ifndef SOLV_OPTS_HPP
#define SOLV_OPTS_HPP

namespace cr {
  struct solv_options{
    uint fast_lower_bound_layers_wait;
    uint slow_lower_bound_layers_wait;
    bool use_Bbridge_rule;
    bool elaborate_branch_selection;
    float keep_searching_if_bnum_above;
    uint max_size_for_Y_lookahead;
  };
  const solv_options default_opts = {
    1, // fast_lower_bound_layers_wait
    8, // slow_lower_bound_layers_wait
    true, // use_Bbridge_rule
    false, // elaborate branch selection
    2.5, // keep searching for branching applications if bnum is above this number
    30, // maximum size of G to allow performing Y_lookahead
  };

};


inline ostream& operator<<(ostream& os, const cr::solv_options& opts){
  return os << " 1. slow lower bound computation each " << opts.slow_lower_bound_layers_wait << " layers"<<endl;
}

#endif

