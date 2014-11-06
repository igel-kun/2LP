#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "branching.hpp"

namespace cr{

  // compute a lower bound for g, using more elaborate (but also slower) methods if requested
  uint compute_lower_bound(graph g, const bool more_elaborate = false);
  uint compute_lower_bound(graph& g, const solv_options& opts, const uint depth);
  
  // compute an upper bound for g, using more elaborate (but also slower) methods if requested
  // NOTE: we indicate failure by returning the empty solution!
  solution_t upper_bound_simple(const instance& _I);
};
