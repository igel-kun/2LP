#ifndef WORM_HPP
#define WORM_HPP

#include "../util/statistics.hpp"
#include "solv_opts.hpp"
#include "defs.hpp"


namespace cr{
  // run the complete worm-trace and return the number of edge deletions it took
  bool run_worm_trace(instance& I, solution_t& sol, stats_t& stats, const solv_options& opts = default_opts);

}
#endif
