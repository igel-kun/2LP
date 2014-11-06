#ifndef GLOBAL_HPP
#define GLOBAL_HPP

#include "../util/graphs.hpp"
#include "../solv/solv_opts.hpp"
#include "../solv/branching.hpp"
#include "../reduction/prr.hpp" // for BRR7 (path infos)
#include "../util/statistics.hpp"
#include "defs.hpp"


namespace cr {
  // apply the B-bridge rule:
  // Let T' be a proper subtree of the component tree T of G whose unique incident B-bridge in G is {u,v} (with v in T').
  // Let S1 , S2' , and S3 denote optimal solutions of G[T'], G[T'-v], and G[T'+u] and
  // let S2 := S2 cup {{v, w} in E(G) | w = u}.
  // (a) If |S1 | < |S3 |, then apply S1 + {u, v} to G.
  // (b) If not (a), but |S2 | = |S3 |, then apply S2 to G.
  // (c) If neither (a) nor (b), then apply S3 to G.
  solution_t apply_Bbridge_rule(instance& I, stats_t& stats, const solv_options& solv_opts, const uint depth);

  // split rule: find bridge {u,v} such that deg(v)>1 and all w in N[v] have nldeg < 3, then split {u,v} from v
  // return whether a split was done
  bool apply_split_rule(instance& I);

  // apply the final reduction rule when nothing else is applicable anymore
  // CAUTION: final_RR in its current form is incorrect!!!
  //bool apply_final_RR(instance& I, stats_t& stat, solution_t& sol);

  // return whether the generalized PRR4 is applicable
  bool prr4_gen_applicable(const vertex_p& separator);

  bool trr3_gen(instance& I, stats_t& stat, const vertex_p& v, solution_t& sol);

  // return whether there is a path from u to v avoiding x and any generators
  bool exists_gen_free_path(const vertex_p& u, const vertex_p& v, const vertex_p& x);


  // perform Y-lookahead: for a vertex v with Y-pendant vw, check if deg(v)-1 > upper_bound and if so, delete vw
  bool Y_lookahead(instance& I, stats_t& stat, solution_t& sol, const vertex_p& v, const uint upper_bound);
  bool Y_lookahead(instance& I, stats_t& stat, const solv_options& opts, solution_t& sol, const uint upper_bound);
  bool Y_lookahead(instance& I, stats_t& stat, const solv_options& opts, solution_t& sol);
}

#endif
