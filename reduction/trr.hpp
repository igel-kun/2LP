#ifndef TRR_H
#define TRR_H

#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "../util/statistics.hpp"
#include "defs.hpp"

namespace cr{
  // update the TRR infos of all vertices, including the subtree NH
  solution_t update_TRR_infos(instance& I, stats_t& stats);
  // update parents trr_infos from child (e := ( parent -> child) )
  void update_trr_infos_from_child(const edge_p& e);

  // apply Tree Reduction Rule X to v assuming all TRRs have been applied to the subtree below v
  bool trr14_subtree(instance& I, stats_t& stats, const vertex_p& v);
  solution_t trr2_subtree(instance& I, stats_t& stats, const vertex_p& v);
  solution_t trr3_subtree(instance& I, stats_t& stats, const vertex_p& v);
  solution_t trr5_subtree(instance& I, stats_t& stats, const vertex_p& v);

  // perofrm the trrs on vertex v and return whether something changed. Store the deleted edges in local_operations
  solution_t perform_trrs(instance& I, stats_t& stats, const vertex_p& v);

  // Tree Reduction Rule 6 removes caterpillars
  bool trr6(instance& I);

  // applies the tree reduction rules from v upwards after v got thrown out of the cyclic core
  // stops at do_not_cross or a cyclic-core vertex, whichever comes first
  // returns the last found vertex on the cyclic core if there is such a vertex (or any vertex otherwise)
  // cautious: make sure v's subtree_NH are updated before running this!
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const vertex_p& do_not_cross);
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const edgelist&);
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const vertexset&);
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v);

  // apply all trrs to all vertices, return number of deletions done
  // TODO: this is much like graph::update_subtree_NH, try to collapse them into one "bottom-up applicator"
  solution_t apply_trrs(instance& I, stats_t& stats);
  solution_t apply_trrs(graph& g, stats_t& stats);


  // quick pendant addition, without checks or anything
  inline void add_nothing(graph& g, const vertex_p& v, const string& name = ""){};
  inline void add_leaf(graph& g, const vertex_p& v, const string& name = ""){
    vertex_p w(g.add_vertex_fast((name == "" ? v->name : name) + "~"));
    const edge_p e(g.add_edge_fast(w, v));
    // update TRR infos
    update_trr_infos_from_child(e);
  }
  inline void add_P2(graph& g, const vertex_p& v, const string& name = ""){
    const vertex_p w(g.add_vertex_fast((name == "" ? v->name : name) + "~"));
    const edge_p e(g.add_edge_fast(w, v));
    add_leaf(g, w, w->name);
    // update TRR infos
    update_trr_infos_from_child(e);
  }
  inline void add_2P2(graph& g, const vertex_p& v, const string& name = ""){
    add_P2(g, v, (name == "" ? v->name : name));
    add_P2(g, v, (name == "" ? v->name : name) + "~~~");
  }
  inline void add_Y(graph& g, const vertex_p& v, const string& name = ""){
    const vertex_p w(g.add_vertex_fast((name == "" ? v->name : name) + "~"));
    const edge_p e(g.add_edge_fast(w, v));
    add_2P2(g, w, w->name);
    // update TRR infos
    update_trr_infos_from_child(e);
  }
  // for an edge uv, delete uv and add a new Y-graph to v
  inline void Ygraphify(graph& g, const edge_p& e){
    vertex_p v(e->head);
    vertex_p u(e->get_tail());

    g.delete_edge(e);
    add_Y(g, v, u->name);
    u->name += '\'';
  }


}

#endif
