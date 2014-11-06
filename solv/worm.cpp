#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "../reduction/trr.hpp"
#include "../reduction/prr.hpp"
#include "../reduction/global.hpp"
#include "bounds.hpp"
#include "branching.hpp"
#include "worm.hpp"

#include <unordered_map>
#include <unordered_set>


namespace cr{
  
  inline bool has_favorable_pendant(const vertex_pc& v){
    if(!v->trr_infos.leaves.empty()) return true;
    if(!v->trr_infos.ptwos.empty()) return true;
    return false;
  }

  inline bool has_long_deg2_path(const vertex_p& v){
    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
      if(e->head->cyc_core_degree() == 2)
        if(get_next_on_deg2path(e)->head->cyc_core_degree() == 2) return true;
    return false;
  }


  vertex_p find_backbone_vertex(graph& g){
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v){
      switch(v->cyc_core_degree()){
        case 0:
        case 1: continue;
        case 2: 
          if(has_favorable_pendant(v)) return v;
          break;
        default:
          if(has_favorable_pendant(v)) return v;
          if(has_long_deg2_path(v)) return v;
      }
    }
    return g.vertices.end();
  }


  // grow a caterpillar to the next cyc-core-degree->2 vertex, avoiding do_not_cross
  bool grow_to_deg3(edge_p& e, const vertex_p& do_not_cross, instance& I, stats_t& stats, const solv_options& opts, solution_t& sol){

    while(true){
      if(e->head == do_not_cross) return false;
      if(I.k <= 0) return false;

      vertex_p v(e->head);
      vertex_p u(e->get_tail());

      // if v has any Y-graph pendants, delete the connecting edge
      edgelist& ygraphs(v->trr_infos.ygraphs);
      while(!ygraphs.empty()){
        I.delete_edge(ygraphs.front(), sol);
        ygraphs.pop_front();
      }

      // if v has a P2, then all edges except e are Y-graphed and the caterpillar ends here
      if(!v->trr_infos.ptwos.empty()){
        for(edge_p f = v->adj_list.begin(); f != v->adj_list.end(); ++f)
          if(f->head != u) Ygraphify(I.g, f);
        return false;
      }
      
      // if v is on a deg2-path, then continue on the path
      if(v->cyc_core_degree() == 2) e = get_next_on_deg2path(e); else return true;  // otherwise, we've found a deg3 vertex
    }
  }



  bool grow_cat_from(edge_p& e, const vertex_p& do_not_cross, instance& I, stats_t& stats, const solv_options& opts, solution_t& sol){
    vertex_p v(e->head);
    vertex_p u(e->get_tail());
    // protect v
    v->prot = true;
    // 1. apply TRRs
    sol += update_TRR_infos(I, stats);
    // 2. find all deg-2-paths from v and apply PRRs
    list<path_info_t> infos;
    apply_prrs_to_vertex(I, stats, sol, infos, v, I.g.get_dfs_id());

    // quick sanity check: if I have less than 7 vertices, then I cannot have a 2-claw, thus the solution is FES
    if(I.g.vertices.size() < 7) { sol += solv_small_instance(I); return (I.k >= 0); }
    if(I.k <= 0) return false;

    // 3. branch on which remaining path to take
    edgelist choices;
    for(edge_p f = v->adj_list.begin(); f != v->adj_list.end(); ++f)
      if((f->head != u) && (f->head->is_on_cyclic_core()))
        choices.push_back(f);

    solution_t best_sol;
    for(edge_p c : choices){
      // copy instance
      unordered_map<uint, vertex_p> id_to_vertex;
      instance Iprime(I, &id_to_vertex);
      solution_t solprime;
      // Y-graphify 
      for(edge_p f : choices) if(f->head != c->head)
        Ygraphify(Iprime.g, find_edge(id_to_vertex[f->get_tail()->id], id_to_vertex[f->head->id]));

      edge_p cprime(find_edge(id_to_vertex[c->get_tail()->id], id_to_vertex[c->head->id]));
      // if we can successfully grow to the next deg3 vertex, then call grow_cat_from again
      if(grow_to_deg3(cprime, id_to_vertex[do_not_cross->id], Iprime, stats, opts, solprime)){
        if(!grow_cat_from(cprime, id_to_vertex[do_not_cross->id], Iprime, stats, opts, solprime)) continue;
      } else // otherwise, call run_worm_trace
        if(!run_worm_trace(Iprime, solprime, stats, opts)) continue;

      // if we arrive here, we've found a solution in solprime
      best_sol = solprime;
      // now looking for smaller solutions than solprime
      I.k = solprime.size() - 1;
    }
    if(best_sol.empty()) return false; else { sol += best_sol; return true; }
  }

  bool run_worm_trace(instance& I, solution_t& sol, stats_t& stats, const solv_options& opts){
    // quick sanity check: if I have less than 7 vertices, then I cannot have a 2-claw, thus the solution is FES
    if(I.g.vertices.size() < 7) { sol += solv_small_instance(I); return (I.k >= 0); }
    if(I.k <= 0) return false;

    vertex_p v = find_backbone_vertex(I.g);
    
    if(v == I.g.vertices.end()){
    } else {
      // v is a backbone vertex in some optimal solution, start growing the caterpillar here
    }

    return false;
  }

}; // end namespace


