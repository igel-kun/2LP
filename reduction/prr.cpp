#include "prr.hpp"
#include "trr.hpp"
#include <sstream>
#include "global.hpp"

namespace cr{

  // add a leaf to v, keeping trr_infos up to date; returns edge to the new leaf
  edge_p& copy_leaf(graph& g, const vertex_p& v, const vertex_p& leaf){
    // create the vertex
    const vertex_p new_leaf(g.add_vertex_fast(leaf->name + '\''));
    new_leaf->dfs_id = leaf->dfs_id;
    // connect it to v
    const edge_p to_v(g.add_edge_fast(new_leaf, v));

    // register the new leaf in the trr_infos (subtree NH)
    update_trr_infos_from_child(to_v);

    // and return edge to the new leaf
    return to_v->get_reversed();
  }

  // add a P2 to v keeping trr_infos up to date; return edge to the center vertex
  edge_p& copy_P2(graph& g, const vertex_p& v, const vertex_p& center){
    // create the vertex
    const vertex_p new_center(g.add_vertex_fast(center->name + '\''));
    new_center->dfs_id = center->dfs_id;
    // create a new leaf at the center
    for(edge_pp leaf = center->trr_infos.leaves.begin(); leaf != center->trr_infos.leaves.end(); ++leaf)
      copy_leaf(g, new_center, (*leaf)->head);
    // connect the new center to v
    const edge_p to_v(g.add_edge_fast(new_center, v));
    // update the trr_infos of v
    update_trr_infos_from_child(to_v);
    // and return edge to the new center
    return to_v->get_reversed();
  }

  // add a Y graph to v keeping trr_infos up to date; return edge to the center vertex
  edge_p& copy_Y(graph& g, const vertex_p& v, const vertex_p& center){
    const vertex_p new_center(g.add_vertex_fast(center->name + '\''));
    new_center->dfs_id = center->dfs_id;
    // copy a leaf if there is one at center
    for(edge_pp leaf = center->trr_infos.leaves.begin(); leaf != center->trr_infos.leaves.end(); ++leaf)
      copy_leaf(g, new_center, (*leaf)->head);
    // copy all P2's that are at center
    for(edge_pp ptwo = center->trr_infos.ptwos.begin(); ptwo != center->trr_infos.ptwos.end(); ++ptwo)
      copy_P2(g, new_center, (*ptwo)->head);
    // connect the new center to v
    const edge_p to_v(g.add_edge_fast(new_center, v));
    // update the trr_infos of v
    update_trr_infos_from_child(to_v);
    // and return edge to the new center
    return to_v->get_reversed();
  }

  // copy v's complete pendant tree to vprime (including dfs_id's)
  void copy_pendant(graph& g, const vertex_p& v, const vertex_p& vprime){
    // step 1: copy leaves
    for(edge_p e : v->trr_infos.leaves) copy_leaf(g, vprime, e->head);
    
    // step 2: copy P2s
    for(edge_p e : v->trr_infos.ptwos) copy_P2(g, vprime, e->head);

    // step 3: copy Ys
    for(edge_p e : v->trr_infos.ygraphs) copy_Y(g, vprime, e->head);
  }




// ***********************************
// * NEW functions for applying PRRs *
// **********************************f

  // PRR1 applies to v if, for (x,v,w,u), any of the following holds:
  // 1. w is on the backbone
  // 2. w has deg-2 and u is inner vertex without P2
  bool prr1_applicable(const vertex_p& v){
    // v is v_i;  consider each of v's cyclic core neighbors as v_{i+1}
    const edgelist nh(v->get_cyclic_core_neighbors());
    for(auto &e : nh){
      const vertex_p& w(e->head);
      if(!w->is_on_backbone()){
        // if w's pendant is Y, then PRR1 doesn't apply
        if(w->degree() == 2){
          // if w is inner and u is inner without P2, then PRR1 applies
          const vertex_p& u(w->first_cyclic_core_neighbor_except(v)->head);
          if(u->cyc_core_degree() == 2)
             if(u->trr_infos.ptwos.empty()) return true;
        }
      } else return true; // if w is on the backbone, then apply prr1
    }
    return false;
  }

  // PRR2 applies to v if, for (u,v,w,x) v & w are cyc-core-deg-2 and any of the following hold:
  // 1. w has a Y-graph
  // 2. w is inner and u is on the backbone
  // 3. w and u are inner vertices without Y-pendants
  bool prr2_applicable(const vertex_p& v){
     // v is v_i;  consider each of v's cyclic neighbors as v_{i+1}
    uint count_single_neighbors = 0;
    const edgelist nh(v->get_cyclic_core_neighbors());
    // v has to be cyclic_core_degree-2
    if(nh.size() != 2) return false;
    for(auto &e : nh){
      const vertex_p& w(e->head);
      // 1. T^w is Y
      if(w->cyc_core_degree() == 2){
        if(!w->pendant_is_Y()){
          DEBUG5(if(!w->pendant_is_single()) FAIL("got a backbone vertex next to a Y-graph in PRR2, this shouldn't happen"));
          ++count_single_neighbors;
          // 3. the vertex x after w is on the backbone
          const vertex_p& x(w->first_cyclic_core_neighbor_except(v)->head);
          if(x->cyc_core_degree() == 2)
            if(x->is_on_backbone()) return true;
        } else return true;
      }
    }
    // 2. if w and u are inner vertices with singletons, then PRR2 applies
    if(count_single_neighbors == 2) return true;
    return false;
  }

  bool prr3_applicable(const path_info_t& info){
    // if there are no tokens, pp3 should be applicable
    if(!info.generators.empty()) return false;
    // if we have more than one separator, we can contract everything between them
    if(!info.separators.empty()) return true;
    // if we have no separators (and no Y-graphs) but the path length is 3, it's cool too :)
    if(!info.pendantYs.empty()) return false;
    if(info.length == 3) return true;
    return false;
  }

  bool prr4_applicable(const path_info_t& info){
    // prr4 is applicable if there is a token separator and a generator on the path
    // use the generalized PRR4 instead
    if(!info.separators.empty()) return prr4_gen_applicable(*info.separators.begin());
    return false;
  }

  bool prr5_applicable(const path_info_t& info){
    // prr5 is applicable if there are more than 2 generators and no separator
    return (info.separators.empty() && info.generators.size() > 2);
  }

  bool prr6_applicable(const path_info_t& info){
    const vertex_p& v(info.start->get_tail());
    // prr6 is applicable if...
    // ... the path is a cycle AND
    if(info.end->head != v) return false;
    // it has at most one generator
    if(info.generators.size() > 1) return false;
    return true;
  }

  // NOTE: prr7 assumes that prr6 was applied
  bool prr7_applicable(const path_info_t& info){
    const vertex_p& v(info.start->get_tail());
    // prr 7 s applicable if...
    // ... the path is a cycle AND
    if(info.end->head != v) return false;
    // v's pendant is a P2
    if(v->trr_infos.ptwos.empty()) return false;
    return true;
  }

  // do the edge deletion and handle the path_infos
  solution_t perform_prr1(instance& I, const vertex_p& v, path_info_t& info){
    DEBUG2(cout << "this is PRR1 for "<<*v<<endl);
    solution_t sol;
    
    const edge_p e(v->trr_infos.ygraphs.front());
    const vertex_p w(e->head);

    // remove the Ygraph from v's trr_infos
    v->trr_infos.ygraphs.pop_front();

    // kill the edge to the center of v's Y-graph
    I.delete_edge(e, sol);

    // and delete the split off caterpiller
    I.g.delete_component(w);

    return sol;
  }

  solution_t perform_prr2(instance& I, stats_t& stat, const vertex_p& v, path_info_t& info){
    DEBUG2(cout << "this is PRR2 for "<<*v<<endl);
    solution_t sol;


    // get the two neighbors of v
    const edgelist nh(v->get_cyclic_core_neighbors());
    vertex_p u(nh.front()->head);
    vertex_p w(nh.back()->head);
    
    DEBUG2(cout << "got cyclic neighbors "<<u<<" and "<<w<<endl);

    // if the path is actually a triangle,
    if(find_edge(u, w) != u->adj_list.end()){
      // save info.start in case we delete it
      const vertex_p path_start(info.start->get_tail());
      // then delete both edges incident to v
      I.delete_edges(nh, sol);
      
      // let the TRRs gnaw on the newfound edge
      // the gnawing should start at the lower one (where degree() - subtree_NH == 1)
      // also, we don't want to gnaw if it's the anchor of the path
      if(u->degree() - u->subtree_NH() == 1){
        if(u != path_start)
          sol += apply_trrs_upwards_after_cut(I, stat, u, w);
      } else {
        if(w != path_start)
          sol += apply_trrs_upwards_after_cut(I, stat, w, u);
      }
    } else {
      // else skip over v by {u,w} and update u and w (could have become separators)
      sol += v->name + "->?";
      I.k--;

      // remember if we need to update start/end of the path
      const bool update_start(info.start->head == v);
      const bool update_end(info.end->get_tail() == v);
      // order u and w correctly
      if(update_start)
        if(w == info.start->get_tail()) swap(u, w);
      if(update_end)
        if(u == info.end->head) swap(u, w);

      // delete the two incident edges and
      I.g.delete_edges(nh);
      
      // add the skip edge
      const edge_p skip_edge(I.g.add_edge_fast(u,w));

      if(update_start) info.start = skip_edge;
      if(update_end) info.end = skip_edge;

      // obstrusify u & w if they are on the deg-2 path, so verification is possible
      // because if (3) matched, then the solution actually has to be transformed,
      // although the size is the same
      if(u->cyc_core_degree() == 2) u->name += '*';
      if(w->cyc_core_degree() == 2) w->name += '*';

      // set the separator if any of the two has become one
      if(u->is_separator()) info.separators.insert(u); else info.separators.erase(u);
      if(w->is_separator()) info.separators.insert(w); else info.separators.erase(w);

      // also, the length of the path just decreased
      info.length--;

    }
    // either way, mark info invalid
    info.valid = false;

    // and delete the component of v
    I.g.delete_component(v);
    return sol;
  }

  bool perform_prr3(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "this is PRR3 for "<<info<<endl);

    if(info.separators.empty()){
      assert(info.length == 3);
      const vertex_p& u(info.start->get_tail());
      const vertex_p& v(info.end->head);
      if(!u->is_on_backbone()) { add_leaf(I.g, u); sol += perform_trrs(I, stat, u); }
      if(!v->is_on_backbone()) { add_leaf(I.g, v); sol += perform_trrs(I, stat, v); }
    } else{
      if(info.start->get_tail() == info.end->head){
        const vertex_p v(info.end->head);
        // if we're on a cycle and there are only separators, delete info.start
        vertex_p x(info.start->head);
        I.delete_edge(info.start, sol);
        // apply TRR
        sol += apply_trrs_upwards_after_cut(I, stat, x, v);
        info.valid = false;
        return true;
      } else {
        edge_p e(info.start);
        edge_p f(info.end->get_reversed());
        bool change = false;

        const vertex_p u(e->head);
        const vertex_p v(f->head);
        const vertex_p x(info.start->get_tail());
        const vertex_p y(info.end->head);
        // here, we have   x --e--> u ... v <--f-- y
        
        if(u != v){
          // u != v since length > 4, so delete everything in between and insert the edge
          e = get_next_on_deg2path(e);
          // we have ...  u --e--> ... v <--f-- y
          // so, delete e and f
          I.g.delete_edge(e);
          I.g.delete_edge(f);
          // delete the component of v
          I.g.delete_component(v);
          // and connect u and y
          const edge_p f2 = I.g.add_edge_fast(u, y);

          // finally add a leaf to our new separator u
          if(u->trr_infos.leaves.empty()) add_leaf(I.g, u);

          // note that the path stays valid if we mofify its length and the separators
          info.separators.clear();
          info.separators.insert(u);
          info.length = 2;
          info.end = f2;
          change = true;
        }
        // finally, give both endpoints leaves if they don't already have one (or a P2)
        if(!x->is_on_backbone()) { add_leaf(I.g, x); sol += perform_trrs(I, stat, x); }
        if(!y->is_on_backbone()) { add_leaf(I.g, y); sol += perform_trrs(I, stat, y); }
        return change;
      }
    }
  }

  solution_t perform_prr4(instance& I, stats_t& stat, path_info_t& info){
    vertex_p v(*info.separators.begin());
    DEBUG2(cout << "this is PRR4 for "<<v << " whose P2 are "<<v->trr_infos.ptwos<<endl);
    solution_t sol;

    const vertex_p do_not_cross(info.start->get_tail());
    const edge_p& e(v->first_cyclic_core_neighbor());
    const vertex_p w(e->head);

    I.g.delete_edge(e);
    // create a copy of v
    vertex_p vprime(I.g.add_vertex_fast((string)*v + "'"));
    // and hang it onto w
    I.g.add_edge_fast(vprime, w);
    // then copy v's pendant to vprime
    copy_pendant(I.g, v, vprime);
    // and, finally, gnaw the 2path away
    DEBUG2(cout << "first, gnaw away from "<< v << " to " << do_not_cross << endl);
    sol += apply_trrs_upwards_after_cut(I, stat, v, do_not_cross);
    DEBUG2(cout << "second, gnaw away from "<< vprime << " to " << do_not_cross<< endl);
    sol += apply_trrs_upwards_after_cut(I, stat, vprime, do_not_cross);

    // oh, right, don't forget to invalidate the path info
    info.valid = false;
    return sol;
  }

  solution_t perform_prr5(instance& I, stats_t& stat, path_info_t& info){
    DEBUG2(cout << "this is PRR5 for "<<info << endl);
    solution_t sol;
    // if the number of generators is odd, then start erasing BEFORE the first generator, otherwise, BEHIND it
    const edge_p first_to_del(((info.generators.size() & 1) == 1) ?
          info.generators.front()  :   get_next_on_deg2path(info.generators.front())  );
    // the last edge to delete is always the one pointing to the last generator
    const edge_p last_to_del(info.generators.back());
    
    // decrease k and register edges in solution
    // if generators is even, decrease k by |generators|/2 - 1, otherwise by |generators|/2 - 1/2
    // which is the same as always decreasing by floor( |generators|/2 - 1/2 )
    const uint deletes((info.generators.size() - 1) >> 1);
    I.k -= deletes;
    // we don't yet know which edges are going to be deleted..
    for(uint i = 0; i < deletes; ++i)
      sol += "[some edge between " + (string)(*info.generators.front()->head) + " and " + (string)(*info.generators.back()->head) + "]";

    // do the graph modification:
    // 0. save the first and last vertex in order to be able to reconnect
    const vertex_p first_vertex(first_to_del->get_tail());
    vertex_p last_vertex(last_to_del->head);
    // and save one vertex to be able to delete the component in between
    const vertex_p comp_to_del(first_to_del->head);
    // (save the name of first_to_del in case we cut a cycle)
    const string first_to_del_name(*first_to_del);

    // 1. delete the first and last edge,
    I.g.delete_edge(first_to_del);
    I.g.delete_edge(last_to_del);
    // 2. delete the connected component and
    I.g.delete_component(comp_to_del);
    // 3. add a new edge, unless it is already there, in which case we broke a cycle
    edge_p e(find_edge(first_vertex, last_vertex));
    if(e != first_vertex->adj_list.end()){
      // if we broke a cycle, then delete the edge (first_vertex, last_vertex) (destroy the 2-cycle)
      sol += first_to_del_name;
      I.k--;
      // update the TRRs just one edge, and let the parent vertex update the rest
      sol += apply_trrs_upwards_after_cut(I, stat, last_vertex, first_vertex);
    } else I.g.add_edge_fast(first_vertex, last_vertex);
    // mark the info invalid, since no other rule appies to this path any more
    info.valid = false;

    return sol;
  }

  solution_t perform_prr6(instance& I, stats_t& stat, path_info_t& info){
    DEBUG2(cout << "this is PRR6 for "<<info << endl);
    solution_t sol;
    
    const vertex_p v(info.end->head);
    if(!info.separators.empty()) FAIL("encountered cycle with separators in PRR6! This should have been reduced by PRR3!");
    if(info.generators.empty()){
      // if there are neither separators nor generators, then it's a triangle of singletons and at most 1 Ygraph
      // so delete the second edge
      edge_p e(get_next_on_deg2path(info.start));
      vertex_p y(e->head);
      I.delete_edge(e, sol);
      // apply TRRs
      vertex_p x(info.start->head);
      sol += apply_trrs_upwards_after_cut(I, stat, x, v);
      sol += apply_trrs_upwards_after_cut(I, stat, y, v);
    } else {
      // if there is exactly one generator, then delete an edge that's far from the generator
      edge_p e(info.generators.front()->get_reversed());
      // if the generator is not next to v go a step further
      if(e->head != v) e = get_next_on_deg2path(e);
      // then, delete e
      vertex_p x(e->get_tail());
      vertex_p y(e->head);
      I.delete_edge(e, sol);
      // apply TRRs
      sol += apply_trrs_upwards_after_cut(I, stat, x, v);
      sol += apply_trrs_upwards_after_cut(I, stat, y, v);
    }
    info.valid = false;
    return sol;
  }

  solution_t perform_prr7(instance& I, stats_t& stat, path_info_t& info){
    DEBUG2(cout << "this is PRR7 for " << info << endl);
    solution_t sol;
    const vertex_p v(info.end->head);

    const edge_p e =
      info.start->head->is_generator() ? info.start : get_next_on_deg2path(info.start);
    // if v has a P2, then delete an edge towards one of them and let the TRRs do the rest
    vertex_p x(e->head);
    vertex_p y(e->get_tail());
    I.delete_edge(e, sol);
    // apply TRRs
    sol += apply_trrs_upwards_after_cut(I, stat, x, v);
    sol += apply_trrs_upwards_after_cut(I, stat, y, v);
    // mark this path invalid
    info.valid = false;

    return sol;
  }


// ********************************************
// * NEW overhead functions for applying PRRs *
// ********************************************

#define BUDGET_EXCEEDED(x) {I.k = -1; info.valid = false; return x;}
  bool prr12_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    // find a vertex with a Y-pendant and apply PRR1&2
    bool change = false;
    while(info.valid && !info.pendantYs.empty()){
      const vertex_p v(info.pendantYs.front());
      DEBUG2(cout << "checking prr12 for "<<*v<<endl);

      // in both PRR1 and PRR2, we need at least one edge deletion, so we'll return failure if the budget ran out
      if(prr1_applicable(v)){
        change = true;
        if(I.k <= 0) BUDGET_EXCEEDED(false) else{
          DO_STAT(stat.reduct_application[PRR1]++);
          // remove v from the pendantYs (and the separator, if necessary, and start/end)
//          info.separators.erase(v);
          info.pendantYs.pop_front();
          sol += perform_prr1(I, v, info);
        }
      } else if(prr2_applicable(v)){
        change = true;
        if(I.k <= 0) BUDGET_EXCEEDED(false) else {
          DO_STAT(stat.reduct_application[PRR2]++);
//          info.separators.erase(v);
          info.pendantYs.pop_front();
          sol += perform_prr2(I, stat, v, info);
        }
        // if neither prr1 nor prr2 applies, we face a short path with just a Y-graph -> nothing to do
      } else break;
    }
    return change;
  }
  
  // if there is no token generator on the degree-2-path, then replace the path by a simple P4 (given its not already shorter)
  bool prr3_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "checking prr3"<<endl);
    if(prr3_applicable(info)){
      DO_STAT(stat.reduct_application[PRR3]++);
      return perform_prr3(I, stat, info, sol);
    }
    else return false;
  }

  // if there is a token separator v on a path with tokens, split the path open at v
  bool prr4_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "checking prr4"<<endl);
    if(prr4_applicable(info))
      if(I.k <= 0) BUDGET_EXCEEDED(false) else{
        DO_STAT(stat.reduct_application[PRR4]++);
        sol += perform_prr4(I, stat, info);
        return true;
      }
    else return false;
  }

  // if we have a chained degree-2 path of even (odd) length, compress it to 2 (to 3) generators
  bool prr5_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "checking prr5"<<endl);
    if(prr5_applicable(info)){
      if(I.k <= 0) BUDGET_EXCEEDED(false) else{
        DO_STAT(stat.reduct_application[PRR5]++);
        sol += perform_prr5(I, stat, info);
        return true;
      }
    }
    return false;
  }

  bool prr6_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "checking prr6"<<endl);
    if(prr6_applicable(info)){
      if(I.k <= 0) BUDGET_EXCEEDED(false) else{
        DO_STAT(stat.reduct_application[PRR6]++);
        sol += perform_prr6(I, stat, info);
        return true;
      }
    }
    return false;
  }

  bool prr7_from_infos(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    DEBUG2(cout << "checking prr7"<<endl);
    if(prr7_applicable(info)){
      if(I.k <= 0) BUDGET_EXCEEDED(false) else{
        DO_STAT(stat.reduct_application[PRR7]++);
        sol += perform_prr7(I, stat, info);
        return true;
      }
    }
    return false;
  }


  // delete the second edge (unless delete_first_edge is set) on the path info_delete and update paths and infos
  bool PRR8_delete_second_edge(instance& I, stats_t& stat, solution_t& sol,
      const list<path_info_t>::iterator& info_delete,
      const list<path_info_t>::iterator& info_remain,
      list<path_info_t>& infos,
      unordered_map<vertex_p, list<path_info_t>::iterator, vertex_hasher>& paths,
      bool delete_first_edge = false){ 

    DEBUG2(cout << "This is PRR8 for the paths (gonna keep the second one):"<<endl<<*info_delete<<endl<<*info_remain<<endl);
    // get the second edge, or the unique edge if info_delete has length 1
    edge_p e = info_delete->start;
    if(!delete_first_edge && (info_delete->length > 1)) e = get_next_on_deg2path(e);

    const vertex_p v(info_delete->start->get_tail());
    const vertex_p w(info_delete->end->head);
    vertex_p x(e->head);
    vertex_p y(e->get_tail());


    DO_STAT(stat.reduct_application[PRR8]++);
    I.delete_edge(e, sol);

    // update the path associated with the end vertex of both info_delete and info_remain
    paths[info_remain->end->head] = info_remain;
    // update the list of paths
    infos.erase(info_delete);

    sol += apply_trrs_upwards_after_cut(I, stat, x, v);
    sol += apply_trrs_upwards_after_cut(I, stat, y, v);
    DEBUG2(cout << "done applying PRR8 to this pair"<<endl);
    return true;
  }
  
  // apply PRR8 to a new path in a list of paths
  // anything - separator  ==>  delete the second edge on the non-separator path
  // singletons or Y - Y  ==>  delete at Y
  //
  // NOTE: Prr8 can destroy the complete graph, so after any application has been found, return to caller!!!
  bool apply_PRR8(instance& I, stats_t& stat, solution_t& sol,
      const list<path_info_t>::iterator& path,
      list<path_info_t>& infos,
      unordered_map<vertex_p, list<path_info_t>::iterator, vertex_hasher>& paths){

    assert(!infos.empty());
    const vertex_p& v(path->end->head);

    // path is only interesting if it has no generator
    if(path->generators.empty()){
      // from here, neither path nor old_path has generators
      DEBUG2(cout << "PRR8: checking "<<*path<<endl);

      // if another deg2-path to v exists, then apply PRR8
      if(paths.find(v) != paths.end()){
        auto old_path = paths[v];
        if(!old_path->pendantYs.empty())
          return PRR8_delete_second_edge(I, stat, sol, old_path, path, infos, paths);
        if(!path->pendantYs.empty())
          return PRR8_delete_second_edge(I, stat, sol, path, old_path, infos, paths);

        // from here, neither path nor old_path has a Ygraph
        if(path->separators.empty() || old_path->separators.empty()){
          // if one of the two paths contains separators, then delete the second edge of the other
          if(!path->separators.empty())
            return PRR8_delete_second_edge(I, stat, sol, old_path, path, infos, paths);
          if(!old_path->separators.empty())
            return PRR8_delete_second_edge(I, stat, sol, path, old_path, infos, paths);

          // from here, neither path nor old_path contains separators (both are singlton-paths)
          // if either endpoint is on the backbone (has leaf or P2), delete an edge away from that endpoint
          if(path->start->get_tail()->is_on_backbone())
            return PRR8_delete_second_edge(I, stat, sol, path, old_path, infos, paths);
          if(path->end->head->is_on_backbone())
            return PRR8_delete_second_edge(I, stat, sol, path, old_path, infos, paths, (path->length < 3));
        }
      } else paths[v] = path;
    }
    return false; 
  }

  // get infos along the deg2-path initiated by e, assumes that e is on a cyc-core-deg2-path
  path_info_t get_path_infos(const edge_p& e, const uint dfs_id){
    path_info_t infos;
    const vertex_p& v(e->get_tail());
    edge_p next(e);

    // also count end-vertices of the path as generators
    if(v->is_generator()) infos.end_generators.push_back(next->get_reversed());
    // set startpoint of path
    infos.start = e;
    // while the next vertex is valid and on the deg2-path and we're not back where we started (in case we're on a cycle)
    while((next->head->cyc_core_degree() == 2) && ((infos.length == 0) ? true : (next->head != v))){
      const vertex_p& v(next->head);

      DEBUG2(cout << "exploring path along " << *next << " with trr("<<*v<<") = "<< v->trr_infos << endl);
      DEBUG2(cout << "deg("<<*v<<")="<<v->degree()<<", subtree("<<*v<<")="<<v->subtree_NH()<< endl);

      // mark this vertex seen
      v->dfs_id = dfs_id;
      
      // count the length
      ++infos.length;

      // analyze the next vertex
      if(v->is_separator()) infos.separators.insert(v);
      if(v->is_generator()) infos.generators.push_back(next);
      
      // mark vertices with Y-graph pendants
      if(!v->trr_infos.ygraphs.empty()) infos.pendantYs.push_back(v);

      // advance along the path
      next = get_next_on_deg2path(next);
      DEBUG2(cout << "next, going along "<< next << " with cyc core deg "<< next->head->cyc_core_degree() << endl);
    }
    // also count end-vertices of the path as generators
    if(next->head->is_generator() && (next->head != v)) infos.end_generators.push_back(next);

    DEBUG2(cout << "finished exploring at "<< *next << endl);
    // set the end of the path and mark it valid
    infos.end = next;
    infos.valid = true;

    return infos;
  }

  // act on path information
  // CAUTION: this can destroy alot - use safe iterators!
  bool act_on_path_info(instance& I, stats_t& stat, path_info_t& info, solution_t& sol){
    // first look for Y-graphs and apply PRR1 & PRR2
    DEBUG2(cout << "acting on " << info << endl);
    bool change = false;

    change |= prr12_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR12? "<<info.valid<<endl);
    if(!info.valid) return change;

    // changes made by PRR3 should not cause a re-run
    change |= prr3_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR3? "<<info.valid<<endl);
    if(!info.valid) return change;

    change |= prr4_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR4? "<<info.valid<<endl);
    if(!info.valid) return change;

    change |= prr5_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR5? "<<info.valid<<endl);
    if(!info.valid) return change;

    change |= prr6_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR6? "<<info.valid<<endl);
    if(!info.valid) return change;

    change |= prr7_from_infos(I, stat, info, sol);
    DEBUG2(cout << "G changed? "<<change<<", is our path still valid after PRR7? "<<info.valid<<endl);

    return change;
  }

  // find the first path of v
  edge_p find_first_path(const vertex_p v, const uint dfs_id){
    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
      if(e->head->is_on_cyclic_core() && (e->head->dfs_id != dfs_id))
        return e;
    return v->adj_list.end();
  }
  // find the next path, ignoring info.end and any paths we have seen or are not cyc-core-deg2-paths
  edge_p find_next_path(const path_info_t info, const uint dfs_id){
    // get the current edge and vertex
    edge_p result(info.start);
    const vertex_p& v(result->get_tail());
    // advance result until it's either at the end of v's neighbor list or it's valid
    do{
      ++result;
      if(result == v->adj_list.end()) break;
    } while((result->head->dfs_id == dfs_id) || (!result->head->is_on_cyclic_core()) || (result == info.end->get_reversed()));
    DEBUG2(if(result != v->adj_list.end()) cout << "next is going to be "<< result << endl; else cout << "this is the end"<<endl);
    return result;
  }

  // for each path, get its info and act on it
  bool apply_prrs_to_vertex(instance& I, stats_t& stat, solution_t& sol,
      list<path_info_t>& infos,
      const vertex_p& v,
      const uint dfs_id){ 
    // first, apply trr3_gen
    if(trr3_gen(I, stat, v, sol)) return true;

    // for PRR8, save a map mapping end-vertices to paths
    unordered_map<vertex_p, list<path_info_t>::iterator, vertex_hasher> paths;
   
    // mark v as considered now
    v->dfs_id = dfs_id;

    bool change = false;
    edge_p path_via = find_first_path(v, dfs_id);

    while(path_via != v->adj_list.end()){
      // collect the path infos
      path_info_t path_info(get_path_infos(path_via, dfs_id));
      DEBUG2(cout << "found path:"<< path_info << endl);
      // increase the iterator
      path_via = find_next_path(path_info, dfs_id);
      // finally, act on the path
      change |= act_on_path_info(I, stat, path_info, sol);
      // if it's still valid, safe the path for later use by BRR78
      if(path_info.valid) {
        infos.push_front(path_info);
        // Note that PRR8 can destroy the next path, so start over
        if(apply_PRR8(I, stat, sol, infos.begin(), infos, paths)) return true;
      }
      DEBUG2(cout << "done acting on this path, solution now " << sol << ", budget left: " << I.k << ", changed? "<<change<<endl);
      // break if v's cyclic core degree has become 2 (or if it was 2 before and we were on a cycle (unique path))
      if(v->cyc_core_degree() < 3) break;
    }
    return change;
  }

  // apply PRRs and TRRs (if v got kicked off the cyclic core) to v, updating v if it is lost
  // return whether the graph changed
  bool apply_prrs_and_trrs_to_vertex(instance& I, stats_t& stat, solution_t& sol,
      list<path_info_t>& infos,
      vertex_p& v,
      const uint dfs_id){
    if(apply_prrs_to_vertex(I, stat, sol, infos, v, dfs_id)){
      DEBUG2(cout << "all paths of "<< *v << " are done, now applying TRRs from it"<<endl);
      // in case all our paths have been cut, apply the TRRs upwards
      sol += apply_trrs_upwards_after_cut(I, stat, v);
      return true;
    } else return false;
  }

  // apply the path reduction rules to the instance I
  solution_t apply_prrs(instance& I, const solv_options& opts, stats_t& stat, list<path_info_t>& infos){
    // make sure the cyclic core is up to date
    solution_t sol(update_TRR_infos(I, stat));
    bool change;

    // exhaustively apply PRRs
    do{
      // clear infos from old iterations, since their deg2-paths are likely destroyed now
      infos.clear();

      // save the last vertex with cyc_core_degree == 2, in case we have only a cycle
      vertex_p last_cyc;
      bool has_cyc_deg3 = false;
      bool has_cyc_deg2 = false;
      change = false;
  
      // get a new dfs_id from I.g
      uint dfs_id(I.g.get_dfs_id());

      for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end();){
        DEBUG2(cout << "considering "<< v->name << " with trr: " << v->trr_infos << " & subs: "<< v->subtree_NH() << endl);
        // also make sure that it's on a cycle
        if(v->cyc_core_degree() > 2){
          DEBUG2(cout << *v << " is cyclic-core-degree >2"<< endl);
          has_cyc_deg3 = true;
          list<path_info_t> v_infos;

          // run PRR1 on v (might kill v's Y-graph)
          if(prr1_applicable(v)) perform_prr1(v);
          // run path detection and PRR application from v
          if(!apply_prrs_and_trrs_to_vertex(I, stat, sol, v_infos, v, dfs_id)) {
            // if none of the PRRs apply to v, then try the Y-lookahead rule and go on (Y-lookahead doesn't change paths)
            if(I.g.vertices.size() < opts.max_size_for_Y_lookahead) Y_lookahead(I, stat, sol, v, I.k);
            infos += v_infos;
            ++v;
          } else change = true;
        } else {
          // why is there no intrinsic invalid value of iterators?
          if(v->cyc_core_degree() == 2) has_cyc_deg2 = true;
          ++v;
        }
      }
      // handle degerate inputs:
      if(!change && !has_cyc_deg3){
        // if we are on a cycle, then apply PPRs from last_cyc
        if(has_cyc_deg2){
          DEBUG2(cout << "okay, we're on a cycle! Applying to all cyc-core vertices..."<<endl);
          for(vertex_p v = I.g.vertices.begin();  v != I.g.vertices.end();){
            // get a new DFS ID 'cause we want to check out the cycle from multiple perspectives
            dfs_id = I.g.get_dfs_id();

            DEBUG2(cout << "considering "<< v <<endl);
            if(v->is_on_cyclic_core()){
              if(apply_prrs_and_trrs_to_vertex(I, stat, sol, infos, v, dfs_id)) {
                change = true;
                break;
              } else ++v;
            } else ++v;
          }
        }
        // finally, apply TRR6
        change |= trr6(I);
      }
      DEBUG2(cout << "PRRs applied, solution so far: "<< sol << " did we change G? "<<change<<endl);
      // clean up and return the solution
    } while(change && (I.k > 0) && (!I.g.vertices.empty()));
    // we may have to apply TRR6 if we finished without creating a single cycle first
    if(I.k == 0) trr6(I);
    
    DEBUG2(cout << "returning from apply_prrs"<<endl);
    return sol;
  }
  
} // end of namespace




