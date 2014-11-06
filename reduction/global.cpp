#include "global.hpp"
#include "../util/statistics.hpp"
#include "../solv/bounds.hpp"
#include "trr.hpp"
#include <algorithm>

namespace cr {


  // returns true iff b can be split at it's head
  bool is_splittable(const edge_p& b){
    const vertex_p& v(b->head);
    // if v is a leaf, then it's not splittable
    if(v->degree() <= 1) return false;
    // if v has nldeg > 2, then v is not splittable
    if(v->nldeg() > 2) return false;
    // if any of v's neighbors has nldeg > 2, then v is not splittable
    for(list<edge>::const_iterator w = v->adj_list.begin(); w != v->adj_list.end(); ++w)
      if(w->head->nldeg() > 2) return false;
    return true;
  }
  
  // apply the split rule to a specific bridge b, which is splittable
  // NOTE: this might invalidate (delete) some bridges
  void apply_split_rule(instance& I, edge_p b){
    const vertex_p& v(b->head);
    const vertex_p& u(b->get_tail());
    // delete the bridge
    // and add a new deg-1 neighbor to u
    const vertex_p vprime(I.g.add_vertex_fast(v->name + '\''));
    DEBUG2(cout << "replacing "<<b<<" by the new edge ("<<vprime<<", "<<u<<")"<<endl);
    I.g.add_edge_fast(vprime, u);
    I.g.delete_edge(b);
    // apply TRRs from vprime upwards
  }
  
  // split rule: find bridge {u,v} such that deg(v)>1 and all w in N[v] have nldeg < 3, then split {u,v} from v
  bool apply_split_rule(instance& I){
    // first, get all bridges
    edgelist bridges = I.g.get_Bbridges();
    // for each bridge, see if the rule applies
    for(edge_pp b = bridges.begin(); b != bridges.end();)
      if(is_splittable(*b)){
        // remember to advance b before deleting anything
        apply_split_rule(I, *(b++));
        return true;
      } else if(is_splittable((*b)->get_reversed())){
        // remember to advance b before deleting anything
        apply_split_rule(I, (*(b++))->get_reversed());
        return true;
      } else ++b;
    // finally, apply the TRRs if something happened
    return false;
  }



// ==== Reduction Rule 2.1 ======

  // verify conditions of Reduction Rule 1:
  // a) G[L] does not contain B-bridges
  // b) the cyclic neighborhood of u is 2 (neighbors v,w)
  // c) no vertex in L\ N[u] has nldeg >2
  bool final_RR_applicable(graph& g, const edge_p& inbridge){
    list<vertex_p> to_consider;
    DEBUG2(cout << "checking if FinRR applies to "<<inbridge<<endl);
    // check b)
    const edgelist u_nh(inbridge->head->get_cyclic_neighbors());
    if(u_nh.size() != 2) return false;
    // don't cross inbridge
    const uint dfs_id(g.get_dfs_id());
    inbridge->head->dfs_id = dfs_id;
    if(u_nh.front()->head->cyc_core_degree() > 3) return false;
    if(u_nh.back()->head->cyc_core_degree() > 3) return false;
    u_nh.front()->head->dfs_id = dfs_id;
    u_nh.back()->head->dfs_id = dfs_id;
    // and start with the two neighbors of u
    to_consider.push_back(u_nh.front()->head);
    to_consider.push_back(u_nh.back()->head);
    
    while(!to_consider.empty()){
      const vertex_p v(to_consider.front());
      to_consider.pop_front();
      for(edge_pc e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
        if((e->head->dfs_id != dfs_id) && (!e->is_Abridge())){
          DEBUG2(cout << "FinRR: checking "<<e<<endl);
          // check a)
          if(e->is_Bbridge()) return false;
          // check c)
          if(e->head->nldeg() > 2) return false;

          e->head->dfs_id = dfs_id;
          to_consider.push_back(e->head);
        }
    }
    return true;
  }
  // return the only cyc-core-deg3 vertex in a leaf of the component tree
  bool get_final_RR_application(graph& g, edge_p& inbridge){
    const edgelist B(g.get_Bbridges());
    DEBUG2(cout << "got B-bridges: "<<B<<endl);
    // if we have a cycle without B-bridges, then just return any edge pointing to a cyclic vertex
    if(B.empty() && (get_FES(g) < 3))
      for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v) if(v->is_on_cyclic_core())
        for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
          if(e->head->is_on_cycle())
            if(final_RR_applicable(g, e)) { inbridge = e; return true;}
    // otherwise, we have B-bridges
    for(edge_ppc e = B.begin(); e != B.end(); ++e){
      if(final_RR_applicable(g, *e)) { inbridge = *e; return true;}
      if(final_RR_applicable(g, (*e)->get_reversed())) { inbridge = (*e)->get_reversed(); return true;}
    }
    DEBUG2(cout << "couldn't find a comptree leaf to apply final_RR"<<endl);
    return false;
  }

  edge_p get_farthest_edge_on_cycle(const vertex_p& u){
    edgelist nh(u->get_cyclic_neighbors());
    if(nh.size() > 2) FAIL("epic fail in get_farthest_edge_on_cycle");
    edge_p left = nh.front();
    edge_p right = nh.back();
    do{
      left = get_next_on_cycle(left);
      right = get_next_on_cycle(right);
    } while((left != right->get_reversed()) && (left->head != right->head));
    return left;
  }
  inline bool is_dirty(const vertex_p& v){
    return !(v->trr_infos.ptwos.empty());
  }

  // apply the final reduction rule when nothing else is applicable anymore
  bool apply_final_RR(instance& I, stats_t& stat, solution_t& sol){
    // make sure bridges are marked correctly
    DEBUG2(cout << "this is FinalRR"<<endl);
    edge_p inbridge;
    if(!get_final_RR_application(I.g, inbridge)) return false;
    vertex_p u(inbridge->head);

    DEBUG2(cout << "got component-tree leaf "<<u<<endl);

    // get all neighbors of u along non-bridges
    edgelist nh(u->get_cyclic_neighbors());
    if(nh.size() != 2) FAIL("epic fail in final_RR");

    // the two cyclic neighbors of u are v and w
    vertex_p w(nh.back()->head);
    vertex_p v(nh.front()->head);

    // make sure that v is dirty whenever w is dirty
    if(is_dirty(w) && !is_dirty(v)) swap(v, w);

    DEBUG2(cout << "after swap, got v = "<<v<<" w = "<<w<<endl);
    // apply the actual reduction:
    edge_p a;
    // if v cyclic degree of v is exactly 2 and w is dirty, then delete farthest edge from u; otherwise, delete {u,v}
    if((v->non_bridge_degree() == 2) && is_dirty(w)){
      a = get_farthest_edge_on_cycle(u);
    } else {
      a = u->first_non_bridge_neighbor_except(w);
    }
    DEBUG3(cout << "RR1: edge deletion: "<<*a<<endl);
    v = a->head;
    w = a->get_tail();
    I.delete_edge(a, sol);

    // apply TRRs upwards after deleting the edge a
    sol += apply_trrs_upwards_after_cut(I, stat, v, w);
    sol += apply_trrs_upwards_after_cut(I, stat, w);

    DO_STAT(stat.reduct_application[Fin]++);
    return true;
  }




  // recurse for a specific pendant to be hung to v
  solution_t recurse_for(instance& Ismall,
                         const vertex_p v,
                         void (*modify)(graph&, const vertex_p&, const string&),
                         stats_t& stat,
                         const solv_options& solv_opts,
                         const uint depth)
  {
    // 1. copy instance
    unordered_map<uint, vertex_p> id_to_vertex;
    instance J(Ismall, &id_to_vertex);
    // 2. modify instance
    modify(J.g, id_to_vertex[v->id], v->name);
    // 3. recurse
    DEBUG2(cout << "recursing for "<<J.g<<endl);
    solution_t S(run_branching_algo(J, stat, solv_opts, depth+1));
    if(J.g.vertices.empty() && (J.k >= 0)) return S; else return solution_t();
  }

  // Note: technically, this is not a reduction rule, but a branching rule. Hence, we'll need the solv_options
  // the following threshold indicates how big the FES must be on both sides in order to apply
#define BBRule_global_FES_threshold 4
  // apply the B-bridge rule:
  // let b = uv be a b-bridge separating C1 from C2 (u in C1, v in C2)
  // S1: solution of C1 with pendant(v) = singleton
  // S2: solution of C1 with pendant(v) = leaf
  // S3: solution of C1 with pendant(v) = P2
  // S4: solution of C1 with pendant(v) = Y
  solution_t apply_Bbridge_rule(instance& I, stats_t& stat, const solv_options& solv_opts, const uint depth){
    // first, get me a good B-bridge: the score of a B-bridge is how evenly it splits the graph G
    // get the list of Bbridges, weighted by the FES in one of the two(!) components they split off
#warning TODO: weight edges correctly (FES instead of just |V| and min of both sides)

    // for correctness of the weighted B-bridges, we need that G is connected
    assert(I.g.cc_number == 1);

    // if the global FES is below the threshold, don't apply us
    const uint big_FES = get_FES(I.g);
    if(big_FES < BBRule_global_FES_threshold) return solution_t();

    weighted_edges B(I.g.get_weighted_Bbridges());
    // if there are no B-bridges, then just retrn here
    if(B.empty()) return solution_t();
    DEBUG2(cout << "found B-bridges: "<<B<<endl);
    // get the best possible B-bridge
    uint best_score = 0;
    edge_p uv = B.begin()->first;
    for(auto e: B){
      if(e.second < I.g.vertices.size()/2){
        if(e.second > best_score){
          best_score = e.second;
          uv = e.first;
        }
      } else { // if e's component is larger than |V|/2, then take e's reverse
        const uint score(I.g.vertices.size() - e.second);
        if(score > best_score){
          best_score = score;
          uv = e.first->get_reversed();
        }
      }
    }
    // if it is, we've found a B-bridge to apply ourselves to
    DEBUG2(cout << "Bbridge rule: found "<<uv<<" with score "<<best_score<<" (at "<<I.g.vertices.size()<<" vertices)"<<endl);
//     // FES on both sides of the bridge should be reasonably big - !! NO, its fine even if the FES is small !!
//    if(best_score < BBRule_application_threshold) return solution_t();
   
    // 0. get G[X]
    const vertex_p u(uv->get_tail());
    vertex_p v(uv->head);
    const uint v_id(v->id);
    const bool uv_was_permanent(uv->is_permanent);


    // delete uv
    I.g.delete_edge(uv);  // from here, the component containing v is smaller than the other one
    // extract the component containing v as Ismall
    instance Ismall;
    unordered_map<uint, vertex_p>* id_to_vertex = new unordered_map<uint, vertex_p>();
    I.g.copy_component(v, Ismall.g, id_to_vertex);
    DEBUG1(cout << "got translation matrix "<<*id_to_vertex<<" & looking for "<<v<<" (id: "<<v->id<<")"<<endl);
    I.g.delete_component(v);
    v = id_to_vertex->at(v_id);
    Ismall.k = I.k;
    delete id_to_vertex;

    // start computing S1,..,S4
    solution_t S1, S2, S3, S4;
    
    DO_STAT(uint fes[5]; fes[0] = get_FES(Ismall.g); fes[1] = fes[2] = fes[3] = fes[4] = big_FES - fes[0] );
    DO_STAT(unsigned char created_instances = 2);
    // 0. compute just any optimal solution
    DEBUG2(cout << "Step 0: getting S4 "<<endl);
    S4 = recurse_for(Ismall, v, &add_leaf, stat, solv_opts, depth);
    DEBUG2(cout << "Step 0: got solution S4 = "<<S4<<endl);
    if(S4.empty()) {I.k = -1; return solution_t();}
    // I'm only interested in solutions matching this bound for Ismall
    Ismall.k = S4.size();

    // 1. see if we can solve G[X]-uv with less edge deletions than needed for G[X] (only if uv isn't permanent)
    if(!uv_was_permanent) {
      DO_STAT(created_instances++);
      Ismall.k = S4.size() - 1;
      DEBUG2(cout << "Step 1: getting S1 "<<endl);
      S1 = recurse_for(Ismall, v, &add_nothing, stat, solv_opts, depth);
      DEBUG2(cout << "Step 1: got solution S1 = "<<S1<<endl);
      if(!S1.empty()){
        S1 += u->name + "->" + v->name;
        I.k -= S1.size();
        DEBUG2(cout << "Step 1: continuing to solve the rest with S1 = "<<S1<<" and remaining k = "<<I.k<<endl);
        S1 += run_branching_algo(I, stat, solv_opts, depth+1);
        DO_STAT(stat.add_BRule(Bbridge, fes, created_instances));
        return S1;
      } else Ismall.k++; // if this branch failed, restore Ismall's k value
    }
    
    // if we're still here, then no optimal solution contains uv
    // 2. see if some optimal solution contains Z (note that this is only possible if there are no permanent edges incident to v, except uv)
    edge_p perm = v->adj_list.begin();
    while(perm != v->adj_list.end()) if(perm->is_permanent) break; else ++perm;
    if(perm == v->adj_list.end()){
      // okay, no permanent edges found
      DO_STAT(created_instances++);
      DEBUG2(cout << "Step 2: getting S2 "<<endl);
      S2 = recurse_for(Ismall, v, &add_Y, stat, solv_opts, depth);
      DEBUG2(cout << "Step 2: got solution S2 = "<<S2<<endl);
      if(!S2.empty()){
        I.k -= S2.size();
        add_leaf(I.g, u);
        DEBUG2(cout << "Step 2: continuing to solve the rest with S2 = "<<S2<<" and remaining k = "<<I.k<<endl);
        S2 += run_branching_algo(I, stat, solv_opts, depth+1);
        DO_STAT(stat.add_BRule(Bbridge, fes, created_instances));
        return S2;
      }
    } else DEBUG2(cout << "Step 2: all edges around "<<v<<" are permanent"<<endl);

    // if we're still here, then no optimal solution contains uv or Z
    // 3. see if some optimal solution contains Z-ux for some x
    DO_STAT(stat.add_BRule(Bbridge, fes, ++created_instances));
    DEBUG2(cout << "Step 3: getting S3 "<<endl);
    S3 = recurse_for(Ismall, v, &add_P2, stat, solv_opts, depth);
    DEBUG2(cout << "Step 3: got solution S3 = "<<S3<<endl);
    if(!S3.empty()){
      I.k -= S3.size();
      add_P2(I.g, u);
      DEBUG2(cout << "Step 3: continuing to solve the rest with S3 = "<<S3<<" and remaining k = "<<I.k<<endl);
      S3 += run_branching_algo(I, stat, solv_opts, depth+1);
      return S3;
    }
    
    // if none of the others finished, S4 is going to be my solution
    // note that, since |S3|>|S4|, we know that there is a Y-graph dangling at u after deleting S4
    I.k -= S4.size();
    add_Y(I.g, u);
    S4 += run_branching_algo(I, stat, solv_opts, depth+1);
    return S4;
  }


  // returns whether there is a generator-free path from u to v avoiding x
  // if u and v are adjacent and we want to forbid going to v, set x = u
  bool exists_gen_free_path(const vertex_p& u, const vertex_p& v, const vertex_p& x){
    // for good measure, check whether we can even reach u and v
    if(!u->trr_infos.ptwos.empty() || !v->trr_infos.ptwos.empty()) return false;
    // can't use graphs DFS-id since we must keep the dfs id's in tact for PRRs
    unordered_set<vertex_p, vertex_hasher> visited;
    list<vertex_p> path_from;
    // if x == u, then we mean that the edge uv should not be taken, so add all neighbors of u except v
    if(x == u){
      visited.insert(u);
      for(auto &e : u->adj_list)
        if(e.head != v) path_from.push_back(e.head);
    } else path_from.push_back(u);
    // do BFS to find v
    while(!path_from.empty()){
      const vertex_p w(path_from.front());
      path_from.pop_front();

      // mark w
      visited.insert(w);
      // see if we're x or a generator
      if((w == x) || !w->trr_infos.ptwos.empty()) continue;
      // if we have reached v, then return success
      if(w == v) return true;

      // add the non-visited neibhbors of w
      for(auto &e : w->adj_list)
        if(visited.find(e.head) == visited.end())
          path_from.push_back(e.head);
    }
    return false;
  }

  // generalized PRR4
  bool prr4_gen_applicable(const vertex_p& separator){
    // 1. if we're at a 2P2, PRR4 is applicable
    if(separator->trr_infos.ptwos.size() > 1) return true;
    // 2. get neighbors of x
    const edgelist el(separator->get_cyclic_core_neighbors());
    assert(el.size() == 2);
    // 3. find a path from one to the other avoiding generators and x
    return !exists_gen_free_path(el.front()->head, el.back()->head, separator);
  }

  // new TRR3:
  // if |P| > 1, then for all vw such that w is not the center of any of the first & last P2,
  // delete vw and add a Y-graph to v
  bool trr3_gen(instance& I, stats_t& stat, const vertex_p& v, solution_t& sol){
    // check applicability
    if(!v->is_on_cyclic_core()) return false;
    if(v->trr_infos.ptwos.size() < 2) return false;

    DO_STAT(stat.reduct_application[TRR3]++);

    // get cyclic neighbors of v
    edgelist cn(v->get_cyclic_core_neighbors());

    // prepare the vertex set that should not be crossed by TRRs
    vertexset do_not_cross;
    for(edge_pc e : cn) do_not_cross.insert(e->head);

    // next, Ygraphify the non-P2 neighbors
    for(edge_p e : cn){
      edge_p f(e);
      vertex_p u(f->head);

      do_not_cross.erase(u);

      Ygraphify(I.g, f);
      sol += apply_trrs_upwards_after_cut(I, stat, u, do_not_cross);
    }
    return true;
  }


  // Y-graph lookahead:
  // if v's pendant vw is a Y-graph:
  // 1. if deg(v) - 1 exceeds the upper bound, then delete vw
  // 2. if deg(v) - 1 matches the upper bound, then try to delete all but vw and apply TRR6.
  //    If that dosen't give the empty graph, delete vw
  bool Y_lookahead(instance& I, stats_t& stat, solution_t& sol, const vertex_p& v, const uint upper_bound){
    if(!v->trr_infos.ygraphs.empty()){
      DEBUG2(cout << "testing "<<v<<" for YL"<<endl);
      const uint cyc_v_deg(v->degree() - 1);
      if(cyc_v_deg < upper_bound)
        return false;
      else if(cyc_v_deg == upper_bound){
        // make a copy I' of I
        unordered_map<uint, vertex_p> id_to_vertex;
        instance Iprime(I, &id_to_vertex);
        const vertex_p vprime(id_to_vertex[v->id]);
        const vertex_p wprime(id_to_vertex[v->trr_infos.ygraphs.front()->head->id]);
        // delete N(v) - Y from I'
        for(edge_p f = vprime->adj_list.begin(); f != vprime->adj_list.end();) if(f->head == wprime) ++f; else
          Iprime.g.delete_edge(f++);
        // see if this solved the instance
        trr6(Iprime);
        // if not, delete vw in I
        if(!Iprime.g.vertices.empty()){
          DEBUG2(cout << "This is YL for "<<v<<" with upper bound "<<upper_bound<<endl);
          I.delete_edge(v->trr_infos.ygraphs.front(), sol);
          v->trr_infos.ygraphs.pop_front();
          DO_STAT(stat.reduct_application[YL]++);
          return true;
        } else return false;
      } else {
        DEBUG2(cout << "This is YL for "<<v<<" with upper bound "<<upper_bound<<endl);
        // here, cyc_v_deg > upper_bound, so delete vw
        I.delete_edge(v->trr_infos.ygraphs.front(), sol);
        v->trr_infos.ygraphs.pop_front();
        DO_STAT(stat.reduct_application[YL]++);
        return true;
      }
    } else return false;
  }

  // apply to all vertices if the graph is small enough to allow success
  bool Y_lookahead(instance& I, stats_t& stat, const solv_options& opts, solution_t& sol, const uint upper_bound){
    if(I.g.vertices.size() > opts.max_size_for_Y_lookahead) return false; else {
      bool result = false;
      for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v){
        if(Y_lookahead(I, stat, sol, v, upper_bound)) result = true;
      }
      return result;
    }
  }
  // use I.k as upper bound if no upper bound is given
  bool Y_lookahead(instance& I, stats_t& stat, const solv_options& opts, solution_t& sol){
    return Y_lookahead(I, stat, opts, sol, I.k);
  }




}

