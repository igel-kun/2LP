#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "../util/b_vector.hpp"
#include "../reduction/trr.hpp"
#include "../reduction/prr.hpp"
#include "../reduction/global.hpp"
#include "bounds.hpp"
#include "branching.hpp"

#include <algorithm> // for sort
#include <unordered_map>
#include <unordered_set>


namespace cr{

  edge_p skip_deg2path(edge_p direction, const vertex_p& do_not_cross, uint *length = NULL){
    DEBUG2(cout << "skipping non-intrusively in direction "<<*direction<<endl);
    while(direction->head != do_not_cross){
      const vertex_p& to(direction->head);
      
      // if we looped, then return
      if(to == do_not_cross) return direction;
      // if a cyc-deg-3 vertex is reached, return it
      if(to->cyc_core_degree() != 2) return direction;

      // if there are no more cyclic neighbors, return the current vertex
      direction = get_next_on_deg2path(direction);

      // advance the length-counter
      if(length) ++(*length);
    }
    return direction;
  }


  // return true if it contains permanent edges
  bool contains_permanent(edgelist& el){
    for(edgelist::const_iterator e = el.begin(); e != el.end(); ++e)
      if((*e)->is_permanent){
        return true;
      }
    return false;
  }

  // return whether deleting e1 makes e2 an A-bridge
  inline bool del_one_makes_other_Abridge(graph& g, const edge_p& e1, const edge_p& e2){
    const edge_p e(skip_deg2path(e1, e2->get_tail()));
    return (e->get_reversed() == e2);
  }


  // a triangle of edges w <--b-- v --a--> u --c--> v is degenerate...
  // if v is cyclic-deg-2 and incident to a Y-graph,
  // OR
  // if u, and w are on an all-singleton v-v-(deg-2 path)
  //
  // in degenerate triangles, a (=uv) & b (=vw) are deleted together
  inline bool triangle_is_degenerate(const vertex_p& v, const vertex_p& u, const vertex_p& w){
    if(v->cyc_core_degree() == 2)
      if(!v->trr_infos.ygraphs.empty())
        return true;
    if(u->cyc_core_degree() == 2)
      if(w->cyc_core_degree() == 2)
        if(v->trr_infos.count() == 0)
          if(u->trr_infos.count() == 0)
            if(w->trr_infos.count() == 0)
              return true;
    return false;
  }

  // get branching ops for BRR1
  // TODO: this is not optimal! maintain a list of triangles in G instead!
  // TODO: modify for triangles with just one deg-3-vertex!
  bool BRR1(const vertex_p v, branchlist& br){
    edgelist non_bridges(v->get_cyclic_neighbors());
    DEBUG2(cout << "this is BRR1 for "<<*v<<" with "<<v->adj_list.size()<< " neighbors"<< endl);
    for(auto i = non_bridges.begin(); i != non_bridges.end(); ++i){
      const edge_p& a(*i);
      for(auto j = i; j != non_bridges.end(); ++j){
        const edge_p& b(*j);
        if(a != b){
          edge_p c = find_edge(a->head, b->head);
          if(c != a->head->adj_list.end()){
            branch_op bo(Triangle);
            if(triangle_is_degenerate(v, a->head, b->head)){
              // register {a,b}-branch
              if(!a->is_permanent && !b->is_permanent) bo.branches.push_back({graph_mod_t(a), graph_mod_t(b)});
              // register {c}-branch
              if(!c->is_permanent) bo.branches.push_back({c});
            } else {
              if(!a->is_permanent) bo.branches.push_back({a});
              if(!b->is_permanent) bo.branches.push_back({b});
              if(!c->is_permanent) bo.branches.push_back({c});
            }
            DEBUG2(cout << "found 3CYC on "<<*a<<", "<<*b<<", "<<*c<<", non-permanent: "<<bo.branches<<endl);
            if(!bo.branches.empty()){
              bo.bnum = bo.branches.size();
              br.push_back(bo);
              // one BRR1 for each vertex is enuff, it's not getting better branching numbers anyway
              return true;
            }
          } // end if
        } // end if
      } // end for
    } // end for
    return false;
  } // end function


  edgelist::const_iterator get_nth(const edgelist& el, uint n){
    edgelist::const_iterator i = el.begin();
    while(n){
      ++i;
      --n;
    }
    return i;
  }
  edgelist::const_iterator operator+(const edgelist::const_iterator i, uint n){
    edgelist::const_iterator j = i;
    while(n){
      ++j;
      --n;
    }
    return j;
  }
  edgelist::iterator operator+(const edgelist::iterator i, uint n){
    edgelist::iterator j = i;
    while(n){
      ++j;
      --n;
    }
    return j;
  }

  // return whether BRR2 applies to the v_i vertices
  bool BRR2(const list<claw_leg>& legs, branchlist& br){
    // BRR2 requires all |E_i|=1
    branch_op bop(Claw0);
    edgelist heads;
    for(auto i = legs.begin(); i != legs.end(); ++i) 
      if(i->E.size() == 1) {
        add_branch(bop, i->E); 
        heads.push_back(i->head);
      } else return false;

    /* new branching deletes {e1,e2,e3}
    // edges {e1,e2}
    if(nonempty_branches > 1)
      add_branch(bop, { *(vis.begin()), *(vis.begin()+1) });
    */
    if(heads.size() > 2) add_branch(bop, heads);


    DEBUG2(cout << "BRR2 computed branch op "<<bop<<endl);
    br.push_back(bop);
    return true;
  }
  // return whether BRR3 applies to the v_i vertices, assuming that BRR2 does not apply
  bool BRR3(const list<claw_leg>& legs, branchlist& br){
    uint big_vertices = 0;
    branch_op bop(Claw1);
    edgelist small_heads;

    for(auto i = legs.begin(); i != legs.end(); ++i)
      if(i->E.size() > 1){
        // if we've already seen a big vertex, then fail
        if(++big_vertices > 1) return false; else {
          add_branch(bop, i->E);
          add_branch(bop, { i->head });
        }
      } else {
        // we're at a deg-2 vertex
        add_branch(bop, i->E);
        small_heads.push_back(i->head);
      }
    // edge {e2,e3}
    if(small_heads.size() > 1) add_branch(bop, small_heads);

    DEBUG2(cout << "BRR3 computed branch op "<<bop<<endl);
    br.push_back(bop);
    return true;
  }
  // return whether BRR4 applies to the v_i vertices, assuming that BRR2-3 does not apply
  bool BRR4(const list<claw_leg>& legs, branchlist& br){
    uint big_vertices = 0;
    branch_op bop(Claw2);

    for(auto i = legs.begin(); i != legs.end(); ++i)
      if(i->E.size() > 1){
        // if we've already seen two big vertex, then fail
        if(++big_vertices > 2) return false; else {
          add_branch(bop, i->E);
          add_branch(bop, { i->head });
        }
      } else add_branch(bop, i->E);

    DEBUG2(cout << "BRR4 computed branch op "<<bop<<endl);
    br.push_back(bop);
    return true;
  }
  // return whether BRR5 applies to the v_i vertices
  bool BRR5(const list<claw_leg>& legs, branchlist& br){
    branch_op bop(Claw3);

    for(auto i = legs.begin(); i != legs.end(); ++i){
      add_branch(bop, i->E);
      add_branch(bop, { i->head });
    }

    DEBUG2(cout << "BRR5 computed branch op "<<bop<<endl);
    br.push_back(bop);
    return true;
  }

  // return whether there are separators on the deg2path starting with e
  bool skip_deg2path_finding_separators(edge_p& e, const vertex_p& do_not_cross){
    bool result = false;
    while((e->head != do_not_cross) && (e->head->cyc_core_degree() < 3)){
      if(e->head->is_separator()) {
        DEBUG2(cout << "skip_deg2path found separator "<<e->head<<" with trr: "<<e->head->trr_infos<<endl);
        result = true;
      }
      e = get_next_on_deg2path(e);
    }
    return result;
  }
  // get clean neighbors, meaning that for each two edges whose deletion makes the other an Abridge, keep just one
  // also, remove leaves and Abridges
  void get_clean_neighbors_and_disallowed(const vertex_p& v, edgelist& clean_nh, edge_p& disallowed){
    vertexset used;
    disallowed = v->adj_list.end();
    const bool has_P2_pendant(!v->trr_infos.ptwos.empty());

    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e) if(!e->is_Abridge()){
      DEBUG2(cout << "considering "<<e<<endl);
      if(e->is_bridge || has_P2_pendant){
        // use B-bridges
        used.insert(e->head);
        clean_nh.push_back(e);
        edge_p f(e);
        if(disallowed == v->adj_list.end()){
          if(skip_deg2path_finding_separators(f, v)) disallowed = e;
        }
      } else {
        // see if we come back to v from e
        edge_p f(e);
        bool has_separator = skip_deg2path_finding_separators(f, v);
        
        if((f->head != v) || (used.find(f->get_tail()) == used.end())){
          // if we don't come back or we come back to v, or its the first time we go along this path, then use e
          used.insert(e->head);
          clean_nh.push_back(e);
          if(has_separator)
            if(disallowed == v->adj_list.end()) disallowed = e;
        }
      }
    }
    DEBUG2(cout << "found "<<clean_nh<<endl);
  }

  // Branching Rule 6:
  // guess two legs of the caterpillar at v & Y-graphify all else
  // if v is not on the backbone in some optimal solution, then v is a leaf, so guess the backbone neighbor of v
  bool BRR6(const vertex_p& v, branchlist& br){
    DEBUG2(cout << "checking "<<v<<" with cyclic_core_degree "<< v->cyc_core_degree()<<" and nldeg "<< v->nldeg()<<endl);
    if(v->is_on_cyclic_core() && (v->nldeg() > 2) ){
      branch_op bop(Token);
      // get the ptwos
      const edgelist& ptwos(v->trr_infos.ptwos);
      const unsigned char has_P2(ptwos.empty() ? 0 : 1);
      // assert by reducedness wrt. trr3_gen that |P2|<2
      DEBUG2(assert(ptwos.size() < 2));
      // get a copy of all incident edges
      list<edgeset> keep_legs;
      edgeset to_keep;

      // get clean neighbors, meaning that for each two edges whose deletion makes the other an Abridge, keep just one
      // also, remove leaves and Abridges
      edgelist clean_nh;
      edge_p disallowed;
      get_clean_neighbors_and_disallowed(v, clean_nh, disallowed);
      const bool has_disallowed(disallowed != v->adj_list.end());
      // if we have just 2 interesting neighbors, then this is not going to work
      if(clean_nh.size() + has_P2 < 3) return false;

      DEBUG2(cout << "this is BRR6 for "<<v<<" with "<<" P2: "<<(uint)has_P2; if(has_disallowed) cout<<" disallowed: "<< disallowed; cout << " clean NH: "<<clean_nh<<endl);

      // build a list of legs to keep into keep_legs, all others are to be deleted
      if(!(has_P2 && has_disallowed)){ // if there is just 1 p2, then we need another leg
        for(edge_pp keep1 = clean_nh.begin(); keep1 != clean_nh.end(); ++keep1) if(*keep1 != disallowed){
          to_keep.insert(*keep1);
          if(!has_P2 && !has_disallowed){ // if there is no p2 at all, we need a second leg
            for(edge_pp keep2 = keep1; ++keep2 != clean_nh.end();){
              to_keep.insert(*keep2);
              keep_legs.push_back(to_keep);
              to_keep.erase(*keep2);
            }
          } else keep_legs.push_back(to_keep);
          to_keep.erase(*keep1);
        } 
      } else keep_legs.push_back(to_keep);
      DEBUG2(cout << "Token rule: keeping legs in "<<keep_legs<<endl);

      // for each leg, invert it to get the branching op
      for(auto leg = keep_legs.begin(); leg != keep_legs.end(); ++leg){
        // add 'disallowed' to each leg
        if(has_disallowed) leg->insert(disallowed);

        edgelist branch;
        for(edge_pp e = clean_nh.begin(); e != clean_nh.end(); ++e)
          if(leg->find(*e) == leg->end()) branch.push_back(*e);
        
        DEBUG2(cout << "creating branch for "<<*leg<<": "<<branch<<endl);
        add_branch(bop, branch, Yify);
      }
      // if v has no leaves or ptwos, then v might not be on the backbone but instead become a leaf!
      if(!v->is_on_backbone()){
        // if v becomes a leaf, then guess its neighbor and delete the rest
        // first, check if there is a permanent edge incident to v
        edge_p perm = v->adj_list.end();
        for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
          if(e->is_permanent) {
            if(perm != v->adj_list.end()) {
              // v is incident to TWO permanent edges, then it cannot become a leaf!
              DEBUG2(cout <<"BRR6 created branch operation: "<<bop<<" with "<<bop.branches.size()<<" branches"<<endl);
              br.push_back(bop);
              return true;
            } else perm = e; // if we found the first permanent edge, remember it
          }
        if(perm != v->adj_list.end()){
          // if there is a permanent edge, this edge is the one to the backbone, so keep it
          edgelist el;
          for(edge_p f = v->adj_list.begin(); f != v->adj_list.end(); ++f) if(f != perm) el.push_back(f);
          add_branch(bop, el, Del);
        } else {
          for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e){
            // get a list of all incident edges of v except e
            edgelist el;
            for(edge_p f = v->adj_list.begin(); f != v->adj_list.end(); ++f) if(f != e) el.push_back(f);
            // add el as new branch
            add_branch(bop, el, Del);
          }
        }
      }
      DEBUG2(cout <<"BRR6 created branch operation: "<<bop<<" with "<<bop.branches.size()<<" branches"<<endl);
      br.push_back(bop);
      return true;
    } else return false;
  }


  bool BB_branching_applicable(const path_info_t& info){
    // there cannot be generators on the path AND
    if(!info.generators.empty()) return false;
    // the end vertices should be on the backbone
    const vertex_p u(info.start->get_tail());
    const vertex_p v(info.end->head);
    if(!u->is_on_backbone()) return false;
    if(!v->is_on_backbone()) return false;
    // see if u and v are weakly separated in G-p, where p is the deg2-path we're considering
    const vertex_p x(info.start->head);
    // note that, if |p|=1, then x == v but exists_gen_free_path can handle that if we call it with (v, u, x)
    return !exists_gen_free_path(v, u, x);
  }

  bool BRR78(const path_info_t& info, branchlist& br){
    DEBUG2(cout<<"this is BRR78 for the path "<<info<<endl);
    if(BB_branching_applicable(info)){
      DEBUG2(cout<<"we're applicable"<<endl);
      branch_op bop(Deg2Path);

      if(info.length > 1){
        // advance info.start by one if info.start->head is not a separator
        const edge_p to_del_left( info.start->head->is_separator() ? info.start : get_next_on_deg2path(info.start));
        add_branch(bop, edgelist(1, to_del_left));
        // note that we need to reverse info->end in order for the special branch to delete the deg-2 path
        if(info.length > 2){
          const edge_p rev_end(info.end->get_reversed());
          if(rev_end->head->is_separator()){
            add_branch(bop, edgelist(1, rev_end));
          } else if(rev_end->head != to_del_left->head) { // if we have a len-3 path without separators, don't create a second branch
            add_branch(bop, edgelist(1, get_next_on_deg2path(rev_end)));
          }
        }
      } else add_branch(bop, edgelist(1, info.start));
      // add an existing branch to represent the guess that a caterpillar includes this deg2-path.
      // This should work with all branch-number computations
      if(!bop.branches.empty()) bop.branches.push_back(bop.branches.front());

      DEBUG2(cout << "created branch-op "<<bop<<endl);
      br.push_back(bop);
      return true;
    } else return false;
  }

  void apply_BRR78(instance& I, const edgelist& el, solution_t& sol){
    edge_p to_del(el.front());
    DEBUG2(cout << "applying BRR78 to "<<el<<endl);
    if(to_del->is_permanent){
      // if our branch is permanent, then this is the special branch that cuts the deg-2 path
      const vertex_p u(to_del->get_tail());
      while((to_del->head->cyc_core_degree() < 3) && (to_del->head != u)){
        const edge_p next(get_next_on_deg2path(to_del));
        I.g.delete_edge(to_del);
        to_del = next;
      }
      // finally, we arrived at the end of the deg-2 path
      const vertex_p v(to_del->head);
      
      DEBUG2(cout << "BRR78: special branch; deleting deg-2 path between "<< u<<" and "<<v<< endl);
      I.g.delete_edge(to_del);
      // after having cut the deg-2 path, add pendant P2's to u and v
      add_P2(I.g, u);
      add_P2(I.g, v);
      // and that's it, the solution remains untouched
    } else {
      // if we don't have a separator, then advance the edge to delete by one
      if((to_del->head->cyc_core_degree() < 3) && !to_del->head->is_separator())
        to_del = get_next_on_deg2path(to_del);
      I.delete_edges(edgelist(1, to_del), sol);
    }
  }








  // note:
  // in A bridges, one endpoint is the parent of the other!
  // a B-bridge is a bridge that is not an A-bridge

  // return whether e points to an eligible branching head
  bool is_eligible_branching_head(const edge_p& e){
    DEBUG2(cout<<"testing "<<*e<<endl);
    // branching into y-graphs on deg2-paths is not good! (e is non-relevant Abridge)
    const vertex_p& v(e->get_tail());
    if(e->head->degree() == 1) return false;
    DEBUG2(cout<<"passed 1st test: "<<*(e->head)<<" is not a leaf"<<endl);
    if((v->non_bridge_degree() == 2) && (!v->trr_infos.ygraphs.empty()))
      if(v->trr_infos.ygraphs.front() == e) return false;
    DEBUG2(cout<<"passed 2nd test: "<<*e<<" is not a relevant A-bridge"<<endl);
//  note: killing B-bridges is, in fact, okay with the analysis, since it splits components off the graph
//    if(e->head->is_incident_to_Bbridge()) return false;
//    DEBUG2(cout<<"passed 3nd test: "<<*(e->head)<<" is not incident to B-bridges"<<endl);
    return true;      
  }

  void get_non_Abridge_branching_heads(const vertex_p& v, edgelist& el){
    DEBUG2(cout << "getting non-bridge branching heads of "<<*v<<endl);
    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
      if(!e->is_Abridge())
        if(is_eligible_branching_head(e))
          el.push_back(e);

    DEBUG2(cout << "found " << (el.empty()?"none":"these: ") << el << endl);
  }
  // select three branching heads to fulfill the following conditions:
/*
 *  (1) deg(v1) = 2 ==> deg(v2) = 2 and deg(v2) = 2 ==> deg(v3) = 2.
 *  (2) if deg(v1) = 2, then deleting {v, v1} does not make {v, v2} a bridge, and
 *  (3) if deg(v1) > deg(v2) = 2, then deleting {v, v2} does not make {v, v3} a bridge.
 */
  bool first_is_better_branchhead(const edge_p& a, const edge_p& b){
    // sorting by degree is not best. Instead, deg-1 is best but for deg>1, the bigger the better!
    const vertex_p u(b->head);
    const vertex_p v(a->head);
    // 1. prefer deg-2 vertices
    if(v->degree() == 2) return true;
    if(u->degree() == 2) return false;
    // 2. then, prefer vertices with permanent edges around them
    if(!v->trr_infos.leaves.empty() ||  !v->trr_infos.ptwos.empty()) return true;
    if(!u->trr_infos.leaves.empty() ||  !u->trr_infos.ptwos.empty()) return false;
    // 3. then, prefer vertices with high degree
    return v->degree() > b->head->degree();
  }
  bool greater_than(const edge_p& a, const edge_p& b){
    return a->head->degree() > b->head->degree();
  }
  bool first_is_better_branch(const modlist_t& E1, const modlist_t& E2){
    // first, prefer size-1 edgelists
    if(E1.size() == 1) return true;
    if(E2.size() == 1) return false;
    // then, prefer smaller over larger
    return E1.size() < E2.size();
  }


  // if deleting *i makes *j a bridge, then replace *j with *(++j)
  void bring_in_order(edgelist& el, const edge_pp& i, edge_pp j){
    const vertex_p& v((*i)->get_tail());
    edge_p e(*i);
    // find out whether going along (v,v2) returns to v with (v3,v)
    DEBUG2(cout << "checking whether killing "<<**i<<" makes "<<**j<<" a bridge"<<endl);
    skip_deg2path(e, v);
    if(e->get_reversed() == *j){
      DEBUG2(cout << "indeed so! Maybe there is more on the edgelist: "<<el<<" ("<<el.size()<<" items)"<<endl);
      DEBUG2(cout << "i is "<<**i<<endl);
      DEBUG2(cout << "j is "<<**j<<endl);
      // swap the second and third items on el
      // NOTE: if there are no more items on the list, then v4 is incident to a B-bridge,
      // and therefore, not an eligible branching head - in this case, we cannot branch here
      edgelist::iterator k = j;
      ++k;
      if(k == el.end()){
        el.clear();
      } else {
        swap(*j, *k);
      }
    }
  }
  // TODO: preferably get ptwos as branching heads
  bool select_branching_heads(edgelist& el){
    const vertex_p& v(el.front()->get_tail());
    const uint ptwos(v->trr_infos.ptwos.size());

    // if there is not enuff to branch on, return false
    DEBUG2(cout << "totalling "<<el.size()+ptwos<<" possible branching heads"<<endl);
    if(el.size() + ptwos < 3) return false;
    // if we have at least 3 ptwos, then TRR's have not been applied!
    assert(ptwos < 3);

    DEBUG2(cout << "sorting by degree"<<endl);
//    sort(el.begin(), el.end(), greater_than);
    el.sort(first_is_better_branchhead); // sort by decreasing quality (best first)
    DEBUG2(cout << "new order: "<<el<<endl);

    if(v->trr_infos.ptwos.empty()){
      
      edgelist::iterator i = el.begin();
      // ensure (2)
      if(el.front()->head->degree() == 2){
        DEBUG2(cout<<"ensuring (2) for "<<**i<<endl);
        edgelist::iterator j = i;
        bring_in_order(el, i, ++j);
      };
      // bring_in_order can determine that the branching heads are invalid. In this case, return
      if(el.empty()){
        DEBUG2(cout << "branching heads were not valid, returning with empty edgelist"<<endl);
        return false;
      }
      ++i;
      // ensure (3)
      if((el.front()->head->degree() > 2) &&
        ((*i)->head->degree() == 2)){
        DEBUG2(cout<<"ensuring (3) for "<<**i<<endl);
        edgelist::iterator j = i;
        bring_in_order(el, i, ++j);
      }
      // bring_in_order can determine that the branching heads are invalid. In this case, return
      if(el.empty()){
        DEBUG2(cout << "branching heads were not valid, returning with empty edgelist"<<endl);
        return false;
      }
    } else DEBUG2(cout << *v << " has "<<ptwos<<" P2s, using those"<<endl);

    el.erase(el.begin() + (3 - ptwos), el.end());
    DEBUG2(cout<<"everything OK, branching heads are now "<<el<<endl);
    return true;
  }
  // returns true if there is a path from u to v avoiding do_not_cross and returns the last_edge on the path
  bool find_path_avoiding(const vertex_p& u, const vertex_p& v, const edge_p& do_not_cross, const size_t dfs_id, edge_p& last_edge){
    if(u->dfs_id == dfs_id) return false;

    u->dfs_id = dfs_id;
    for(edge_p e = u->adj_list.begin(); e != u->adj_list.end(); ++e) if(e != do_not_cross){
      if(e->head == v){
        last_edge = e;
        DEBUG2(cout << "heureka, found path from "<<*u<<" to "<<*v<<" avoiding "<<*do_not_cross<<". Last edge is "<<*e<<endl);
        return true;
      } else if(find_path_avoiding(e->head, v, do_not_cross, dfs_id, last_edge)) return true;
    }
    return false;
  }

  
  // new version of compute_Ei that returns a set containing _all_ edges incident to vi except {v,vi}
  void compute_Ei(graph& g, const edge_p& ei, edgelist& el){
    const vertex_p& vi(ei->head);
    DEBUG2(cout << "finding E_i for "<<*ei<<" (non-bridge-degree of "<<*vi<<": "<<vi->non_bridge_degree()<<")"<<endl);

    for(edge_p e = vi->adj_list.begin(); e != vi->adj_list.end(); ++e)
      if(e != ei->get_reversed())
        el.push_back(e);
  }

  
  // get branching ops for BRR2
  bool BRR2_to_5(graph& g, const vertex_p& v, branchlist& br){
    edgelist branch_heads;
    if(v->nldeg() < 3) return false;
    DEBUG2(cout << "This is BRR2-5 for "<<*v<<endl);
    // first, get the branching heads
    get_non_Abridge_branching_heads(v, branch_heads);
    // if we can't get any branching heads, then just return
    if(branch_heads.empty()) return false;
    // then, select the ones that are good
    DEBUG2(cout << "selecting correct branching heads from "<<branch_heads<<endl;);
    if(!select_branching_heads(branch_heads)) return false;

    list<claw_leg> legs;
    for(edgelist::iterator e = branch_heads.begin(); e != branch_heads.end(); ++e){
      edgelist el;
      // compute the sets E_i and store em in Eis
      compute_Ei(g, *e, el);
      legs.push_back( (claw_leg){*e, el } );
    }
    DEBUG2(cout << "going into branching rules with "<<legs<<endl;);

    if(BRR2(legs, br)) return true;
    if(BRR3(legs, br)) return true;
    if(BRR4(legs, br)) return true;
    if(BRR5(legs, br)) return true;

    FAIL("epic fail in BRR2_to_5");
  }


  // return the number of single-delete-branches in a branch-op
  uint single_branches(const branch_op& b){
    uint result = 0;
    for(auto ml = b.branches.begin(); ml != b.branches.end(); ++ml)
      if((ml->size() == 1) && (ml->front().type == Del)) ++result;
    return result;
  }

  pair<branch_op, float> select_best_branch_from_list(const branchlist& br){
    DEBUG2(cout << "accumulated "<<br.size()<<" branching ops, now choosing the best one"<<endl);
    branchlist::const_iterator best_op;
    float best_bnum = FLT_MAX;
    uint best_single_branches = UINT_MAX;
    for(auto bop = br.begin(); bop != br.end(); ++bop){
      float bnum = bop->bnum;
      if(bnum == 0) bnum = branch_number(*bop); else DEBUG2(cout << "bnum already known: "<<bnum<<endl);
      DEBUG2(cout<<bnum<< "-branching: "<<*bop<<endl);
      
      if(bnum > best_bnum) continue;
      const uint singles(single_branches(*bop));
      if((bnum < best_bnum) || ((bnum == best_bnum) && (singles < best_single_branches))){
        best_bnum = bnum;
        best_op = bop;
        best_single_branches = singles;
      }
    }
    branch_op bo( *best_op);
    // sort the single branches to the beginning
    bo.branches.sort(first_is_better_branch);
    // and exchange first branch for the largest
    bo.branches.push_front(bo.branches.back());
    bo.branches.pop_back();
    DEBUG2(cout << "best branching op is "<<bo<<" with bnum "<<best_bnum<<endl);
    return make_pair(bo, best_bnum);
  }

  bool get_best_branch_op(graph& g, branch_op& bo, const list<path_info_t>& path_infos, const bool quick_select, const float branch_threshold){
    branchlist br;

    // check BRR6 first, it has the best chance to produce a size-1 branching (and if so, it produces the best size-1 branching)
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v)
      if(BRR6(v, br))
        if(br.back().branches.size() == 1){
          bo = br.back();
          return true; // if we found a reduction, we're happy
        }
    DEBUG2(cout << "done applying BRR6"<<endl);
    if(quick_select && !br.empty()) {
      pair<branch_op, float> best_branch(select_best_branch_from_list(br));
      if(best_branch.second <= branch_threshold){
        // continue trying branching rules if we didn't find something with bnum at most two
        bo = best_branch.first;
        return true;
      } else {
        br.clear();
        br.push_back(best_branch.first);
      }
    }

/*
 * NOTE: BRR78 is utter bullshit and I cannot tell what the fuck itactually does
    // try BRR 7 & 8
    DEBUG2(cout << "we got "<<path_infos.size() << " paths to check BRR78 on:"<<endl);
    DEBUG2(cout << path_infos<<endl);
    for(auto p = path_infos.begin(); p != path_infos.end(); ++p)
      if(BRR78(*p, br)){
        if(br.back().branches.size() == 1){
          bo = br.back();
          return true;
        }
      }
    DEBUG2(cout << "done applying BRR7 & 8"<<endl);
    if(quick_select && !br.empty()) {
      pair<branch_op, float> best_branch(select_best_branch_from_list(br));
      if(best_branch.second <= branch_threshold){
        // continue trying branching rules if we didn't find something with bnum at most two
        bo = best_branch.first;
        return true;
      } else {
        br.clear();
        br.push_back(best_branch.first);
      }
    }
    */
  
    // check BRR1
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v) if(v->is_on_cycle()){
      if(BRR1(v, br))
        if(br.back().branches.size() == 1){
          bo = br.back();
          return true; // if we found a reduction, we're happy
        }
    } else DEBUG1(cout << "not applying BRR1 to "<<v<<" since "<<v->degree()<<"-"<<v->incident_bridges<<"=0"<<endl);
    DEBUG2(cout << "done applying BRR1"<<endl);
    // try BRR 2-5
/*    for(vertex_p v = g.vertices.begin(); v != g.vertices.end(); ++v) if(v->is_on_cycle()){
      if(BRR2_to_5(g, v, br))
        if(br.back().branches.size() == 1){
          bo = br.back();
          return true; // if we found a reduction, we're happy
        }
    }
    DEBUG2(cout << "done applying BRR2-5"<<endl);
*/    
    // if we found no branching applications, then final_RR should be applied
    if(br.empty()){
      DEBUG2(cout << "couldn't branch: no branchable vertices"<<endl);
      return false;
    }

    bo = select_best_branch_from_list(br).first;
    return true;
  }



  // apply a Deletion operation
  inline void apply_Del(instance& I, const edge_p& e, solution_t& sol){
    I.delete_edge(e, sol);
  }

  // apply a Ygraphify operation
  void apply_Yify(instance& I, const edge_p& e, solution_t& sol, const string name_add = ""){
    const vertex_p v(e->get_tail());
    const vertex_p u(e->head);
    const string uname(u->name);
    string name(v->name + name_add);

    if(!e->is_permanent){
      // delete the edge
      I.g.delete_edge(e);
      // and add a new Y-graph
      add_Y(I.g, u, name);
    } else {
      // if the edge we want to delete is permanent, we will delete all but e at the head
      const edge_p e_rev(e->get_reversed());
      for(edge_p f = u->adj_list.begin(); f != u->adj_list.end();)
        if(f == e_rev) ++f; else apply_Del(I, f++, sol);
    }
    // add a leaf with some name if v is not already recognizable as on backbone
    if(!v->is_on_backbone()) add_leaf(I.g, v, uname + '*');
  }

  void apply_one_branch(instance& I, const branch_type& t, const modlist_t& ml, solution_t& sol){
    for(auto gmod : ml){
      switch(gmod.type){
        case Del:
          apply_Del(I, gmod.e, sol);
          break;
        case Yify:
          apply_Yify(I, gmod.e, sol);
          break;
      }
    }
  }



  // process a branching operation: for each of its branches, mke a copy of the graph, apply its changes, and recurse
  solution_t apply_branch_op(branch_op& bo, instance& I, stats_t& stat, const solv_options& opts, const uint depth){
    solution_t min_sol;
    int known_solution = I.k + 1;
    // for each branch in the branch list
    for(auto ml = bo.branches.begin(); ml != bo.branches.end(); ++ml){
      // if the branch exceeds the budget (recall that empty branches mean size-1), then don't do it
      if((bo.type != Token) && (bo.type != Deg2Path))
        if((int)ml->size() > min(I.k, known_solution - 1))
          continue;
      DEBUG2(cout << "depth " << depth << " branch: "<<*ml<<endl);
      solution_t solprime;
      // save the first edge of ml in case we need to mark it permanent
      graph_mod_t to_be_permanent(ml->front());
      // make a copy Iprime of I, translating the modlist ml to the new graph
      unordered_map<uint, vertex_p> id_to_vertex;
      instance Iprime(I, &id_to_vertex);
      // convert the mod list for the new graph
      for(auto &gmod : *ml) gmod.e = convert_edge(gmod.e, id_to_vertex);
      // we only need to find solutions that are better than what we have
      Iprime.k = min(I.k, known_solution - 1);
      // delete the edges of this branch
      apply_one_branch(Iprime, bo.type, *ml, solprime);
      // and recurse
      solprime += run_branching_algo(Iprime, stat, opts, depth+1);

      // if we were successfull, save the solution
      if(Iprime.g.vertices.empty() && !(Iprime.k < 0)){
        DEBUG2(cout << solprime << " is indeed a valid solution and its size is " << solprime.size()<<endl);
        min_sol = solprime;
        // we've found a solution that should be smaller than known_solution
        known_solution = solprime.size();
      }
      // mark edges permanent in I (for the next branch)
      // recheck Sud05, but I think we can only mark edgesets of size one permanent!
      if((ml->size() == 1) && (ml->front().type == Del)) {
        DEBUG2(cout << "marking size-1 branch "<<to_be_permanent<<" permanent"<<endl);
        to_be_permanent.e->mark_permanent();
//#error why is this edge not permanent in the next branch???
      }

      DEBUG2(cout << "depth "<<depth<<": current min solution is "<<min_sol << " most recent: "<<solprime << " now searching for solutions of size " << min(I.k, known_solution - 1) <<endl);
    }
    return min_sol;
  }


  // main function solving the problem!
  // TODO: we copy the graph, this is inefficient! improve!
  solution_t run_branching_algo(instance& I, stats_t& stat, const solv_options& opts, const uint depth){
    // keep track of the search tree size
    DO_STAT(stat.searchtree_nodes++);
    DO_STAT(stat.searchtree_depth = max(stat.searchtree_depth, depth));
    DEBUG5(if(stat.searchtree_nodes % 10000 == 0) cout << "currently at "<< stat.searchtree_nodes<<" nodes"<<endl;);
    
    // quick sanity check: if I have less than 7 vertices, then I cannot have a 2-claw, thus the solution is FES
    if(I.g.vertices.size() < 7) return solv_small_instance(I);

    solution_t sol;
    // [1.] apply preprocessing
    DEBUG4(cout << "=== Phase 1 (depth "<<depth<<"): TRRs ===== (k = "<<I.k<<")"<<endl);
    sol += apply_trrs(I, stat);
    DEBUG2(I.g.write_to_stream(std::cout));
    DEBUG4(cout << "=== Phase 2 (depth "<<depth<<"): PRRs ("<< sol.size() <<" dels, k = "<<I.k<<") ====="<<endl);
    // first, apply the split rule
    apply_split_rule(I);
    // the PRRs can give us a list of deg-2 paths which we use for BRR78
    list<path_info_t> deg2paths;
    sol += apply_prrs(I, opts, stat, deg2paths);
    DEBUG2(cout << "welcome back to run_branching_algo (depth "<<depth<<") - sol: "<<sol<<endl);

    // if preprocessing solved I, then return success
    if(I.g.vertices.empty() && !(I.k < 0)) return sol;
   
    // return failure if preprocessing already took all the operations
    // we do another check later, but here, we can avoid computing the lower bound
    if(I.k <= 0) return solution_t();

    // quick sanity check: if I am reduced with respect to the PRRs and TRR and I have less than 8 vertices, then any FES is a solution
    if(I.g.vertices.size() < 8) { sol += solv_small_instance(I); return sol; }

    DEBUG2(cout << "got "<<deg2paths.size()<<" deg2paths:"<<endl; for(auto i = deg2paths.begin(); i != deg2paths.end(); ++i) cout << *i << endl;);
    // [2.] get a lower bound
    // prepare the graph (it shouldn't be required for finding connected components, but we'll need it later anyways)
    // the (better) lower bound by star packing is slower, so just use it every once in a while
    int lower_bnd = compute_lower_bound(I.g, opts, depth);

    DEBUG4(cout << "=== Phase 3 (depth "<<depth<<"): compare budget (" << I.k <<") to lower bound ("<< lower_bnd <<") ====="<<endl);

    // if the lower bound already exceeds our budget then give up
    if(lower_bnd > I.k) {I.k = -1; return solution_t();}

    // [3.] split off connected components if possible
    DEBUG4(cout << "=== Phase 4 (depth "<<depth<<"): global reduction rules & detect connected components ===" << endl);

    // compute the cc number
    I.g.mark_bridges();

    if(I.g.cc_number > 1){
      instance Iprime;
      DEBUG2(cout << "detected " << I.g.cc_number << " components, splitting g"<<endl;);
      
      // split off a connected component from I
      split_off_component(I.g, Iprime.g);
      Iprime.k = I.k;
      
      solution_t rec_sol;
      DEBUG2(cout << "split into components of size "<<I.g.vertices.size()<< " and "<<Iprime.g.vertices.size()<<endl);
      // recurse for both components, summing up solutions (smaller one first)
      if(I.g.vertices.size() < Iprime.g.vertices.size()){
        rec_sol = run_branching_algo(I, stat, opts, depth+1);
        // if there was not enough budget to solve I, return failure
        if(!I.g.vertices.empty() || (I.k < 0)) {I.k = -1; return solution_t();}
        // otherwise, use the remaining budget for the split-off component
        Iprime.k -= rec_sol.size();
        rec_sol += run_branching_algo(Iprime, stat, opts, depth+1);
        // if there was not enough budget to solve this component, then return failure
        if(!Iprime.g.vertices.empty() || (Iprime.k < 0)) {I.k = -1; return solution_t();}
        // otherwise, copy the recursive solution into the overall solution
      } else {
        rec_sol = run_branching_algo(Iprime, stat, opts, depth+1);
        // if there was not enough budget to solve Iprime, return failure
        if(!Iprime.g.vertices.empty() || (Iprime.k < 0)) {I.k = -1; return solution_t();}
        // otherwise, use the remaining budget for the split-off component
        I.k -= rec_sol.size();
        rec_sol += run_branching_algo(I, stat, opts, depth+1);
        // if there was not enough budget to solve this component, return failure
        if(!I.g.vertices.empty() || (I.k < 0)) {I.k = -1; return solution_t();}
        // otherwise, copy the recursive solution into the overall solution
      }
      // if all went well, return success
      sol += rec_sol;
      return sol;
    }

    // try to apply the B-bridge rule, which is possible if and only if G has B-bridgs
    if(opts.use_Bbridge_rule) {
      solution_t bb_sol(apply_Bbridge_rule(I, stat, opts, depth));
      // if we produced edges, recurse (otherwise, our path_infos become invalid)
      if(!bb_sol.empty()){
        DEBUG2(cout << "Bbridge rule got partial solution "<<bb_sol<<", recursing now"<<endl);
        sol += bb_sol;
        sol += run_branching_algo(I, stat, opts, depth);
        return sol;
      }
    }
    // CAUTION: final_RR is incorrect!!!
    //if(apply_final_RR(I, stat, sol)){
    //  sol += run_branching_algo(I, stat, opts, depth);
    //  return sol;
    //}


    // [4.] start branching
    DEBUG4(cout << "=== Phase 5 (depth "<<depth<<"): branchings ("<< sol.size() <<" dels) ====="<<endl);
    DEBUG3(I.g.write_to_stream(std::cout));

    // do the actual branching: first, get a good (the BEST! ^^) branching operation
    branch_op bo;
    if(get_best_branch_op(I.g, bo, deg2paths, !opts.elaborate_branch_selection, opts.keep_searching_if_bnum_above)){

//      if(branch_number(bo) > 2.01){
//        cout << I.g;
//        FAIL("branch number > 2 for this");
//      }

      switch(bo.branches.size()){
        case 0:
          // if there is an empty branch, this means that all
          // branch-edges are permanent and, hence, we have already
          // seen an optimal solution, so no need to continue
          return solution_t();

        case 1:
          // if we have only one branch, then its a reduction rule
          // and we don't need to copy the graph
          DO_STAT(stat.add_BRule(bo));
          
          DEBUG3(cout << "only one branch: "<<bo.branches<<" - wont copy the graph"<<endl);
          apply_one_branch(I, bo.type, bo.branches.front(), sol);
          // and recurse
          sol += run_branching_algo(I, stat, opts, depth+1);
          return sol;

        default:
          DO_STAT(stat.add_BRule(bo));

          solution_t min_sol(apply_branch_op(bo, I, stat, opts, depth));
          DEBUG2(cout << "depth "<<depth<<": done applying branching, we have "; if(min_sol.empty()) cout << "no solution"; else cout << " a solution: "<<min_sol; cout<< endl);
          // if no solution was found, return failure, otherwise merge the minimum solution into sol and clear the graph
          if(min_sol.empty()) return solution_t(); else{
            I.g.clear();
            sol += min_sol;
          }
          // everything went well, so return the solution
          DEBUG3(cout << "returning with "<<stat.searchtree_nodes<<" nodes and solution "<< sol<<endl);
          return sol;
      }
    } else {
      cout << I.g << endl;
      FAIL("no reduction and no branching applies! This shouldn't happen!");
      // no branchable vertices have been found,
      sol += run_branching_algo(I, stat, opts, depth+1);

      DEBUG3(cout << "returning with "<<stat.searchtree_nodes<<" nodes "<<endl);
      return sol;
    }
  }




}; // end namespace


