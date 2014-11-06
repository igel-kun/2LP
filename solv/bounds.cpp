#include "bounds.hpp"
#include "../reduction/trr.hpp"
#include <unordered_map>
#include <unordered_set>

namespace cr{

  inline bool edge_p_equal(const edge_p& a, const edge_p& b){
    return (a == b) || (a->get_reversed() == b);
  }
  bool edge_p_cmp(const edge_p& a, const edge_p& b){

    edge_p e1 = (a->head->id > a->get_reversed()->head->id?a:a->get_reversed());
    edge_p e2 = (b->head->id > b->get_reversed()->head->id?b:b->get_reversed());
    
    const vertex_p& v(e1->get_reversed()->head);

    if(v->id == e2->get_reversed()->head->id){
      return distance(e1, v->adj_list.end()) > distance(e2, v->adj_list.end());
    } else return v->id < e2->get_reversed()->head->id;
    
  }

  bool degree_is_less(const vertex& v1, const vertex& v2){
    return v1.degree() < v2.degree();
  }
  // compute a lower bound using a packing of 2-stars
  uint star_packing(graph g){
    int k = 0;
    // first, copy the graph and mark all bridges
    // second, sort the vertices by their degree
    g.vertices.sort(degree_is_less);
    // then check each vertex for being the center of a 2-star
    for(vertex_p v = g.vertices.begin(); v != g.vertices.end();++v){
      // save the edges we want to delete, we have to be able to check whether we already marked one for deletion
      edgeset to_delete;
      // mark this vertex and all its neighbors as in the 2-star
      // for each neighbor, try to find another edge
      for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e){
        // get the first neighbor of e->head
        edge_p a(e->head->adj_list.begin());
        // while a either leads to v or is already in to_delete, advance a
        do{
          if(a->head != v)
            if(to_delete.find(a) == to_delete.end())
              break;
          ++a;
        } while(a != e->head->adj_list.end());
        // if we found someone with the properties above, then mark the two edges of the branch for deletion
        if(a != e->head->adj_list.end()){
          to_delete.insert(e);
          // put a in reversed, so we can check for the unreversed edge from the other side later on
          to_delete.insert(a->get_reversed());
          DEBUG2(cout << "inserting "<<e <<" and "<<a->get_reversed()<<endl);
        }
      }
      // to_delete.size() / 2  is the number of rays in this 2-star
      uint rays = to_delete.size() >> 1;
      if(rays > 2 ){
        // all but two rays have to be destroyed
        k += rays - 2;
        DEBUG2(cout << "found a star with "<<rays<<" rays (now, k = "<<k<<"):" <<endl);
        DEBUG2(for(edgeset::iterator i = to_delete.begin(); i != to_delete.end(); ++i) cout << *i << endl);

        for(edgeset::iterator i = to_delete.begin(); i != to_delete.end(); ++i)
          g.delete_edge(*i);
      }
    }
    // finally, add an FES of the remaining graph
    k += get_FES(g);
    return k;
  }

  // compute lower bound for k
  uint compute_lower_bound(graph& g, const bool more_elaborate){
    uint FES_lb = get_FES(g);
    uint STAR_lb = (more_elaborate ? star_packing(g) : 0);
    DEBUG2(cout << "lower bounds: FES: " << FES_lb << " STAR: "<<STAR_lb<<endl);
    return max(FES_lb, STAR_lb);
//    if(more_elaborate) return STAR_lb; else return FES_lb;
  }

  uint compute_lower_bound(graph& g, const solv_options& opts, const uint depth){
    uint lower_bound = 0;
    // accumulate lower bounds if depth devides by their respective layer_wait's
    if((depth % opts.fast_lower_bound_layers_wait) == 0) lower_bound = max(lower_bound, get_FES(g));
    if((depth % opts.slow_lower_bound_layers_wait) == 0) lower_bound = max(lower_bound, star_packing(g));
    return lower_bound;
  }


  // make the given vertex nldeg2 by deleting more or less random edges
  solution_t make_nldeg2(instance& I, const vertex_p& v){
    solution_t sol;
    edge_p e(v->adj_list.begin());
    uint nldeg = v->nldeg();

    while((nldeg > 2) && (e != v->adj_list.end())) 
      if(!e->is_Abridge() && !e->is_permanent){ // don't delete A-bridges or permanent edges
        const vertex_p& w(e->head);
        edge_p to_del;

        if(w->degree() == 2){ // if w is deg-2, then return the edge that's not {v,w}, unless this one is an Abridge or permanent
          to_del = w->adj_list.begin();
          if(to_del->head == v){
            ++to_del;
            if(to_del->is_permanent || to_del->is_Abridge()) --to_del;
          }
        } else // otherwise, return {v,w}
          to_del = e;

        // increase e before it might get invalidated by deleting to_del
        ++e;
        // delete to_del
        I.delete_edge(to_del, sol);
        // if we cannot do better than the budget, then return the empty solution
        if(I.k <= 0) return solution_t();
        // no matter what we did, the nldeg decreased
        --nldeg;
      } else ++e;
    return sol;
  }

  // find a 2-claw and destroy it by deleting random edges
  solution_t upper_bound_simple(const instance& _I){
    solution_t sol;
    instance I(_I);
    // mark bridges, since we'll want to use is_Abridge()
    I.g.mark_bridges();

    // make each vertex nldeg2
    for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v)
      sol += make_nldeg2(I, v);

    // if the budget is gone, just return the empty solution
    if(I.k <= 0) return solution_t();
    // finally, get rid of all 2-claws that are surrounded by A-bridges
    stats_t tmp;
    sol += apply_trrs(I, tmp);
    // add an FES of the rest of the graph, unless it's over the budget
    const uint FES(get_FES(I.g));
    DEBUG2(cout << "upper bound computed by heuristic is "<<sol.size()<<"+"<<FES<<", while k = "<<I.k<<endl);
    if(I.k <= (int)FES) return solution_t();

    for(uint i = 0; i < FES; ++i)
      sol += "[a non-bridge]";
    return sol;
  }

};
