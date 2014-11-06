#include "trr.hpp"

namespace cr{

  // update the subtree neighborhood for each vertex (also mark all type-A bridges (except Y-"heads") permanent)
  solution_t update_TRR_infos(instance& I, stats_t& stats){
    solution_t sol;
    // don't do anything if subtrees are up-to-date
    if(I.g.subtrees_marked) return sol;

    DEBUG2(cout << "updating subtree NH"<<endl);

    list<vertex_p> subtree_done;
    // first, set all subtree_NH's to zero and save all leaves, since their (empty) subtrees are done
    for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v){
      v->invalidate_parent();
      v->trr_infos.clear();
      if(v->degree() == 1) // save all leaves
        if(!v->prot) // except if they are protected
          subtree_done.push_back(v);
    }
    DEBUG2(cout << "found leaves: " << subtree_done << endl);

    // this is more of a leaves-to-root BFS then a dfs, ah well...
    const uint dfs_id = I.g.get_dfs_id();
    // As long as there are unconsidered vertices who have their subtree done, update their subtree NH
    while(!subtree_done.empty()){
      vertex_p v(subtree_done.front());
      subtree_done.pop_front();
    
      DEBUG2(cout << "next in line: " << *v << " with deg = "<< v->degree() << " & subs = " << v->subtree_NH()<<endl);
      // mark v as considered
      v->dfs_id = dfs_id;

      // apply Tree Reductions to the parent
      sol += perform_trrs(I, stats, v);

      // get my parent, we cannot use get_parent() here, since subtrees are not up to date!!
      // but since our subtrees are all done, our parent is the only one whose dfs_id is old
      edge_p to_parent(v->adj_list.begin());
      while(to_parent == v->adj_list.end() ? false : (to_parent->head->dfs_id == dfs_id)) ++to_parent;
      // for the next part, we need the parent
      // (note that, if all neighbors have been considered, then we are the root of a tree and don't need further processing)
      if(to_parent != v->adj_list.end()){
        const vertex_p& parent(to_parent->head);
        DEBUG2(cout << "found parent " << parent << " with deg: "<< parent->degree() << " & subs: " << parent->subtree_NH()<<endl);

        // update trr_infos (including the subtree NH) and the parent of v
        update_trr_infos_from_child(to_parent);

        // mark the parent to be considered if all but one of its neighbors have been deemed in the subtree
        if(! (parent->degree() > parent->subtree_NH() + 1) ){
          DEBUG2(cout << "adding "<<*parent<<" to be considered"<<endl);
          subtree_done.push_back(parent);
        }
      }
    }
    // finally, apply TRRs to the vertices on the cyclic core, since we never did that in the loop
    for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v) 
      if(v->is_on_cyclic_core())
        sol += perform_trrs(I, stats, v);
    

    DEBUG2(cout << "done updating subtree_NH" << endl);
    I.g.subtrees_marked = true;
    return sol;
  }



  // update parents trr_infos from child (e := ( parent <--- child) )
  void update_trr_infos_from_child(const edge_p& e){
    const edge_p& to_child(e->get_reversed());
    const vertex_p& child(to_child->head);
    const vertex_p& parent(e->head);
    const trr_info_t& child_infos(child->trr_infos);

    // always push the newest item to the front
    // otherwise TRR3 might delete a permanent edge after PRR4 cut a de2path
    child->set_parent(e);
    switch(child->degree()){
      case 1:
        // I'm the only neighbor of my child, so add it as leaf
        parent->trr_infos.leaves.push_front(to_child);
        to_child->mark_permanent();
        break;
      case 2:
        // my child has one neighbor aside from me
        if(!child_infos.ygraphs.empty()) // if it's a Y-graph, then we have a 2claw here
          parent->trr_infos.tclaws.push_front(to_child);
        else if(!child_infos.leaves.empty()){ // otherwise, it should be a leaf
          parent->trr_infos.ptwos.push_front(to_child);
          to_child->mark_permanent();
        } else FAIL("epic fail in update_tr_infos_from_child (deg("<<child<<")=2, but no Y or leaf); forgot to apply TRRs?");
        break;
      default:
        if(child_infos.ptwos.size() > 1)
          parent->trr_infos.ygraphs.push_front(to_child);
        else FAIL("epic fail in update_tr_infos_from_child (deg("<<child<<")>2, but no 2 ptwos); forgot to apply TRRs?");
//        else // don't allow Ys with a leaf to register as P2 in the parent
//          if(!child_infos.leaves.empty())
//            parent->trr_infos.ptwos.push_back(to_child); // never use the same vertex for P2 and Y
    }
    DEBUG2(cout << "updated trr_infos of "<<*parent<<" with " << *child<< " ("<<child_infos<<") Now: "<<parent->trr_infos<<" (sub: "<<parent->subtree_NH()<<" deg: "<<parent->degree()<<")"<<endl);
  }


  // apply Tree Reduction Rule X to v assuming all TRRs have been applied to the subtree below v

  // if nldeg(v) < 3   and   there is a ptwo, then delete the leaf of the ptwo
  // then, delete all but one leaf
  bool trr14_subtree(instance& I, stats_t& stats, const vertex_p& v){
    edgelist& leaves(v->trr_infos.leaves);
    edgelist& ptwos(v->trr_infos.ptwos);
    bool result = false;

    // first, if there is a P_2, then
    if((ptwos.size() == 1) && v->trr_infos.ygraphs.empty() && v->trr_infos.tclaws.empty() && !v->is_on_cyclic_core()){
      DO_STAT(stats.reduct_application[TRR4]++);
      DEBUG2(cout << "applying TRR4 to " << *v << " L: "<< leaves << " P2: "<<ptwos<< ")" <<endl);
      // delete its leaf - it's the unique leaf in (ptwos.begin())->head's trr_infos :)
      const vertex_p middle(ptwos.front()->head); // middle vertex is unique
      I.g.delete_vertex(middle->trr_infos.leaves.front()->head);
      middle->trr_infos.leaves.clear();
      // and register the middle vertex as new leaf
      leaves.push_back(ptwos.front());
      ptwos.pop_front();
      result = true;
    }

    // next, delete all but one leaf, unless there's a P2 in which case delete all leaves
    if((leaves.size() > 1) || (!ptwos.empty())){
      DO_STAT(stats.reduct_application[TRR1]++);
      DEBUG2(cout << "TRR1: " << *v << " has " << leaves.size() << " leaves" << endl);
      // delete all leaves but the one with the smallest id
      edge_pp to_del(leaves.begin());
      // if we have no P2, then keep a leaf
      if(ptwos.empty()) ++to_del;
      while(to_del != leaves.end()) {
        const vertex_p u((*to_del)->head);
        to_del = leaves.erase(to_del);
        I.g.delete_vertex(u);
      }
      DEBUG2(cout << "cleared all leaves of " << *v);
      DEBUG2(if(ptwos.empty()) cout << " except " << leaves.front(); cout << endl);
      result = true;
    }
    return result;
  }

  // while |Y|>0 and ( |Y|+|L|+|P|>1 or there is some permanent edge) delete some y of Y & reduce k by one
  solution_t trr2_subtree(instance& I, stats_t& stat, const vertex_p& v){
    const uint leaves_and_ptwos = v->trr_infos.leaves.size() + v->trr_infos.ptwos.size();
    edgelist& ygraphs(v->trr_infos.ygraphs);
    list<vertex_p> to_del;
    solution_t sol;

    // find whether there is a permanent edge around v, if v is on the cyclic core
    bool has_permanent = false;
    if(v->is_on_cyclic_core()){
      const edgelist cyc_nh(v->get_cyclic_core_neighbors());
      for(edge_ppc e = cyc_nh.begin(); e != cyc_nh.end(); ++e)
        if((*e)->is_permanent) { has_permanent = true; break; }
    }

    while(!ygraphs.empty() && ((ygraphs.size() + leaves_and_ptwos > 1) || has_permanent)){
      DO_STAT(stat.reduct_application[TRR2]++);
      DEBUG2(cout << "applying TRR2 to " << *v << " with " << v->subtree_NH() << " of " << v->degree() << " neighbors in subtree" << endl);

      DEBUG3(cout << "TRR2: edge deletion: "<<*(ygraphs.front())<<endl);
      // delete the edge to the y-graph and delete the component containing the ygraph
      to_del.push_back(ygraphs.front()->head);
      I.delete_edge(ygraphs.front(), sol);
      ygraphs.pop_front();
    }
    // delete the y-graphs that were cut off
    for(vertex_p u : to_del) I.g.delete_component(u);
    return sol;
  }

  //  while |P|>2 delete a p in P and decrease k by one; give an edge deletion hint, otherwise, we might
  //  delete the wrong edge after splitting the graph with PRR4
  solution_t trr3_subtree(instance& I, stats_t& stat, const vertex_p& v){
    edgelist& ptwos(v->trr_infos.ptwos);
    list<vertex_p> to_del;
    solution_t sol;

    while(ptwos.size() > 2){
      DO_STAT(stat.reduct_application[TRR3]++);
      DEBUG3(cout << "TRR3: edge deletion: "<<*(ptwos.front())<<endl);
      // remove the connection edge and let delete_component do the rest
      to_del.push_back(ptwos.front()->head);
      I.g.delete_edge(ptwos.front());
      I.k--;
      sol += (string)*v + "->?";
      ptwos.pop_front();
    }
    // delete the edges that were cut off
    for(vertex_p u : to_del) I.g.delete_component(u);
    return sol;
  }



  // if deg(w)=2 and |Y|=1 then delete w, delete y in Y and decrease k by one
  // note: we use the tclaws trr_info at the parent v of w to do this
  // note: please make sure trr2 has been applied to w before (that's what "_subtree" says...)
  solution_t trr5_subtree(instance& I, stats_t& stat, const vertex_p& v){
    if(!v->trr_infos.tclaws.empty()){
      DO_STAT(stat.reduct_application[TRR5]++);
      DEBUG2(cout << "applying TRR5 to " << *v << "(subtree "<<v->subtree_NH() << ", deg "<< v->degree() << ")" << endl);
      const edge_p& tc(v->trr_infos.tclaws.front());
      const vertex_p w(tc->head);
      solution_t sol;

      // delete the edge {v,w} and throw away the component of w
      I.delete_edge(tc, sol);
      I.g.delete_component(w);
      v->trr_infos.tclaws.pop_front();

      // note that some other TRR could now be applicable to v (for instance, TRR4)!
      sol += perform_trrs(I, stat, v);
      return sol;
    }
    return solution_t();
  }


  // actualy do the application on a single vertex whose trr_infos are accurate!
  // return true if something happend and set local_operations to the number of edges deleted
  solution_t perform_trrs(instance& I, stats_t& stats, const vertex_p& v){
    // apply the tree reduction rules!
    // TODO: here, we ignore whether or not trr14 changed the graph since it has no impact outside of the TRRs (I hope)
    // if we are a leaf, then thre is nothing to be done
    if(v->trr_infos.empty()) return solution_t();

    solution_t sol;
    uint old_size;
    do{
      // we remember the size of the solution before applying trrs, so we can perform as long as there are changes
      old_size = sol.size();
      trr14_subtree(I, stats, v);
      sol += trr2_subtree(I, stats, v);
      sol += trr3_subtree(I, stats, v);
      sol += trr5_subtree(I, stats, v);
    } while(old_size != sol.size());
    DEBUG2(cout << "done performing TRRs to " << *v << ", sol = " << sol << " k = "<<I.k<<endl);
    return sol;
  }


 
  // return whether there is a vertex with nldeg>2 or a cycle in the component of v
  bool dfs_discover_nldegthree(const size_t dfs_id, const vertex_p& v, const vertex_p& parent){
    if(v->dfs_id != dfs_id){
      bool akku = false;
      v->dfs_id = dfs_id;
      // be sure to visit the whole connected component, even if a nldeg>2 vertex is found!
      // Otherwise, later searches in this component might miss this vertex
      // trr_infos cannot be trusted since a vertex with subtree_NH == degree() may not have some leaves registered
      if((v->nldeg() > 2) || v->prot){
        DEBUG2(cout << "found that "<<*v<<" has nldeg "<<v->nldeg()<< " (or enjoys devine protection)"<< endl);
        akku = true;
      }
      for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
        if(e->head != parent)
          if(dfs_discover_nldegthree(dfs_id, e->head, v))
            akku = true;// don't break because the cc must be marked
      return akku;
    } else return true;// discovered a cycle, which is equally as good
  }

  
  
  // Tree Reduction Rule 6 removes caterpillars
  bool trr6(instance& I){
    bool result = false;
    size_t dfs_id = I.g.get_dfs_id(); // we're slightly abusing dfs_ids to mark vertices
    list<vertex_p> to_del;

    for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v)
      if((v->degree() <= 1) && (v->dfs_id != dfs_id))
        if(!dfs_discover_nldegthree(dfs_id, v, v)){
          to_del.push_back(v);
          result = true;
        }

    for(list<vertex_p>::iterator v = to_del.begin(); v != to_del.end(); ++v)
      I.g.delete_component(*v);

    DEBUG2(cout << "done with TRR6 - " <<I.g.vertices.size() << " verts and "<<I.g.num_edges()<<" edges left"<< endl);
    return result;
  }



  // applies the tree reduction rules from v upwards after v got thrown out of the cyclic core
  // returns with v being the last found vertex (usually do_not_cross or a cyclic-core vertex)
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const vertexset& do_not_cross){
    solution_t sol;

    DEBUG2(cout << "applying TRRs upwards from " << *v <<endl<< " with: "<<v->trr_infos<<" (sub: "<<v->subtree_NH()<<" deg: "<<v->degree()<<"),"<<endl<<" wont cross "<<do_not_cross<<endl);
    while(true){
      // apply the TRRs
      sol += perform_trrs(I, stats, v);

      // if we meet do_not_cross, then return
      if(do_not_cross.find(v) != do_not_cross.end()) return sol;
      // if v is protected then return
      if(v->prot) return sol;

      // if we're on the cyclic core already, then stop processing
      // if we're at a new root, then stop
      if(v->degree() - v->subtree_NH() != 1) return sol;
   
      // my parent is my only neighbor who does not have me as parent
      DEBUG2(cout << "looking for parent of "<<*v << endl);
      bool success = false;
      const edge_p to_parent(v->get_parent(&success));
      DEBUG2(if(success) cout << "found parent "<<*to_parent->head<< endl; else cout << "no parent found"<<endl;);
      // if none of my neighbors was on the cyclic core, then we're all done
      if(!success) return sol;


      // update the subtree_NH and trr_infos of parent
      update_trr_infos_from_child(to_parent);
      
      // get next vertex
      v = to_parent->head;
    }
  }
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v){
    vertexset l;
    return apply_trrs_upwards_after_cut(I, stats, v, l);
  }

  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const vertex_p& do_not_cross){
    vertexset l({ do_not_cross });
    return apply_trrs_upwards_after_cut(I, stats, v, l);
  }
  
  solution_t apply_trrs_upwards_after_cut(instance& I, stats_t& stats, vertex_p& v, const edgelist& do_not_cross){
    vertexset l;
    for(edge_pc e : do_not_cross) l.insert(e->head);
    return apply_trrs_upwards_after_cut(I, stats, v, l);
  }

  // apply all trrs to all vertices, return number of deletions done
  solution_t apply_trrs(instance& I, stats_t& stats){
    solution_t sol(update_TRR_infos(I, stats));
    // don't forget TRR6
    if(trr6(I)) 
      DO_STAT(stats.reduct_application[TRR6]++);

    return sol;
  }

};
