#include "graphs.hpp"
#include <unordered_map>
#include <sstream>

namespace cr{

  uint hash_edge::operator()(const edge_p& e) const{
    uint id1(e->head->id);
    uint id2(e->get_tail()->id);
    return id1 * max(id1, id2) + id2;
  }

  // for all j in {0..sizeof(result)-1}:
  //    let V[j] be the set of vertices whose ID is j modulo sizeof(result)
  //    let D[j] be the set of last bits of the degrees of all vertices in V[j]
  //    set j'th most significant bit of result to the parity of D[j]
  uint hash_graph::operator()(const graph& g) const{
    uint result = 0;
    for(auto v = g.vertices.begin(); v != g.vertices.end(); ++v)
      result ^= ( (v->degree() & 1) << (v->id % (8*sizeof(result))) );
  // since the sum over all degrees is always even, the parity of result is always even
  // to compensate for the lost image space, shift it to the left and set the least significant
  // bit of result to the parity of the vertex set
    result = ( result << 1 ) | ( g.vertices.size() & 1 );
    return result;
  }


    
  bool graph::operator==(const graph& g) const{
    // make sure *this and g have the same # of vertices
    if(vertices.size() != g.vertices.size()) return false;
    auto v(vertices.begin());
    auto vprime(g.vertices.begin());
    while(v != vertices.end()){
      if(v->id != vprime->id) return false;
      if(v->degree() != vprime->degree()) return false;
      // test the sets of incident edges by the ids of the heads
      set<uint> el;
      for(auto e = v->adj_list.begin(); e != v->adj_list.end(); ++v)
        el.insert(e->head->id);
      for(auto e = v->adj_list.begin(); e != v->adj_list.end(); ++v)
        if(el.find(e->head->id) == el.end()) return false;

      ++v; ++vprime;
    }
    return true;
  }


  graph::graph(const graph& g, unordered_map<uint, vertex_p>* id_to_vertex):
    // copy graph infos
    current_dfs_id(g.current_dfs_id),
    current_id(0),
    bridges_marked(g.bridges_marked),
    subtrees_marked(false),
    edgenum(0),
    cc_number(g.cc_number)
  {
    DEBUG2(cout << "copy constructing a new graph with "<<g.vertices.size()<<" vertices and "<<g.num_edges()<<" edges"<<endl);
    add_disjointly(g, id_to_vertex);
  }

  // copy constructor - NOTE THAT trr_infos ARE NOT up to date for the copy
  graph::graph(const graph& g, edgelist* const el):
    // copy graph infos
    current_dfs_id(g.current_dfs_id),
    current_id(0),
    bridges_marked(g.bridges_marked),
    subtrees_marked(false),
    edgenum(0),
    cc_number(g.cc_number)
  {
    unordered_map<uint, vertex_p> id_to_vertex;
    add_disjointly(g, &id_to_vertex);
    // if we are also tasked with translating the edgelist el, then do so using id_to_vertex
    if(!el->empty()){
      DEBUG2(cout << "translating edgelist "<< *el << " using "<<id_to_vertex<<endl);
      edgelist::iterator eprime = el->begin();
      bool reached_end;
      do{
        // advance the iterator, since we're going to delete e (and add the corresponding edge to the front of el)
        edgelist::iterator e(eprime++);
        reached_end = (eprime == el->end());

        // push the translated edge to the front of el
        DEBUG2(cout << "translating ("<<id_to_vertex[(*e)->get_tail()->id]<<","<<id_to_vertex[(*e)->head->id]<<")"<<" with adjacency list "<<id_to_vertex[(*e)->head->id]->adj_list<<endl);
        el->push_front(find_edge(id_to_vertex[(*e)->get_tail()->id], id_to_vertex[(*e)->head->id]));
        // and delete e
        el->erase(e);
      } while(!reached_end);
      DEBUG2(cout << "done"<<endl);
    }
  }


    vertex::vertex(const vertex& v): id(v.id), prot(v.prot), dfs_id(v.dfs_id),
    incident_bridges(v.incident_bridges),
    name(v.name),
    parent_valid(false){}


    edgelist vertex::get_cyclic_neighbors(){
      edgelist result;
      // TODO: this is not optimal, maybe maintain a list of cyclic neighbors!
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(!(e->is_bridge)) result.push_back(e);
      return result;
    }

    edgelist vertex::get_cyclic_core_neighbors(){
      edgelist result;
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->head->is_on_cyclic_core()) result.push_back(e);
      return result;
    }

    edge_p vertex::first_cyclic_core_neighbor_except(const vertex_p& except, bool *success){
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->head->is_on_cyclic_core() && (e->head != except)) {
          if(success) *success = true;
          return e;
        }
      if(success) *success = false;

      return adj_list.end();
    }

    edge_p vertex::first_cyclic_core_neighbor(bool *success){
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->head->is_on_cyclic_core()) {
          if(success) *success = true;
          return e;
        }
      if(success) *success = false;

      return adj_list.end();
    }

    edge_p vertex::first_non_bridge_neighbor_except(const vertex_p& except, bool *success){
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(!e->is_bridge && (e->head != except)) {
          if(success) *success = true;
          return e;
        }
      if(success) *success = false;

      return adj_list.end();
    }

    edge_p vertex::first_non_bridge_neighbor(bool *success){
      for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e)
        if(!e->is_bridge) {
          if(success) *success = true;
          return e;
        }
      if(success) *success = false;

      return adj_list.end();
    }

    const edge_p& vertex::get_parent(bool *success) {
      // if we are on the cyclic core, we don't have a parent
      if(!is_on_cyclic_core()){
        DEBUG1(cout << *this << " is not on cylic core since deg = " << degree() << " & subs = " << subtree_NH() << endl);
        // if we already know the parent, just use it
        if(has_parent()) {
          DEBUG1(cout << "I already have a parent" << endl);
          DEBUG1(cout << "I'll just return "<<*parent->head<<endl);
          if(success) *success = true;
        } else {
          // otherwise, find out the parent
          for(edge_p e = adj_list.begin(); e != adj_list.end(); ++e){
            DEBUG1(cout << "is " << e->head << " my parent?"<<endl);
            DEBUG1(cout << e->head<<" infos: subtree: "<<e->head->subtree_NH()<<" deg: "<<e->head->degree()<<")");
            if(e->head->has_parent() ? // if the head has no parent, or a parent that is not us, then this is our parent
              (e->head->get_parent() != e->get_reversed()) : true){
              DEBUG1(cout << " YES!" << endl);
              // if we found a parent, register it for later use
              set_parent(e);
              if(success) *success = true;
              return parent;
            } else DEBUG1(cout << " no..." << endl);
          }
          if(success) *success = false;
        }
      } else {
        DEBUG1(cout << "I'm on the core, no parent for me"<<endl);
        if(success) *success = false;
      }
      return parent;
    }

    uint vertex::nldeg() const{
      uint non_leaves = 0;
      // note that the trr_infos cannot be trusted because this is being called from TRR6, where trr_infos may be broken
      for(list<edge>::const_iterator e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->head->degree() > 1) non_leaves++;
      return non_leaves;
    }


    bool vertex::is_separator() const{
      if(cyc_core_degree() != 2) return false;
      if(trr_infos.ptwos.size() > 1) return true;
      if(!trr_infos.ptwos.empty()) return false;
      if(!trr_infos.ygraphs.empty()) return false;
      if(!trr_infos.leaves.empty()) return true;
      // only degree-two vertices remain from here on

      // find any token-generating neighbors of v
      for(list<edge>::const_iterator e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->head->is_generator() || (e->head->cyc_core_degree() > 2))
          return false;
      // if noone generates a token, then v is a separator
      return true;
    }

    bool vertex::is_incident_to_Bbridge() const{
      for(list<edge>::const_iterator e = adj_list.begin(); e != adj_list.end(); ++e)
        if(e->is_Bbridge()) return true;
      return false;
    }




  // return number of vertices present in the graph
  uint graph::num_vertices() const{
    return vertices.size();
  }

  // return number of edges present in the graph
  uint graph::num_edges() const{
    return edgenum;
  }

  // get an id to perform a dfs
  uint graph::get_dfs_id(){
    ++current_dfs_id;
    // id current_dfs_id got rolled back (carry flag :P) roll back all vertices
    if(!current_dfs_id){
      for(vertex_p v = vertices.begin(); v != vertices.end(); ++v)
        v->dfs_id = 0;
      ++current_dfs_id;
    }
    return current_dfs_id;
  }
  // clear the graph (remove all vertices and edges)
  void graph::clear(){
    vertices.clear();
    bridges_marked = true;
    subtrees_marked = true;
    edgenum = 0;
    cc_number = 0;
  }

  // find a vertex by specifying its id, return vertices.end() if its just not there
  vertex_p graph::find_vertex_by_id(const uint id){
    for(vertex_p i = vertices.begin(); i != vertices.end(); ++i)
      if(i->id == id) return i;
    return vertices.end();
  }

  // add a vertex with a brand new id to the graph and return a fresh iterator to it
  // for id, just take the last vertex in the list and increase his id
  vertex_p graph::add_vertex_fast(){
    return add_vertex_fast(++current_id);
  }

  vertex_p graph::add_vertex_fast(const string& s){
    vertex_p v(add_vertex_fast());
    v->name = s;
    return v;
  }

  // add a vertex given as id to the graph and return a fresh iterator to it
  vertex_p graph::add_vertex_fast(const uint id){
    return vertices.insert(vertices.end(), vertex(id));
  }

  vertex_p graph::add_vertex_fast(const uint id, const string& s){
    vertex_p v(add_vertex_fast(id));
    v->name = s;
    return v;
  }


  // add a vertex given as id to the graph and return a fresh iterator to it
  // if it was there already, return an iterator to this one!
  vertex_p graph::add_vertex_secure(const uint id){
    vertex_p duplicate = find_vertex_by_id(id);
    if(duplicate == vertices.end())
      return add_vertex_fast(id);
    else
      return duplicate;
  }


  // add an edge to the graph - modify adjacency lists
  // this is the _fast_ variant: no check is done whether this edge already exists!
  edge_p graph::add_edge_fast(const vertex_p& u, const vertex_p& w){
    edge eu(u);
    edge ew(w);

    DEBUG2(cout << "adding edge "<<*u<<"-"<<*w<<endl);
    edge_p uadj_pos = u->adj_list.insert(u->adj_list.end(), ew);
    edge_p wadj_pos = w->adj_list.insert(w->adj_list.end(), eu);

    uadj_pos->head_adj_pos = wadj_pos;
    wadj_pos->head_adj_pos = uadj_pos;

    bridges_marked = false;
    subtrees_marked = false;
    edgenum++;

    return uadj_pos;
  }

  edge_p graph::add_edge_fast(const vertex_p& u, const vertex_p& v, const edge_pc& copy_from){
    edge_p result(add_edge_fast(u, v));
    if(copy_from->is_permanent) result->mark_permanent();
    if(copy_from->is_bridge) result->mark_bridge();
    return result;
  }
  // add an edge to the graph - modify adjacency lists
  // this is the _secure_ variant: all sanity checks are performed
  edge_p graph::add_edge_secure(const vertex_p& u, const vertex_p& v){
    // first, disallow loops
    if(u == v) return u->adj_list.end();

    // next, disallow parallel edges
    if(adjacent(u, v)) return u->adj_list.end();

    // if the edge is not yet there, insert it
    return add_edge_fast(u,v);
  }

    // delete a vertex and all incident edge
    void graph::delete_vertex(const vertex_p& v){
      DEBUG1(cout << "deleting "<< *v << endl);
      // delete incident edges
      while(!v->adj_list.empty()) delete_edge(v->adj_list.begin());
      // and remove it from the vertex list
      vertices.erase(v);
    }
    void graph::delete_vertices(list<vertex_p>& vl){
      for(list<vertex_p>::iterator v = vl.begin(); v != vl.end(); ++v)
        delete_vertex(*v);
    }

    // delete an edge
    edge_p graph::delete_edge(const edge_p& e){
      const edge_p& mirror_e(e->head_adj_pos);
      
      DEBUG1(cout << "deleting edge "<< *e << endl);

      vertex_p w(e->head);
      vertex_p u(mirror_e->head);

      if(e->is_bridge){
        u->incident_bridges--;
        w->incident_bridges--;
        cc_number++;
      }

      list<edge>& wadj(w->adj_list);
      list<edge>& uadj(u->adj_list);

      // take care of parent
      if(u->has_parent()) if(u->parent->head == w) u->invalidate_parent();
      if(w->has_parent()) if(w->parent->head == u) w->invalidate_parent();

      // don't forget to update graph variables
      edgenum--;
      bridges_marked = false;
      subtrees_marked = false;

      // perform the delete and return the next edge_p "in line"
      uadj.erase(mirror_e);
      return wadj.erase(e);
    }

    void graph::delete_edges(const edgelist& l){
      for(edgelist::const_iterator e = l.begin(); e != l.end(); ++e)
        delete_edge(*e);
    }

    // delete the connected component containing vertex v
    void graph::delete_component(const vertex_p& v){
      // first, copy the list of neighbors of v
      DEBUG1(cout << "destroying component of " << *v << endl);
      uint dfs_id = get_dfs_id();

      list<vertex_p> to_destroy;
      to_destroy.push_back(v);
      v->dfs_id = dfs_id;

      // destroy the component via BFS
      while(!to_destroy.empty()){
        vertex_p u(to_destroy.front());
        to_destroy.pop_front();

        // add all neighbors...
        for(edge_p e = u->adj_list.begin(); e != u->adj_list.end(); ++e)
          if(e->head->dfs_id != dfs_id){ // ... that are not already in the 'to_destroy' list...
            // put them in the list
            to_destroy.push_back(e->head);
            // and mark them as being DOOOOOOOMED!!
            e->head->dfs_id = dfs_id;
          }
        // and destroy u
        delete_vertex(u);
      }
      cc_number--;
    }


  // simple output,
  // Prints the edgelist of g (plus the number of vertices/edges in verbose mode)
  // g: a graph
  // P: Property map to access the information that should be print for the vertices 
  // out: output stream.
  void graph::write_to_stream(ostream& out, const bool verbose) const{
    
    if(verbose){
      out<<"number of vertices: "<<  num_vertices()  <<endl;
      out<<"number of edges: "<<    num_edges() <<endl;
    }
  
    set<vertex> seen;
    for(list<vertex>::const_iterator v = vertices.begin(); v != vertices.end(); ++v){
      // mark v 'seen' so no edges involving v are printed later
      //std::pair<set<vertex>::iterator, bool> seen_entry = 
        seen.insert(*v);
      for(list<edge>::const_iterator e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
        // if v's current head has not been 'seen' yet, print the edge
        if(seen.find(*(e->head)) == seen.end()){
          out<< *e;
          if(verbose){
            if(e->is_permanent) out << " (P)";
            if(e->is_bridge) out << " (B)";
          }
          out << endl;
        }
    }
  } // end write_graph


  // simple input
  // reads the edgelist into g
  // in: input stream. 
  void graph::read_from_stream(istream& in){
    typedef unordered_map<string, vertex_p>::iterator map_iter;
    
    unordered_map<string, vertex_p> id2vertex;
    string name0,name1;
    pair<map_iter, bool> i,j;

    // clear the graph
    clear();

    while(in){
      // read the two endpoints
      in >> name0;
      in >> name1;

      // get their respective vertices (or create if they don't exist yet)
      i = id2vertex.insert(pair<string, vertex_p>(name0,vertex_p()));
      if(i.second) { // new vertex
        i.first->second = add_vertex_fast();
        i.first->second->name = name0;
      }
      j = id2vertex.insert(make_pair(name1,vertex_p()));
      if(j.second) { // new vertex
        j.first->second = add_vertex_fast();
        j.first->second->name = name1;
      }
      add_edge_secure(i.first->second,j.first->second);
    }
  } // end of read_graph


  // simple input
  // reads the edgelist into f
  // infile: input file namf
  void graph::read_from_file(const char* infile)
  {
    ifstream f(infile);
    read_from_stream(f);
  }


// my own bridge finder, takes linear time and less space per vertex
// this also outputs for each bridge b the number of vertices in the one of the two(!)
// components split off by deleting b
  uint my_bridge_finder(const vertex_p& v, const vertex_p& parent, uint dfs_id, edgelist& bridgelist, list<uint>& comp_sizes){
    v->tarjan_infos.number = parent->tarjan_infos.number + 1;
    v->incident_bridges = 0;
    v->dfs_id = dfs_id;
    uint subtree_size = 1;
    edge_p to_parent;
    uint smallest_num = INT_MAX;

    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e) if(e->head != parent){
      const vertex_p& w(e->head);
      if(w->dfs_id != dfs_id){ // this is a forward edge
        // recurse, saving the number of vertices in the sub-DFS-tree
        subtree_size += my_bridge_finder(w, v, dfs_id, bridgelist, comp_sizes);
        // after recursion, we know whether e is a bridge or not
        if(!e->is_bridge) smallest_num = min(smallest_num, w->tarjan_infos.number);
      } else smallest_num = min(smallest_num, w->tarjan_infos.number);
    } else to_parent = e; // save the edge to the parent in case we want to mark it as bridge
    v->tarjan_infos.number = smallest_num;
    // if parent's number is still smaller than v's, then {v,parent} is indeed a bridge
    if(parent->tarjan_infos.number < smallest_num) {
      DEBUG1(cout << "found bridge: " << v<<"->"<<parent<< " splitting away "<< subtree_size <<" vertices"<< endl);
      to_parent->mark_bridge();
      bridgelist.push_back(to_parent);
      comp_sizes.push_back(subtree_size);
    }
    // in any case, return the number of vertices in this DFS subtree
    return subtree_size;
  }


// Tarjan Bridgefinder follows, lucky for us, C++ doesn't allow declaring functions inside functions so I cannot make it more readable...
// theory largely from wikipedia
// I extended it a bit to also return the size of the component at the head of each bridge (when deleting the bridge)

  // phase 1: give every vertex a number in preorder and compute the tarjan infos
  void tarjan_dfs(vertex_p v, vertex_p parent, uint &current_number, edgelist &bridgelist, list<uint>& comp_sizes){
    // initialize tarjan infos
    v->tarjan_infos.number = current_number;
    v->tarjan_infos.L = current_number;
    v->tarjan_infos.H = current_number;
    v->tarjan_infos.ND = 1;

    // next vertex gets next number
    current_number++;
    
    // go through adjacent vertices and dive whenever we haven't met it before
    for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); e++){
      if(e->head->tarjan_infos.number == 0){ // head is uninitialized
        tarjan_dfs(e->head, v, current_number, bridgelist, comp_sizes);

        const tarjan_info_t& child_infos(e->head->tarjan_infos);
        // update tarjan_infos in v with "tree edges"
        v->tarjan_infos.ND += child_infos.ND;
        v->tarjan_infos.L = min(v->tarjan_infos.L, child_infos.L);
        v->tarjan_infos.H = max(v->tarjan_infos.H, child_infos.H);
        // check whether e is a bridge and mark, don't forget to mark the other direction as well
        if((child_infos.L == child_infos.number) && (child_infos.H < child_infos.number + child_infos.ND)) {
          DEBUG1(cout << "found bridge: " << e<< " splitting away "<< e->head->tarjan_infos.ND <<" vertices"<< endl);
          e->mark_bridge();
          bridgelist.push_back(e);
          comp_sizes.push_back(child_infos.ND);
        }
      } else if(e->head != parent){
        const tarjan_info_t& child_infos(e->head->tarjan_infos);
        // update tarjan_infos in v with "non-tree edges"
        v->tarjan_infos.L = min(v->tarjan_infos.L, child_infos.L);
        v->tarjan_infos.H = max(v->tarjan_infos.H, child_infos.H);
      } // end if
    } // end for
//    DEBUG1(cout << "tarjan infos for " << *v << ": num = " << v->tarjan_infos.number << " L = " << v->tarjan_infos.L << " H = " << v->tarjan_infos.H << " ND = " << v->tarjan_infos.ND << endl);
  }

  // return all bridges in a graph
  void graph::compute_bridges(edgelist& bridgelist, list<uint>& split_off_sizes){
    cc_number = 0;
    if(num_edges() == 0) {bridges_marked = true; return;}

    // prepare tarjan_infos and incident_bridges
    for(vertex_p v = vertices.begin(); v != vertices.end(); ++v) {
      v->incident_bridges = 0;
      v->tarjan_infos = tarjan_info_t();
      v->tarjan_infos.number = 0;
    }

    uint number = 1;
    // compute the connected components
    for(vertex_p v = vertices.begin(); v != vertices.end(); ++v) // need to iterate in case there are multiple components
      if(v->tarjan_infos.number == 0){
        cc_number++;
        tarjan_dfs(v, vertices.end(), number, bridgelist, split_off_sizes);
      }
//    uint dfs_id = get_dfs_id();
//    for(vertex_p v = vertices.begin(); v != vertices.end(); ++v) // need to iterate in case there are multiple components
//      if(v->dfs_id != dfs_id){
//        cc_number++;
//        my_bridge_finder(v, v, dfs_id, bridgelist, split_off_sizes);
//      }

    DEBUG2(cout << "found "<<bridgelist.size()<<" bridges"<<endl);
    bridges_marked = true;
  }

  edgelist graph::get_bridges(){
    edgelist bridges;
    list<uint> split_off_sizes;
    // compute both bridges and split_off_sizes, but ignore the latter
    compute_bridges(bridges, split_off_sizes);
    return bridges;
  }

  weighted_edges graph::get_weighted_bridges(){
    weighted_edges result;
    edgelist bridges;
    list<uint> split_off_sizes;
    // compute both bridges and split_off_sizes and bring them into the right format
    compute_bridges(bridges, split_off_sizes);
    // convert the two lists into a map
    edgelist::iterator e = bridges.begin();
    list<uint>::iterator i = split_off_sizes.begin();
    while(i != split_off_sizes.end()) // both should have the same size
      result.insert(make_pair(*e++, *i++));
    return result;
  }


  // mark all bridges in a graph (also update cyclic core!), UNLESS the bridges are still marked from a previous run
  void graph::mark_bridges(){
    if(!bridges_marked) get_bridges();
  }

  edgelist graph::get_Bbridges(){
    edgelist el(get_bridges());
    // remove non-Bbridges from the bridge-list
    for(edge_pp e = el.begin(); e != el.end();)
      if(!(*e)->is_Bbridge()) el.erase(e++); else ++e;
    return el;
  }

  weighted_edges graph::get_weighted_Bbridges(){
    weighted_edges bridges(get_weighted_bridges());
    // remove non-Bbridges from the bridge-list
    for(weighted_edges::iterator e = bridges.begin(); e != bridges.end();)
      if(!e->first->is_Bbridge()) bridges.erase(e++); else ++e;
    return bridges;
  }


  // get all degree-one vertices
  list<vertex_p> graph::get_leaves(){
    list<vertex_p> tmp;
    for(vertex_p v = vertices.begin(); v != vertices.end(); ++v)
      if(v->degree() == 1) tmp.push_back(v);
    return tmp;
  }

  // copy the connected component containing v from *this to gto
  void graph::copy_component(const vertex_p& v, graph& gto, unordered_map<uint, vertex_p>* id_to_vertex){
    const uint dfs_id(get_dfs_id());
    list<vertex_p> to_consider;

    // prepare to_consider
    to_consider.push_back(v);

    const bool destroy_map(id_to_vertex == NULL);
    if(destroy_map) id_to_vertex = new unordered_map<uint, vertex_p>();

    // and do the BFS
    while(!to_consider.empty()){
      vertex_p u(to_consider.front());
      to_consider.pop_front();

      if(u->dfs_id != dfs_id){
        // copy v into comp, saving the map from v's ID to v's copy
        DEBUG2(cout << "copying "<<*u<<"(id: "<<u->id<<")"<<endl;);
        vertex_p my_new_vertex(gto.add_vertex_fast(u->name));
        id_to_vertex->insert(make_pair(u->id, my_new_vertex));
        u->dfs_id = dfs_id;
  
        // add the neighbors of u to be considered
        for(edge e : u->adj_list)
          if(e.head->dfs_id != dfs_id)
            to_consider.push_back(e.head);
          else // or add the edge if they already have been considered
            gto.add_edge_fast(my_new_vertex, id_to_vertex->at(e.head->id));
      }
    }
    if(destroy_map) delete id_to_vertex;
  }

  // return the size of the connected component of v
  // this is a much lighter version of the copy_component function above
  uint graph::component_size(const vertex_p& v){
    uint count = 0;
    uint dfs_id = get_dfs_id();
    list<vertex_p> to_consider;

    // prepare to_consider
    to_consider.push_back(v);

    // and do the BFS
    while(!to_consider.empty()){
      vertex_p u(to_consider.front());
      to_consider.pop_front();

      if(u->dfs_id != dfs_id){
        count++;
        u->dfs_id = dfs_id;
  
        // add the neighbors of u to be considered
        for(edge_p e = u->adj_list.begin(); e != u->adj_list.end(); ++e)
          if(e->head->dfs_id != dfs_id)
            to_consider.push_back(e->head);
      }
    }
    return count;
  }

  // split a component off of a graph
  void split_off_component(graph& g, graph& comp, unordered_map<uint, vertex_p>* id_to_vertex){
    if(g.cc_number < 2) return;
    if(g.vertices.empty()) return;

    g.copy_component(g.vertices.begin(), comp, id_to_vertex);
    g.delete_component(g.vertices.begin());
  }


  void graph::add_disjointly(const graph& Gfrom, unordered_map<uint, vertex_p>* id_to_vertex){
    if(Gfrom.vertices.empty()) return;
    const bool destroy_map(id_to_vertex == NULL);
    if(destroy_map) id_to_vertex = new unordered_map<uint, vertex_p>();
    // copy vertices
    for(vertex_pc x = Gfrom.vertices.begin(); x != Gfrom.vertices.end(); ++x){
      vertex_p y(add_vertex_fast(x->name));
      y->prot = x->prot;
      (*id_to_vertex)[x->id] = y;
      // copy edges involving v
      for(edge_pc e = x->adj_list.begin(); e != x->adj_list.end(); ++e){
        const unordered_map<uint, vertex_p>::const_iterator i(id_to_vertex->find(e->head->id));
        if(i != id_to_vertex->end()) add_edge_fast(i->second, y, e);
      }
    }
    if(destroy_map) delete id_to_vertex;
  }


  // get an FES via spanning tree
  edgelist get_a_FES(graph& g){
    // make sure we have vertices
    if(g.vertices.empty()) return edgelist();

    edgelist fes;
    const uint dfs_id(g.get_dfs_id());
    unordered_set<vertex_p, vertex_hasher> to_consider;
    to_consider.insert(g.vertices.begin());
    g.vertices.front().dfs_id = dfs_id;

    while(!to_consider.empty()){
      const vertex_p v(*to_consider.begin());
      to_consider.erase(v);

      for(edge_p e = v->adj_list.begin(); e != v->adj_list.end(); ++e)
        if(e->head->dfs_id != dfs_id){
          e->head->dfs_id = dfs_id;
          to_consider.insert(e->head);
        } else if(to_consider.find(e->head) != to_consider.end()) fes.push_back(e);
    }
    DEBUG2(cout << "computed FES "<< fes << endl);
    return fes;
  }

};





