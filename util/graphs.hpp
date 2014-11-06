#ifndef GRAPHS_HPP
#define GRAPHS_HPP

#include <cstdlib>
#include <set>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>

#include <fstream>
#include <iostream>
#include <sstream>

#include "defs.hpp"


using namespace std;


 
namespace cr {

  // my graph will be a list of adjacency lists
  //  for each vertex, its adjacency list contains pointers to an edge
  //  edges contain pointers to both its endpoints _and_ to the appropriate positions in the respective adjacency lists

  class vertex;
  class edge;
  class graph;
  class instance;


  typedef list<vertex>::iterator vertex_p;
  typedef list<vertex>::const_iterator vertex_pc;
  typedef list<edge>::iterator edge_p;
  typedef list<edge>::const_iterator edge_pc;

  typedef list<vertex_p> vertexlist;
  typedef list<edge_p> edgelist;
  
  typedef vertexlist::iterator vertex_pp;
  typedef vertexlist::const_iterator vertex_ppc;
  typedef edgelist::iterator edge_pp;
  typedef edgelist::const_iterator edge_ppc;

  // hash an edge using the IDs of its vertices
  class hash_edge {
    public:
      uint operator()(const edge_p& e) const;
  };
  typedef unordered_map<edge_p, uint, hash_edge> weighted_edges;
  typedef pair<edge_p, uint> weighted_edge;

  // hash a graph using the ids of its vertices and their degrees
  class hash_graph {
    public:
      uint operator()(const graph& e) const;
  };

    // can't take a list of edges here because some reduction rules replace parts without knowing the exact edge to delete; they will just register "some edge incident to v"
  typedef list<string> solution_t;

  struct tarjan_info_t{
    uint number, L, H, ND;

    tarjan_info_t():number(0),L(0),H(0),ND(0){}
  };
  struct trr_info_t{
    edgelist leaves,ptwos,ygraphs,tclaws;
    // edge FROM child TO parent!

    inline void clear(){leaves.clear();ptwos.clear();ygraphs.clear();tclaws.clear();}
    inline uint count() const{
      return leaves.size() + ptwos.size() + ygraphs.size() + tclaws.size();
    }
    inline bool empty() const{
      return leaves.empty() && ptwos.empty() && ygraphs.empty() && tclaws.empty();
    }

  };
};



namespace cr{
  class vertex {
    friend class edge;
    friend class graph;
  public:
    uint id;

    // vertices can be protected from removal, we need this for the Bbridge_rule
    bool prot;

    // used to detect whether he was visited in the current DFS
    // since c++ stores bools in uint's anyways, we might as well use this to perform less vertex-updates using "dfs id's":
    // each time a dfs is to be perfomed, get a new id from the graph and use this id to check whether we already visited the vertex
    uint dfs_id;

    // number of neighbors that are not in the cyclic core
    uint incident_bridges;
    string name;

    // the adjacency list of the vertex
    list<edge> adj_list;

    // some infos related to TRRs
    trr_info_t trr_infos;

    // finally, tarjan infos
    tarjan_info_t tarjan_infos;
  private:
    bool parent_valid; // DAMN YOU C++ for not having a "NULL"-iterator to represent "uninitialized"
    edge_p parent;

  public:

    /****************
     * constructors
     ****************/

    vertex(const uint new_id):id(new_id),prot(false),dfs_id(0),incident_bridges(0),parent_valid(false){}
    vertex(const vertex& v);
    
    // empty destructor for testing purposes, TODO: remove
    ~vertex(){};

    inline uint subtree_NH() const {return trr_infos.count();}
    inline bool operator<(const vertex& v) const {return id < v.id;}
    inline bool operator==(const vertex& v) const {return id == v.id;}
    inline operator string() const{return name;}

    // This is preferred to is_on_cyclic_core, because if we only have a tree, the root doesn't have a parent!
    inline uint degree() const {
      return adj_list.size();
    }

    void init(const uint _dfs_id, const uint _incident_bridges, const string& _name){
      dfs_id = _dfs_id;
      incident_bridges = _incident_bridges;
      name = _name;
      parent_valid = false;
    }
    void init_from(const vertex_p& v){
      init(v->dfs_id, v->incident_bridges, v->name);
    }
    vertex& operator=(const vertex& v) {
      init(v.dfs_id, v.incident_bridges, v.name);
      return *this;
    }

    inline uint cyc_core_degree() const{
      const uint d = degree() - subtree_NH();
      if(d == 1) return 0; else return d;
    }
    // a vertex is on the cyclic core if its cyc_core_degree() is not zero
#define is_on_cyclic_core() cyc_core_degree()
//    inline bool is_on_cyclic_core() const{
//      return cyc_core_degree();
//    }

    edgelist get_cyclic_neighbors();
    edgelist get_cyclic_core_neighbors();

    // TODO: the next functions are very similar, maybe unify them with a test functional
    // TODO: also all 4 functions are suboptimal, maybe maintain a list of cyclic neighbors!

    // returns edge to the first cyclic core neighbor of v (*success = true) or any vertex (*success = false) if none exists
    edge_p first_cyclic_core_neighbor_except(const vertex_p& except, bool *success = NULL);
    edge_p first_cyclic_core_neighbor(bool *success = NULL);
    
    // returns edge to the first non bridge neighbor of v (*success = true) or any vertex (*success = false) if there is none
    edge_p first_non_bridge_neighbor_except(const vertex_p& except, bool *success = NULL);
    edge_p first_non_bridge_neighbor(bool *success = NULL);

    // get the parent of the vertex, assuming all subtrees have their parent set to us
    inline bool has_parent() const {
      return parent_valid;
    }
    inline void invalidate_parent() {
      parent_valid = false;
    }
   
    const edge_p& get_parent(bool *success = NULL);
    
    inline void set_parent(const edge_p& e){
      parent_valid = true;
      parent = e;
    }



    uint nldeg() const;

    inline uint non_bridge_degree() const{
      return degree() - incident_bridges;
    }
#define is_on_cycle() non_bridge_degree()
//    inline bool is_on_cycle() const{
//      return non_bridge_degree();
//    }


    inline bool pendant_is_single() const{
      return trr_infos.empty();
    }
    inline bool pendant_is_Y() const{
      return !trr_infos.ygraphs.empty();
    }

    inline bool is_generator() const{
      return !trr_infos.ptwos.empty();
    }
    // decide whether v is a token separator
    bool is_separator() const;
    bool is_incident_to_Bbridge() const;
    
    inline bool is_on_backbone() const{
      return !(trr_infos.leaves.empty() && trr_infos.ptwos.empty());
    }
  };

  class vertex_hasher{
  public:
    uint operator()(const vertex_p& x) const{
      return x->id;
    }
  };
  typedef unordered_set<vertex_p, vertex_hasher> vertexset;


  class edge{
    friend class vertex;
    friend class graph;
  public:
    bool is_bridge;
    bool is_permanent;

    const vertex_p head;
    edge_p head_adj_pos;


    /******************
     * constructors
     ******************/
    edge(const vertex_p& _head):is_bridge(false),is_permanent(false),head(_head){}
    // empty destructor for testing purposes, TODO: remove
    ~edge(){};
        
    edge_p& get_reversed(){return head_adj_pos;}
    const edge_p& get_reversed() const{return head_adj_pos;}

    const vertex_p& get_head() const {return head;}
    const vertex_p& get_tail() const {return get_reversed()->head;}


    bool operator==(const cr::edge& e) const{
      return (head == e.head) && (head_adj_pos->head == e.head_adj_pos->head);
    }
    operator string() const {return (string)*(head_adj_pos->head) + "->" + (string)*(head);}
    operator int() const {return head_adj_pos->head->id * 1000 + head->id;}

    void init(const bool _is_permanent, const bool _is_bridge){
      is_permanent = _is_permanent;
      is_bridge = _is_bridge;
      head_adj_pos->is_permanent = _is_permanent;
      head_adj_pos->is_bridge = _is_bridge;
    }
    void init_from(const edge_p& e){
      init(e->is_permanent, e->is_bridge);
    }
    edge& operator=(const edge& e){
      init(e.is_permanent, e.is_bridge);
      return *this;
    }

    inline void mark_permanent(const bool mark = true) {
      is_permanent = mark;
      get_reversed()->is_permanent = mark;
    }
    inline void mark_bridge(const bool mark = true){
      is_bridge = mark;
      head->incident_bridges += (mark ? 1 : -1);

      get_reversed()->is_bridge = mark;
      get_tail()->incident_bridges += (mark ? 1 : -1);
    }

    inline bool is_Bbridge() const{
      if(!is_bridge) return false;
      const vertex_p& v(get_tail());
      return (v->is_on_cyclic_core() && head->is_on_cyclic_core());
    }

    inline bool is_Abridge() const{
      return (is_bridge && !is_Bbridge());
    }

    bool is_relevant_Abridge() const{
      // I am a bridge
      if(!is_bridge) return false;

      // I am not a B bridge
      const vertex_p& v(get_tail());
      if((v->is_on_cyclic_core() && head->is_on_cyclic_core())) return false; // == is_Bbridge(), but no overhead

      // I am relevant
      if(v->cyc_core_degree() == 2) return false;
      if(!v->is_incident_to_Bbridge()) return false;
      return true;
    }


  };

  class edge_hasher{
  public:
    uint operator()(const edge_p& x) const{
      const vertex_hasher h;
      return h(x->head) + h(x->get_tail());
    }
  };
  typedef unordered_set<edge_p, edge_hasher> edgeset;

  class graph {
  private:
    void compute_bridges(edgelist& bridgelist, list<uint>& split_off_sizes);
  public:
    uint current_dfs_id;
    uint current_id;
    // are the bridges up to date?
    bool bridges_marked;
    bool subtrees_marked;
    // number of edges in the graph
    uint edgenum;
    uint cc_number;

    list<vertex> vertices;

    /****************************
     * constructors
     ***************************/
    graph():current_dfs_id(1),current_id(0),bridges_marked(false),subtrees_marked(false),edgenum(0),cc_number(0),vertices(){}
    // initialize while translating the edge list el
    // note that the edgelist may change, but the pointer wont
    graph(const graph& g, edgelist * const el);
    graph(const graph& g, unordered_map<uint, vertex_p>* id_to_vertex = NULL);
    graph(const graph& g, edge_p& e){
      edgelist el; el.push_back(e);
      graph(g, &el);
    }

    /**************************
     * read-only informative functions
     **************************/

    uint num_vertices() const;
    uint num_edges() const;

    // get a new id to perform dfs or whatever, the new id is guaranteed to not occur in the graph (unless the programmer did something stupid)
    // setting dfs_ids to 0 in the graph corresponds to clearing the dfs_id, since this function will never return 0
    uint get_dfs_id();

    // find a vertex by specifying its id, return vertices.end() if its just not there
    vertex_p find_vertex_by_id(const uint id);

    // simple output,
    // Prints the edgelist of g (plus the number of vertices/edges in verbose mode)
    // out: output stream.
    void write_to_stream(ostream& out, const bool verbose = true) const;

    // test whether two graphs are equal (including vertex id's and _the_order_in_the_vertex_list)
    // TODO: the second condition can be dropped if we keep a map:id->vertex in the graph instead of a list<vertex>
    bool operator==(const graph& g) const;

    /************************
     * graph modifications
     ************************/

    // clear the graph (remove all vertices and edges)
    void clear();

    // add a vertex given as id to the graph and return a fresh iterator to it
    vertex_p add_vertex_fast();
    vertex_p add_vertex_fast(const string& s);
    vertex_p add_vertex_fast(const uint id);
    vertex_p add_vertex_fast(const uint id, const string& s);
    vertex_p add_vertex_secure(const uint id);

    // add an edge to the graph - modify adjacency lists
    // this is the _fast_ variant: no check is done whether this edge already exists!
    edge_p add_edge_fast(const vertex_p& u, const vertex_p& v);
    edge_p add_edge_fast(const vertex_p& u, const vertex_p& v, const edge_pc& copy_from);

    // add an edge to the graph - modify adjacency lists
    // this is the _secure_ variant: all sanity checks are performed
    edge_p add_edge_secure(const vertex_p& u, const vertex_p& v);

    // delete a vertex and all incident edge
    void delete_vertex(const vertex_p& v);
    void delete_vertices(list<vertex_p>& vl);

    // delete an edge, return next edge_p in "current" adjacency list
    edge_p delete_edge(const edge_p& e);
    void delete_edges(const edgelist& l);

    // delete the connected component containing vertex v
    void delete_component(const vertex_p& v);

    // delete the connected component containing vertex v
    void copy_component(const vertex_p& v, graph& gto, unordered_map<uint,vertex_p>* id_to_vertex = NULL);
    void copy_graph(const vertex_p& v, graph& gto, unordered_map<uint,vertex_p>* id_to_vertex = NULL);

    // count the number of vertices in the component of v
    uint component_size(const vertex_p& v);

    // simple input
    // reads the edgelist into g
    // in: input stream.
    void read_from_stream(istream& in);

    // simple input
    // reads the edgelist into f
    // infile: input file namf
    void read_from_file(const char* infile);

    // mark all bridges
    void mark_bridges();
    // mark & return all bridges
    edgelist get_bridges();
    // mark & return all bridges, weighted by the size of a cc they split off the graph
    weighted_edges get_weighted_bridges();

    // mark all bridges & return all B-bridges
    edgelist get_Bbridges();
    // mark all bridges & return all B-bridges, weighted by the size of a cc they split off the graph
    weighted_edges get_weighted_Bbridges();

    // add a graph to this one
    void add_disjointly(const graph& Gfrom, unordered_map<uint, vertex_p>* id_to_vertex = NULL);


//    void update_subtree_NH();
    // get all degree-one vertices
    list<vertex_p> get_leaves();
  };




  inline edge_p find_edge(const vertex_p u, const vertex_p v){
    for(edge_p e = u->adj_list.begin(); e != u->adj_list.end(); ++e)
      if(e->head == v) return e;
    return u->adj_list.end();
  }

  inline bool adjacent(const vertex_p u, const vertex_p v){
    return (find_edge(u, v) != u->adj_list.end());
  }

  inline bool contains_Abridge(const edgelist& E){
    for(edgelist::const_iterator i = E.begin(); i != E.end(); ++i)
      if((*i)->is_Abridge()) return true;
    return false;
  }
  inline bool contains_non_relevant_Abridge(const edgelist& E){
    for(edgelist::const_iterator i = E.begin(); i != E.end(); ++i)
      if((*i)->is_Abridge() && !(*i)->is_relevant_Abridge()) return true;
    return false;
  }



  // continue along the given edge on a deg2path
  inline edge_p get_next_on_deg2path(const edge_p& e, bool* success = NULL){
    const vertex_p& v(e->head);
    if(v->cyc_core_degree() != 2){
      FAIL("epic fail in get_next_on_deg2path(" << *e << ")");
    } else
      return v->first_cyclic_core_neighbor_except(e->get_tail(), success);
  }
  // continue along the given edge on a cycle (non_bridges)
  inline edge_p get_next_on_cycle(const edge_p& e){
    return e->head->first_non_bridge_neighbor_except(e->get_tail());
  }

  inline uint get_FES(graph& g){
    g.mark_bridges();
    DEBUG1(cout << "get_FES: g has " << g.edgenum << "(="<<g.num_edges()<<") edges, "<<g.cc_number<<" components, and "<<g.vertices.size()<<" verts"<<endl);
    return g.edgenum + g.cc_number - g.vertices.size();
  }

  // get a feedback edge set via spanning tree
  edgelist get_a_FES(graph& g);

  // split a component off of a graph
  void split_off_component(graph& g, graph& comp, unordered_map<uint, vertex_p>* id_to_vertex = NULL);
  

  inline edge_p first_neighbor_non_dfs_id(const vertex_p& v, const uint dfs_id){
    list<edge>& adj(v->adj_list);
    for(edge_p e = adj.begin(); e != adj.end(); ++e)
      if(e->head->dfs_id != dfs_id) return e;
  }

  // TODO: make this more efficient!!
  inline edge_p convert_edge(const edge_p& e, const unordered_map<uint, vertex_p>& id_to_vertex){
    return find_edge(id_to_vertex.at(e->get_tail()->id), id_to_vertex.at(e->head->id));
  }

  class instance {
  public:
    graph g;
    int k;


    instance():g(),k(0){}
    instance(const graph& _g):g(_g),k(0){}
    instance(const graph& _g, const uint _k):g(_g),k(_k){}
    // instanciate, translating the edgelist el to the new graph
    // note that the edgelist may change, but the pointer wont
    instance(const instance& I, edgelist * const el):g(I.g,el),k(I.k){}
    instance(const instance& I, unordered_map<uint, vertex_p>* id_to_vertex = NULL):g(I.g,id_to_vertex),k(I.k){}
    instance(const instance& I, edge_p& e):g(I.g,e),k(I.k){}

    // delete an edge, registering it in the solution
    edge_p delete_edge(const edge_p& e){
      // decrement k
      k--;
      // and do the actual deleting
      return g.delete_edge(e);
    }
    
    edge_p delete_edge(const edge_p& e, solution_t& solution){
      // register e in the solution
      solution += (string)*e;

      // do the edge deletion
      return delete_edge(e);
    }

    void delete_edges(const edgelist& l){
      for(edgelist::const_iterator e = l.begin(); e != l.end(); ++e)
        delete_edge(*e);
    }
    void delete_edges(const edgelist& l, solution_t& solution){
      for(edgelist::const_iterator e = l.begin(); e != l.end(); ++e)
        delete_edge(*e, solution);
    }

  };

  // convert an edgelist to a solution_t
  inline solution_t edgelist_to_solution(const edgelist& el){
    solution_t result;
    for(edge_ppc e = el.begin(); e != el.end(); ++e) result += (string)**e;
    return result;
  }




};

//ostream& operator<<(ostream& os, const cr::solution_t& s) {return os << s.i; }

inline ostream& operator<<(ostream& os, const cr::vertex& v) {return os << (string)v; }
inline ostream& operator<<(ostream& os, const cr::vertex_p& v) {return os << (const cr::vertex&)(*v); }
inline ostream& operator<<(ostream& os, const cr::vertex_pc& v) {return os << (const cr::vertex&)(*v); }

inline ostream& operator<<(ostream& os, const cr::edge& e) {return os << (string)e;}
inline ostream& operator<<(ostream& os, const cr::edge_p& e) {return os << (cr::edge)*e;}
inline ostream& operator<<(ostream& os, const cr::edge_pc& e) {return os << (cr::edge)*e;}

inline ostream& operator<<(ostream& os, const cr::trr_info_t& t){
  return os << "L="<<t.leaves<<" P="<<t.ptwos<<" Y="<<t.ygraphs<<" 2C="<<t.tclaws;
}
inline ostream& operator<<(ostream& os, const cr::graph& g){
  g.write_to_stream(os);
  return os;
}




#endif
