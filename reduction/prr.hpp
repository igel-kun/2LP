#ifndef PRR_H
#define PRR_H

#include "../util/defs.hpp"
#include "../util/graphs.hpp"
#include "../util/statistics.hpp"
#include "../solv/solv_opts.hpp"
#include "defs.hpp"

namespace cr{

  // add a leaf to v, keeping trr_infos up to date; returns edge to the new leaf
  edge_p& copy_leaf(graph& g, const vertex_p& v, const vertex_p& leaf);
  // add a P2 to v keeping trr_infos up to date; return edge to the center vertex
  edge_p& copy_P2(graph& g, const vertex_p& v, const vertex_p& center);
  // add a Y graph to v keeping trr_infos up to date; return edge to the center vertex
  edge_p& copy_Y(graph& g, const vertex_p& v, const vertex_p& center);


  struct path_info_t{
    // start and end of the path
    edge_p start, end;
    // list of generators: edge points TO the generator
    list<edge_p> generators;
    list<edge_p> end_generators;
    // list of vertices with pendant Y-graphs on this path
    list<vertex_p> pendantYs;

    // do we have a separator
    vertexset separators;

    // length is the number of inner edges, that is the number of inner vertices + 1
    uint length;
    // used to invalidate the path, for example when destroying it
    bool valid;

    path_info_t():length(1),valid(false){}
    bool is_loop(){return start->head == end->get_reversed()->head;}
  };




  // apply the path reduction rules to the instance I, producing a list of deg2paths (represented by their infos)
  solution_t apply_prrs(instance& I, const solv_options& opts, stats_t& stat, list<path_info_t>& infos);

  bool apply_prrs_to_vertex(instance& I, stats_t& stat, solution_t& sol, list<path_info_t>& infos, const vertex_p& v, const uint dfs_id);
}


inline ostream& operator<<(ostream& os, const cr::path_info_t& t){
  if(!t.valid)
    return os << "INVALID ";
  else
    return os << "path: "<< *t.start << ".." << *t.end << " (len " << t.length << ") gen: " << t.generators << " sep: " << t.separators<< " pY: " << t.pendantYs;
}

#endif
