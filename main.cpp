
#include "util/graphs.hpp"
#include "util/statistics.hpp"
#include "reduction/trr.hpp"
#include "reduction/prr.hpp"
#include "solv/branching.hpp"
#include "solv/bounds.hpp"
#include "solv/verify.hpp"
#include "solv/solv_opts.hpp"
#include "math.h"

void usage(const char* progname, std::ostream& o){
  o << "usage: " << progname << " file <file to read> [more opts]" << std::endl;
  o << "       " << progname << " rand <vertices> <additional edges> [more opts]"<< std::endl;
  o << "more opts: " << " -lbmod x\t <int>\t apply slower (more powerful) lower bound each x layers (def: "<<cr::default_opts.slow_lower_bound_layers_wait<<")"<< std::endl;
  o << "           " << " -BB x\t {0,1}\t control application of Bbridge branching rule (0=no, 1=yes) (def: "<<cr::default_opts.use_Bbridge_rule <<")"<<std::endl;
  o << "           " << " -YL x\t <int>\t perform Y-lookahead if G has fewer than x vertices (def: "<< cr::default_opts.max_size_for_Y_lookahead<<")"<< std::endl;
  exit(1);
}

void get_random_graph(cr::graph& g, const size_t num_vertices, const size_t num_additional_edges){
  srand(time(NULL));
  cr::vertex_p verts[num_vertices];
  uint j = 0;
  char tmp[20];

  verts[0] = g.add_vertex_fast(std::string("0"));
  for(uint i = 1; i < num_vertices; ++i){
    sprintf(tmp,"%d",i);
    verts[i] = g.add_vertex_fast(std::string(tmp));
    j = round(((double)(i-1) * random()) / RAND_MAX);
    g.add_edge_fast(verts[i], verts[j]);
  }
  for(uint i = 0; i < num_additional_edges; ++i){
    bool success = false;
    do{
      const cr::vertex_p v(verts[(size_t)round(((double)(num_vertices-1) * random()) / RAND_MAX)]);
      success = (g.add_edge_secure(v, verts[(size_t)round(((double)(num_vertices-1) * random()) / RAND_MAX)]) == v->adj_list.end());
    } while(!success);
  }
}

const std::pair<string, int> _requires_params[] = {
  { "file", 1 },
  { "rand",  2 },
  { "-lbmod", 1 },
  { "-BB", 1 },
  { "-YL", 1 }
};
// global arguments with their parameters
std::map<string, std::vector<string> > arguments;

void parse_args(int argc, char** argv, cr::instance& I, cr::solv_options& opts){
  int arg_ptr;
  std::map<std::string, int>  requires_params(std::begin(_requires_params), std::end(_requires_params));
  arg_ptr = 1;
  while(arg_ptr < argc){
    const std::string arg(argv[arg_ptr++]);
    // if the argument is not registered in requires_args, then exit with usage
    if(requires_params.find(arg) == requires_params.end()) usage(argv[0], std::cerr);
    // if there are not enough parameters for this argument
    if(argc < arg_ptr + requires_params[arg]) usage(argv[0], std::cerr);
    // otherwise fill the argument map
    std::vector<string> params(requires_params[arg]);
    for(int i = 0; i < requires_params[arg]; ++i)
      params[i] = argv[arg_ptr++];
    arguments.insert(make_pair(arg, params));
  }
}

int main(int argc, char** argv)
{
  if(argc<2) usage(argv[0], std::cerr);

  cr::instance I;
  cr::solution_t sol;
  cr::solv_options opts(cr::default_opts);

  // parse the arguments, filling 'arguments'
  parse_args(argc, argv, I, opts);

  if(arguments.find("rand") != arguments.end())
    get_random_graph(I.g, stoi(arguments["rand"][0]), stoi(arguments["rand"][1]));
  else if(arguments.find("file") != arguments.end()) I.g.read_from_file(arguments["file"][0].c_str());
  else usage(argv[0], std::cerr);
  if(arguments.find("-lbmod") != arguments.end()) opts.slow_lower_bound_layers_wait = stoi(arguments["-lbmod"][0]);
  if(arguments.find("-BB") != arguments.end()) opts.use_Bbridge_rule = stoi(arguments["-BB"][0]);
  if(arguments.find("-YL") != arguments.end()) opts.max_size_for_Y_lookahead = stoi(arguments["-YL"][0]);
  I.k = INT_MAX;

  //cout <<"Options: " << endl << opts << endl;


  // get an upper bound on k to start with
  cr::solution_t upper_bound(upper_bound_simple(I));
  I.k = upper_bound.size();

  cr::instance Iprime(I);

/*
  std::cout << "prior to  TRRs: " << I.g.vertices.size() << " vertices and " << I.g.edgenum << " edges" << std::endl;
  sol += apply_trrs(I);

  std::cout << "after TRRs: " << std::endl;
  I.g.write_to_stream(std::cout);

  sol += apply_prrs(I);

  std::cout << "after PRRs: " << std::endl;
  I.g.write_to_stream(std::cout);

  std::list<cr::edge_p> bridgelist = I.g.get_bridges();

  uint verts = I.g.vertices.size();
  uint edges = I.g.edgenum;
  uint ccs = I.g.cc_number;
  uint lower_bound = cr::compute_lower_bound(I.g);

  std::cout << "stats: |V|: "<<verts<<" |E|: "<<edges<<" #cc: "<<ccs<<" FES: "<<ccs+edges-verts<<" bridges: "<<bridgelist.size()<<" lowerbound: "<<lower_bound<<std::endl;
*/
  cr::stats_t stats;
  stats.input_FES=cr::get_FES(I.g);
  sol += cr::run_branching_algo(I, stats, opts);

  //std::cout << "verifying size-"<<sol.size()<<" solution " << sol << endl;
  if(! verify_solution(Iprime, sol)) {cout << "======= EPIC FAIL: VERIFICATION FAILED ======" << endl; exit(1);}

  std::cout << "solution: "<< sol << " size: "<<sol.size()<<std::endl;
  std::cerr << stats <<std::endl;
  output_parser_friendly(cout, stats);
}
