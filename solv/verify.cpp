#include "verify.hpp"
#include "../reduction/trr.hpp"
#include "../util/statistics.hpp"
#include "../solv/branching.hpp"
#include <map>
#include <string>
#include <algorithm>

namespace cr{

  vector<string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
      std::stringstream ss(s);
      std::string item;
      while (std::getline(ss, item, delim)) {
          elems.push_back(item);
      }
      return elems;
  }


  bool verify_solution(instance I, solution_t sol){
    map<string, vertex_p> name_to_vertex;
    DEBUG3(uint sol_size = sol.size());
    DEBUG2(cout << "======== final phase: verification of "<<sol<<endl);
    // start with applying TRR6 to get rid of caterpillars
    
    // first, get the ID to vertex map
    for(vertex_p v = I.g.vertices.begin(); v != I.g.vertices.end(); ++v)
      name_to_vertex.insert(make_pair(v->name, v));

    for(solution_t::iterator s = sol.begin(); s != sol.end();){
      uint delim_pos = s->find("->");
      if(delim_pos != s->npos){
        // parse the two vertex IDs
        string name1 = s->substr(0, delim_pos);
        string name2 = s->substr(delim_pos+2);
        // remove primes
        name1.erase(remove(name1.begin(), name1.end(), '\''), name1.end());
        name2.erase(remove(name2.begin(), name2.end(), '\''), name2.end());
        // delete edge
        if(name_to_vertex.find(name1) != name_to_vertex.end() && name_to_vertex.find(name2) != name_to_vertex.end()){
          edge_p to_delete = find_edge(name_to_vertex[name1], name_to_vertex[name2]);
          if(to_delete != name_to_vertex[name1]->adj_list.end()){
            I.g.delete_edge(to_delete);
            // remove from solution
            sol.erase(s++);
            // remove caterpillars from the remaining instance, this way, we also notice non-optimality
          } else {DEBUG2(cout << "couldn't find "<<*s<<" in the graph"<<endl); ++s;}
        } else {DEBUG2(cout << *s << " is a special edge, I cannot delete it"<<endl); ++s;}
      } else {DEBUG2(cout << *s << " is a special edge, I cannot delete it"<<endl); ++s;}
    }
    DEBUG3(cout << " successfully deleted "<< sol_size - sol.size()<< " edges, "<< sol.size() << " to go"<<endl);
    trr6(I);

    I.k = sol.size();
    // use the branching algorithm to solve the remaining graph
    stats_t stats;
    solution_t new_sol(run_branching_algo(I, stats));
    // if I got solved within the budget, then the solution is verified, otherwise, it's not
    if( (new_sol.size() == sol.size()) && I.g.vertices.empty()){
      DEBUG4(cout << "verified " << sol << " by " << new_sol<<endl);
      return true;
    } else {
      DEBUG4(cout << "failed to verify! verified " << new_sol<< " (size "<<new_sol.size()<<") but tested "<<sol<<" (size "<<sol.size()<<")"<<endl);
      DEBUG4(if(!I.g.vertices.empty()) cout << "remaining graph: "<<endl<<I.g<<endl);
      return false;
    }
  }

}
