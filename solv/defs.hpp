#ifndef SOLV_DEFS_HPP
#define SOLV_DEFS_HPP

#include "../util/graphs.hpp"

namespace cr{
  enum mod_type {Del, Yify}; // either delete the edge or Y-graphify it
  enum branch_type {Triangle, Claw0, Claw1, Claw2, Claw3, Deg2Path, Token, Bbridge};


  struct claw_leg {
    edge_p head;
    edgelist E;
  };

  struct graph_mod_t {
    mod_type type;
    edge_p e;

    graph_mod_t(const edge_p& _e, const mod_type& _type = Del): type(_type), e(_e){}
  };
  typedef list<graph_mod_t> modlist_t;

  struct branch_op {
    branch_type type;
    // the actual branches
    list<modlist_t> branches;
    float bnum;

    // constructors
    branch_op(const branch_type t):type(t),branches(),bnum(0){};
    branch_op() {};
  };
  typedef list<branch_op> branchlist;
}

inline ostream& operator<<(ostream& os, const cr::branch_type& t){
  switch(t){
    case cr::Triangle: os <<"Triangle"; break;
    case cr::Claw0:    os <<"Claw0"; break;
    case cr::Claw1:    os <<"Claw1"; break;
    case cr::Claw2:    os <<"Claw2"; break;
    case cr::Claw3:    os <<"Claw3"; break;
    case cr::Deg2Path: os <<"Deg2Path"; break;
    case cr::Token:    os <<"Token"; break;
    case cr::Bbridge:  os <<"B-bridge"; break;
  }
  return os;
}

inline ostream& operator<<(ostream& os, const cr::mod_type& t){
  switch(t){
    case cr::Del:    os <<"del"; break;
    case cr::Yify:  os <<"Yy"; break;
  }
  return os;
}

inline ostream& operator<<(ostream& os, const cr::graph_mod_t& m){
  return os << m.type << ' ' << m.e;
}


#endif
