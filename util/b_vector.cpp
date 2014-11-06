#include "b_vector.hpp"
#include "../solv/defs.hpp"
#include <cmath>

#define PRECISION 5


// branching number code by Joseph, Chuang-Chieh Lin (lincc@cs.ccu.edu.tw or josephcclin@gmail.com)
inline float CH_POLY(const uint *A, const float var, const uint n){
	float s = 0;
 	for(uint i=0; i != n; i++) s += pow(var,*(A+i));
  return 1-s; 
}


float branch_number(const uint *BV, const uint n){
  float temp = 0;
  float poly_result = 1;
  for(uint d=1; d<=PRECISION; d++) { 
    poly_result = 1;
    while (poly_result > 0) {
      temp += pow(0.1, d);
      poly_result = CH_POLY(BV, temp, n);
    } 
    temp -= pow(0.1, d);
  }
  return (1/temp);
}

float branch_number(const cr::branch_op& bop){
  const std::list<cr::modlist_t>& branches(bop.branches);

  if(branches.empty()) return FLT_MAX;
  if(branches.size() == 1) return 1/(float)branches.size();

  uint *b_vector = new uint[branches.size()];
  uint index = 0;
  // convert the bop.branches into a vector of sizes
  // ATTENTION: note that we use empty branches to represent some branchings (BRR78), so count empty branches as size-1
  for(auto i = branches.begin(); i != branches.end(); ++i){
    b_vector[index] = 0;
    for(auto b = i->begin(); b != i->end(); ++b)
      if((b->type == cr::Del) || (b->type == cr::Yify))
        b_vector[index]++;
    ++index;
  }

  float result = branch_number(b_vector, index);
  delete[] b_vector;
  return result;
}


