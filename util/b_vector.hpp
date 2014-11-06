#ifndef B_VECTOR_HPP
#define B_VECTOR_HPP

#include "../util/defs.hpp"
#include "../solv/defs.hpp"

float branch_number(const uint *BV, const uint n);
float branch_number(const cr::branch_op& bop);

#endif
