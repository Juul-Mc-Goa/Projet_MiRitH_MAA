#ifndef FIELD_ARITHMETICS_H_
#define FIELD_ARITHMETICS_H_

#include "types.h"
typedef unsigned int uint;

uint scalar_add(uint a, uint b);
uint scalar_neg(uint a);
uint scalar_mul(uint a, uint b);

uint compute_double(uint a);
uint compute_mul(uint a, uint b);
void print_gf_16_addition_table();
void print_gf_16_mul_table();
#endif // FIELD_ARITHMETICS_H_
