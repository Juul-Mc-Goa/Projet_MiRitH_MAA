#include "field_arithmetics.h"
#include "constants.h"
#include <stdio.h>

typedef unsigned int uint;

uint scalar_add(uint a, uint b) { return GF_16_ADD_TABLE[a][b]; }
uint scalar_neg(uint a) { return a; }
uint scalar_mul(uint a, uint b) { return GF_16_MUL_TABLE[a][b]; }

// used for the `compute_mul` function
uint compute_double(uint a) {
  if (a & (1 << 3)) {
    return (a << 1) ^ GF_16.polynomial;
  } else {
    return (a << 1);
  }
}

// used for the `print_gf_16_mul_table` function
uint compute_mul(uint a, uint b) {
  uint result = 0;
  uint shifted_a = a;

  for (uint i = 0; i < 4; i++) {
    if (b & (1 << i)) {
      result ^= shifted_a;
    }
    shifted_a = compute_double(shifted_a);
  }

  return result;
}

// Generate addition table in GF(16). The output is then used in the
// `GF_16_ADD_TABLE` constant.
void print_gf_16_addition_table() {
  for (uint i = 0; i < 16; i++) {
    printf("{");
    for (uint j = 0; j < 15; j++) {
      printf("%u, ", i ^ j);
    }
    printf("%u},\n", i ^ 15);
  }
}

// Generate multiplication table in GF(16). The output is then used in the
// `GF_16_MUL_TABLE` constant.
void print_gf_16_mul_table() {
  for (uint i = 0; i < 16; i++) {
    printf("{");
    for (uint j = 0; j < 15; j++) {
      printf("%u, ", compute_mul(i, j));
    }
    printf("%u},\n", compute_mul(i, 15));
  }
}
