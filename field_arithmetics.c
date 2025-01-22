#include "field_arithmetics.h"
#include "constants.h"
#include <stdio.h>

typedef unsigned int uint;

uint scalar_add(uint a, uint b) { return GF_16_ADD_TABLE[a][b]; }
uint scalar_neg(uint a) { return a; }
uint scalar_mul(uint a, uint b, FiniteField field) {
  return GF_16_MUL_TABLE[a][b];
}

// used for the `compute_mul` function
uint compute_double(uint a, FiniteField field) {
  if (a & (1 << (field.log_field_size - 1))) {
    return (a << 1) ^ field.polynomial;
  } else {
    return (a << 1);
  }
}

// used for the `print_gf_16_mul_table` function
uint compute_mul(uint a, uint b, FiniteField field) {
  uint result = 0;
  uint shifted_a = a;

  for (uint i = 0; i < field.log_field_size; i++) {
    if (b & (1 << i)) {
      result ^= shifted_a;
    }
    shifted_a = scalar_double(shifted_a, field);
  }

  return result;
}

void print_gf_16_addition_table() {
  for (uint i = 0; i < 16; i++) {
    for (uint j = 0; j < 16; j++) {
      printf("%u, ", i ^ j);
    }
    printf("\n");
  }
}

void print_gf_16_mul_table() {
  for (uint i = 0; i < 16; i++) {
    for (uint j = 0; j < 16; j++) {
      printf("%u, ", scalar_mul(i, j, GF_16));
    }
    printf("\n");
  }
}
