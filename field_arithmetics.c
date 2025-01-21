#include "field_arithmetics.h"

typedef unsigned int uint;

uint scalar_add(uint a, uint b) { return a ^ b; }
uint scalar_neg(uint a) { return a; }

uint scalar_double(uint a, FiniteField field) {
  if (a & (1 << (field.log_field_size - 1))) {
    return (a << 1) ^ field.polynomial;
  } else {
    return (a << 1);
  }
}

uint scalar_mul(uint a, uint b, FiniteField field) {
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

FiniteField gf_16() {
  FiniteField field;
  field.field_size = 16;
  field.log_field_size = 4;
  field.polynomial = 19; // X^4 + X + 1

  return field;
}
