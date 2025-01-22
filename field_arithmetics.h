#ifndef FIELD_ARITHMETICS_H_
#define FIELD_ARITHMETICS_H_

typedef unsigned int uint;

typedef struct FiniteField {
  uint field_size;
  uint log_field_size;
  uint polynomial;
} FiniteField;

uint scalar_add(uint a, uint b);
uint scalar_neg(uint a);
uint scalar_double(uint a, FiniteField field);
uint scalar_mul(uint a, uint b, FiniteField field);

void print_gf_16_addition_table();
void print_gf_16_mul_table();
#endif // FIELD_ARITHMETICS_H_
