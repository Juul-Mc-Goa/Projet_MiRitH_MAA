#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "field_arithmetics.h"
#include "key_generation.h"

extern const char HEX_CHAR_TABLE[16];
extern const FiniteField GF_16;
extern const SignatureParameters PARAMS_1_A_FAST;
extern const SignatureParameters PARAMS_1_A_SHORT;
extern const uint GF_16_ADD_TABLE[16][16];
extern const uint GF_16_MUL_TABLE[16][16];
#endif // CONSTANTS_H_
