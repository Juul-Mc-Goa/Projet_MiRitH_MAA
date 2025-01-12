#ifndef KEY_GENERATION_H_
#define KEY_GENERATION_H_

#include "matrix.h"
#include <gmp.h>
#include <stdbool.h>

typedef unsigned int uint;

// A struct holding all the parameters for the signature
typedef struct SignatureParameters {
  uint lambda;
  MatrixSize matrix_dimension;
  mpz_t moduli;
} SignatureParameters;

// A struct holding
// 1. a `seed`, in a array of `lambda` booleans
// 2. a matrix `m0` containing elements of F_q (stored as gnu MP integers)
typedef struct PublicKey {
  uint lambda;
  bool *seed;
  Matrix m0;
} PublicKey;

// A container for a public/private key pair. Returned by the `key_gen`
// function.
typedef struct PublicPrivateKeyPair {
  uint lambda;
  PublicKey public_key;
  bool *private_key;
} PublicPrivateKeyPair;

bool *allocate_seed(uint lambda);
void generate_seed(bool *seed, uint lambda);

void generate_prime(mpz_t result, uint lambda);

PublicPrivateKeyPair allocate_key_pair(SignatureParameters parameters);
void clear_key_pair(PublicPrivateKeyPair key_pair);
PublicPrivateKeyPair key_gen(uint lambda, uint prime_length, MatrixSize size);

#endif // KEY_GENERATION_H_
