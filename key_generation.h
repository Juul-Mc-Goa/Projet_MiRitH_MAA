#ifndef KEY_GENERATION_H_
#define KEY_GENERATION_H_

#include "matrix.h"
#include <gmp.h>
#include <stdbool.h>

typedef unsigned int uint;

// A struct holding all the parameters for the signature
typedef struct SignatureParameters {
  uint lambda;                 // security parameter
  MatrixSize matrix_dimension; // dimension (m, n) of the full matrices
  uint prime_bit_count; // the size of the prime defining the finite field
  uint target_rank;     // rank `r` of the solution
  uint solution_size;   // size `k` of the solution vector
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
void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            uint bit_count);

PublicPrivateKeyPair allocate_key_pair(SignatureParameters parameters);
void clear_key_pair(PublicPrivateKeyPair key_pair);
PublicPrivateKeyPair key_gen(SignatureParameters params);

#endif // KEY_GENERATION_H_
