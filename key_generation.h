#ifndef KEY_GENERATION_H_
#define KEY_GENERATION_H_

#include "field_arithmetics.h"
#include "matrix.h"
#include <gmp.h>
#include <stdbool.h>

typedef unsigned int uint;

// A struct holding all the parameters for the signature
typedef struct SignatureParameters {
  uint lambda;                 // security parameter
  MatrixSize matrix_dimension; // dimension (m, n) of the full matrices
  FiniteField field;           // the field used for computations
  uint target_rank;            // the rank `r` of the solution
  uint solution_size;          // the size `k` of the solution vector
  uint first_challenge_size; // the number of columns `s` in the first challenge
  uint number_of_parties;    // the number of parties `N`
  uint tau;                  // the number of rounds
} SignatureParameters;

// A struct holding
// 1. a `seed`, in a array of `lambda` booleans
// 2. a matrix `m0` containing elements of F_q (stored as gnu MP integers)
typedef struct PublicKey {
  uint lambda;
  bool *seed;
  Matrix m0;
} PublicKey;

// A struct holding a `seed` in an array of `lambda` booleans.
typedef struct PrivateKey {
  uint lambda;
  bool *seed;
} PrivateKey;

// A container for a public/private key pair. Returned by the `key_gen`
// function.
typedef struct PublicPrivateKeyPair {
  PublicKey public_key;
  PrivateKey private_key;
} PublicPrivateKeyPair;

bool *allocate_seed(uint lambda);
void seed_random_state(bool *seed, uint lambda, gmp_randstate_t random_state);

void allocate_key_pair(PublicPrivateKeyPair *key_pair,
                       SignatureParameters parameters);
void clear_key_pair(PublicPrivateKeyPair key_pair);
PublicPrivateKeyPair key_gen(SignatureParameters params);

#endif // KEY_GENERATION_H_
