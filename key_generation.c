#include "gmp.h"
#include <stdbool.h>
#include <stdlib.h>

typedef unsigned int uint;

// A struct holding
// 1. a `seed`, in a array of `lambda` booleans
// 2. a matrix `m0` of dimension `number_of_rows` * `number_of_columns`,
//    containing elements of F_q (stored as gnu MP integers)
typedef struct PublicKey {
  uint lambda;
  bool *seed;
  uint number_of_rows;
  uint number_of_columns;
  mpz_t **m0;
} PublicKey;

// A container for a public/private key pair. Returned by the `key_gen`
// function.
typedef struct PublicPrivateKeyPair {
  uint lambda;
  PublicKey public_key;
  bool *private_key;
} PublicPrivateKeyPair;

bool *allocate_seed(uint lambda) {
  bool *result = malloc(lambda * sizeof(bool));
  return result;
}

void generate_seed(bool *seed, uint lambda) {
  // todo
}

PublicPrivateKeyPair key_gen(uint lambda) {
  // todo
}
