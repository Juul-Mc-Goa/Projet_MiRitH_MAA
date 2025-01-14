#include "key_generation.h"
#include "matrix.h"
#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

typedef unsigned int uint;

bool *allocate_seed(uint lambda) {
  bool *result = malloc(lambda * sizeof(bool));
  return result;
}

/* generate the seed, using the `getrandom` system call (available on Linux) */
void generate_seed(bool *seed, uint lambda) {
  uint byte_size = 8;
  uint size = lambda / byte_size;
  if ((lambda % byte_size) != 0) {
    size += 1;
  }

  // the getrandom function outputs a random array of bytes, that we package
  // in a char* variable. then we unpack it into a boolean array.
  char *raw_seed = malloc(size);
  ssize_t result_size = getrandom(raw_seed, size, 0);

  if (result_size != size) {
    printf("did not fill the required number of random booleans.\n");
    printf("requested: %u bytes, filled: %ld bytes\n", size, result_size);
    return;
  }

  for (uint i = 0; i < lambda; i++) {
    uint raw_index = i / byte_size;
    uint bit_index = i % byte_size;
    uint random_val = raw_seed[raw_index];

    seed[i] = (random_val & (1 << bit_index)) != 0 ? true : false;
  }

  free(raw_seed);
  raw_seed = NULL;
}

void seed_as_string(char *result, uint lambda, bool *seed) {
  for (uint i = 0; i < lambda; i++) {
    result[i] = seed[i] ? '1' : '0';
  }
}

void generate_prime(mpz_t result, uint bit_length) {
  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  mpz_urandomb(result, random_state, bit_length);
  mpz_nextprime(result, result);

  gmp_randclear(random_state);
}

PublicPrivateKeyPair allocate_key_pair(SignatureParameters parameters) {
  PublicPrivateKeyPair result;

  result.lambda = parameters.lambda;
  result.private_key = malloc(parameters.lambda * sizeof(bool));

  result.public_key.lambda = parameters.lambda;
  result.public_key.seed = malloc(parameters.lambda * sizeof(bool));
  allocate_matrix(&result.public_key.m0, parameters.matrix_dimension);

  return result;
}

void clear_key_pair(PublicPrivateKeyPair key_pair) {
  free(key_pair.private_key);
  free(key_pair.public_key.seed);
  clear_matrix(&key_pair.public_key.m0);
}

void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            uint bit_count) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      mpz_urandomb(m->data[i][j], random_state, bit_count);
      mpz_mod(m->data[i][j], m->data[i][j], m->moduli);
    }
  }
}

PublicPrivateKeyPair key_gen(SignatureParameters params) {
  uint lambda = params.lambda;
  PublicPrivateKeyPair result;
  result.lambda = lambda;

  // private key generation
  bool *private_seed = allocate_seed(lambda);
  generate_seed(private_seed, lambda);
  result.private_key = private_seed;

  // public key generation
  // 1. seed generation
  bool *public_seed = allocate_seed(lambda);
  generate_seed(public_seed, lambda);

  // 2. random matrix generation
  Matrix m0;
  allocate_matrix(&m0, params.matrix_dimension);
  generate_prime(m0.moduli, params.prime_bit_count);

  // 2.1. gmp random state initialization
  gmp_randstate_t public_random_state;
  gmp_randinit_default(public_random_state);

  gmp_randstate_t private_random_state;
  gmp_randinit_default(private_random_state);

  // 2.2. gmp random state seeding
  char *seed_str = malloc(lambda * sizeof(char));

  mpz_t seed;
  mpz_init(seed);

  // public seed
  seed_as_string(seed_str, lambda, public_seed);
  mpz_set_str(seed, seed_str, 2);
  gmp_randseed(public_random_state, seed);

  // private seed
  seed_as_string(seed_str, lambda, private_seed);
  mpz_set_str(seed, seed_str, 2);
  gmp_randseed(private_random_state, seed);

  free(seed_str);

  // 2.3. gmp random integer generation
  // TODO: - generate with `public_seed`:
  //         - (DONE) k matrices M_1, ... M_k of size (m, n)
  //       - generate with `private seed`:
  //         - (DONE) a vector alpha of size k
  //         - a matrix K of size (r, n - r)
  //         - a matrix E' of size (m, r)

  mpz_t *alpha = malloc(params.solution_size * sizeof(mpz_t));
  mpz_init_set_ui(alpha[0], 1);

  Matrix sum;
  allocate_matrix(&sum, params.matrix_dimension);
  fill_matrix_with_zero(&sum);

  Matrix m_i;
  allocate_matrix(&m_i, params.matrix_dimension);
  for (uint i = 1; i <= params.solution_size; i++) {
    // generate M_i
    generate_random_matrix(&m_i, public_random_state, params.prime_bit_count);

    // generate alpha_i
    mpz_init(alpha[i]);
    mpz_urandomb(alpha[i], private_random_state, params.prime_bit_count);

    // compute sum -= alpha_i * M_i
    scalar_product(&m_i, alpha[i], m_i);
    matrix_opposite(&m_i);
    matrix_sum(&sum, sum, m_i);
  }

  clear_matrix(&sum);
  clear_matrix(&m_i);
  gmp_randclear(public_random_state);
  gmp_randclear(private_random_state);

  return result;
}
