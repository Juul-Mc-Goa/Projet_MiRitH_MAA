#include "key_generation.h"
#include "field_arithmetics.h"
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
    printf("\nDid not fill the required number of random booleans.\n");
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

void seed_random_state(bool *seed, uint lambda, gmp_randstate_t random_state) {
  char *seed_str = malloc(lambda * sizeof(char));
  mpz_t gmp_seed;
  mpz_init(gmp_seed);

  // public seed
  seed_as_string(seed_str, lambda, seed);
  mpz_set_str(gmp_seed, seed_str, 2);
  gmp_randseed(random_state, gmp_seed);

  free(seed_str);
}

void allocate_key_pair(PublicPrivateKeyPair *key_pair,
                       SignatureParameters params) {
  key_pair->private_key.lambda = params.lambda;
  key_pair->private_key.seed = malloc(params.lambda * sizeof(bool));

  key_pair->public_key.lambda = params.lambda;
  key_pair->public_key.seed = malloc(params.lambda * sizeof(bool));
  allocate_matrix(&key_pair->public_key.m0, params.field,
                  params.matrix_dimension);
}

void clear_key_pair(PublicPrivateKeyPair key_pair) {
  free(key_pair.private_key.seed);
  free(key_pair.public_key.seed);
  clear_matrix(&key_pair.public_key.m0);
}

uint generate_random_element(gmp_randstate_t random_state, FiniteField field) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, random_state, field.log_field_size);

  uint result = mpz_get_ui(temp);
  mpz_clear(temp);

  return result;
}

void generate_random_vector(uint *vector, uint vector_size,
                            gmp_randstate_t random_state, FiniteField field) {
  for (uint i = 0; i < vector_size; i++) {
    vector[i] = generate_random_element(random_state, field);
  }
}

void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            FiniteField field) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = generate_random_element(random_state, field);
    }
  }
}

PublicPrivateKeyPair key_gen(SignatureParameters params) {
  uint lambda = params.lambda;
  PublicPrivateKeyPair result;
  allocate_key_pair(&result, params);

  // private key generation
  result.private_key.lambda = lambda;
  generate_seed(result.private_key.seed, lambda);

  // public key generation
  // 1. seed generation
  result.public_key.lambda = lambda;
  generate_seed(result.public_key.seed, lambda);

  // 2. random matrix generation
  // 2.1. gmp random state initialization
  gmp_randstate_t public_random_state;
  gmp_randinit_default(public_random_state);

  gmp_randstate_t private_random_state;
  gmp_randinit_default(private_random_state);

  // 2.2. gmp random state seeding
  // public seed
  seed_random_state(result.public_key.seed, lambda, public_random_state);
  // private seed
  seed_random_state(result.private_key.seed, lambda, private_random_state);

  // 2.3. gmp random integer generation
  // 2.3.1 generate M_1, ..., M_k, and field elements alpha_1, ..., alpha_k,
  //       then store sum = alpha_1 * M_1 + ... + alpha_k * M_k
  uint *alpha = malloc(params.solution_size * sizeof(uint));
  alpha[0] = 1;

  Matrix sum;
  allocate_matrix(&sum, params.field, params.matrix_dimension);
  fill_matrix_with_zero(&sum);

  Matrix m_i;
  allocate_matrix(&m_i, params.field, params.matrix_dimension);
  for (uint i = 1; i <= params.solution_size; i++) {
    // generate M_i
    generate_random_matrix(&m_i, public_random_state, params.field);

    // generate alpha_i
    alpha[i] = generate_random_element(private_random_state, params.field);

    // compute sum += alpha_i * M_i
    scalar_product(&m_i, alpha[i], m_i);
    matrix_sum(&sum, sum, m_i);
  }

  // 2.3.2 generate a random matrix K
  Matrix K;
  MatrixSize K_size = {params.target_rank,
                       params.matrix_dimension.n - params.target_rank};
  allocate_matrix(&K, params.field, K_size);
  generate_random_matrix(&K, private_random_state, params.field);

  // 2.3.3 generate a random matrix E_R
  Matrix E_R;
  MatrixSize E_R_size = {params.matrix_dimension.m, params.target_rank};
  allocate_matrix(&E_R, params.field, E_R_size);
  generate_random_matrix(&E_R, private_random_state, params.field);

  // 2.3.4 compute E from K and E_R
  Matrix E;
  allocate_matrix(&E, params.field, params.matrix_dimension);

  // left side: compute E_L = E_R * K
  // E_L is just a view of E so it does not need to be allocated or freed
  Matrix E_L;
  MatrixSize left_size = {E.size.m, K_size.n};

  E_L.field = params.field;
  E_L.size = left_size;
  E_L.data = E.data;
  matrix_product(&E_L, E_R, K);

  // right side: copy E_R into E
  for (uint i = 0; i < E.size.m; i++) {
    for (uint j = K_size.n; j < E.size.n; j++) {
      E.data[i][j] = E_R.data[i][j - K_size.n];
    }
  }

  // 2.3.5 compute M_0 = E - sum: here field_size is a power of two, so -1 = 1,
  // and M_0 = E + sum
  matrix_opposite(&sum);
  matrix_sum(&result.public_key.m0, E, sum);

  clear_matrix(&sum);
  clear_matrix(&m_i);
  clear_matrix(&K);
  clear_matrix(&E_R);
  clear_matrix(&E);
  gmp_randclear(public_random_state);
  gmp_randclear(private_random_state);

  return result;
}
