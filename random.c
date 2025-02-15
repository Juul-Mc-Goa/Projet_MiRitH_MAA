#include "random.h"
#include "constants.h"
#include "gmp.h"
#include "key_generation.h"
#include "matrix.h"
#include "mpc.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

typedef unsigned char uchar;

/* generate the seed, using the `getrandom` system call (available on Linux) */
void generate_seed(bool *seed, uint lambda) {
  uint byte_size = 8;
  uint size = ceil(lambda / 8.0);

  // the getrandom function outputs a random array of bytes, that we package
  // in a char* variable. then we unpack it into a boolean array.
  char *raw_seed = malloc(size * sizeof(uchar));
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
}

void generate_uchar_seed(uchar *seed, uint lambda) {
  uint size = ceil(lambda / 8.0);
  ssize_t result_size = getrandom(seed, size * sizeof(uchar), 0);
  if (result_size != size) {
    printf("\nDid not fill the required number of random bytes.\n");
    printf("requested: %u bytes, filled: %ld bytes\n", size, result_size);
    return;
  }
}

bool *allocate_seed(uint lambda) {
  bool *result = malloc(lambda * sizeof(bool));
  return result;
}

void allocate_uchar_seed(uchar **result, uint lambda) {
  uint size = ceil(lambda / 8.0);
  *result = malloc(sizeof(uchar) * size);
}

/* Pack a boolean array into a string (in hexadecimal). */
void seed_as_string(char *result, uint lambda, bool *seed) {
  uint result_size = ceil((double)lambda / 4.0);

  uint current_int = 0;
  uint current_char_index = 0;

  for (uint i = 0; i < (4 * result_size); i++) {
    // shift by one to the left to make room for the new boolean
    current_int <<= 1;
    // assign LSB to the new boolean
    if (i < lambda) {
      current_int |= seed[i];
    }
    // if we finished processing a 4bit chunk, we update result
    if ((i & 3) == 3) {
      result[current_char_index] = HEX_CHAR_TABLE[current_int];
      current_int = 0;
      current_char_index += 1;
    }
  }
}

void seed_random_state(bool *seed, uint lambda, gmp_randstate_t random_state) {
  uint str_size = ceil((double)lambda / 4.0);
  char *seed_str = malloc(str_size * sizeof(char));

  mpz_t gmp_seed;
  mpz_init(gmp_seed);

  seed_as_string(seed_str, lambda, seed);
  mpz_set_str(gmp_seed, seed_str, 2);
  gmp_randseed(random_state, gmp_seed);

  free(seed_str);
}

/* Generate a random element in `field`. */
uint generate_random_element(gmp_randstate_t random_state, FiniteField field) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, random_state, field.log_field_size);

  uint result = mpz_get_ui(temp);
  mpz_clear(temp);

  return result % field.field_size;
}

/* Generate a random matrix with coefficients in `field`. */
void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            FiniteField field) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = generate_random_element(random_state, field);
    }
  }
}

/* Uses `random_state` to generate an array of random matrices. */
void generate_random_instance(MinRankInstance *instance, uint matrix_count,
                              gmp_randstate_t random_state, FiniteField field) {
  for (uint i = 0; i < matrix_count; i++) {
    generate_random_matrix(&instance->matrix_array[i], random_state, field);
  }
}

/* Unpack a `MinRankSolution` from the given private key. The `solution`
 * must be initialized (by `init_solution` for example). This uses
 * `priv_key.seed` to seed a random state, then generates random `alpha` and
 * `K`. */
void unpack_solution_from_private_key(MinRankSolution *solution,
                                      gmp_randstate_t random_state,
                                      SignatureParameters params,
                                      PrivateKey priv_key) {
  solution->alpha.data[0][0] = 1;
  for (uint i = 1; i <= params.solution_size; i++) {
    solution->alpha.data[0][i] =
        generate_random_element(random_state, params.field);
  }

  generate_random_matrix(&solution->K, random_state, params.field);
}

void unpack_instance_from_public_key(MinRankInstance *instance,
                                     gmp_randstate_t random_state,
                                     SignatureParameters params,
                                     PublicKey pub_key) {
  instance->matrix_count = params.solution_size + 1;

  // copy `m0` from the public key
  matrix_init_set(&instance->matrix_array[0], pub_key.m0);

  for (uint i = 1; i <= params.solution_size; i++) {
    generate_random_matrix(&instance->matrix_array[i], random_state,
                           params.field);
  }
}
