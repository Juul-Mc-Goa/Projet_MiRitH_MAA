#include "random.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

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

uint generate_random_element(gmp_randstate_t random_state, FiniteField field) {
  mpz_t temp;
  mpz_init(temp);
  mpz_urandomb(temp, random_state, field.log_field_size);

  uint result = mpz_get_ui(temp);
  mpz_clear(temp);

  return result;
}

void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            FiniteField field) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = generate_random_element(random_state, field);
    }
  }
}

void generate_random_instance(MinRankInstance *instance, uint solution_size, gmp_randstate_t random_state,
                              FiniteField field) {
  for (uint i = 0; i < solution_size; i++) {
    generate_random_matrix(&instance->matrix_array[i], random_state, field);
  }
}
