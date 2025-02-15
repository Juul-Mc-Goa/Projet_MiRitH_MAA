#include "sign.h"
#include "random.h"
#include "constants.h"
#include "gmp.h"
#include "key_generation.h"
#include "matrix.h"
#include "mpc.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/random.h>

typedef unsigned char uchar;

/* generate the seed, using the `getrandom` system call (available on Linux) */
void generate_seed(seed_t seed) {
  ssize_t result_size = getrandom(seed.data, seed.size * sizeof(uchar), 0);
  if (result_size != seed.size) {
    printf("\nDid not fill the required number of random bytes.\n");
    printf("requested: %u bytes, filled: %ld bytes\n", seed.size, result_size);
    return;
  }
}

void allocate_seed(seed_t *seed, uint lambda) {
  uint size = ceil(lambda / 8.0);
  seed->data = malloc(sizeof(uchar) * size);
  seed->size = size;
}

void clear_seed(seed_t *seed) { free(seed->data); }

void print_seed(seed_t seed) {
  for (uint i = 0; i < seed.size; i++) {
    char left_char = HEX_CHAR_TABLE[(uint8_t)seed.data[i] >> 4];
    char right_char = HEX_CHAR_TABLE[(uint8_t)seed.data[i] & 15];
    printf("%c%c ", left_char, right_char);
  }
  printf("\n");
}

/* Convert a byte array into a GMP integer. */
void seed_to_mpz(uchar *string, size_t string_size, mpz_t *big_int) {
  // each uchar is converted to 2 hex digits
  char *hex_string = malloc(sizeof(char) * (2 * string_size + 1));

  for (uint i = 0; i < string_size; i++) {
    hex_string[2 * i] = HEX_CHAR_TABLE[(uint8_t)string[i] >> 4];
    hex_string[2 * i + 1] = HEX_CHAR_TABLE[(uint8_t)string[i] & 15];
  }
  hex_string[2 * string_size] = (char)0;

  mpz_set_str(*big_int, hex_string, 16);

  free(hex_string);
}

void seed_random_state(seed_t seed, gmp_randstate_t random_state) {
  mpz_t gmp_seed;
  mpz_init(gmp_seed);

  seed_to_mpz(seed.data, seed.size, &gmp_seed);
  gmp_randseed(random_state, gmp_seed);
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
