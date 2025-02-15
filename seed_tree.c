#include "constants.h"
#include "random.h"
#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned char uchar;

/* Convert a byte array into a GMP integer. */
void string_to_mpz(uchar *string, size_t string_size, mpz_t *big_int) {
  // each uchar is converted to 2 hex digits
  char *hex_string = malloc(sizeof(char) * 2 * string_size);

  for (uint i = 0; i < string_size; i++) {
    hex_string[2 * i] = HEX_CHAR_TABLE[(uint8_t)string[i] >> 4];
    hex_string[2 * i + 1] = HEX_CHAR_TABLE[(uint8_t)string[i] & 15];
  }

  mpz_set_str(*big_int, hex_string, 16);

  free(hex_string);
}

// Function to initialize the PRG state
void PRG_init(const uchar *salt, const uchar *seed, uint lambda,
              gmp_randstate_t *prg_state) {
  size_t size = ceil(lambda / 8.0);
  uchar *salt_and_seed = malloc(sizeof(uchar) * 3 * size);

  // copy salt
  memcpy(salt_and_seed, salt, 2 * size);
  // copy seed
  memcpy(&salt_and_seed[2 * size], salt, size);
  // use the result to seed prg_state
  mpz_t gmp_seed;
  mpz_init(gmp_seed);
  string_to_mpz(salt_and_seed, 3 * size, &gmp_seed);
  gmp_randseed(*prg_state, gmp_seed);

  mpz_clear(gmp_seed);
}

// Function to generate bytes from the PRG state
void PRG_bytes(gmp_randstate_t prg_state, size_t length,
               unsigned char *output) {
  // the GMP integer holding the random byte
  mpz_t temp;
  mpz_init(temp);

  for (uint i = 0; i < length; i++) {
    // generate a random byte
    mpz_urandomb(temp, prg_state, 8);
    int8_t random_byte = mpz_get_ui(temp);
    output[i] = (uchar)random_byte;
  }
  mpz_clear(temp);
}

// Function to generate the seeds of the parties for every round
void TreePRG(const unsigned char *salt, const unsigned char *seed,
             size_t lambda, size_t N, unsigned char **output_seeds) {
  size_t N_bitlength = (size_t)ceil(log2(N));

  // tree_size is (almost) 2 * 2^{N_bitlength}
  size_t tree_size = (size_t)(1 << (1 + N_bitlength)) - 1;

  // allocate tree
  unsigned char **tree = (uchar **)malloc(tree_size * sizeof(unsigned char *));
  for (size_t i = 0; i < tree_size; i++) {
    tree[i] = (uchar *)malloc(lambda * sizeof(unsigned char));
  }

  // initialize prg_state
  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  // Initialize the tree with the initial seed
  memcpy(tree[0], seed, lambda);

  // Generate the tree
  for (size_t i = 0; i < tree_size; i++) {
    PRG_init(salt, tree[i], lambda, &prg_state);

    unsigned char child_seed_left[lambda];
    unsigned char child_seed_right[lambda];
    PRG_bytes(prg_state, lambda, child_seed_left);
    PRG_bytes(prg_state, lambda, child_seed_right);

    if (2 * i + 1 < tree_size) {
      memcpy(tree[2 * i + 1], child_seed_left, lambda);
    }
    if (2 * i + 2 < tree_size) {
      memcpy(tree[2 * i + 2], child_seed_right, lambda);
    }
  }

  // Extract the seeds for the parties
  size_t j = (size_t)(1 << N_bitlength);
  for (size_t i = 0; i < N; i++) {
    memcpy(output_seeds[i], tree[j + i], lambda);
  }

  // Free the allocated memory
  for (size_t i = 0; i < tree_size; i++) {
    free(tree[i]);
  }
  free(tree);
}

// Function to pack the seeds used by all the parties but the i_star in a
// particular round
void Tree_pack(const unsigned char **tree, size_t i_star, size_t N,
               size_t lambda, unsigned char *packed) {
  size_t initial_index = 1;
  size_t packed_length = 0;
  size_t logN = (size_t)ceil(log2(N));

  for (size_t j = 1; j <= logN; j++) {
    size_t t = (size_t)floor(i_star / pow(2, logN - j));
    if (t % 2 == 1) {
      memcpy(packed + packed_length * lambda, tree[2 * initial_index], lambda);
      packed_length++;
      initial_index = 2 * initial_index + 1;
    } else {
      memcpy(packed + packed_length * lambda, tree[2 * initial_index + 1],
             lambda);
      packed_length++;
      initial_index = 2 * initial_index;
    }
  }
}

// Function to unpack the seeds used by all the parties but the i_star in a
// particular round
void Tree_unpack(const unsigned char *packed, size_t i_star, size_t N,
                 size_t lambda, unsigned char **unpacked) {
  size_t logN = (size_t)ceil(log2(N));
  size_t temp_size = (size_t)pow(2, logN);
  unsigned char **temp =
      (unsigned char **)malloc(temp_size * sizeof(unsigned char *));
  for (size_t i = 0; i < temp_size; i++) {
    temp[i] = (unsigned char *)malloc(lambda * sizeof(unsigned char));
    if (i < N) {
      unpacked[i] = NULL;
    }
  }

  size_t j = logN - 1;
  size_t gap = 0;
  size_t packed_index = 0;

  while (j != (size_t)-1) {
    size_t j_0 = (size_t)floor(i_star / pow(2, j)) % 2;
    size_t aux_initial_index = (1 - j_0) * (size_t)pow(2, j);
    size_t initial_index = aux_initial_index + gap;
    gap += j_0 * (size_t)pow(2, j);

    unsigned char **tree_0 =
        (unsigned char **)malloc((size_t)pow(2, j) * sizeof(unsigned char *));
    for (size_t i = 0; i < (size_t)pow(2, j); i++) {
      tree_0[i] = (unsigned char *)malloc(lambda * sizeof(unsigned char));
    }

    TreePRG(NULL, packed + packed_index * lambda, lambda, (size_t)pow(2, j),
            tree_0);
    packed_index++;

    for (size_t i = 0; i < (size_t)pow(2, j); i++) {
      memcpy(temp[initial_index + i], tree_0[i], lambda);
    }

    for (size_t i = 0; i < (size_t)pow(2, j); i++) {
      free(tree_0[i]);
    }
    free(tree_0);

    if (j == 0) {
      break;
    }
    j--;
  }

  for (size_t i = 0; i < N; i++) {
    if (i != i_star) {
      memcpy(unpacked[i], temp[i], lambda);
    }
  }

  for (size_t i = 0; i < temp_size; i++) {
    free(temp[i]);
  }
  free(temp);
}
