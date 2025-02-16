#include "packing.h"
#include "gmp.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void allocate_packed_last_state(uchar *result, uint lambda,
                                MatrixSize alpha_size, MatrixSize K_size,
                                MatrixSize C_size) {
  uint alpha_length = alpha_size.m * alpha_size.n;
  uint K_length = K_size.m * K_size.n;
  uint C_length = C_size.m * C_size.n;
  uint total = alpha_length + K_length + C_length;
  result = malloc(sizeof(uint8_t) * (lambda + total / 2));
}

void pack_boolean(uchar *result, bool element, uint *bit_index) {
  uint current_index = *bit_index >> 3, current_offset = *bit_index & 7;

  uint8_t new_value = (uint8_t)result[current_index];
  new_value ^= ((uint8_t)element << (7 - current_offset));
  result[current_index] = (uchar)new_value;

  *bit_index += 1;
  if (current_offset == 7) {
    result[current_index + 1] = (uchar)0;
  }
}

void pack_four_bit_value(uchar *result, uint element, uint *bit_index) {
  uint current_index = *bit_index >> 3;
  uint current_offset = *bit_index & 7;

  if (current_offset == 0) {
    uint8_t new_value = ((uint8_t)element) << 4;
    result[current_index] = (uchar)new_value;
  } else {
    uint8_t new_value = (uint8_t)result[current_index] ^ (uint8_t)element;
    result[current_index] = (uchar)new_value;
    result[current_index + 1] = (uchar)0;
  }
  *bit_index += 4;
}

void pack_32bit_value(uchar *result, uint element, uint *bit_index) {
  for (uint i = 0; i < 8; i++) {
    pack_four_bit_value(result, element & 15, bit_index);
    element >>= 4;
  }
}

void pack_last_state(uchar *result, seed_t salt, uint lambda, uint l, uint i,
                     Matrix alpha, Matrix K, Matrix C) {

  // pack `salt`
  memcpy(result, salt.data, salt.size);

  // pack `l, i`
  result[salt.size] = (uchar)((l >> 24) & 255);
  result[salt.size + 1] = (uchar)((l >> 16) & 255);
  result[salt.size + 2] = (uchar)((l >> 8) & 255);
  result[salt.size + 3] = (uchar)(l & 255);

  result[salt.size + 4] = (uchar)((i >> 24) & 255);
  result[salt.size + 5] = (uchar)((i >> 16) & 255);
  result[salt.size + 6] = (uchar)((i >> 8) & 255);
  result[salt.size + 7] = (uchar)(i & 255);

  uint bit_index = (salt.size << 3) + 64;

  // pack `alpha`
  for (uint i = 0; i < alpha.size.n; i++) {
    pack_four_bit_value(result, alpha.data[0][i], &bit_index);
  }
  // pack `K`
  for (uint i = 0; i < K.size.m; i++) {
    for (uint j = 0; j < K.size.n; j++) {
      pack_four_bit_value(result, K.data[i][j], &bit_index);
    }
  }
  // pack `C`
  for (uint i = 0; i < C.size.m; i++) {
    for (uint j = 0; j < C.size.n; j++) {
      pack_four_bit_value(result, C.data[i][j], &bit_index);
    }
  }
}

/* Unpack a `MinRankSolution` from the given private key. The `solution`
 * must be initialized (by `init_solution` for example). This uses
 * `priv_key.seed` to seed a random state, then generates random `alpha` and
 * `K`. */
void unpack_solution_from_private_key(MinRankSolution *solution,
                                      gmp_randstate_t prg_state,
                                      SignatureParameters params,
                                      PrivateKey priv_key) {
  solution->alpha.data[0][0] = 1;
  for (uint i = 1; i <= params.solution_size; i++) {
    solution->alpha.data[0][i] =
        generate_random_element(prg_state, params.field);
  }

  generate_random_matrix(&solution->K, prg_state, params.field);
}

void unpack_instance_from_public_key(MinRankInstance *instance,
                                     gmp_randstate_t prg_state,
                                     SignatureParameters params,
                                     PublicKey pub_key) {
  instance->matrix_count = params.solution_size + 1;

  // copy `m0` from the public key
  matrix_init_set(&instance->matrix_array[0], pub_key.m0);

  for (uint i = 1; i <= params.solution_size; i++) {
    generate_random_matrix(&instance->matrix_array[i], prg_state, params.field);
  }
}
