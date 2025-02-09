/* Check that the sums of the additive sharings of `M_left, M_right` really are
 * `M_left, M_right`. */
#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"
#include "../mpc.h"
#include "../random.h"
#include "mpc_common.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("---------------------------------------------- beginning the share "
         "alpha test...\n");

  uint number_of_parties = 3, target_rank = 1, s = 2, solution_size = 4;
  uint matrix_count = solution_size + 1;

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate a random instance
  MinRankInstance instance;
  MatrixSize input_matrix_size = {3, 3};
  init_instance(&instance, matrix_count, input_matrix_size);
  generate_random_instance(&instance, matrix_count, random_state, GF_16);

  // generate a random vector `alpha`
  Matrix alpha;
  MatrixSize alpha_size = {1, matrix_count};
  allocate_matrix(&alpha, GF_16, alpha_size);
  generate_random_matrix(&alpha, random_state, GF_16);

  // compute M using the global `alpha`, and split the result into `M_left,
  // M_right`. Check that it is indeed equal to the sum of each party's `M`.
  Matrix M, shared_M;
  MatrixSize M_size = {input_matrix_size.m, input_matrix_size.n};
  allocate_matrix(&M, GF_16, M_size);
  allocate_matrix(&shared_M, GF_16, M_size);

  PartyState *parties = malloc(sizeof(PartyState) * number_of_parties);

  init_parties_and_sum_M(parties, instance, M, shared_M, alpha, s,
                         number_of_parties, target_rank, random_state);

  clear_matrix(&alpha);
  clear_matrix(&M);
  clear_matrix(&shared_M);
  clear_parties(parties, number_of_parties);
  clear_instance(&instance);
  gmp_randclear(random_state);
}
