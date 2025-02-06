/* Check that the sum of the additive sharing of `S` really is * `S`. */
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
  printf("-------------------------------------------- beginning the share"
         " alpha/A test...\n");
  uint solution_size = 4;
  uint number_of_parties = 3;
  uint target_rank = 1;
  uint s = 2;

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate a random instance
  MinRankInstance instance;
  MatrixSize input_matrix_size = {3, 3};
  init_instance(&instance, solution_size, input_matrix_size);
  generate_random_instance(&instance, solution_size, random_state, GF_16);

  // generate a random vector `alpha`
  Matrix alpha;
  MatrixSize alpha_size = {1, solution_size};
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

  Matrix M_left, M_right;
  M_right.data = malloc(M.size.m * sizeof(uint *));
  split_matrix(&M_left, &M_right, M, target_rank);

  // generate two random matrices `A, R`
  Matrix A, R, S;
  MatrixSize A_size = {s, target_rank}, R_size = {s, input_matrix_size.m};
  allocate_matrix(&A, GF_16, A_size);
  allocate_matrix(&R, GF_16, R_size);
  allocate_matrix(&S, GF_16, A_size);

  generate_A_R_and_sum_S(A, R, S, M_right, parties, number_of_parties,
                         random_state);

  // manually free the right part
  free(M_right.data);

  clear_matrix(&alpha);
  clear_matrix(&M);
  clear_matrix(&S);
  clear_matrix(&shared_M);
  clear_parties(parties, number_of_parties);
  clear_instance(&instance);
  gmp_randclear(random_state);
}
