/* Check that the sum of the additive sharing of `V` really is * `V`. */
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
  printf("---------------------------------------- beginning the share "
         "alpha/A/C/K test...\n");
  uint solution_size = 4;
  uint matrix_count = solution_size + 1;
  uint number_of_parties = 3;
  uint target_rank = 1;
  uint s = 2;

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
  allocate_matrix(&M, GF_16, input_matrix_size);
  allocate_matrix(&shared_M, GF_16, input_matrix_size);

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

  // generate K, C and sum V
  Matrix K, C, V;
  MatrixSize C_size = {s, input_matrix_size.n - target_rank};
  allocate_matrix(&K, GF_16, C_size);
  allocate_matrix(&C, GF_16, C_size);
  allocate_matrix(&V, GF_16, C_size);

  generate_random_matrix(&K, random_state, GF_16);
  matrix_product(&C, A, K);

  // generate a share of `C` and `K`, and use them to compute `V`
  share_c_k_and_update(C, K, R, S, random_state, number_of_parties, parties);

  // compute S*K
  matrix_product(&V, S, K);

  // compute R*M_left
  Matrix temp;
  allocate_matrix(&temp, GF_16, V.size);
  matrix_product(&temp, R, M_left);

  // V -= R*M_left
  matrix_opposite(&temp);
  matrix_sum(&V, V, temp);

  // V -= C
  matrix_opposite(&C);
  matrix_sum(&V, V, C);

  Matrix shared_V;
  allocate_matrix(&shared_V, GF_16, C_size);
  compute_global_v(&shared_V, parties, number_of_parties);

  printf("\nV (not shared):\n");
  print_matrix(&V);
  printf("\nV (shared):\n");
  print_matrix(&shared_V);

  // manually free the right part
  free(M_right.data);

  clear_matrix(&alpha);
  clear_matrix(&A);
  clear_matrix(&C);
  clear_matrix(&M);
  clear_matrix(&S);
  clear_matrix(&shared_M);
  clear_matrix(&shared_V);
  clear_matrix(&temp);
  clear_matrix(&V);
  clear_parties(parties, number_of_parties);
  clear_instance(&instance);
  gmp_randclear(random_state);
}
