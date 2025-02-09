#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"
#include "../mpc.h"
#include "../random.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* Initialize `parties`, then share `alpha` and check that the shared M is equal
 * to the unshared one. */
void init_parties_and_sum_M(PartyState *parties, MinRankInstance instance,
                            Matrix M, Matrix shared_M, Matrix alpha, uint s,
                            uint number_of_parties, uint target_rank,
                            gmp_randstate_t random_state) {
  MatrixSize input_matrix_size = instance.matrix_array[0].size;

  // compute M using the global `alpha`
  matrix_big_weighted_sum(&M, alpha.data[0], instance.matrix_array,
                          instance.matrix_count);

  // share `alpha` and compute each party's `M_left, M_right`
  init_parties(parties, number_of_parties, s, input_matrix_size, target_rank);
  share_alpha_and_update(alpha, random_state, number_of_parties, instance,
                         parties);

  // sum each party's part
  Matrix shared_M_left, shared_M_right;
  fill_matrix_with_zero(&shared_M);

  // we have to allocate `shared_M_right`
  shared_M_right.data = malloc(shared_M.size.m * sizeof(uint *));
  split_matrix(&shared_M_left, &shared_M_right, shared_M, target_rank);

  for (uint i = 0; i < number_of_parties; i++) {
    matrix_sum(&shared_M_left, parties[i].M_left, shared_M_left);
    matrix_sum(&shared_M_right, parties[i].M_right, shared_M_right);
  }

  printf("\nM (not shared): \n");
  print_matrix(&M);

  printf("\nM (shared): \n");
  print_matrix(&shared_M);

  free(shared_M_right.data);
}

void generate_A_R_and_sum_S(Matrix A, Matrix R, Matrix S, Matrix M_right,
                            PartyState *parties, uint number_of_parties,
                            gmp_randstate_t random_state) {
  generate_random_matrix(&A, random_state, GF_16);
  generate_random_matrix(&R, random_state, GF_16);

  printf("\nR:\n");
  print_matrix(&R);
  printf("\nA:\n");
  print_matrix(&A);

  // compute `S = R * M_right + A`
  matrix_product(&S, R, M_right);
  matrix_sum(&S, S, A);

  // generate a share of `A`, and use it to compute `S`
  share_a_and_update(A, R, random_state, number_of_parties, parties);

  // compute the sum of `parties[i].S`
  Matrix shared_S;
  allocate_matrix(&shared_S, GF_16, A.size);
  compute_global_s(&shared_S, parties, number_of_parties);

  printf("\nS (not shared):\n");
  print_matrix(&S);
  printf("\nS (shared):\n");
  print_matrix(&shared_S);
}
