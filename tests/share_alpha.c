/* Check that the sums of the additive sharings of `M_left, M_right` really are
 * `M_left, M_right`. */
#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"
#include "../mpc.h"
#include "../random.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("---------------------------------------------- beginning the share "
         "alpha test...\n");
  uint solution_size = 4;
  uint number_of_parties = 3;
  uint target_rank = 1;
  uint s = 2;

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate random matrices for the input
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
  // M_right`
  Matrix M, M_left, M_right;
  MatrixSize M_size = {input_matrix_size.m, input_matrix_size.n};
  allocate_matrix(&M, GF_16, M_size);

  matrix_big_weighted_sum(&M, alpha.data[0], instance.matrix_array,
                          instance.solution_size);

  // we have to allocate `M_right.data`: it's a (as of yet not allocated) (uint
  // **) where the data at each index `i` is a pointer to data that belongs to
  // `M`
  M_right.data = malloc(M.size.m * sizeof(uint *));

  split_matrix(&M_left, &M_right, M, target_rank);

  // share `alpha` and compute each party's `M_left, M_right`
  PartyState *parties = malloc(sizeof(PartyState) * number_of_parties);
  init_parties(parties, number_of_parties, s, input_matrix_size, target_rank);

  share_alpha_and_update(alpha, random_state, number_of_parties, instance,
                         parties);

  // sum each party's part
  Matrix shared_M, shared_M_left, shared_M_right;
  allocate_matrix(&shared_M, GF_16, M_size);

  fill_matrix_with_zero(&shared_M);

  // we have to allocate `shared_M_right`
  shared_M_right.data = malloc(shared_M.size.m * sizeof(uint *));

  split_matrix(&shared_M_left, &shared_M_right, shared_M,
               target_rank);
  for (uint i = 0; i < number_of_parties; i++) {
    matrix_sum(&shared_M_left, parties[i].M_left, shared_M_left);
    matrix_sum(&shared_M_right, parties[i].M_right, shared_M_right);
  }

  printf("weighted sum without sharing: \n");
  print_matrix(&M);

  printf("weighted sum with sharing: \n");
  print_matrix(&shared_M);

  // manually free the right parts
  free(shared_M_right.data);
  free(M_right.data);

  clear_matrix(&alpha);
  clear_matrix(&M);
  clear_matrix(&shared_M);
  clear_parties(parties, number_of_parties);
  clear_instance(&instance);
  gmp_randclear(random_state);
}
