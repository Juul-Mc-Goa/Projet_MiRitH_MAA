#include "mpc.h"
#include "constants.h"
#include "key_generation.h"
#include "matrix.h"
#include <stdbool.h>
#include <stdlib.h>

/* Share the vector `alpha` into `number_of_parties` parts, then use
 * it to update the array `parties` of `PartyState` with their share
 * of `M_left` and `M_right`. */
void share_alpha_and_update(Matrix alpha, gmp_randstate_t random_state,
                            uint number_of_parties, MinRankInstance instance,
                            PartyState *parties) {
  uint solution_size = instance.solution_size; // number of input matrices
  Matrix local_alpha, alpha_sum;
  MatrixSize alpha_size = {1, solution_size};
  allocate_matrix(&local_alpha, GF_16, alpha_size);
  allocate_matrix(&alpha_sum, GF_16, alpha_size);

  for (uint i = 0; i < number_of_parties - 1; i++) {
    generate_random_matrix(&local_alpha, random_state, GF_16);
    compute_local_m(&parties[i], instance, local_alpha, solution_size);
    matrix_sum(&alpha_sum, alpha_sum, local_alpha);
  }

  // last share must be `alpha - sum`
  // (here this is `alpha + sum` as we are in characteristic 2)
  matrix_opposite(&alpha_sum);
  matrix_sum(&local_alpha, alpha_sum, alpha);
  compute_local_m(&parties[number_of_parties - 1], instance, local_alpha,
                  solution_size);

  clear_matrix(&local_alpha);
  clear_matrix(&alpha_sum);
}

/* Share the matrix `A` into `number_of_parties` parts, then use it to
 * update `parties` with their share of `S`. */
void share_a_and_update(Matrix A, Matrix R, gmp_randstate_t random_state,
                        uint number_of_parties, PartyState *parties) {
  Matrix A_sum, local_A;
  allocate_matrix(&local_A, GF_16, A.size);
  allocate_matrix(&A_sum, GF_16, A.size);
  fill_matrix_with_zero(&A_sum);

  for (uint i = 0; i < number_of_parties - 1; i++) {
    generate_random_matrix(&local_A, random_state, GF_16);
    matrix_sum(&A_sum, A_sum, local_A);

    compute_local_s(&parties[i], R, local_A);
  }

  matrix_opposite(&A_sum);
  matrix_sum(&local_A, A, A_sum);
  compute_local_s(&parties[number_of_parties - 1], R, local_A);

  clear_matrix(&A_sum);
  clear_matrix(&local_A);
}

/* Share the two matrices `C, K` into `number_of_parties` parts, then use them
 * to update `parties` with their share of `V`. */
void share_c_k_and_update(Matrix C, Matrix K, Matrix R, Matrix S,
                          gmp_randstate_t random_state, uint number_of_parties,
                          PartyState *parties) {
  Matrix C_sum, local_C;
  allocate_matrix(&local_C, GF_16, C.size);
  allocate_matrix(&C_sum, GF_16, C.size);
  fill_matrix_with_zero(&C_sum);

  Matrix K_sum, local_K;
  allocate_matrix(&local_K, GF_16, K.size);
  allocate_matrix(&K_sum, GF_16, K.size);
  fill_matrix_with_zero(&K_sum);

  for (uint i = 0; i < number_of_parties - 1; i++) {
    generate_random_matrix(&local_C, random_state, GF_16);
    generate_random_matrix(&local_K, random_state, GF_16);
    matrix_sum(&C_sum, C_sum, local_C);
    matrix_sum(&K_sum, K_sum, local_K);

    compute_local_v(&parties[i], S, local_K, R, local_C);
  }

  matrix_opposite(&C_sum);
  matrix_opposite(&K_sum);
  matrix_sum(&local_C, C, C_sum);
  matrix_sum(&local_K, K, K_sum);
  compute_local_v(&parties[number_of_parties - 1], S, local_K, R, local_C);

  clear_matrix(&C_sum);
  clear_matrix(&K_sum);
  clear_matrix(&local_C);
  clear_matrix(&local_K);
}

/* Check that the given `solution` to the MinRank `instance` is correct, with
 * multi-party computation. Arguments:
 * - `gmp_randstate_t` is needed to generate shares,
 * - `number_of_parties` specifies how many parties are to be created,
 * - `R` is a random matrix typically send by the verifier,
 * - `instance` is essentially an array of matrices,
 * - `solution` contains a vector `alpha` and a matrix `K`.
 * */
bool mpc_check_solution(gmp_randstate_t random_state, uint number_of_parties,
                        Matrix R, MinRankInstance instance,
                        MinRankSolution solution) {
  MatrixSize size = instance.matrix_array[0].size; // size of each matrix `M_i`
  uint n = size.n,
       s = R.size.m,            // intermediate matrix size
      r = solution.target_rank; // the rank of the solution

  PartyState *parties;
  parties = malloc(sizeof(PartyState) * number_of_parties);
  init_parties(parties, number_of_parties, s, size, r);

  // generate a share of `solution.alpha` and let each party use its share
  // to compute `M_left, M_right`
  share_alpha_and_update(solution.alpha, random_state, number_of_parties,
                         instance, parties);

  // generate a random matrix `A`
  Matrix A;
  MatrixSize A_size = {s, r};
  allocate_matrix(&A, GF_16, A_size);
  generate_random_matrix(&A, random_state, GF_16);

  // generate a share of `A`, and use it to compute `S`
  share_a_and_update(A, R, random_state, number_of_parties, parties);

  // compute the sum of `parties[i].S`
  Matrix S;
  allocate_matrix(&S, GF_16, A_size);
  compute_global_s(&S, parties, number_of_parties);

  // compute `C = AK`
  Matrix C;
  MatrixSize C_size = {s, n - r};
  allocate_matrix(&C, GF_16, C_size);
  matrix_product(&C, A, solution.K);

  // generate a share of `C` and `K`, and use them to compute `V`
  share_c_k_and_update(C, solution.K, R, S, random_state, number_of_parties,
                       parties);

  Matrix V;
  allocate_matrix(&V, GF_16, C_size);
  compute_global_v(&V, parties, number_of_parties);

  bool is_valid = matrix_is_zero(V);

  clear_parties(parties, number_of_parties);
  clear_matrix(&A);
  clear_matrix(&S);
  clear_matrix(&C);
  clear_matrix(&V);

  return is_valid;
}

void init_party_state(PartyState *state, uint s, MatrixSize size,
                      uint target_rank) {

  MatrixSize size_left = {size.m, size.n - target_rank},
             size_right = {size.m, target_rank}, size_S = {s, target_rank},
             size_V = {s, size.n - target_rank};

  allocate_matrix(&state->M_left, GF_16, size_left);
  allocate_matrix(&state->M_right, GF_16, size_right);
  allocate_matrix(&state->S, GF_16, size_S);
  allocate_matrix(&state->V, GF_16, size_V);
}

void clear_party_state(PartyState *s) {
  clear_matrix(&s->M_left);
  clear_matrix(&s->M_right);
  clear_matrix(&s->S);
  clear_matrix(&s->V);
}

void init_parties(PartyState *parties, uint number_of_parties, uint s,
                  MatrixSize size, uint target_rank) {
  for (uint i = 0; i < number_of_parties; i++) {
    init_party_state(&parties[i], s, size, target_rank);
  }
}

void clear_parties(PartyState *parties, uint size) {
  for (uint i = 0; i < size; i++) {
    clear_party_state(&parties[i]);
  }

  free(parties);
}

void init_instance(MinRankInstance *instance, uint solution_size,
                   MatrixSize input_matrix_size) {
  instance->solution_size = solution_size;
  instance->matrix_array = malloc(sizeof(Matrix) * solution_size);

  for (uint i = 0; i < solution_size; i++) {
    allocate_matrix(&instance->matrix_array[i], GF_16, input_matrix_size);
  }
}

void clear_instance(MinRankInstance *instance) {
  for (uint i = 0; i < instance->solution_size; i++) {
    clear_matrix(&instance->matrix_array[i]);
  }
  free(instance->matrix_array);
}

/* compute the weighted sum of `instance.matrix_array` with `local_alpha`. */
void compute_local_m(PartyState *state, MinRankInstance instance,
                     Matrix local_alpha, uint solution_size) {
  // split each `instance.matrix_array[i]` into a left and right part
  Matrix *left_part = malloc(sizeof(Matrix) * solution_size),
         *right_part = malloc(sizeof(Matrix) * solution_size);

  // we have to allocate each uint** in the right part
  for (uint i = 0; i < solution_size; i++) {
    right_part[i].data = malloc(state->M_right.size.m * sizeof(uint *));
  }

  split_each_matrix(left_part, right_part, instance.matrix_array,
                    state->M_left.size.n, instance.solution_size);

  // compute the weighted sum on each part
  matrix_big_weighted_sum(&state->M_left, local_alpha.data[0], left_part,
                          instance.solution_size);
  matrix_big_weighted_sum(&state->M_right, local_alpha.data[0], right_part,
                          instance.solution_size);
}

void compute_local_s(PartyState *state, Matrix R, Matrix local_A) {
  matrix_product(&state->S, R, state->M_right);
  matrix_sum(&state->S, state->S, local_A);
}

/* Just sum each local `S`. */
void compute_global_s(Matrix *S, PartyState *parties, uint length) {
  fill_matrix_with_zero(S);
  for (uint i = 0; i < length; i++) {
    matrix_sum(S, *S, parties[i].S);
  }
}

/* Update `state` by assigning `state.V = S*K - R*M_left - C`. */
void compute_local_v(PartyState *state, Matrix global_S, Matrix local_K,
                     Matrix R, Matrix local_C) {
  // compute S*K
  matrix_product(&state->V, global_S, local_K);

  // compute R*M_left
  Matrix temp;
  allocate_matrix(&temp, GF_16, state->V.size);
  matrix_product(&temp, R, state->M_left);

  // V -= R*M_left
  matrix_opposite(&temp);
  matrix_sum(&state->V, state->V, temp);

  // V -= C
  matrix_opposite(&local_C);
  matrix_sum(&state->V, state->V, local_C);

  clear_matrix(&temp);
}

/* Just sum each local `V`. */
void compute_global_v(Matrix *V, PartyState *parties, uint length) {
  fill_matrix_with_zero(V);
  for (uint i = 0; i < length; i++) {
    matrix_sum(V, *V, parties[i].V);
  }
}
