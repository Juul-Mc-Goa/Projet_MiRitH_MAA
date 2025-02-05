#ifndef MPC_H_
#define MPC_H_

#include "matrix.h"

typedef unsigned int uint;

typedef struct MinRankInstance {
  uint solution_size;
  Matrix *matrix_array;
} MinRankInstance;

typedef struct MinRankSolution {
  uint solution_size;
  uint target_rank;
  Matrix alpha;
  Matrix K;
} MinRankSolution;

typedef struct PartyState {
  Matrix M_left;
  Matrix M_right;
  Matrix S;
  Matrix V;
} PartyState;

void clear_party_state(PartyState *s);
void init_party_state(PartyState *state, uint s, MatrixSize size,
                      uint target_rank);
void init_parties(PartyState *parties, uint number_of_parties, uint s,
                  MatrixSize size, uint target_rank);
void clear_parties(PartyState *parties, uint size);

void split_each_matrix(Matrix *result_left, Matrix *result_right, uint mid,
                       Matrix *input, uint length);

bool mpc_check_solution(gmp_randstate_t random_state, uint number_of_parties,
                        Matrix R, MinRankInstance instance,
                        MinRankSolution solution);

void share_alpha_and_update(Matrix alpha, gmp_randstate_t random_state,
                            uint number_of_parties, MinRankInstance instance,
                            PartyState *parties);
void share_a_and_update(Matrix A, Matrix R, gmp_randstate_t random_state,
                        uint number_of_parties, PartyState *parties);
void share_c_k_and_update(Matrix C, Matrix K, Matrix R, Matrix S,
                          gmp_randstate_t random_state, uint number_of_parties,
                          PartyState *parties);

void compute_local_m(PartyState *state, MinRankInstance instance,
                     Matrix local_alpha, uint solution_size);
void compute_local_s(PartyState *state, Matrix R, Matrix local_A);
void compute_global_s(Matrix *S, PartyState *parties, uint length);
void compute_local_v(PartyState *state, Matrix global_S, Matrix local_K,
                     Matrix R, Matrix local_C);
void compute_global_v(Matrix *V, PartyState *parties, uint length);

#endif // MPC_H_
