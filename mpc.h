#ifndef MPC_H_
#define MPC_H_

#include "matrix.h"
#include "types.h"

/* Main function */
bool mpc_check_solution(gmp_randstate_t random_state, uint number_of_parties,
                        Matrix R, MinRankInstance instance,
                        MinRankSolution solution);

/* Utilities used by `mpc_check_solution` */
void share_alpha_and_update(Matrix alpha, gmp_randstate_t random_state,
                            uint number_of_parties, MinRankInstance instance,
                            PartyState *parties);
void share_a_and_update(Matrix A, Matrix R, gmp_randstate_t random_state,
                        uint number_of_parties, PartyState *parties);
void share_c_k_and_update(Matrix C, Matrix K, Matrix R, Matrix S,
                          gmp_randstate_t random_state, uint number_of_parties,
                          PartyState *parties);
void compute_local_m(PartyState *state, MinRankInstance instance,
                     Matrix local_alpha, uint matrix_count);
void compute_local_s(PartyState *state, Matrix R, Matrix local_A);
void compute_global_s(Matrix *S, PartyState *parties, uint length);
void compute_local_v(PartyState *state, Matrix global_S, Matrix local_K,
                     Matrix R, Matrix local_C);
void compute_global_v(Matrix *V, PartyState *parties, uint length);

/* Memory management */
void init_party_state(PartyState *state, uint s, MatrixSize size,
                      uint target_rank);
void clear_party_state(PartyState *s);

void init_parties(PartyState *parties, uint number_of_parties, uint s,
                  MatrixSize size, uint target_rank);
void clear_parties(PartyState *parties, uint size);

void init_instance(MinRankInstance *instance, uint matrix_count,
                   MatrixSize input_matrix_size);
void clear_instance(MinRankInstance *instance);

void init_solution(MinRankSolution *solution, uint matrix_count,
                   uint target_rank, MatrixSize input_matrix_size);
void clear_solution(MinRankSolution *solution);

#endif // MPC_H_
