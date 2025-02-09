#ifndef MPC_COMMON_H_
#define MPC_COMMON_H_

#include "../matrix.h"
#include "../mpc.h"

void init_parties_and_sum_M(PartyState *parties, MinRankInstance instance,
                            Matrix M, Matrix shared_M, Matrix alpha, uint s,
                            uint number_of_parties, uint target_rank,
                            gmp_randstate_t random_state);

void generate_A_R_and_sum_S(Matrix A, Matrix R, Matrix S, Matrix M_right,
                            PartyState *parties, uint number_of_parties,
                            gmp_randstate_t random_state);
#endif // MPC_COMMON_H_
