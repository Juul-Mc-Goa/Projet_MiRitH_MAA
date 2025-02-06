#ifndef RANDOM_H_
#define RANDOM_H_

#include "field_arithmetics.h"
#include "matrix.h"
#include "mpc.h"

#include <gmp.h>
#include <stdbool.h>

void generate_seed(bool *seed, uint lambda);
uint generate_random_element(gmp_randstate_t random_state, FiniteField field);
void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            FiniteField field);
void generate_random_instance(MinRankInstance *instance, uint solution_size, gmp_randstate_t random_state,
                              FiniteField field);

#endif // RANDOM_H_
