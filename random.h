#ifndef RANDOM_H_
#define RANDOM_H_

#include "field_arithmetics.h"
#include "key_generation.h"
#include "matrix.h"
#include "mpc.h"

#include <gmp.h>
#include <stdbool.h>

typedef unsigned char uchar;

void allocate_seed(seed_t *seed, uint lambda);
void clear_seed(seed_t *seed);
void print_seed(seed_t seed);
void generate_seed(seed_t seed);
void seed_to_mpz(uchar *string, size_t string_size, mpz_t *big_int);
void seed_random_state(seed_t seed, gmp_randstate_t random_state);
uint generate_random_element(gmp_randstate_t random_state, FiniteField field);
void generate_random_matrix(Matrix *m, gmp_randstate_t random_state,
                            FiniteField field);
void generate_random_instance(MinRankInstance *instance, uint solution_size,
                              gmp_randstate_t random_state, FiniteField field);
#endif // RANDOM_H_
