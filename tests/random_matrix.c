#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("------------------------------ beginning random matrix test...\n");
  Matrix m;
  m.size.m = 4;
  m.size.n = 4;
  FiniteField field = GF_16;
  allocate_matrix(&m, field, m.size);

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);
  uint lambda = 5;
  bool *seed = allocate_seed(lambda);
  generate_seed(seed, lambda);
  seed_random_state(seed, lambda, random_state);

  generate_random_matrix(&m, random_state, field);
  print_matrix(&m);

  clear_matrix(&m);
  gmp_randclear(random_state);
  free(seed);
}
