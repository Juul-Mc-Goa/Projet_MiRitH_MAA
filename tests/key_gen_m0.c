/* Perform the same computations as in `key_gen` and check that we get the same
 * result. */
#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"
#include "../random.h"

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("--------------------------------------------------- beginning key "
         "gen m0 test...\n");
  SignatureParameters params;
  params.lambda = 4;
  params.matrix_dimension.m = 3;
  params.matrix_dimension.n = 3;
  params.field = GF_16;
  params.target_rank = 1;
  params.solution_size = 4;
  params.first_challenge_size = 2;
  params.number_of_parties = 2;
  params.tau = 2;

  // call `key_gen`
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);
  printf("finished key_gen.\n");

  // perform the same computations
  gmp_randstate_t public_random_state;
  gmp_randstate_t private_random_state;

  // public seed
  gmp_randinit_default(public_random_state);
  seed_random_state(key_pair.public_key.seed, params.lambda,
                    public_random_state);
  // private seed
  gmp_randinit_default(private_random_state);
  seed_random_state(key_pair.private_key.seed, params.lambda,
                    private_random_state);
  uint *alpha = malloc(params.solution_size * sizeof(uint));
  alpha[0] = 1;

  Matrix sum;
  allocate_matrix(&sum, params.field, params.matrix_dimension);
  fill_matrix_with_zero(&sum);

  Matrix m_i;
  allocate_matrix(&m_i, params.field, params.matrix_dimension);
  for (uint i = 1; i <= params.solution_size; i++) {
    // generate M_i
    generate_random_matrix(&m_i, public_random_state, params.field);

    // generate alpha_i
    alpha[i] = generate_random_element(private_random_state, params.field);

    // compute sum += alpha_i * M_i
    scalar_product(&m_i, alpha[i], m_i);
    matrix_sum(&sum, sum, m_i);
  }

  // generate a random matrix K
  Matrix K;
  MatrixSize K_size = {params.target_rank,
                       params.matrix_dimension.n - params.target_rank};
  allocate_matrix(&K, params.field, K_size);
  generate_random_matrix(&K, private_random_state, params.field);

  // 2.3.3 compute E from K and E_R
  Matrix E, E_L, E_R;
  allocate_matrix(&E, params.field, params.matrix_dimension);
  fill_matrix_with_zero(&E);

  E_R.data = malloc(params.matrix_dimension.m * sizeof(uint *));
  split_matrix(&E_L, &E_R, E, params.target_rank);

  // generate a random matrix E_R
  generate_random_matrix(&E_R, private_random_state, params.field);

  // compute E_L = E_R * K
  matrix_product(&E_L, E_R, K);

  // 2.3.4 compute M_0 = E - sum
  Matrix m0;
  allocate_matrix(&m0, GF_16, key_pair.public_key.m0.size);

  matrix_opposite(&sum);
  matrix_sum(&m0, E, sum);

  printf("key_pair m0:\n");
  print_matrix(&key_pair.public_key.m0);
  printf("\n");

  printf("local m0:\n");
  print_matrix(&m0);
  printf("\n");

  free(E_R.data);
  clear_matrix(&sum);
  clear_matrix(&m_i);
  clear_matrix(&K);
  clear_matrix(&E);
  clear_matrix(&m0);
  gmp_randclear(public_random_state);
  gmp_randclear(private_random_state);
  clear_key_pair(key_pair);
}
