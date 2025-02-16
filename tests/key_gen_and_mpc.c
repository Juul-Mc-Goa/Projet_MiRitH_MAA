/* Generates a `PubliPrivateKeyPair`, then:
 * 1. unpack the public key to get `MinRankInstance`,
 * 2. unpack the private key to get a `MinRankSolution`,
 * 3. check that the solution is correct with `mpc_check_solution`. */
#include "../constants.h"
#include "../key_generation.h"
#include "../mpc.h"
#include "../random.h"
#include "../packing.h"

#include <gmp.h>
#include <stdio.h>

int main(int argc, char **argv) {
  printf("--------------------------------------- beginning the keygen and mpc"
         " test...\n");
  // define the signature parameters
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

  uint matrix_count = params.solution_size + 1;

  printf("parameters:\n");
  printf("- lambda: %u\n- matrix size: (%u, %u)\n- target rank: %u\n"
         "- solution size: %u\n- first challenge size: %u\n- number of "
         "parties: %u\n- tau: %u\n\n",
         params.lambda, params.matrix_dimension.m, params.matrix_dimension.n,
         params.target_rank, params.solution_size, params.first_challenge_size,
         params.number_of_parties, params.tau);

  // generate the keys
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);
  printf("\npublic key:\n- lambda: %u\n- seed: ", key_pair.public_key.lambda);
  print_seed(key_pair.public_key.seed);
  printf("\n- m0:\n");
  print_matrix(&key_pair.public_key.m0);

  printf("\nprivate key:\n- lambda: %u\n- seed: ", key_pair.private_key.lambda);
  print_seed(key_pair.private_key.seed);

  gmp_randstate_t public_random_state, private_random_state;
  gmp_randinit_default(public_random_state);
  gmp_randinit_default(private_random_state);

  seed_random_state(key_pair.private_key.seed, private_random_state);
  seed_random_state(key_pair.public_key.seed, public_random_state);

  // unpack MinRank instance
  MinRankInstance instance;
  init_instance(&instance, matrix_count, params.matrix_dimension);
  unpack_instance_from_public_key(&instance, public_random_state, params,
                                  key_pair.public_key);

  // unpack MinRank solution
  MinRankSolution solution;
  init_solution(&solution, matrix_count, params.target_rank,
                params.matrix_dimension);
  unpack_solution_from_private_key(&solution, private_random_state, params,
                                   key_pair.private_key);
  printf("\nsolution:\nalpha: ");
  print_matrix(&solution.alpha);
  printf("K:\n");
  print_matrix(&solution.K);

  // check that `solution` is correct
  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate a random matrix `R`
  Matrix R;
  MatrixSize R_size = {params.first_challenge_size, params.matrix_dimension.m};
  allocate_matrix(&R, params.field, R_size);
  generate_random_matrix(&R, random_state, params.field);

  bool result = mpc_check_solution(random_state, params.number_of_parties, R,
                                   instance, solution);
  printf("is solution correct?  %u\n", result);

  clear_key_pair(key_pair);
  clear_instance(&instance);
  clear_solution(&solution);
  clear_matrix(&R);
  gmp_randclear(random_state);
  gmp_randclear(public_random_state);
  gmp_randclear(private_random_state);
}
