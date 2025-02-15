#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"
#include "../random.h"

#include <stdint.h>
#include <stdio.h>

int main(int argc, char **argv) {
  printf("----------------------------------------------- beginning key "
         "generation test...\n");
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

  printf("parameters:\n");
  printf("- lambda: %u\n- matrix size: (%u, %u)\n- target rank: %u\n"
         "- solution size: %u\n- first challenge size: %u\n- number of "
         "parties: %u\n- tau: %u\n\n",
         params.lambda, params.matrix_dimension.m, params.matrix_dimension.n,
         params.target_rank, params.solution_size, params.first_challenge_size,
         params.number_of_parties, params.tau);

  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);

  printf("public key:\n");
  printf("- lambda: %u\n", key_pair.public_key.lambda);
  printf("- seed: ");
  print_seed(key_pair.public_key.seed);
  printf("- matrix m0:\n");
  print_matrix(&key_pair.public_key.m0);
  printf("\n");

  printf("private key:\n");
  printf("- lambda: %u\n", key_pair.private_key.lambda);
  printf("- seed: ");
  print_seed(key_pair.private_key.seed);
  printf("\n");

  printf("finished key generation.\n");
}
