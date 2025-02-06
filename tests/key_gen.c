#include "../constants.h"
#include "../key_generation.h"
#include "../matrix.h"

#include <stdio.h>

void print_seed(bool *seed, uint lambda) {
  for (uint i = 0; i < lambda; i++) {
    printf("%u", seed[i]);
  }
  printf("\n");
}

int main(int argc, char **argv) {
  printf("----------------------------------------------- beginning key generation test...\n");
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

  PublicPrivateKeyPair key_pair = key_gen(params);

  printf("public key:\n");
  printf("- lambda: %u\n", key_pair.public_key.lambda);
  printf("- seed: ");
  print_seed(key_pair.public_key.seed, key_pair.public_key.lambda);
  printf("- matrix m0:\n");
  print_matrix(&key_pair.public_key.m0);
  printf("\n");

  printf("private key:\n");
  printf("- lambda: %u\n", key_pair.private_key.lambda);
  printf("- seed: ");
  print_seed(key_pair.private_key.seed, key_pair.private_key.lambda);
  printf("\n");

  printf("finished key generation.\n");
}
