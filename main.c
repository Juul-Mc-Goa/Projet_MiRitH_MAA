#include "constants.h"
#include "key_generation.h"
#include "mpc.h"
#include "packing.h"
#include "random.h"
#include "sign.h"
#include "types.h"

#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <stdio.h>

int main(int argc, char **argv) {
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

  uchar *message = (uchar *)"hello";
  uint msg_size = 5;
  uchar ***commits;
  allocate_all_commits(&commits, params);

  // generate the keys
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  MinRankInstance instance;
  init_instance(&instance, params.solution_size + 1, params.matrix_dimension);
  unpack_instance_from_public_key(&instance, prg_state, params,
                                  key_pair.public_key);

  MinRankSolution solution;
  init_solution(&solution, params.solution_size + 1, params.target_rank,
                params.matrix_dimension);
  unpack_solution_from_private_key(&solution, prg_state, params,
                                   key_pair.private_key);

  PartyData **data;
  allocate_all_party_data(&data, params);

  // generate salt
  seed_t salt;
  allocate_seed(&salt, 2 * params.lambda);
  generate_seed(salt);

  phase_one(commits, salt, data, params, instance, solution);

  Matrix *challenges;
  allocate_challenges(&challenges, params);

  PartyState **parties;
  allocate_all_parties(&parties, params);
  uchar *h1;
  allocate_hash_digest(h1, params.lambda);

  printf("allocated challenges, h1, parties\n");

  phase_two(challenges, h1, message, msg_size, params, salt.data, commits);
  printf("completed phase two\n");
  phase_three(challenges, instance, parties, data, params);
  printf("completed phase three\n");
  // TODO: Apply phase_two, phase_three, phase_four
}
