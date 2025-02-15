#include "constants.h"
#include "field_arithmetics.h"
#include "key_generation.h"
#include "matrix.h"
#include "packing.h"
#include "random.h"
#include "sign.h"
#include "types.h"

#include <math.h>
#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <stdint.h>

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

  uint seed_size = (uint)ceil(params.lambda / 8.0);

  uchar ***commits = malloc(sizeof(uchar **) * params.tau);
  for (uint round = 0; round < params.tau; round++) {
    commits[round] = malloc(sizeof(uchar *) * params.number_of_parties);
    for (uint party = 0; party < params.number_of_parties - 1; party++) {
      commits[round][party] = malloc(sizeof(uchar) * seed_size);
    }
  }

  // generate the keys
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);
  MinRankInstance instance;
  generate_random_instance(&instance, params.solution_size, prg_state, GF_16);

  unpack_instance_from_public_key(&instance, prg_state, params,
                                  key_pair.public_key);

  MinRankSolution solution;
  unpack_solution_from_private_key(&solution, prg_state, params,
                                   key_pair.private_key);

  // TODO: Apply phase_one, phase_two, phase_three, phase_four
}
