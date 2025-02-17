#include "verif.h"
#include "random.h"
#include "sign.h"
#include "types.h"

#include <string.h>

void verif(PublicKey key, uchar *message, uchar *digest, uint digest_size,
           SignatureParameters params) {
  // unpack `params`
  uint N = params.number_of_parties;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint lambda = params.lambda;
  uint hash_size = commit_size(lambda);

  uint current_index = 0;
  uint bit_index;

  seed_t salt;
  allocate_seed(&salt, 2 * params.lambda);
  generate_seed(salt);
  memcpy(salt.data, digest, salt.size);
  current_index += salt.size;

  uchar *h1 = malloc(sizeof(uchar) * hash_size);
  uchar *h2 = malloc(sizeof(uchar) * hash_size);

  memcpy(h1, &digest[current_index], hash_size);
  memcpy(h2, &digest[current_index + hash_size], hash_size);
  current_index += 2 * hash_size;

  // unpack challenges
  Matrix *first_challenges;
  allocate_first_challenge(&first_challenges, params);
  prg_first_challenge(first_challenges, h1, params);

  uint *second_challenges;
  allocate_second_challenge(&second_challenges, params);
  prg_second_challenge(second_challenges, h2, params);
}
