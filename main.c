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

void print_digest(uchar *digest, uint digest_size) {
  printf("digest size: %u\ndigest:\n", digest_size);
  for (uint i = 0; i < digest_size; i++) {
    printf("%02x", (uint)digest[i]);
  }
  printf("\n");
}
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

  // generate the keys
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);

  uchar *digest1, *digest2;
  uint digest_size1 = sign(&digest1, message, msg_size, key_pair, params);
  print_digest(digest1, digest_size1);

  printf("test\n");
  uint digest_size2 = sign(&digest2, message, msg_size, key_pair, params);
  print_digest(digest2, digest_size2);
}
