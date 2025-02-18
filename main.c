#include "constants.h"
#include "key_generation.h"
#include "random.h"
#include "sign.h"
#include "types.h"

#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void print_digest(uchar *digest, uint digest_size) {
  for (uint i = 0; i < digest_size; i++) {
    printf("%02x", (uint)digest[i]);
  }
  printf("\n");
}
int main(int argc, char **argv) {
  /* SignatureParameters params = PARAMS_1_A_SHORT; */
  SignatureParameters params = PARAMS_1_A_FAST;

  uchar *message = (uchar *)"test a b c";
  uint msg_size = strlen(message);
  printf("message: %s, length: %u\n", message, msg_size);

  clock_t start_time = clock();
  // generate the keys
  PublicPrivateKeyPair key_pair;
  allocate_key_pair(&key_pair, params);
  key_gen(&key_pair, params);

  uchar *digest;
  uint digest_size = sign(&digest, message, msg_size, key_pair, params);

  clock_t end_time = clock();
  printf("---------------------------------------------------------------------"
         "-----------\n");
  print_digest(digest, digest_size);
  double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  printf("\ndigest size: %u\n", digest_size);
  printf("time: %lf s\n", total_time);
  printf("clock cycles: %ld\n", end_time - start_time);
}
