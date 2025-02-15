#include "constants.h"
#include "field_arithmetics.h"
#include "key_generation.h"
#include "matrix.h"
#include "random.h"
#include "seed_tree.h"

#include <gmp.h>
#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/evp.h>

int test_keccak_digest() {
  EVP_MD_CTX *ctx = NULL; // Message Digest context
  EVP_MD *keccak = NULL;
  const unsigned char msg[] = {0x00, 0x01, 0x02, 0x03};
  unsigned int len = 0;
  unsigned char *outdigest = NULL;
  int ret = 1;

  /* Create a context for the digest operation */
  ctx = EVP_MD_CTX_new();
  if (ctx == NULL)
    goto err;

  /*
   * Fetch the KECCAK algorithm implementation for doing the digest.
   * - first NULL param: use the "default" library context
   * - second NULL param: no particular criteria for the implementation
   */
  keccak = EVP_MD_fetch(NULL, "KECCAK-224", NULL);
  if (keccak == NULL)
    goto err;

  /* Initialise the digest operation */
  if (!EVP_DigestInit_ex(ctx, keccak, NULL))
    goto err;

  /*
   * Pass the message to be digested. This can be passed in over multiple
   * EVP_DigestUpdate calls if necessary
   */
  if (!EVP_DigestUpdate(ctx, msg, sizeof(msg)))
    goto err;

  /* Allocate the output buffer */
  outdigest = OPENSSL_malloc(EVP_MD_get_size(keccak));
  if (outdigest == NULL)
    goto err;

  /* Now calculate the digest itself */
  if (!EVP_DigestFinal_ex(ctx, outdigest, &len))
    goto err;

  /* Print out the digest result */
  BIO_dump_fp(stdout, outdigest, len);

  ret = 0;

err:
  /* Clean up all the resources we allocated */
  OPENSSL_free(outdigest);
  EVP_MD_free(keccak);
  EVP_MD_CTX_free(ctx);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

int main(int argc, char **argv) {
  uint lambda = 4;
  // seed generation
  /* bool *seed = allocate_seed(lambda); */
  /* generate_seed(seed, lambda); */

  /* printf("generated seed:\n"); */
  /* for (uint i = 0; i < 4; i++) { */
  /*   printf("%u: %u\n", i, seed[i]); */
  /* } */

  /* // print addition table in GF(16) */
  /* printf("\n\nAddition Table for GF(16):\n"); */
  /* print_gf_16_addition_table(); */

  /* // print multiplication table in GF(16) */
  /* printf("\n\nMultiplication Table for GF(16):\n"); */
  /* print_gf_16_mul_table(); */

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);
  uchar *str_salt, *str_seed;

  allocate_uchar_seed(&str_salt, 2 * lambda);
  generate_uchar_seed(str_salt, 2 * lambda);
  printf("generated str_salt\n");

  allocate_uchar_seed(&str_seed, lambda);
  generate_uchar_seed(str_seed, lambda);
  printf("generated str_seed\n");

  uint n = 7;
  uchar **output_seeds = malloc(sizeof(uchar *) * n);
  for (uint i = 0; i < n; i++) {
    output_seeds[i] = (uchar *)malloc(lambda * sizeof(uchar));
  }

  TreePRG(str_salt, str_seed, lambda, n, output_seeds);

  for (uint i = 0; i < 7; i++) {
    for (uint j = 0; j < lambda; j++) {
      uint8_t first_half = (uint8_t)output_seeds[i][j] >> 4;
      uint8_t second_half = (uint8_t)output_seeds[i][j] & 15;
      printf("%c%c", HEX_CHAR_TABLE[first_half], HEX_CHAR_TABLE[second_half]);
    }
    printf("\n");
  }
}
