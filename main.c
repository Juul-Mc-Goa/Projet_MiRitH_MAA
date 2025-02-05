#include "constants.h"
#include "field_arithmetics.h"
#include "key_generation.h"
#include "matrix.h"
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

void test_key_gen() {
  printf("------------------------------ beginning key generation test...\n");
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

  PublicPrivateKeyPair key_pair = key_gen(params);
  printf("finished key generation.\n");
}

void test_random_matrix() {
  printf("------------------------------ beginning random matrix test...\n");
  Matrix m;
  m.size.m = 4;
  m.size.n = 4;
  FiniteField field = GF_16;
  allocate_matrix(&m, field, m.size);

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);
  uint lambda = 5;
  bool *seed = allocate_seed(lambda);
  generate_seed(seed, lambda);
  seed_random_state(seed, lambda, random_state);

  generate_random_matrix(&m, random_state, field);
  print_matrix(&m);

  clear_matrix(&m);
  gmp_randclear(random_state);
  free(seed);
}

int main(int argc, char **argv) {
  // seed generation
  uint lambda = 4;
  bool *seed = allocate_seed(lambda);
  generate_seed(seed, lambda);

  printf("generated seed:\n");
  for (uint i = 0; i < 4; i++) {
    printf("%u: %u\n", i, seed[i]);
  }

  // matrices
  test_matrix_sum();
  test_matrix_prod();
  test_random_matrix();

  // key generation
  test_key_gen();

  // print addition table in GF(16)
  printf("\n\nAddition Table for GF(16):\n");
  print_gf_16_addition_table();

  // print multiplication table in GF(16)
  printf("\n\nMultiplication Table for GF(16):\n");
  print_gf_16_mul_table();
}
