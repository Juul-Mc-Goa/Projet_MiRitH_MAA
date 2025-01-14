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

int main(int argc, char **argv) {
  uint lambda = 4;
  bool *seed = allocate_seed(lambda);
  generate_seed(seed, lambda);

  printf("generated seed:\n");
  for (uint i = 0; i < 4; i++) {
    printf("%u: %u\n", i, seed[i]);
  }

  MatrixSize size;
  size.m = 3;
  size.n = 3;
  uint **matrix = malloc(size.m * sizeof(uint *));
  matrix[0] = malloc(size.n * sizeof(uint));
  matrix[1] = malloc(size.n * sizeof(uint));
  matrix[2] = malloc(size.n * sizeof(uint));

  matrix[0][0] = 1;
  matrix[0][1] = 0;
  matrix[0][2] = 0;

  matrix[1][0] = 0;
  matrix[1][1] = 1;
  matrix[1][2] = 0;

  matrix[2][0] = 0;
  matrix[2][1] = 0;
  matrix[2][2] = 1;

  Matrix m;
  allocate_matrix(&m, size);
  matrix_set_ui(&m, matrix);

  print_matrix(&m);
  clear_matrix(&m);
  free(matrix);
}
