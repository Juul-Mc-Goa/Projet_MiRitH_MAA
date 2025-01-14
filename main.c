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

uint **identity_matrix(MatrixSize size) {
  uint **result = malloc(size.m * sizeof(uint *));
  for (uint i = 0; i < size.m; i++) {
    result[i] = malloc(size.n * sizeof(uint));
    for (uint j = 0; j < size.n; j++) {
      result[i][j] = (i == j) ? 1 : 0;
    }
  }

  return result;
}

uint **rotation_matrix(MatrixSize size) {
  uint **result = malloc(size.m * sizeof(uint *));
  for (uint i = 0; i < size.m; i++) {
    result[i] = malloc(size.n * sizeof(uint));
    for (uint j = 0; j < size.n; j++) {
      if ((j + 1) % size.n == i) {
        result[i][j] = 1;

      } else {
        result[i][j] = 0;
      }
    }
  }

  return result;
}

void test_matrix_sum() {
  printf("------------------------------ beginning matrix sum test...\n");
  MatrixSize size;
  size.m = 5;
  size.n = 5;

  uint **identity = identity_matrix(size);
  uint **rotation = rotation_matrix(size);

  Matrix m_id;
  allocate_matrix(&m_id, size);
  matrix_init_set_ui(&m_id, identity);

  printf("Id matrix:\n");
  print_matrix(&m_id);

  Matrix m_rot;
  allocate_matrix(&m_rot, size);
  matrix_init_set_ui(&m_rot, rotation);

  printf("circular permutation matrix:\n");
  print_matrix(&m_rot);

  Matrix m_sum;
  allocate_matrix(&m_sum, size);
  matrix_sum(&m_sum, m_id, m_rot);
  printf("matrix sum (Id + Rot):\n");
  print_matrix(&m_sum);

  clear_matrix(&m_id);
  clear_matrix(&m_rot);
  clear_matrix(&m_sum);
  free(identity);
  free(rotation);
}

void test_matrix_prod() {
  printf("------------------------------ beginning matrix product test...\n");
  MatrixSize size;
  size.m = 5;
  size.n = 5;

  uint **identity = identity_matrix(size);
  uint **rotation = rotation_matrix(size);

  Matrix m_id;
  allocate_matrix(&m_id, size);
  matrix_init_set_ui(&m_id, identity);

  printf("Id matrix:\n");
  print_matrix(&m_id);

  Matrix m_rot;
  allocate_matrix(&m_rot, size);
  matrix_init_set_ui(&m_rot, rotation);

  printf("circular permutation matrix:\n");
  print_matrix(&m_rot);

  Matrix m_prod;
  allocate_matrix(&m_prod, size);
  matrix_product(&m_prod, m_rot, m_rot);
  printf("matrix product (Rot * Rot):\n");
  print_matrix(&m_prod);

  clear_matrix(&m_id);
  clear_matrix(&m_rot);
  clear_matrix(&m_prod);
  free(identity);
  free(rotation);
}

void test_random_matrix() {
  printf("------------------------------ beginning random matrix test...\n");
  Matrix m;
  m.size.m = 4;
  m.size.n = 4;
  allocate_matrix(&m, m.size);

  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);
  generate_prime(m.moduli, 5);

  generate_random_matrix(&m, random_state, 4);
  print_matrix(&m);

  clear_matrix(&m);
  gmp_randclear(random_state);
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

  // prime generation
  mpz_t prime;
  mpz_init(prime);
  generate_prime(prime, 5);
  gmp_printf("generated prime: %Zd\n", prime);

  // matrices
  test_matrix_sum();
  test_matrix_prod();
  test_random_matrix();
}
