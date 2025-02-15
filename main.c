#include "constants.h"
#include "field_arithmetics.h"
#include "key_generation.h"
#include "matrix.h"
#include "packing.h"
#include "random.h"
#include "seed_tree.h"
#include "sign.h"
#include "types.h"

#include <math.h>
#include <openssl/bio.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <stdint.h>

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

  // TODO: Apply phase_one, phase_two, phase_three
}
