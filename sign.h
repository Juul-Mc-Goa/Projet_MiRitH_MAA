#ifndef SIGN_H_
#define SIGN_H_

#include "types.h"
#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>

void allocate_commit(uchar *commit, uint lambda);
void allocate_hash_digest(uchar *digest, uint lambda);
void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda);

void initialize_sha3(EVP_MD_CTX *ctx, uint lambda);
int hash0(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda, uint l,
          uint i, seed_t state);
int hash0_last(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda,
               uint l, uint i, seed_t seed, Matrix alpha, Matrix K, Matrix C);
int hash1(uchar *digest, EVP_MD_CTX *context, uchar *message, uint message_size,
          uchar *salt, uint lambda, uint number_of_rounds,
          uint number_of_parties, uchar ***commits);

int phase_one(uchar ***commits, seed_t *seeds,SignatureParameters params,
              MinRankInstance instance, MinRankSolution solution);
#endif // SIGN_H_
