#ifndef SIGN_H_
#define SIGN_H_

#include "types.h"
#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>

void allocate_commit(uchar **commit, uint lambda);
void allocate_all_commits(uchar ****commits, SignatureParameters params);
void allocate_all_party_data(PartyData ***data, SignatureParameters params);
void allocate_hash_digest(uchar **digest, uint lambda);
void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda);
void allocate_challenges(Matrix **challenges, SignatureParameters params);
void allocate_all_parties(PartyState ***parties, SignatureParameters params);
void initialize_sha3(EVP_MD_CTX **ctx, uint lambda);
int hash0(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda, uint l,
          uint i, seed_t state);
int hash0_last(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda,
               uint l, uint i, seed_t seed, Matrix alpha, Matrix K, Matrix C);
int hash1(uchar *digest, EVP_MD_CTX *context, uchar *message, uint message_size,
          uchar *salt, uint lambda, uint number_of_rounds,
          uint number_of_parties, uchar ***commits);

void phase_one(uchar ***commits, seed_t salt, PartyData **data,
               SignatureParameters params, MinRankInstance instance,
               MinRankSolution solution);
void phase_two(Matrix *challenges, uchar *h1, uchar *message, uint message_size,
               SignatureParameters params, uchar *salt, uchar ***commits);
void phase_three(Matrix *challenges, MinRankInstance instance,
                 PartyState **parties, PartyData **data,
                 SignatureParameters params);
void phase_four(uchar *h2, uint *second_challenges, uchar *message,
                uint message_size, seed_t salt, uchar *h1, PartyState **parties,
                SignatureParameters params) ;

#endif // SIGN_H_
