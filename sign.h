#ifndef SIGN_H_
#define SIGN_H_

#include "types.h"
#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>

uint commit_size(uint lambda);
void allocate_commit(uchar **commit, uint lambda);
void allocate_all_commits(uchar ****commits, SignatureParameters params);
void allocate_all_party_seeds(uchar ****party_seeds,
                              SignatureParameters params);
void allocate_all_party_data(PartyData ***data, SignatureParameters params);
void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda);
uint allocate_signature_digest(uchar **digest, uint *second_challenges,
                               SignatureParameters params);
void allocate_first_challenge(Matrix **challenges, SignatureParameters params);
void allocate_second_challenge(uint **challenges, SignatureParameters params);
void allocate_all_parties(PartyState ***parties, SignatureParameters params);
void initialize_sha3(EVP_MD_CTX **ctx, uint lambda);

int hash0(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda, uint l,
          uint i, seed_t state);
int hash0_last(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda,
               uint l, uint i, seed_t seed, Matrix alpha, Matrix K, Matrix C);
int hash1(uchar *digest, EVP_MD_CTX *context, uchar *message, uint message_size,
          uchar *salt, uint lambda, uint number_of_rounds,
          uint number_of_parties, uchar ***commits);

void prg_first_challenge(Matrix *challenges, uchar *h1,
                         SignatureParameters params) ;
void prg_second_challenge(uint *challenges, uchar *h2,
                          SignatureParameters params);

void phase_one(uchar ***commits, uchar ***party_seeds, seed_t salt,
               PartyData **data, SignatureParameters params,
               MinRankInstance instance, MinRankSolution solution);
void phase_two(Matrix *challenges, uchar *h1, uchar *message, uint message_size,
               SignatureParameters params, uchar *salt, uchar ***commits);
void phase_three(Matrix *challenges, MinRankInstance instance,
                 PartyState **parties, PartyData **data,
                 SignatureParameters params);
void phase_four(uchar *h2, uint *second_challenges, uchar *message,
                uint message_size, seed_t salt, uchar *h1, PartyState **parties,
                SignatureParameters params);
void phase_five(uchar *digest, uint digest_size, uint *second_challenges,
                seed_t salt, uchar *h1, uchar *h2, uchar ***commits,
                uchar ***party_seeds, PartyState **parties, PartyData **data,
                SignatureParameters params);

uint sign(uchar **digest, uchar *message, uint msg_size,
          PublicPrivateKeyPair key_pair, SignatureParameters params);
#endif // SIGN_H_
