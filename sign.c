#include "constants.h"
#include "gmp.h"
#include "matrix.h"
#include "mpc.h"
#include "packing.h"
#include "random.h"
#include "seed_tree.h"
#include "types.h"

#include <math.h>
#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>

void allocate_commit(uchar *commit, uint lambda) {
  commit = malloc(sizeof(uchar) * (uint)ceil(lambda / 8.0));
}

void allocate_all_commits(uchar ***commits, SignatureParameters params) {
  uint lambda = params.lambda;
  uint rounds = params.tau;
  uint N = params.number_of_parties;

  commits = malloc(sizeof(uchar **) * rounds);
  for (uint round = 0; round < rounds; round++) {
    commits[round] = malloc(sizeof(uchar *) * N);
    for (uint party; party < N; party++) {
      allocate_commit(commits[round][party], lambda);
    }
  }
}

void allocate_all_party_data(PartyData **data, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint m = params.matrix_dimension.m;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint lambda = params.lambda;
  uint rounds = params.tau;
  uint matrix_count = params.solution_size + 1;
  size_t seed_size = (size_t)ceil(lambda / 8.0);

  MatrixSize alpha_size = {1, matrix_count};
  MatrixSize A_size = {s, r};
  MatrixSize K_size = {r, n - r};
  MatrixSize C_size = {s, n - r};

  data = malloc(sizeof(PartyData *) * rounds);

  for (uint round = 0; round < rounds; round++) {
    data[round] = malloc(sizeof(PartyData) * N);
    for (uint party; party < N; party++) {
      allocate_matrix(&data[round][party].alpha, GF_16, alpha_size);
      allocate_matrix(&data[round][party].A, GF_16, A_size);
      allocate_matrix(&data[round][party].C, GF_16, C_size);
      allocate_matrix(&data[round][party].K, GF_16, K_size);
    }
  }
}

void allocate_challenges(Matrix **challenges, SignatureParameters params) {
  uint s = params.first_challenge_size;
  uint m = params.matrix_dimension.m;
  uint rounds = params.tau;
  MatrixSize R_size = {s, m};

  *challenges = malloc(sizeof(Matrix) * rounds);

  for (uint round = 0; round < rounds; round++) {
    allocate_matrix(&(*challenges)[round], GF_16, R_size);
  }
}

void allocate_all_parties(PartyState ***parties, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint m = params.matrix_dimension.m;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint rounds = params.tau;

  *parties = malloc(sizeof(PartyState *) * rounds);

  for (uint round = 0; round < rounds; round++) {
    init_parties(*parties[round], N, s, params.matrix_dimension, r);
  }
}

void allocate_hash_digest(uchar *digest, uint lambda) {
  size_t hash_size = (size_t)ceil(lambda / 8.0);
  digest = malloc(sizeof(uchar) * hash_size);
}

void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda) {
  size_t seed_size = (size_t)ceil(lambda / 8.0);
  party_seeds = malloc(number_of_parties * sizeof(uchar *));

  for (uint party = 0; party < number_of_parties; party++) {
    party_seeds[party] = malloc(seed_size * sizeof(uchar));
  }
}

void initialize_sha3(EVP_MD_CTX *ctx, uint lambda) {
  EVP_MD *sha3 = NULL;
  ctx = EVP_MD_CTX_new();
  if (ctx == NULL)
    goto err;

  /*
   * Fetch the KECCAK algorithm implementation for doing the digest.
   * - first NULL param: use the "default" library context
   * - second NULL param: no particular criteria for the implementation
   */
  sha3 = EVP_MD_fetch(NULL, "SHA3-256", NULL);
  if (sha3 == NULL)
    goto err;

  /* Initialise the digest operation */
  if (!EVP_DigestInit_ex(ctx, sha3, NULL))
    goto err;

err:
  EVP_MD_free(sha3);
}

/* applies SHA3 to `(salt, l, i, state)` */
int hash0(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda, uint l,
          uint i, seed_t state) {
  int ret = 1;
  // allocate `msg`
  uint salt_size = lambda >> 2, int_size = 32, state_size = lambda >> 3;
  uint total_size = salt_size + 2 * int_size + state_size;
  uchar *msg = malloc(sizeof(uchar) * total_size);

  uint bit_index = 0;
  msg[0] = (uchar)0;
  // pack `salt`
  for (uint j = 0; j < 2 * lambda; j++) {
    pack_boolean(msg, salt.data[j], &bit_index);
  }
  // pack `l, i`
  pack_32bit_value(msg, l, &bit_index);
  pack_32bit_value(msg, i, &bit_index);
  // pack `state`
  for (uint j = 0; j < lambda; j++) {
    pack_boolean(msg, state.data[j], &bit_index);
  }
  // pass the message to be digested
  if (!EVP_DigestUpdate(context, msg, sizeof(msg)))
    goto err;

  // calculate the digest
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;

  ret = 0;

err:
  /* Clean up all the resources we allocated */
  free(msg);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

/* applies SHA3 to `(salt, l, i, alpha, K, C)` */
int hash0_last(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda,
               uint l, uint i, seed_t seed, Matrix alpha, Matrix K, Matrix C) {
  int ret = 1;

  // compute various sizes to obtain the length of `msg`
  uint salt_size = lambda >> 2;
  uint seed_size = lambda >> 3;
  uint int_size = 32;

  uint alpha_size = alpha.size.n;
  uint K_length = K.size.m * K.size.n;
  uint C_length = C.size.m * C.size.n;
  uint total_matrices_length = alpha_size + K_length + C_length;

  uint total_size = salt_size + 2 * int_size + seed_size;
  total_size += (total_matrices_length >> 1);
  if ((total_matrices_length & 1) != 0) {
    total_size += 1;
  }

  // allocate `msg`
  uchar *msg = malloc(sizeof(uchar) * total_size);
  pack_last_state(msg, salt, lambda, l, i, alpha, K, C);

  // pass the message to be digested
  if (!EVP_DigestUpdate(context, msg, sizeof(msg)))
    goto err;

  // calculate the digest
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;

  ret = 0;

err:
  /* Clean up all the resources we allocated */
  free(msg);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

/* Concatenate `message, salt` and each `commits[i][j]` before applying SHA3. */
int hash1(uchar *digest, EVP_MD_CTX *context, uchar *message, uint message_size,
          uchar *salt, uint lambda, uint number_of_rounds,
          uint number_of_parties, uchar ***commits) {
  int ret = 1;
  // `salt` has 2 * lambda booleans, packed into lambda / 4 `uchar`s
  uint salt_size = lambda >> 2;
  uchar *msg = malloc(sizeof(uchar) * salt_size);

  // pass `message` to the context
  if (!EVP_DigestUpdate(context, message, sizeof(message)))
    goto err;

  // pack `salt` into `msg`
  uint bit_index = 0;
  for (uint i = 0; i < message_size; i++) {
    pack_boolean(msg, salt[i], &bit_index);
  }

  // pass `msg` to the context
  if (!EVP_DigestUpdate(context, msg, sizeof(msg)))
    goto err;

  // pass `commit` to the context
  for (uint i = 0; i < number_of_rounds; i++) {
    for (uint j = 0; j < number_of_parties; j++) {
      if (!EVP_DigestUpdate(context, commits[i][j], sizeof(commits[i][j])))
        goto err;
    }
  }

  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;

err:
  free(msg);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

void hash2(uchar *msg, uint msg_size) {
  // TODO
}

/* Phase one: generate `params.tau` instances of the MPC protocol. The
 * result is stored in `commits`: a 2D array of SHA3 digests. */
void phase_one(uchar ***commits, seed_t salt, PartyData **data,
               SignatureParameters params, MinRankInstance instance,
               MinRankSolution solution) {
  // unpack `params`
  uint k = params.solution_size;
  uint N = params.number_of_parties;
  uint m = params.matrix_dimension.m;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint lambda = params.lambda;
  uint tau = params.tau;
  uint matrix_count = params.solution_size + 1;
  size_t seed_size = (size_t)ceil(lambda / 8.0);

  seed_t round_seed, party_seed;
  allocate_seed(&round_seed, lambda);
  allocate_seed(&party_seed, lambda);

  // initialize context for hashes (used in `commits`)
  EVP_MD_CTX *hash_context;
  EVP_MD *sha3;
  initialize_sha3(hash_context, lambda);

  // create various matrices
  Matrix alpha_sum, A_sum, C, C_sum, K_sum;
  MatrixSize A_size = {s, r};
  MatrixSize C_size = {s, r};

  Matrix alpha = solution.alpha;
  Matrix K = solution.K;

  allocate_matrix(&alpha_sum, GF_16, alpha.size);
  allocate_matrix(&A_sum, GF_16, A_size);
  allocate_matrix(&K_sum, GF_16, K.size);
  allocate_matrix(&C, GF_16, C_size);
  allocate_matrix(&C_sum, GF_16, C_size);

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  Matrix A;
  allocate_matrix(&A, GF_16, A_size);
  generate_random_matrix(&A, prg_state, GF_16);

  uchar **party_seeds = malloc(sizeof(uchar *) * N);
  for (uint i = 0; i < N; i++) {
    party_seeds[i] = malloc(sizeof(uchar) * seed_size);
  }

  for (uint round = 0; round < params.tau; round++) {
    fill_matrix_with_zero(&A_sum);
    fill_matrix_with_zero(&alpha_sum);
    fill_matrix_with_zero(&K_sum);
    fill_matrix_with_zero(&C_sum);

    generate_seed(round_seed);
    TreePRG(&salt, &round_seed, (size_t)lambda, N, party_seeds);

    for (uint party = 0; party < N - 1; party++) {
      party_seed.data = party_seeds[party];
      PRG_init(&salt, &party_seed, lambda, &prg_state);
      generate_random_matrix(&data[round][party].A, prg_state, GF_16);
      generate_random_matrix(&data[round][party].alpha, prg_state, GF_16);
      generate_random_matrix(&data[round][party].C, prg_state, GF_16);
      generate_random_matrix(&data[round][party].K, prg_state, GF_16);

      // update `{matrix}_sum`
      matrix_sum(&alpha_sum, alpha_sum, alpha);
      matrix_sum(&A_sum, A_sum, A);
      matrix_sum(&K_sum, K_sum, K);
      matrix_sum(&C_sum, C_sum, C);

      hash0(commits[round][party], hash_context, salt, lambda, round, party,
            party_seed);
    }
    // last party is special
    party_seed.data = party_seeds[N - 1];
    PRG_init(&salt, &party_seed, lambda, &prg_state);
    generate_random_matrix(&data[round][N - 1].A, prg_state, GF_16);
    generate_random_matrix(&data[round][N - 1].alpha, prg_state, GF_16);
    generate_random_matrix(&data[round][N - 1].C, prg_state, GF_16);
    generate_random_matrix(&data[round][N - 1].K, prg_state, GF_16);
    // compute last alpha, K, C
    matrix_sum(&data[round][N - 1].alpha, solution.alpha, alpha_sum);
    matrix_sum(&data[round][N - 1].K, solution.K, K_sum);
    matrix_product(&C, data[round][N - 1].A, data[round][N - 1].K);
    matrix_sum(&data[round][N - 1].C, C, C_sum);
    hash0_last(commits[round][N - 1], hash_context, salt, lambda, round, N - 1,
               party_seed, alpha, K, C);
  }

  clear_seed(&round_seed);
  clear_seed(&party_seed);
  clear_matrix(&alpha_sum);
  clear_matrix(&A_sum);
  clear_matrix(&C);
  clear_matrix(&C_sum);
  clear_matrix(&K_sum);
  EVP_MD_CTX_free(hash_context);
  gmp_randclear(prg_state);
}

/* Compute a hash of `message || salt || commits`, then use it to create
 * the challenges for each round. */
void phase_two(Matrix *challenges, uchar *h1, uchar *message, uint message_size,
               SignatureParameters params, uchar *salt, uchar ***commits) {
  uint lambda = params.lambda;
  uint tau = params.tau;

  seed_t salt_wrap;
  salt_wrap.size = (size_t)ceil(lambda / 4.0);
  salt_wrap.data = salt;

  // compute the hash
  EVP_MD_CTX *context = NULL;
  initialize_sha3(context, lambda);

  hash1(h1, context, message, message_size, salt, lambda, tau, tau, commits);

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  // create the challenges
  for (uint round = 0; round < tau; round++) {
    PRG_init(&salt_wrap, NULL, lambda, &prg_state);
    generate_random_matrix(&challenges[round], prg_state, GF_16);
  }

  gmp_randclear(prg_state);
}

void phase_three(MinRankInstance instance, PartyState **parties,
                 PartyData **data, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint s = params.first_challenge_size;
  uint r = params.target_rank; // rank of the solution
  uint matrix_count = instance.matrix_count;

  MatrixSize size = params.matrix_dimension;
  uint n = size.n;
  MatrixSize C_size = {s, n - r};
  Matrix S;
  allocate_matrix(&S, GF_16, C_size);

  for (uint round = 0; round < params.tau; round++) {
    for (uint party = 0; party < N; party++) {
      compute_local_m(&parties[round][party], instance,
                      data[round][party].alpha, matrix_count);
      compute_local_s(&parties[party], data[round][parties].R,
                      data[round][parties].A);
    }
    compute_global_s(&S, parties[round], N);
    for (uint party = 0; party < N; party++) {
      compute_local_v(&parties[round][party], S, data[round][party].K,
                      data[round][party].R, data[round][party].C);
    }
  }
}

void phase_four() {}

void sign(uchar *digest, MinRankInstance instance, MinRankSolution solution,
          uchar *message, uint msg_size, SignatureParameters params) {
  uchar ***commits;
  allocate_all_commits(commits, params);

  PartyData **data;
  allocate_all_party_data(data, params);

  // generate salt
  seed_t salt;
  allocate_seed(&salt, 2 * params.lambda);
  generate_seed(salt);

  phase_one(commits, salt, data, params, instance, solution);

  Matrix *challenges;
  allocate_challenges(&challenges, params);
  PartyState **parties;
  allocate_all_parties(&parties, params);
  uchar *h1;
  allocate_hash_digest(h1, params.lambda);

  phase_two(challenges, h1, message, msg_size, params, salt.data, commits);
  phase_three(instance, parties, data, params);
  /* phase_four(message, salt, h1); */
}

/* Phase 3 of signing
 *
 * Given as input:
 * - the public key 'M';
 * - the secret key 'E';
 * - the matrix 'A';
 * - the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr';
 * computes and hashes the additive shares 'S_shr', 'V_shr'.
 * The matrix 'S_shr[i]' is written over 'A_shr'. */
/* void sign_phase3(uint8_t M[K + 1][matrix_bytes_size(M, N)], */
/*                  uint8_t E[matrix_bytes_size(M, N)], */
/*                  uint8_t A[matrix_bytes_size(S, R)], */
/*                  uint8_t A_shr[N_PARTIES][matrix_bytes_size(S, R)], */
/*                  uint8_t C_shr[N_PARTIES][matrix_bytes_size(S, N - R)], */
/*                  uint8_t a_shr[N_PARTIES][matrix_bytes_size(K, 1)], */
/*                  uint8_t K_shr[N_PARTIES][matrix_bytes_size(R, N - R)], */
/*                  uint8_t R[matrix_bytes_size(S, M)]) { */
/*   // MPC (Multi-Party Computation) implementation goes here */
/* } */

/* Phase 4 of signing
 *
 * Given as input:
 * - the hash hash2;
 * computes the challenge i_star */
/* void sign_phase4(uint32_t i_star[TAU], const hash_t hash2) { */
/*   size_t l; */
/*   prng_t prng; */

/*   // Initialize PRNG from 'hash2' */
/*   prng_init(&prng, hash2, NULL); */

/*   for (l = 0; l < TAU; l++) { */
/*     uint32_t r; */
/*     // Generate random number and compute i_star[l] */
/*     prng_bytes(&prng, &r, sizeof(r)); */
/*     i_star[l] = r % N_PARTIES; */
/*   } */
/* } */
