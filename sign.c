#include "constants.h"
#include "gmp.h"
#include "matrix.h"
#include "mpc.h"
#include "openssl/cryptoerr_legacy.h"
#include "packing.h"
#include "random.h"
#include "seed_tree.h"
#include "types.h"

#include <math.h>
#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>
#include <stdio.h>
#include <string.h>

// TODO: actually implement the dependance on lambda
uint commit_size(uint lambda) { return EVP_MD_get_size(EVP_sha3_256()); }

void allocate_commit(uchar **commit, uint lambda) {
  *commit = OPENSSL_malloc(commit_size(lambda));
}

void clear_commit(uchar *commit) { OPENSSL_free(commit); }

void allocate_all_commits(uchar ****commits, SignatureParameters params) {
  uint lambda = params.lambda;
  uint rounds = params.tau;
  uint N = params.number_of_parties;

  *commits = malloc(sizeof(uchar **) * rounds);
  uchar ***commits_value = *commits;
  for (uint round = 0; round < rounds; round++) {
    commits_value[round] = malloc(sizeof(uchar *) * N);
    for (uint party = 0; party < N; party++) {
      allocate_commit(&commits_value[round][party], lambda);
    }
  }
}

void clear_all_commits(uchar ***commits, SignatureParameters params) {
  uint rounds = params.tau;
  uint N = params.number_of_parties;

  for (uint round = 0; round < rounds; round++) {
    for (uint party = 0; party < N; party++) {
      clear_commit(commits[round][party]);
    }
    free(commits[round]);
  }
  free(commits);
}

void allocate_all_party_seeds(uchar ****party_seeds,
                              SignatureParameters params) {
  uint N = params.number_of_parties;
  uint tau = params.tau;
  uint lambda = params.lambda;
  size_t seed_size = (size_t)ceil(lambda / 8.0);

  *party_seeds = malloc(sizeof(uchar ***) * tau);
  uchar ***seeds_val = *party_seeds;
  for (uint round = 0; round < tau; round++) {
    seeds_val[round] = malloc(sizeof(uchar *) * N);
    for (uint party = 0; party < N; party++) {
      seeds_val[round][party] = malloc(sizeof(uchar) * seed_size);
    }
  }
}

void clear_all_party_seeds(uchar ***party_seeds, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint tau = params.tau;
  uint lambda = params.lambda;
  size_t seed_size = (size_t)ceil(lambda / 8.0);

  for (uint round = 0; round < tau; round++) {
    party_seeds[round] = malloc(sizeof(uchar *) * N);
    for (uint party = 0; party < N; party++) {
      party_seeds[round][party] = malloc(sizeof(uchar) * seed_size);
      free(party_seeds[round][party]);
    }
    free(party_seeds[round]);
  }
  free(party_seeds);
}

void allocate_all_party_data(PartyData ***data, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint rounds = params.tau;
  uint matrix_count = params.solution_size + 1;

  MatrixSize alpha_size = {1, matrix_count};
  MatrixSize A_size = {s, r};
  MatrixSize K_size = {r, n - r};
  MatrixSize C_size = {s, n - r};

  *data = malloc(sizeof(PartyData *) * rounds);
  PartyData **data_value = *data;

  for (uint round = 0; round < rounds; round++) {
    data_value[round] = malloc(sizeof(PartyData) * N);
    for (uint party = 0; party < N; party++) {
      allocate_matrix(&data_value[round][party].alpha, GF_16, alpha_size);
      allocate_matrix(&data_value[round][party].A, GF_16, A_size);
      allocate_matrix(&data_value[round][party].C, GF_16, C_size);
      allocate_matrix(&data_value[round][party].K, GF_16, K_size);
    }
  }
}

void clear_all_party_data(PartyData **data, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint rounds = params.tau;

  for (uint round = 0; round < rounds; round++) {
    for (uint party = 0; party < N; party++) {
      clear_matrix(&data[round][party].alpha);
      clear_matrix(&data[round][party].A);
      clear_matrix(&data[round][party].C);
      clear_matrix(&data[round][party].K);
    }
    free(data[round]);
  }
}

void allocate_first_challenge(Matrix **challenges, SignatureParameters params) {
  uint s = params.first_challenge_size;
  uint m = params.matrix_dimension.m;
  uint rounds = params.tau;
  MatrixSize R_size = {s, m};

  *challenges = malloc(sizeof(Matrix) * rounds);

  for (uint round = 0; round < rounds; round++) {
    allocate_matrix(&(*challenges)[round], GF_16, R_size);
  }
}

void clear_first_challenge(Matrix *challenges, SignatureParameters params) {
  uint rounds = params.tau;
  for (uint round = 0; round < rounds; round++) {
    clear_matrix(&challenges[round]);
  }
}

void allocate_second_challenge(uint **challenges, SignatureParameters params) {
  *challenges = malloc(sizeof(uint) * params.tau);
}

void allocate_all_parties(PartyState ***parties, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint rounds = params.tau;

  *parties = malloc(sizeof(PartyState *) * rounds);
  PartyState **parties_value = *parties;

  for (uint round = 0; round < rounds; round++) {
    parties_value[round] = malloc(sizeof(PartyState) * N);
    init_parties(parties_value[round], N, s, params.matrix_dimension, r);
  }
}

void clear_all_parties(PartyState **parties, SignatureParameters params) {
  uint N = params.number_of_parties;
  uint rounds = params.tau;

  for (uint round = 0; round < rounds; round++) {
    clear_parties(parties[round], N);
  }
  free(parties);
}

void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda) {
  size_t seed_size = (size_t)ceil(lambda / 8.0);
  party_seeds = malloc(number_of_parties * sizeof(uchar *));

  for (uint party = 0; party < number_of_parties; party++) {
    party_seeds[party] = malloc(seed_size * sizeof(uchar));
  }
}

uint allocate_signature_digest(uchar **digest, uint *second_challenges,
                               SignatureParameters params) {
  uint lambda = params.lambda;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint tau = params.tau;
  uint N = params.number_of_parties;

  uint seed_size = (uint)ceil(lambda / 8.0);
  uint salt_size = seed_size << 1;
  uint hash_size = commit_size(lambda);
  uint alpha_length = params.solution_size + 1;
  uint K_length = r * (n - r);
  uint C_length = s * (n - r);
  uint S_length = s * r;

  uint total_size = salt_size + 2 * hash_size;

  for (uint round = 0; round < tau; round++) {
    // party states
    total_size += (N - 1) * seed_size;
    // last party has a different state
    uint i_star = second_challenges[round] % N;
    if (i_star != N)
      total_size += (uint)ceil((alpha_length + K_length + C_length) / 2.0);

    // commit for i_star
    total_size += hash_size;

    // S matrix for i_star
    total_size += (uint)ceil(S_length / 2.0);
  }

  // one null byte to end the digest
  total_size += 1;

  *digest = malloc(sizeof(uchar) * total_size);
  return total_size;
}

void initialize_sha3(EVP_MD_CTX **ctx, uint lambda) {
  *ctx = EVP_MD_CTX_new();
  EVP_MD_CTX *ctx_value = *ctx;
  if (ctx_value == NULL)
    goto err;

  /* Initialise the digest operation */
  if (!EVP_DigestInit_ex(ctx_value, EVP_sha3_256(), NULL))
    goto err;

err:
  ERR_print_errors_fp(stderr);
}

/* applies SHA3 to `(salt, l, i, state)` */
int hash0(uchar *digest, EVP_MD_CTX *context, seed_t salt, uint lambda, uint l,
          uint i, seed_t state) {
  int ret = 1;
  // allocate `msg`
  uint int_size = 4;
  uint total_size = salt.size + 2 * int_size + state.size;
  uchar *msg = malloc(sizeof(uchar) * total_size);

  // pack `salt`
  memcpy(msg, salt.data, salt.size);
  // pack `l, i`
  msg[salt.size] = (uchar)((l >> 24) & 255);
  msg[salt.size + 1] = (uchar)((l >> 16) & 255);
  msg[salt.size + 2] = (uchar)((l >> 8) & 255);
  msg[salt.size + 3] = (uchar)(l & 255);

  msg[salt.size + 4] = (uchar)((i >> 24) & 255);
  msg[salt.size + 5] = (uchar)((i >> 16) & 255);
  msg[salt.size + 6] = (uchar)((i >> 8) & 255);
  msg[salt.size + 7] = (uchar)(i & 255);
  // pack `state`
  memcpy(&msg[salt.size + 8], state.data, state.size);

  // pass the message to be digested
  if (!EVP_DigestUpdate(context, msg, total_size))
    goto err;

  // calculate the digest
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;

  // reinit the hash context
  if (!EVP_DigestInit_ex2(context, EVP_sha3_256(), NULL))
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
  uint salt_size = salt.size;
  uint seed_size = seed.size;
  uint int_size = 4;

  uint alpha_size = alpha.size.n;
  uint K_length = K.size.m * K.size.n;
  uint C_length = C.size.m * C.size.n;
  uint total_matrices_length = alpha_size + K_length + C_length;

  uint total_size = salt_size + 2 * int_size + seed_size;
  total_size += ceil(total_matrices_length / 2.0);

  // allocate `msg`
  uchar *msg = malloc(sizeof(uchar) * total_size);
  pack_last_state(msg, salt, lambda, l, i, alpha, K, C);

  // pass the message to be digested
  if (!EVP_DigestUpdate(context, msg, total_size))
    goto err;

  // calculate the digest
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;

  // reinit the hash context
  if (!EVP_DigestInit_ex2(context, EVP_sha3_256(), NULL))
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
  uint salt_size = ceil(lambda / 4.0);

  // pass `message` to the context
  if (!EVP_DigestUpdate(context, message, message_size))
    goto err;
  // pass `msg` to the context
  if (!EVP_DigestUpdate(context, salt, salt_size))
    goto err;
  //
  // pass `commit` to the context
  for (uint i = 0; i < number_of_rounds; i++) {
    for (uint j = 0; j < number_of_parties; j++) {
      if (!EVP_DigestUpdate(context, commits[i][j], commit_size(lambda)))
        goto err;
    }
  }
  // compute the hash
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;
  // reinit the hash context
  if (!EVP_DigestInit_ex2(context, EVP_sha3_256(), NULL))
    goto err;

  ret = 0;

err:
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

void hash2(uchar *digest, uchar *message, uint message_size, seed_t salt,
           uchar *h1, PartyState **parties, EVP_MD_CTX *context,
           SignatureParameters params) {
  uint lambda = params.lambda;
  uint tau = params.tau;
  uint N = params.number_of_parties;
  uint s = params.first_challenge_size;
  uint r = params.target_rank;
  uint n = params.matrix_dimension.n;

  int ret = 1;
  // `salt` has 2 * lambda booleans, packed into lambda / 4 `uchar`s
  uint salt_size = ceil(lambda / 4.0);

  // pass `message` to the context
  if (!EVP_DigestUpdate(context, message, message_size))
    goto err;

  // pass `salt` to the context
  if (!EVP_DigestUpdate(context, salt.data, salt_size))
    goto err;

  // pass `h1` to the context
  if (!EVP_DigestUpdate(context, h1, salt_size))
    goto err;

  // pass `S, V` to the context
  // size of S: (s, r)
  // size of V: (s, n - r)
  uint total_size = tau * N * (s * r) * (s * (n - r));
  uchar *packed_matrices = malloc(sizeof(uchar) * total_size);
  pack_all_S_and_V(packed_matrices, parties, tau, N);

  // compute the hash
  if (!EVP_DigestFinal_ex(context, digest, NULL))
    goto err;
  // reinit the hash context
  if (!EVP_DigestInit_ex2(context, EVP_sha3_256(), NULL))
    goto err;

err:
  if (ret != 0)
    ERR_print_errors_fp(stderr);
}

/* Phase one: generate `params.tau` instances of the MPC protocol. The
 * result is stored in `commits`: a 2D array of SHA3 digests. */
void phase_one(uchar ***commits, uchar ***party_seeds, seed_t salt,
               PartyData **data, SignatureParameters params,
               MinRankInstance instance, MinRankSolution solution) {
  // unpack `params`
  uint N = params.number_of_parties;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint lambda = params.lambda;

  seed_t round_seed, party_seed;
  allocate_seed(&round_seed, lambda);
  allocate_seed(&party_seed, lambda);

  // initialize context for hashes (used in `commits`)
  EVP_MD_CTX *hash_context;
  initialize_sha3(&hash_context, lambda);

  // create various matrices
  Matrix alpha_sum, A_sum, C, C_sum, K_sum;
  MatrixSize A_size = {s, r};
  MatrixSize C_size = {s, n - r};

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

  for (uint round = 0; round < params.tau; round++) {
    fill_matrix_with_zero(&A_sum);
    fill_matrix_with_zero(&alpha_sum);
    fill_matrix_with_zero(&K_sum);
    fill_matrix_with_zero(&C_sum);

    generate_seed(round_seed);
    TreePRG(&salt, &round_seed, (size_t)lambda, N, party_seeds[round]);

    for (uint party = 0; party < N - 1; party++) {
      party_seed.data = party_seeds[round][party];
      PRG_init(&salt, &party_seed, lambda, &prg_state);
      generate_random_matrix(&data[round][party].A, prg_state, GF_16);
      generate_random_matrix(&data[round][party].alpha, prg_state, GF_16);
      generate_random_matrix(&data[round][party].C, prg_state, GF_16);
      generate_random_matrix(&data[round][party].K, prg_state, GF_16);

      // update `{matrix}_sum`
      matrix_sum(&alpha_sum, alpha_sum, data[round][party].alpha);
      matrix_sum(&A_sum, A_sum, data[round][party].A);
      matrix_sum(&K_sum, K_sum, data[round][party].K);
      matrix_sum(&C_sum, C_sum, data[round][party].C);

      hash0(commits[round][party], hash_context, salt, lambda, round, party,
            party_seed);
    }
    // last party is special
    party_seed.data = party_seeds[round][N - 1];
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

void prg_first_challenge(Matrix *challenges, uchar *h1,
                         SignatureParameters params) {
  uint lambda = params.lambda;
  uint tau = params.tau;

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  seed_t h1_wrap = {commit_size(lambda), h1};

  // create the challenges
  for (uint round = 0; round < tau; round++) {
    PRG_init(&h1_wrap, NULL, lambda, &prg_state);
    generate_random_matrix(&challenges[round], prg_state, GF_16);
  }

  gmp_randclear(prg_state);
}

/* Compute a hash of `message || salt || commits`, then use it to create
 * the challenges for each round. */
void phase_two(Matrix *challenges, uchar *h1, uchar *message, uint message_size,
               SignatureParameters params, uchar *salt, uchar ***commits) {
  uint lambda = params.lambda;
  uint tau = params.tau;
  uint N = params.number_of_parties;

  // compute the hash
  EVP_MD_CTX *context;
  initialize_sha3(&context, lambda);

  hash1(h1, context, message, message_size, salt, lambda, tau, N, commits);
  prg_first_challenge(challenges, h1, params);

  EVP_MD_CTX_free(context);
}

void phase_three(Matrix *challenges, MinRankInstance instance,
                 PartyState **parties, PartyData **data,
                 SignatureParameters params) {
  uint N = params.number_of_parties;
  uint s = params.first_challenge_size;
  uint r = params.target_rank; // rank of the solution
  uint matrix_count = instance.matrix_count;

  MatrixSize S_size = {s, r};
  Matrix S;
  allocate_matrix(&S, GF_16, S_size);

  for (uint round = 0; round < params.tau; round++) {
    for (uint party = 0; party < N; party++) {
      compute_local_m(&parties[round][party], instance,
                      data[round][party].alpha, matrix_count);
      compute_local_s(&parties[round][party], challenges[round],
                      data[round][party].A);
    }
    compute_global_s(&S, parties[round], N);
    for (uint party = 0; party < N; party++) {
      compute_local_v(&parties[round][party], S, data[round][party].K,
                      challenges[round], data[round][party].C);
    }
  }

  clear_matrix(&S);
}

void prg_second_challenge(uint *challenges, uchar *h2,
                          SignatureParameters params) {
  uint lambda = params.lambda;
  uint tau = params.tau;

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  seed_t h2_wrap = {commit_size(lambda), h2};

  PRG_init(&h2_wrap, NULL, lambda, &prg_state);

  uchar *random_bytes = malloc(sizeof(uchar) * 4 * tau);
  PRG_bytes(prg_state, 4 * tau, random_bytes);

  for (uint round = 0; round < tau; round++) {
    challenges[round] = ((uint)random_bytes[4 * round] << 24);
    challenges[round] += ((uint)random_bytes[4 * round + 1] << 16);
    challenges[round] += ((uint)random_bytes[4 * round + 2] << 8);
    challenges[round] += (uint)random_bytes[4 * round + 3];
  }

  gmp_randclear(prg_state);
}

void phase_four(uchar *h2, uint *second_challenges, uchar *message,
                uint message_size, seed_t salt, uchar *h1, PartyState **parties,
                SignatureParameters params) {
  uint lambda = params.lambda;

  // compute the hash
  EVP_MD_CTX *context;
  initialize_sha3(&context, lambda);
  hash2(h2, message, message_size, salt, h1, parties, context, params);

  // generate challenges
  prg_second_challenge(second_challenges, h2, params);

  EVP_MD_CTX_free(context);
}

void phase_five(uchar *digest, uint digest_size, uint *second_challenges,
                seed_t salt, uchar *h1, uchar *h2, uchar ***commits,
                uchar ***party_seeds, PartyState **parties, PartyData **data,
                SignatureParameters params) {
  uint tau = params.tau;
  uint N = params.number_of_parties;

  uint hash_size = commit_size(params.lambda);
  uint seed_size = (uint)ceil(params.lambda / 8.0);

  uint current_index = 0;
  uint bit_index;
  memcpy(digest, salt.data, salt.size);
  current_index += salt.size;

  memcpy(&digest[current_index], h1, hash_size);
  memcpy(&digest[current_index + hash_size], h1, hash_size);
  current_index += 2 * hash_size;

  for (uint round = 0; round < tau; round++) {
    // pack all states
    uint i_star = second_challenges[round] % N;
    for (uint party = 0; party < N; party++) {
      if (party != i_star) {
        memcpy(&digest[current_index], party_seeds[round][party], seed_size);
        current_index += seed_size;
        if (party == N) {
          bit_index = current_index << 3;
          pack_matrix(&digest[current_index], data[round][party].alpha,
                      &bit_index);
          pack_matrix(&digest[current_index], data[round][party].K, &bit_index);
          pack_matrix(&digest[current_index], data[round][party].C, &bit_index);
          current_index = (uint)ceil(bit_index / 8.0);
        }
      }
    }
    // pack commit of i_star
    memcpy(&digest[current_index], commits[round][i_star], hash_size);
    current_index += hash_size;

    // pack S matrix of i_star
    bit_index = current_index << 3;
    pack_matrix(&digest[current_index], parties[round][i_star].S, &bit_index);
    current_index = (uint)ceil(bit_index / 8.0);
  }

  digest[digest_size] = (uchar)0;
}

uint sign(uchar **digest, uchar *message, uint msg_size,
          PublicPrivateKeyPair key_pair, SignatureParameters params) {
  uint matrix_count = params.solution_size + 1;
  uint target_rank = params.target_rank;

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  // unpack the MinRank instance and solution from the key pair
  MinRankInstance instance;
  init_instance(&instance, matrix_count, params.matrix_dimension);
  unpack_instance_from_public_key(&instance, prg_state, params,
                                  key_pair.public_key);

  MinRankSolution solution;
  init_solution(&solution, matrix_count, target_rank, params.matrix_dimension);
  unpack_solution_from_private_key(&solution, prg_state, params,
                                   key_pair.private_key);

  // generate salt
  seed_t salt;
  allocate_seed(&salt, 2 * params.lambda);
  generate_seed(salt);

  // phase one
  uchar ***commits;
  allocate_all_commits(&commits, params);
  PartyData **data;
  allocate_all_party_data(&data, params);
  uchar ***party_seeds;
  allocate_all_party_seeds(&party_seeds, params);

  phase_one(commits, party_seeds, salt, data, params, instance, solution);

  // phase two
  Matrix *challenges;
  allocate_first_challenge(&challenges, params);
  PartyState **parties;
  allocate_all_parties(&parties, params);
  uchar *h1;
  allocate_commit(&h1, params.lambda);

  phase_two(challenges, h1, message, msg_size, params, salt.data, commits);

  // phase three
  phase_three(challenges, instance, parties, data, params);

  // phase four
  uchar *h2;
  allocate_commit(&h2, params.lambda);
  uint *second_challenges;
  allocate_second_challenge(&second_challenges, params);

  phase_four(h2, second_challenges, message, msg_size, salt, h1, parties,
             params);

  // phase five
  uint digest_size =
      allocate_signature_digest(digest, second_challenges, params);

  phase_five(*digest, digest_size, second_challenges, salt, h1, h2, commits,
             party_seeds, parties, data, params);

  clear_seed(&salt);
  clear_all_commits(commits, params);
  clear_all_party_data(data, params);
  clear_all_party_seeds(party_seeds, params);
  clear_first_challenge(challenges, params);
  clear_all_parties(parties, params);
  OPENSSL_free(h1);
  OPENSSL_free(h2);
  free(second_challenges);

  return digest_size;
}
