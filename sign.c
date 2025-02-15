#include "constants.h"
#include "matrix.h"
#include "mpc.h"
#include "packing.h"
#include "random.h"

#include <openssl/bio.h>
#include <openssl/crypto.h>
#include <openssl/err.h>
#include <openssl/evp.h>


void allocate_commit(uchar *commit, uint lambda) {
  commit = malloc(sizeof(uchar) * lambda);
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

/* Phase one: generate `params.tau` instances of the MPC protocol. The
 * result is stored in `commits`: a 2D array of SHA3 digests. */
int phase_one(uchar ***commits, SignatureParameters params,
              MinRankInstance instance, MinRankSolution solution) {
  // unpack `params`
  uint number_of_parties = params.number_of_parties, lambda = params.lambda,
       tau = params.tau, matrix_count = params.solution_size + 1;
  uint ret = 1;

  EVP_MD_CTX *random_state;
  initialize_shake256(&random_state);

  // generate salt
  seed_t salt;
  allocate_seed(&salt, 2 * lambda);
  generate_seed(salt);

  seed_t current_round_seed, current_party_seed;
  allocate_seed(&current_round_seed, lambda);
  allocate_seed(&current_party_seed, lambda);

  // initialize context for hashes (used in `commits`)
  EVP_MD_CTX *hash_context = NULL;
  ret = initialize_sha3(&hash_context, lambda);
  if (ret != 0)
    goto err;

  // create various matrices
  Matrix alpha, alpha_sum, A, A_sum, C, C_sum, K, K_sum;
  MatrixSize alpha_size = {1, matrix_count},
             A_size = {params.first_challenge_size, params.target_rank},
             K_size = {params.target_rank,
                       params.matrix_dimension.n - params.target_rank},
             C_size = {A_size.m, K_size.n};

  allocate_matrix(&alpha, GF_16, alpha_size);
  allocate_matrix(&alpha_sum, GF_16, alpha_size);

  allocate_matrix(&A, GF_16, A_size);
  allocate_matrix(&A_sum, GF_16, A_size);

  allocate_matrix(&K, GF_16, K_size);
  allocate_matrix(&K_sum, GF_16, K_size);

  allocate_matrix(&C, GF_16, C_size);
  allocate_matrix(&C_sum, GF_16, C_size);

  // for each round
  for (uint l = 0; l < tau; l++) {
    fill_matrix_with_zero(&A_sum);
    fill_matrix_with_zero(&alpha_sum);
    fill_matrix_with_zero(&K_sum);
    fill_matrix_with_zero(&C_sum);

    // generate a seed
    generate_seed(current_round_seed);

    // use this seed to generate `number_of_parties` seeds
    for (uint i = 0; i < number_of_parties - 1; i++) {
      generate_seed_from_randstate(current_party_seed, salt, current_round_seed,
                                   lambda, &random_state);
      // use this fresh seed to generate A, alpha, K, C
      seed_random_state(current_party_seed, &hash_context);
      generate_random_matrix(&A, random_state, GF_16);
      generate_random_matrix(&alpha, random_state, GF_16);
      generate_random_matrix(&K, random_state, GF_16);
      generate_random_matrix(&C, random_state, GF_16);

      // update `{matrix}_sum`
      matrix_sum(&alpha_sum, alpha_sum, alpha);
      matrix_sum(&A_sum, A_sum, A);
      matrix_sum(&K_sum, K_sum, K);
      matrix_sum(&C_sum, C_sum, C);

      // store the commit
      allocate_hash_digest(&commits[l][i], hash_context, lambda);

      if (commits[l][i] == NULL)
        goto err;

      ret = hash0(commits[i][l], hash_context, salt, lambda, l, i,
                  current_party_seed);

      if (ret != 0)
        goto err;
    }
    // last party is special: we need to enforce various relations
    generate_seed_from_randstate(current_party_seed, salt, current_round_seed,
                                 lambda, &random_state);
    seed_random_state(current_party_seed, &random_state);
    // generate last A
    generate_random_matrix(&A, random_state, GF_16);
    // compute last alpha, K, C
    matrix_sum(&alpha, solution.alpha, alpha_sum);
    matrix_sum(&K, solution.K, K_sum);
    matrix_product(&C, A, K);
    matrix_sum(&C, C, C_sum);
    allocate_hash_digest(&commits[l][number_of_parties - 1], hash_context,
                         lambda);

    if (commits[l][number_of_parties - 1] == NULL)
      goto err;

    ret = hash0_last(commits[l][number_of_parties - 1], hash_context, salt,
                     lambda, l, number_of_parties - 1, current_party_seed,
                     alpha, K, C);

    if (ret != 0)
      goto err;
  }

  // TODO: execute MPC protocol, store the shares of S and V
  // TODO: compute h2: hash of `msg || salt || (S, V)_i`
  //       use it to generate second challenges
  // TODO: output `h1 || h2 || (almost each (state_i, com_i, S_i))`

err:
  clear_seed(&current_round_seed);
  clear_seed(&current_party_seed);
  clear_matrix(&alpha);
  clear_matrix(&alpha_sum);
  clear_matrix(&A);
  clear_matrix(&A_sum);
  clear_matrix(&C);
  clear_matrix(&C_sum);
  clear_matrix(&K);
  clear_matrix(&K_sum);

  EVP_MD_CTX_free(hash_context);
  EVP_MD_CTX_free(random_state);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}

/* Compute a hash of `message || salt || commits`, then use it to create
 * the challenges for each round. */
int phase_two(Matrix *challenges, uchar *message, uint message_size,
              uint number_of_rounds, uint number_of_parties, uint lambda,
              uchar *salt, uchar ***commits) {
  uint ret = 1;

  // compute the hash
  EVP_MD_CTX *context = NULL;
  ret = initialize_sha3(&context, lambda);
  if (ret != 0)
    goto err;

  uchar *digest;
  allocate_hash_digest(&digest, context, lambda);
  hash1(digest, context, message, message_size, salt, lambda, number_of_rounds,
        number_of_parties, commits);

  // create the challenges

err:
  EVP_MD_CTX_free(context);
  if (ret != 0)
    ERR_print_errors_fp(stderr);
  return ret;
}
void allocate_party_seeds(uchar **party_seeds, uint number_of_parties,
                          uint lambda) {
  size_t seed_size = (size_t)ceil(lambda / 8.0);
  party_seeds = malloc(number_of_parties * sizeof(uchar *));

  for (uint party = 0; party < number_of_parties; party++) {
    party_seeds[party] = malloc(seed_size * sizeof(uchar));
  }
}

void sign_phase1(Matrix *M_share, Matrix alpha, Matrix K, uchar *message,
                 SignatureParameters params) {
  uint k = params.solution_size;
  uint N = params.number_of_parties;
  uint m = params.matrix_dimension.m;
  uint n = params.matrix_dimension.n;
  uint r = params.target_rank;
  uint s = params.first_challenge_size;
  uint lambda = params.lambda;

  uchar **party_seeds;
  allocate_party_seeds(party_seeds, N, lambda);

  seed_t salt, round_seed, party_seed;
  allocate_seed(&salt, 2 * lambda);
  allocate_seed(&round_seed, lambda);
  allocate_seed(&party_seed, lambda);

  generate_seed(salt, 2 * lambda);

  gmp_randstate_t prg_state;
  gmp_randinit_default(prg_state);

  Matrix A_share, alpha_share, C_share, K_share;
  MatrixSize A_size = {s, r};
  MatrixSize alpha_size = {1, k + 1};
  MatrixSize C_size = {s, r};
  MatrixSize K_size = {r, n - r};

  allocate_matrix(&A_share, GF_16, A_size);
  allocate_matrix(&alpha_share, GF_16, alpha_size);
  allocate_matrix(&C_share, GF_16, C_size);
  allocate_matrix(&K_share, GF_16, K_size);

  uchar ***commits = malloc(sizeof(uchar **) * params.tau);

  for (uint round = 0; round < params.tau; round++) {
    generate_seed(round_seed, lambda);
    TreePRG(&salt, &round_seed, (size_t)lambda, N, party_seeds);

    commits[round] = malloc(sizeof(uchar *) * N);

    for (uint party = 0; party < N - 1; party++) {
      party_seed.data = party_seeds[party];
      PRG_init(&salt, &party_seed, lambda, &prg_state);
      generate_random_matrix(&A_share, prg_state, GF_16);
      generate_random_matrix(&alpha_share, prg_state, GF_16);
      generate_random_matrix(&C_share, prg_state, GF_16);
      generate_random_matrix(&K_share, prg_state, GF_16);
    }
  }
}

/* Phase 1 of signing
 *
 * Given as input:
 * - the salt 'salt';
 * - the round 'l';
 * - the seeds 'seed' of the round;
 * - the secret key 'a', 'K';
 * compute the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr',
 * the matrix 'A', and the commitments 'com'. */
/* void sign_phase1(Matrix A, Matrix A_shr, Matrix C_shr, Matrix a_shr, */
/*                  Matrix K_shr, hash_t com[N_PARTIES], */

/*                  hash_t salt, uint32_t l, seed_t seed[N_PARTIES], */
/*                  const uint8_t a[matrix_bytes_size(K, 1)], */
/*                  const uint8_t K[matrix_bytes_size(R, N - R)]) { */
/*   uint32_t i; */

/*   for (i = 0; i < N_PARTIES; i++) { */
/*     prng_t prng; */

/*     // Initialize PRNG from 'seed[i]' */
/*     prng_init(&prng, salt, seed[i]); */

/*     // Generate the random matrix 'A_shr[i]' */
/*     matrix_init_random(A_shr[i], S, R, &prng); */

/*     if (i != N_PARTIES - 1) { */
/*       // Generate random matrices 'a_shr[i]', 'C_shr[i]', 'K_shr[i]' */
/*       matrix_init_random(a_shr[i], K, 1, &prng); */
/*       matrix_init_random(C_shr[i], S, N - R, &prng); */
/*       matrix_init_random(K_shr[i], R, N - R, &prng); */

/*       // Compute the commitment */
/*       hash_digest0(com[i], salt, l, i, seed[i]); */
/*     } else { */
/*       uint32_t j; */

/*       // Set a_shr[N_PARTIES - 1] = a - sum_{j < N_PARTIES - 1} a_shr[j] */
/*       matrix_copy(a_shr[N_PARTIES - 1], a, K, 1); */

/*       for (j = 0; j < N_PARTIES - 1; j++) { */
/*         matrix_subtract(a_shr[N_PARTIES - 1], a_shr[j], K, 1); */
/*       } */

/*       // Set K_shr[N_PARTIES - 1] = K - sum_{j < N_PARTIES - 1} K_shr[j] */
/*       matrix_copy(K_shr[N_PARTIES - 1], K, R, N - R); */

/*       for (j = 0; j < N_PARTIES - 1; j++) { */
/*         matrix_subtract(K_shr[N_PARTIES - 1], K_shr[j], R, N - R); */
/*       } */

/*       // Open 'A' by summing all A_shr */
/*       matrix_init_zero(A, S, R); */

/*       for (j = 0; j < N_PARTIES; j++) { */
/*         matrix_add(A, A_shr[j], S, R); */
/*       } */

/*       // Compute C_shr[N_PARTIES - 1] = A * K - sum_{j < N_PARTIES - 1}
 * C_shr[j] */
/*       matrix_product(C_shr[N_PARTIES - 1], A, K, S, R, N - R); */

/*       for (j = 0; j < N_PARTIES - 1; j++) { */
/*         matrix_subtract(C_shr[N_PARTIES - 1], C_shr[j], S, N - R); */
/*       } */

/*       // Compute the commitment with auxiliary data */
/*       hash_digest0_aux(com[i], salt, l, i, seed[N_PARTIES - 1], */
/*                        a_shr[N_PARTIES - 1], K_shr[N_PARTIES - 1], */
/*                        C_shr[N_PARTIES - 1]); */
/*     } */
/*   } */
/* } */

/* Phase 3 of signing
 *
 * Given as input:
 * - the public key 'M';
 * - the secret key 'E';
 * - the matrix 'A';
 * - the additive shares 'A_shr', 'C_shr', 'a_shr', 'K_shr';
 * computes and hashes the additive shares 'S_shr', 'V_shr'.
 * The matrix 'S_shr[i]' is written over 'A_shr'. */
void sign_phase3(uint8_t M[K + 1][matrix_bytes_size(M, N)],
                 uint8_t E[matrix_bytes_size(M, N)],
                 uint8_t A[matrix_bytes_size(S, R)],
                 uint8_t A_shr[N_PARTIES][matrix_bytes_size(S, R)],
                 uint8_t C_shr[N_PARTIES][matrix_bytes_size(S, N - R)],
                 uint8_t a_shr[N_PARTIES][matrix_bytes_size(K, 1)],
                 uint8_t K_shr[N_PARTIES][matrix_bytes_size(R, N - R)],
                 uint8_t R[matrix_bytes_size(S, M)]) {
  // MPC (Multi-Party Computation) implementation goes here
}

/* Phase 4 of signing
 *
 * Given as input:
 * - the hash hash2;
 * computes the challenge i_star */
void sign_phase4(uint32_t i_star[TAU], const hash_t hash2) {
  size_t l;
  prng_t prng;

  // Initialize PRNG from 'hash2'
  prng_init(&prng, hash2, NULL);

  for (l = 0; l < TAU; l++) {
    uint32_t r;
    // Generate random number and compute i_star[l]
    prng_bytes(&prng, &r, sizeof(r));
    i_star[l] = r % N_PARTIES;
  }
}
