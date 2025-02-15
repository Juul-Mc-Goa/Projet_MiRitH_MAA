#include "constants.h"
#include "key_generation.h"
#include "matrix.h"
#include "mpc.h"
#include "random.h"
#include "seed_tree.h"
#include "types.h"

#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

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
