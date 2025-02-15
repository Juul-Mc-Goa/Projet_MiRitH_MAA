#ifndef TYPES_H_
#define TYPES_H_

typedef unsigned char uchar;
typedef unsigned int uint;

/* A seed is stored as an array of `unsigned char`, along with its length. */
typedef struct seed_t {
  uint size;
  uchar *data;
} seed_t;

/* To store the field F_{2^k}, we store `2^k`, `k`, and the polynomial used for
 * the field definition. */
typedef struct FiniteField {
  uint field_size;
  uint log_field_size;
  uint polynomial;
} FiniteField;

/* Size of a matrix: `m` is the number of lines, `n` is the number of rows. */
typedef struct MatrixSize {
  uint m;
  uint n;
} MatrixSize;

/* A matrix, stored as a 2D array of `unsigned int`, plus a field for the finite
 * field, and another one for the matrix size. */
typedef struct Matrix {
  uint **data;
  FiniteField field;
  MatrixSize size;
} Matrix;

// A struct holding all the parameters for the signature
typedef struct SignatureParameters {
  uint lambda;                 // security parameter
  MatrixSize matrix_dimension; // dimension (m, n) of the full matrices
  FiniteField field;           // the field used for computations
  uint target_rank;            // the rank `r` of the solution
  uint solution_size;          // the size `k` of the solution vector
  uint first_challenge_size;   // the number of rows `s` in the first challenge
  uint number_of_parties;      // the number of parties `N`
  uint tau;                    // the number of rounds
} SignatureParameters;

// A struct holding
// 1. a `seed`, in a array of `lambda` booleans
// 2. a matrix `m0` containing elements of F_q (stored as gnu MP integers)
typedef struct PublicKey {
  uint lambda;
  seed_t seed;
  Matrix m0;
} PublicKey;

// A struct holding a `seed` in an array of `lambda` booleans.
typedef struct PrivateKey {
  uint lambda;
  seed_t seed;
} PrivateKey;

// A container for a public/private key pair. Returned by the `key_gen`
// function.
typedef struct PublicPrivateKeyPair {
  PublicKey public_key;
  PrivateKey private_key;
} PublicPrivateKeyPair;

typedef struct MinRankInstance {
  uint matrix_count;
  Matrix *matrix_array;
} MinRankInstance;

typedef struct MinRankSolution {
  uint matrix_count;
  uint target_rank;
  Matrix alpha;
  Matrix K;
} MinRankSolution;

typedef struct PartyState {
  Matrix M_left;
  Matrix M_right;
  Matrix S;
  Matrix V;
} PartyState;

typedef struct PartyData {
  Matrix alpha;
  Matrix A;
  Matrix R;
  Matrix C;
  Matrix K;
} PartyData;

#endif // TYPES_H_
