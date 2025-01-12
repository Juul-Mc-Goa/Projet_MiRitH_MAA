#include "matrix.h"
#include <gmp.h>
#include <stdlib.h>

typedef unsigned int uint;

void allocate_matrix(Matrix m, MatrixSize size) {
  Matrix result;
  result.size = size;
  result.data = malloc(size.m * sizeof(mpz_t *));

  for (uint i = 0; i < size.m; i++) {
    result.data[i] = malloc(size.n * sizeof(mpz_t));
    for (uint j = 0; i < size.n; j++) {
      mpz_init(result.data[i][j]);
    }
  }

  mpz_init(result.moduli);
}

void clear_matrix(Matrix m) {
  for (uint i = 0; i < m.size.m; i++) {
    for (uint j = 0; i < m.size.n; j++) {
      mpz_clear(m.data[i][j]);
    }
    free(m.data[i]);
  }

  mpz_clear(m.moduli);
  free(m.data);
}

/* Multiply a matrix `m` by a gmp integer `scalar`, store the output in
 * `result`. */
void scalar_product(Matrix result, mpz_t scalar, Matrix m) {
  for (uint i = 0; i < m.size.m; i++) {
    for (uint j = 0; i < m.size.n; j++) {
      mpz_mul(m.data[i][j], m.data[i][j], scalar);
    }
  }
}

/* Sum a list of `k` matrices, store the output in `result`. */
void matrix_sum(Matrix result, Matrix *summands, uint k) {
  uint m = summands[0].size.m, n = summands[0].size.n;

  for (uint i = 0; i < m; i++) {
    for (uint j = 0; i < n; j++) {
      mpz_set_ui(result.data[i][j], 0);
      for (uint l = 0; l < k; l++) {
        mpz_add(result.data[i][j], result.data[i][j], summands[l].data[i][j]);
      }
    }
  }
}

/* Multiply two matrices. No check is done on the matrices dimensions.
 * So it is assumed that
 * - `m_left.size.n == m_right.size.m`,
 * - `result.size.m == m_left.size.m`,
 * - `result.size.n == m_right.size.n`. */
void matrix_product(Matrix result, Matrix m_left, Matrix m_right) {
  uint m = result.size.m, n = result.size.n, mid_dimension = m_left.size.n;

  for (uint i = 0; i < m; i++) {
    for (uint j = 0; i < n; j++) {
      // compute `result[i][j] = 0`
      mpz_set_ui(result.data[i][j], 0);
      for (uint l = 0; l < mid_dimension; l++) {
        // compute `result[i][j] += m_left[i][l] * m_right[l][j]`
        mpz_addmul(result.data[i][j], m_left.data[i][l], m_right.data[l][j]);
      }
    }
  }
}
