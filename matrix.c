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

  free(m.moduli);
  free(m.data);
}

void scalar_product(Matrix result, mpz_t scalar, Matrix m) {}
void matrix_sum(Matrix result, Matrix *summands) {}
void matrix_product(Matrix result, Matrix m_left, Matrix m_right) {}
