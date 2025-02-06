#include "matrix.h"
#include "constants.h"
#include "field_arithmetics.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned int uint;

void allocate_matrix(Matrix *result, FiniteField field, MatrixSize size) {
  result->size.m = size.m;
  result->size.n = size.n;
  result->data = malloc(size.m * sizeof(uint *));

  for (uint i = 0; i < size.m; i++) {
    result->data[i] = malloc(size.n * sizeof(uint));
  }

  result->field = field;
}

void clear_matrix(Matrix *m) {
  for (uint i = 0; i < m->size.m; i++) {
    free(m->data[i]);
  }
  free(m->data);
}

void copy_into_matrix(Matrix *m, uint **array) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = array[i][j];
    }
  }
}

void matrix_init_set(Matrix *m, Matrix src_matrix) {
  m->field = src_matrix.field;
  m->size = src_matrix.size;

  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = src_matrix.data[i][j];
    }
  }
}

void fill_matrix_with_zero(Matrix *m) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = 0;
    }
  }
}

/* Split `input` into a left part and a right part. */
void split_matrix(Matrix *left, Matrix *right, Matrix input, uint target_rank) {
  left->size.m = input.size.m;
  left->size.n = input.size.n - target_rank;
  left->data = input.data;

  right->size.m = input.size.m;
  right->size.n = target_rank;
  // store the pointer to element of index `mid`
  for (uint i = 0; i < input.size.m; i++) {
    right->data[i] = &input.data[i][left->size.n];
  }
}

/* Split each matrix in `input`. */
void split_each_matrix(Matrix *left, Matrix *right, Matrix *input,
                       uint target_rank, uint length) {
  for (uint i = 0; i < length; i++) {
    split_matrix(&left[i], &right[i], input[i], target_rank);
  }
}

bool matrix_is_zero(Matrix m) {
  for (uint i = 0; i < m.size.m; i++) {
    for (uint j = 0; j < m.size.n; j++) {
      if (m.data[i][j] != 0) {
        return false;
      }
    }
  }
  return true;
}

void print_matrix(Matrix *m) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      printf("%u ", m->data[i][j]);
    }
    printf("\n");
  }
}

void vector_sum(uint *result, uint length, uint *left, uint *right) {
  for (uint i = 0; i < length; i++) {
    result[i] = scalar_add(left[i], right[i]);
  }
}

/* Sum 2 matrices, store the output in `result`. */
void matrix_sum(Matrix *result, Matrix left, Matrix right) {
  uint m = left.size.m, n = left.size.n;

  for (uint i = 0; i < m; i++) {
    for (uint j = 0; j < n; j++) {
      result->data[i][j] = scalar_add(left.data[i][j], right.data[i][j]);
    }
  }
}

/* Sum a list of `k` matrices, store the output in `result`. */
void matrix_big_sum(Matrix *result, Matrix *summands, uint k) {
  uint m = summands[0].size.m, n = summands[0].size.n;

  for (uint i = 0; i < m; i++) {
    for (uint j = 0; j < n; j++) {
      result->data[i][j] = 0;
      for (uint l = 0; l < k; l++) {
        result->data[i][j] =
            scalar_add(result->data[i][j], summands[l].data[i][j]);
      }
    }
  }
}

void matrix_big_weighted_sum(Matrix *result, uint *weights, Matrix *summands,
                             uint k) {
  fill_matrix_with_zero(result);

  Matrix temp;
  allocate_matrix(&temp, GF_16, result->size);

  for (uint i = 0; i < k; i++) {
    scalar_product(&temp, weights[i], summands[i]);
    matrix_sum(result, *result, temp);
  }
}

void vector_opposite(uint *vec, uint length) {
  for (uint i = 0; i < length; i++) {
    vec[i] = scalar_neg(vec[i]);
  }
}

void matrix_opposite(Matrix *m) {
  for (uint i = 0; i < m->size.m; i++) {
    for (uint j = 0; j < m->size.n; j++) {
      m->data[i][j] = scalar_neg(m->data[i][j]);
    }
  }
}

/* Multiply a matrix `m` by a gmp integer `scalar`, store the output in
 * `result`. */
void scalar_product(Matrix *result, uint scalar, Matrix m) {
  for (uint i = 0; i < m.size.m; i++) {
    for (uint j = 0; j < m.size.n; j++) {
      result->data[i][j] = scalar_mul(scalar, m.data[i][j]);
    }
  }
}

/* Multiply two matrices. No check is done on the matrices dimensions.
 * So it is assumed that
 * - `m_left.size.n == m_right.size.m`,
 * - `result.size.m == m_left.size.m`,
 * - `result.size.n == m_right.size.n`. */
void matrix_product(Matrix *result, Matrix m_left, Matrix m_right) {
  uint mid_dimension = m_left.size.n;

  for (uint i = 0; i < result->size.m; i++) {
    for (uint j = 0; j < result->size.n; j++) {
      // compute `result[i][j] = 0`
      result->data[i][j] = 0;
      for (uint l = 0; l < mid_dimension; l++) {
        // compute `result[i][j] += m_left[i][l] * m_right[l][j]`
        result->data[i][j] ^= scalar_mul(m_left.data[i][l], m_right.data[l][j]);
      }
    }
  }
}
