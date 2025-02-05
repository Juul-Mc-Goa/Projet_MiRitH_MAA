#ifndef MATRIX_H_
#define MATRIX_H_

#include "field_arithmetics.h"
#include <gmp.h>
#include <stdbool.h>

typedef unsigned int uint;

typedef struct MatrixSize {
  uint m;
  uint n;
} MatrixSize;

typedef struct Matrix {
  uint **data;
  FiniteField field;
  MatrixSize size;
} Matrix;

void allocate_matrix(Matrix *result, FiniteField field, MatrixSize size);
void copy_into_matrix(Matrix *m, uint **array);
void matrix_init_set(Matrix *m, Matrix src_matrix);
void fill_matrix_with_zero(Matrix *m);
bool matrix_is_zero(Matrix m);
void print_matrix(Matrix *m);
void clear_matrix(Matrix *m);

void vector_sum(uint *result, uint length, uint *left, uint *right);
void vector_opposite(uint *vec, uint length);

void matrix_sum(Matrix *result, Matrix left, Matrix right);
void matrix_big_sum(Matrix *result, Matrix *summands, uint k);
void matrix_big_weighted_sum(Matrix *result, uint *weights, Matrix *summands,
                             uint k);
void matrix_opposite(Matrix *m);

void scalar_product(Matrix *result, uint scalar, Matrix m);
void matrix_product(Matrix *result, Matrix m_left, Matrix m_right);

#endif // MATRIX_H_
