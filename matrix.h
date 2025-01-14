#ifndef MATRIX_H_
#define MATRIX_H_

#include <gmp.h>

typedef unsigned int uint;

typedef struct MatrixSize {
  uint m;
  uint n;
} MatrixSize;

typedef struct Matrix {
  mpz_t **data;
  mpz_t moduli;
  MatrixSize size;
} Matrix;

void allocate_matrix(Matrix *m, MatrixSize size);
void matrix_set_ui(Matrix *m, uint **uint_matrix);
void fill_matrix_with_zero(Matrix *m);
void print_matrix(Matrix *m);
void clear_matrix(Matrix *m);

void matrix_sum(Matrix *result, Matrix left, Matrix right);
void matrix_big_sum(Matrix *result, Matrix *summands, uint k);
void matrix_opposite(Matrix *m);

void scalar_product(Matrix *result, mpz_t scalar, Matrix m);
void matrix_product(Matrix *result, Matrix m_left, Matrix m_right);

#endif // MATRIX_H_
