#ifndef MATRIX_COMMON_H_
#define MATRIX_COMMON_H_

#include "../matrix.h"

typedef unsigned int uint;

uint **identity_matrix(MatrixSize size);
uint **rotation_matrix(MatrixSize size);

#endif // MATRIX_COMMON_H_
