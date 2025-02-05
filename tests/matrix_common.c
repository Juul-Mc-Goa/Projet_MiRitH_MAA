#include "../field_arithmetics.h"
#include "../matrix.h"

#include <stdlib.h>

uint **identity_matrix(MatrixSize size) {
  uint **result = malloc(size.m * sizeof(uint *));
  for (uint i = 0; i < size.m; i++) {
    result[i] = malloc(size.n * sizeof(uint));
    for (uint j = 0; j < size.n; j++) {
      result[i][j] = (i == j) ? 1 : 0;
    }
  }

  return result;
}

uint **rotation_matrix(MatrixSize size) {
  uint **result = malloc(size.m * sizeof(uint *));
  for (uint i = 0; i < size.m; i++) {
    result[i] = malloc(size.n * sizeof(uint));
    for (uint j = 0; j < size.n; j++) {
      if ((j + 1) % size.n == i) {
        result[i][j] = 1;

      } else {
        result[i][j] = 0;
      }
    }
  }

  return result;
}
