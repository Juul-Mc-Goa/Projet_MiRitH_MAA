#include "../constants.h"
#include "../matrix.h"

#include "matrix_common.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("----------------------------------------------- beginning matrix product test...\n");
  MatrixSize size;
  size.m = 5;
  size.n = 5;
  FiniteField field = GF_16;

  uint **rotation = rotation_matrix(size);

  Matrix m_rot;
  allocate_matrix(&m_rot, field, size);
  copy_into_matrix(&m_rot, rotation);

  printf("circular permutation matrix:\n");
  print_matrix(&m_rot);

  Matrix m_prod;
  allocate_matrix(&m_prod, field, size);
  matrix_product(&m_prod, m_rot, m_rot);
  printf("matrix product (Rot * Rot):\n");
  print_matrix(&m_prod);

  clear_matrix(&m_rot);
  clear_matrix(&m_prod);
  free(rotation);
}
