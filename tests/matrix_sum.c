#include "../constants.h"
#include "../matrix.h"

#include "matrix_common.h"

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  printf("------------------------------ beginning matrix sum test...\n");
  MatrixSize size = {5, 5};
  FiniteField field = GF_16;

  uint **identity = identity_matrix(size);
  uint **rotation = rotation_matrix(size);

  Matrix m_id;
  allocate_matrix(&m_id, field, size);
  copy_into_matrix(&m_id, identity);

  printf("Id matrix:\n");
  print_matrix(&m_id);

  Matrix m_rot;
  allocate_matrix(&m_rot, field, size);
  copy_into_matrix(&m_rot, rotation);

  printf("circular permutation matrix:\n");
  print_matrix(&m_rot);

  Matrix m_sum;
  allocate_matrix(&m_sum, field, size);
  matrix_sum(&m_sum, m_id, m_rot);
  printf("matrix sum (Id + Rot):\n");
  print_matrix(&m_sum);

  clear_matrix(&m_id);
  clear_matrix(&m_rot);
  clear_matrix(&m_sum);
  free(identity);
  free(rotation);
}
