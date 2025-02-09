/* Compare the result of `A(B + C)` and `AB + AC` */
#include "../constants.h"
#include "../matrix.h"
#include "../random.h"

#include <stdio.h>

int main(int argc, char **argv) {
  printf("------------------------------------- beginning the matrix "
         "distributivity test...\n");
  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate 3 random matrices
  Matrix A, B, C;
  MatrixSize size = {7, 7};
  allocate_matrix(&A, GF_16, size);
  allocate_matrix(&B, GF_16, size);
  allocate_matrix(&C, GF_16, size);

  generate_random_matrix(&A, random_state, GF_16);
  generate_random_matrix(&B, random_state, GF_16);
  generate_random_matrix(&C, random_state, GF_16);

  // compute intermediate matricse
  Matrix AB, AC, B_plus_C;
  allocate_matrix(&AB, GF_16, size);
  allocate_matrix(&AC, GF_16, size);
  allocate_matrix(&B_plus_C, GF_16, size);

  matrix_product(&AB, A, B);
  matrix_product(&AC, A, C);
  matrix_sum(&B_plus_C, B, C);

  // compute the two results
  Matrix result1, result2;
  allocate_matrix(&result1, GF_16, size);
  allocate_matrix(&result2, GF_16, size);

  matrix_sum(&result1, AB, AC);
  matrix_product(&result2, A, B_plus_C);

  printf("AB + AC:\n");
  print_matrix(&result1);
  printf("\nA(B + C):\n");
  print_matrix(&result2);

  clear_matrix(&A);
  clear_matrix(&B);
  clear_matrix(&C);
  clear_matrix(&AB);
  clear_matrix(&AC);
  clear_matrix(&B_plus_C);
  clear_matrix(&result1);
  clear_matrix(&result2);
}
