/* Compare the result of `A(BC)` and `(AB)C` */
#include "../constants.h"
#include "../matrix.h"
#include "../random.h"

#include <stdio.h>

int main(int argc, char **argv) {
  printf("------------------------------------- beginning the matrix "
         "associativity test...\n");
  gmp_randstate_t random_state;
  gmp_randinit_default(random_state);

  // generate 3 random matrices
  Matrix A, B, C;
  MatrixSize size = {2, 2};
  allocate_matrix(&A, GF_16, size);
  allocate_matrix(&B, GF_16, size);
  allocate_matrix(&C, GF_16, size);

  generate_random_matrix(&A, random_state, GF_16);
  generate_random_matrix(&B, random_state, GF_16);
  generate_random_matrix(&C, random_state, GF_16);

  printf("A:\n");
  print_matrix(&A);
  printf("\nB:\n");
  print_matrix(&B);
  printf("\nC:\n");
  print_matrix(&C);
  printf("\n");

  // compute intermediate matricse
  Matrix AB, BC;
  allocate_matrix(&AB, GF_16, size);
  allocate_matrix(&BC, GF_16, size);

  matrix_product(&AB, A, B);
  matrix_product(&BC, B, C);

  printf("AB:\n");
  print_matrix(&AB);
  printf("\nBC:\n");
  print_matrix(&BC);
  printf("\n");

  // compute the two results
  Matrix result1, result2;
  allocate_matrix(&result1, GF_16, size);
  allocate_matrix(&result2, GF_16, size);

  matrix_product(&result1, A, BC);
  matrix_product(&result2, AB, C);

  printf("A(BC):\n");
  print_matrix(&result1);
  printf("\n(AB)C:\n");
  print_matrix(&result2);
  printf("\n");

  clear_matrix(&A);
  clear_matrix(&B);
  clear_matrix(&C);
  clear_matrix(&AB);
  clear_matrix(&BC);
  clear_matrix(&result1);
  clear_matrix(&result2);
}
