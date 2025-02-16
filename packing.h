#ifndef PACKING_H_
#define PACKING_H_

#include "key_generation.h"
#include "random.h"

void allocate_packed_last_state(uchar *result, uint lambda,
                                MatrixSize alpha_size, MatrixSize K_size,
                                MatrixSize C_size);
void pack_boolean(uchar *result, bool element, uint *bit_index);
void pack_four_bit_value(uchar *result, uint element, uint *bit_index);
void pack_32bit_value(uchar *result, uint element, uint *bit_index);
void pack_last_state(uchar *result, seed_t salt, uint lambda, uint l, uint i,
                     Matrix alpha, Matrix K, Matrix C);
void pack_all_S_and_V(uchar *result, PartyState **parties, uint tau, uint N);
void unpack_solution_from_private_key(MinRankSolution *solution,
                                      gmp_randstate_t prg_state,
                                      SignatureParameters params,
                                      PrivateKey pk);
void unpack_instance_from_public_key(MinRankInstance *instance,
                                     gmp_randstate_t prg_state,
                                     SignatureParameters params,
                                     PublicKey pub_key);

#endif // PACKING_H_
