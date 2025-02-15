#ifndef SEED_TREE_H_
#define SEED_TREE_H_

#include <gmp.h>

typedef unsigned int uint;
typedef unsigned char uchar;

void PRG_init(const uchar *salt, const uchar *seed, uint lambda,
              gmp_randstate_t *prg_state);
void PRG_bytes(gmp_randstate_t prg_state, size_t length, unsigned char *output);
void TreePRG(const uchar *salt, const uchar *seed, size_t lambda, size_t N,
             uchar **output_seeds);

#endif // SEED_TREE_H_
