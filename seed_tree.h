#ifndef SEED_TREE_H_
#define SEED_TREE_H_

#include "types.h"

#include <gmp.h>

void PRG_init(const seed_t *salt, const seed_t *seed, uint lambda,
              gmp_randstate_t *prg_state);
void PRG_bytes(gmp_randstate_t prg_state, size_t length, unsigned char *output);
void TreePRG(const seed_t *salt, const seed_t *seed, size_t lambda, size_t N,
             uchar **output_seeds);

#endif // SEED_TREE_H_
