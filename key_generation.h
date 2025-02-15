#ifndef KEY_GENERATION_H_
#define KEY_GENERATION_H_

#include "field_arithmetics.h"
#include "matrix.h"
#include "types.h"
#include <gmp.h>
#include <stdbool.h>

typedef unsigned int uint;


void allocate_key_pair(PublicPrivateKeyPair *key_pair,
                       SignatureParameters parameters);
void clear_key_pair(PublicPrivateKeyPair key_pair);
void key_gen(PublicPrivateKeyPair *result,
                             SignatureParameters params);

#endif // KEY_GENERATION_H_
