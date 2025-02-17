#ifndef VERIF_H_
#define VERIF_H_

#include "types.h"

void verif(PublicKey key, uchar *message, uchar *digest, uint digest_size,
           SignatureParameters params) ;

#endif // VERIF_H_
