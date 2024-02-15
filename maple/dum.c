
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include "maplec.h"

ALGEB retry(MKernelVector kv, ALGEB *args){

    mpz_ptr ptr = MapleToGMPInteger(kv, args[1]);

    return ToMapleInteger(kv,44);
    //return ToMapleNULL(kv);
}





/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
