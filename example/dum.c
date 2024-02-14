#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <stdlib.h>


#include "nmod_poly_mat_extra.h"
#include "sagemath_extra.h"

#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>
#include <flint/profiler.h>

// verify determinant
int verify_determinant(const nmod_poly_t det, const nmod_poly_mat_t mat, flint_rand_t state)
{
    // checking dimensions
    if (mat->r != mat->c)
    {
        printf("~~~ verify determinant ~~~ INCORRECT: dimension mismatch\n");
        return 0;
    }

    // compute determinant via Flint's native method
    nmod_poly_t det_correct;
    nmod_poly_init(det_correct, mat->modulus);
    nmod_poly_mat_det(det_correct, mat);

    if (! nmod_poly_equal(det, det_correct))
    {
        printf("~~~ verify determinant ~~~ INCORRECT: determinant is wrong\n");
        nmod_poly_print_pretty(det, "X"); printf("\n");
        nmod_poly_print_pretty(det_correct, "X"); printf("\n");
        nmod_poly_clear(det_correct);
        return 0;
    }
    else
    {
        nmod_poly_clear(det_correct);
        return 1;
    }
}



int main(int argc, char ** argv)
{


    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());


    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, 4, 4 , 2);
    nmod_poly_mat_rand(mat, state, 2);


    nmod_poly_t det;
    nmod_poly_init(det, mat->modulus);

    nmod_poly_mat_det(det, mat);

    nmod_poly_print_pretty(det, "x"); printf("\n");


    nmod_poly_mat_print_pretty(mat,"x"); printf("\n");

    nmod_poly_mat_clear(mat);




}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
