
#include <time.h>
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include "maplec.h"

#include "nmod_poly_mat_extra.h"

// ******** Mettre des const **********

// ******** VOIR POUR LE MODULO *******


// Todo use a list for converting 
//    --> Cf linbox  lb-maple.C. lbConvertPolynomial

// Todo variable x 

ALGEB to_maple_nmod_poly(MKernelVector kv, nmod_poly_t p){

    ALGEB maple_p; 

    slong d = nmod_poly_degree(p);

    maple_p = ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, 0));

    ALGEB f = EvalMapleStatement(kv,"proc(pol,cc,dd) return (pol +cc*x^dd); end proc;");

    for (slong i = 1; i <= d; i++) 
        maple_p = EvalMapleProc(kv, f, 3, maple_p, ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, i)),ToMapleInteger(kv,i));

    return maple_p;   
}

ALGEB to_maple_nmod_poly_mat(MKernelVector kv, nmod_poly_mat_t A){

    ALGEB maple_A; 

    slong d = nmod_poly_degree(p);

    maple_p = ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, 0));

    ALGEB f = EvalMapleStatement(kv,"proc(pol,cc,dd) return (pol +cc*x^dd); end proc;");

    for (slong i = 1; i <= d; i++) 
        maple_p = EvalMapleProc(kv, f, 3, maple_p, ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, i)),ToMapleInteger(kv,i));

    return maple_p;   
}


// MapleToInteger64 vs slong and ulong 
// From a maple list args[1], [... , degree, coeff...], args[2] = length/2, nb of monomials 
// modulus assigned somewhere else 

void get_nmod_poly_b(nmod_poly_t p, MKernelVector kv, ALGEB *args){

    ALGEB coeffs = args[1]; // Work with coeffs in input? 

    M_INT l = MapleToInteger64(kv,args[2]);

    MapleALGEB_Printf(kv, " Length %a", ToMapleInteger(kv,l));

    for (M_INT i=0; i<l; i++) {

        nmod_poly_set_coeff_ui(p, MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+1)), 
                               MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+2)));

    }
}

// From a string at the flint format 

void get_nmod_poly(nmod_poly_t p, MKernelVector kv, ALGEB *args){ 

    nmod_poly_set_str(p, MapleToString(kv,args[1]));

}


ALGEB pget(MKernelVector kv, ALGEB *args){

    nmod_poly_t p;

    nmod_poly_init(p, 13);

    get_nmod_poly(p, kv, args);
    
    return to_maple_nmod_poly(kv,p);

}





ALGEB getmatrix(MKernelVector kv, ALGEB *args){


    ALGEB A = args[1]; // Doesn't work with P in input? 

    M_INT m,n;   // slong flint or M_INT ? 

    m = RTableUpperBound(kv, A, 1);

    n = RTableUpperBound(kv, A, 2);

    MapleALGEB_Printf(kv, " Rows %a columns %a",ToMapleInteger(kv,m),ToMapleInteger(kv,n));


    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    nmod_poly_mat_t mat;

    nmod_poly_mat_init(mat, m, n , 13);

    nmod_poly_mat_rand(mat, state, 4);

    slong i=1,j=2;

    nmod_poly_set_coeff_ui(mat->rows[i]+j, 0, 56);
    nmod_poly_set_coeff_ui(mat->rows[i]+j, 3, 17);


    nmod_poly_t res;
    nmod_poly_init(res,mat->modulus);
    nmod_poly_set(res,mat->rows[i]+j);

    // ALGEB maplepol; 

    // //--------------

    // slong d = nmod_poly_degree(res);


    // maplepol = ToMapleInteger(kv,nmod_poly_get_coeff_ui(res, 0));


    // ALGEB f = EvalMapleStatement(kv,"proc(pol,cc,dd) return (pol +cc*x^dd); end proc;");

    // for (slong i = 1; i <= d; i++) {

    // maplepol = EvalMapleProc(kv, f, 3, maplepol, ToMapleInteger(kv,nmod_poly_get_coeff_ui(res, i)),ToMapleInteger(kv,i));


    // }

    return to_maple_nmod_poly(kv,res); 


}






/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
