
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

// Voir où les initialisations 

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

    slong m = A->r;
    slong n = A->c;

    RTableSettings setting;
    RTableGetDefaults(kv, &setting);
    setting.num_dimensions=2;
    setting.subtype=RTABLE_MATRIX;
    
    M_INT bounds[4];
    bounds[0] = 1;
    bounds[1] = m;
    bounds[2] = 1;
    bounds[3] = n;

    maple_A = RTableCreate(kv,&setting,NULL,bounds);


    M_INT index[2];
    RTableData tmp;

    for (slong i=1; i<m+1; i++)
        for (slong j=1; j<n+1; j++){

            index[0]=i;
            index[1]=j;

            tmp.dag =  to_maple_nmod_poly(kv, nmod_poly_mat_entry(A,i-1,j-1));

            RTableAssign(kv, maple_A, index, tmp);
        }

        return maple_A;

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

void get_nmod_poly(nmod_poly_t p, MKernelVector kv, ALGEB pol){ 

    nmod_poly_set_str(p, MapleToString(kv,pol));

}


ALGEB pget(MKernelVector kv, ALGEB *args){

    nmod_poly_t p;

    nmod_poly_init(p, 13);

    get_nmod_poly(p, kv, args[1]);
    
    return to_maple_nmod_poly(kv,p);

}



// !!!!  Modulus

void get_nmod_poly_mat(nmod_poly_mat_t A,   MKernelVector kv, ALGEB maple_A){


    //ALGEB maple_A = args[1]; // Doesn't work with P in input? 

    M_INT m,n;   // slong flint or M_INT ? 

    m = RTableUpperBound(kv, maple_A, 1);

    n = RTableUpperBound(kv, maple_A, 2);

    nmod_poly_mat_init(A, m, n , 13);

    M_INT index[2];

    RTableData tmp;             

    for (slong i=1; i<m+1; i++)
        for (slong j=1; j<n+1; j++){

            index[0]=i;
            index[1]=j;
            tmp = RTableSelect(kv,maple_A,index);

            nmod_poly_set_str(nmod_poly_mat_entry(A,i-1,j-1), MapleToString(kv,tmp.dag));

        }


}


ALGEB mget(MKernelVector kv, ALGEB *args){

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, kv, args[1]);
    
    return to_maple_nmod_poly_mat(kv,A);

}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
