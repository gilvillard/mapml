/* 
   Copyright (C) 2024 Gilles Villard

   This file is part of mapml. mapml is free software: you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation,
   either version 2.1 of the License, or (at your option) any later version.

   mapml is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with mapml. If not, see <http://www.gnu.org/licenses/>. */

// *******************************************************


#ifndef MAPML_CONVERSION_C
#define MAPML_CONVERSION_C

#include "mapml.h"
#include "conversion.h"


/**********************************************************
 * 
 * Converts an nmod_poly_t to its maple equivalent 
 *   !! no modulus is transmitted
 * 
 * Todo: use a list instead, cf linbox  lb-maple.C. lbConvertPolynomial
 * 
 *********************************************************/


ALGEB nmod_poly_to_algeb(MKernelVector kv, const nmod_poly_t p){

    ALGEB maple_p; 

    slong d = nmod_poly_degree(p);

    maple_p = ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, 0));

    ALGEB f = EvalMapleStatement(kv,"proc(pol,c,d) return (pol +c*x^d); end proc;");

    for (slong i = 1; i <= d; i++) 
        maple_p = EvalMapleProc(kv, f, 3, maple_p, ToMapleInteger(kv,nmod_poly_get_coeff_ui(p, i)),ToMapleInteger(kv,i));

    return maple_p;   
}


/**********************************************************
 * 
 * Converts an nmod_mat_poly to its maple equivalent 
 *   !! no modulus is transmitted
 * 
 * Todo: convert matrix coefficients directly ? 
 *    not by global matrix conversion ? 
 * 
 ***********************************************************/

ALGEB nmod_mat_poly_to_algeb(MKernelVector kv, const nmod_mat_poly_t Ain){

    ALGEB maple_A; 

    slong m = Ain->r;
    slong n = Ain->c;

    nmod_poly_mat_t A;

    // Initialization necessary? 
    nmod_poly_mat_init(A, m, n , Ain->mod.n);

    nmod_poly_mat_set_from_mat_poly(A,Ain);

    //------------------------------
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

            tmp.dag =  nmod_poly_to_algeb(kv, nmod_poly_mat_entry(A,i-1,j-1));

            RTableAssign(kv, maple_A, index, tmp);
        }

    // Put this ?
    nmod_poly_mat_clear(A);

    return maple_A;

}


/**********************************************************
 * 
 * Converts an nmod_poly_mat to its maple equivalent 
 *   !! no modulus is transmitted
 * 
 ***********************************************************/

ALGEB nmod_poly_mat_to_algeb(MKernelVector kv, const nmod_poly_mat_t A){

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

            tmp.dag =  nmod_poly_to_algeb(kv, nmod_poly_mat_entry(A,i-1,j-1));

            RTableAssign(kv, maple_A, index, tmp);
        }

    return maple_A;

}


/**********************************************************
 * 
 * Converts a maple string polynomial representation 
 *  [deg+1  modulus coefficients]
 *  to an nmod_poly_t t  
 * 
 *  ALGEB stringpol: a string
 * 
 ***********************************************************/

void get_nmod_poly(nmod_poly_t p, const mp_limb_t modulus, MKernelVector kv, ALGEB stringpol){ 

    nmod_poly_init(p, modulus);

    nmod_poly_set_str(p, MapleToString(kv,stringpol));

}


// MapleToInteger64 vs slong and ulong 
// From a maple list args[1], [... , degree, coeff...], args[2] = length/2, nb of monomials 
// modulus assigned somewhere else 

// void get_nmod_poly_b(nmod_poly_t p, MKernelVector kv, ALGEB *args){

//     ALGEB coeffs = args[1]; // Work with coeffs in input? 

//     M_INT l = MapleToInteger64(kv,args[2]);

//     MapleALGEB_Printf(kv, " Length %a", ToMapleInteger(kv,l));

//     for (M_INT i=0; i<l; i++) {

//         nmod_poly_set_coeff_ui(p, MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+1)), 
//            MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+2)));

//     }
// }



/**********************************************************
 * 
 * Converts a maple string matrix polynomial representation 
 *  [i,j][deg+1  modulus coefficients]
 *  to an nmod_poly_mat_t   
 * 
 *  ALGEB string_A: a matrix of strings
 * 
 ***********************************************************/


void get_nmod_mat_poly(nmod_mat_poly_t Aout,   const mp_limb_t modulus, MKernelVector kv, ALGEB string_A){


    //ALGEB maple_A = args[1]; // Doesn't work with P in input? 

    M_INT m,n;   // slong flint or M_INT ? 

    m = RTableUpperBound(kv, string_A, 1);

    n = RTableUpperBound(kv, string_A, 2);

    nmod_poly_mat_t A;

    nmod_poly_mat_init(A, m, n , modulus); // Initializes the polynomials 

    M_INT index[2];

    RTableData tmp;             

    for (slong i=1; i<m+1; i++)
        for (slong j=1; j<n+1; j++){

            index[0]=i;
            index[1]=j;
            tmp = RTableSelect(kv,string_A,index);

            // polynomials have been initialized 
            nmod_poly_set_str(nmod_poly_mat_entry(A,i-1,j-1), MapleToString(kv,tmp.dag));

        }

    nmod_mat_poly_init(Aout, m, n, modulus);

    slong len = nmod_poly_mat_max_length(A);

    nmod_mat_poly_set_trunc_from_poly_mat(Aout,A,len);

    nmod_poly_mat_clear(A);

}


/**********************************************************
 * 
 * Converts a maple string polynomial matrix representation 
 *  [i,j][deg+1  modulus coefficients]
 *  to an nmod_poly_mat_t   
 * 
 *  ALGEB string_A: a matrix of strings
 * 
 ***********************************************************/


    void get_nmod_poly_mat(nmod_poly_mat_t A,   const mp_limb_t modulus, MKernelVector kv, ALGEB string_A){


    //ALGEB maple_A = args[1]; // Doesn't work with P in input? 

    M_INT m,n;   // slong flint or M_INT ? 

    m = RTableUpperBound(kv, string_A, 1);

    n = RTableUpperBound(kv, string_A, 2);

    nmod_poly_mat_init(A, m, n , modulus); // Initializes the polynomials 

    M_INT index[2];

    RTableData tmp;             

    for (slong i=1; i<m+1; i++)
        for (slong j=1; j<n+1; j++){

            index[0]=i;
            index[1]=j;
            tmp = RTableSelect(kv,string_A,index);

            // polynomials have been initialized 
            nmod_poly_set_str(nmod_poly_mat_entry(A,i-1,j-1), MapleToString(kv,tmp.dag));

        }


    }


#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
