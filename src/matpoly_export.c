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


#ifndef MAPML_MATPOLY_EXPORT_C
#define MAPML_MATPOLY_EXPORT_C


#include "mapml.h"
#include "conversion.h"



/**********************************************************
 * 
 * modulo matrix polynomial determinant  
 * 
 *  ALGEB args[1]: matrix polynomial, vector entries 
 *        args[2]: modulus 
 * 
 * Returns the determinant as a list of coefficients 
 * 
 ***********************************************************/


ALGEB pm_determinant(MKernelVector kv, ALGEB *args){

   ALGEB stringmat=args[1];

   mp_limb_t modulus = MapleToInteger64(kv,args[2]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, stringmat);

   nmod_poly_t p;  

    nmod_poly_init(p, modulus); // Inittialization probably not done by flint below ? 

    nmod_poly_mat_det(p, A);
    
    return nmod_poly_to_algeb(kv,p);

 }


/**********************************************************
 * 
 * modulo matrix polynomial mbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns a polynomial matrix, list entries 
 * 
 *  Calls pml poly_mat i.e. matrix of polynomials 
 *     for divide and conquer w.r.t the order 
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

 ALGEB pm_matrix_mbasis(MKernelVector kv, ALGEB *args){

   ALGEB stringmat=args[2];

   mp_limb_t modulus = MapleToInteger64(kv,args[4]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, stringmat);

   slong m = A ->r;

   nmod_poly_mat_t M;

   nmod_poly_mat_init(M, m, m, modulus); 

   //slong res_shift[m];

   ALGEB maple_shift = args[1];

   slong shift[m];
   for (ulong i = 0; i < m; i++) 
      shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

   ulong order = MapleToInteger64(kv,args[3]);

   nmod_poly_mat_mbasis(M, shift, A, order);
   //mbasis(M, res_shift, A, order, shift);

   return nmod_poly_mat_to_algeb(kv,M);

}


//ALGEB coeffs = args[1]; // Work with coeffs in input? 

//     M_INT l = MapleToInteger64(kv,args[2]);

//     MapleALGEB_Printf(kv, " Length %a", ToMapleInteger(kv,l));

//     for (M_INT i=0; i<l; i++) {

//         nmod_poly_set_coeff_ui(p, ), 
//            MapleToInteger64(kv,MapleListSelect(kv,coeffs,2*i+2)));








#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
