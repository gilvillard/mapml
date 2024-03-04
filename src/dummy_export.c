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


#ifndef MAPML_DUMMY_EXPORT_C
#define MAPML_DUMMY_EXPORT_C


#include "mapml.h"
#include "conversion.h"



/**********************************************************
 * 
 * maple polynomial round trip (for cheks) 
 * 
 * Converts a maple string polynomial representation 
 *  [deg+1  modulus coefficients]
 *  back to a maple polynomial   
 * 
 *  ALGEB args[1]: polynomial string 
 *        args[2]: modulus 
 * 
 ***********************************************************/


ALGEB polynomial_rt(MKernelVector kv, ALGEB *args){

    ALGEB stringpol=args[1];

    mp_limb_t modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_t p;  
    // The polynomial is initialized by nmod_poly_to_algeb

    get_nmod_poly(p, modulus, kv, stringpol);
    
    return nmod_poly_to_algeb(kv,p);

}

/**********************************************************
 * 
 * maple matrix polynomial round trip (for cheks) 
 * 
 * Converts a maple string polynomial matrix representation 
 *  [deg+1  modulus coefficients]
 *  back to a maple polynomial  matrix  
 * Internal: mat_poly
 *  
 *  ALGEB args[1]: polynomial matrix of string 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/


ALGEB matpoly_rt(MKernelVector kv, ALGEB *args){

    ALGEB stringmat=args[1];

    mp_limb_t modulus = MapleToInteger64(kv,args[2]);

    nmod_mat_poly_t A;

    get_nmod_mat_poly(A, modulus, kv, stringmat);
    
    return nmod_mat_poly_to_algeb(kv,A);

}


/**********************************************************
 * 
 * maple polynomial matrix round trip (for cheks) 
 * 
 * Converts a maple string polynomial matrix representation 
 *  [deg+1  modulus coefficients]
 *  back to a maple polynomial matrix 
 * Internal: poly_mat
 * 
 *  
 *  ALGEB args[1]: polynomial matrix of string 
 *        args[2]: modulus 
 * 
 * 
 ***********************************************************/


ALGEB polymat_rt(MKernelVector kv, ALGEB *args){

    ALGEB stringmat=args[1];

    mp_limb_t modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, modulus, kv, stringmat);
    
    return nmod_poly_mat_to_algeb(kv,A);

}

ALGEB polymat_rt2(MKernelVector kv, ALGEB *args){

    ALGEB stringmat=args[1];

    mp_limb_t modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat2(A, modulus, kv, stringmat);
    
    return nmod_poly_mat_to_algeb(kv,A);

}



#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
