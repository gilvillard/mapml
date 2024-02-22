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
 *  ALGEB args[1]: matrix polynomial string 
 *        args[2]: modulus 
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




#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
