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

#include <string.h>
#include "mapml.h"
#include "conversion.h"


ALGEB coeffs(MKernelVector kv, ALGEB *args){

    ALGEB pol=args[1];

    fmpq_poly_t p;

    //get_fmpq_poly(p, kv, pol); 

    //MapleALGEB_Printf(kv, " ++++  \n");

    fmpq_t res;

    fmpq_init(res);

    //fmpq_set_str(q, MapleToString(kv,rat), 10);

    MapleALGEB_Printf(kv, " \n");
    MapleALGEB_Printf(kv, " ++++  \n");


    fmpq_poly_get_coeff_fmpq(res, p, 0);
    
    MapleALGEB_Printf(kv, " \n");
    MapleALGEB_Printf(kv, " ++++  \n");

    ALGEB s;
    char *str;
    str = fmpq_get_str(str, 10, res);

    //*s=ToMapleString(kv,str);

    MapleALGEB_Printf(kv, " ---- \n");
    MapleALGEB_Printf(kv, str);

    return ToMapleString(kv,str);

}


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

    ALGEB vectmat=args[1];

    mp_limb_t modulus = MapleToInteger64(kv,args[2]);

    nmod_poly_mat_t A;

    get_nmod_poly_mat(A, modulus, kv, vectmat);

    nmod_poly_t p;  

    nmod_poly_init(p, modulus); // Inittialization probably not done by flint below ? 

    nmod_poly_mat_det(p, A);
    
    return nmod_poly_to_algeb(kv,p);

}


/*******************************************************************
 * 
 * Approximant row basis for a vector (m x 1 matrix) of derivatives 
 * 
 *  ALGEB args[1]: shift
 *        args[2]: polynomial (vector of coefficients)  
 *        args[3]: order of approximation 
 *        args[4]: modulus 
 *        args[5]: method, "MBasis" or "PMBasis" 
 * 
 *  Returns ALGEB  res:=[M,dct] 
 *    M: a polynomial matrix, list entries, m x m 
 *    dct: list of out defects  
 * 
 *  Eventually, be careful with the sign between either 
 *      the defect (e.g. in gfun) or shifts = -dct in pml
 * 
 *******************************************************************/


ALGEB pm_diff_mbasis(MKernelVector kv, ALGEB *args){


    //double t = 0.0;
    //clock_t tt;

    mp_limb_t modulus = MapleToInteger64(kv,args[4]);

    nmod_poly_t p;    // The polynomial that will be differentiated 
    nmod_poly_init(p,modulus);
    get_nmod_poly(p, modulus, kv, args[2]);


    //The dimension is deduced from the size of the shift 
    M_INT m = MapleNumArgs(kv,args[1]);


    // the vector of derivatives 
    nmod_poly_mat_t Vdiff;
    nmod_poly_mat_init(Vdiff, m, 1, modulus);  // No use of vectors for the moment, matrices 

    nmod_poly_set(nmod_poly_mat_entry(Vdiff,0,0),p);

    for (ulong i = 1; i < m; i++)
     nmod_poly_derivative(nmod_poly_mat_entry(Vdiff,i,0),nmod_poly_mat_entry(Vdiff,i-1,0));


    // Resulting basis  
    nmod_poly_mat_t M;
    nmod_poly_mat_init(M, m, m, modulus); 

    ALGEB maple_shift = args[1];
    slong shift[m];
    for (ulong i = 0; i < m; i++) 
        shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

    ulong order = MapleToInteger64(kv,args[3]);


    // appromimant computation using PML 

    //tt = clock();


    if (strcmp(MapleToString(kv,args[5]),"MBasis") ==0) 
        nmod_poly_mat_mbasis(M, shift, Vdiff, order); 
    else 
        nmod_poly_mat_pmbasis(M, shift, Vdiff, order); 
   

   //t = (double)(clock()-tt) / CLOCKS_PER_SEC;


   //MapleALGEB_Printf(kv, MapleToString(kv,args[5]));
   //MapleALGEB_Printf(kv, " ++++  Time diffbasis %f ms\n", ToMapleFloat(kv,t*1000));


   // Construction of the result [M, dct]
   ALGEB outdct;

   outdct= MapleListAlloc(kv,m);

   for (slong i = 1; i <= m; i++) 
        MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));
    
   ALGEB res= MapleListAlloc(kv,2);
   MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
   MapleListAssign(kv,res,2,outdct);

   return res;

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
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

 ALGEB pm_matrix_mbasis(MKernelVector kv, ALGEB *args){


   double t = 0.0;
   clock_t tt;
  

   ALGEB vectmat=args[2];

   mp_limb_t modulus = MapleToInteger64(kv,args[4]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, vectmat);

   slong m = A ->r;

   nmod_poly_mat_t M;

   nmod_poly_mat_init(M, m, m, modulus); 

   //slong res_shift[m];

   ALGEB maple_shift = args[1];

   slong shift[m];
   for (ulong i = 0; i < m; i++) 
      shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

   ulong order = MapleToInteger64(kv,args[3]);


   tt = clock();

   nmod_poly_mat_mbasis(M, shift, A, order);
   //mbasis(M, res_shift, A, order, shift);

   t = (double)(clock()-tt) / CLOCKS_PER_SEC;
   MapleALGEB_Printf(kv, " ++++  Time mbasis %f ms\n", ToMapleFloat(kv,t*1000));


   ALGEB outdct;

   outdct= MapleListAlloc(kv,m);

   for (slong i = 1; i <= m; i++) 
        MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));

    
   ALGEB res= MapleListAlloc(kv,2);
   MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
   MapleListAssign(kv,res,2,outdct);

   return res;

}


/**********************************************************
 * 
 * modulo matrix polynomial pmbasis  
 * 
 *  ALGEB args[1]: shift
 *        args[2]: matrix polynomial, vector entries 
 *        args[3]: order
 *        args[4]: modulus 
 * 
 *  Returns M,dct 
 *    M: a polynomial matrix, list entries 
 *    dct: the out defects  
 * 
 *     !!! Be careful with the sign either 
 *         defect (e.g. in gfun) or shifts = -dct in pml
 * 
 * 
 ***********************************************************/


// Row basis, how to specify ? 

 ALGEB pm_matrix_pmbasis(MKernelVector kv, ALGEB *args){


   double t = 0.0;
   clock_t tt;
  

   ALGEB vectmat=args[2];

   mp_limb_t modulus = MapleToInteger64(kv,args[4]);

   nmod_poly_mat_t A;

   get_nmod_poly_mat(A, modulus, kv, vectmat);

   slong m = A ->r;

   nmod_poly_mat_t M;

   nmod_poly_mat_init(M, m, m, modulus); 

   //slong res_shift[m];

   ALGEB maple_shift = args[1];

   slong shift[m];
   for (ulong i = 0; i < m; i++) 
      shift[i]=MapleToInteger64(kv,MapleListSelect(kv,maple_shift,i+1));

   ulong order = MapleToInteger64(kv,args[3]);


   tt = clock();

   nmod_poly_mat_pmbasis(M, shift, A, order);
   //mbasis(M, res_shift, A, order, shift);

   t = (double)(clock()-tt) / CLOCKS_PER_SEC;
   MapleALGEB_Printf(kv, " ++++  Time pmbasis %f ms\n", ToMapleFloat(kv,t*1000));


   ALGEB outdct;

   outdct= MapleListAlloc(kv,m);

   for (slong i = 1; i <= m; i++) 
        MapleListAssign(kv,outdct,i, ToMapleInteger(kv,shift[i-1]));

    
   ALGEB res= MapleListAlloc(kv,2);
   MapleListAssign(kv,res,1,nmod_poly_mat_to_algeb(kv,M));
   MapleListAssign(kv,res,2,outdct);

   return res;

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
