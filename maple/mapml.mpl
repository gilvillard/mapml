


mapml:= module()

	local ModuleLoad, ModuleUnload;

    description "XXX";

	local flint_path, mapml_path,ff:

	export pm_check_io,Dev:

	local convertpol,matpoly_rt;


	option package, load=Init, unload=End: 
	
	
	flint_path := "/usr/local/lib/libpmlflint.dylib";
	mapml_path := "/usr/local/lib/libmapml.dylib";

	ModuleLoad:= proc( )
        printf("Loading mapml\n");
    end proc:

    # Explicitly call ModuleLoad here so the type is registered when this
    # code is cut&pasted in.  ModuleLoad gets called when the module is
    # read from the repository, but the code used to define a module
    # (like the command below) is not called at library read time.
    ModuleLoad();

    #ModuleUnload:= proc( )
    #printf("");
    #end proc:

    # GV rather put in ModuleLoad ? 
    
	ff:=define_external("nmod_poly_mat_init", MAPLE, LIB = flint_path):
		

	convertpol:=proc(p,modulus) local d,i,ll; if (degree(p) <0) then d:=-1: else d:=degree(p) fi:   
    	ll:=sprintf("%a %a", d+1,modulus); 
    	for i from 0 to d do ll:=cat(ll," ", sprintf(%a,coeff(p,x,i))); od: 
	return(ll):
	end proc:



	#  Definition of the low level submodule 
	#  -------------------------------------

	Dev:= module()

	    description "YYY";

		export matpoly_rt;					

		matpoly_rt := define_external("matpoly_rt",MAPLE, LIB = mapml_path);


	end module: 	


	#  Definition of the high level submodule 
	#  --------------------------------------




	#  Test functions  
	#  --------------

	matpoly_rt:=proc(A,modulus) local stringA;

   		stringA:=map(t->convertpol(t,modulus),A);

	   return Dev:-matpoly_rt(stringA,modulus);

	end proc:


	pm_check_io :=proc() local rr,t,A,B,Z,val,q;

		q:=11;
		m:=2; n:=4; 
		Z:=Matrix(m,n,0);
		rr:=t->randpoly(x,degree=2) mod q:
		A:=RandomMatrix(m,n,generator=rr);
		B:=matpoly_rt(A,q);

		val:=Equal(map(t-> expand(t) mod q,A-B),Z)
		return();



	end proc: 

end module: 	


