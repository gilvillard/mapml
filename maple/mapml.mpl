


mapml:= module()

	local ModuleLoad, ModuleUnload;

        description "XXX";

	local pmlflintpath, mapmlpath, pmlflint_lib, mapml_lib, ff, 
	      tottest, succeeded;

	

	# Local functions 
	# ---------------

	local convertpol; 

	local matpoly_rt;

	# High level submodule functions 
	# ------------------------------

	export pmCheck_Io, pmDeterminant, pmCheck_Determinant, Dev:



	option package, load=Init, unload=End: 
	
	
	pmlflintpath := "/Users/gvillard/univ-linalg/flint-extras"; 

	mapmlpath := "/usr/local";

	pmlflint_lib := cat(pmlflintpath,"/libpmlflint.dylib");
        mapml_lib := cat(mapmlpath,"/lib/libmapml.dylib");

	ModuleLoad:= proc( )
        printf("Loading mapml\n");
        pmCheck_Io();
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
    
	ff:=define_external("nmod_poly_mat_init", MAPLE, LIB = pmlflint_lib):
		

	# Local functions for conversions
	# -------------------------------


	convertpol:=proc(p,modulus) local d,i,ll; if (degree(p) <0) then d:=-1: else d:=degree(p) fi:   
    	ll:=sprintf("%a %a", d+1,modulus); 
    	for i from 0 to d do ll:=cat(ll," ", sprintf(%a,coeff(p,x,i))); od: 
	return(ll):
	end proc:



	#  Definition of the low level submodule 
	#  -------------------------------------

	Dev:= module()

	    description "YYY";

		export matpoly_rt, pm_determinant;					

		matpoly_rt := define_external("matpoly_rt",MAPLE, LIB = mapml_lib);

		pm_determinant := define_external("pm_determinant",MAPLE, LIB = mapml_lib);

	end module: 	




	#  Definition of the high level submodule 
	#  --------------------------------------

	pmDeterminant :=proc(A,q) local t,stringA; 

		stringA:=map(t->convertpol(t,q),A);

		return Dev:-pm_determinant(stringA,q);

	end proc: 

		

	#  Test functions  
	#  --------------

	tottest:=0: succeeded:=0:


	# export
	pmCheck_Determinant :=proc() local rr,A,val,q,n,d1,d2;

		tottest:=0: succeeded:=0: # GV here or global ? 

		q:=11;
		n:=4; 
		
		rr:=t->randpoly(x,degree=2) mod q:
		A:=RandomMatrix(n,n,generator=rr);
		
		d1:=pmDeterminant(A,q);

		d2:=Determinant(A) mod q;

		val:= expand(d1-d1) mod q;

		if (val=0) then tottest:=tottest+1: succeeded:=succeeded+1: else tottest:=tottest+1: fi:


        printf("%d tests passed, over %d",succeeded,tottest);

		return(val);

	end proc: 

	# local
	matpoly_rt:=proc(A,modulus) local stringA;

   		stringA:=map(t->convertpol(t,modulus),A);

	   return Dev:-matpoly_rt(stringA,modulus);

	end proc:

	# export
	pmCheck_Io :=proc() local rr,t,A,B,Z,val,q,m,n;

		tottest:=0: succeeded:=0: # GV here or global ? 

		q:=11;
		m:=2; n:=4; 
		Z:=Matrix(m,n,0);
		rr:=t->randpoly(x,degree=2) mod q:
		A:=RandomMatrix(m,n,generator=rr);
		B:=matpoly_rt(A,q);

		val:=Equal(map(t-> expand(t) mod q,A-B),Z);

		if (val) then tottest:=tottest+1: succeeded:=succeeded+1: else tottest:=tottest+1: fi:

		q:=2;
		m:=40; n:=44; 
		Z:=Matrix(m,n,0);
		rr:=t->randpoly(x,degree=8) mod q:
		A:=RandomMatrix(m,n,generator=rr);
		B:=matpoly_rt(A,q);

		val:=Equal(map(t-> expand(t) mod q,A-B),Z);

		if (val) then tottest:=tottest+1: succeeded:=succeeded+1: else tottest:=tottest+1: fi:



		q:=179424673;
		m:=20; n:=16; 
		Z:=Matrix(m,n,0);
		rr:=t->randpoly(x,degree=8) mod q:
		A:=RandomMatrix(m,n,generator=rr);
		B:=matpoly_rt(A,q);

		val:=Equal(map(t-> expand(t) mod q,A-B),Z);

		if (val) then tottest:=tottest+1: succeeded:=succeeded+1: else tottest:=tottest+1: fi:



        printf("%d tests passed, over %d",succeeded,tottest);

		return(val);

	end proc: 

end module: 	


