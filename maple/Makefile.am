
EXTRA_DIST = mapml.mpl mapml.mpl.bak maple-path.sh 

include_mapmldir= $(includedir)/mapml

# GV how to make mapml.mpl depend on mapml.mpl.bak ? Did not succeed for that 

install-data-local:
	sh maple-path.sh "$(PML_LIBDIR)"  "$(prefix)"
	$(INSTALL_DATA) mapml.mpl $(includedir)/mapml/mapml.mpl
	


