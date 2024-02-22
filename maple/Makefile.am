
EXTRA_DIST = mapml.mpl mapml.mpl.bak maple-path.sh 

include_mapmldir= $(includedir)/mapml

include_mapml_HEADERS=  mapml.mpl

install-data-local:
	sh maple-path.sh "$(PML_LIBDIR)"  "$(prefix)"
	


