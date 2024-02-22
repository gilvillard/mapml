
EXTRA_DIST = mapml.mpl mapml.mpl.bak maple-path.sh 

include_mapmldir= $(includedir)/mapml


install-data-local:
	sh maple-path.sh "$(PML_LIBDIR)"  "$(prefix)"
	$(INSTALL_DATA) mapml.mpl $(includedir)/mapml/mapml.mpl
	


