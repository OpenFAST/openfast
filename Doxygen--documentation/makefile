REVNUM := `svn info ../.. |grep Revision: | grep -o [0-9]*`


#doxygen bits
default:
	sed s/"PROJECT_NUMBER.*"/"PROJECT_NUMBER = \"Revision: $(REVNUM) (last commit)\""/ -i Doxyfile;
	doxygen;
	$(MAKE) -C latex/;
	cp latex/refman.pdf InflowWind--doxygen.pdf;

