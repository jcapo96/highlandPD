# Check the pressence of cvmfs to enable compilation of tensorflow and related libraries

IF(EXISTS /cvmfs/larsoft.opensciencegrid.org/)
	  set(COMPILETF 1)
endif()
