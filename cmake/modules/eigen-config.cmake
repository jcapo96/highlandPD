# Locates the eigen directories.

IF(EXISTS /cvmfs/larsoft.opensciencegrid.org/)
          set(EIGEN_INCLUDE_DIRS /cvmfs/larsoft.opensciencegrid.org/products/tensorflow/v1_12_0c/Linux64bit+2.6-2.12-e19-py2-prof/include/eigen/)
endif()