# Locates the protobuf library and include directories.

if(DEFINED ENV{COMPILETF})
          set(PROTOBUF_LIBRARIES /cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_11_2a/Linux64bit+2.6-2.12-e19/lib64/libprotobuf.so.3.11.2.0)
          set(PROTOBUF_INCLUDE_DIRS /cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_11_2a/Linux64bit+2.6-2.12-e19/include/)
endif()