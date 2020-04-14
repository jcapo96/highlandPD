#---To be used from the binary tree--------------------------------------------------------------------------
set(HIGHLAND_PD_INCLUDE_DIR_SETUP "
# HIGHLAND_PD configured for use from the build tree - absolute paths are used.
set(HIGHLAND_PD_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)
")
set(HIGHLAND_PD_LIBRARY_DIR_SETUP "
# HIGHLAND_PD configured for use from the build tree - absolute paths are used.
set(HIGHLAND_PD_LIBRARY_DIR ${CMAKE_SOURCE_DIR}/lib)
")
set(HIGHLAND_PD_BINARY_DIR_SETUP "
# HIGHLAND_PD configured for use from the build tree - absolute paths are used.
set(HIGHLAND_PD_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
")
set(HIGHLAND_PD_LIBRARIES_SETUP "
# HIGHLAND_PD configured for use from the build tree - absolute paths are used.
set(HIGHLAND_PD_LIBRARIES \"-L${CMAKE_SOURCE_DIR}/lib -lhighlandPD\")
")

configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/HIGHLAND_PDConfig.cmake.in
               ${CMAKE_SOURCE_DIR}/HIGHLAND_PDConfig.cmake @ONLY NEWLINE_STYLE UNIX)

