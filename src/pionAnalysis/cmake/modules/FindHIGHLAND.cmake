MESSAGE(STATUS "Looking for HIGHLAND...")

# HIGHLAND configured for use from the build tree - absolute paths are used.
set(HIGHLAND_INCLUDE_DIRS $ENV{HIGHLANDPATH}/include)


# HIGHLAND configured for use from the build tree - absolute paths are used.
set(HIGHLAND_LIBRARY_DIR $ENV{HIGHLANDPATH}/lib)


# HIGHLAND configured for use from the build tree - absolute paths are used.
set(HIGHLAND_BINARY_DIR $ENV{HIGHLANDPATH}/bin)


# HIGHLAND configured for use from the build tree - absolute paths are used.
set(HIGHLAND_LIBRARIES "-L$ENV{HIGHLANDPATH}/lib -lhighland")

set(HIGHLAND_FOUND true)

