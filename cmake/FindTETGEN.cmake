# Find tetgen
#
# Find the tetgen includes and library
# 
# if you need to add a custom library search path, do it via CMAKE_PREFIX_PATH 
# 
# This module defines
#  TETGEN_INCLUDE_DIRS, where to find header, etc.
#  TETGEN_LIBRARIES, the libraries needed to use tetgen.
#  TETGEN_FOUND, If false, do not try to use tetgen.
#  TETGEN_INCLUDE_PREFIX, include prefix for tetgen.

# Additional modules
include(FindPackageHandleStandardArgs)

IF(NOT TETGEN_ROOT_DIR)
 	SET(TETGEN_ROOT_DIR "$ENV{TETGEN_ROOT_DIR}")
ENDIF()

if (WIN32)
	# Find include files
	find_path(
		TETGEN_INCLUDE_DIR
		NAMES tetgen.h
		PATHS
		$ENV{PROGRAMFILES}/include/
		${TETGEN_ROOT_DIR}/
		${TETGEN_ROOT_DIR}/include/
		DOC "The directory where tetgen.h resides")

	# Find library files
	find_library(
		TETGEN_LIBRARY
		NAMES tet
		PATHS
		$ENV{PROGRAMFILES}/lib/
		${TETGEN_ROOT_DIR}/
		${TETGEN_ROOT_DIR}/lib/
		${TETGEN_ROOT_DIR}/src/
        ${TETGEN_ROOT_DIR}/src/Release/
        ${TETGEN_ROOT_DIR}/src/Debug/
        ${TETGEN_ROOT_DIR}/x64/Release/
        ${TETGEN_ROOT_DIR}/x64/Debug/)
else()
	# Find include files
	find_path(
		TETGEN_INCLUDE_DIR
		NAMES tetgen.h
		PATHS
		/usr/include/
		/usr/local/include/
		/sw/include/
		/opt/local/include/
		${TETGEN_ROOT_DIR}/
		${TETGEN_ROOT_DIR}/include/
		DOC "The directory where tetgen.h resides")

	# Find library files
	# Try to use static libraries
	find_library(
		TETGEN_LIBRARY
		NAMES tet
		PATHS
		/usr/lib64/
		/usr/lib/
		/usr/local/lib64/
		/usr/local/lib/
		/sw/lib/
		/opt/local/lib/
		${TETGEN_ROOT_DIR}/
		${TETGEN_ROOT_DIR}/lib/
		${TETGEN_ROOT_DIR}/src/
		DOC "The TET library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(tetgen DEFAULT_MSG TETGEN_INCLUDE_DIR TETGEN_LIBRARY)
if (TETGEN_FOUND)
	set(TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})
	set(TETGEN_LIBRARIES ${TETGEN_LIBRARY})	
endif()

# Hide some variables
mark_as_advanced (TETGEN_INCLUDE_DIR TETGEN_LIBRARY)