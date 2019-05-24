# Find STB
#
# Find the STB includes
# 
# if you need to add a custom library search path, do it via CMAKE_PREFIX_PATH 
# 
# This module defines
#  STB_INCLUDE_DIRS, where to find header, etc.
#  STB_FOUND, If false, do not try to use STB.
#  STB_INCLUDE_PREFIX, include prefix for STB.

# Additional modules
include(FindPackageHandleStandardArgs)

IF(NOT STB_ROOT_DIR)
 	SET(STB_ROOT_DIR "$ENV{STB_ROOT_DIR}")
ENDIF()

if (WIN32)
	# Find include files
	find_path(
		STB_INCLUDE_DIR
		NAMES stb.h
		PATHS
		$ENV{PROGRAMFILES}/include/
		${STB_ROOT_DIR}/
		${STB_ROOT_DIR}/include/
		DOC "The directory where stb.h resides")
else()
	# Find include files
	find_path(
		STB_INCLUDE_DIR
		NAMES stb.h
		PATHS
		/usr/include/
		/usr/local/include/
		/sw/include/
		/opt/local/include/
		${STB_ROOT_DIR}/
		${STB_ROOT_DIR}/include/
		DOC "The directory where stb.h resides")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(stb DEFAULT_MSG STB_INCLUDE_DIR)
if (STB_FOUND)
	set(STB_INCLUDE_DIRS ${STB_INCLUDE_DIR})
endif()

# Hide some variables
mark_as_advanced (STB_INCLUDE_DIR)