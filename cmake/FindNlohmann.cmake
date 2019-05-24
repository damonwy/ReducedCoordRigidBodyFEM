# Find nlohmannjson
#
# Find the nlohmannjson includes
# 
# if you need to add a custom library search path, do it via CMAKE_PREFIX_PATH 
# 
# This module defines
#  NLOHMANN_INCLUDE_DIRS, where to find header, etc.
#  NLOHMANN_FOUND, If false, do not try to use nlohmannjson.
#  NLOHMANN_INCLUDE_PREFIX, include prefix for nlohmannjson.

# Additional modules
include(FindPackageHandleStandardArgs)

IF(NOT NLOHMANN_ROOT_DIR)
 	SET(NLOHMANN_ROOT_DIR "$ENV{NLOHMANN_ROOT_DIR}")
ENDIF()

if (WIN32)
	# Find include files
	find_path(
		NLOHMANN_INCLUDE_DIR
		NAMES nlohmann/json.hpp
		PATHS
		$ENV{PROGRAMFILES}/include/
		${NLOHMANN_ROOT_DIR}/
		${NLOHMANN_ROOT_DIR}/include/
		DOC "The directory where nlohmann/json.h resides")
else()
	# Find include files
	find_path(
		NLOHMANN_INCLUDE_DIR
		NAMES nlohmann/json.hpp
		PATHS
		/usr/include/
		/usr/local/include/
		/sw/include/
		/opt/local/include/
		${NLOHMANN_ROOT_DIR}/
		${NLOHMANN_ROOT_DIR}/include/
		DOC "The directory where nlohmann/json.h resides")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(nlohmann DEFAULT_MSG NLOHMANN_INCLUDE_DIR)
if (NLOHMANN_FOUND)
	set(NLOHMANN_INCLUDE_DIRS ${NLOHMANN_INCLUDE_DIR})
endif()

# Hide some variables
mark_as_advanced (NLOHMANN_INCLUDE_DIR)
