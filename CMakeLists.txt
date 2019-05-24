CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Name of the project
PROJECT(RedMax)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

################################################################################
### Custom setting for lab machine ###
IF(WIN32)
  SET(GLM_ROOT_DIR "D:\\Libraries\\glm-0.9.8.5\\glm")
  SET(GLFW_ROOT_DIR "D:\\Libraries\\glfw-3.2.1")
  SET(GLEW_ROOT_DIR "D:\\Libraries\\glew-2.1.0")
  SET(EIGEN3_ROOT_DIR "D:\\Libraries\\eigen-eigen-6c3592226be0")
  SET(TETGEN_ROOT_DIR "D:\\Libraries\\tetgen1.5.1-beta1")
  SET(MOSEK_ROOT_DIR "D:\\Libraries\\mosek_8.1_x64\\8")
  SET(PARDISO_ROOT_DIR "D:\\Libraries\\PARDISO")
  SET(MKL_ROOT_DIR "C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries_2018\\windows\\mkl")
  SET(JSONCPP_ROOT_DIR "D:\\Libraries\\jsoncpp\\jsoncpp_x64-windows")
  SET(STB_ROOT_DIR "D:\\Libraries\\stb-master")
  SET(LIBIGL_ROOT_DIR "D:\\Libraries\\libigl-master\\libigl")
#  SET(DEF_DIR_CMAES "D:\\Libraries\\CMA-ESpp-master")
#  SET(DEF_DIR_LIBLBFGS "D:\\Libraries\\liblbfgs-master")
#  SET(DEF_DIR_GUROBI "D:\\Libraries\\gurobi\\win64")
#  SET(DEF_DIR_YIXUAN "D:\\Libraries\\yixuan-LBFGSpp-ad499fe")
#  SET(DEF_DIR_DLIB "D:\\Libraries\\dlib-19.13")
#  SET(DEF_DIR_QUADPROGPP "D:\\Libraries\\quadprogpp")
#  SET(DEF_DIR_NLOPT "D:\\Libraries\\nlopt")
ENDIF()

################################################################################
### Custom setting for mac ###
IF(APPLE)
  SET(GLM_ROOT_DIR "/Users/yingwang/Libraries/glm-0.9.9.3")
  SET(GLFW_ROOT_DIR "/Users/yingwang/Libraries/glfw-master")
  SET(GLEW_ROOT_DIR "/Users/yingwang/Libraries/glew-2.1.0")
  SET(EIGEN3_ROOT_DIR "/Users/yingwang/Libraries/eigen-eigen-08a1f23d4f46")
  SET(TETGEN_ROOT_DIR "/Users/yingwang/Libraries/tetgen1.5.1")
  SET(MOSEK_ROOT_DIR "/Users/yingwang/Libraries/mosek/8")

  SET(PARDISO_ROOT_DIR "/Users/yingwang/Libraries/pardiso")
  SET(MKL_ROOT_DIR "/opt/intel/mkl")
  SET(JSONCPP_ROOT_DIR "/Users/yingwang/Libraries/jsoncpp-0.10.7")
  SET(NLOHMANN_ROOT_DIR "/Users/yingwang/Libraries/json-3.5.0")
  SET(STB_ROOT_DIR "/Users/yingwang/Libraries/stb-master")
  SET(LIBIGL_ROOT_DIR "/Users/yingwang/Libraries/libigl-2.0.0")
ENDIF()

################################################################################
### Available options ###

# Is this the solution?
# Override with `cmake -DSOL=ON ..`
OPTION(SOL "Solution" OFF)
OPTION(REDMAX_WITH_MOSEK    "Use MOSEK"    ON)
OPTION(REDMAX_WITH_PARDISO  "Use PARDISO"  ON)
OPTION(REDMAX_WITH_MKL      "Use MKL"      ON)
OPTION(REDMAX_WITH_JSONCPP  "Use JSONCPP"  ON)
OPTION(REDMAX_WITH_NLOHMANN "Use NlOHMANN" ON)
OPTION(REDMAX_WITH_STB      "Use STB"      ON)

################################################################################

# Use glob to get the list of all source files.
# We don't really need to include header and resource files to build, but it's
# nice to have them also show up in IDEs.
IF(${SOL})
  FILE(GLOB_RECURSE SOURCES "src0/*.cpp")
  FILE(GLOB_RECURSE HEADERS "src0/*.h")
ELSE()
  FILE(GLOB_RECURSE SOURCES "src/*.cpp")
  FILE(GLOB_RECURSE HEADERS "src/*.h")
ENDIF()
FILE(GLOB_RECURSE GLSL "resources/*.glsl")

# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

################################################################################
### Compile the Eigen3 part ###
find_package (Eigen3 3.3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

################################################################################
### Compile the GLM part ###
find_package (GLM REQUIRED)
INCLUDE_DIRECTORIES(${GLM_INCLUDE_DIRS})

################################################################################
### Compile the GLFW part ###
find_package (GLFW REQUIRED)
INCLUDE_DIRECTORIES(${GLFW_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLFW_LIBRARIES})

################################################################################
### Compile the GLEW part ###
find_package (GLEW REQUIRED)
INCLUDE_DIRECTORIES(${GLEW_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_LIBRARIES})

################################################################################
### Compile the TETGEN part ###
find_package (TETGEN REQUIRED)
INCLUDE_DIRECTORIES(${TETGEN_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${TETGEN_LIBRARIES})

################################################################################
### Compile the JSONCPP part ###
IF(REDMAX_WITH_JSONCPP)
  find_package (JSONCPP REQUIRED)
  INCLUDE_DIRECTORIES(${JSONCPP_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${JSONCPP_LIBRARIES})
  ADD_DEFINITIONS(-DREDMAX_JSONCPP)
ENDIF()

################################################################################
### Compile the MKL part ###
IF(REDMAX_WITH_MKL)
  find_package (MKL REQUIRED)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MKL_LIBRARIES})
  ADD_DEFINITIONS(-DREDMAX_MKL)
ENDIF()

################################################################################
### Compile the PARDISO part ###
IF(REDMAX_WITH_PARDISO)
  find_package (PARDISO REQUIRED)
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${PARDISO_LIBRARIES})
  ADD_DEFINITIONS(-DREDMAX_PARDISO)
ENDIF()

################################################################################
### Compile the NLOHMANN part ###
IF(REDMAX_WITH_NLOHMANN)
  find_package (NLOHMANN REQUIRED)
  INCLUDE_DIRECTORIES(${NLOHMANN_INCLUDE_DIRS})
  ADD_DEFINITIONS(-DREDMAX_NLOHMANN)
ENDIF()

################################################################################
### Compile the MOSEK part ###
IF(REDMAX_WITH_MOSEK)
  find_package(MOSEK REQUIRED)
  INCLUDE_DIRECTORIES(${MOSEK_DIRS})
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MOSEK_LIBRARIES})
  ADD_DEFINITIONS(-DREDMAX_MOSEK)
ENDIF()

################################################################################
### Compile the STB part ###
IF(REDMAX_WITH_STB)
  find_package(STB REQUIRED)
  INCLUDE_DIRECTORIES(${STB_DIRS})
  ADD_DEFINITIONS(-DREDMAX_STB)
ENDIF()

################################################################################
### OS specific options and libraries ###
IF(WIN32)
  # c++11 is enabled by default.
  # -Wall produces way too many warnings.
  # -pedantic is not supported.
  # Disable warning 4996.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} opengl32.lib)
ELSE()
# Enable all pedantic warnings.
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
    # Add required frameworks for GLFW.
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -pthread")
    #Link the Linux OpenGL library
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "GL")
  ENDIF()
ENDIF()
