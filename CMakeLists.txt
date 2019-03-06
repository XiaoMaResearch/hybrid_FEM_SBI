cmake_minimum_required(VERSION 2.8)



project(hybrid_fem_bie)
enable_language(CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${hybrid_fem_bie_SOURCE_DIR}/cmake/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${hybrid_fem_bie_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "-std=c++0x" CACHE STRING "" FORCE)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
mark_as_advanced(BUILD_SHARED_LIBS)

#set(CMAKE_CXX_COMPILER clang-omp++ CACHE STRING "C++ compiler" FORCE)
find_package(GSL REQUIRED)
find_package(FFTW REQUIRED)
include_directories("/usr/local/include")
include_directories("src")
include_directories("src/bie")
include_directories("src/bie/interface_laws")
include_directories("src/fem")
include_directories("src/tools")

include_directories("~/Cpp_Libarary/libigl/include")
include_directories("~/Cpp_Libarary/Eigen")
include_directories("~/Cpp_Libarary/Eigen/unsupported")
include_directories(${GSL_INCLUDE_PATH})
include_directories(${FFTW_INCLUDE_PATH})
# include exodus directories

add_subdirectory(src)


# ------------------------------------------------------
# tests
option(WI_TESTS "tests" ON)

if(WI_TESTS)
  add_subdirectory(tests)
endif()

# junk
option(WI_JUNK "junk" ON)

if(WI_JUNK)
  add_subdirectory(junk)
endif()
# parallel _simualtion

# ------------------------------------------------------
# add simulation folder
# this folder does not exist
# so that you can add a version controlled simulation folder (similar to tests)
option(WI_SIMULATIONS "simulations" ON)
#
if (WI_SIMULATIONS)
   add_subdirectory(simulations)
 endif()
