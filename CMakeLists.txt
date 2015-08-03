CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(rnavier CXX)

set(CMAKE_MODULE_PATH ${rnavier_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

#------------------------------------------------------------------------
# Turn on GSL
#------------------------------------------------------------------------
find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})

#------------------------------------------------------------------------
# Turn on boost and ndarray
#------------------------------------------------------------------------
find_package(Boost)
if (${Boost_FOUND}) 
    message("Found Boost Headers in ${Boost_INCLUDE_DIRS}")
else()
    message("Could not find Boost headers")
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(ndarray)

#------------------------------------------------------------------------
# Turn on root
#------------------------------------------------------------------------
option(RNAVIER_ROOT "Link with ROOT Librariers" ON)
if (RNAVIER_ROOT)
   message("RNAVIER_ROOT is on") 
   add_definitions("-DRNAVIER_ROOT") 
   find_package(ROOT)
   include_directories(AFTER ${ROOT_INCLUDE_DIR})
   message("RNAVIER_ROOT links with ${ROOT_LIBRARIES}")
endif()


#------------------------------------------------------------------------
# Set options
#------------------------------------------------------------------------
include(CMakeDependentOption)

# These options are switched on by default
option (EOSS95P "Compile the EOS s95p" ON)

cmake_dependent_option(PHOBOSIC "Add source code to include phobos glauber" ON RNAVIER_ROOT OFF)

add_subdirectory(src)

