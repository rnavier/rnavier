# Try to find GTS library in this project
#
# It defines the following variables:
#  GTS_FOUND - system has GTS lib
#  GTS_INCLUDE_DIRS - where to find headers 
#  GTS_LIBRARIES - full path to the libraries
#  GTS_LIBRARY_DIRS, the directory where library is found
# 
#  CMAKE_GTS_FLAGS  = Unix compiler flags for GTS, essentially "`gts-config --cxxflags`"
#  GTS_LINK_DIRECTORIES = link directories, useful for rpath on Unix
 
set( GTS_FOUND OFF )
 
set(GTS_CONFIG_EXECUTABLE /bin/gts-config)
#set(GTS_CONFIG_EXECUTABLE ${PROJECT_SOURCE_DIR}/gts/bin/gts-config)


if( GTS_CONFIG_EXECUTABLE ) 
  set( GTS_FOUND ON )
   
  # run the gsl-config program to get cxxflags
  execute_process(
    COMMAND sh "${GTS_CONFIG_EXECUTABLE}" --cflags
    OUTPUT_VARIABLE GTS_CFLAGS
    RESULT_VARIABLE RET
    ERROR_QUIET
    )

  if( RET EQUAL 0 )
    string( STRIP "${GTS_CFLAGS}" GTS_CFLAGS )
    separate_arguments( GTS_CFLAGS )

    # parse definitions from cflags; drop -D* from CFLAGS
    string( REGEX MATCHALL "-D[^;]+"
      GTS_DEFINITIONS  "${GTS_CFLAGS}" )
    string( REGEX REPLACE "-D[^;]+;" ""
      GTS_CFLAGS "${GTS_CFLAGS}" )

    # parse include dirs from cflags; drop -I prefix
    string( REGEX MATCHALL "-I[^;]+"
      GTS_INCLUDE_DIRS "${GTS_CFLAGS}" )
    string( REPLACE "-I" ""
      GTS_INCLUDE_DIRS "${GTS_INCLUDE_DIRS}")
    string( REGEX REPLACE "-I[^;]+;" ""
      GTS_CFLAGS "${GTS_CFLAGS}")

    message("GTS_DEFINITIONS=${GTS_DEFINITIONS}")
    message("GTS_INCLUDE_DIRS=${GTS_INCLUDE_DIRS}")
    message("GTS_CFLAGS=${GTS_CFLAGS}")
  else( RET EQUAL 0 )
    set( GTS_FOUND FALSE )
  endif( RET EQUAL 0 )

  # run the gsl-config program to get the libs
  execute_process(
    COMMAND sh "${GTS_CONFIG_EXECUTABLE}" --libs
    OUTPUT_VARIABLE GTS_LIBRARIES
    RESULT_VARIABLE RET
    ERROR_QUIET
    )
  if( RET EQUAL 0 )
    string(STRIP "${GTS_LIBRARIES}" GTS_LIBRARIES )
    separate_arguments( GTS_LIBRARIES )

    # extract linkdirs (-L) for rpath (i.e., LINK_DIRECTORIES)
    string( REGEX MATCHALL "-L[^;]+"
      GTS_LIBRARY_DIRS "${GTS_LIBRARIES}" )
    string( REPLACE "-L" ""
      GTS_LIBRARY_DIRS "${GTS_LIBRARY_DIRS}" )
  else( RET EQUAL 0 )
    set( GTS_FOUND FALSE )
  endif( RET EQUAL 0 )
   
  MARK_AS_ADVANCED(
    GTS_CFLAGS
  )

else( GTS_CONFIG_EXECUTABLE )
    message( STATUS "FindGTS: gts-config not found.")
endif( GTS_CONFIG_EXECUTABLE )
 
if( GTS_FOUND )
  if( NOT GTS_FIND_QUIETLY )
    message( STATUS "FindGTS: Found both GTS headers and library" )
  endif( NOT GTS_FIND_QUIETLY )
else( GTS_FOUND )
  if( GTS_FIND_REQUIRED )
    message( FATAL_ERROR "FindGTS: Could not find GTS headers or library" )
  endif( GTS_FIND_REQUIRED )
endif( GTS_FOUND )
