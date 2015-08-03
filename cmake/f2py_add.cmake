# Compile a module based on sources:
#
# Method 1: Just sources
#
#   F2PY_ADD(MODULE foo SRC file1.f90 file2.f90 ...)
#
# Method 2: Just fortran sources and exernal libraries
#
#   F2PY_ADD(MODULE foo SRC foo1.f90 foo2.f90 ... 
#       [EXTERNAL_LIBS externlibs])
#
#   if EXTERNAL_LIBS is defined, this will be passed on to the f2py command
#   line in the appropriate place. externlibs might be for example:
#  
#    -DF2PY_REPORT_ON_ARRAY_COPY=1 -L/opt/local/lib -lgsl -lgslcblas
#
# Method 3:  (Recommended) Link with a shared library
#
#    F2PY_ADD(MODULE foo 
#       SRC file1.f90  file2.f90 ....
#       LIBS lib1 lib2  ...
#       [DONT_COMPILE_SRC] 
#       [EXTERNAL_LIBS externlibs] )
#
#    Constructs and compiles the wrapper code, and links with lib1, lib2 etc.
#    This does NOT actually compile the fortran code it is assumed that the
#    necessary objects are in lib1 lib2 etc
#     
#    (Recommended) If the [DONT_COMPILE_SRC] flag is present, then it is
#    assumed that the fortran objects are in the libraries and should not be
#    compiled by f2py.  This is actually the recomended way to build python
#    modules. Then the signature file is created (unless the SIGNATURE flag is
#    supplied), and the c-wrapper is compiled and linked against the
#    libraries.
#
# Custom Signature Files:
#
#    In order to  custom signature file (as opposed to the automatic one)
#    use the SIGNATURE flag. Then call for example
# 
#    F2PY_ADD(MODULE foo 
#        SIGNATURE SRC myfoo.pyf file1.f90 file2.f90 ...
#        [EXTERNAL_LIBS] ${foobar})
#  
#    Other valid calling sequences are for example
# 
#    F2PY_ADD(MODULE foo 
#        SIGNATURE SRC myfoo.pyf 
#        LIB lib1 DONT_COMPILE_SRC)
#
# Installation and DESTINATION:
#
#    All of the forms of the command can add a DESTINATION directory,
# 
#    F2PY_ADD(MODULE foo SRC bar.f90 DESTINATION lib)
#
function(F2PY_ADD)

  include(cmake_parse_arguments)
   
  cmake_parse_arguments(F2PY_ADD "DONT_COMPILE_SRC;SIGNATURE" "MODULE;DESTINATION" "SRC;LIBS;EXTERNAL_LIBS" ${ARGN})

  # check the arguements:
  # - check that MODULE is defined
  # - check that SRC is defined
  # - check that if DONT_COMPILE_SRC and SIGNATURE is defined that SRC 
  if (NOT F2PY_ADD_MODULE) 
     message(FATAL_ERROR "[F2PY] Undefined module in f2py_add. Add the module name to the arguements f2py_add(MODULE foo ..)")
  endif()
  if (NOT F2PY_ADD_SRC) 
     message(FATAL_ERROR "[F2PY] Undefined source in f2py_add. Add the source to the arguments f2py_add(... SRC foo.f90..)")
  endif()
  if (F2PY_ADD_DONT_COMPILE_SRC AND NOT F2PY_ADD_LIBS)
     message(WARNING "[F2PY] f2py_add is called with DONT_COMPILE_SRC but with out any LIBS")
  endif()

  # give information
  message(STATUS "Building f2py module ${F2PY_ADD_MODULE} with ${F2PY_ADD_SRC}")
  if (F2PY_ADD_EXTERNAL_LIBS)
     message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with external lib ${F2PY_ADD_EXTERNAL_LIBS}")
  endif()
  if (F2PY_ADD_LIBS)
     message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with libraries ${F2PY_ADD_LIBS}")
  endif()
  if (F2PY_ADD_SIGNATURE)
     message(STATUS "F2PY: Building f2py module ${F2PY_ADD_MODULE} with custom signature")
  endif()
   
  # Get a list of the include directories.
  # The f2py --include-paths option used when generating a signature file
  get_directory_property(_inc_dirs INCLUDE_DIRECTORIES)
  string(REPLACE ";" ":" _inc_paths "${_inc_dirs}")

  # Make the source filenames absolute. The result is in _abs_srcs
  set(_abs_srcs)
  foreach(_src ${F2PY_ADD_SRC})
    get_filename_component(_abs_src ${_src} ABSOLUTE)
    list(APPEND _abs_srcs ${_abs_src})
  endforeach(_src ${F2PY_ADD_SRC})

  # Make a line like -llib1 -l1ib2 -llib3 etc. The result is in _flibs
  if (F2PY_ADD_LIBS) 
     set(_flibs)
     foreach(_lib ${F2PY_ADD_LIBS}) 
        list(APPEND _flibs -l${_lib})
        get_target_property(_liblocation ${_lib} LOCATION)
        get_filename_component(_libpath ${_liblocation} PATH)
        find_library(_libpath ${_lib})
        message("Found the library path ${_liblocation}")
        message("Found the path ${_libpath}")
     endforeach()
  endif()

  # A custom signature file has not been provided. Generate signature file
  # module.pyf using f2py -m module fortran source
  if (NOT F2PY_ADD_SIGNATURE) 
     add_custom_command(
        OUTPUT ${F2PY_ADD_MODULE}.pyf
        COMMAND f2py --quiet --include-paths ${_inc_paths}  -m ${F2PY_ADD_MODULE} --overwrite-signature -h ${F2PY_ADD_MODULE}.pyf ${_abs_srcs} 
        DEPENDS ${F2PY_ADD_SRC}
        COMMENT "[F2PY] Generating fortran signature files ${F2PY_ADD_MODULE}.pyf"
     )
     list(INSERT 0 ${F2PY_ADD_MODULE}.pyf ${_abs_srcs})
  endif()

  # If a signature file, module.pyf, is included
  # The commands are identical, but the DEPENDS are different. 
  # If a signature is provided then the DEPENDS is just the source
  # otherwise the DEPENDS is just the source and signature.pyf
  if (F2PY_ADD_SIGNATURE) 
     set(_depends ${F2PY_ADD_SRC})
  else()
     set(_depends ${F2PY_ADD_MODULE}.pyf ${F2PY_ADD_SRC})
  endif()

  add_custom_command(
      OUTPUT ${F2PY_ADD_MODULE}.so
      COMMAND f2py --quiet --include-paths ${_inc_paths} -m ${F2PY_ADD_MODULE} -c ${_abs_srcs} -L${CMAKE_CURRENT_BINARY_DIR} ${_flibs} ${F2PY_ADD_EXTERNAL_LIBS} 
      DEPENDS ${_depends}
      COMMENT "[F2PY] Building Fortran to Python interface module ${F2PY_ADD_MODULE}"
  )

  # Add a custom target <name> to trigger the generation of the python module.
  add_custom_target(${F2PY_ADD_MODULE} ALL DEPENDS ${F2PY_ADD_MODULE}.so)

  # Install the python module
#  if (F2PY_ADD_DESTINATION)
#     install(TARGETS ${F2PY_ADD_MODULE}.so RUNTIME DESTINATION ${F2PY_ADD_DESTINATION})
#  endif()

endfunction(F2PY_ADD)




### -----------------------------------------------------------------------------
### Macro to generate a Python interface module from one or more Fortran sources
###
### Usage: add_f2py_module(<module-name> <src1>..<srcN> DESTINATION <install-dir>
###
#macro (add_f2py_module _name)

#  # Precondition check.
#  if(NOT F2PY_EXECUTABLE)
#    message(FATAL_ERROR "add_f2py_module: f2py executable is not available!")
#  endif(NOT F2PY_EXECUTABLE)

#  # Parse arguments.
#  string(REGEX REPLACE ";?DESTINATION.*" "" _srcs "${ARGN}")
#  string(REGEX MATCH "DESTINATION;.*" _dest_dir "${ARGN}")
#  string(REGEX REPLACE "^DESTINATION;" "" _dest_dir "${_dest_dir}")

#  # Sanity checks.
#  if(_srcs MATCHES "^$")
#    message(FATAL_ERROR "add_f2py_module: no source files specified")
#  endif(_srcs MATCHES "^$")
#  if(_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";")
#    message(FATAL_ERROR "add_f2py_module: destination directory invalid")
#  endif(_dest_dir MATCHES "^$" OR _dest_dir MATCHES ";")

#  # Get the compiler-id and map it to compiler vendor as used by f2py.
#  # Currently, we only check for GNU, but this can easily be extended. 
#  # Cache the result, so that we only need to check once.
#  if(NOT F2PY_FCOMPILER)
#    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#      if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#        set(_fcompiler "gnu95")
#      else(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#        set(_fcompiler "gnu")
#      endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    else(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#      set(_fcompiler "F2PY_FCOMPILER-NOTFOUND")
#    endif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(F2PY_FCOMPILER ${_fcompiler} CACHE STRING
#      "F2PY: Fortran compiler type by vendor" FORCE)
#    if(NOT F2PY_FCOMPILER)
#      message(STATUS "[F2PY]: Could not determine Fortran compiler type. "
#                     "Troubles ahead!")
#    endif(NOT F2PY_FCOMPILER)
#  endif(NOT F2PY_FCOMPILER)

#  # Set f2py compiler options: compiler vendor and path to Fortran77/90 compiler.
#  if(F2PY_FCOMPILER)
#    set(_fcompiler_opts "--fcompiler=${F2PY_FCOMPILER}")
#    list(APPEND _fcompiler_opts "--f77exec=${CMAKE_Fortran_COMPILER}")
#    if(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#      list(APPEND _fcompiler_opts "--f90exec=${CMAKE_Fortran_COMPILER}")
#    endif(CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#  endif(F2PY_FCOMPILER)

#  # Make the source filenames absolute.
#  set(_abs_srcs)
#  foreach(_src ${_srcs})
#    get_filename_component(_abs_src ${_src} ABSOLUTE)
#    list(APPEND _abs_srcs ${_abs_src})
#  endforeach(_src ${_srcs})

#  # Get a list of the include directories.
#  # The f2py --include_paths option, used when generating a signature file,
#  # needs a colon-separated list. The f2py -I option, used when compiling
#  # the sources, must be repeated for every include directory.
#  get_directory_property(_inc_dirs INCLUDE_DIRECTORIES)
#  string(REPLACE ";" ":" _inc_paths "${_inc_dirs}")
#  set(_inc_opts)
#  foreach(_dir ${_inc_dirs})
#    list(APPEND _inc_opts "-I${_dir}")
#  endforeach(_dir)

#  # Define the command to generate the Fortran to Python interface module. The
#  # output will be a shared library that can be imported by python.
#  add_custom_command(OUTPUT ${_name}.so
#    COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} -h ${_name}.pyf
#            --include_paths ${_inc_paths} --overwrite-signature ${_abs_srcs}
#    COMMAND ${F2PY_EXECUTABLE} --quiet -m ${_name} -c ${_name}.pyf
#            ${_fcompiler_opts} ${_inc_opts} ${_abs_srcs}
#    DEPENDS ${_srcs}
#    COMMENT "[F2PY] Building Fortran to Python interface module ${_name}")

#  # Add a custom target <name> to trigger the generation of the python module.
#  add_custom_target(${_name} ALL DEPENDS ${_name}.so)

#  # Install the python module
#  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_name}.so
#    DESTINATION ${_dest_dir})

#endmacro (add_f2py_module)

