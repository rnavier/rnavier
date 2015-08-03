

set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

# defines target library hydro
include_directories(hydro)
add_subdirectory(hydro)

#defines target object library eoss95p
if (EOSS95P)
    include_directories(eoss95p)
    add_subdirectory(eoss95p)
    set(objs ${objs} $<TARGET_OBJECTS:eoss95p>)
endif()

#defines target object library for phobosic 
if (PHOBOSIC)
    include_directories(phobosic)
    add_subdirectory(phobosic)
    set(objs ${objs} $<TARGET_OBJECTS:phobosic>)
endif()

# put together the hiextras library
add_library(hiextras ${objs})
if (RNAVIER_ROOT)
    target_link_libraries(hiextras ${ROOT_LIBRARIES})
endif()
target_link_libraries(hiextras hydro)

# Start of applications codes
add_subdirectory(test_core)
add_subdirectory(test_code2d)
add_subdirectory(sodic)
add_subdirectory(bjperturbation)
add_subdirectory(gubseric)




