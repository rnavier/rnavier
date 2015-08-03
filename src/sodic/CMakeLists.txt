
set(src
    TEOSIdeal.cxx
    TICSodTest.cxx)

add_executable(test_sod.exe test_sod.cxx ${src})
target_link_libraries(test_sod.exe hydro) 

add_executable(test_sod_2d.exe test_sod_2d.cxx ${src})
target_link_libraries(test_sod_2d.exe hydro) 
