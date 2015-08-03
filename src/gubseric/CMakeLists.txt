

add_executable(test_icgb_hydro.exe test_icgb_hydro.cxx TICGbTest.cxx)
add_executable(test_fo.exe test_fo.cxx TICGbTest.cxx)
target_link_libraries(test_icgb_hydro.exe hydro) 
target_link_libraries(test_fo.exe hydro) 
