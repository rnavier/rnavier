
add_executable(test_code2d.exe test_code2d.cxx)
target_link_libraries(test_code2d.exe hydro hiextras)

add_executable(test_code3d.exe test_code3d.cxx)
target_link_libraries(test_code3d.exe hydro hiextras)
