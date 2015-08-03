add_executable(test_icbj_hydro.exe test_icbj_hydro.cxx)
target_link_libraries(test_icbj_hydro.exe hydro)

add_executable(test_icbj.exe test_icbj.cxx)
target_link_libraries(test_icbj.exe hydro)

add_executable(test_gaus.exe test_gaus.cxx)
target_link_libraries(test_gaus.exe hydro)

add_executable(test_solveSource.exe test_solveSource.cxx)
target_link_libraries(test_solveSource.exe hydro)

add_executable(test_ndarray.exe test_ndarray.cxx)
