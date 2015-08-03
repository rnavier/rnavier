include_directories(makecode)

set(hydro_source
    TBRSSSStep.cxx
    TBRSSSSolverS2.cxx
    TIOBuffer.cxx
    TEOS.cxx
    TBRSSS.cxx
    TRNavier3DBj.cxx
    RNavierPath.cxx
    TFOMaker.cxx
    TDomain3D.cxx
    TGrid3D.cxx
    TIC3D.cxx
    TICBjTest.cxx
    THydro3DBj.cxx
    THYAnalysis.cxx
    TBjRegulator.cxx
    )

add_library(hydro ${hydro_source}) 
target_link_libraries(hydro ${GSL_LIBRARIES})
