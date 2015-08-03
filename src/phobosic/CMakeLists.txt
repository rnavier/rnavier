ROOT_GENERATE_DICTIONARY(
    ${CMAKE_CURRENT_BINARY_DIR}/phobos_dict
    TPhobosMC.h
    TGlbEvent.h
    LINKDEF phobosic_linkdef.h
    )

set(phobosic_src 
    ${CMAKE_CURRENT_BINARY_DIR}/phobos_dict.cxx
    TPhobosMC.cxx
    TICPhobosMC.cxx
    TGlbEvent.cxx
    )

add_library(phobosic OBJECT ${phobosic_src})
