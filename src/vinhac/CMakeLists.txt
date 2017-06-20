
aux_source_directory(. WINHAC_SOURCES)
aux_source_directory(./interfaces WINHAC_SOURCES)
aux_source_directory(./pdf WINHAC_SOURCES)
aux_source_directory(./core WINHAC_SOURCES)
aux_source_directory(./demo WINHAC_SOURCES)
aux_source_directory(./handlers WINHAC_SOURCES)
aux_source_directory(./input WINHAC_SOURCES)
aux_source_directory(./model/ WINHAC_SOURCES)
aux_source_directory(./util/ WINHAC_SOURCES)
aux_source_directory(./model/amplitudes WINHAC_SOURCES)
aux_source_directory(./model/util WINHAC_SOURCES)

add_library(winhac SHARED ${WINHAC_SOURCES} )
target_link_libraries(winhac ${LHAPDF_LIB} ${LIBS})
IF(ROOT_FOUND)
	target_link_libraries(winhac ${ROOT_LIBRARIES})
ENDIF (ROOT_FOUND)

IF(PYTHIA8_FOUND AND BOOST_FOUND)
	target_link_libraries(winhac ${PYTHIA8_LIB} ${BOOST_LIB} )
ENDIF(PYTHIA8_FOUND AND BOOST_FOUND)


install(TARGETS winhac
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)


