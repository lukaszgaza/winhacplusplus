
aux_source_directory(. WINHAC_SOURCES)

add_library(winhac SHARED ${WINHAC_SOURCES} )
target_link_libraries(winhac ${LHAPDF} ${LIBS})
IF(ROOT_FOUND)
	target_link_libraries(winhac ${ROOT_LIBRARIES})
ENDIF (ROOT_FOUND)


install(TARGETS winhac
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)


