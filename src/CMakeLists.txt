
add_executable (Demo demo/Demo.cpp)
target_link_libraries(Demo ${LIBS} winhac winhacUtils winhacUtilsSanc)

IF(PYTHIA8_FOUND AND BOOST_FOUND)
	include_directories(${PYTHIA8_INCLUDE})
	include_directories(${BOOST_INCLUDE})
	add_executable (DemoPythia demo/DemoPythia.cpp)
	target_link_libraries(DemoPythia ${LIBS} ${PYTHIA8_LIB} ${BOOST_LIB} winhac winhacUtils winhacUtilsSanc)
	add_custom_target(runDemoPythia ${EXECUTABLE_OUTPUT_PATH}/DemoPythia ${PYTHIA8_PREFIX}/xmldoc ../input ../output 10 ../input/UserFile.xml | tee runDemoPythia.log)
	add_custom_target(runDemoPythia1000 ${EXECUTABLE_OUTPUT_PATH}/DemoPythia ${PYTHIA8_PREFIX}/xmldoc ../input ../output 1000 ../input/UserFile.xml | tee runDemoPythia.log)
	add_custom_target(runDemoPythia1000000 ${EXECUTABLE_OUTPUT_PATH}/DemoPythia ${PYTHIA8_PREFIX}/xmldoc ../input ../output 1000000 ../input/UserFile.xml | tee runDemoPythia.log)
	add_custom_target(runDemoPythia10000 ${EXECUTABLE_OUTPUT_PATH}/DemoPythia ${PYTHIA8_PREFIX}/xmldoc ../input ../output 10000 ../input/UserFile.xml | tee runDemoPythia.log)
	add_custom_target(runDemoPythia100000 ${EXECUTABLE_OUTPUT_PATH}/DemoPythia ${PYTHIA8_PREFIX}/xmldoc ../input ../output 100000 ../input/UserFile.xml | tee runDemoPythia.log)
	add_dependencies(runDemoPythia DemoPythia)
	add_dependencies(runDemoPythia1000 DemoPythia)
	add_dependencies(runDemoPythia10000 DemoPythia)
	add_dependencies(runDemoPythia100000 DemoPythia)
	add_dependencies(runDemoPythia1000000 DemoPythia)
ENDIF(PYTHIA8_FOUND AND BOOST_FOUND)


add_subdirectory(vinhac)
add_subdirectory(utils)
