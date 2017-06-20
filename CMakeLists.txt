


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(VINHAC)


enable_language(Fortran)
set(CMAKE_CXX_FLAGS "-Wall -O2")
set(CMAKE_Fortran_FLAGS "-O2 ${CMAKE_Fortran_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# reading manual build properties
include(ReadBuildProperties)

# finding and configuring modules
include(FindROOT)
include(FindXerces)
include(UseDoxygen)
include(FindLHAPDF)
include(FindPythia8)
include(FindHepMC)
include(FindBoost)


IF(HEPMC_FOUND)
	add_definitions(-DHEPMC_ENABLED)
	include_directories(${HEPMC_INCLUDE})
	set(LIBS ${LIBS} ${HEPMC_LIB})
ENDIF(HEPMC_FOUND)


IF(ROOT_FOUND)
	add_definitions(-DROOT_HIST)
	include_directories(${ROOT_INCLUDE_DIR})
ENDIF(ROOT_FOUND)


SET(EXECUTABLE_OUTPUT_PATH ${VINHAC_BINARY_DIR}/build/exec)
SET(LIBRARY_OUTPUT_PATH ${VINHAC_BINARY_DIR}/build/lib)
set(CMAKE_BINARY_DIR ${CMAKE_BINARY_DIR}/build/bin)
set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(PROJECT_BUILD_DIR ${PROJECT_BINARY_DIR}/build)
set(PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR}/build/bin)
set(PROJECT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

set(LIBS ${LIBS} ${XERCES_CPP_LIB})


add_subdirectory(src)

INSTALL(
DIRECTORY ${PROJECT_SOURCE_DIR}/ 
DESTINATION ${CMAKE_INSTALL_PREFIX}/include
          FILES_MATCHING PATTERN "*.h" PATTERN "*svn*" EXCLUDE
          PATTERN "*CMakeFiles*" EXCLUDE
)



# CUSTOM TARGETS

add_custom_target(run1 ${EXECUTABLE_OUTPUT_PATH}/Demo input output 1 input/UserFile.xml | tee run.log)
add_dependencies(run1 Demo)

add_custom_target(run10 ${EXECUTABLE_OUTPUT_PATH}/Demo input output 10 input/UserFile.xml | tee run.log)
add_dependencies(run10 Demo)

add_custom_target(run1000 ${EXECUTABLE_OUTPUT_PATH}/Demo input output 1000 input/UserFile.xml  | tee run.log)
add_dependencies(run1000 Demo)

add_custom_target(run100000 ${EXECUTABLE_OUTPUT_PATH}/Demo input output 100000 input/UserFile.xml | tee run.log)
add_dependencies(run100000 Demo)

add_custom_target(run1000000 ${EXECUTABLE_OUTPUT_PATH}/Demo input output 1000000 input/UserFile.xml | tee run.log)
add_dependencies(run1000000 Demo)

add_custom_target(restore_input_from_backup cp -rf input.backup/* input )

add_custom_target(test_zero_mass_quarks cp -rf test/_zero_mass_quarks/* input )

add_custom_target(test_W+-_born_e cp -rf test/born/W+-/e/* input )
add_custom_target(test_W+_born_e cp -rf test/born/W+/e/* input )
add_custom_target(test_W-_born_e cp -rf test/born/W-/e/* input )
add_custom_target(test_W+-_born_e_mu cp -rf test/born/W+-/e-mu/* input )
add_custom_target(test_W+_born_e_mu cp -rf test/born/W+/e-mu/* input )
add_custom_target(test_W-_born_e_mu cp -rf test/born/W-/e-mu/* input )
add_custom_target(test_W+-_born_mu cp -rf test/born/W+-/mu/* input )
add_custom_target(test_W+_born_mu cp -rf test/born/W+/mu/* input )
add_custom_target(test_W-_born_mu cp -rf test/born/W-/mu/* input )

add_custom_target(test_W+-_born_ps_e cp -rf test/born-ps/W+-/e/* input )
add_custom_target(test_W+_born_ps_e cp -rf test/born-ps/W+/e/* input )
add_custom_target(test_W-_born_ps_e cp -rf test/born-ps/W-/e/* input )
add_custom_target(test_W+-_born_ps_e_mu cp -rf test/born-ps/W+-/e-mu/* input )
add_custom_target(test_W+_born_ps_e_mu cp -rf test/born-ps/W+/e-mu/* input )
add_custom_target(test_W-_born_ps_e_mu cp -rf test/born-ps/W-/e-mu/* input )
add_custom_target(test_W+-_born_ps_mu cp -rf test/born-ps/W+-/mu/* input )
add_custom_target(test_W+_born_ps_mu cp -rf test/born-ps/W+/mu/* input )
add_custom_target(test_W-_born_ps_mu cp -rf test/born-ps/W-/mu/* input )

add_custom_target(test_W+-_born_ps_had_e cp -rf test/born-ps-had/W+-/e/* input )
add_custom_target(test_W+_born_ps_had_e cp -rf test/born-ps-had/W+/e/* input )
add_custom_target(test_W-_born_ps_had_e cp -rf test/born-ps-had/W-/e/* input )
add_custom_target(test_W+-_born_ps_had_e_mu cp -rf test/born-ps-had/W+-/e-mu/* input )
add_custom_target(test_W+_born_ps_had_e_mu cp -rf test/born-ps-had/W+/e-mu/* input )
add_custom_target(test_W-_born_ps_had_e_mu cp -rf test/born-ps-had/W-/e-mu/* input )
add_custom_target(test_W+-_born_ps_had_mu cp -rf test/born-ps-had/W+-/mu/* input )
add_custom_target(test_W+_born_ps_had_mu cp -rf test/born-ps-had/W+/mu/* input )
add_custom_target(test_W-_born_ps_had_mu cp -rf test/born-ps-had/W-/mu/* input )

add_custom_target(test_W+-_rad_e cp -rf test/rad/W+-/e/* input )
add_custom_target(test_W+_rad_e cp -rf test/rad/W+/e/* input )
add_custom_target(test_W-_rad_e cp -rf test/rad/W-/e/* input )
add_custom_target(test_W+-_rad_e_mu cp -rf test/rad/W+-/e-mu/* input )
add_custom_target(test_W+_rad_e_mu cp -rf test/rad/W+/e-mu/* input )
add_custom_target(test_W-_rad_e_mu cp -rf test/rad/W-/e-mu/* input )
add_custom_target(test_W+-_rad_mu cp -rf test/rad/W+-/mu/* input )
add_custom_target(test_W+_rad_mu cp -rf test/rad/W+/mu/* input )
add_custom_target(test_W-_rad_mu cp -rf test/rad/W-/mu/* input )

add_custom_target(test_W+-_rad_ps_e cp -rf test/rad-ps/W+-/e/* input )
add_custom_target(test_W+_rad_ps_e cp -rf test/rad-ps/W+/e/* input )
add_custom_target(test_W-_rad_ps_e cp -rf test/rad-ps/W-/e/* input )
add_custom_target(test_W+-_rad_ps_e_mu cp -rf test/rad-ps/W+-/e-mu/* input )
add_custom_target(test_W+_rad_ps_e_mu cp -rf test/rad-ps/W+/e-mu/* input )
add_custom_target(test_W-_rad_ps_e_mu cp -rf test/rad-ps/W-/e-mu/* input )
add_custom_target(test_W+-_rad_ps_mu cp -rf test/rad-ps/W+-/mu/* input )
add_custom_target(test_W+_rad_ps_mu cp -rf test/rad-ps/W+/mu/* input )
add_custom_target(test_W-_rad_ps_mu cp -rf test/rad-ps/W-/mu/* input )

add_custom_target(test_W+-_rad_ps_had_e cp -rf test/rad-ps-had/W+-/e/* input )
add_custom_target(test_W+_rad_ps_had_e cp -rf test/rad-ps-had/W+/e/* input )
add_custom_target(test_W-_rad_ps_had_e cp -rf test/rad-ps-had/W-/e/* input )
add_custom_target(test_W+-_rad_ps_had_e_mu cp -rf test/rad-ps-had/W+-/e-mu/* input )
add_custom_target(test_W+_rad_ps_had_e_mu cp -rf test/rad-ps-had/W+/e-mu/* input )
add_custom_target(test_W-_rad_ps_had_e_mu cp -rf test/rad-ps-had/W-/e-mu/* input )
add_custom_target(test_W+-_rad_ps_had_mu cp -rf test/rad-ps-had/W+-/mu/* input )
add_custom_target(test_W+_rad_ps_had_mu cp -rf test/rad-ps-had/W+/mu/* input )
add_custom_target(test_W-_rad_ps_had_mu cp -rf test/rad-ps-had/W-/mu/* input )

################



SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "VINHAC")
SET(CPACK_PACKAGE_VENDOR "Kamil Sobol")
SET(CPACK_GENERATOR "TGZ")
SET(CPACK_PACKAGE_EXECUTABLES "Demo" "Combine")
INCLUDE(CPack)


#Some variables printout.
#include(DebugInfo)


