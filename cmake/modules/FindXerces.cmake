

MESSAGE(STATUS "Looking for Xerces...")

FIND_LIBRARY(XERCES_CPP_LIB
NAMES
 xerces-c
PATHS
 /usr/lib
 /usr/local/lib
 ${VINHAC_BINARY_DIR}/lib/xerces-c-3.0.1/build/lib
)


FIND_PATH(XERCES_CPP_INCLUDE
 xercesc/util/PlatformUtils.hpp
 /usr/include
 /usr/local/include
 ${VINHAC_BINARY_DIR}/lib/xerces-c-3.0.1/build/include
)
 MESSAGE("Xerces from sources: ${XERCES_CPP_LIB} ${XERCES_CPP_INCLUDE}")
IF (NOT EXISTS ${XERCES_CPP_LIB} OR NOT EXISTS ${XERCES_CPP_INCLUDE})
 MESSAGE("Building Xerces from sources: ${XERCES_CPP_LIB} ${XERCES_CPP_INCLUDE}")
  EXEC_PROGRAM("cd lib && ./build_xerces.sh")
  find_library(XERCES xerces-c ${VINHAC_BINARY_DIR}/lib/xerces-c-3.0.1/build/lib)
  link_directories(${VINHAC_BINARY_DIR}/lib/xerces-c-3.0.1/build/lib)
  
ELSE (NOT EXISTS ${XERCES_CPP_LIB} OR NOT EXISTS ${XERCES_CPP_INCLUDE})
 SET(XERCES ${XERCES_CPP_LIB})
ENDIF (NOT EXISTS ${XERCES_CPP_LIB} OR NOT EXISTS ${XERCES_CPP_INCLUDE})

include_directories(${XERCES_CPP_INCLUDE})

