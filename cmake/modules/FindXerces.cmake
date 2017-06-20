

MESSAGE(STATUS "Looking for Xerces...")

FIND_LIBRARY(XERCES_CPP_LIB
NAMES
 xerces-c
PATHS
 /usr/lib
 /usr/local/lib
 ${XERCES_PREFIX}/lib
)


FIND_PATH(XERCES_CPP_INCLUDE
 xercesc/util/PlatformUtils.hpp
 /usr/include
 /usr/local/include
 ${XERCES_PREFIX}/include
)


IF ((NOT EXISTS ${XERCES_CPP_LIB}) OR (NOT EXISTS ${XERCES_CPP_INCLUDE}))
	MESSAGE( FATAL_ERROR "!!!!! XERCES-C library not found.        !!!!!\n"
			     "!!!!! Please go to the lib directory     !!!!!\n"
			     "!!!!! and install it, then edit          !!!!!\n"
                             "!!!!! config/build.properties file       !!!!!")
ELSE ((NOT EXISTS ${XERCES_CPP_LIB}) OR (NOT EXISTS ${XERCES_CPP_INCLUDE}))
 SET(XERCES ${XERCES_CPP_LIB})
ENDIF ((NOT EXISTS ${XERCES_CPP_LIB}) OR (NOT EXISTS ${XERCES_CPP_INCLUDE}))

MESSAGE(STATUS "Looking for Xerces... - found " ${XERCES_CPP_LIB} )
MESSAGE(STATUS "Looking for Xerces... - found " ${XERCES_CPP_INCLUDE} )
include_directories(${XERCES_CPP_INCLUDE})

