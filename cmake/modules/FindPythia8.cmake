#Author: Kamil Sobol
# This module tries to find the Pythia8 installation on your system.

MESSAGE(STATUS "Looking for Pythia8...")


# try to find Pythia8 in user defined path
FIND_LIBRARY(PYTHIA8_LIB
	NAMES
		pythia8
	PATHS
		${PYTHIA8_PREFIX}/lib
	NO_DEFAULT_PATH
)


# if not try to find Pythia8 in standard instalation paths
IF(${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND")
	FIND_LIBRARY(PYTHIA8_LIB
		NAMES
			pythia8
		PATHS
			/usr/lib
			/usr/local/lib
	)
	FIND_PATH(PYTHIA8_INCLUDE
	 	Pythia.h
	 PATHS
		 /usr/include
		 /usr/local/include
	)
ELSE(${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND")
	FIND_PATH(PYTHIA8_INCLUDE
	 	Pythia.h
	 PATHS
		 ${PYTHIA8_PREFIX}/include
	 NO_DEFAULT_PATH
	)
ENDIF(${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND")


# final printout.
IF((${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND") OR (${PYTHIA8_INCLUDE} MATCHES "PYTHIA8_INCLUDE-NOTFOUND"))
	SET(PYTHIA8_FOUND FALSE)
	MESSAGE( STATUS	 "\n"
			     "!!!!! Pythia8                            !!!!!\n"
			     "!!!!! shared library or includes         !!!!!\n"
			     "!!!!! not found.                         !!!!!\n"
			     "!!!!! DemoPythia will not be built       !!!!!\n"
			     "!!!!! If you have it installed           !!!!!\n"
			     "!!!!! in custom localisation please edit !!!!!\n"
                             "!!!!! config/build.properties file       !!!!!")
ELSE((${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND") OR (${PYTHIA8_INCLUDE} MATCHES "PYTHIA8_INCLUDE-NOTFOUND"))
	SET(PYTHIA8_FOUND TRUE)
	MESSAGE(STATUS "Looking for Pythia8... - found " ${PYTHIA8_LIB} )
	MESSAGE(STATUS "Looking for Pythia8... - found " ${PYTHIA8_INCLUDE} )
ENDIF((${PYTHIA8_LIB} MATCHES "PYTHIA8_LIB-NOTFOUND") OR (${PYTHIA8_INCLUDE} MATCHES "PYTHIA8_INCLUDE-NOTFOUND"))


