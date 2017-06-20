#Author: Kamil Sobol
# This module tries to find the Boost installation on your system.

MESSAGE(STATUS "Looking for Boost...")


# try to find Boost in user defined path
FIND_LIBRARY(BOOST_LIB
	NAMES
		boost_thread-mt
	PATHS
		${BOOST_PREFIX}/stage/lib
	NO_DEFAULT_PATH
)


# if not try to find Boost in standard instalation paths
IF(${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND")
	FIND_LIBRARY(BOOST_LIB
		NAMES
			boost_thread-mt
		PATHS
			/usr/lib
			/usr/lib64
			/usr/local/lib
	)
	FIND_PATH(BOOST_INCLUDE
	 	boost/version.hpp
	 PATHS
	 	/usr/include
	 	/usr/local/include
	)
ELSE(${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND")
	FIND_PATH(BOOST_INCLUDE
	 	boost/version.hpp
	 PATHS
		 ${BOOST_PREFIX}
	 NO_DEFAULT_PATH
	)
ENDIF(${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND")


# final printout.
IF((${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND") OR (${BOOST_INCLUDE} MATCHES "BOOST_INCLUDE-NOTFOUND"))
	SET(BOOST_FOUND FALSE)
	MESSAGE( STATUS	 "\n"
			     "!!!!! Boost                            !!!!!\n"
			     "!!!!! shared library or includes         !!!!!\n"
			     "!!!!! not found.                         !!!!!\n"
			     "!!!!! DemoPythia will not be built       !!!!!\n"
			     "!!!!! If you have it installed           !!!!!\n"
			     "!!!!! in custom localisation please edit !!!!!\n"
                             "!!!!! config/build.properties file       !!!!!")
ELSE((${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND") OR (${BOOST_INCLUDE} MATCHES "BOOST_INCLUDE-NOTFOUND"))
	SET(BOOST_FOUND TRUE)
	MESSAGE(STATUS "Looking for Boost... - found " ${BOOST_LIB} )
	MESSAGE(STATUS "Looking for Boost... - found " ${BOOST_INCLUDE} )
ENDIF((${BOOST_LIB} MATCHES "BOOST_LIB-NOTFOUND") OR (${BOOST_INCLUDE} MATCHES "BOOST_INCLUDE-NOTFOUND"))


