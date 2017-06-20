#Author: Kamil Sobol
#This reads Java-like properties build file

MESSAGE(STATUS "Reading config/build.properties file ...")

#read a file
FILE(READ "config/build.properties" BUILD_PROPERTIES)

#construct list of lines
STRING(REGEX REPLACE ";" "\\\\;" BUILD_PROPERTIES "${BUILD_PROPERTIES}")
STRING(REGEX REPLACE "\n" ";" BUILD_PROPERTIES "${BUILD_PROPERTIES}")

#iterate through the uncommented lines
FOREACH (PROPERTY ${BUILD_PROPERTIES})
	IF(NOT ${PROPERTY} MATCHES "^#.*$")
		STRING(REGEX REPLACE "=" ";" KEY_VALUE "${PROPERTY}")
			LIST(GET KEY_VALUE 0 KEY)
			LIST(GET KEY_VALUE 1 VALUE)
			SET(${KEY} ${VALUE})
	ENDIF(NOT ${PROPERTY} MATCHES "^#.*$")
ENDFOREACH(PROPERTY)

