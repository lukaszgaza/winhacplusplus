
#add_executable (testcmake testcmake.cpp)
add_executable (Demo Demo.cpp)

target_link_libraries(Demo ${LIBS} winhac winhacUtils winhacUtilsSanc)


add_subdirectory(vinhac)
add_subdirectory(utils)
