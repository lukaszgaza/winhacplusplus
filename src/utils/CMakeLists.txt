
aux_source_directory(. UTILS_SOURCES)
aux_source_directory(./FOAM FOAM_SOURCES)
aux_source_directory(./CLHEP/Random/src CLHEP_SOURCES)
aux_source_directory(./CLHEP/Vector/src CLHEP_SOURCES)
set(SANC_SOURCES 
./SANC/ewcdec/B0D_gen.f  ./SANC/ewcdec/BFF_gen.f  
./SANC/ewcdec/C0_gen.f  ./SANC/ewcdec/delta_r.f   
./SANC/ewcdec/fermse.f  ./SANC/ewcdec/W_decay_main.f 
./SANC/ewcdec/B0F_gen.f  ./SANC/ewcdec/bse_bc.f   
./SANC/ewcdec/C0IR.f    ./SANC/ewcdec/dhardc.f    
./SANC/ewcdec/ineett.f  ./SANC/ewcdec/W_dec.f 
./SANC/ewcdec/B0P_gen.f  ./SANC/ewcdec/bse_fc.f   
./SANC/ewcdec/C0IRF.f   ./SANC/ewcdec/eett_aux.f  
./SANC/ewcdec/JintW.f   ./SANC/ewcdec/xspenz.f
 
)

add_library(winhacUtils SHARED ${UTILS_SOURCES} ${FOAM_SOURCES} ${CLHEP_SOURCES} ./SANC/wh_ewcdec.f )
add_library(winhacUtilsSanc SHARED ${SANC_SOURCES})

add_executable (Combine farm/Combine.cpp)


install(TARGETS winhacUtils
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)
            
install(TARGETS winhacUtilsSanc
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)


