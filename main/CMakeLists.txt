FIND_PACKAGE( ROOT REQUIRED)





AUX_SOURCE_DIRECTORY( src library_sources )
SET(LandauGauss_name LandauGaussDict.cxx)

file(GLOB Demo_HEADERS *.h)
source_group("Headers" FILES ${Demo_HEADERS})



INCLUDE_DIRECTORIES(  ${ROOT_INCLUDE_DIR} )

# fill list of header files 
set(LandauGauss_HEADERS 
 C:/projects/git/LandauGauss/main/include/landgausFit.h
  )

SET(MY_INCLUDE_DIRECTORIES ./main/inc ./inc ../../main/inc ${ROOT_INCLUDE_DIR} )

# set everything needed for the root dictonary and create the
# dictionary
set(LandauGauss_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/${LandauGauss_name}) 

ROOT_GENERATE_DICTIONARY_OLD("${LandauGauss_HEADERS}" "${LandauGauss_DICTIONARY}" "${MY_INCLUDE_DIRECTORIES}")
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${library_sources}  ${LandauGauss_name} )



SOURCE_GROUP(Headers  inc/*.h)

target_link_libraries( ${PROJECT_NAME} ${EUDAQ_THREADS_LIB} ${ROOT_LIBRARIES} ${ADDITIONAL_LIBRARIES})


INSTALL(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
