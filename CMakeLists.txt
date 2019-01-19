#############################################################################
#
# CMAKE file for building tkLayout software
#
#  Author: Z. Drasal, CERN
#
#############################################################################

#----------------------------------------------------------------------------
# CMAKE compatibility issues
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)

#----------------------------------------------------------------------------
# Project name and version
PROJECT(tklayout)

#----------------------------------------------------------------------------
# Project options
OPTION(INSTALL_DOC "Set to ON to build/install Documentation" OFF )


SET(CXX_STANDARD 17)

# Output directories
SET( EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}/bin" CACHE PATH
    "EXECUTABLE_OUTPUT_PATH" FORCE )
SET( LIBRARY_OUTPUT_PATH "${PROJECT_BINARY_DIR}/lib" CACHE PATH
    "LIBRARY_OUTPUT_PATH" FORCE )
MARK_AS_ADVANCED( EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH )


#----------------------------------------------------------------------------
# Project dependencies

FIND_PACKAGE( Boost COMPONENTS program_options filesystem system REQUIRED )

find_package(ROOT COMPONENTS Geom RIO HistPainter)
include_directories(${ROOT_INCLUDE_DIRS})
link_libraries(${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Include, source, dirs 
#
# Compiler flags
#SET ( CMAKE_CXX_COMPILER "g++" )
ADD_DEFINITIONS( "-Wl,--copy-dt-needed-entries" )
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17  -g -pedantic -Wall -Wno-reorder")

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/source/Geometry/include 
                     ${PROJECT_SOURCE_DIR}/source/Analysis/include
                     ${PROJECT_SOURCE_DIR}/source/Extractors/include
                     ${PROJECT_SOURCE_DIR}/source/TinyXML/include
                     ${PROJECT_SOURCE_DIR}/source/Tools/include
                     ${PROJECT_SOURCE_DIR}/source/Tracking/include
                     ${PROJECT_SOURCE_DIR}/source/Web/include
                     ${Boost_INCLUDE_DIR}
                     ${ROOT_INCLUDE_DIRS} )
FILE( GLOB exe_sources ${PROJECT_SOURCE_DIR}/source/*.cc )
FILE( GLOB obj_sources ${PROJECT_SOURCE_DIR}/source/Geometry/src/*.cc
                       ${PROJECT_SOURCE_DIR}/source/Analysis/src/*.cc 
                       ${PROJECT_SOURCE_DIR}/source/Extractors/src/*.cc
                       ${PROJECT_SOURCE_DIR}/source/Tools/src/*.cc 
                       ${PROJECT_SOURCE_DIR}/source/Tracking/src/*.cc
                       ${PROJECT_SOURCE_DIR}/source/Web/src/*.cc )
FILE( GLOB obj_headers ${PROJECT_SOURCE_DIR}/source/Geometry/include/*.h
                       ${PROJECT_SOURCE_DIR}/source/Analysis/include/*.h
                       ${PROJECT_SOURCE_DIR}/source/Extractors/include/*.h 
                       ${PROJECT_SOURCE_DIR}/source/Tools/include/*.h  
                       ${PROJECT_SOURCE_DIR}/source/Tracking/include/*.h 
                       ${PROJECT_SOURCE_DIR}/source/Web/include/*.h ) 
FILE( GLOB tinyxml_sources ${PROJECT_SOURCE_DIR}/source/TinyXML/src/*.cc )

# Get source files with main function defined
SET ( sources "" )
FOREACH( file ${exe_sources} )
 IF ( ${file} MATCHES "tklayout.cc" ) 
  SET( source_tklayout ${file} )
 ENDIF()
 IF ( ${file} MATCHES "setup.cc" ) 
  SET( source_setup ${file} )
 ENDIF()
 IF ( ${file} MATCHES "delphize.cc" )
  SET( source_delphize ${file} )
 ENDIF()
ENDFOREACH()

# Get files on which setup or delphize depend
FOREACH( file ${obj_sources} )
 IF( ${file} MATCHES "MainConfigHandler.cc" )
  SET( source_mainhandler ${file} )
 ENDIF()
 IF( ${file} MATCHES "GraphVizCreator.cc" )
  SET( source_vizgraphcreator ${file} )
 ENDIF()
 IF( ${file} MATCHES "MessageLogger.cc" )
  SET( source_messagelogger ${file} )
 ENDIF()
 IF( ${file} MATCHES "string_functions.cc" )
  SET( source_stringfunctions ${file} )
 ENDIF()
  IF( ${file} MATCHES "StringConverter.cc" )
  SET( source_stringconverter ${file} )
 ENDIF()
ENDFOREACH()

#----------------------------------------------------------------------------
# Set revision
# a custom target that is always built
ADD_CUSTOM_TARGET(revisiontag ALL)

# creates svnversion.h using cmake script
ADD_CUSTOM_COMMAND(TARGET revisiontag COMMAND ${CMAKE_COMMAND}
   -DSOURCE_DIR=${CMAKE_SOURCE_DIR} 
   -P ${CMAKE_SOURCE_DIR}/GetGitVersion.cmake)

#--------------------------------------------------------
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)            # Set automatically rpath for dynamic linking of external libraries
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") # Set automatically rpath for dynamic linking of internally created & installed libraries

ADD_EXECUTABLE(tklayout ${source_tklayout} ${obj_sources} ${obj_headers} )
ADD_EXECUTABLE(setup.bin ${source_setup} ${source_vizgraphcreator} ${source_messagelogger} ${source_mainhandler} ${source_stringfunctions} ${source_stringconverter} ${obj_headers} )
ADD_EXECUTABLE(delphize ${source_delphize} )

ADD_LIBRARY(tinyXML SHARED ${tinyxml_sources})

# explicitly say that the executable depends on custom target
ADD_DEPENDENCIES(tklayout revisiontag)

TARGET_LINK_LIBRARIES(tklayout ${Boost_LIBRARIES} ${ROOT_LIBRARIES} tinyXML )
TARGET_LINK_LIBRARIES(setup.bin ${Boost_LIBRARIES} )
TARGET_LINK_LIBRARIES(delphize ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
INSTALL(TARGETS tklayout  RUNTIME DESTINATION bin)
INSTALL(TARGETS setup.bin RUNTIME DESTINATION bin)
INSTALL(TARGETS delphize  RUNTIME DESTINATION bin)
INSTALL(TARGETS tinyXML   ARCHIVE DESTINATION lib
                          LIBRARY DESTINATION lib
                          RUNTIME DESTINATION bin)
INSTALL(DIRECTORY style DESTINATION share/tkLayout)
INSTALL(DIRECTORY xml DESTINATION share/tkLayout)
INSTALL(DIRECTORY geometries DESTINATION share/tkLayout)
INSTALL(DIRECTORY config DESTINATION share/tkLayout)

IF(CMAKE_HOST_UNIX)
    
    # Set the default version of GitRevision.cc file
    INSTALL( CODE "MESSAGE(STATUS \"Setting GitRevision.cc to default\" )
                   EXECUTE_PROCESS(COMMAND bash -c \"if [ -e GitRevision.orig.cc ]; then rm -f ${PROJECT_SOURCE_DIR}/source/Tools/src/GitRevision.cc; mv GitRevision.orig.cc ${PROJECT_SOURCE_DIR}/source/Tools/src/GitRevision.cc; fi\")" )
    

    # Configure tkLayout
    #INSTALL( CODE "MESSAGE(STATUS \"Configuring tkLayout software...\" )
    #               MESSAGE(STATUS \" Target directories setup:\" )
    #		   EXECUTE_PROCESS(COMMAND bash -c \"if [ -e ${PROJECT_BINARY_DIR}/bin/setup.bin ]; then ${PROJECT_BINARY_DIR}/bin/setup.bin --dirNames; fi\")" ) 
    #INSTALL( CODE "MESSAGE(STATUS \" Installation status:\" )
    #               EXECUTE_PROCESS(COMMAND bash -c \"if ! ${PROJECT_SOURCE_DIR}/bin/setup.bin --checkDir; then echo \\\"Problem during installation\\\"; else echo \\\"Installation successful\\\"; fi\")" )

ENDIF()

#----------------------------------------------------------------------------
# Documentation
#
# Find Doxygen
IF( INSTALL_DOC )
FIND_PACKAGE( Doxygen )
IF( DOXYGEN_FOUND )

    # Add doc target
    ADD_CUSTOM_TARGET( doc ) 
    ADD_CUSTOM_COMMAND( TARGET doc
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${PROJECT_SOURCE_DIR}/README.md" "${PROJECT_SOURCE_DIR}/doc/README.md"    
        COMMAND "${DOXYGEN_EXECUTABLE}" ARGS "Doxyfile"
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
        COMMENT "Building API Documentation in ${PROJECT_SOURCE_DIR}/doc... "
        VERBATIM )

ELSE()
    MESSAGE( STATUS "Doxygen not found in your system!!" )
    IF( INSTALL_DOC )
        MESSAGE( STATUS "INSTALL_DOC forced to OFF" )
        SET( INSTALL_DOC OFF )
    ENDIF()
ENDIF()

# Install documentation
    # make sure doxygen is executed (make doc) before make install
    INSTALL( CODE "EXEC_PROGRAM(${CMAKE_BUILD_TOOL} ${PROJECT_BINARY_DIR} ARGS doc)" )
    # install documentation
    INSTALL( DIRECTORY "${PROJECT_SOURCE_DIR}/doc"
            DESTINATION .
            PATTERN "*CVS*" EXCLUDE )
ENDIF()


#----------------------------------------------------------------------------
# Uninstall 
#
CONFIGURE_FILE(
    "${CMAKE_CURRENT_SOURCE_DIR}/CMakeUninstall.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
    
