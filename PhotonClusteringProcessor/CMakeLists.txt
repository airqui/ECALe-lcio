########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 3.12 FATAL_ERROR)
########################################################

# fix for lcio: all c++ sources need to be compiled with -Wno-long-long
SET( ENV{CXXFLAGS} "-Wno-long-long $ENV{CXXFLAGS}" )

# require proper c and c++
# -pedantic cannot be set globally because it cannot be removed for certain
# source files only: grep pedantic src/CMakeLists.txt
#SET( ENV{CXXFLAGS} "-Wall -ansi $ENV{CXXFLAGS}" ) # FIXME -pedantic
SET( ENV{CFLAGS} "-w -ansi -pedantic $ENV{CFLAGS}" ) # -w suppress all warning messages

# project name
PROJECT( TwoPhotonClustering )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )
#SET (CMAKE_CXX_STANDARD 11)


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )

#FIND_PACKAGE( Boost REQUIRED )
#include_directories( SYSTEM ${Boost_INCLUDE_DIRS} )

FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( ROOT REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( GSL REQUIRED )

FIND_PACKAGE ( FastJet )
IF ( NOT ${FastJet_VERSION} VERSION_LESS 3 )
  FIND_PACKAGE( FastJet 3 REQUIRED COMPONENTS siscone siscone_spherical fastjetplugins fastjetcontribfragile fastjettools)
ELSEIF ( NOT ${FastJet_VERSION} VERSION_LESS 2.4.2 )
  FIND_PACKAGE( FastJet 2 REQUIRED COMPONENTS SISConePlugin )
ENDIF()

INCLUDE_DIRECTORIES( SYSTEM ${FastJet_INCLUDE_DIRS} )
LINK_LIBRARIES( ${FastJet_LIBRARIES} ${FastJet_COMPONENT_LIBRARIES} )

FOREACH( pkg Marlin MarlinUtil LCIO ROOT GSL )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( BEFORE SYSTEM ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

#INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${Marlin_LIBRARIES} )
#LINK_LIBRARIES( ${ROOT_LIBRARIES} )
#LINK_LIBRARIES( ${LCIO_LIBRARIES} )
#ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()

FIND_PACKAGE( aidaTT )
ADD_DEFINITIONS( "-D AIDATT_USE_DD4HEP -DUSE_GBL -DUSE_LCIO" )
INCLUDE_DIRECTORIES( SYSTEM ${aidaTT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${aidaTT_LIBRARIES} )
ADD_DEFINITIONS( ${aidaTT_DEFINITIONS} )


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++17" )
#ADD_DEFINITIONS( "-Wall -ansi -pedantic -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs" )
#ADD_DEFINITIONS( "-Wall -ansi -pedantic")
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
