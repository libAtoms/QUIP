PROJECT(FoX Fortran)

SET(CMAKE_VERBOSE_MAKEFILE FALSE)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#======================#
# Enable testing
#======================#

ENABLE_TESTING()
INCLUDE(CTest)

#======================#
# Build paths
#======================#

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})


#======================#
# Enable/disable libraries
#======================#

OPTION(${PROJECT_NAME}_ENABLE_WXML "Enable/disable general Fortran XML output library" ON)
OPTION(${PROJECT_NAME}_ENABLE_WCML "Enable/disable Chemical Markup Language output library" ON)
OPTION(${PROJECT_NAME}_ENABLE_WKML "Enable/disable Keyhole Markup Language output library" ON)
OPTION(${PROJECT_NAME}_ENABLE_DOM "Enable/disable Document Object Model input library" ON)
OPTION(${PROJECT_NAME}_ENABLE_EXAMPLES "Build examples" ON)

#======================#
# set another compiler #
#======================#
# to set another compiler than the default, 
# do the usual from within the build directory:
# $ FC=g95 
# $ rm -rf * && cmake ../
# $ make 

#=====================#
# verbose make output #
#=====================#
# to see out all commands used during building,
# do as follows:
# $ rm -rf * && cmake ../
# $ make VERBOSE=1

#======================#
# speed up compilation #
#======================#
# to speed up compilation by using multiple cores/processors,
# do the following:
# $ rm -rf * && cmake ../
# $ make -j

#==============================#
# get the machine architecture #
#==============================#

SET(OS ${CMAKE_SYSTEM})
SET(MACHINE_ARCH ${CMAKE_SYSTEM_PROCESSOR})

MESSAGE(STATUS "determining operating system and architecture:")
MESSAGE("   -> your operating system is : ${OS}")
MESSAGE("   -> your architecture is     : ${MACHINE_ARCH}")


#==========================#
# check if m4 is available #
#==========================#

MESSAGE(STATUS "Searching for m4 scripting language")
FIND_PROGRAM(M4_PROGRAM m4)
IF(${M4_PROGRAM} MATCHES "NOTFOUND")
  MESSAGE("   -> WARNING : could not find m4, cannot generate fortran files from m4 sources")
ELSE(${M4_PROGRAM} MATCHES "NOTFOUND")
  MESSAGE("   -> ${M4_PROGRAM}")
ENDIF(${M4_PROGRAM} MATCHES "NOTFOUND")

#============================#
# set some other definitions #
#============================#

# file name extensions
SET(m4_ext .m4)
SET(f90pp_ext .F90)
SET(f90_ext .f90)

###===================================###
###===================================###
### determine some preprocessor flags ###
###===================================###
###===================================###

#===========================================================#
# determine end-of-line character based on operating system #
#===========================================================#

INCLUDE(${${PROJECT_NAME}_SOURCE_DIR}/cmake/DetermineEOL.cmake)

#============================================#
# determine if and how FLUSH intrinsic works #
#============================================#

INCLUDE(${${PROJECT_NAME}_SOURCE_DIR}/cmake/CheckFlushIntrinsic.cmake)

#=====================================================#
# check for 'associated in restricted expression' bug #
#=====================================================#

INCLUDE(${${PROJECT_NAME}_SOURCE_DIR}/cmake/CheckAssociatedBug.cmake)

#============================================#
# determine if and how ABORT intrinsic works #
#============================================#

INCLUDE(${${PROJECT_NAME}_SOURCE_DIR}/cmake/CheckAbortIntrinsic.cmake)

SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FPPFLAGS}")

# Base library directories
ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/fsys)
ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/utils)
ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/common)

# wxml directory
IF(${PROJECT_NAME}_ENABLE_WXML OR ${PROJECT_NAME}_ENABLE_WCML OR ${PROJECT_NAME}_ENABLE_WKML OR ${PROJECT_NAME}_ENABLE_DOM OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    SET(${PROJECT_NAME}_ENABLE_WXML ON)
    ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/wxml)
ELSE()
    GET_PROPERTY(THE_INCLUDE_DIRS GLOBAL PROPERTY INCLUDE_DIRECTORIES)
    STRING(REPLACE ${${PROJECT_NAME}_SOURCE_DIR}/wxml "" NEW_INCLUDE_DIRS "${THE_INCLUDE_DIRS}")
    SET_PROPERTY(GLOBAL PROPERTY INCLUDE_DIRECTORIES ${NEW_INCLUDE_DIRS})
ENDIF()

# wkml directory
IF(${PROJECT_NAME}_ENABLE_WKML OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    SET(${PROJECT_NAME}_ENABLE_WKML ON)
    ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/wkml)
ELSE()
    GET_PROPERTY(THE_INCLUDE_DIRS GLOBAL PROPERTY INCLUDE_DIRECTORIES)
    STRING(REPLACE ${${PROJECT_NAME}_SOURCE_DIR}/wkml "" NEW_INCLUDE_DIRS "${THE_INCLUDE_DIRS}")
    SET_PROPERTY(GLOBAL PROPERTY INCLUDE_DIRECTORIES ${NEW_INCLUDE_DIRS})
ENDIF()

# wcml directory
IF(${PROJECT_NAME}_ENABLE_WCML OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    SET(${PROJECT_NAME}_ENABLE_WCML ON)
    ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/wcml)
ELSE()
    GET_PROPERTY(THE_INCLUDE_DIRS GLOBAL PROPERTY INCLUDE_DIRECTORIES)
    STRING(REPLACE ${${PROJECT_NAME}_SOURCE_DIR}/wcml "" NEW_INCLUDE_DIRS "${THE_INCLUDE_DIRS}")
    SET_PROPERTY(GLOBAL PROPERTY INCLUDE_DIRECTORIES ${NEW_INCLUDE_DIRS})
ENDIF()

# dom and sax directories
IF(${PROJECT_NAME}_ENABLE_DOM OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    SET(${PROJECT_NAME}_ENABLE_DOM ON)
    ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/sax)
    ADD_SUBDIRECTORY(${${PROJECT_NAME}_SOURCE_DIR}/dom)
ELSE()
    GET_PROPERTY(THE_INCLUDE_DIRS GLOBAL PROPERTY INCLUDE_DIRECTORIES)
    STRING(REPLACE ${${PROJECT_NAME}_SOURCE_DIR}/sax "" TEMP_DIRS "${THE_INCLUDE_DIRS}")
    STRING(REPLACE ${${PROJECT_NAME}_SOURCE_DIR}/dom "" NEW_INCLUDE_DIRS "${TEMP_DIRS}")
    SET_PROPERTY(GLOBAL PROPERTY INCLUDE_DIRECTORIES ${NEW_INCLUDE_DIRS})
ENDIF()

# examples directory
IF(${PROJECT_NAME}_ENABLE_EXAMPLES OR ${PROJECT_NAME}_ENABLE_EXAMPLES)
    SET(${PROJECT_NAME}_ENABLE_WCML ON)
    ADD_SUBDIRECTORY(examples)
ELSE()
    GET_PROPERTY(THE_INCLUDE_DIRS GLOBAL PROPERTY INCLUDE_DIRECTORIES)
    STRING(REPLACE "${${PROJECT_NAME}_SOURCE_DIR}/examples" "" NEW_INCLUDE_DIRS "${THE_INCLUDE_DIRS}")
    SET_PROPERTY(GLOBAL PROPERTY INCLUDE_DIRECTORIES ${NEW_INCLUDE_DIRS})

ENDIF()


