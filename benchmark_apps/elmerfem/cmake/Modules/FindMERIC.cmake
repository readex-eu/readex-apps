# CMake script for finding MERIC for Elmer
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

message("FindMERIC.cmake started...")

SET(_mericIfaceF90 "meric_mod.F90")
SET(_mericLibName "libmericmpi.so")

# If MKL_LIBRARIES libraries are already defined, do nothing
SET(MERIC_FOUND FALSE)

IF (NOT MERIC_ROOT)
  SET(MERIC_ROOT "$ENV{MERIC_ROOT}")
ENDIF()

#SET(_mericIncludePaths
#  "$ENV{MERIC_ROOT}/include"
#  "${MERIC_ROOT}/include"
#  INTERNAL
#)

SET(_mericLibPaths
  "$ENV{MERIC_ROOT}"
  "${MERIC_ROOT}"
  INTERNAL
)

#TODO different subdir
SET(_mericInterfaceSrcPaths
  "$ENV{MERIC_ROOT}/${_mericIfaceF90}"
  "${MERIC_ROOT}/${_mericIfaceF90}"
  INTERNAL
)

# Find Feti4i library
#FIND_LIBRARY(MERIC_LIBRARIES ${_mericLibName}${SHL_EXTENSION} HINTS ${_mericLibPaths})
#SET(MERIC_LIBRARIES ${PROJECT_SOURCE_DIR}/meric/meric_fortran_test/lib/${_mericLibName} CACHE FILE "")
SET(MERIC_LIBRARIES "${MERIC_ROOT}lib/${_mericLibName}" CACHE FILE "")

# Find the actual interface file
#FIND_FILE(MERIC_INTERFACE_SOURCE NAMES ${_mericIfaceF90} PATHS ${_mericInterfaceSrcPaths})
#SET(MERIC_INTERFACE_SOURCE ${PROJECT_SOURCE_DIR}/meric/${_mericIfaceF90} CACHE FILE "")
SET(MERIC_INTERFACE_SOURCE "${MERIC_ROOT}include/${_mericIfaceF90}" CACHE FILE "")

message(STATUS "MERIC_LIBRARIES=${MERIC_LIBRARIES}")
message(STATUS "MERIC_INTERFACE_SOURCE=${MERIC_INTERFACE_SOURCE}")


IF(MERIC_LIBRARIES AND MERIC_INTERFACE_SOURCE)
  SET(MERIC_FOUND TRUE)
ENDIF()

IF(MERIC_FOUND)
  IF (NOT MERIC_FIND_QUIETLY)
    MESSAGE(STATUS "A library with MERIC API found.")
  ENDIF()
ELSE()
  IF (MERIC_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${MERIC_FAILMSG})
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  MERIC_LIBRARIES
  MERIC_INTERFACE_SOURCE
)
