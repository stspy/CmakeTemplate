cmake_minimum_required(VERSION 2.6)
#
# Set the project name.
#

PROJECT (CmakeTemplate)

SET (EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin CACHE INTERNAL
  "Where to put the executables")

SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE INTERNAL
  "Where to put the libraries")

SET (DEBUG 1)

#IF (DEBUG)
#	SET (ENV{VTK_DIR} D:/ProgramFiles/VTK6.0/VTK_Win64_Debug)
#ELSE ()
#	SET (ENV{VTK_DIR} D:/ProgramFiles/VTK6.0/VTK_Win64_Release)
#ENDIF ()


SUBDIRS ( 
  Apps 
  ExtLibs
  Tempclass
  )


#FIND_PACKAGE(VTK REQUIRED)
#INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR}/TempClass/wcdMinPath)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})




 
  
