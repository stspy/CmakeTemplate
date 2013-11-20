#
# Source files
#
# Here is where you can add the name of your local Transform classes.
#

FILE (GLOB ALL_SOURCES "*.cpp" "*.c")
FILE (GLOB ALL_INCLUDES "*.hpp" "*.h")

SET (ALL_SRCS 
	 ${ALL_SOURCES}
	 ${ALL_INCLUDES}
)
 
# Create C++ library
ADD_LIBRARY (miiMinPath2D STATIC ${ALL_SRCS})


TARGET_LINK_LIBRARIES(miiMinPath2D  ${OpenCV_LIBS})

# The installation path for the link-library
INSTALL_TARGETS (/lib miiMinPath2D)

# The installation path for the include-files
INSTALL (FILES ${ALL_INCLUDES} DESTINATION include)
