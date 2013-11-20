FILE(GLOB src_files "*.*")


SET (ALL_SRCS ${src_files})


ADD_EXECUTABLE (ImageReader ${ALL_SRCS})


TARGET_LINK_LIBRARIES (ImageReader ${OpenCV_LIBS} miiMinPath2D)



INSTALL_TARGETS (/bin ImageReader)