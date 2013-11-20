FILE(GLOB src_files "*.*")


SET (ALL_SRCS ${src_files})


ADD_EXECUTABLE (VideoCap ${ALL_SRCS})


TARGET_LINK_LIBRARIES(VideoCap ${OpenCV_LIBS} miiMinPath2D)



INSTALL_TARGETS(/bin VideoCap)