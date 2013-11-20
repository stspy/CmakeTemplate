FILE (GLOB src_files "*.*")

SET (ALL_SRCS ${src_files})

ADD_EXECUTABLE (VtkLineDisp ${ALL_SRCS})

TARGET_LINK_LIBRARIES (VtkLineDisp ${VTK_LIBRARIES})


INSTALL_TARGETS (/bin VtkLineDisp)