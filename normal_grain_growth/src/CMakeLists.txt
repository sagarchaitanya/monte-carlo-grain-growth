project(grain_growth)

set(CMAKE_BUILD_TYPE Release)

set(HEADERS
	${HEADER_DIR}/main.h
)

set(SOURCES
#	*..cpp
)

add_executable(a.out main.cpp ${SOURCES} ${HEADERS})
target_link_libraries(a.out
	${MPI_LIBRARIES}
	${VTK_LIBRARIES}
)
