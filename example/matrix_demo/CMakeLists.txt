file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)
add_executable(matrix_demo ${srcs})
target_include_directories(matrix_demo PUBLIC include)

target_link_libraries(matrix_demo PUBLIC matrix)