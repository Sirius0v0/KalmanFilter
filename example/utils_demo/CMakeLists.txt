file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)
add_executable(utils_demo ${srcs})
target_include_directories(utils_demo PUBLIC include)
target_link_libraries(utils_demo PUBLIC utils)
