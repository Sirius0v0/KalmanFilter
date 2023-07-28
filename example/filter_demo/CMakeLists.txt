file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h)
add_executable(filter_demo ${srcs})
target_include_directories(filter_demo PUBLIC include)

target_link_libraries(filter_demo PUBLIC Filters)