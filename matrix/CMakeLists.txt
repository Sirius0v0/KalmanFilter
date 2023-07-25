file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h include/*.hpp)
add_library(matrix STATIC ${srcs})
target_include_directories(matrix PUBLIC include)
