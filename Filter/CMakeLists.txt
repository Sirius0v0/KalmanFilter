file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h include/*.hpp)
add_library(Filters STATIC ${srcs})
target_include_directories(Filters PUBLIC include)
target_link_libraries(Filters PUBLIC Matrix)
