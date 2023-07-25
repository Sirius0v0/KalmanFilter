file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h include/*.hpp)
add_library(KFExceptions STATIC ${srcs})
target_include_directories(KFExceptions PUBLIC include)
