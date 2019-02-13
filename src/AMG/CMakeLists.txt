

include_directories("${CMAKE_SOURCE_DIR}/include/AMG")

file(GLOB_RECURSE AMG_SOURCES "${PROJECT_SOURCE_DIR}/src/AMG/*.c")

add_library(amg STATIC ${AMG_SOURCES} )

