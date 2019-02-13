

include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")

file(GLOB_RECURSE QUAD_SOURCES "${PROJECT_SOURCE_DIR}/src/QuadFormulas/*.C")
 
 
add_library(quad STATIC ${QUAD_SOURCES})

