

include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")

file(GLOB_RECURSE REF_SOURCES "${PROJECT_SOURCE_DIR}/src/Refinement/*.C")
 
 
add_library(ref STATIC ${REF_SOURCES})

