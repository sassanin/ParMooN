

include_directories("${CMAKE_SOURCE_DIR}/include/General")

file(GLOB_RECURSE GEN_SOURCES "${PROJECT_SOURCE_DIR}/src/General/*.C")
 
 
add_library(gen STATIC ${GEN_SOURCES})

