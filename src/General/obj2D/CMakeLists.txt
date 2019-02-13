

include(${CMAKE_SOURCE_DIR}/config.cmake)

include_directories("${CMAKE_SOURCE_DIR}/include/General")

file(GLOB_RECURSE GEO_SOURCES "${PROJECT_SOURCE_DIR}/src/General/*.C")
 
 
add_library(gen STATIC ${GEO_SOURCES})

