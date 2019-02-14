# ========================================================================
#     This file is part of ParMooN.
# 
#     ParMooN (cmg.cds.iisc.ac.in/parmoon) is free software:
#     you can redistribute it and/or modify it under the terms of the 
#     GNU General Public License as published by the Free Software 
#     Foundation, either version 3 of the License, or (at your option)
#     any later version. 
#
#     ParMooN is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
# =======================================================================
# PURPOSE :  configuration file for the ParMooN Version 1.1
# AUTHOR  : Sashikumaar Ganesan, CDS, IISc Bangalore, India
# DATE    : 05 June 2015
# ========================================================================

# normally no need to alter any thing in this file
# please DO NOT modify unless if it is necessary
# beginners are advised NOT to alter anything the below contents

cmake_minimum_required (VERSION 2.8)

project (ParMooN)

include(${CMAKE_SOURCE_DIR}/UserConfig.cmake)

set (ParMooN_VERSION_MAJOR 1)
set (ParMooN_VERSION_MINOR 0)

SET(BUILD_SHARED_LIBS ON)

# set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/intel/mkl/lib/intel64)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/intel/mkl/lib/intel64)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/software/lib)
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{HOME}/Software/lib)

set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/local/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} /opt/intel/mkl/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/intel/mkl/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/software/include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} $ENV{HOME}/Software/include)

# set(ENV{PKG_CONFIG_PATH} $ENV{PKG_CONFIG_PATH}:/opt/local/lib/pkgconfig )
set(PARMOON_EXTLIB_PATH ${PROJECT_SOURCE_DIR}/EXT_LIB)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


if(NOT "${AParMooN_ARCH}" STREQUAL "CRAY64")
# search the blas in the system, if available it will be default
message("==================Searching BLAS in the System - begin===================================================")

if(AParMooN_USE_SYSTEM_BLAS)
  find_package(BLAS)
  include_directories(${BLAS_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    if(BLAS_FOUND)
      message("Available system BLAS_LIB =${BLAS_LIBRARIES}")
    elseif(NOT BLAS_FOUND)
      message("BLAS_NOT_FOUND in your system, Searching in ParMooN")
    endif()
  endif()
endif(AParMooN_USE_SYSTEM_BLAS)

if(NOT BLAS_FOUND)
    message("BLAS_NOT_FOUND in your system, Searching in ParMooN")
    find_package(ACMLBLAS REQUIRED)
    include_directories(${BLAS_INCLUDE_DIRS})
    if(CMAKE_VERBOSE_MAKEFILE)
      if(BLAS_FOUND)
        message("Available system BLAS_LIB =${BLAS_LIBRARIES}")
      elseif(NOT BLAS_FOUND)
        message("BLAS_NOT_FOUND in ParMooN")
      endif(BLAS_FOUND)
    endif(CMAKE_VERBOSE_MAKEFILE)        
endif(NOT BLAS_FOUND)

message("==================Searching BLAS in the System - end=====================================================")

message("BLAS_LIBRARIES  = ${BLAS_LIBRARIES}")


if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")


#   message("==================Searching MPI in the System - begin===================================================")
find_package(MPI)

   if("${AParMooN_MPI_IMPLEMENTATION}" STREQUAL "MPICH")
       set(MPI_LIBS -lmpich -lmpichf90 -lmpichcxx -lgfortran )
   elseif("${AParMooN_MPI_IMPLEMENTATION}" STREQUAL "OPENMPI")
      set(MPI_LIBS -lmpi -lmpifort -lmpicxx)
   elseif("${AParMooN_MPI_IMPLEMENTATION}" STREQUAL "INTELMPI")
      set(MPI_LIBS -lifcore)    
   elseif("${AParMooN_MPI_IMPLEMENTATION}" STREQUAL "MACMPI")
     set(MPI_LIBS -lmpich -lmpichf90 -lmpichcxx)             
   endif("${AParMooN_MPI_IMPLEMENTATION}" STREQUAL "MPICH")   

message(" MPI_LIBS " = ${MPI_LIBS})

 # search the MUMPS in the system, if available it wil be default
 message("==================Searching MUMPS in the System - begin=====================================================")
 find_package(MUMPS)
 include_directories(${MUMPS_INCLUDE_DIRS})
 if(CMAKE_VERBOSE_MAKEFILE)
  if(MUMPS_FOUND)
    message("Available system MUMPS_LIB =${MUMPS_LIBRARIES}")
  elseif(NOT MUMPS_FOUND)
    message(FATAL_ERROR "MUMPS_NOT_FOUND in your system or in the ParMooN")
  endif()
 endif()
 message("==================Searching MUMPS in the System - end=======================================================")

endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")

endif(NOT "${AParMooN_ARCH}" STREQUAL "CRAY64")

# search the UMFPACK in the system, if available it wil be default
message("==================Searching UMFPACK in the System - begin===================================================")
find_package(UMFPACK)
include_directories(${UMFPACK_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(UMFPACK_FOUND)
    message("Available system UMFPACK_LIB =${UMFPACK_LIBRARIES}")
  elseif(NOT UMFPACK_FOUND)
    message(FATAL_ERROR "UMFPACK_NOT_FOUND in your system or in the ParMooN")
  endif(UMFPACK_FOUND)
endif(CMAKE_VERBOSE_MAKEFILE)
message("==================Searching UMFPACK in the System - end=====================================================")

# search the GRIDGEN in the system, if available it wil be default
message("==================Searching GRIDGEN in the System - begin==================================================")
find_package(GRIDGEN)
include_directories(${GRIDGEN_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(GRIDGEN_FOUND)
    message("Available system GRIDGEN_LIB =${GRIDGEN_LIBRARIES}")
  elseif(NOT GRIDGEN_FOUND)
    message(FATAL_ERROR "GRIDGEN_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching GRIDGEN in the System - end====================================================")

# search the TETGEN in the system, if available it wil be default
if("${AParMooN_GEO}" STREQUAL "3D")
  message("==================Searching TETGEN in the System - begin==================================================")
  find_package(TETGEN)
  include_directories(${TETGEN_INCLUDE_DIRS})
  if(CMAKE_VERBOSE_MAKEFILE)
    if(TETGEN_FOUND)
      message("Available system TETGEN_LIB =${TETGEN_LIBRARIES}")
    elseif(NOT TETGEN_FOUND)
      message(FATAL_ERROR "TETGEN_NOT_FOUND in your system or in the ParMooN")
    endif()
  endif()
message("==================Searching TETGEN in the System - end======================================================")
endif("${AParMooN_GEO}" STREQUAL "3D")

if (NOT CONFIGURED_ONCE)
    set(CMAKE_CXX_FLAGS "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
    set(CMAKE_C_FLAGS   "${warnings}"
        CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

# search the TECPLOT in the system, if available it wil be default
message("==================Searching TECPLOT in the System - begin==================================================")
find_package(TECPLOT)
include_directories(${TECPLOT_INCLUDE_DIRS})
if(CMAKE_VERBOSE_MAKEFILE)
  if(TECPLOT_FOUND)
    message("Available system TECPLOT_LIB =${TECPLOT_LIBRARIES}")
  elseif(NOT TECPLOT_FOUND)
    message(FATAL_ERROR "TECPLOT_NOT_FOUND in your system or in the ParMooN")
  endif()
endif()
message("==================Searching TECPLOT in the System - end====================================================")


# select the compiler type,
# CMAKE_BUILD_TYPE [DEBUG|RELEASE|RELWITHDEBINFO|MINSIZEREL]
set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "options")
set(PARMOON_INCLUDE )

include (CMakeForceCompiler)
  set(CMAKE_SYSTEM_NAME ${AParMooN_ARCH})
  
if("${AParMooN_ARCH}" STREQUAL "MAC64")

  if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
    find_path(MPI_CPATH   mpicc-mpich-mp)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc-mpich-mp)
    find_path(MPI_CXXPATH   mpicxx-mpich-mp)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx-mpich-mp)      
  elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
    find_path(CPATH  clang)
    set(CMAKE_C_COMPILER ${CPATH}/clang  )
    find_path(CXXPATH clang++)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/clang++)  
  endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
 
elseif("${AParMooN_ARCH}" STREQUAL "LINUX64")

  if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
    find_path(MPI_CPATH   mpicc)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpicc)
    find_path(MPI_CXXPATH   mpicxx)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpicxx)      
  elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
    find_path(CPATH  gcc)
    set(CMAKE_C_COMPILER ${CPATH}/gcc)
    find_path(CXXPATH g++)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/g++)  
  endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
  
elseif("${AParMooN_ARCH}" STREQUAL "INTEL64")

  if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
    find_path(MPI_CPATH   mpiicc)
    set(CMAKE_C_COMPILER ${MPI_CPATH}/mpiicc)
    find_path(MPI_CXXPATH   mpiicpc)  
    set(CMAKE_CXX_COMPILER ${MPI_CXXPATH}/mpiicpc)      
  elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
  
    include_directories("/opt/intel/include")
  
    find_path(CPATH  icc)
    set(CMAKE_C_COMPILER ${CPATH}/icc)
    find_path(CXXPATH icpc)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/icpc)  
  endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")  
  
# message("MPI_CXXPATH  = ${CXXPATH}")

elseif("${AParMooN_ARCH}" STREQUAL "CRAY64")

    find_path(CPATH cc)
    set(CMAKE_C_COMPILER ${CPATH}/cc)
    find_path(CXXPATH CC)  
    set(CMAKE_CXX_COMPILER ${CXXPATH}/CC)  

endif()


# set the define for compiler  options
  set(CMAKE_C_FLAGS "${PARMOON_C_DEF}")
  set(CMAKE_C_FLAGS_DEBUG "-g ${PARMOON_C_DEF}")
  set(CMAKE_C_FLAGS_RELEASE "O3 -s  ${PARMOON_C_DEF}")
  set(CMAKE_CXX_FLAGS "${PARMOON_CXX_DEF} " )
  set(CMAKE_CXX_FLAGS_DEBUG "${PARMOON_CXX_DEF} ")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -s ${PARMOON_CXX_DEF}") 
 
include_directories("${CMAKE_SOURCE_DIR}/include/General")
include_directories("${CMAKE_SOURCE_DIR}/include/FE")
include_directories("${CMAKE_SOURCE_DIR}/include/Geometry")
include_directories("${CMAKE_SOURCE_DIR}/include/QuadFormulas")
include_directories("${CMAKE_SOURCE_DIR}/include/FE1D")
include_directories("${CMAKE_SOURCE_DIR}/include/FE2D")
include_directories("${CMAKE_SOURCE_DIR}/include/FE3D")
include_directories("${CMAKE_SOURCE_DIR}/include/AMG")
include_directories("${CMAKE_SOURCE_DIR}/include/PBE")
include_directories("${CMAKE_SOURCE_DIR}/include/Refinement")
include_directories("${CMAKE_SOURCE_DIR}/include/System")
include_directories("${CMAKE_SOURCE_DIR}/include/Parallel")
include_directories("${CMAKE_SOURCE_DIR}/Examples")

add_subdirectory(src/AMG) 
add_subdirectory(src/Refinement) 
add_subdirectory(src/QuadFormulas)  
add_subdirectory(src/Geometry)          
add_subdirectory(src/General) 
add_subdirectory(src/FE) 
add_subdirectory(src/System) 
# add_subdirectory(src/PBE) 


message("CMAKE_SOURCE_DIR  = ${CMAKE_SOURCE_DIR}")


if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
   add_subdirectory(src/Parallel) 
endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")

add_executable(parmoon_${AParMooN_GEO}_${AParMooN_PARALLEL_TYPE}.exe  ${AParMooN_MODEL})

if( "${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI" )
   target_link_libraries(parmoon_${AParMooN_GEO}_${AParMooN_PARALLEL_TYPE}.exe "sys" "par" "geo"  "fe" "gen" "geo" "fe" "gen"  
                      "amg"  "quad" "ref" "geo" "fe" "quad" "geo"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${MPIBLACS_LIBRARIES}  ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES}  ${MUMPS_LIBRARIES}   ${BLAS_LIBRARIES} ${MPI_LIBS})  
elseif( "${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
   target_link_libraries(parmoon_${AParMooN_GEO}_${AParMooN_PARALLEL_TYPE}.exe "sys" "geo"  "fe" "gen" "geo" "fe" "gen"  
                      "amg"  "quad" "ref" "geo" "fe" "quad" "geo"  ${UMFPACK_LIBRARIES} ${TECPLOT_LIBRARIES}
                      ${BLAS_LIBRARIES} ${TETGEN_LIBRARIES} ${GRIDGEN_LIBRARIES} )
endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")

message("============================================================================================================"  )  
message("AParMooN_GEO=${AParMooN_GEO}")
message("AParMooN_MODEL=${AParMooN_MODEL}")
message("AParMooN_ARCH=${AParMooN_ARCH}")
message("AParMooN_PARALLEL_TYPE=${AParMooN_PARALLEL_TYPE}")
message("C_COMPILER=${CMAKE_C_COMPILER}")
message("CXX_COMPILER=${CMAKE_CXX_COMPILER}")
message("PARMOON_CXX_DEF=${PARMOON_CXX_DEF}")
message("EXE_FILE=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/parmoon_${AParMooN_GEO}_${AParMooN_PARALLEL_TYPE}.exe")
message("============================================================================================================"  )        
message("Configuration for the ParMooN Version 1.1")  
message("         Sashikumaar Ganesan")
message("         Dept. of Computational & Data Sciences")
message("         Indian Institute of Science, Bangalore, India")  
message("mailt0 : sashi-AT-iisc.ac.in") 
message("============================================================================================================"  ) 




