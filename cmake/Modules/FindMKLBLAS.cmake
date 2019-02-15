# ===================================================================
# This is FindBLAS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, CDS, IISc Bangalore, India
# date: 14 Feb 2018
# searching for MKL BLAS lib in the system 
# if found, this will define
#  BLAS_FOUND - System has BLAS
#  BLAS_INCLUDE_DIRS - The BLAS include directories
#  BLAS_LIBRARIES - The libraries needed to use BLAS
# ===================================================================
if(BLAS_INCLUDES AND BLAS_LIBRARIES)
  set(BLAS_FIND_QUIETLY TRUE)
endif(BLAS_INCLUDES AND BLAS_LIBRARIES)

if(NOT BLAS_FOUND)
 
  find_path(BLAS_INCLUDE_DIR   mkl.h PATHS $ENV{ACMBLASDIR}/include ${CMAKE_INCLUDE_PATH})
  
  if(NOT BLAS_INCLUDE_DIR)
    find_path(BLAS_INCLUDE_DIR   mkl.fi PATHS $ENV{ACMBLASDIR}/include ${CMAKE_INCLUDE_PATH})
  endif(NOT BLAS_INCLUDE_DIR)
  
  find_library(BLAS_LIBRARY NAMES mkl_core PATHS $ENV{BLASDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(BLAS_LIBDIR ${BLAS_LIBRARY} PATH)
  find_library(BLAS_LIBRARY_MP NAMES mkl_intel_lp64 PATHS ${BLAS_LIBDIR})
  find_library(BLAS_LIBRARY_SQ NAMES mkl_sequential PATHS ${BLAS_LIBDIR})
  
  if(BLAS_LIBRARY)  
    # combine mumps and its deps 
    if(BLAS_LIBRARY_SQ)
      set(BLAS_LIBRARY   ${BLAS_LIBRARY_SQ} ${BLAS_LIBRARY}  ) 
    else(BLAS_LIBRARY_SQ)   
      set(BLAS_LIBRARY FALSE)
    endif(BLAS_LIBRARY_SQ)    
    
    if(BLAS_LIBRARY_MP)
      set(BLAS_LIBRARY   ${BLAS_LIBRARY_MP} ${BLAS_LIBRARY}  ) 
    else(BLAS_LIBRARY_MP)   
      set(BLAS_LIBRARY FALSE)
    endif(BLAS_LIBRARY_MP)
    
    # set BLAS
    if(BLAS_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(BLAS_LIBRARIES ${BLAS_LIBRARY})
      set(BLAS_INCLUDE_DIRS ${BLAS_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set BLAS_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(BLAS  DEFAULT_MSG
                                        BLAS_LIBRARY BLAS_INCLUDE_DIR)

      mark_as_advanced(BLAS_INCLUDE_DIR BLAS_LIBRARY)
    endif(BLAS_LIBRARY)  
  endif()

endif()


