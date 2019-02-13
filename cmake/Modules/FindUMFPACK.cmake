# ===================================================================
# This is FindUMFPACK.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a UMFPACK lib in the system 
# if found, this will define
#  UMFPACK_FOUND - System has UMFPACK
#  UMFPACK_INCLUDE_DIRS - The UMFPACK include directories
#  UMFPACK_LIBRARIES - The libraries needed to use UMFPACK
# ===================================================================
if(UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif(UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)

if(NOT UMFPACK_FOUND)
 if(USE_SYSTEM_UMFPACK)   
  find_path(UMFPACK_INCLUDE_DIR   umfpack.h PATHS $ENV{UMFPACKDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(UMFPACK_LIBRARY NAMES umfpack PATHS $ENV{UMFPACKDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
  find_library(UMFPACK_LIBRARY_SUITESE NAMES suitesparseconfig PATHS ${UMFPACK_LIBDIR})
  find_library(UMFPACK_LIBRARY_AMD NAMES amd PATHS ${UMFPACK_LIBDIR})
  find_library(UMFPACK_LIBRARY_CHOLMOD NAMES cholmod PATHS ${UMFPACK_LIBDIR}) 
  find_library(UMFPACK_LIBRARY_CCOLMOD NAMES ccholmod PATHS ${UMFPACK_LIBDIR})   
  find_library(UMFPACK_LIBRARY_CAMD NAMES camd PATHS ${UMFPACK_LIBDIR})
  find_library(UMFPACK_LIBRARY_COLAMD NAMES colamd PATHS ${UMFPACK_LIBDIR})
 endif(USE_SYSTEM_UMFPACK)      
  if(NOT UMFPACK_LIBRARY)
    message("UMFPACK not found in the system, so checking the availability in ParMooN for the selected AParMooN_ARCH=${AParMooN_ARCH}")
    find_path(UMFPACK_INCLUDE_DIR  umfpack.h PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Include)
    find_library(UMFPACK_LIBRARY NAMES umfpack_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/UMFPACK/Lib) 
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
    find_library(UMFPACK_LIBRARY_SUITESE NAMES suitesparseconfig_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR})
    find_library(UMFPACK_LIBRARY_AMD NAMES amd_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR})
    find_library(UMFPACK_LIBRARY_CHOLMOD NAMES cholmod_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR}) 
    find_library(UMFPACK_LIBRARY_CCOLMOD NAMES ccolamd_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR})   
    find_library(UMFPACK_LIBRARY_CAMD NAMES camd_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR})    
    find_library(UMFPACK_LIBRARY_COLAMD NAMES colamd_${AParMooN_ARCH} PATHS ${UMFPACK_LIBDIR})    
  endif(NOT UMFPACK_LIBRARY)
  
  if(UMFPACK_LIBRARY)  
    # deps for mumps
    find_library(A_LIB NAMES scalapack PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH})
    find_library(PARMETIS_LIB NAMES parmetis PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH} ${PARMOON_EXTLIB_PATH}/Metis)
    find_library(METIS_LIB NAMES metis PATHS ${UMFPACK_LIBDIR} ${CMAKE_LIBRARY_PATH})     
    if(NOT SCALAPACK_LIB)
       find_library(SCALAPACK_LIB NAMES scalapack_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS)
    endif()
     if(NOT PARMETIS_LIB)
         find_library(PARMETIS_LIB NAMES parmetis_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis)
    endif()
    if(NOT METIS_LIB)
      find_library(METIS_LIB NAMES metis_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/Metis)   
    endif()      
    
    # combine mumps and its deps    
#     if(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD AND UMFPACK_LIBRARY_CHOLMOD AND UMFPACK_LIBRARY_CCOLMOD AND UMFPACK_LIBRARY_CAMD AND UMFPACK_LIBRARY_COLAMD)
#       set(UMFPACK_LIBRARY  ${UMFPACK_LIBRARY}  ${UMFPACK_LIBRARY_CHOLMOD}  ${UMFPACK_LIBRARY_CCOLMOD}  ${UMFPACK_LIBRARY_CAMD}  ${UMFPACK_LIBRARY_AMD}  ${UMFPACK_LIBRARY_COLAMD}   ${UMFPACK_LIBRARY_SUITESE}) 
#     endif(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD AND UMFPACK_LIBRARY_CHOLMOD AND UMFPACK_LIBRARY_CCOLMOD AND UMFPACK_LIBRARY_CAMD)
#     
    if(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)
      set(UMFPACK_LIBRARY  ${UMFPACK_LIBRARY} ${UMFPACK_LIBRARY_AMD} ${UMFPACK_LIBRARY_SUITESE} ) 
    else(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)   
      set(UMFPACK_LIBRARY FALSE)
    endif(UMFPACK_LIBRARY_SUITESE AND UMFPACK_LIBRARY_AMD)
    
    # set UMFPACK
    if(UMFPACK_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY})
      set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set UMFPACK_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(UMFPACK  DEFAULT_MSG
                                        UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR)

      mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY)
    endif(UMFPACK_LIBRARY)  
  endif()

endif()


