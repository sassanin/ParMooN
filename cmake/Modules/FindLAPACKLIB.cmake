# ===================================================================
# This is FindLAPACK.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 17 June 2015
# searching for a LAPACK lib in the system 
# if found, this will define
#  LAPACK_FOUND - System has LAPACK
#  LAPACK_INCLUDE_DIRS - The LAPACK include directories
#  LAPACK_LIBRARIES - The libraries needed to use LAPACK
# ===================================================================
if(LAPACK_INCLUDES AND LAPACK_LIBRARIES)
  set(LAPACK_FIND_QUIETLY TRUE)
endif(LAPACK_INCLUDES AND LAPACK_LIBRARIES)

if(NOT LAPACK_FOUND)
    find_library(SCALAPACK_LIBRARY NAMES scalapack_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS) 
    find_library(BLACS_LIBRARY NAMES blacs_MPI_${AParMooN_ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS)     
    find_library(BLACS_F77LIBRARY NAMES blacsF77init_MPI_${AParMooN_ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS) 
    
#         message("SCALAPACK_LIBRARY: ${SCALAPACK_LIBRARY}")
#          message("BLACS_LIBRARY: ${BLACS_LIBRARY}")       
#         message("BLACS_F77LIBRARY: ${BLACS_F77LIBRARY}")
  
    if(SCALAPACK_LIBRARY)
      set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${SCALAPACK_LIBRARY} ) 
    else(SCALAPACK_LIBRARY)   
      set(LAPACK_LIBRARY FALSE)
       message(FATAL_ERROR "SCALAPACK_LIBRARY: ${SCALAPACK_LIBRARY} not found !") 
    endif(SCALAPACK_LIBRARY)
    
    if(BLACS_LIBRARY)
      set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${BLACS_LIBRARY} ) 
    else(BLACS_LIBRARY)   
      set(LAPACK_LIBRARY FALSE)
    endif(BLACS_LIBRARY)    
    
     if(BLACS_F77LIBRARY)
      set(LAPACK_LIBRARY  ${LAPACK_LIBRARY}  ${BLACS_F77LIBRARY} ) 
    else(BLACS_F77LIBRARY)   
      set(LAPACK_LIBRARY FALSE)
    endif(BLACS_F77LIBRARY)   
    

    
    # set LAPACK
        
    if(LAPACK_LIBRARY)    
       include(FindPackageHandleStandardArgs)
    
      set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
      set(LAPACK_INCLUDE_DIR ${LAPACK_LIBRARY}) 
      set(LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set LAPACK_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(LAPACK  DEFAULT_MSG
                                        LAPACK_LIBRARY LAPACK_INCLUDE_DIR)

      mark_as_advanced(LAPACK_INCLUDE_DIR LAPACK_LIBRARY)
    endif(LAPACK_LIBRARY)  

endif(NOT LAPACK_FOUND)



