# ===================================================================
# This is FindMPIBLACS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 17 June 2015
# searching for a MPIBLACS lib in the system 
# if found, this will define
#  MPIBLACS_FOUND - System has MPIBLACS
#  MPIBLACS_INCLUDE_DIRS - The MPIBLACS include directories
#  MPIBLACS_LIBRARIES - The libraries needed to use MPIBLACS
# ===================================================================
if(MPIBLACS_INCLUDES AND MPIBLACS_LIBRARIES)
  set(MPIBLACS_FIND_QUIETLY TRUE)
endif(MPIBLACS_INCLUDES AND MPIBLACS_LIBRARIES)

if(NOT MPIBLACS_FOUND)
    find_library(MPIBLACS_LIBRARY NAMES scalapack_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS) 
    find_library(BLACS_LIBRARY NAMES blacs_MPI_${AParMooN_ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS)     
    find_library(BLACS_F77LIBRARY NAMES blacsF77init_MPI_${AParMooN_ARCH}-0 PATHS ${PARMOON_EXTLIB_PATH}/MPIBLACS) 
    
#         message("MPIBLACS_LIBRARY: ${MPIBLACS_LIBRARY}")
#          message("BLACS_LIBRARY: ${BLACS_LIBRARY}")       
#         message("BLACS_F77LIBRARY: ${BLACS_F77LIBRARY}")
  
    if(MPIBLACS_LIBRARY)
      set(MPIBLACS_LIBRARY  ${MPIBLACS_LIBRARY}  ${MPIBLACS_LIBRARY} ) 
    else(MPIBLACS_LIBRARY)   
      set(MPIBLACS_LIBRARY FALSE)
       message(FATAL_ERROR "MPIBLACS_LIBRARY: ${MPIBLACS_LIBRARY} not found !") 
    endif(MPIBLACS_LIBRARY)
    
    if(BLACS_LIBRARY)
      set(MPIBLACS_LIBRARY  ${MPIBLACS_LIBRARY}  ${BLACS_LIBRARY} ) 
    else(BLACS_LIBRARY)   
      set(MPIBLACS_LIBRARY FALSE)
    endif(BLACS_LIBRARY)    
    
     if(BLACS_F77LIBRARY)
      set(MPIBLACS_LIBRARY  ${MPIBLACS_LIBRARY}  ${BLACS_F77LIBRARY} ) 
    else(BLACS_F77LIBRARY)   
      set(MPIBLACS_LIBRARY FALSE)
    endif(BLACS_F77LIBRARY)   
    

    
    # set MPIBLACS
        
    if(MPIBLACS_LIBRARY)    
       include(FindPackageHandleStandardArgs)
    
      set(MPIBLACS_LIBRARIES ${MPIBLACS_LIBRARY})
      set(MPIBLACS_INCLUDE_DIR ${MPIBLACS_LIBRARY}) 
      set(MPIBLACS_INCLUDE_DIRS ${MPIBLACS_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set MPIBLACS_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(MPIBLACS  DEFAULT_MSG
                                        MPIBLACS_LIBRARY MPIBLACS_INCLUDE_DIR)

      mark_as_advanced(MPIBLACS_INCLUDE_DIR MPIBLACS_LIBRARY)
    endif(MPIBLACS_LIBRARY)  

endif(NOT MPIBLACS_FOUND)



