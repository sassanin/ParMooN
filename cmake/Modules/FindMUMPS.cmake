# ===================================================================
# This is FindMUMPS.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a MUMPS lib in the system 
# if found, this will define
#  MUMPS_FOUND - System has MUMPS
#  MUMPS_INCLUDE_DIRS - The MUMPS include directories
#  MUMPS_LIBRARIES - The libraries needed to use MUMPS
# ===================================================================
if(MUMPS_INCLUDES AND MUMPS_LIBRARIES)
  set(MUMPS_FIND_QUIETLY TRUE)
endif(MUMPS_INCLUDES AND MUMPS_LIBRARIES)

if(NOT MUMPS_FOUND)
 if(USE_SYSTEM_MUMPS)  
  find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS $ENV{MUMPSDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(MUMPS_LIBRARY NAMES dmumps PATHS $ENV{MUMPSDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(MUMPS_LIBDIR ${MUMPS_LIBRARY} PATH)
  find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common PATHS ${MUMPS_LIBDIR})
 endif(USE_SYSTEM_MUMPS)
 
  if(NOT MUMPS_LIBRARY)
    message("MUMPS not found in the system, so checking the availability in ParMooN for the selected AParMooN_ARCH=${AParMooN_ARCH}")
    find_path(MUMPS_INCLUDE_DIR  mumps_compat.h PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/include)
    find_library(MUMPS_LIBRARY NAMES dmumps_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/MUMPS/lib) 
    get_filename_component(MUMPS_LIBDIR ${MUMPS_LIBRARY} PATH)
    find_library(MUMPS_LIBRARY_COMMON NAMES mumps_common_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR})
  endif(NOT MUMPS_LIBRARY)
  
  if(MUMPS_LIBRARY)  

    # deps for mumps
    find_library(SCALAPACK_LIB NAMES scalapack PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})
    find_library(PARMETIS_LIB NAMES parmetis PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH} ${PARMOON_EXTLIB_PATH}/Metis)
    find_library(METIS_LIB NAMES metis PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
   
    find_library(PORD_LIB NAMES pord_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
    find_library(PTESMUMPS_LIB NAMES ptesmumps_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
    find_library(PTSCOTCH_LIB NAMES ptscotch_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
    find_library(SCOTCH_LIB NAMES scotch_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
    find_library(PTSCOTCHERREXIT_LIB NAMES ptscotcherrexit_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
    find_library(PTSCOTCHPARMETIS_LIB NAMES ptscotchparmetis_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  

#     find_library(SCOTCHERR_LIB NAMES scotcherr_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
#     find_library(SCOTCHERREXIT_LIB NAMES scotcherrexit_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})  
#     find_library(PTSCOTCHERR_LIB NAMES ptscotcherr_${AParMooN_ARCH} PATHS ${MUMPS_LIBDIR} ${CMAKE_LIBRARY_PATH})          
#                         
            
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
    if(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
      set(MUMPS_LIBRARY  ${MUMPS_LIBRARY} ${MUMPS_LIBRARY_COMMON} ${SCALAPACK_LIB} ${PARMETIS_LIB} ${METIS_LIB} ${PORD_LIB} ${PTESMUMPS_LIB} ${PTSCOTCH_LIB} ${SCOTCH_LIB} ${PTSCOTCHERREXIT_LIB} ${PTSCOTCHPARMETIS_LIB} ${SCOTCHERR_LIB} ${SCOTCHERREXIT_LIB} ${PTSCOTCHERR_LIB}) 
    else(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
      if(NOT SCALAPACK_LIB)
       message("MUMPS found in system but scalapack not found in the search path")
      endif()
       if(NOT PARMETIS_LIB)
       message("MUMPS found in system but parmetis not found in the search path")
      endif()
      if(NOT METIS_LIB)
       message("MUMPS found in system but metis not found in the search path")
      endif()
      set(MUMPS_LIBRARY FALSE)
    endif(MUMPS_LIBRARY_COMMON AND SCALAPACK_LIB AND METIS_LIB)
    
    # set MUMPS
    if(MUMPS_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(MUMPS_LIBRARIES ${MUMPS_LIBRARY})
      set(MUMPS_INCLUDE_DIRS ${MUMPS_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set MUMPS_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(MUMPS  DEFAULT_MSG
                                        MUMPS_LIBRARY MUMPS_INCLUDE_DIR)

      mark_as_advanced(MUMPS_INCLUDE_DIR MUMPS_LIBRARY)
    endif(MUMPS_LIBRARY)  
  endif()

endif()


