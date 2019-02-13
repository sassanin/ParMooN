# ===================================================================
# This is FindTETGEN.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a TETGEN lib in the system 
# if found, this will define
#  TETGEN_FOUND - System has TETGEN
#  TETGEN_INCLUDE_DIRS - The TETGEN include directories
#  TETGEN_LIBRARIES - The libraries needed to use TETGEN
# ===================================================================
if(TETGEN_INCLUDES AND TETGEN_LIBRARIES)
  set(TETGEN_FIND_QUIETLY TRUE)
endif(TETGEN_INCLUDES AND TETGEN_LIBRARIES)

if(NOT TETGEN_FOUND)
 if(USE_SYSTEM_TETGEN)  
  find_path(TETGEN_INCLUDE_DIR   tetgen.h PATHS $ENV{TETGENDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(TETGEN_LIBRARY NAMES tet PATHS $ENV{TETGENDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(TETGEN_LIBDIR ${TETGEN_LIBRARY} PATH)
 endif(USE_SYSTEM_TETGEN)   
 
  if(NOT TETGEN_LIBRARY)
    message("TETGEN not found in the system, so checking the availability in ParMooN for the selected AParMooN_ARCH=${AParMooN_ARCH}")
    find_path(TETGEN_INCLUDE_DIR  tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
    find_library(TETGEN_LIBRARY NAMES tet_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
  endif(NOT TETGEN_LIBRARY)
  
  if(TETGEN_LIBRARY)      
    # set TETGEN
    if(TETGEN_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(TETGEN_LIBRARIES ${TETGEN_LIBRARY})
      set(TETGEN_INCLUDE_DIRS ${TETGEN_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set TETGEN_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(TETGEN  DEFAULT_MSG
                                        TETGEN_LIBRARY TETGEN_INCLUDE_DIR)

      mark_as_advanced(TETGEN_INCLUDE_DIR TETGEN_LIBRARY)
    endif(TETGEN_LIBRARY)  
  endif()

endif()


