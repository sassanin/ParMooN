# ===================================================================
# This is FindGRIDGEN.cmake file for the ParMooN Version 1.1
# written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# date: 07 June 2015
# searching for a GRIDGEN lib in the system 
# if found, this will define
#  GRIDGEN_FOUND - System has GRIDGEN
#  GRIDGEN_INCLUDE_DIRS - The GRIDGEN include directories
#  GRIDGEN_LIBRARIES - The libraries needed to use GRIDGEN
# ===================================================================
if(GRIDGEN_INCLUDES AND GRIDGEN_LIBRARIES)
  set(GRIDGEN_FIND_QUIETLY TRUE)
endif(GRIDGEN_INCLUDES AND GRIDGEN_LIBRARIES)

if(NOT GRIDGEN_FOUND)
 if(USE_SYSTEM_GRIDGEN) 
  find_path(GRIDGEN_INCLUDE_DIR   gridgen.h PATHS $ENV{GRIDGENDIR}/include ${CMAKE_INCLUDE_PATH})
  find_library(GRIDGEN_LIBRARY NAMES gridgen PATHS $ENV{GRIDGENDIR}/lib ${CMAKE_LIBRARY_PATH})
  get_filename_component(GRIDGEN_LIBDIR ${GRIDGEN_LIBRARY} PATH)
 endif(USE_SYSTEM_GRIDGEN)    
  if(NOT GRIDGEN_LIBRARY)
    message("GRIDGEN not found in the system, so checking the availability in ParMooN for the selected AParMooN_ARCH=${AParMooN_ARCH}")
    find_path(GRIDGEN_INCLUDE_DIR  gridgen.h PATHS ${PARMOON_EXTLIB_PATH}/GridGen)
    find_library(GRIDGEN_LIBRARY NAMES gridgen_${AParMooN_ARCH} PATHS ${PARMOON_EXTLIB_PATH}/GridGen)
  endif(NOT GRIDGEN_LIBRARY)
  
  if(GRIDGEN_LIBRARY)      
    # set GRIDGEN
    if(GRIDGEN_LIBRARY)
      include(FindPackageHandleStandardArgs)
    
      set(GRIDGEN_LIBRARIES ${GRIDGEN_LIBRARY})
      set(GRIDGEN_INCLUDE_DIRS ${GRIDGEN_INCLUDE_DIR})

      # handle the QUIETLY and REQUIRED arguments and set GRIDGEN_FOUND to TRUE
      # if all listed variables are TRUE
      find_package_handle_standard_args(GRIDGEN  DEFAULT_MSG
                                        GRIDGEN_LIBRARY GRIDGEN_INCLUDE_DIR)

      mark_as_advanced(GRIDGEN_INCLUDE_DIR GRIDGEN_LIBRARY)
    endif(GRIDGEN_LIBRARY)  
  endif()

endif()


