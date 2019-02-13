# # ===================================================================
# # This is FindPARDISO.cmake file for the ParMooN Version 1.1
# # written by Sashikumaar Ganesan, SERC, IISc Bangalore, India
# # date: 17 June 2015
# # searching for a PARDISO lib in the system 
# # if found, this will define
# #  PARDISO_FOUND - System has PARDISO
# #  PARDISO_INCLUDE_DIRS - The PARDISO include directories
# #  PARDISO_LIBRARIES - The libraries needed to use PARDISO
# # ===================================================================
# if(PARDISO_INCLUDES AND PARDISO_LIBRARIES)
#   set(PARDISO_FIND_QUIETLY TRUE)
# endif(PARDISO_INCLUDES AND PARDISO_LIBRARIES)
# 
# if(NOT PARDISO_FOUND)
#  if(USE_SYSTEM_PARDISO)  
# #   find_path(PARDISO_INCLUDE_DIR   tetgen.h PATHS $ENV{PARDISODIR}/include ${CMAKE_INCLUDE_PATH})
#   find_library(PARDISO_LIBRARY NAMES tet PATHS $ENV{PARDISODIR}/lib ${CMAKE_LIBRARY_PATH})
#   get_filename_component(PARDISO_LIBDIR ${PARDISO_LIBRARY} PATH)
#  endif(USE_SYSTEM_PARDISO)   
#  
#   if(NOT PARDISO_LIBRARY)
#     message("PARDISO not found in the system, so checking the availability in ParMooN for the selected ARCH=${ARCH}")
# #     find_path(PARDISO_INCLUDE_DIR  tetgen.h PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
#     find_library(PARDISO_LIBRARY NAMES tet_${ARCH} PATHS ${PARMOON_EXTLIB_PATH}/tetgen)
#   endif(NOT PARDISO_LIBRARY)
#   
#   if(PARDISO_LIBRARY)      
#     # set PARDISO
#     if(PARDISO_LIBRARY)
#       include(FindPackageHandleStandardArgs)
#     
#       set(PARDISO_LIBRARIES ${PARDISO_LIBRARY})
#       set(PARDISO_INCLUDE_DIRS ${PARDISO_INCLUDE_DIR})
# 
#       # handle the QUIETLY and REQUIRED arguments and set PARDISO_FOUND to TRUE
#       # if all listed variables are TRUE
#       find_package_handle_standard_args(PARDISO  DEFAULT_MSG
#                                         PARDISO_LIBRARY PARDISO_INCLUDE_DIR)
# 
#       mark_as_advanced(PARDISO_INCLUDE_DIR PARDISO_LIBRARY)
#     endif(PARDISO_LIBRARY)  
#   endif()
# 
# endif()
# 
# 
