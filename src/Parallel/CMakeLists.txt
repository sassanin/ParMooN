
include_directories("${CMAKE_SOURCE_DIR}/include/Parallel")
include_directories("${CMAKE_SOURCE_DIR}/EXT_LIB/Metis/Lib")


if("${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
 
     list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MumpsSolver.C")
     list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/SeqParDirectSolver.C")

     add_library(par STATIC ${PAR_SOURCES})

elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI"  OR "${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID" )

   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MeshPartition.C") 
   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/MumpsSolver.C")
   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFEMapper3D.C")
   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParFECommunicator3D.C")
   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDirectSolver.C")
   list(APPEND PAR_SOURCES "${PROJECT_SOURCE_DIR}/src/Parallel/ParDiso.C")

   add_library(par STATIC ${PAR_SOURCES})

endif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
 
