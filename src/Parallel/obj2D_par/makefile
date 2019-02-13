###################################
#
#  Makefile for MooNMD
#  author: Sashikumaar Ganesan
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = \
         -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE1D -I../../../include/FE2D \
         -I../../../include/AMG  -I../../../include/Parallel  \
         -I../../../EXT_LIB/Metis -I../../../EXT_LIB/Metis/MetisInclude\
         -I../../../EXT_LIB/Metis/Lib/ -I../../../EXT_LIB/MUMPS/include 

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_par.a
   MEMBERS = \
    $(SRCPATH)MeshPartition.o \
    $(SRCPATH)ParFECommunicator2D.o 
# \
# 	$(SRCPATH)MumpsSolver.o  \
# 	$(SRCPATH)ParFECommunicator2D.o  \
# 	$(SRCPATH)ParDirectSolver.o \
# 	$(SRCPATH)ParMooNBLAS.o \
# 	$(SRCPATH)ParVector.o \
# 	$(SRCPATH)ParVectorNSE.o 

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache
