###################################
#
#  Makefile for ParMooN
#  author: Sashikumaar Ganesan
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib3D/
INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE2D -I../../../include/FE1D \
         -I../../../include/FE3D -I../../../include/AMG \
         -I../../../EXT_LIB/Metis/MetisInclude -I../../../EXT_LIB/Metis/Lib  \
	  -I../../../EXT_LIB/MUMPS/include -I../../../include/Parallel 

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_par.a

MEMBERS = \
	$(SRCPATH)MeshPartition.o \
	$(SRCPATH)MumpsSolver.o  \
	$(SRCPATH)ParFEMapper3D.o \
	$(SRCPATH)ParFECommunicator3D.o \
	$(SRCPATH)ParDirectSolver.o \
	$(SRCPATH)ParDiso.o

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files 
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

