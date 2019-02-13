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
         -I../../../include/PBE -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../include/System \
         -I../../../include/Parallel  -I../../../EXT_LIB/MUMPS/include \
         -I../../../EXT_LIB/tecplot/include

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_System_par.a
   MEMBERS = \
        $(SRCPATH)Assemble.o \
        $(SRCPATH)AssembleMat3D.o \
        $(SRCPATH)SystemCD3D.o \
        $(SRCPATH)SystemTCD3D.o \
        $(SRCPATH)SystemNSE3D.o \
        $(SRCPATH)SystemTNSE3D.o \
        $(SRCPATH)SystemTNSE3D_ALE.o 
###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

