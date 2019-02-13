###################################
#
#  Makefile for ParMooN
#  author: Sashikumaar Ganesan
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE1D -I../../../include/FE2D \
         -I../../../include/AMG \
         -I../../../include/System \
         -I../../../include/Parallel \
         -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../EXT_LIB/tecplot/include


include ../../../Make.inc

DEFINEPARAMS = ../../../include/General/DefineParams.h

LIBRARY = $(LIBPATH)lib_System.a
   MEMBERS = \
	$(SRCPATH)SystemNSE2D.o \
	$(SRCPATH)SystemCD2D.o \
	$(SRCPATH)SystemTCD2D.o \
	$(SRCPATH)SystemTCD2D_ALE.o \
	$(SRCPATH)SystemTNSE2D.o \
	$(SRCPATH)SystemTNSE2D_ALE.o 
###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

