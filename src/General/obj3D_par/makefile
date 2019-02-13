###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
LIBPATH = ../../../lib3D/
SRCPATH = 
INC    = -I../../../include/General -I../../../include/Geometry \
         -I../../../include/Refinement -I../../../include/FE \
         -I../../../include/QuadFormulas -I../../../include/AMG\
         -I../../../include/Parallel 

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_general_par.a
   MEMBERS = \
        $(SRCPATH)Database.o \
        $(SRCPATH)MainUtilities.o \
        $(SRCPATH)MooNMD_Io.o \
        $(SRCPATH)ReadParam.o \
        $(SRCPATH)Utilities.o \
        $(SRCPATH)Vector.o

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

ReadParam.o : ../ReadParam.C
	@echo $<
	@$(CPP) $(CPPFLAG_$(MACH)_FAST) $(DEFINE) -c $< -o $@

