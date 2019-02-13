###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

include ../../../Make.inc

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/General -I../../../include/Geometry \
         -I../../../include/Refinement -I../../../include/FE \
         -I../../../include/QuadFormulas -I../../../include/AMG \
	 -I../../../include/PBE

LIBRARY = $(LIBPATH)lib_pbe.a
   MEMBERS = \
        $(SRCPATH)aggregation.o \
        $(SRCPATH)basictools.o \
        $(SRCPATH)faltung.o \
        $(SRCPATH)discont.o \
        $(SRCPATH)faltung_hilfe.o \
        $(SRCPATH)fepc_easy.o \
        $(SRCPATH)folge.o \
        $(SRCPATH)folgen_vektor.o \
        $(SRCPATH)funktion.o \
        $(SRCPATH)koeffizienten.o \
        $(SRCPATH)piecewise_lin_cont.o

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache


