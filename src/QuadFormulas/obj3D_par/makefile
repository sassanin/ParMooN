###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
LIBPATH = ../../../lib3D/
SRCPATH = 
INC    = -I../../../include/QuadFormulas -I../../../include/General \
         -I../../../include/FE -I../../../include/Geometry

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_QF_par.a
   MEMBERS = \
        $(SRCPATH)QuadFormula.o \
        $(SRCPATH)QuadFormula1D.o \
        $(SRCPATH)QuadFormula2D.o \
        $(SRCPATH)QuadFormula3D.o \
        $(SRCPATH)QuadFormulaHexa.o \
        $(SRCPATH)QuadFormulaQuad.o \
        $(SRCPATH)QuadFormulaTetra.o \
        $(SRCPATH)QuadFormulaTria.o


###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache
