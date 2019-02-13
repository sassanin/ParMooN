###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/AMG

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_amg.a
   MEMBERS = \
	$(SRCPATH)amg_1d_prec.o \
	$(SRCPATH)amg_2d_prec.o \
	$(SRCPATH)amg_blas.o \
	$(SRCPATH)amg_blas_nse2d.o \
	$(SRCPATH)amg_blas_nse3d.o \
	$(SRCPATH)amg_coarsen.o \
	$(SRCPATH)amg_iter.o \
	$(SRCPATH)amg_low.o \
	$(SRCPATH)amg_solve_main.o \
	$(SRCPATH)amg_solvers.o \
	$(SRCPATH)amg_sp.o \
	$(SRCPATH)amg_variables.o 

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY)
