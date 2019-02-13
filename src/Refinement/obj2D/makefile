###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/Refinement -I../../../include/Geometry \
         -I../../../include/General

include ../../../Make.inc

DEFINEPARAMS = ../../../include/General/DefineParams.h

LIBRARY = $(LIBPATH)lib_refinement.a
   MEMBERS = \
	$(SRCPATH)RefLineDesc.o \
	$(SRCPATH)RefMortar0Desc.o \
	$(SRCPATH)RefMortar1Desc.o \
	$(SRCPATH)RefMortarLineDesc.o \
	$(SRCPATH)RefNoRef.o \
	$(SRCPATH)RefQuadBis0Desc.o \
	$(SRCPATH)RefQuadBis1Desc.o \
	$(SRCPATH)RefQuad1Conf0Desc.o \
	$(SRCPATH)RefQuad1Conf1Desc.o \
	$(SRCPATH)RefQuad1Conf2Desc.o \
	$(SRCPATH)RefQuad1Conf3Desc.o \
	$(SRCPATH)RefQuad2Conf0Desc.o \
	$(SRCPATH)RefQuad2Conf1Desc.o \
	$(SRCPATH)RefQuad2Conf2Desc.o \
	$(SRCPATH)RefQuad2Conf3Desc.o \
	$(SRCPATH)RefQuadRegDesc.o \
	$(SRCPATH)RefQuadToTri0Desc.o \
	$(SRCPATH)RefQuadToTri1Desc.o \
	$(SRCPATH)RefTriBis0Desc.o \
	$(SRCPATH)RefTriBis1Desc.o \
	$(SRCPATH)RefTriBis2Desc.o \
	$(SRCPATH)RefTriBis01Desc.o \
	$(SRCPATH)RefTriBis02Desc.o \
	$(SRCPATH)RefTriBis10Desc.o \
	$(SRCPATH)RefTriBis12Desc.o \
	$(SRCPATH)RefTriBis20Desc.o \
	$(SRCPATH)RefTriBis21Desc.o \
	$(SRCPATH)RefTriRegDesc.o \
	$(SRCPATH)RefineMortar.o \
	$(SRCPATH)Refinement.o 

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

Refinement.o: ../Refinement.C $(DEFINEPARAMS)
