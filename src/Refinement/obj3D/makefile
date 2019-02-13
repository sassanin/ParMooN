###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib3D/
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
	$(SRCPATH)Refinement.o \
	$(SRCPATH)RefHexaRegDesc.o \
	$(SRCPATH)RefTetraRegDesc.o \
	$(SRCPATH)RefTetraReg0Desc.o \
	$(SRCPATH)RefTetraReg1Desc.o \
	$(SRCPATH)RefTetraReg2Desc.o \
	$(SRCPATH)RefTetraBis0Desc.o \
	$(SRCPATH)RefTetraBis1Desc.o \
	$(SRCPATH)RefTetraBis2Desc.o \
	$(SRCPATH)RefTetraBis3Desc.o \
	$(SRCPATH)RefTetraBis4Desc.o \
	$(SRCPATH)RefTetraBis5Desc.o \
	$(SRCPATH)RefTetraBis01Desc.o \
	$(SRCPATH)RefTetraBis02Desc.o \
	$(SRCPATH)RefTetraBis03Desc.o \
	$(SRCPATH)RefTetraBis04Desc.o \
	$(SRCPATH)RefTetraBis05Desc.o \
	$(SRCPATH)RefTetraBis10Desc.o \
	$(SRCPATH)RefTetraBis12Desc.o \
	$(SRCPATH)RefTetraBis13Desc.o \
	$(SRCPATH)RefTetraBis14Desc.o \
	$(SRCPATH)RefTetraBis15Desc.o \
	$(SRCPATH)RefTetraBis20Desc.o \
	$(SRCPATH)RefTetraBis21Desc.o \
	$(SRCPATH)RefTetraBis23Desc.o \
	$(SRCPATH)RefTetraBis24Desc.o \
	$(SRCPATH)RefTetraBis25Desc.o \
	$(SRCPATH)RefTetraBis30Desc.o \
	$(SRCPATH)RefTetraBis32Desc.o \
	$(SRCPATH)RefTetraBis34Desc.o \
	$(SRCPATH)RefTetraBis35Desc.o \
	$(SRCPATH)RefTetraBis40Desc.o \
	$(SRCPATH)RefTetraBis41Desc.o \
	$(SRCPATH)RefTetraBis43Desc.o \
	$(SRCPATH)RefTetraBis45Desc.o \
	$(SRCPATH)RefTetraBis51Desc.o \
	$(SRCPATH)RefTetraBis52Desc.o \
	$(SRCPATH)RefTetraBis53Desc.o \
	$(SRCPATH)RefTetraBis54Desc.o \
	$(SRCPATH)RefTetraQuad0Desc.o \
	$(SRCPATH)RefTetraQuad1Desc.o \
	$(SRCPATH)RefTetraQuad2Desc.o \
	$(SRCPATH)RefTetraQuad3Desc.o 

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

Refinement.o: ../Refinement.C $(DEFINEPARAMS)
