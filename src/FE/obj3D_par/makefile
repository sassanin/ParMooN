###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
LIBPATH = ../../../lib3D/
SRCPATH = 

INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE2D -I../../../include/FE1D \
         -I../../../include/FE3D -I../../../include/AMG \
         -I../../../include/PBE -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../include/Parallel \
         -I../../../EXT_LIB/tecplot/include -I../../../EXT_LIB/MUMPS/include

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_FE_par.a
   MEMBERS = \
	$(SRCPATH)ADICell.o \
	$(SRCPATH)ADICell1D.o \
	$(SRCPATH)ADISystem.o \
	$(SRCPATH)ADISystem1D.o \
	$(SRCPATH)ALE_3D.o \
        $(SRCPATH)Assemble2D.o \
        $(SRCPATH)Assemble3D.o \
        $(SRCPATH)Aux2D3D.o \
        $(SRCPATH)AuxParam2D.o \
        $(SRCPATH)AuxParam3D.o \
        $(SRCPATH)BaseFunct1D.o \
        $(SRCPATH)BaseFunct2D.o \
        $(SRCPATH)BaseFunct3D.o \
        $(SRCPATH)Bcgs.o \
        $(SRCPATH)Blas1.o \
        $(SRCPATH)Blas2.o \
        $(SRCPATH)Bulk.o \
        $(SRCPATH)Bulk_3d4d.o \
        $(SRCPATH)Cg.o \
        $(SRCPATH)ConvDiff.o \
        $(SRCPATH)ConvDiff3D.o \
        $(SRCPATH)Convolution.o \
        $(SRCPATH)DirectSolver.o \
        $(SRCPATH)DiscreteForm2D.o \
        $(SRCPATH)DiscreteForm3D.o \
        $(SRCPATH)FE1D.o \
        $(SRCPATH)FE2D.o \
        $(SRCPATH)FE3D.o \
        $(SRCPATH)FE2DMapper.o \
        $(SRCPATH)FE2DMapper1Reg.o \
        $(SRCPATH)FE3DMapper.o \
        $(SRCPATH)FE3DMapper1Reg.o \
        $(SRCPATH)FEDatabase2D.o \
        $(SRCPATH)FEDatabase3D.o \
        $(SRCPATH)FEDesc1D.o \
        $(SRCPATH)FEDesc2D.o \
        $(SRCPATH)FEDesc3D.o \
        $(SRCPATH)FEFunction1D.o \
        $(SRCPATH)FEFunction2D.o \
        $(SRCPATH)FEFunction3D.o \
        $(SRCPATH)FEM_TVD_FCT.o \
        $(SRCPATH)FESpace.o \
        $(SRCPATH)FESpace1D.o \
        $(SRCPATH)FESpace2D.o \
        $(SRCPATH)FESpace3D.o \
        $(SRCPATH)FEVectFunct3D.o \
        $(SRCPATH)FgmresIte.o \
        $(SRCPATH)FixedPointIte.o \
        $(SRCPATH)FreeSurface3D.o \
        $(SRCPATH)HNDesc.o \
        $(SRCPATH)HangingNode.o \
        $(SRCPATH)HexaAffin.o \
        $(SRCPATH)HexaIsoparametric.o \
        $(SRCPATH)HexaTrilinear.o \
        $(SRCPATH)ItMethod.o \
        $(SRCPATH)JacobiIte.o \
        $(SRCPATH)LinAlg.o \
        $(SRCPATH)LineAffin.o \
        $(SRCPATH)LocalProjection.o \
        $(SRCPATH)MGComponents2D.o \
        $(SRCPATH)MGComponents3D.o \
        $(SRCPATH)MGLevel3D.o \
        $(SRCPATH)Matrix.o \
        $(SRCPATH)Matrix3D.o \
        $(SRCPATH)MultiGrid3D.o \
        $(SRCPATH)MultiGridIte.o \
        $(SRCPATH)MultiGridScaIte.o \
        $(SRCPATH)NSE3D_FixPo.o \
        $(SRCPATH)NSE3D_FixPoSkew.o \
        $(SRCPATH)NSE3D_Friction_FixPo.o \
        $(SRCPATH)NSE3D_Newton.o \
        $(SRCPATH)NSE_MGLevel.o \
        $(SRCPATH)NSE_MGLevel1.o \
        $(SRCPATH)NSE_MGLevel2.o \
        $(SRCPATH)NSE_MGLevel3.o \
        $(SRCPATH)NSE_MGLevel4.o \
	$(SRCPATH)NSE_MGLevel5.o \
        $(SRCPATH)NSE_MGLevel14.o \
	$(SRCPATH)NSE_MultiGrid.o \
	$(SRCPATH)NodalFunctional1D.o \
        $(SRCPATH)NodalFunctional2D.o \
        $(SRCPATH)NodalFunctional3D.o \
        $(SRCPATH)Output2D.o \
        $(SRCPATH)Output3D.o \
	$(SRCPATH)PardisoSolver.o \
        $(SRCPATH)QuadAffin.o \
        $(SRCPATH)QuadIsoparametric.o \
        $(SRCPATH)QuadBilinear.o \
        $(SRCPATH)Queue.o \
        $(SRCPATH)RationalLES.o \
        $(SRCPATH)RefTrans3D.o \
        $(SRCPATH)RFB.o \
        $(SRCPATH)Solver.o \
        $(SRCPATH)SquareMatrix.o \
        $(SRCPATH)SquareMatrix1D.o \
        $(SRCPATH)SquareMatrix2D.o \
        $(SRCPATH)SquareMatrix3D.o \
        $(SRCPATH)SquareStructure.o \
        $(SRCPATH)SquareStructure1D.o \
        $(SRCPATH)SquareStructure2D.o \
        $(SRCPATH)SquareStructure3D.o \
        $(SRCPATH)SSORIte.o \
        $(SRCPATH)Structure.o \
        $(SRCPATH)Structure3D.o \
        $(SRCPATH)Superconvergence3D.o \
        $(SRCPATH)TCD3D.o \
	$(SRCPATH)TimeDiscRout.o \
        $(SRCPATH)TNSE3D_FixPo.o \
        $(SRCPATH)TNSE3D_Newton.o \
        $(SRCPATH)TetraAffin.o \
        $(SRCPATH)TetraIsoparametric.o \
        $(SRCPATH)TriaAffin.o \
        $(SRCPATH)TriaIsoparametric.o \
        $(SRCPATH)Upwind3D.o \
        $(SRCPATH)Urea_3d4d.o \
        $(SRCPATH)VMS.o \
	$(SRCPATH)Windtunnel_3d4d.o 

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

