###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
SRCPATH = 
LIBPATH = ../../../lib2D/
INC    = -I../../../include/General -I../../../include/FE \
         -I../../../include/Geometry -I../../../include/QuadFormulas \
         -I../../../include/FE1D -I../../../include/FE2D \
         -I../../../include/AMG \
         -I../../../include/Parallel \
         -I../../../EXT_LIB/UMFPACK/Include \
         -I../../../EXT_LIB/tecplot/include


include ../../../Make.inc

DEFINEPARAMS = ../../../include/General/DefineParams.h

LIBRARY = $(LIBPATH)lib_FE_par.a
   MEMBERS = \
	$(SRCPATH)ADICell.o \
	$(SRCPATH)ADICell1D.o \
	$(SRCPATH)ADISystem.o \
	$(SRCPATH)ADISystem1D.o \
	$(SRCPATH)Assemble2D.o \
	$(SRCPATH)Assemble2D_edge_Oseen.o \
        $(SRCPATH)AuxParam2D.o \
	$(SRCPATH)BaseFunct1D.o \
	$(SRCPATH)BaseFunct2D.o \
	$(SRCPATH)Blas1.o \
	$(SRCPATH)Blas2.o \
        $(SRCPATH)Bulk.o \
	$(SRCPATH)Bulk_2d3d.o \
	$(SRCPATH)CD2DErrorEstimator.o \
	$(SRCPATH)ConvDiff2D.o \
	$(SRCPATH)Convolution.o \
	$(SRCPATH)Darcy2D.o \
	$(SRCPATH)DirectSolver.o \
	$(SRCPATH)DiscreteForm2D.o \
	$(SRCPATH)Drop.o \
	$(SRCPATH)FE1D.o \
	$(SRCPATH)FE2D.o \
	$(SRCPATH)FE2DMapper.o \
	$(SRCPATH)FE2DMapper1Reg.o \
	$(SRCPATH)FEDatabase2D.o \
	$(SRCPATH)FEDesc1D.o \
	$(SRCPATH)FEDesc2D.o \
	$(SRCPATH)FEFunction1D.o \
	$(SRCPATH)FEFunction2D.o \
	$(SRCPATH)FEM_TVD_FCT.o \
	$(SRCPATH)FreeSurface2D.o \
	$(SRCPATH)FESpace.o \
	$(SRCPATH)FESpace1D.o \
	$(SRCPATH)FESpace2D.o \
	$(SRCPATH)FEVectFunct1D.o \
	$(SRCPATH)FEVectFunct2D.o \
	$(SRCPATH)FgmresIte.o \
	$(SRCPATH)FixedPointIte.o \
	$(SRCPATH)HNDesc.o \
	$(SRCPATH)HangingNode.o \
	$(SRCPATH)ItMethod.o \
	$(SRCPATH)JacobiIte.o \
	$(SRCPATH)LinAlg.o \
	$(SRCPATH)LineAffin.o \
	$(SRCPATH)LocalProjection.o \
	$(SRCPATH)MainRoutines2D.o \
	$(SRCPATH)MGComponents2D.o \
	$(SRCPATH)MGLevel2D.o \
	$(SRCPATH)Matrix.o \
	$(SRCPATH)Matrix2D.o \
	$(SRCPATH)MovingNavierStokes.o \
	$(SRCPATH)MultiGrid2D.o \
	$(SRCPATH)MultiGridIte.o \
	$(SRCPATH)MultiGridScaIte.o \
	$(SRCPATH)NSE2DErrorEstimator.o \
	$(SRCPATH)NSE2D_AxialSymm3D_FixPo.o \
	$(SRCPATH)NSE2D_EquOrd_FixPo.o \
	$(SRCPATH)NSE2D_FixPo.o \
	$(SRCPATH)NSE2D_FixPoSkew.o \
	$(SRCPATH)NSE2D_FixPoRot.o \
	$(SRCPATH)NSE2D_Friction_FixPo.o \
	$(SRCPATH)NSE2D_Newton.o \
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
	$(SRCPATH)Output2D.o \
	$(SRCPATH)PDAE_Index2D_2.o \
	$(SRCPATH)QuadAffin.o \
	$(SRCPATH)QuadBilinear.o \
	$(SRCPATH)QuadIsoparametric.o \
	$(SRCPATH)RefTrans1D.o \
	$(SRCPATH)RefTrans2D.o \
	$(SRCPATH)RFB.o \
	$(SRCPATH)Solver.o \
	$(SRCPATH)SquareMatrix.o \
	$(SRCPATH)SquareMatrix1D.o \
	$(SRCPATH)SquareMatrix2D.o \
	$(SRCPATH)SquareMatrixNSE2D.o \
	$(SRCPATH)StructureNSE2D.o \
	$(SRCPATH)SquareStructure.o \
	$(SRCPATH)SquareStructure1D.o \
	$(SRCPATH)SquareStructure2D.o \
	$(SRCPATH)SSMUM.o \
	$(SRCPATH)Structure.o \
	$(SRCPATH)Structure2D.o \
	$(SRCPATH)Superconvergence2D.o \
	$(SRCPATH)TCD2D.o \
	$(SRCPATH)TimeDiscRout.o \
	$(SRCPATH)TNSE2D_FixPo.o \
	$(SRCPATH)TNSE2D_FixPoRot.o \
	$(SRCPATH)TNSE2D_FixPo_SSMUM.o \
	$(SRCPATH)TriaAffin.o \
	$(SRCPATH)TriaIsoparametric.o \
	$(SRCPATH)Upwind.o\
	$(SRCPATH)VMS.o
###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache

Assemble.o: ../Assemble.C $(DEFINEPARAMS)
Assemble2D_edge_convdiv.o: ../Assemble2D_edge_convdiv.C $(DEFINEPARAMS)
Assemble2D_edge_Oseen.o: ../Assemble2D_edge_Oseen.C $(DEFINEPARAMS)
FESpace2D.o: ../FESpace2D.C $(DEFINEPARAMS)
SquareStructure.o: ../SquareStructure.C $(DEFINEPARAMS)
Structure.o: ../Structure.C $(DEFINEPARAMS)

