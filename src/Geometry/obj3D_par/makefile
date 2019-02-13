###################################
#
#  Makefile for MooNMD
#  author: Volker Behns
#
###################################

RELPATH = ../../..
LIBPATH = ../../../lib3D/
SRCPATH = 
INC    = -I../../../include/Geometry -I../../../include/FE \
         -I../../../include/General -I../../../include/QuadFormulas \
         -I../../../EXT_LIB/tetgen

include ../../../Make.inc

LIBRARY = $(LIBPATH)lib_geo_par.a
   MEMBERS = \
	$(SRCPATH)BaseCell.o \
	$(SRCPATH)BdCircle.o \
	$(SRCPATH)BdCylinder.o \
	$(SRCPATH)BDEdge3D.o \
	$(SRCPATH)BdLine.o \
	$(SRCPATH)BdPlane.o \
	$(SRCPATH)BdPolygon.o \
	$(SRCPATH)BdSphere.o \
	$(SRCPATH)BdSpline.o \
	$(SRCPATH)BdNonUniformSpline.o \
	$(SRCPATH)BdNoPRM.o \
	$(SRCPATH)BdWall.o \
	$(SRCPATH)Brick.o \
	$(SRCPATH)BoundComp.o \
	$(SRCPATH)BoundComp2D.o \
	$(SRCPATH)BoundComp3D.o \
	$(SRCPATH)BoundEdge.o \
	$(SRCPATH)BoundFace.o \
	$(SRCPATH)BoundPart.o \
	$(SRCPATH)Collection.o \
	$(SRCPATH)Domain.o \
	$(SRCPATH)Edge.o \
	$(SRCPATH)GridCell.o \
	$(SRCPATH)Hexahedron.o \
	$(SRCPATH)InnerEdge.o \
	$(SRCPATH)InterfaceJoint.o \
	$(SRCPATH)InterfaceJoint3D.o \
	$(SRCPATH)IsoBoundEdge.o \
	$(SRCPATH)IsoBoundFace.o \
	$(SRCPATH)IsoEdge3D.o \
	$(SRCPATH)IsoInterfaceJoint.o \
	$(SRCPATH)IsoInterfaceJoint3D.o \
	$(SRCPATH)IsoJointEqN.o \
	$(SRCPATH)It_Between.o \
	$(SRCPATH)It_EQ.o \
	$(SRCPATH)It_EQLevel.o \
	$(SRCPATH)It_Finest.o \
	$(SRCPATH)It_LE.o \
	$(SRCPATH)It_LELevel.o \
	$(SRCPATH)It_Mortar.o \
	$(SRCPATH)It_OCAF.o \
	$(SRCPATH)It_Search.o \
	$(SRCPATH)Iterator.o \
	$(SRCPATH)Joint.o \
	$(SRCPATH)JointCollection.o \
	$(SRCPATH)JointEqN.o \
	$(SRCPATH)Line.o \
	$(SRCPATH)MacroCell.o \
	$(SRCPATH)Mapper.o \
	$(SRCPATH)MortarBaseJoint.o \
	$(SRCPATH)MortarJoint.o \
	$(SRCPATH)Parallelogram.o \
	$(SRCPATH)PeriodicJoint.o \
	$(SRCPATH)Quadrangle.o \
	$(SRCPATH)ReadGeo.o \
	$(SRCPATH)Rectangle.o \
	$(SRCPATH)RefDesc.o \
	$(SRCPATH)ShapeDesc.o \
	$(SRCPATH)SubDomainEdge3D.o \
	$(SRCPATH)SubDomainJoint.o \
	$(SRCPATH)SubDomainHaloJoint.o \
	$(SRCPATH)Tests.o \
	$(SRCPATH)Tetrahedron.o \
	$(SRCPATH)TetGenMeshLoader.o \
	$(SRCPATH)Triangle.o \
	$(SRCPATH)Vertex.o

###############

$(LIBRARY) : $(MEMBERS)
	@echo "ar "$(LIBRARY)
	@$(AR) $@ $(MEMBERS)

clean:
	@echo removing object files
	@$(RM) $(MEMBERS) $(LIBRARY) ii_files SunWS_cache
