/** ==========================================================================
#    This file is part of the finite element software ParMooN.
# 
#    ParMooN (cmg.cds.iisc.ac.in/parmoon) is a free finite element software  
#    developed by the research groups of Prof. Sashikumaar Ganesan (IISc, Bangalore),
#    Prof. Volker John (WIAS Berlin) and Prof. Gunar Matthies (TU-Dresden):
#
#    ParMooN is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    ParMooN is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with ParMooN.  If not, see <http://www.gnu.org/licenses/>.
#
#    If your company is selling a software using ParMooN, please consider 
#    the option to obtain a commercial license for a fee. Please send 
#    corresponding requests to sashi@iisc.ac.in

# =========================================================================*/ 
   
#ifndef __ALLCLASSES__
#define __ALLCLASSES__

class TADICell;
class TADICell1D;
class TADISystem;
class TADISystem1D;
class TAuxParam2D;
class TAuxParam3D;
class TAssemble;
class TAssemble3D;
class TBaseCell;
class TBaseFunct1D;
class TBaseFunct2D;
class TBaseFunct3D;
class TBdCircle;
class TBDEdge3D;
class TBdLine;
class TBdPlane;
class TBdPolygon;
class TBdSpline;
class TBrick;
class TBoundComp;
class TBoundComp2D;
class TBoundComp3D;
class TBoundEdge;
class TBoundFace;
class TBoundPart;
class TBoundPoint;
class TCD2DErrorEstimator;
class TCollection;
class TDatabase;
class TDiscreteForm2D;
class TDiscreteForm3D;
class TDomain;
class TEdge;
class TFE1D;
class TFE2D;
class TFE3D;
class TFEDatabase2D;
class TFEDatabase3D;
class TFEDesc1D;
class TFEDesc2D;
class TFEDesc3D;
class TFEFunction2D;
class TFEFunction3D;
class TFE2DMapper;
class TFE2DMapper1Reg;
class TFE3DMapper;
class TFESpace;
class TFESpace1D;
class TFESpace2D;
class TFESpace3D;
class TFEVectFunct1D;
class TFEVectFunct2D;
class TGridCell;
class THNDesc;
class THangingNode;
class THexaAffin;
class THexaTrilinear;
class THexahedron;
class TInnerEdge;
class TInterfaceJoint;
class TInterfaceJoint3D;
class TIsoBoundEdge;
class TIsoEdge3D;
class TIsoInterfaceJoint;
class TIsoInterfaceJoint3D;
class TIsoJointEqN;
class TIt_Between;
class TIt_EQ;
class TIt_EQLevel;
class TIt_Finest;
class TIt_LE;
class TIt_LELevel;
class TIt_Mortar;
class TIt_OCAF;
class TIt_Search;
class TIterator;
class TJoint;
class TJointCollection;
class TJointEqN;
class TJoint_2to1;
class TLine;
class TMGLevel2D;
class TMGLevel3D;
class TMacroCell;
class TMapper;
class TMatrix;
class TMatrix2D;
class TMatrix3D;
class TMortarBaseJoint;
class TMortarJoint;
class TMultiGrid2D;
class TMultiGrid3D;
class TNS2DErrorEstimator;
class TNodalFunctional1D;
class TNodalFunctional2D;
class TNodalFunctional3D;
class TOutput2D;
class TParallelogram;
class TParFECommunicator2D;
class TParVector;
class TParVectorNSE;
class TParVectorNSE3D;
class TPeriodicJoint;
class TQuadAffin;
class TQuadBilinear;
class TQuadFormula;
class TQuadFormula1D;
class TQuadFormula2D;
class TQuadFormula3D;
class TQuadFormulaHexa;
class TQuadFormulaQuad;
class TQuadFormulaTetra;
class TQuadFormulaTria;
class TQuadIsoparametric;
class TQuadrangle;
class TRectangle;
class TRefDesc;
class TRefHexaRegDesc;
class TRefLineDesc;
class TRefMortar0Desc;
class TRefMortar1Desc;
class TRefMortarLineDesc;
class TRefNoRef;
class TRefQuad1Conf0Desc;
class TRefQuad1Conf1Desc;
class TRefQuad1Conf2Desc;
class TRefQuad1Conf3Desc;
class TRefQuad2Conf0Desc;
class TRefQuad2Conf1Desc;
class TRefQuad2Conf2Desc;
class TRefQuad2Conf3Desc;
class TRefQuadBis0Desc;
class TRefQuadBis1Desc;
class TRefQuadRegDesc;
class TRefQuadToTri0Desc;
class TRefQuadToTri1Desc;
class TRefTetraRegDesc;
class TRefTetraReg0Desc;
class TRefTetraReg1Desc;
class TRefTetraReg2Desc;
class TRefTetraRegDesc;
class TRefTrans1D;
class TRefTrans2D;
class TRefTrans3D;
class TRefTriBis0Desc;
class TRefTriBis1Desc;
class TRefTriBis2Desc;
class TRefTriRegDesc;
class TShapeDesc;
class TSquareMatrix;
class TSquareMatrix1D;
class TSquareMatrix2D;
class TSquareMatrix3D;
class TSquareStructure;
class TSquareStructure1D;
class TSquareStructure2D;
class TSquareStructure3D;
class TSubDomainHaloJoint;
class TSubDomainEdge3D;
class TSubDomainJoint;
class TNSE2DMGLevel;
class TNSE2DMGLevel1;
class TNSE2DMGLevel2;
class TNSE2DMGLevel3;
class TNSE2DMGLevel4;
class TNSE2DMultiGrid;
class TStructure;
class TStructure2D;
class TStructure3D;
class TTetraAffin;
class TTetrahedron;
class TTriaAffin;
class TTriaIsoparametric;
class TLineAffin;
class TTriangle;
class TVertex;
class TVirtCell;

#endif
