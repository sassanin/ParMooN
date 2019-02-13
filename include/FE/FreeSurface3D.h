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
   
#ifdef __3D__

#ifndef __FREESURFACE3D__
#define __FREESURFACE3D__

#include <Vertex.h>
#include <SquareMatrix3D.h>
#include <FEVectFunct3D.h>

// ========================================================================
// declaration for grid handling
// ========================================================================
void GridBoundCondition(double x, double y, double z, BoundCond &cond);

// ========================================================================
// auxiliary routines
// ========================================================================
void SortVertices(TVertex **Array, int length);
void SortVerticesXYZ(TVertex **Array, int length);
double GetVolume(TCollection *Coll);

// ========================================================================
// calculate all parameters which are needed for calculation
// ========================================================================
void CalculateAllParameters();

// ========================================================================
// calculate all field parameters
// ========================================================================
void CalculateFields();

// ========================================================================
// find normal and tangential vectors for slip d.o.f.
// ========================================================================
void FindVectorsForSlipDOF(TFESpace3D *fespace,
        int &N_FaceDOF, int &N_EdgeDOF,
        int* &FaceDOF, int* &EdgeDOF, 
        double* &FaceVectors, double* &EdgeVectors);

// ========================================================================
// manipulate square matrices due to u.n=0 constraint
// ========================================================================
void ManipulateSquareMatrices(TSquareStructure3D *sqstructure,
        double *a11, double *a12, double *a13,
        double *a21, double *a22, double *a23,
        double *a31, double *a32, double *a33,
        int N_FaceDOF, int N_EdgeDOF,
        int *FaceDOF, int *EdgeDOF,
        double* FaceVectors, double* EdgeVectors);

// ========================================================================
// manipulate square matrices due to u.n=0 constraint
// ========================================================================
void ManipulateMatricesAndRhs(TStructure3D *structure,
        double *b1t, double *b2t, double *b3t,
        double *f1, double *f2, double *f3,
        int N_FaceDOF, int N_EdgeDOF,
        int *FaceDOF, int *EdgeDOF,
        double* FaceVectors, double* EdgeVectors);

// find all joints which form the free surface
void FindFreeSurface(TCollection *Coll,
                     int &N_SurfaceJoints, int* &CellNumbers,
                     int* &JointNumbers);
		     
void FindFreeSurfaceFromJointType(TCollection *Coll, JointType type,
				  int &N_SurfaceJoints, int* &CellNumbers,
				  int* &JointNumbers);

// calculate normal vectors on free surface
void CalculateNormals(TCollection *Coll,
                     int N_SurfaceJoints, int *CellNumbers,
                     int *JointNumbers,
                     TFESpace3D *fespace,
                     double* &n1, double* &n2, double* &n3,
                     double* &len);

// ========================================================================
// ATTENTION !!!
// the current implementation does not check the boundary condition
// it is assumed that all given faces are at the free boundary !!!
// ========================================================================
void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
                 int *CellNumbers, int *JointNumbers,
                 TFEFunction3D *potential, double dt,
                 TSquareMatrix3D *Aii,
                 double *rhs1, double *rhs2, double *rhs3);
		 
void FreeSurfInt(TCollection *Coll, int N_BoundFaces,
		 int *CellNumbers, int *JointNumbers,
		 TAuxParam3D *aux, double dt,
		 TSquareMatrix3D **Aii,
		 double *rhs1, double *rhs2, double *rhs3);
		 
void FreeSurfInt_new(TCollection *Coll, int N_BoundFaces, int *CellNumbers, int *JointNumbers,
		     double dt, TSquareMatrix3D **Aii, double *rhs1, double *rhs2, double *rhs3);
		     
void FreeSurfInt_Sphere(TFESpace3D *fespace, double dt,
			double *rhs1, double *rhs2, double *rhs3);

#endif // __FREESURFACE3D__

#endif // __3D__
