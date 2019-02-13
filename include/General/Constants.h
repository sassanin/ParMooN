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
   
// =======================================================================
// @(#)Constants.h        1.15 04/13/00
// 
// Purpose:     constains important constants which are used by several 
//              other classes
//
// Author:      Volker Behns  09.07.97
//
// =======================================================================

#ifndef __CONSTANTS__
#define __CONSTANTS__

#ifdef NULL
  #undef NULL
#endif

#define NULL 0

#include <math.h>

#ifdef __2D__
#define GEO_DIM 2
#else
#define GEO_DIM 3
#endif


#define Pi M_PI
#define ln2 log((double) 2.0)

#ifndef __LINUX__
  #define Abs(x) abs(x)
#else
  #define Abs(x) fabs(x)
#endif

#define SizeOfPointer  ((long) 4)
#define SizeOfInt      ((long) 4)
#define SizeOfFloat    ((long) 4)
#define SizeOfDouble   ((long) 8)

#ifdef __COMPAQ__
  #undef SizeOfPointer
  #define SizeOfPointer ((long) 8)
#endif

#ifdef __HPaCC64__
  #undef SizeOfPointer
  #define SizeOfPointer ((long) 8)
#endif

#ifdef FALSE
#undef FALSE
#endif
#ifdef TRUE
#undef TRUE
#endif

#ifdef _MPI
  #define BYADD 0
  #define BYMASTER 1
  #define BYOWN 2
#endif

enum boolean {FALSE, TRUE};

#define N_BOUNDCOND 10
#define Magnetic 1
#define conductivity 3.46e6
#define Liquid_Density 6360
#define Kinematic_Viscosity 3.4e-7
#define a_duct .04
#define d_duct .01
#define Ha a_duct*Magnetic*sqrt(conductivity/(Liquid_Density*Kinematic_Viscosity))
#define THERMAL_CONDUCTIVITY 16.5
#define SPECIFIC_HEAT 152.6
#define Prandtl_Number Liquid_Density*Kinematic_Viscosity*SPECIFIC_HEAT/(THERMAL_CONDUCTIVITY)
#define MHD_K 2*Ha*(d_duct*d_duct)/(a_duct*a_duct) 

enum BoundCond { DIRICHLET, NEUMANN, ROBIN, SLIP, FREESURF, 
                 SLIP_FRICTION_PENETRATION_RESISTANCE, 
                 INTERFACE, SUBDOMAIN_INTERFACE, SUBDOMAIN_HALOBOUND, DIRICHLET_WEAK};

enum JointType {Joint, JointEqN, MortarBaseJoint, MortarJoint,
                BoundaryPoint, BoundaryEdge, BoundaryFace,
                InterfaceJoint, PeriodicJoint, IsoInterfaceJoint,
                IsoJointEqN, IsoBoundEdge, IsoBoundFace,
                Joint_2to1, 
                InterfaceJoint3D, IsoInterfaceJoint3D,
                SubDomainJoint, SubDomainHaloJoint, InnerInterfaceJoint,
                InnerEdge, IsoEdge3D, BDEdge3D};

typedef void DoubleFunct1D(double, double *);
typedef void DoubleFunct2D(double, double, double *);
typedef void DoubleFunct3D(double, double, double, double *);
typedef void DoubleFunctND(int, double *, double *);
typedef void DoubleFunctVect(double *, double *);
typedef int IntFunct2D(double, double);
typedef double DoubleFunct2Param(double, double);

typedef void BoundCondFunct3D(int, double, double, double, BoundCond &);
typedef void BoundValueFunct3D(int, double, double, double, double &);
typedef void BoundCondFunct2D(int, double, BoundCond &);
typedef void BoundValueFunct2D(int, double, double &);

typedef void ErrorMethod2D(int, double *, double *, 
                           double *, double *, double,
                           double **, double **,
                           double **, double *);

typedef void ErrorMethod3D(int, double *, double *, double *,
                           double *, double *, double,
                           double **, double **,
                           double **, double *);

typedef DoubleFunctVect ParamFct;

typedef void CoeffFct2D(int n_points, double *X, double *Y,
                        double **param, double **coeffs);

typedef void CoeffFct3D(int n_points, double *X, double *Y, double *Z,
                        double **param, double **coeffs);
                        
typedef void CoeffFctND(int n_points, int N_Dim, double **Coords,
                        double **param, double **coeffs);  

typedef void AssembleFct2D(double, double *, double, double **, 
                           int *, double ***, double **);

typedef void AssembleFctParam2D(double, double *, double *,
                                double, double **, 
                                int *, double ***, double **);

typedef void AssembleFct3D(double, double *, double, double **, 
                           int *, double ***, double **);

typedef void AssembleFctParam3D(double, double *, double *,
                                double, double **, 
                                int *, double ***, double **);

class TBaseCell;
typedef void ManipulateFct2D(int, double **, double **, TBaseCell *);

typedef void ManipulateFct3D(int, double **, double **, TBaseCell *);

class TCollection;
typedef void EvalAllNF(TCollection *, TBaseCell *, double *, double *);
typedef void EvalJointNF(TCollection *, TBaseCell *, int, double *, double *);

class TFESpace2D;
typedef void CheckWrongNeumannNodesFunct2D(TCollection *, TFESpace2D *,
					   int &, int* &,
					   int* &, 
					   double* &);

class TSquareMatrix;
class TMatrix;
typedef void MatVecProc(TSquareMatrix **, TMatrix **, double *, double *);
typedef void DefectProc(TSquareMatrix **, TMatrix **, double *, double *,
                        double *);

class TVertex;
class TIsoBoundEdge;
typedef void ModifyMeshCoords(double , double , double &, double &, double );
typedef void ModifyBoundCoords(int *, TVertex **, TIsoBoundEdge **,  double *, double);

typedef void ModifyMeshCoords_3D(double , double ,double , double &, double &,double &, double );
typedef void ModifyBoundCoords_3D(int, TVertex **, int, int *, TBaseCell **, double);
class TFEVectFunct3D;
typedef void MoveBound_3D(TFEVectFunct3D *, TFEVectFunct3D *, double *, int, TVertex **, int, int *, TBaseCell **, double);

#ifdef _MPI
class TParVectorNSE3D;
typedef void ParDefectProc(TSquareMatrix **, TMatrix **, TParVectorNSE3D  *, TParVectorNSE3D *,
                        TParVectorNSE3D *);
#endif

typedef int TypeBoundSwitchFunct2D(int, double );

class TFEFunction2D;
class TFEVectFunct2D;
typedef void EvaluateSolutionFunct2D(TFEFunction2D **, TFEVectFunct2D **, 
                                     double *, int *);


#define GALERKIN         1
#define SDFEM            2
#define SUPG             2
#define UPWIND           3
#define SMAGORINSKY      4
#define SMAGORINSKY_EXPL 5
#define CLASSICAL_LES    6
#define GL00_CONVOLUTION 7
#define GL00_AUX_PROBLEM 8
#define VMS_PROJECTION   9
#define VMS_PROJECTION_EXPL 10
#define VMS_RFB_EXPL     11
#define VMS_PROJECTION_SD     12
#define VMS_RFB_EXPL_COUPLED     13
#define LERAY_ALPHA_EXPL   16
#define LERAY_ALPHA_IMPL   17
#define NSE_RFB          21
#define SDFEM_DIVDIV     5
#define LOCAL_PROJECTION 14
#define LOCAL_PROJECTION_2_LEVEL 15
#define CIP              4
#define DG               5
#define GLS              6
#define FD               0
#define ENO_3           13
#define WENO_5          14
#define HEATLINE        100

#define OSEEN_PROBLEM    13

#define HMM86            0 
#define TP86_1           1
#define TP86_2           2 
#define JSW87            3 
#define GdC88            4 
#define dCG91            5
#define dCA03            6 
#define AS97             7
#define C93              8 
#define KLR02_1          9
#define KLR02_2          10
#define KLR02_3          11
#define KLR02_4          12
#define J90              13
#define BE02_1           14 
#define BE02_2           15 
#define BH04             16 
#define BE05_1           17 
#define BE05_2           18 
#define LP96             19
#define CS99             20
#define MH_Kno06         21
#define BE02_3           22 
//#define ALG_FLUX_CORR    23 
#define Y_Z_BETA         24 
#define JSW87_1          25 
#define FEM_FCT          50
#define FEM_FCT_LIN      51
#define FEM_TVD          52
#define GENERAL_SOLD    200 
 
#define SD_POWER_H  1.0
#define SD_FACTOR_H 1.0

#define BULK_FWE_FDM_UPWIND    0
#define BULK_BWE_FDM_UPWIND    1
#define BULK_BWE_FEM_SUPG      2
#define BULK_FEM_FCT       3

#define WINDTUNNEL_LAYER_NUMBER_X_CONST 4
#define WINDTUNNEL_DIM_Y_CONST 47 
// be careful number of grid points plus 1 right dimension 46
#define WINDTUNNEL_DIM_Z_CONST 20
// be careful number of grid points plus 1 right dimension 19
#define WINDTUNNEL_DIM_R_CONST 47
//#define WINDTUNNEL_DIM_R_CONST 51
// be careful number of grid points plus 1 right dimension 45

#define WINDTUNNEL_FWE_FDM_UPWIND    0
#define WINDTUNNEL_BWE_FDM_UPWIND    1
#define WINDTUNNEL_FEM_FCT           2

#define UREA_FWE_FDM_UPWIND    0
#define UREA_BWE_FDM_UPWIND    1
#define UREA_FEM_FCT           2

#define PB_RKV_ENO             4

#define ABS(i)        (((i)<0) ? (-(i)) : (i))

#ifdef MIN
  #undef MIN
#endif

#ifdef MAX
  #undef MAX
#endif

#ifndef __HPaCC__
#ifndef __HPaCC64__
#define MIN(x,y)    (((x)<(y)) ? (x) : (y))
#define MAX(x,y)    (((x)>(y)) ? (x) : (y))
#endif
#endif

#define POW2(i)       (1<<(i))
#define POW(a,b)      ((b==1)?(a):((b==2)?((a)*(a)):((b==0)?(1):(pow(a,b)))))
#define ABSDIFF(a,b)  (fabs((a)-(b)))

#define SIGN(a)       ((a>0)?(1):((a<0)?(-1):(0)))

double GetTime();
int GetMemory();

#define AMG_SOLVE 0
#define GMG 1
#define DIRECT 2

// maybe we should define more types here (like CD, TCD, TNSE, TSTOKES,...)
#define NSE    0
#define STOKES 1
#define OSEEN  2

#define SCALAR 100

#define MAXN_LEVELS 100

#ifdef _MPI
  #define BEGIN_SEQ \
    MPI_Barrier(comm); \
    for(int _i=0;_i<size;++_i) \
    { \
      if ( rank == _i ) \
      { 
      
  #define END_SEQ \
      } \
      MPI_Barrier(comm); \
    }
#endif


#endif
