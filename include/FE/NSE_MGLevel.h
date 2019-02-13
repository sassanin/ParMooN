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
// @(#)NSE_MGLevel.h        1.10 07/03/00
//
// Class:       TNSE_MGLevel
// Purpose:     abstract super class for a Stokes/Navier-Stokes system
//
// Author:      Volker John 25.07.2000
//
// History:     25.07.2000 start of implementation
//
// =======================================================================

#ifndef __NSE_MGLEVEL__
#define __NSE_MGLEVEL__

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

#include <Matrix.h>
#ifdef __2D__
   #include <SquareMatrix2D.h>
#endif  
#ifdef __3D__
   #include <SquareMatrix3D.h>
#endif  

class TNSE_MGLevel
{
  protected:
    /** level number */
    int Level;

    /** type number */
    int Type;

#ifdef __2D__
    /** FE space for velocity u */
    TFESpace2D *USpace;

   /** FE space for pressure p */
    TFESpace2D *PSpace;
#endif  

#ifdef __3D__
    /** FE space for velocity u */
    TFESpace3D *USpace;

   /** FE space for pressure p */
    TFESpace3D *PSpace;
#endif  

    /** number of Dirichlet velocity degrees */
    int N_Dirichlet;

    /** number of active velocity degrees */
    int N_Active;

    /** hanging node bound for velocity space */
    int HangingNodeBound;

    /** number of ALL degrees of freedom */
    int N_DOF;

    /** number of all velocity degrees of freedom */
    int N_UDOF;

    /** number of all pressure degrees of freedom */
    int N_PDOF;

    /** velocity space */
    int VelocitySpace;

    /** velocity space */
    int PressureSpace;

    /** array with right-hand sides f1 */
    double *Rhs1;

    /** array with right-hand sides f2 */
    double *Rhs2;

#ifdef __3D__
    /** array with right-hand sides f3 */
    double *Rhs3;
#endif  

    /** array with right-hand sides g */
    double *RhsP;

    /** array with approximate solution u1 */
    double *U1;

    /** array with approximate solution u2 */
    double *U2;

#ifdef __3D__
    /** array with approximate solution u2 */
    double *U3;
#endif  

    /** array with approximate solution p */
    double *P;

    /** number of auxiliary vectors */
    int N_Aux;

    /** array of auxiliary vectors */
    double **Aux;

    /** array for additional data, e.g. ILU decomposition */
    double *Additional;

    /** damping coefficient */
    double alpha;

    /** update coefficient */
    double beta;

    /** collection for Cell-Vanka */
    TCollection *VankaColl;

    /** array for downwind numbering */
    int *downwind;

    /** array for finding local velocity-velocity couplings */
    int *velo_velo_local_coupling;

    /** array for finding local velocity-pressure couplings in gradient */
    int *gradient_local_coupling;

    /** array for finding local velocity-pressure couplings in divergence */
    int *divergence_local_coupling;
    
#ifdef _MPI
    
    TParFECommunicator3D *ParCommU,*ParCommP;

#endif 
  public:
    /** constructor */
    TNSE_MGLevel(int level, double *f1, double *u1,
                 int n_aux, double *al, int VelocitySpace, 
                 int PressureSpace,
                 TCollection *coll);

    /** destructor */
    ~TNSE_MGLevel();

    /** return i-th auxiliary vector */
    double *GetAuxVector(int i);

#ifdef __2D__
    /** return FunctionVectors */
    void GetSolution(double* &u1, double* &u2, double* &p)
      { u1 = U1; u2 = U2; p = P; }

    /** return FunctionVectors */
    void GetRhs(double* &f1, double* &f2, double* &f3)
    { f1 = Rhs1; f2 = Rhs2; f3 = RhsP;}
#endif  
#ifdef __3D__
    /** return FunctionVectors */
    void GetSolution(double* &u1, double* &u2, double* &u3, double* &p)
      { u1 = U1; u2 = U2; u3 = U3; p = P; }

    /** return FunctionVectors */
    void GetRhs(double* &f1, double* &f2, double* &f3, double* &f4)
    { f1 = Rhs1; f2 = Rhs2; f3 = Rhs3; f4 = RhsP;}
#endif  

    /** return FunctionVectors */
    void GetSolutionVector(double* &u1)
      { u1 = U1;}

    /** return FunctionVectors */
    void GetRhsVector(double* &f1)
    { f1 = Rhs1;}

    /** return AuxVectors */
    double **GetAuxVectors()
    { return Aux; }

    /** return number of velocity degrees of freedom */
    int GetN_UDOF()
    { return N_UDOF; }

    /** return number of pressure degrees of freedom */
    int GetN_PDOF()
    { return N_PDOF; }

    /** return velo space */
    int GetVelocitySpace()
    { return VelocitySpace; }

    /** return velo space */
    int GetPressureSpace()
    { return PressureSpace; }

    /** return HangingNodeBound for this level */
    int GetHangingNodeBound()
    { return HangingNodeBound; }

    /** return number of Dirichlet nodes */
    int GetN_Dirichlet()
    { return N_Dirichlet; }

    /** calculate defect */
    virtual void Defect(double *u1, double *f1, double *d1, double &res);

    /** update solution */
    void Update(double *u1, double *v1);

    /** correct Dirichlet and hanging nodes */
    virtual void CorrectNodes(double *u1);

    /** set hanging nodes in order to satisfy coupling condition */
    void SetHangingNodes(double *u1);

    /** correct defect */
    void CorrectDefect(double *v1);

    /** reset vector to zero */
    void Reset(double *v1);

#ifdef __2D__
    /** return FE space for velocity u*/
    TFESpace2D *GetUSpace()
    { return USpace; }

    /** return FE space for pressure p*/
    TFESpace2D *GetPSpace()
    { return PSpace; }
#endif  

#ifdef __3D__
    /** return FE space for velocity u*/
    TFESpace3D *GetUSpace()
    { return USpace; }

    /** return FE space for pressure p*/
    TFESpace3D *GetPSpace()
    { return PSpace; }
#endif  

    TCollection *GetCollection()
    { return VankaColl; }

    void SetCollection(TCollection *coll)
    { VankaColl = coll; }

    /** Vanka smoother */
    virtual void CellVanka(double *u1, double *rhs1, double *aux, 
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

     /** nodal Vanka smoother */
    virtual void NodalVanka(double *u1, double *rhs1, double *aux, 
        int N_Parameters, double *Parameters, int smoother, int N_Levels);

    /** solve exact on this level */
    virtual void SolveExact(double *u1, double *rhs1);

    /** solve exact on this level */
    virtual void SolveExactUMFPACK(double *u1, double *rhs1, int &umfpack_flag);

    /** Braess Sarazin smoother */
    virtual void BraessSarazin(double *u1, double *rhs1, double *aux,
        int N_Parameters, double *Parameters,int N_Levels);

    /** step length control for Vanka */
    virtual double StepLengthControl(double *u1, double *u1old, double *def1,
				     int N_Parameters, double *Parameter);

    /** print all matrices and both right hand sides */
    virtual void PrintAll();

    /** return Type */
    int GetType()
    { return Type; }
    
#ifdef _MPI
    TParFECommunicator3D* GetParCommU(){    
      return ParCommU;      
    }
    TParFECommunicator3D* GetParCommP(){
      return ParCommP;      
    }

    virtual void UpdateHaloRhs(double*, double*); 


    
#endif

};

#endif
