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
// %W% %G%
//
// Class:       TMGLevel3D
// Purpose:     store all data for one level in a multi grid method in 3d
//
// Author:      Gunar Matthies 26.06.2000
//
// History:     26.06.2000 start of implementation
//
// =======================================================================

#ifndef __MGLEVEL3D__
#define __MGLEVEL3D__

#include <SquareMatrix3D.h>
#ifdef _MPI   
   #ifdef __3D__
    #include <ParFECommunicator3D.h>
   #else
    #include <ParFECommunicator2D.h>
   #endif
#endif 

// #include <omp.h>
class TMGLevel3D
{
  protected:
    /** level number */
    int Level;

    /** FE space */
    TFESpace3D *FESpace;

    /** permutation vector */
    int *Permutation;

    /** number of active nodes */
    int N_Active;

    /** upper bound for hanging node number */
    int HangingNodeBound;

    /** number of Dirichlet nodes */
    int N_Dirichlet;

    /** number of all degrees of freedom */
    int N_DOF;

    /** used matrix */
    TSquareMatrix3D *A;

    /** structure of used matrix */
    TSquareStructure3D *MatrixStructure;

    /** row pointer for matrix */
    int *RowPtr;

    /** column number vector */
    int *KCol;

    /** matrix entries */
    double *Entries;

    /** array with right-hand sides */
    double *Rhs;

    /** array with approximate solution */
    double *X;

    /** number of auxiliary vectors */
    int N_Aux;

    /** array of auxiliary vectors */
    double **Aux;

    /** array for additional data, e.g. ILU decomposition */
    double *Additional;

    /** generate ILU decomposition */
    void ILUDecomposition();
    
#ifdef _MPI
     /** number of all degrees of freedom in own cells*/
     int OwnN_DOF;  
     
     double *Temp_arr;
    
     TParFECommunicator3D *ParComm; 
     TParFEMapper3D *ParMapper;
     
     /** Reorder of sol array */
     int *Reorder,N_Master,N_Int,N_Dept;
     
     int *Reorder_M,*Reorder_D1,*Reorder_D2,*Reorder_I;
     
     int N_InterfaceM, N_Dept1, N_Dept2;// N_Dept3;
     
#endif
     
#ifdef _HYBRID     
     /** Coloring variables */
     int N_CMaster, N_CDept1, N_CDept2, N_CInt, *ptrCMaster, *ptrCDept1, *ptrCDept2, *ptrCInt;
#endif

  public:

    
    /** constructor */
    TMGLevel3D(int level, TSquareMatrix3D *A,
             double *rhs, double *sol, int n_aux,
             int *permutation);
    
#ifdef _MPI
/** constructor for parallel */
    TMGLevel3D(int level, TSquareMatrix3D *a, double *rhs, double *sol, 
                       TParFECommunicator3D *parComm, TParFEMapper3D *parMapper, int n_aux,
                       int *permutation);
    double *GetTemp_arr()
    { return Temp_arr; }
#endif  

    /** destructor */
    ~TMGLevel3D();

    /** return i-th auxiliary vector */
    double *GetAuxVector(int i);

    /** return FunctionVectors */
    double *GetSolution()
    { return X; }

    /** return Rhs */
    double *GetRhs()
    { return Rhs; }
    
    

    /** return AuxVectors */
    double **GetAuxVectors()
    { return Aux; }

    /** return system matrix */
    TSquareMatrix3D *GetMatrix()
    { return A; }

    /** return number of degrees of freedom */
    int GetN_DOF()
    { return N_DOF; }

    /** get HangingNodeBound */
    int GetHangingNodeBound()
    { return HangingNodeBound; }

    /** get number of Dirichlet nodes */
    int GetN_Dirichlet()
    { return N_Dirichlet; }

    /** calculate defect */
    void Defect(double *sol, double *f, double *d, double &res);

    /** update solution */
    void Update(double *sol, double *upd);

    /** correct Dirichlet and hanging nodes */
    void CorrectNodes(double *vect);

    /** correct defect */
    void CorrectDefect(double *vect);

    /** reset vector to zero */
    void Reset(double *vect);

    /** return FE space */
    TFESpace3D *GetFESpace()
    { return FESpace; }

    /** smoother */
    void ILU(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void SSOR(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** smoother */
    void Jacobi(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);
     

    
    /** smoother */
    void Block2x2(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);

    /** solve exact on this level */
    void SolveExact(double *u1, double *rhs1);

    /** step length control */
    double StepLengthControl(double *u, 
                         double *uold, 
                         double *def,
                         int N_Parameters, 
                         double *Parameters);
    
    
        
#ifdef _MPI       

    TParFECommunicator3D *GetParComm()
      { return ParComm; }
      
       void SOR_Re(double *sol, double *f, double *aux,
        int N_Parameters, double *Parameters);   
    
#endif

#ifdef _HYBRID
        void SOR_Re_Color(double *sol, double *f, double *aux,
			  int N_Parameters, double *Parameters, bool firstTime, bool lastTime);
	void SOR_Re_Color(double *sol, double *f, double *aux, int N_Parameters, double *Parameters,int smooth);
	
	void SOR_Re_Color_Coarse(double *sol, double *f, double *aux,
				 int N_Parameters, double *Parameters);
#endif
       
};

#endif
