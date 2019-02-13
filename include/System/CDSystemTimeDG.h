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
   
/** ************************************************************************ 
*
* @class     TCDSystemTimeDG
* @brief     stores the information of a 2D scalar dG in time discretization
* @author    Sashikumaar Ganesan, 
* @date      23.12.15
* @History    
 ************************************************************************  */


#ifndef __CDSYSTEMTIMEDG__
#define __CDSYSTEMTIMEDG__

#include <SquareMatrix2D.h>
#include <string.h>

/**class for 2D scalar system  dG in time discretization */
class TCDSystemTimeDG
{
  protected:
    /** total number of dof in dG solution */
    int N_Eqn;   
    
    /** block mat size in dG system */
    int N_U, N_UActive;   
    
    /** dG solution array  */
    double *Sol_dG, *rhs_dG;
    
    /** number of block matrices  */
    int N_RowBlockMat, N_ColBlockMat;  
    
    /** no-nonzers, rowptr and column index of dG system matrix */
    int Sys_N_Entries, *Sys_rowptr, *Sys_colindex;    

    /** mass and stiffness matrices  */
    TSquareMatrix2D *Mat_M, *Mat_A, *Sys_Mat;

    /** mass and stiffness matrices in ALE approach  */
    TSquareMatrix2D *Mat_M_Qp1, *Mat_A_Qp1; // dG(1)
    double *Rhs_Qp1;
    bool CONSERVATIVEALE, QpMatricsAdded;
     
    /** system mat structure */
    TSquareStructure2D *Sys_structure;
    
    /** rhs fespace for assemblling RHS*/
    TFESpace2D *RhsSpace;
    
    /** BilinearCoeffs for assemblling RHS*/
    CoeffFct2D *BilinearCoeffs;
    
  public:
    /** constructor */
     TCDSystemTimeDG(int n_RowBlockMat, int n_ColBlockMat, TSquareMatrix2D *mat_M, TSquareMatrix2D *mat_A);

    /** destrcutor */
    ~TCDSystemTimeDG();
    
    void AddBilinear(CoeffFct2D *bilinearCoeffs)
          {BilinearCoeffs =  bilinearCoeffs;}

    /** add additinal matrics on quadpts in ALE approach */
    void AddQp1Matrices(TSquareMatrix2D *mat_M_Qp1,TSquareMatrix2D *mat_A_Qp1, double *rhs_Qp1)
    {
     Mat_M_Qp1 = mat_M_Qp1; 
     Mat_A_Qp1 = mat_A_Qp1;
     Rhs_Qp1 = rhs_Qp1;
     QpMatricsAdded = TRUE;
    }
          
    /** assemble the system matrix */
    virtual void AssembleSysMat(double *Mu_old, double *Rhs);
    
    /** assemble the system matrix */   
    virtual void AssembleALESysMat_Qp1(double *Mu_old, double *Rhs);
    
    /** solve the system matrix */
    void AssembleRhs(int N_Rhs, double *T, double *Rhs);

    /** all sys mat entries will be reset to zero */
    void ResetSysMat()
    {  memset(Sys_Mat->GetEntries(), 0, Sys_N_Entries*SizeOfDouble);  }
    
    /** solve the system matrix */
    void SovedGSystem();
 
    /** solve dG system and return the Sol at end t^n system matrix */
    virtual void SoveTimedG(double *Sol);   
    
    /** set the ALE form CONSERVATIVEALE T/F */
    void SetALEForm(bool conservative)
    {  CONSERVATIVEALE = conservative; }
    
};

#endif
