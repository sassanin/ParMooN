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
// @(#)ADISystem1D.h        4.1 13.11.09
// 
// Class:       TADISystem1D
// Purpose:     class for  ADISystem1D

//
// Author:      Sashikumaar Ganesan (13.11.09)
//
// History:     start of implementation 13.11.09 (Sashikumaar Ganesan)
//
// =======================================================================

#ifndef __ADISYSTEM1D__
#define __ADISYSTEM1D__

#include <ADISystem.h>
#include <FESpace1D.h>
#include <FEFunction1D.h>

class TADISystem1D : public TADISystem
{
  protected:

    /** x value of physical space */
    double X;

    /** y value of physical space */
    double Y;

    /** z value of physical space */
    double Z;   

    /** Fe space of configuration space */
    TFESpace1D *FESpace1D_Intl;

    /** fe function, needed ofr interpolation */
    TFEFunction1D *FEFunction_Intl;

    /** mass matrices for all QuadPt will not change */
    TSquareMatrix1D *M_Intl;

    /** stiffness matrix for all QuadPts*/
    TSquareMatrix1D *A_Intl;

    /** supg mass matrices for all QuadPt will not change */
    TSquareMatrix1D *S_Intl;

    /** supg stiffness matrix for all QuadPts*/
    TSquareMatrix1D *K_Intl;


    /** Initial condition*/
     DoubleFunctND *Initial;

    /** Boundary values*/
     DoubleFunct2D  *BDValue_LMin;
     DoubleFunct2D  *BDVal_LMax;

  private:
   int ConstructAllInfo();

  public:
    /** constructor */
    TADISystem1D(TFEFunction1D *FEFunction_intl, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
                           TSquareMatrix1D *S_internal, TSquareMatrix1D *K_internal,
                           double x, double y, DoubleFunctND *initial, 
                           DoubleFunct2D *bDValue_LMin, DoubleFunct2D *bDVal_LMax,
                           double *Sol,  double *OldSol, double *b, 
                           double *Defect, double *Intlposl, DoubleFunctND *growthAndNuc);

    /** constructor */
    TADISystem1D(TFEFunction1D *FEFunction_intl, TSquareMatrix1D *M_internal, TSquareMatrix1D *A_internal,
                 double x, double y, double z, double *Sol, double *OldSol, double *b, double *Defect,
                 DoubleFunctND *growthfunct);
  
    //void Interpolate(double *Sol, DoubleFunct3D *Exact);
    void SetDirichletBc(BoundCond cond_Lmin, BoundCond cond_Lmax, double BDValue_Lmin, double BDValue_Lmax);

    void Interpolate(double *Sol, DoubleFunctND *Exact);
    
    void SolveAllLevels(double *MatValues_orig, double *SolFromX, double *AggrRhs, double &C, double C_Sat, 
                        double &T, CoeffFctND *BilinearCoeffs, 
                        double tau, BoundCond cond_Lmin, BoundCond cond_Lmax,  DoubleFunct3D *Exact );

    void AssembleInitRhs(double C, double T, CoeffFctND *BilinearCoeffs, 
                         BoundCond cond_Lmin, BoundCond cond_Lmax);
 
    void AssembleARhs_FD(double Conv, CoeffFctND *Bilinear);
	
    void AssembleARhs(double Conv, CoeffFctND *Bilinear);
    
    void AssembleARhs_SUPG(double Conv, CoeffFctND *Bilinear);

    void AssembleARhs_DG(double Conv, double Bnuc, CoeffFctND *Bilinear, BoundCond cond_Lmin, BoundCond cond_Lmax);
    
    double GetWeightedF();

   int GetN_InternalLevels()
    {
     return N_V;
    }

    double GetQ3Max(double *currsol);

    /** destrcutor */
    ~TADISystem1D();

};

#endif
