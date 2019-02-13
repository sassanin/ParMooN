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
// @(#)NSE_MultiGrid.h        1.4 02/08/00
//
// Class:       TNSE_MultiGrid
// Purpose:     store all data for a multi grid method for 
//              Stokes or Navier-Stokes problems
//
// Author:      Volker John 25.07.2000
//
// History:     25.07.2000 start of implementation
//
// =======================================================================

#ifndef __NSE_MULTIGRID__
#define __NSE_MULTIGRID__

#include <NSE_MGLevel.h>

class TNSE_MultiGrid
{
  protected:
    /** number of levels */
    int N_Levels;

    /** number of problems */
    int N_Problems;

    /** number of parameters */
    int N_Parameters;

    /** array of double parameters */
    double *Parameters;

    /** array of multi grid levels */
    TNSE_MGLevel *MultiGridLevels[MAXN_LEVELS];

#ifdef __2D__
    /** array of FE spaces */
    TFESpace2D *USpaces[MAXN_LEVELS];

    /** array of FE spaces */
    TFESpace2D *PSpaces[MAXN_LEVELS];
#endif  
#ifdef __3D__
    /** array of FE spaces */
    TFESpace3D *USpaces[MAXN_LEVELS];

    /** array of FE spaces */
    TFESpace3D *PSpaces[MAXN_LEVELS];
#endif  

    /** array of u1 vectors on each level */
    double **U1Vectors[MAXN_LEVELS];

    /** array of u2 vectors on each level */
    double **U2Vectors[MAXN_LEVELS];

#ifdef __3D__
    /** array of u3 vectors on each level */
    double **U3Vectors[MAXN_LEVELS];
#endif  

    /** array of p vectors on each level */
    double **PVectors[MAXN_LEVELS];

    /** right-hand side vectors */
    double **Rhs1Vectors[MAXN_LEVELS];

    /** right-hand side vectors */
    double **Rhs2Vectors[MAXN_LEVELS];

#ifdef __3D__
    /** right-hand side vectors */
    double **Rhs3Vectors[MAXN_LEVELS];
#endif  

    /** right-hand side vectors */
    double **RhsPVectors[MAXN_LEVELS];

    /** auxiliary vectors */
    double **AuxVectors[MAXN_LEVELS];

    /** number of recursions */
    int mg_recursions[MAXN_LEVELS];

  public:
    /** constructor */
    TNSE_MultiGrid(int n_problems, int n_parameters, double *parameters);

    /** return number of multi grid levels */
    int GetN_Levels()
    { return N_Levels; }

    /** add new level as finest */
    void AddLevel(TNSE_MGLevel *MGLevel);

    /** replace level i by given MGLevel and return old level */
    TNSE_MGLevel *ReplaceLevel(int i, TNSE_MGLevel *MGLevel);

    /** return i-th level as TMGLevel object */
    TNSE_MGLevel *GetLevel(int i)
    { return MultiGridLevels[i]; }

    /** restrict u1, u2 from finest grid to all coarser grids */
    void RestrictToAllGrids();

    /** set correct values for Dirichlet nodes on grid i */
    void SetDirichletNodes(int i);

    /** return residual on grid i */
    double GetResidual(int i);

    /** cycle on level i */
    void Cycle(int i, double &res);

    /** set recursion for multigrid */ 
    void SetRecursion(int levels);

    /** get parameter */
    double GetParam(int i);

    /** get parameter */
    void SetParam(int i, double a);

};

#endif
