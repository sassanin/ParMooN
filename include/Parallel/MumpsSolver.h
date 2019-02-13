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
// @(#)MumpsSolver.h
//
// Class:      TMumpsSolver
// Purpose:    Solve equation system by  MPI based MUMPS parallel direct solvers
//
// Author:     Sashikumaar Ganesan (30.09.09)
//
// History:    Start of implementation 30.09.09 (Sashikumaar Ganesan)
//
// =======================================================================
#if defined(_MPI) || defined(_SMPI)
 extern "C"
 {
  #  include "dmumps_c.h"
 }

#ifndef __MUMPSSOLVER__
#define __MUMPSSOLVER__

 #include <SquareMatrix.h>
 #include <Matrix.h>
 #include <SquareMatrix3D.h>
 #include <Matrix3D.h>
 #include <ParFECommunicator3D.h>

/** general class for all parallel direct solvers */

class TMumpsSolver
{
  protected:

  /** MPI_Comm for which the fespace communications are needed */
  MPI_Comm Comm;

  /**  internal pointers of MUMPS solver*/
  DMUMPS_STRUC_C id;
 
  /** global rhs */
  double *MumpsRhs;
 
  /** */
  bool FactorFlag;
  
  
  public:
   /** constructor */
   TMumpsSolver(int N_Eqns, int M_dist_Nz, int *M_dist_Irn, int *M_dist_Jcn, int N_Rhs);

   void FactorizeAndSolve(double *Mat_loc, double *rhs);

   void Solve(double *Mat_loc, double *rhs);
   
   void Clean();
  
    /** destructor */
    ~TMumpsSolver();
};
#endif


#endif
