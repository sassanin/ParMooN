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
* @brief     source file for TAssemble
* @author    Sashikumaar Ganesan
* @date      30.05.15
* @History 
 ************************************************************************  */

#include <Assemble.h>



TAssemble::TAssemble(int n_fespaces, int n_sqmatrices, int n_matrices, int n_rhs, double **rhs)
{
  N_SqMatrices = n_sqmatrices;
  N_Matrices = n_matrices;
  N_FeSpaces = n_fespaces;    
  N_Rhs = n_rhs;
  N_AllMatrices = n_sqmatrices+n_matrices;
  
  
  if(N_SqMatrices)
   {  
    GlobalNumbers = new int* [n_sqmatrices];
    BeginIndex = new int* [n_sqmatrices];
    HangingEntries = new double* [n_sqmatrices];

   } //if(N_SqMatric 
   
  if(n_matrices)
  {
    TestGlobalNumbers = new int* [n_matrices];
    AnsatzGlobalNumbers = new int* [n_matrices];
    TestBeginIndex = new int* [n_matrices];
    AnsatzBeginIndex = new int* [n_matrices];
  } // endif n_matrices
  
  if(n_rhs)
   {   
    Rhs = new double*[n_rhs];     
    RhsBeginIndex = new int* [n_rhs];
    RhsGlobalNumbers = new int* [n_rhs];
    LocRhs = new double* [n_rhs]; 
   }

 if(N_AllMatrices)
  LocMatrices = new double** [N_AllMatrices]; 
  
}

TAssemble::~TAssemble()
{
  
  
}
