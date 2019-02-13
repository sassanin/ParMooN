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
   
#include <Convolution.h>
#include <Database.h>
#include <TNSE2D_Routines.h>
#include <stdlib.h>


/* ======================================================================
   depending on the global orientation of the normals at inner edges the sign 
   of the basis functions need to be inversed. This is only important for 
   Raviart-Thomas finite elements. (called from DarcyRaviartThomas(...))
   
   This information is stored in the FE-Descriptor. However the FE-Descriptor 
   is not accessible in the assembling routine DarcyRaviartThomas. Therefore we
   use global parameters in TDatabase. These need to be set for every cell.
   
   Here it is assumed that Raviart-Thomas elements of order 0,1,2, or 3 on
   triangles or quadrilaterals are used
*/
int GetSignOfThisDOF(int N_DOF, int DOF)
{
  switch (N_DOF)
  {
  case 3:// Raviart-Thomas zeroth order, triangles
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[DOF];
    break;
    // note that RT2 on quadrilaterals and RT3 on triangles have 24 basis functions
  case 4:// Raviart-Thomas zeroth order, quadrilaterals
    return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[DOF];
    break;
  case 6: // BDM first order, triangles
    // degree of freedom on an edge, no inner degree of freedom
    return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
    break;
  case 8:
    // Raviart-Thomas first order, triangles
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
      if(DOF<6) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%2)/2];
      else // inner degree of freedom
        return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1011) //BDM first order, quadrilaterals
    {
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    }
    break;
  case 12:
    // Raviart-Thomas first order, quadrilaterals
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1001)
    {
    if(DOF<8) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%2)/2];
    else // inner degree of freedom
      return 1;
    }
    else if (TDatabase::ParamDB->VELOCITY_SPACE == 1012) //BDM second order, triangles
    {
      if(DOF<9) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    break;
  case 14://BDM second order, quadrilaterals
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 15:// Raviart-Thomas second order, triangles
    if(DOF<9) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%3)/3];
    else // inner degree of freedom
      return 1;
    break;
  case 20:// BDM third order. triangles
    if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 22://BDM third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else // inner degree of freedom
      return 1;
    break;
  case 24:
    if(TDatabase::ParamDB->VELOCITY_SPACE == 1002)
    {
      // Raviart-Thomas second order, quadrilaterals
      if(DOF<12) // degree of freedom on an edge
        return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%3)/3];
      else // inner degree of freedom
        return 1;
    }
    else if(TDatabase::ParamDB->VELOCITY_SPACE == 1003)
    {
      // Raviart-Thomas third order, triangles
      if(DOF<12) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_TRIA[(DOF-DOF%4)/4];
      else // inner degree of freedom
        return 1;
    }
    else
    {
      Error("VELOCITY_SPACE has to be set to either 1002 or 1003\n");
      exit(0);
    }
    break;
  case 40:// Raviart-Thomas third order, quadrilaterals
    if(DOF<16) // degree of freedom on an edge
      return TDatabase::ParamDB->NORMAL_ORIENTATION_QUAD[(DOF-DOF%4)/4];
    else
      return 1;
    break;
  default:
       OutPut("N_DOF" << N_DOF << endl);
       OutPut("WARNING: Unknown Raviart-Thomas or BDM element !" << endl);
       return 1;
       break;
  }
}

// ======================================================================
// (DarcyType 1)
// Standard Galerkin with Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
// elements
// ======================================================================
void BilinearAssembleDarcyGalerkin(double Mult, double *coeff, double *param,
                                   double hK, double **OrigValues,
                                   int *N_BaseFuncts, double ***LocMatrices,
                                   double **LocRhs)
{
  double val;
  double ansatz, ansatz_x_00, ansatz_y_00;
  double test00, test_x_00, test_y_00, test_x_10, test_y_01;
  double test_div;
  
  // ( A  B1 )   ( 0 2 )
  // ( B2 C  )   ( 3 1 )
  
  double **MatrixA = LocMatrices[0];
  double **MatrixC = LocMatrices[1];
  double **MatrixB1 = LocMatrices[2];
  double **MatrixB2 = LocMatrices[3];
  
  double *Rhs0 = LocRhs[0];
  double *Rhs1 = LocRhs[1];
  
  int N_U = N_BaseFuncts[0];
  int N_P = N_BaseFuncts[1];

  double *Orig0 = OrigValues[0];   // u
  double *Orig1 = OrigValues[1];   // p
  double *Orig2 = OrigValues[2];   // u_x
  double *Orig3 = OrigValues[3];   // u_y

  double c0 = coeff[0];            // sigma
  double c1 = coeff[1];          // f1
  double c2 = coeff[2];          // f2
  double c3 = coeff[3];          // g(x,y)

  int *newsign = new int[N_U];
  for(int i=0;i<N_U;i++)
  {
    // here check whether signs should be inverted
    newsign[i] = GetSignOfThisDOF(N_U,i); // in Darcy2D.C
  }
  // A, B1, B2
  for(int i=0;i<N_U;i++)
  {
    // A:
    test_x_00 = newsign[i]*Orig0[i];
    test_y_00 = newsign[i]*Orig0[N_U+i];
    
    Rhs0[i] += Mult*(c1*test_x_00 + c2*test_y_00);
    
    for(int j=0;j<N_U;j++)
    {
      ansatz_x_00 = newsign[j]*Orig0[j];
      ansatz_y_00 = newsign[j]*Orig0[N_U+j];

      // A: u_x v_x + u_y v_y
      val  = c0*(test_x_00*ansatz_x_00 + test_y_00*ansatz_y_00);
      MatrixA[i][j] += Mult * val;
    }
    // B1, B2:
    test_x_10 = newsign[i]*Orig2[i];
    test_y_01 = newsign[i]*Orig3[N_U+i];
    test_div = test_x_10 + test_y_01;
    for(int j=0;j<N_P;j++)
    {
      ansatz = Orig1[j];
      val = Mult*test_div*ansatz;
      // (p div v)
      MatrixB1[i][j] -= val;
      // (q, div u)
      MatrixB2[j][i] -= val; // slow (consider moving this into another loop)
    }
  }
  
  for(int i=0;i<N_P;i++)
  {
    test00 = Orig1[i];
    // assemble rhs: div u = g
    // rhs: -(g,q)
    // c3 = g(x,y)
    Rhs1[i] -= Mult*test00*c3;
    // C:
    /*
    for(int j=0;j<N_P;j++)
    {
      ansatz = Orig1[j];
      val = Mult*test*ansatz;
      // (p,q)
      MatrixC[i][j] += val;
    }
    */
  }
  delete [] newsign;
}
