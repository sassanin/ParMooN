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
// Class:       TFE3D
// Purpose:     store all information for one finite element class
//
// Author:      Gunar Matthies  19.11.99
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
# endif

#include <FE3D.h>
#include <FEDatabase3D.h>
#include <Database.h>

/** constructor */
TFE3D::TFE3D()
{
}

/** constructor with data */
TFE3D::TFE3D(BaseFunct3D basefunct_id, NodalFunctional3D nodalfunctional_id,
         RefTrans3D reftransid, FEDesc3D fedesc_id, int n_info)
{
  BaseFunct_ID = basefunct_id;
  BaseFunct = TFEDatabase3D::GetBaseFunct3D(BaseFunct_ID);

  NodalFunctional_ID = nodalfunctional_id;
  NodalFunctional  = TFEDatabase3D::GetNodalFunctional3D(NodalFunctional_ID);

  RefTransID = reftransid;

  FEDesc_ID = fedesc_id;
  FEDesc = TFEDatabase3D::GetFEDesc3D(FEDesc_ID);

  N_Info = n_info;
  N_DOF = BaseFunct->GetDimension();

  Size = N_Info + N_DOF;
}

/** check N[i](b[j]) = delta[ij] */
void TFE3D::CheckNFandBF()
{
  int i,j,k,l, N_Points;
  double *xi, *eta, *zeta;
  
  NodalFunctional->GetPointsForAll(N_Points, xi, eta, zeta);
  
  int n_basis_functions = BaseFunct->GetDimension();
  int baseVectDim = BaseFunct->GetBaseVectDim();
  double PointValues[N_Points*baseVectDim];
  double FunctionalValues[n_basis_functions];
  double AllPointValues[N_Points][n_basis_functions*baseVectDim];

  for(k=0;k<N_Points;k++)
    BaseFunct->GetDerivatives(D000, xi[k], eta[k], zeta[k],
                              AllPointValues[k]);
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==TDatabase::ParamDB->Par_P0)
#endif
  if(TDatabase::ParamDB->SC_VERBOSE > 0)
  {
   cout << "CheckNFandBF: " << "BaseFunct_ID: " << BaseFunct_ID << " ";
   cout << "NodalFunctional_ID: " << NodalFunctional_ID << endl;
  }
  for(k=0;k<N_DOF;k++)
  {
    for(l=0;l<N_Points;l++)
      for(i = 0; i< baseVectDim; i++)
        PointValues[l+N_Points*i] = AllPointValues[l][k+n_basis_functions*i];

    NodalFunctional->GetAllFunctionals(NULL, NULL, PointValues,
                          FunctionalValues);

    for(i=0;i<N_DOF;i++)
    {
      if(fabs(FunctionalValues[i])<1e-10)
        FunctionalValues[i] = 0;
      if( i == k )
        if( fabs(FunctionalValues[i]-1) > 1e-8 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
      if( i != k )
        if( fabs(FunctionalValues[i]-0) > 1e-8 )
          cout << "BF: " << k << " NF: " << i << " " << FunctionalValues[i] << endl;
    }
  }

#ifdef _MPI
  if(rank==TDatabase::ParamDB->Par_P0)
#endif 
  if(TDatabase::ParamDB->SC_VERBOSE > 0)
  cout << endl;
}
