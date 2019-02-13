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
// @(#)Upwind.C        1.6 11/24/99
//
// Purpose:     upwind stabilization
//
// Authors:     Volker Behns / Gunar Matthies  18.10.99
//              Volker John, bug fixes 18. - 22.10.99
//
// =======================================================================

#include <Database.h>
#include <DiscreteForm3D.h>
#include <FEFunction3D.h>
#include <FE3D.h>
#include <SquareMatrix3D.h>
#include <Upwind3D.h>
#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

void GetUpwindValue(double UPWIND_ORDER,double t,double *phi)
{
  double s;

  if (UPWIND_ORDER==1) // Samarskij upwind  
    {
      if (t<0)
        *phi = 1/(2*(1-t));
      else
        *phi = (1+2*t)/(2*(1+t));
    }
  else 
    {
      if (UPWIND_ORDER==0) // sharp upwind 
        *phi = (SIGN(t)+1)/2;
      else
        {
          s = fabs(t);
          *phi = 0.5+0.5*SIGN(t)*POW(s/(s+1), UPWIND_ORDER);
        }
    }
  return;
}


// =======================================================================
// do upwind assembling for first order nonconforming elements
// =======================================================================
void UpwindForNavierStokes3D(TSquareMatrix3D *sqmatrix, TFEFunction3D *u1,
                             TFEFunction3D *u2, TFEFunction3D *u3)
{
  double RE=TDatabase::ParamDB->RE_NR;
  double UPWIND_ORDER=TDatabase::ParamDB->UPWIND_ORDER;
  double UPWIND_FLUX_DAMP=TDatabase::ParamDB->UPWIND_FLUX_DAMP;
  int  UPW_APPL = TDatabase::ParamDB->UPWIND_APPLICATION;
  TBaseCell *cell;
  TCollection *coll;
  int i,j,k,l,m,n,m1, N_Cells, N_Edges, N_Faces, N_Rows, N_Vertices;
  int n_mat, m_mat;
  double *u1vect, *u2vect, *u3vect;
  int *GlobalNumbers, *BeginIndex, *DOF;
  TFESpace3D *fespace;
  double matrix[8][8];
  double val1, val2, val, val3, flux, RE2=RE*0.5;
  double xcoords[8], ycoords[8], zcoords[8];
  double v1, v2, v3, v4;
  double sx, sy, sz, t, phi, s;
  double nx, ny, nz, Area;
  double x,y,z;
  int *ColInd, *RowPtr;
  double *Entries;
  int ActiveBound, DirichletBound, end;
  TShapeDesc *desc;
  FE3D CurrentElement;
  const int *TmpEV, *TmpEF;
  int MaxLen;
  // values in the barycentres of the dual planes
  // weights with respect to the values in the faces
/*  double FctValNTetra[6][6] = { {  0.666666666666666666667, 
                                   -0.333333333333333333333,
                                    0.666666666666666666667 },
                                 {  0.666666666666666666667, 
                                    0.666666666666666666667,
                                   -0.333333333333333333333 },
                                 { -0.333333333333333333333,
                                    0.666666666666666666667, 
                                    0.666666666666666666667 } };*/

  // get pointers to columns, rows and entries of matrix A
  ColInd = sqmatrix->GetKCol();
  RowPtr = sqmatrix->GetRowPtr();
  Entries = sqmatrix->GetEntries();

  // get finite element space for velocity
  fespace = u1->GetFESpace3D();
  
  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get number of free dof
  ActiveBound = fespace->GetActiveBound();
  // get start of dirichlet nodes in dof array
  DirichletBound = fespace->GetHangingBound();
  
  // get vectors of current velocity
  u1vect = u1->GetValues();
  u2vect = u2->GetValues();
  u3vect = u3->GetValues();
  
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // loop over all cells for computing the upwind stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];
    // # of vertices
    N_Vertices = cell->GetN_Vertices();
    // # of faces
    N_Faces = cell->GetN_Faces();
    // get coordinates of the vertices
    for(k=0;k<N_Vertices;k++)
      cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k], zcoords[k]);
    
    // calculate the coordinates of the center of gravity
    sx = sy = sz = 0.0;
    for(k=0;k<N_Vertices;k++)
    {
      sx += xcoords[k];
      sy += ycoords[k];
      sz += zcoords[k];
    }     
    sx /= N_Vertices;
    sy /= N_Vertices;
    sz /= N_Vertices;

    // reset local matrix
    memset(matrix, 0, SizeOfDouble*64);
    CurrentElement = fespace->GetFE3D(i,cell);
    switch(CurrentElement)
    {
      // ********** nonconforming discretization *******************
      case N_P1_3D_T_A:
      case N_Q1_3D_H_A:
      case N_Q1_3D_H_M:                
	// # of edges
	N_Edges = cell->GetN_Edges();
        
	desc = cell->GetShapeDesc();
	desc->GetEdgeFace(TmpEF, MaxLen);
	desc->GetEdgeVertex(TmpEV);

        // for all edges
        for(n=0;n<N_Edges;n++)
        {
	  // vertices if edge n
	  l=TmpEV[2*n];
	  m=TmpEV[2*n+1];
	  // face on edge n 
	  k=TmpEF[MaxLen*n];
	  j=TmpEF[MaxLen*n+1];
          
	  // cout << "l,m,j,k:" << l << m << j << k << endl;
	  //cout << "Koordinaten[l]:" << xcoords[l] << " " << ycoords[l] << " " << zcoords[l] << endl;
	  //cout << "Koordinaten[m]:" << xcoords[m] << " " << ycoords[m] << " " << zcoords[m] << endl;
	  //cout << "Schwerpunkt:" << sx << " " << sy << " " << sz << endl;
          
          // normal to face of dual domain in direction of face k
          nx = 0.5*((ycoords[l]-sy)*(zcoords[m]-sz)-(zcoords[l]-sz)*(ycoords[m]-sy));
          ny = -0.5*((xcoords[l]-sx)*(zcoords[m]-sz)-(zcoords[l]-sz)*(xcoords[m]-sx));
          nz = 0.5*((xcoords[l]-sx)*(ycoords[m]-sy)-(ycoords[l]-sy)*(xcoords[m]-sx));

	  Area=sqrt(nx*nx+ny*ny+nz*nz);

	  nx /= Area;
	  ny /= Area;
	  nz /= Area;

	  //cout << "nx,ny,nz:" << nx << " " << ny << " " << nz << endl;
          
          // compute current velocity in center of edge of dual
          // boundary          
          switch(N_Edges)
          {
            // tetrahedron
            case 6: 
              // dofs of the faces at the current edge
              m=DOF[k];
              m1=DOF[j];
              val1 = 0.75*(u1vect[m]+u1vect[m1]);
              val2 = 0.75*(u2vect[m]+u2vect[m1]);
              val3 = 0.75*(u3vect[m]+u3vect[m1]);
              // dofs of the faces not at the current edge
              for (l=0;l<4;l++)
              {
                if (l==j) continue; 
                if (l==k) continue; 
                m = DOF[l];
                val1 -= 0.25*u1vect[m];
                val2 -= 0.25*u2vect[m];
                val3 -= 0.25*u3vect[m];
              }
              break;
              // hexahedron       
            case 12: 
              // dofs of the faces at the current edge
              m=DOF[k];
              m1=DOF[j];
              val1 = 0.5*(u1vect[m]+u1vect[m1]);
              val2 = 0.5*(u2vect[m]+u2vect[m1]);
              val3 = 0.5*(u3vect[m]+u3vect[m1]);
              //cout << i << " " << m << "  " << m1 << endl;
              //cout << u1vect[m] << " " << u1vect[m1] << " " << val1 << endl;
              //cout << u2vect[m] << " " << u2vect[m1] << " " << val2 << endl;
              //OutPut("val3 " << m << " " << u3vect[m] << " "  << m1 << " " << u3vect[m1] << " " << val3 << endl);
              break;
          } // endswitch
          // compute flux across face of dual boundary from face j to face k
          flux = val1*nx + val2*ny + val3*nz;
          // parameter for upwind function
          //OutPut(CurrentElement << " flux " <<  flux << " " << val1 <<  " " << val2 << " " << val3);
          //OutPut(" " << nx << " "  << ny << " " << nz << endl);
          t = UPWIND_FLUX_DAMP*flux*RE2;
          // compute upwind value
          GetUpwindValue(UPWIND_ORDER,t,&phi);
	  flux *= Area;
          // compute matrix entries
          switch(UPW_APPL)
          {
          case 0:
            //OutPut("mat "<< flux << " " << phi << endl); 
	    matrix[k][k] += flux*phi;
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(1-phi);
            break;
          case 1:  
            matrix[k][k] += flux*(-0.5+phi);
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(0.5-phi);
            break;
          default:
            OutPut("Wrong Upwind Application !!! Default is used !!!"<< endl);
            matrix[k][k] += flux*phi;
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(1-phi);
            break;
          }
        } // endfor k
        break;
     // ********** conforming discretization *******************	
      default:
        OutPut("Upwind works only for nonconforming elements!" << endl);
        exit(-1);
    } // endswitch

    for(m=0;m<N_Faces;m++)
    {
      l = DOF[m];
      //cout << "DOF: " << l << endl;
      if (l < ActiveBound)
      {
        // node l is inner or Neumann node
        end = RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Faces;k++)
          {
            if(DOF[k] == ColInd[n])
            {
              //cout << m << "   " << k << endl << n << endl;
              //OutPut("E " << Entries[n] << " " << matrix[m][k] << endl);
              Entries[n] += matrix[m][k];
            } // endif
          } // endfor k
        } // endfor n
      } // endif l
      else
      {
        if (l < DirichletBound)
        {
          // hanging node
          /*

          // implement later
          l -= ActiveBound;
          end = HangingRowPtr[l+1];
          for(n=HangingRowPtr[l];n<end;n++)
          {
            for(k=0;k<N_;k++)
            {
              if(DOF[k] == HangingColInd[n])
              {
                CurrentHangingEntries[n] += MatrixRow[k];
              } // endif
            } // endfor k
          } // endfor n 
          */
        }
        else
        {
          // Dirichlet node => do nothing
        }
      }
    } // endfor m
  } // endfor i
}


void UpwindForConvDiff(TSquareMatrix3D *sqmatrix, double *RHS,
                       TFESpace3D *fespace, TDiscreteForm3D *DiscreteForm)
{
  static double RE=TDatabase::ParamDB->RE_NR;
  double UPWIND_ORDER=TDatabase::ParamDB->UPWIND_ORDER;
  double UPWIND_FLUX_DAMP=TDatabase::ParamDB->UPWIND_FLUX_DAMP;
  int  UPW_APPL = TDatabase::ParamDB->UPWIND_APPLICATION;
  TBaseCell *cell;
  TCollection *coll;
  TShapeDesc *desc;
  int i,j,k,l,m,n, N_Edges, N_Cells, N_Vertices, N_Rows, N_Faces;
  int n_mat, m_mat;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double matrix[8][8];
  double Area, val1, val2, val, flux;
  double xcoords[8], ycoords[8], zcoords[8];
  double v1, v2, v3, v4;
  double sx, sy, sz, t, phi, s;
  double nx, ny, nz, x, y, z;
  int *ColInd, *RowPtr;
  double *Entries;
  int ActiveBound, DirichletBound, end;
  double *coeff = new double[20], *param = new double[1], RE2=RE*0.5;
  const int *TmpEV, *TmpEF;
  int MaxLen;
  FE3D CurrentElement;

  // get pointers to columns, rows and entries of matrix A
  ColInd = sqmatrix->GetKCol();
  RowPtr = sqmatrix->GetRowPtr();
  Entries = sqmatrix->GetEntries();

  // get arrays with the numbering of the dof
  GlobalNumbers = fespace->GetGlobalNumbers();
  BeginIndex = fespace->GetBeginIndex();

  // get number of free dof
  ActiveBound = fespace->GetActiveBound();
  // get start of dirichlet nodes in dof array
  DirichletBound = fespace->GetHangingBound();
  
  // get collection and number of cells
  coll = fespace->GetCollection();
  N_Cells = coll->GetN_Cells();

  // loop over all cells for computing the upwind stabilization
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    // pointer to global indices of dof connected with this cell
    DOF = GlobalNumbers + BeginIndex[i];
    // # of vertices
    N_Vertices = cell->GetN_Vertices();
    // # of faces
    N_Faces = cell->GetN_Faces();
    // get coordinates of the vertices
    for(k=0;k<N_Vertices;k++)
      cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k], zcoords[k]);

    // calculate the coordinates of the center of gravity
    sx = sy = sz = 0.0;
    for(k=0;k<N_Vertices;k++)
    {
      sx += xcoords[k];
      sy += ycoords[k];
      sz += zcoords[k];
    }     
    sx /= N_Vertices;
    sy /= N_Vertices;
    sz /= N_Vertices;

    // reset local matrix
    memset(matrix, 0, SizeOfDouble*64);

    CurrentElement = fespace->GetFE3D(i,cell);
    switch(CurrentElement)
    {
      // ********** nonconforming discretization *******************
      case N_P1_3D_T_A:
      case N_Q1_3D_H_A:
      case N_Q1_3D_H_M:                
	// # of edges
	N_Edges = cell->GetN_Edges();

	desc = cell->GetShapeDesc();
	desc->GetEdgeFace(TmpEF, MaxLen);
	desc->GetEdgeVertex(TmpEV);

        // for all edges
        for(n=0;n<N_Edges;n++)
        {
	  // vertices if edge n
	  l=TmpEV[2*n];
	  m=TmpEV[2*n+1];
	  // face on edge n 
	  k=TmpEF[MaxLen*n];
	  j=TmpEF[MaxLen*n+1];

	  // cout << "l,m,j,k:" << l << m << j << k << endl;
	  //cout << "Koordinaten[l]:" << xcoords[l] << " " << ycoords[l] << " " << zcoords[l] << endl;
	  //cout << "Koordinaten[m]:" << xcoords[m] << " " << ycoords[m] << " " << zcoords[m] << endl;
	  //cout << "Schwerpunkt:" << sx << " " << sy << " " << sz << endl;

          // normal to face of dual domain in direction of face k
          nx = 0.5*((ycoords[l]-sy)*(zcoords[m]-sz)-(zcoords[l]-sz)*(ycoords[m]-sy));
          ny = -0.5*((xcoords[l]-sx)*(zcoords[m]-sz)-(zcoords[l]-sz)*(xcoords[m]-sx));
          nz = 0.5*((xcoords[l]-sx)*(ycoords[m]-sy)-(ycoords[l]-sy)*(xcoords[m]-sx));

	  Area=sqrt(nx*nx+ny*ny+nz*nz);

	  nx /= Area;
	  ny /= Area;
	  nz /= Area;

	  //cout << "nx,ny,nz:" << nx << " " << ny << " " << nz << endl;

          // midpoints of face of dual boundary
          x = (xcoords[l]+xcoords[m]+sx)/3;
          y = (ycoords[l]+ycoords[m]+sy)/3;
          z = (zcoords[l]+zcoords[m]+sz)/3;

	  // cout << "midpoint: x,y,z:" << x << " " << y << " " << z << endl;
	  
          // compute convection in center of face of dual
          // boundary
          DiscreteForm->GetCoeffFct()(1, &x, &y, &z, &param, &coeff);

	  //cout << "Koeffizienten:" << coeff[1] << " " << coeff[2] << " " << coeff[3] << endl;

          // compute flux across face of dual boundary from face j to face k
          flux = coeff[1]*nx + coeff[2]*ny + coeff[3]*nz;
          // parameter for upwind function
          t = UPWIND_FLUX_DAMP*flux*RE2;
          // compute upwind value
          GetUpwindValue(UPWIND_ORDER,t,&phi);
	  flux *= Area;
          // compute matrix entries
          switch(UPW_APPL)
          {
          case 0:
	    matrix[k][k] += flux*phi;
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(1-phi);
            /*matrix[j][j] += flux*phi;
            matrix[j][k] -= flux*phi;
            matrix[k][j] += flux*(1-phi);
            matrix[k][k] -= flux*(1-phi);*/
            break;
          case 1:  
            matrix[k][k] += flux*(-0.5+phi);
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(0.5-phi);
            break;
          default:
            OutPut("Wrong Upwind Application !!! Default is used !!!"<< endl);
            matrix[k][k] += flux*phi;
            matrix[k][j] -= flux*phi;
            matrix[j][k] += flux*(1-phi);
            matrix[j][j] -= flux*(1-phi);
            break;
          } 
        } // endfor k
        break;
     // ********** conforming discretization *******************
	
      default:
	OutPut("Upwind is not implemented for conforming discretization !" << endl);
        OutPut("Upwind works only for first order elements!" << endl);
        exit(-1);
    } // endswitch

    for(m=0;m<N_Faces;m++)
    {
      l = DOF[m];
      //cout << "DOF: " << l << endl;
      if (l < ActiveBound)
      {
        // node l is inner or Neumann node
        end = RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Faces;k++)
          {
            if(DOF[k] == ColInd[n])
            {
              //cout << m << "   " << k << endl << n << endl;
              Entries[n] += matrix[m][k];
            } // endif
          } // endfor k
        } // endfor n
      } // endif l
      else
      {
        if (l < DirichletBound)
        {
          // hanging node
          /*

          // implement later
          l -= ActiveBound;
          end = HangingRowPtr[l+1];
          for(n=HangingRowPtr[l];n<end;n++)
          {
            for(k=0;k<N_;k++)
            {
              if(DOF[k] == HangingColInd[n])
              {
                CurrentHangingEntries[n] += MatrixRow[k];
              } // endif
            } // endfor k
          } // endfor n 
          */
        }
        else
        {
          // Dirichlet node => do nothing
        }
      }
    } // endfor m
  } // endfor i
}
