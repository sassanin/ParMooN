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
#ifdef __2D__

#include <Database.h>
#include <DiscreteForm2D.h>
#include <FEFunction2D.h>
#include <FEDatabase2D.h>
#include <FE2D.h>
#include <LinAlg.h>
#include <SquareMatrix2D.h>
#include <Upwind.h>
#include <IsoBoundEdge.h>
#include <BoundComp.h>

#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>

void GetUpwindValue(double UPWIND_ORDER,double t,double *phi)
{
  double s;

  if (UPWIND_ORDER==1)                            // Samarskij upwind
  {
    if (t<0)
      *phi = 1/(2*(1-t));
    else
      *phi = (1+2*t)/(2*(1+t));
  }
  else
  {
    if (UPWIND_ORDER==0)                          // sharp upwind
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
void UpwindForNavierStokes(CoeffFct2D *Coeff, TSquareMatrix2D *sqmatrix,
TFEFunction2D *u1, TFEFunction2D *u2)
{
  double RE;
  double UPWIND_ORDER=TDatabase::ParamDB->UPWIND_ORDER;
  double UPWIND_FLUX_DAMP=TDatabase::ParamDB->UPWIND_FLUX_DAMP;
  int  UPW_APPL = TDatabase::ParamDB->UPWIND_APPLICATION;
  TBaseCell *cell;
  TCollection *coll;
  int i,j,k,l,m,n, N_Cells, N_Vertices, N_Rows;
  int n_mat, m_mat;
  double *u1vect, *u2vect;
  int *GlobalNumbers, *BeginIndex, *DOF;
  TFESpace2D *fespace;
  double FctValNQuad[4][4] =
  {
    {
      0.5, 0, 0, 0.5
    }
    ,
    {
      0.5, 0.5, 0, 0
    },
    { 0, 0.5, 0.5, 0 },
    {
      0, 0, 0.5, 0.5
    }
  };
  double FctValNTria[3][3] =
  {
    {
      0.666666666666666666667,
      -0.333333333333333333333,
      0.666666666666666666667
    },
    {
      0.666666666666666666667,
      0.666666666666666666667,
      -0.333333333333333333333
    },
    {
      -0.333333333333333333333,
      0.666666666666666666667,
      0.666666666666666666667
    }
  };
  double FctValCQuad[4][4] =
  {
    {
      0.375, 0.375, 0.125, 0.125
    },
    { 0.125, 0.375, 0.125, 0.375 },
    { 0.125, 0.125, 0.375, 0.375 },
    { 0.375, 0.125, 0.375, 0.125 }
  };
  double FctValCTria[3][3] =
  {
    {
      0.416666666666666666667,
      0.416666666666666666667,
      0.166666666666666666667
    },
    {
      0.166666666666666666667,
      0.416666666666666666667,
      0.416666666666666666667
    },
    {
      0.416666666666666666667,
      0.166666666666666666667,
      0.416666666666666666667
    }
  };
  double matrix[4][4];
  double val1, val2, val, flux;
  double xcoords[4], ycoords[4];
  double v1, v2, v3, v4;
  double sx, sy, t, phi, s;
  double nx, ny;
  double nxs[4], nys[4];
  double x,y;
  double u1old, u2old;
  double *coeffs, *params;
  int *ColInd, *RowPtr;
  double *Entries;
  int ActiveBound, DirichletBound, end;
  FE2D CurrentElement;

  coeffs = new double[6];
  params = new double[9];
  memset(params, 0, 9 * SizeOfDouble);

  // get pointers to columns, rows and entries of matrix A
  ColInd = sqmatrix->GetKCol();
  RowPtr = sqmatrix->GetRowPtr();
  Entries = sqmatrix->GetEntries();

  // get finite element space for velocity
  fespace = u1->GetFESpace2D();
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
    // get coordinates of the vertices
    for(k=0;k<N_Vertices;k++)
      cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);

    // calculate the coordinates of the center of gravity
    sx = sy = 0.0;
    for(k=0;k<N_Vertices;k++)
    {
      sx += xcoords[k];
      sy += ycoords[k];
    }
    sx /= N_Vertices;
    sy /= N_Vertices;

    // reset local matrix
    memset(matrix, 0, SizeOfDouble*16);
    CurrentElement = fespace->GetFE2D(i,cell);

    switch(CurrentElement)
    {
      // ********** nonconforming discretization *******************
      case N_P1_2D_T_A:
      case N_Q1_2D_Q_A:
      case N_Q1_2D_Q_M:
        // for all edges (nodes)
        for(k=0;k<N_Vertices;k++)
        {
          // j : edge which is 'left'
          j = (k+N_Vertices-1) % N_Vertices;
          // normal to edge of dual domain
          nx = sy-ycoords[k];
          ny = xcoords[k]-sx;
          // compute current velocity in center of edge of dual
          // boundary
          switch(N_Vertices)
          {
            // triangle
            case 3: val1 = 0.0;
            val2 = 0.0;
            for(l=0;l<N_Vertices;l++)
            {
              m = DOF[l];
              val = FctValNTria[k][l];
              val1 += u1vect[m]*val;
              val2 += u2vect[m]*val;
            }                                     // endfor l
            break;
            // quadrangle
            case 4: m=DOF[k];
            n=DOF[j];
            val1 = 0.5*(u1vect[m]+u1vect[n]);
            val2 = 0.5*(u2vect[m]+u2vect[n]);
            break;
          }                                       // endswitch
          // compute flux across edge of dual boundary
          flux = val1*nx + val2*ny;
          // get Reynolds number
          Coeff(1, &sx ,&sy, &params, &coeffs);
          RE = 1.0/coeffs[0];
          // parameter for upwind function
          t = UPWIND_FLUX_DAMP*flux*RE*0.5;
          // compute upwind value
          GetUpwindValue(UPWIND_ORDER,t,&phi);
          // compute matrix entries
          switch(UPW_APPL)
          {
            case 0:
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
        }                                         // endfor k
        break;

      default:
        cout << "Upwind works only for first order nonconforming elements!" << endl;
        exit(-1);
    }                                             // endswitch

    for(m=0;m<N_Vertices;m++)
    {
      l=DOF[m];
      // cout << "DOF: " << l << endl;
      if(l<ActiveBound)
      {
        // node l is inner or Neumann node
        end=RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Vertices;k++)
          {
            if(DOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += matrix[m][k];
            }                                     // endif
          }                                       // endfor k
        }                                         // endfor n
      }                                           // endif l
      else
      {
        if(l<DirichletBound)
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
    }                                             // endfor m
  }                                               // endfor i
  delete params;
  delete coeffs;
}                                                 // end of UpwindForNavierStokes


void UpwindForConvDiff(CoeffFct2D *Coeff,
TSquareMatrix2D *sqmatrix, double *RHS,
TFESpace2D *fespace, TDiscreteForm2D *DiscreteForm,
TFEFunction2D *u1, TFEFunction2D *u2,
int ConvIsVelo)
{
  double RE;
  double UPWIND_ORDER=TDatabase::ParamDB->UPWIND_ORDER;
  double UPWIND_FLUX_DAMP=TDatabase::ParamDB->UPWIND_FLUX_DAMP;
  int  UPW_APPL = TDatabase::ParamDB->UPWIND_APPLICATION;
  TBaseCell *cell;
  TCollection *coll;
  int i,j,k,l,m,n, N_Cells, N_Vertices, N_Rows;
  int n_mat, m_mat;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double matrix[4][4];
  double Area, val1, val2, val, flux;
  double xcoords[4], ycoords[4];
  double v1, v2, v3, v4;
  double sx, sy, t, phi, s;
  double nx, ny, x, y;
  double nxs[4], nys[4];
  int *ColInd, *RowPtr;
  double *Entries;
  int ActiveBound, DirichletBound, end;
  double *coeff = new double[6], *param = new double[9];
  FE2D CurrentElement;

  memset(param, 0, 9 * SizeOfDouble);

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
    // get coordinates of the vertices
    for(k=0;k<N_Vertices;k++)
      cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);

    CurrentElement = fespace->GetFE2D(i, cell);

    // calculate the coordinates of the center of gravity
    sx = sy = 0.0;
    for(k=0;k<N_Vertices;k++)
    {
      sx += xcoords[k];
      sy += ycoords[k];
    }
    sx /= N_Vertices;
    sy /= N_Vertices;

    // reset local matrix
    memset(matrix, 0, SizeOfDouble*16);

    switch(CurrentElement)
    {
      // ********** nonconforming discretization *******************
      case N_P1_2D_T_A:
      case N_Q1_2D_Q_A:
      case N_Q1_2D_Q_M:
        // for all edges (nodes)
        for(k=0;k<N_Vertices;k++)
        {
          // j : edge which is 'left'
          j = (k+N_Vertices-1) % N_Vertices;
          // normal to edge of dual domain
          nx = sy-ycoords[k];
          ny = xcoords[k]-sx;
          // midpoints of edge of dual boundary
          x = (xcoords[k]+sx)*0.5;
          y = (ycoords[k]+sy)*0.5;
          // compute convection in center of edge of dual
          // boundary
          if (!ConvIsVelo)
            DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
          else
          {
            u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
            u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
          }
          // compute flux across edge of dual boundary
          flux = coeff[1]*nx + coeff[2]*ny;
          // get diffusion coefficient
          Coeff(1, &sx ,&sy, &param, &coeff);
          RE = 1.0/coeff[0];
          // parameter for upwind function
          t = UPWIND_FLUX_DAMP*flux*RE*0.5;
          // compute upwind value
          GetUpwindValue(UPWIND_ORDER,t,&phi);
          // compute matrix entries
          switch(UPW_APPL)
          {
            case 0:
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
        }                                         // endfor k
        break;
        // ********** conforming discretization *******************
      case C_P1_2D_T_A:
      case C_Q1_2D_Q_A:
      case C_Q1_2D_Q_M:
        if (N_Vertices==3)
        {
          for (k=0;k<3;k++)
          {
            // centre of face (P[k]P[(k+1)%3])
            x=(xcoords[k]+xcoords[(k+1)%3])*0.5;
            y=(ycoords[k]+ycoords[(k+1)%3])*0.5;
            // normal of dual face directed to
            nxs[(2*k+1)%3] = sy-y;
            // P[(k+1)%3]
            nys[(2*k+1)%3] = x-sx;
          }
          // assemble matrix, compute integral
          for(m=0;m<N_Vertices;m++)
          {
            for(n=0;n<N_Vertices;n++)
            {
              if (n==m)
              {                                   // upwind to point P[(m-1)%3]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[(m+2)%3])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[(m+2)%3])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                // the underlying pde
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = -(coeff[1]*nxs[(2*m+2)%3]+coeff[2]*nys[(2*m+2)%3]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                if (UPW_APPL==1)
                  phi += 0.5;
                matrix[m][n] += -flux*(1.0-phi);
                // upwind to point P[(m+1)%3]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[(m+1)%3])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[(m+1)%3])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = (coeff[1]*nxs[(2*m+1)%3]+coeff[2]*nys[(2*m+1)%3]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                if (UPW_APPL==1)
                  phi += 0.5;
                matrix[m][n] += -flux*(1.0-phi);
              }
              else
              {                                   // upwind to point P[n]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[n])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[n])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = (coeff[1]*nxs[(m+n)%3]+coeff[2]*nys[(m+n)%3]);
                if (m==(n+1)%3) flux = -flux;
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                matrix[m][n] += flux*(1.0-phi);
              }
            }                                     // end for n
          }                                       // end for m
        }
        else
        {
          for (k=0;k<4;k++)
          {
            // centre of face (P[k]P[(k+1)%4])
            x=(xcoords[k]+xcoords[(k+1)%4])*0.5;
            y=(ycoords[k]+ycoords[(k+1)%4])*0.5;
            nxs[k] = sy-y;                        // normal of dual face directed to
            nys[k] = x-sx;                        // P[(k+1)%4]
          }
          // assemble matrix, compute integral
          for(m=0;m<N_Vertices;m++)
          {
            if (m<2) m_mat=m;
            if (m==2) m_mat=3;
            if (m==3) m_mat=2;
            for(n=0;n<N_Vertices;n++)
            {
              if (n<2) n_mat=n;
              if (n==2) n_mat=3;
              if (n==3) n_mat=2;
              if (n==m)
              {                                   // upwind to point P[(m+3)%4]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[(m+3)%4])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[(m+3)%4])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                //OutPut(coeff[1] << " " << coeff[2] << " ");
                flux = -(coeff[1]*nxs[(m+3)%4]+coeff[2]*nys[(m+3)%4]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                matrix[m_mat][n_mat] += -flux*(1.0-phi);
                // upwind to point P[(m+1)%4]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[(m+1)%4])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[(m+1)%4])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = (coeff[1]*nxs[m%4]+coeff[2]*nys[m%4]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                matrix[m_mat][n_mat] += -flux*(1.0-phi);
              }
              if(n==(m+1)%4)
              {                                   // upwind to point P[(m+1)%4]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[n])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[n])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = (coeff[1]*nxs[m%4]+coeff[2]*nys[m%4]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                matrix[m_mat][n_mat] += flux*(1.0-phi);
              }
              if(n==(m+3)%4)
              {                                   // upwind to point P[(m+3)%4]
                // coordinates on dual boundary
                x= (xcoords[m]+xcoords[n])*0.25+sx*0.5;
                // coordinates on dual boundary
                y= (ycoords[m]+ycoords[n])*0.25+sy*0.5;
                if (!ConvIsVelo)
                  // get parameter function of the underlying pde
                  DiscreteForm->GetCoeffFct()(1, &x, &y, &param, &coeff);
                else
                {
                  u1->FindGradientLocal(cell,i,x,y,&coeff[1]);
                  u2->FindGradientLocal(cell,i,x,y,&coeff[2]);
                }
                flux = -(coeff[1]*nxs[n%4]+coeff[2]*nys[n%4]);
                // get diffusion coefficient
                Coeff(1, &sx ,&sy, &param, &coeff);
                RE = 1.0/coeff[0];
                // parameter for upwind function
                t= UPWIND_FLUX_DAMP*flux*RE*0.5;
                GetUpwindValue(UPWIND_ORDER,t,&phi);
                matrix[m_mat][n_mat] += flux*(1.0-phi);
              }
            }                                     // end for n
          }                                       // end for m
        }
        break;

      default:
        OutPut("Upwind works only for first order elements!" << endl);
        exit(-1);
    }                                             // endswitch
    for(m=0;m<N_Vertices;m++)
    {
      l = DOF[m];
      // cout << "DOF: " << l << endl;
      if (l < ActiveBound)
      {
        // node l is inner or Neumann node
        end = RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Vertices;k++)
          {
            if(DOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += matrix[m][k];
            }                                     // endif
          }                                       // endfor k
        }                                         // endfor n
      }                                           // endif l
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
    }                                             // endfor m
  }                                               // endfor i
  delete param;
  delete coeff;
}                                                 // UpwindForConvDiff

/******************************************************************************/
//
// IMPROVED MIZUKAMI-HUGHES METHOD (Knobloch, CMAME 2007)
//
/******************************************************************************/

void ComputeParametersMizukamiHughes(TBaseCell *cell, int cell_no,
TFEFunction2D *u, CoeffFct2D *Coeffs,
BoundCondFunct2D *BoundaryCondition,
int *dof, int ActiveBound,
double *c_mh)
{
  int i, N_Edges, found, comp;
  double *coeff, sx, sy, x[3], y[3], b1, b2, val[3], ux, uy, norm_b;
  double phi_x[3], phi_y[3], phi_z, v1, v2, w1, w2, s1, s2;
  double t0, t1, wx, wy, hx, hy, h21x, h21y, h31x, h31y, deter, alpha;
  double beta, h21, h31, v, vx, vy, vx_orth, vy_orth, r2, phi;
  TJoint *joint;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  BoundCond Cond;
  TBoundComp *BoundComp;

  coeff = new double[6];

  // compute coordinates of vertices
  N_Edges = cell->GetN_Edges();
  if (N_Edges != 3)
  {
    OutPut("Improved Mizukami-Hughes Method implemented only for triangles !!!"
      << endl);
    exit(4711);
  }
  sx = sy = 0;
  // compute derivatives for basis functions
  for (i=0;i<N_Edges; i++)
  {
    x[i] = cell->GetVertex(i)->GetX();
    y[i] = cell->GetVertex(i)->GetY();
    sx += x[i];
    sy += y[i];
  }
  sx /= N_Edges;
  sy /= N_Edges;
  // get coefficients, are assumed to be constant on the triangle
  Coeffs(1, &sx, &sy, NULL, &coeff);
  b1 = coeff[1];
  b2 = coeff[2];

  norm_b = b1*b1+b2*b2;
  norm_b = sqrt(norm_b);
  // locally no convection
  if (norm_b < 1e-9)
  {
    c_mh[0] = c_mh[1] = c_mh[2] = 0;
    delete coeff;
    return;
  }

  // find (constant) gradient of discrete solution
  u->FindGradientLocal(cell, cell_no, sx, sy, val);
  ux = val[1];
  uy = val[2];

  // found gives the zone where b point into
  // 0 - vertex zone 0;        1 - edge zone 0
  // 2 - vertex zone 1;        3 - edge zone 1
  // 4 - vertex zone 2;        5 - edge zone 2
  found = -1;
  // find the vertex where b points into the vertex zone or edge zone
  for (i=0; i<N_Edges; i++)
  {
    // gradient of the Lagrange fe function = normal vector of the plane
    v1 = x[(i+1)%N_Edges] - x[i];
    v2 = y[(i+1)%N_Edges] - y[i];

    w1 = x[(i+2)%N_Edges] - x[i];
    w2 = y[(i+2)%N_Edges] - y[i];

    phi_z = v1*w2 - w1*v2;
    phi_x[i] = v2 - w2;
    phi_x[i] /= phi_z;
    phi_y[i] = w1 - v1;
    phi_y[i] /= phi_z;
    val[i] = b1*phi_x[i] + b2*phi_y[i];
  }

  for (i=0;i<N_Edges;i++)
  {
    // vertex zone
    if ((val[i] > 0) && (val[(i+1)%N_Edges] <= 0) && (val[(i+2)%N_Edges] <= 0))
    {
      if (found > -1)
      {
        OutPut(cell_no << " already something found " << found <<
          " found new in vertex zone " << i << endl);
        exit(4711);
      }
      found = 2*i;
    }
    // edge zone
    if ((val[i] < 0) && (val[(i+1)%N_Edges] > 0) && (val[(i+2)%N_Edges] > 0))
    {
      if (found > -1)
      {
        OutPut(cell_no << " already something found " << found <<
          " found new in edge zone " << i << endl);
        exit(4711);
      }
      found = 2*i+1;
    }
  }
  if (found==-1)
  {
    OutPut("Mizukami-Hughes: nothing found !!!"<<endl);
    exit(4711);
  }

  // vertex zone
  switch(found)
  {
    case 0:
      c_mh[0] = 2.0/3.0;
      c_mh[1] = -1.0/3.0;
      c_mh[2] = -1.0/3.0;
      delete coeff;
      return;
    case 2:
      c_mh[0] = -1.0/3.0;
      c_mh[1] = 2.0/3.0;
      c_mh[2] = -1.0/3.0;
      delete coeff;
      return;
    case 4:
      c_mh[0] = -1.0/3.0;
      c_mh[1] = -1.0/3.0;
      c_mh[2] = 2.0/3.0;
      delete coeff;
      return;
  }

  // check if there is a common node with the Dirichlet boundary
  // THIS IS ONLY FOR THE UNIT SQUARE
  for (i=0;i<N_Edges;i++)
  {
    if (dof[i] >= ActiveBound)
    {
      c_mh[0] = -1.0/3.0;
      c_mh[1] = -1.0/3.0;
      c_mh[2] = -1.0/3.0;
      delete coeff;
      return;
    }
  }

  // here should be another condition on the whole triangulation !!!!
  // right angles, iso-sceles triangles are excluded

  // STILL TO INCLUDE A CONDITION !!!

  // compute vector of the edges
  switch(found)
  {
    case 1:
      hx = x[2] - x[1];
      hy = y[2] - y[1];
      // edge which determines VZ_2
      h21x = x[1] - x[0];
      h21y = y[1] - y[0];
      // edge which determines VZ_3
      h31x = x[2] - x[0];
      h31y = y[2] - y[0];
      break;
    case 3:
      hx = x[0] - x[2];
      hy = y[0] - y[2];
      // edge which determines VZ_2
      h21x = x[2] - x[1];
      h21y = y[2] - y[1];
      // edge which determines VZ_3
      h31x = x[0] - x[1];
      h31y = y[0] - y[1];
      break;
    case 5:
      hx = x[1] - x[0];
      hy = y[1] - y[0];
      // edge which determines VZ_2
      h21x = x[0] - x[2];
      h21y = y[0] - y[2];
      // edge which determines VZ_3
      h31x = x[1] - x[2];
      h31y = y[1] - y[2];
      break;
  }

  // next condition
  if (fabs(b1*ux + b2*uy) < 1e-16)
  {
    switch(found)
    {
      case 1:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = 1.0/6.0;
        c_mh[2] = 1.0/6.0;
        delete coeff;
        return;
      case 3:
        c_mh[0] = 1.0/6.0;
        c_mh[1] = -1.0/3.0;
        c_mh[2] = 1.0/6.0;
        delete coeff;
        return;
      case 5:
        c_mh[0] = 1.0/6.0;
        c_mh[1] = 1.0/6.0;
        c_mh[2] = -1.0/3.0;
        delete coeff;
        return;
    }
  }

  // vector orthogonal to \nabla u^h
  wx = -uy;
  wy = ux;

  // check if V_2 is empty, compute parameter for intersection
  deter = -wx*h21y+wy*h21x;
  if (fabs(deter)<1e-20)
  {
    // parallel lines, no intersection
    //OutPut("Determinant vanishes !!!"<<endl);
    alpha = -4711;
  }
  else
  {
    alpha = -wx*b2+wy*b1;
    alpha /= deter;
  }

  deter = -wx*h31y+wy*h31x;
  if (fabs(deter)<1e-20)
  {
    // parallel lines, no intersection
    //OutPut("Determinant vanishes !!!"<<endl);
    beta = -4711;
  }
  else
  {
    beta = -wx*b2+wy*b1;
    beta /= deter;
    //	alpha_a = -(sx-b1)*h31y+(sy-b2)*h31x;
    //alpha_a /= deter;
  }
  // OutPut(b1+alpha_a*wx << " " << sx+beta*h31x << "   " <<
  //	   b2+alpha_a*wy << " " << sy+beta*h31y << "   ");
  //OutPut("a " << alpha << " b " << beta << endl);
  if ((alpha<=0)&&(beta<=0))
  {
    OutPut("negative " << found << " S " << sx << " " << sy <<
      " h21 " << h21x << " " << h21y <<  " h31 " << h31x << " " << h31y <<
      " w " << wx << " " << wy << " b " << b1 << " " << b2 <<
      " a " << alpha << " " << beta << endl);
    exit(4711);
  }
  // V_2 not empty, V_3 empty
  if ((alpha > 0) && (beta <= 0))
  {
    switch(found)
    {
      case 1:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = 2.0/3.0;
        c_mh[2] = -1.0/3.0;
        delete coeff;
        return;
      case 3:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = -1.0/3.0;
        c_mh[2] = 2.0/3.0;
        delete coeff;
        return;
      case 5:
        c_mh[0] = 2.0/3.0;
        c_mh[1] = -1.0/3.0;
        c_mh[2] = -1.0/3.0;
        delete coeff;
        return;
    }
  }
  // V_2 empty, V_3 not empty
  if ((alpha <= 0) && (beta >  0))
  {
    switch(found)
    {
      case 1:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = -1.0/3.0;
        c_mh[2] = 2.0/3.0;
        delete coeff;
        return;
      case 3:
        c_mh[0] = 2.0/3.0;
        c_mh[1] = -1.0/3.0;
        c_mh[2] = -1.0/3.0;
        delete coeff;
        return;
      case 5:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = 2.0/3.0;
        c_mh[2] = -1.0/3.0;
        delete coeff;
        return;
    }
  }

  // V_2 and V_3 not empty
  // define vectors of the method
  s1 = b1/norm_b;
  s2 = b2/norm_b;
  h21 = sqrt(h21x*h21x+h21y*h21y);
  h21x /= h21;
  h21y /= h21;
  h31 = sqrt(h31x*h31x+h31y*h31y);
  h31x /= h31;
  h31y /= h31;
  vx = h21x + h31x;
  vy = h21y + h31y;
  v = sqrt(vx*vx + vy*vy);
  vx /= v;
  vy /= v;

  // change sign of w if necessary
  if (wx*vx+wy*vy < 0)
  {
    wx = -wx;
    wy = -wy;
  }

  vx_orth = -vy;
  vy_orth = vx;
  // change sign if necessary
  if (vx_orth*h31x+vy_orth*h31y < 0)
  {
    vx_orth = -vx_orth;
    vy_orth = -vy_orth;
  }

  // first case
  if (wx*vx_orth + wy*vy_orth < 0)
  {
    if (b1*h21x+b2*h21y < 0)
      r2 = -1;
    else
    {
      if (b1*h21x+b2*h21y > 0)
        r2 = 1;
      else
        r2 = 0;
    }
    t0 = fabs(s1*h21y-s2*h21x);
    t1 = fabs(vx*h21y-vy*h21x);
    r2 = t0/t1 + 1 - r2;
    if (r2>1)
      r2 = 1;
    t0 = fabs(wx*h21y - wy*h21x);
    t1 = vx*h21x + vy*h21y;
    t1 *= r2;
    phi = 2*t0/t1;
    if (phi>1)
      phi = 1;
    t0 = (h21x - h31x)* s1 + (h21y - h31y) * s2;
    t1 = 1 - h21x*h31x - h21y * h31y;
    r2 = -1.0/3.0 + phi *(1 + t0/t1) / 2.0;
    switch(found)
    {
      case 1:
        c_mh[0] = -1.0/3.0;
        c_mh[1] = r2;
        c_mh[2] = 1.0/3.0 - c_mh[1];
        break;
      case 3:
        c_mh[1] = -1.0/3.0;
        c_mh[2] = r2;
        c_mh[0] = 1.0/3.0 - c_mh[2];
        break;
      case 5:
        c_mh[2] = -1.0/3.0;
        c_mh[0] = r2;
        c_mh[1] = 1.0/3.0 - c_mh[0];
        break;
    }
  }
  else
    // second case
  {
    if (b1*h31x+b2*h31y < 0)
      r2 = -1;
    else
    {
      if (b1*h31x+b2*h31y > 0)
        r2 = 1;
      else
        r2 = 0;
    }
    t0 = fabs(s1*h31y-s2*h31x);
    t1 = fabs(vx*h31y-vy*h31x);
    r2 = t0/t1 + 1 - r2;
    if (r2>1)
      r2 = 1;
    t0 = fabs(wx*h31y - wy*h31x);
    t1 = vx*h31x + vy*h31y;
    t1 *= r2;
    phi = 2*t0/t1;
    if (phi>1)
      phi = 1;
    t0 = (h31x - h21x)* s1 + (h31y - h21y) * s2;
    t1 = 1 - h21x*h31x - h21y * h31y;
    r2 = -1.0/3.0 + phi *(1 + t0/t1) / 2.0;
    switch(found)
    {
      case 1:
        c_mh[2] = r2;
        c_mh[1] = 1.0/3.0 - c_mh[2];
        c_mh[0] = -1.0/3.0;
        break;
      case 3:
        c_mh[0] = r2;
        c_mh[2] = 1.0/3.0 - c_mh[0];
        c_mh[1] = -1.0/3.0;
        break;
      case 5:
        c_mh[1] = r2;
        c_mh[0] = 1.0/3.0 - c_mh[1];
        c_mh[2] = -1.0/3.0;
        break;
    }
  }
  delete coeff;
}


/*******************************************************************************/
//
// MIZUKAMI-HUGHES METHOD FOR CONVECTION--DIFFUSION EQUATIONS
//
/*******************************************************************************/

void MizukamiHughes(TSquareMatrix2D *sqmatrix,
double *RHS,
TFESpace2D *fespace,
TFEFunction2D *u,
CoeffFct2D *Coeffs,
BoundCondFunct2D *BoundaryCondition)
{
  TBaseCell *cell;
  TCollection *coll;
  int i,j,k,l,m,n, N_Cells, N_Vertices, N_Rows;
  int n_mat, m_mat;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double matrix[3][3];
  double c_mh[3], b1, b2, val, area;
  double xcoords[3], ycoords[3];
  double v1, v2, v3, v4;
  double sx, sy, t, phi, s;
  double nx, ny, x, y;
  double nxs[4], nys[4];
  int *ColInd, *RowPtr, jj, ii;
  double *Entries;
  int ActiveBound, DirichletBound, end;
  double *coeff = new double[6], *param = new double[1];
  FE2D CurrentElement;
  BaseFunct2D *BaseFuncts, LocBF;
  int *N_BaseFunct, LocN_BF;

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

  // information on the base functions
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

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
    // check if P1 fe
    CurrentElement = fespace->GetFE2D(i, cell);
    if (CurrentElement != C_P1_2D_T_A)
    {
      OutPut("Mizukami-Hughes method implemented only for P_1 and not for fe id "
        << CurrentElement << endl);
      exit(4711);
    }

    // compute local basis functions
    LocN_BF = N_BaseFunct[CurrentElement];
    LocBF = BaseFuncts[CurrentElement];

    // # of vertices
    N_Vertices = cell->GetN_Vertices();
    // get coordinates of the vertices
    for(k=0;k<N_Vertices;k++)
      cell->GetVertex(k)->GetCoords(xcoords[k], ycoords[k]);

    // calculate the coordinates of the center of gravity
    sx = sy = 0.0;
    for(k=0;k<N_Vertices;k++)
    {
      sx += xcoords[k];
      sy += ycoords[k];
    }
    sx /= N_Vertices;
    sy /= N_Vertices;
    // get coefficients
    Coeffs(1, &sx, &sy, NULL, &coeff);
    b1 = coeff[1];
    b2 = coeff[2];

    // compute weights
    c_mh[0] = c_mh[1] = c_mh[2] = -4711;
    ComputeParametersMizukamiHughes(cell, i, u,
      Coeffs, BoundaryCondition, DOF, ActiveBound, c_mh);

    // reset local matrix
    memset(matrix, 0, SizeOfDouble*9);

    // loop over ansatz function
    for (jj=0;jj<3;jj++)
    {
      val = b1*(ycoords[(jj+1)%3] - ycoords[(jj+2)%3])
        + b2 * (xcoords[(jj+2)%3] - xcoords[(jj+1)%3]);
      val /= 2;
      // loop over test function
      for (ii=0;ii<3;ii++)
      {
        matrix[ii][jj] = val * c_mh[ii];
      }
    }

    for(m=0;m<N_Vertices;m++)
    {
      l = DOF[m];
      // cout << "DOF: " << l << endl;
      if (l < ActiveBound)
      {
        // node l is inner or Neumann node
        end = RowPtr[l+1];
        for(n=RowPtr[l];n<end;n++)
        {
          for(k=0;k<N_Vertices;k++)
          {
            if(DOF[k] == ColInd[n])
            {
              // cout << m << "   " << k << endl << n << endl;
              Entries[n] += matrix[m][k];
            }                                     // endif
          }                                       // endfor k
        }                                         // endfor n
      }                                           // endif l
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
    }                                             // endfor m

    // update rhs
    for (ii=0;ii<3;ii++)
    {
      area = (xcoords[1]-xcoords[0])*(ycoords[2]-ycoords[0])
        - (xcoords[2]-xcoords[0])*(ycoords[1]-ycoords[0]);
      area = fabs(area)/2;
      l = DOF[ii];
      // only for rhs not on Dirichlet boundary ???
      if (l < ActiveBound)
      {
        // edge midpoint rule
        val = 0;
        for (jj=0;jj<3;jj++)
        {
          x = (xcoords[jj]+xcoords[(jj+1)%3])/2;
          y = (ycoords[jj]+ycoords[(jj+1)%3])/2;
          Coeffs(1, &x, &y, NULL, &coeff);
          val += coeff[4];
        }
        val *= area * c_mh[ii] / 3.0;
        RHS[l] += val;
      }
    }
  }
}

#endif // #ifdef __2D__
