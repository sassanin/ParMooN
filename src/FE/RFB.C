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
// RFB.C
//
// Purpose:     approximate solution of RFB problems in two-level method
//              using Q_1 finite elements
//
// Author:      Volker John       start 2005/11/28
//
// =======================================================================

#include <Constants.h>
#include <Database.h>
#include <LinAlg.h>
#include <Convolution.h>
#include <stdlib.h>
#include <string.h>

#ifdef __2D__
#include <FEDatabase2D.h>
#include <FEFunction2D.h>
#endif

#ifdef __3D__
#include <FEDatabase3D.h>
#include <FEFunction3D.h>
#include <TNSE3D_Routines.h>
#endif

#ifdef __2D__
void ApproximateRFBSolutionQuadNSE2DOri(TCollection *Coll, TFEFunction2D *u1,
TFEFunction2D *u2, CoeffFct2D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, dof[4], N_U, N_;
  int N_sub2, N_nodes, N_nodes2, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[4],y[4], a[16], b[16], *a_sub, *rhs_sub, test_bary;
  double xc[4], yc[4], area, val, xs, ys, u1_values[3], u2_values[3];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff, *xlocal, *ylocal;
  double *integral_1, *integral_2, xi, eta, *test_value;
  double rhs1_x, rhs2_x, rhs1_y, rhs2_y;
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  TFESpace2D *fespace;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;

  N_sub2 = N_sub * N_sub;
  N_nodes = N_sub-1;
  N_nodes2 = N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes2*N_nodes2];
  memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);

  rhs_sub = new double[4 * N_nodes2];
  memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

  coeff = new double[3];

  xlocal = new double[N_nodes2];
  ylocal = new double[N_nodes2];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace2D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 4)
    {
      OutPut("RFB stabilization only for quadrilateral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      //OutPut("coords " << x[j] << " " << y[j] << endl);
    }

    // area of parallelogramm with vector product
    area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub2 ;
    //OutPut(" " << area << endl);
    memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);
    memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        //OutPut("cell " << ic*N_sub+jc << " : ");
        // compute coordinates of submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);
        //OutPut(xc[0] << " " << yc[0] << " ");
        xc[1] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[1] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[1] /= (N_sub*N_sub);
        yc[1] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[1] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[1] /= (N_sub*N_sub);
        //OutPut(xc[1] << " " << yc[1] << " ");
        xc[2] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[2] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[2] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[2] /= (N_sub*N_sub);
        //OutPut(xc[2] << " " << yc[2] << " ");
        xc[3] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[3] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[3] /=  (N_sub*N_sub);
        yc[3] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[3] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[3] /=  (N_sub*N_sub);
        //OutPut(xc[3] << " " << yc[3] << endl);

        // initialize rhs
        memset(b,0,16*SizeOfDouble);

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<4;jj++)
        {
          a[4*jj] = 1;
          a[4*jj+1] = xc[jj];
          a[4*jj+2] = yc[jj];
          a[4*jj+3] = xc[jj]*yc[jj];
          b[5*jj] = 1;
        }
        //for (j=0;j<16;j++)
        //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

        // solve system for the coefficients of the bilinear function
        SolveMultipleSystems(a,b,4,4,4,4);
        /*for (j=0;j<N_V;j++)
        OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
        b[4*j+2] << " " <<  b[4*j+3] << endl);
        */

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        // enumeration is lexicographically since the indizes of the dof
        // are the same as in the enumeration of the vertices
        dof[0] = (ic-1)*N_nodes+jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + N_nodes;
        dof[3] = dof[0] + N_nodes;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[2] = -1;
          dof[3] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[1] = -1;
          dof[2] = -1;
        }
        //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << endl);
        // just for graphical output
        /*for (j=0;j<4;j++)
        {
          if (dof[j] != -1)
          {
            xlocal[dof[j]] = xc[j];
            ylocal[dof[j]] = yc[j];
          }
        }*/
        // barycenter
        xs = xc[0] + xc[1] + xc[2] + xc[3];
        xs /= 4;
        ys = yc[0] + yc[1] + yc[2] + yc[3];
        ys /= 4;

        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);
        // just for comparison with a known shape of solution
        //u1_values[0] = 1;
        //u2_values[0] = 0;

        // compute value of rhs in (xs,ys)
        Coeffs(1, &xs, &ys, NULL, &coeff);
        //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);
        // just for comparison with a known shape of solution
        //coeff[1] = 1;

        // update matrix entries
        // the matrix is stored column wise
        // ii -- test function
        for(ii=0;ii<4;ii++)
        {
          if (dof[ii] == -1)
            continue;

          // value of test function in bary center
          test_bary = b[4*ii] +  b[4*ii+1] * xs + b[4*ii+2] * ys + b[4*ii+3] * xs * ys;
          //OutPut("test bary " << test_bary << endl);

          // compute first rhs
          val = u1_values[0]*u1_values[1] +  u2_values[0]*u1_values[2];
          val *= test_bary;
          rhs_sub[dof[ii]] -= area*val;
          // compute second rhs
          val = u1_values[0]*u2_values[1] +  u2_values[0]*u2_values[2];
          val *= test_bary;
          rhs_sub[dof[ii]+ N_nodes2] -= area*val;

          // compute third and fourth rhs
          rhs_sub[dof[ii]+ 2*N_nodes2] += area*coeff[1]*test_bary;
          rhs_sub[dof[ii]+ 3*N_nodes2] += area*coeff[2]*test_bary;

          // jj -- ansatz function
          for (jj=0;jj<4;jj++)
          {
            if (dof[jj] == -1)
              continue;
            // entry (jj,ii) exists
            // compute index in array containing the entry

            index = dof[ii]*N_nodes2+dof[jj];
            //OutPut(" ii " << ii << " jj " << jj << " ind " << index )
            // add Laplacian
            // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
            // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
            // approximate integral by midpoint rule
            // x part of the gradient
            val = (b[4*jj+1]+ b[4*jj+3]*ys)*(b[4*ii+1]+ b[4*ii+3]*ys);
            // y part of the gradient
            val += (b[4*jj+2]+ b[4*jj+3]*xs)*(b[4*ii+2]+ b[4*ii+3]*xs);
            val *= eps * area;
            //OutPut(" val " << val << endl);
            a_sub[index] += val;

            // add convective term
            // convection times gradient of ansatz fct.
            val = u1_values[0] * (b[4*jj+1]+ b[4*jj+3]*ys)
              + u2_values[0] * (b[4*jj+2]+ b[4*jj+3]*xs);
            // multiplied with test fct.
            val *= test_bary;
            val *= area;
            a_sub[index] += val;

          }                      // end jj
        }                        // end ii
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*for (jj=0;jj<N_nodes2;jj++)
    {
    for (ii=0;ii<N_nodes2;ii++)
    {
      if (fabs(a_sub[jj+ii*N_nodes2])<1e-10)
    a_sub[jj+ii*N_nodes2] = 0;
      OutPut(a_sub[jj+ii*N_nodes2] << " ");
    }
    OutPut(endl);
    }*/

    // just for comparison with a known shape of solution
    /*for (j=0;j< N_nodes2 ; j++)
      OutPut("rhs " << rhs_sub[j+2*N_nodes2] << endl); */

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    SolveMultipleSystems(a_sub, rhs_sub, N_nodes2, N_nodes2, N_nodes2 , 4);

    // just for comparison with a known shape of solution
    /*    for (j=0;j< N_nodes2 ; j++)
      OutPut(rhs_sub[j] << " " << rhs_sub[j+N_nodes2] << " " <<
      rhs_sub[j+2*N_nodes2] << " " << rhs_sub[j+3*N_nodes2]  << endl);*/
    /*	OutPut(" sol " << xlocal[j] << " " << ylocal[j] << " " <<
      rhs_sub[j+2*N_nodes2] << endl);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    // add contributions from the global function and the rhs
    //OutPut("norm vor " << sqrt(Ddot(4*N_nodes2,rhs_sub,rhs_sub)) );
    for (j=0;j< 2 *N_nodes2 ; j++)
      rhs_sub[j] += rhs_sub[j + 2*N_nodes2];
    //rhs_sub[j] = rhs_sub[j + 2*N_nodes2];

    //OutPut("norm nach " << sqrt(Ddot(2*N_nodes2,rhs_sub,rhs_sub)) << endl);

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase2D::GetFE2D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    test_value = new double[N_];
    integral_1 = new double[N_];
    integral_2 = new double[N_];
    memset(integral_1,0,N_*SizeOfDouble);
    memset(integral_2,0,N_*SizeOfDouble);

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        // compute coordinates of submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);

        xc[1] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[1] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[1] /= (N_sub*N_sub);
        yc[1] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[1] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[1] /= (N_sub*N_sub);

        xc[2] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[2] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[2] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[2] /= (N_sub*N_sub);

        xc[3] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[3] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[3] /=  (N_sub*N_sub);
        yc[3] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[3] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[3] /=  (N_sub*N_sub);

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        dof[0] = (ic-1)*N_nodes+jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + N_nodes;
        dof[3] = dof[0] + N_nodes;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[2] = -1;
          dof[3] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[1] = -1;
          dof[2] = -1;
        }
        // barycenter
        xs = xc[0] + xc[1] + xc[2] + xc[3];
        xs /= 4;
        ys = yc[0] + yc[1] + yc[2] + yc[3];
        ys /= 4;

        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);

        // compute the value of the test functions in the bary center
        // of the submesh cell
        TFEDatabase2D::GetRefFromOrig(RefTrans, xs, ys, xi, eta);
        bf->GetDerivatives(D00, xi, eta, test_value);
        //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
        //       test_value[2]  << " " << test_value[3] << endl);

        // compute the derivative of the solution of the RFB problem
        // in the barycenter of the submesh cell
        // compute for this purpose the bilinear function on the mesh cell

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<4;jj++)
        {
          a[4*jj] = 1;
          a[4*jj+1] = xc[jj];
          a[4*jj+2] = yc[jj];
          a[4*jj+3] = xc[jj]*yc[jj];
          if (dof[jj]!=-1)
          {
            b[jj] =  rhs_sub[dof[jj]];
            b[jj+4] = rhs_sub[dof[jj]+N_nodes2];
          }
          else
            b[jj] = b[jj+4] = 0;
        }

        // solve system for the coefficients of the bilinear function
        SolveMultipleSystems(a,b,4,4,4,2);
        /*for (jj=0;jj<2;jj++)
          OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
          b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/
        rhs1_x = b[1] + b[3]*ys;
        rhs1_y = b[2] + b[3]*xs;
        rhs2_x = b[5] + b[7]*ys;
        rhs2_y = b[6] + b[7]*xs;

        //OutPut("grad " << xs << " " <<  ys << " " << rhs1_x << " " << rhs1_y << endl);

        // update the integrals
        for (j=0;j<N_;j++)
        {
          val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y)*test_value[j]*area;
          integral_1[j] += val;
          val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y)*test_value[j]*area;
          integral_2[j] += val;
        }
      }                          // jc
    }                            // ic
    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
    }
    delete test_value;
    delete integral_1;
    delete integral_2;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;
  delete xlocal;
  delete ylocal;

}


void ApproximateRFBSolutionQuadNSE2D(TCollection *Coll, TFEFunction2D *u1,
TFEFunction2D *u2, CoeffFct2D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, dof[4], N_U, N_, ij;
  int N_sub2, N_nodes, N_nodes2, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[4],y[4], a[16], b[16], *a_sub, *rhs_sub, test[4];
  double xc[4], yc[4], area, val, xs, ys, u1_values[3], u2_values[3];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff, *xlocal, *ylocal;
  double *integral_1, *integral_2, xi, eta, *test_value;
  double rhs1_x, rhs2_x, rhs1_y, rhs2_y;
  double gradx[4], grady[4], gradx_test[4], grady_test[4];
  double x_quad[4], y_quad[4], val_test, val1, val2;
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  TFESpace2D *fespace;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;

  N_sub2 = N_sub * N_sub;
  N_nodes = N_sub-1;
  N_nodes2 = N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes2*N_nodes2];
  memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);

  rhs_sub = new double[4 * N_nodes2];
  memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

  coeff = new double[3];

  xlocal = new double[N_nodes2];
  ylocal = new double[N_nodes2];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace2D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 4)
    {
      OutPut("RFB stabilization only for quadrilateral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      //OutPut("coords " << x[j] << " " << y[j] << endl);
    }

    // area of parallelogramm with vector product
    area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub2 ;
    //OutPut(" " << area << endl);
    memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);
    memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        //OutPut("cell " << ic*N_sub+jc << " : ");
        // compute coordinates of submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);
        //OutPut(xc[0] << " " << yc[0] << " ");
        xc[1] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[1] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[1] /= (N_sub*N_sub);
        yc[1] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[1] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[1] /= (N_sub*N_sub);
        //OutPut(xc[1] << " " << yc[1] << " ");
        xc[2] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[2] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[2] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[2] /= (N_sub*N_sub);
        //OutPut(xc[2] << " " << yc[2] << " ");
        xc[3] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[3] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[3] /=  (N_sub*N_sub);
        yc[3] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[3] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[3] /=  (N_sub*N_sub);
        //OutPut(xc[3] << " " << yc[3] << endl);

        // initialize rhs
        memset(b,0,16*SizeOfDouble);

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<4;jj++)
        {
          a[4*jj] = 1;
          a[4*jj+1] = xc[jj];
          a[4*jj+2] = yc[jj];
          a[4*jj+3] = xc[jj]*yc[jj];
          b[5*jj] = 1;
        }
        //for (j=0;j<16;j++)
        //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

        // solve system for the coefficients of the bilinear function
        SolveMultipleSystems(a,b,4,4,4,4);
        /*for (j=0;j<N_V;j++)
        OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
        b[4*j+2] << " " <<  b[4*j+3] << endl);
        */

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        // enumeration is lexicographically since the indizes of the dof
        // are the same as in the enumeration of the vertices
        dof[0] = (ic-1)*N_nodes+jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + N_nodes;
        dof[3] = dof[0] + N_nodes;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[2] = -1;
          dof[3] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[1] = -1;
          dof[2] = -1;
        }
        //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << endl);
        // just for graphical output
        /*for (j=0;j<4;j++)
        {
          if (dof[j] != -1)
          {
            xlocal[dof[j]] = xc[j];
            ylocal[dof[j]] = yc[j];
          }
        }*/
        // barycenter
        xs = xc[0] + xc[1] + xc[2] + xc[3];
        xs /= 4;
        ys = yc[0] + yc[1] + yc[2] + yc[3];
        ys /= 4;
        // compute quadrature points
        for (ij=0;ij<4;ij++)
        {
          x_quad[ij] = (xc[1]-xc[0])* gauss2_x[ij] + (xc[3]-xc[0])*gauss2_y[ij] + xc[1]+xc[3];
          x_quad[ij] /= 2;
          y_quad[ij] = (yc[1]-yc[0])* gauss2_x[ij] + (yc[3]-yc[0])*gauss2_y[ij] + yc[1]+yc[3];
          y_quad[ij] /= 2;
        }
        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);
        // just for comparison with a known shape of solution
        //u1_values[0] = 1;
        //u2_values[0] = 0;

        // compute value of rhs in (xs,ys)
        Coeffs(1, &xs, &ys, NULL, &coeff);
        //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);
        // just for comparison with a known shape of solution
        //coeff[1] = 1;

        // update matrix entries
        // the matrix is stored column wise
        // ii -- test function
        for(ii=0;ii<4;ii++)
        {
          if (dof[ii] == -1)
            continue;
          // loop over the quadrature points
          val_test = 0;
          for(ij=0;ij<4;ij++)
          {
            // value of test function in bary center
            test[ij] = b[4*ii] +  b[4*ii+1] * x_quad[ij]
              + b[4*ii+2] * y_quad[ij] + b[4*ii+3] * x_quad[ij] * y_quad[ij];
            val_test += test[ij];
            gradx_test[ij] = b[4*ii+1]+ b[4*ii+3]*y_quad[ij];
            grady_test[ij] = b[4*ii+2]+ b[4*ii+3]*x_quad[ij];
          }
          val_test *= area/4.0;

          // compute first rhs
          val = u1_values[0]*u1_values[1] +  u2_values[0]*u1_values[2];
          val *= val_test;
          rhs_sub[dof[ii]] -= val;
          // compute second rhs
          val = u1_values[0]*u2_values[1] +  u2_values[0]*u2_values[2];
          val *= val_test;
          rhs_sub[dof[ii]+ N_nodes2] -= val;

          // compute third and fourth rhs
          rhs_sub[dof[ii]+ 2*N_nodes2] += val_test*coeff[1];
          rhs_sub[dof[ii]+ 3*N_nodes2] += val_test*coeff[2];

          // jj -- ansatz function
          for (jj=0;jj<4;jj++)
          {
            if (dof[jj] == -1)
              continue;
            // entry (jj,ii) exists
            // compute index in array containing the entry

            index = dof[ii]*N_nodes2+dof[jj];
            //OutPut(" ii " << ii << " jj " << jj << " ind " << index )

            val = 0;
            // loop over the quadrature points
            for(ij=0;ij<4;ij++)
            {
              // add Laplacian
              // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
              // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
              // x part of the gradient
              gradx[ij] = b[4*jj+1]+ b[4*jj+3]*y_quad[ij];
              grady[ij] = b[4*jj+2]+ b[4*jj+3]*x_quad[ij];
              val += gradx[ij] * gradx_test[ij] + grady[ij] * grady_test[ij];
            }
            val *= eps * area/ 4;
            //OutPut(" val " << val << endl);
            a_sub[index] += val;

            val = 0;
            // loop over the quadrature points
            for(ij=0;ij<4;ij++)
            {
              // add convective term
              // convection times gradient of ansatz fct.
              val += (u1_values[0] * gradx[ij] + u2_values[0] * grady[ij])
                * test[ij];
            }
            // multiplied with test fct.
            a_sub[index] += val*area/4.0;

          }                      // end jj
        }                        // end ii
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*for (jj=0;jj<N_nodes2;jj++)
    {
    for (ii=0;ii<N_nodes2;ii++)
    {
      if (fabs(a_sub[jj+ii*N_nodes2])<1e-10)
    a_sub[jj+ii*N_nodes2] = 0;
      OutPut(a_sub[jj+ii*N_nodes2] << " ");
    }
    OutPut(endl);
    }*/

    // just for comparison with a known shape of solution
    /*for (j=0;j< N_nodes2 ; j++)
      OutPut("rhs " << rhs_sub[j+2*N_nodes2] << endl); */

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    SolveMultipleSystems(a_sub, rhs_sub, N_nodes2, N_nodes2, N_nodes2 , 4);

    // just for comparison with a known shape of solution
    /*    for (j=0;j< N_nodes2 ; j++)
      OutPut(rhs_sub[j] << " " << rhs_sub[j+N_nodes2] << " " <<
      rhs_sub[j+2*N_nodes2] << " " << rhs_sub[j+3*N_nodes2]  << endl);*/
    /*	OutPut(" sol " << xlocal[j] << " " << ylocal[j] << " " <<
      rhs_sub[j+2*N_nodes2] << endl);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    // add contributions from the global function and the rhs
    //OutPut("norm vor " << sqrt(Ddot(4*N_nodes2,rhs_sub,rhs_sub)) );
    for (j=0;j< 2 *N_nodes2 ; j++)
      rhs_sub[j] += rhs_sub[j + 2*N_nodes2];
    //rhs_sub[j] = rhs_sub[j + 2*N_nodes2];

    //OutPut("norm nach " << sqrt(Ddot(2*N_nodes2,rhs_sub,rhs_sub)) << endl);

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase2D::GetFE2D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    test_value = new double[N_];
    integral_1 = new double[N_];
    integral_2 = new double[N_];
    memset(integral_1,0,N_*SizeOfDouble);
    memset(integral_2,0,N_*SizeOfDouble);

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        // compute coordinates of submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);

        xc[1] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[1] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[1] /= (N_sub*N_sub);
        yc[1] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[1] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[1] /= (N_sub*N_sub);

        xc[2] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[2] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[2] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[2] /= (N_sub*N_sub);

        xc[3] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[3] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[3] /=  (N_sub*N_sub);
        yc[3] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[3] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[3] /=  (N_sub*N_sub);

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        dof[0] = (ic-1)*N_nodes+jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + N_nodes;
        dof[3] = dof[0] + N_nodes;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[2] = -1;
          dof[3] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[1] = -1;
          dof[2] = -1;
        }
        // barycenter
        xs = xc[0] + xc[1] + xc[2] + xc[3];
        xs /= 4;
        ys = yc[0] + yc[1] + yc[2] + yc[3];
        ys /= 4;
        // compute quadrature points
        for (ij=0;ij<4;ij++)
        {
          x_quad[ij] = (xc[1]-xc[0])* gauss2_x[ij] + (xc[3]-xc[0])*gauss2_y[ij] + xc[1]+xc[3];
          x_quad[ij] /= 2;
          y_quad[ij] = (yc[1]-yc[0])* gauss2_x[ij] + (yc[3]-yc[0])*gauss2_y[ij] + yc[1]+yc[3];
          y_quad[ij] /= 2;
        }

        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);

        // compute the bilinear function on the mesh cell
        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<4;jj++)
        {
          a[4*jj] = 1;
          a[4*jj+1] = xc[jj];
          a[4*jj+2] = yc[jj];
          a[4*jj+3] = xc[jj]*yc[jj];
          if (dof[jj]!=-1)
          {
            b[jj] =  rhs_sub[dof[jj]];
            b[jj+4] = rhs_sub[dof[jj]+N_nodes2];
          }
          else
            b[jj] = b[jj+4] = 0;
        }

        // solve system for the coefficients of the bilinear function
        SolveMultipleSystems(a,b,4,4,4,2);
        /*for (jj=0;jj<2;jj++)
          OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
          b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/

        // loop over the quadrature points
        for(ij=0;ij<4;ij++)
        {
          // compute the value of the test functions in the quad point
          TFEDatabase2D::GetRefFromOrig(RefTrans, x_quad[ij], y_quad[ij], xi, eta);
          bf->GetDerivatives(D00, xi, eta, test_value);
          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the quad points
          rhs1_x = b[1] + b[3]*y_quad[ij];
          rhs1_y = b[2] + b[3]*x_quad[ij];
          rhs2_x = b[5] + b[7]*y_quad[ij];
          rhs2_y = b[6] + b[7]*x_quad[ij];

          val1 = u1_values[0] * rhs1_x + u2_values[0] * rhs1_y;
          val2 = u1_values[0] * rhs2_x + u2_values[0] * rhs2_y;

          // update the integrals for all test functions
          for (j=0;j<N_;j++)
          {
            integral_1[j] += val1 * test_value[j];
            integral_2[j] += val2 * test_value[j];
          }
        }                        // ij
        for (j=0;j<N_;j++)
        {
          integral_1[j] *= area/4;
          integral_2[j] *= area/4;
        }
      }                          // jc
    }                            // ic
    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
    }
    delete test_value;
    delete integral_1;
    delete integral_2;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;
  delete xlocal;
  delete ylocal;

}


// =======================================================================
// approximate solution with Q2 finite elements
// =======================================================================

void ApproximateRFBSolutionQuad_Q2_NSE2D(TCollection *Coll, TFEFunction2D *u1,
TFEFunction2D *u2, CoeffFct2D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, dof[9], N_U, N_, ij;
  int N_sub2, N_nodes, N_nodes2, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[4],y[4], a[81], b[81], *a_sub, *rhs_sub, test[8];
  double xc[9], yc[9], area, val, xs, ys, u1_values[3], u2_values[3];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff, *xlocal, *ylocal;
  double *integral_1, *integral_2, xi, eta, *test_value;
  double rhs1_x, rhs2_x, rhs1_y, rhs2_y;
  double gradx[8], grady[8], gradx_test[8], grady_test[8];
  double x_quad[4], y_quad[4], val_test, val1, val2;
  double gauss2_x[4]=
  {
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  double gauss2_y[4]=
  {
    -0.57735026918962576450914878050195746,
    -0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746,
    0.57735026918962576450914878050195746
  };
  TBaseCell *cell;
  TVertex *vertex;
  FE2D CurrentElement;
  TFESpace2D *fespace;
  TFE2D *FE_Obj;
  RefTrans2D RefTrans;
  TBaseFunct2D *bf;

  N_sub2 = N_sub * N_sub;
  N_nodes = 2*N_sub-1;
  N_nodes2 = N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes2*N_nodes2];
  memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);

  rhs_sub = new double[4 * N_nodes2];
  memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

  coeff = new double[3];

  xlocal = new double[N_nodes2];
  ylocal = new double[N_nodes2];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace2D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 4)
    {
      OutPut("RFB stabilization only for quadrilateral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j]);
      //OutPut("coords " << x[j] << " " << y[j] << endl);
    }

    // area of parallelogramm with vector product
    area = fabs((x[1]-x[0])*(y[3]-y[0]) - (x[3]-x[0])*(y[1]-y[0]));
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub2 ;
    //OutPut(" " << area << endl);
    memset(a_sub, 0, N_nodes2*N_nodes2*SizeOfDouble);
    memset(rhs_sub, 0, 4*N_nodes2*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        //OutPut("cell " << ic*N_sub+jc << " : ");
        // compute coordinates of submesh cell
        // first the vertices of the submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);
        //OutPut(xc[0] << " " << yc[0] << " ");
        xc[2] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[2] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[2] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[2] /= (N_sub*N_sub);
        //OutPut(xc[1] << " " << yc[1] << " ");
        xc[8] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[8] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[8] /= (N_sub*N_sub);
        yc[8] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[8] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[8] /= (N_sub*N_sub);
        //OutPut(xc[2] << " " << yc[2] << " ");
        xc[6] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[6] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[6] /=  (N_sub*N_sub);
        yc[6] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[6] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[6] /=  (N_sub*N_sub);
        //OutPut(xc[3] << " " << yc[3] << endl);
        // coordinates on the edges
        xc[1] = (xc[0]+xc[2])/2;
        yc[1] = (yc[0]+yc[2])/2;
        xc[3] = (xc[0]+xc[6])/2;
        yc[3] = (yc[0]+yc[6])/2;
        xc[5] = (xc[2]+xc[8])/2;
        yc[5] = (yc[2]+yc[8])/2;
        xc[7] = (xc[6]+xc[8])/2;
        yc[7] = (yc[6]+yc[8])/2;
        // coordinates in the barycenter
        xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
        yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;

        // for (jj=0;jj<9;jj++)
        //    OutPut(xc[jj] << " " << yc[jj] << endl);

        // initialize rhs
        memset(b,0,81*SizeOfDouble);

        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<9;jj++)
        {
          a[9*jj] = 1;
          a[9*jj+1] = xc[jj];
          a[9*jj+2] = yc[jj];
          a[9*jj+3] = xc[jj]*yc[jj];
          a[9*jj+4] = xc[jj]*xc[jj];
          a[9*jj+5] = yc[jj]*yc[jj];
          a[9*jj+6] = a[9*jj+4]*yc[jj];
          a[9*jj+7] = a[9*jj+5]*xc[jj];
          a[9*jj+8] = a[9*jj+7]*xc[jj];
          b[10*jj] = 1;
        }

        // solve system for the coefficients of the bilinear function
        SolveMultipleSystems(a,b,9,9,9,9);
        //for (jj=0;jj<9;jj++)
        //    OutPut(b[9*jj] << endl);
        //for (j=0;j<N_V;j++)
        //OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
        //b[4*j+2] << " " <<  b[4*j+3] << endl);

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        // enumeration is lexicographically
        // indizes correspond to indizes of xc, yc
        dof[0] = (2*ic-1)*N_nodes+2*jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + 1;
        dof[3] = dof[0] + N_nodes;
        dof[4] = dof[3] + 1;
        dof[5] = dof[4] + 1;
        dof[6] = dof[3] + N_nodes;
        dof[7] = dof[6] + 1;
        dof[8] = dof[7] + 1;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
          dof[2] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[6] = -1;
          dof[7] = -1;
          dof[8] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
          dof[6] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[2] = -1;
          dof[5] = -1;
          dof[8] = -1;
        }
        //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << endl);
        // just for graphical output
        /* for (j=0;j<9;j++)
        {
          if (dof[j] != -1)
          {
            xlocal[dof[j]] = xc[j];
            ylocal[dof[j]] = yc[j];
          }
        }*/
        //for (j=0;j<9;j++)
        //    OutPut("j " << j << " " << dof[j] << endl);
        // barycenter
        xs = xc[4];
        ys = yc[4];
        // compute quadrature points
        for (ij=0;ij<4;ij++)
        {
          x_quad[ij] = (xc[1]-xc[0])* gauss2_x[ij] + (xc[3]-xc[0])*gauss2_y[ij] + xc[1]+xc[3];
          x_quad[ij] /= 2;
          y_quad[ij] = (yc[1]-yc[0])* gauss2_x[ij] + (yc[3]-yc[0])*gauss2_y[ij] + yc[1]+yc[3];
          y_quad[ij] /= 2;
        }
        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);
        // just for comparison with a known shape of solution
        //u1_values[0] = 1;
        //u2_values[0] = 0;

        // compute value of rhs in (xs,ys)
        Coeffs(1, &xs, &ys, NULL, &coeff);
        //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);
        // just for comparison with a known shape of solution
        //coeff[1] = 1;

        // update matrix entries
        // the matrix is stored column wise
        // ii -- test function
        for(ii=0;ii<9;ii++)
        {
          if (dof[ii] == -1)
            continue;
          // loop over the quadrature points
          val_test = 0;
          for(ij=0;ij<4;ij++)
          {
            test[ij] = b[9*ii] +  b[9*ii+1] * x_quad[ij] + b[9*ii+2] * y_quad[ij]
              + b[9*ii+3] * x_quad[ij] * y_quad[ij] + b[9*ii+4] * x_quad[ij] * x_quad[ij]
              + b[9*ii+5] * y_quad[ij] * y_quad[ij] + b[9*ii+6] * x_quad[ij] * x_quad[ij] *y_quad[ij]
              + b[9*ii+7] * x_quad[ij] * y_quad[ij]  * y_quad[ij]
              + b[9*ii+8] * x_quad[ij] * x_quad[ij] * y_quad[ij] * y_quad[ij];
            val_test += test[ij];
            //OutPut("test bary " << test << endl);

            // x part of the gradient
            gradx_test[ij] = b[9*ii+1]+ b[9*ii+3]*y_quad[ij] + 2*b[9*ii+4]*x_quad[ij]
              + 2*b[9*ii+6]*x_quad[ij]*y_quad[ij]+ b[9*ii+7]*y_quad[ij]*y_quad[ij]
              + 2*b[9*ii+8]*x_quad[ij]*y_quad[ij]*y_quad[ij];
            // y part of the gradient
            grady_test[ij] = b[9*ii+2]+ b[9*ii+3]*x_quad[ij] + 2*b[9*ii+5]*y_quad[ij]
              + b[9*ii+6]*x_quad[ij]*x_quad[ij]+ 2*b[9*ii+7]*x_quad[ij]*y_quad[ij]
              + 2*b[9*ii+8]*x_quad[ij]*x_quad[ij]*y_quad[ij];
          }
          val_test *= area/4.0;

          // compute first rhs
          val = u1_values[0]*u1_values[1] +  u2_values[0]*u1_values[2];
          val *= val_test;
          rhs_sub[dof[ii]] -= val;
          // compute second rhs
          val = u1_values[0]*u2_values[1] +  u2_values[0]*u2_values[2];
          val *= val_test;
          rhs_sub[dof[ii]+ N_nodes2] -= val;

          // compute third and fourth rhs
          rhs_sub[dof[ii]+ 2*N_nodes2] += coeff[1]*val_test;
          rhs_sub[dof[ii]+ 3*N_nodes2] += coeff[2]*val_test;

          // jj -- ansatz function
          for (jj=0;jj<9;jj++)
          {
            if (dof[jj] == -1)
              continue;
            // entry (jj,ii) exists
            // compute index in array containing the entry

            index = dof[ii]*N_nodes2+dof[jj];
            //OutPut(" ii " << ii << " jj " << jj << " ind " << index);
            // add Laplacian
            // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
            // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
            // x part of the gradient
            val = 0;
            // loop over the quadrature points
            for(ij=0;ij<4;ij++)
            {
              gradx[ij] = b[9*jj+1]+ b[9*jj+3]*y_quad[ij] + 2*b[9*jj+4]*x_quad[ij]
                + 2*b[9*jj+6]*x_quad[ij]*y_quad[ij]+ b[9*jj+7]*y_quad[ij]*y_quad[ij]
                + 2*b[9*jj+8]*x_quad[ij]*y_quad[ij]*y_quad[ij];
              val += gradx[ij] * gradx_test[ij];
              // y part of the gradient
              grady[ij] = b[9*jj+2]+ b[9*jj+3]*x_quad[ij] + 2*b[9*jj+5]*y_quad[ij]
                + b[9*jj+6]*x_quad[ij]*x_quad[ij]+ 2*b[9*jj+7]*x_quad[ij]*y_quad[ij]
                + 2*b[9*jj+8]*x_quad[ij]*x_quad[ij]*y_quad[ij];
              val += grady[ij] * grady_test[ij];
            }
            val *= eps * area/4.0;
            a_sub[index] += val;

            // add convective term
            // convection times gradient of ansatz fct.
            // loop over the quadrature points
            val = 0;
            for(ij=0;ij<4;ij++)
            {
              // add convective term
              // convection times gradient of ansatz fct.
              val += (u1_values[0] * gradx[ij] + u2_values[0] * grady[ij])
                * test[ij];
            }
            // multiplied with test fct.
            a_sub[index] += val*area/4.0;

          }                      // end jj
        }                        // end ii
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*for (jj=0;jj<N_nodes2;jj++)
    {
    for (ii=0;ii<N_nodes2;ii++)
    {
      if (fabs(a_sub[jj+ii*N_nodes2])<1e-10)
    a_sub[jj+ii*N_nodes2] = 0;
      OutPut(a_sub[jj+ii*N_nodes2] << " ");
    }
    OutPut(endl);
    }*/

    // just for comparison with a known shape of solution
    //for (j=0;j< N_nodes2 ; j++)
    //  OutPut("rhs " << rhs_sub[j+2*N_nodes2] << endl);

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    SolveMultipleSystems(a_sub, rhs_sub, N_nodes2, N_nodes2, N_nodes2 , 4);

    // just for comparison with a known shape of solution
    /*    for (j=0;j< N_nodes2 ; j++)
      OutPut(rhs_sub[j] << " " << rhs_sub[j+N_nodes2] << " " <<
      rhs_sub[j+2*N_nodes2] << " " << rhs_sub[j+3*N_nodes2]  << endl);*/
    /*	OutPut(" sol " << xlocal[j] << " " << ylocal[j] << " " <<
      rhs_sub[j+2*N_nodes2] << endl);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    // add contributions from the global function and the rhs
    //OutPut("norm vor " << sqrt(Ddot(4*N_nodes2,rhs_sub,rhs_sub)) );
    for (j=0;j< 2 *N_nodes2 ; j++)
      rhs_sub[j] += rhs_sub[j + 2*N_nodes2];
    //rhs_sub[j] = rhs_sub[j + 2*N_nodes2];

    //OutPut("norm nach " << sqrt(Ddot(2*N_nodes2,rhs_sub,rhs_sub)) << endl);

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE2D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase2D::GetFE2D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct2D();
    test_value = new double[N_];
    integral_1 = new double[N_];
    integral_2 = new double[N_];
    memset(integral_1,0,N_*SizeOfDouble);
    memset(integral_2,0,N_*SizeOfDouble);

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        // compute coordinates of submesh cell
        // first the vertices of the submesh cell
        xc[0] = (N_sub-jc)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[0] += jc*((N_sub-ic) * x[1] + ic * x[2]);
        xc[0] /= (N_sub*N_sub);
        yc[0] = (N_sub-jc)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[0] += jc*((N_sub-ic) * y[1] + ic * y[2]);
        yc[0] /= (N_sub*N_sub);
        //OutPut(xc[0] << " " << yc[0] << " ");
        xc[2] = (N_sub-jc-1)*((N_sub-ic) * x[0] + ic * x[3]);
        xc[2] += (jc+1)*((N_sub-ic) * x[1] + ic * x[2]);
        xc[2] /= (N_sub*N_sub);
        yc[2] = (N_sub-jc-1)*((N_sub-ic) * y[0] + ic * y[3]);
        yc[2] += (jc+1)*((N_sub-ic) * y[1] + ic * y[2]);
        yc[2] /= (N_sub*N_sub);
        //OutPut(xc[1] << " " << yc[1] << " ");
        xc[8] = (N_sub-jc-1)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[8] += (jc+1)*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[8] /= (N_sub*N_sub);
        yc[8] = (N_sub-jc-1)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[8] += (jc+1)*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[8] /= (N_sub*N_sub);
        //OutPut(xc[2] << " " << yc[2] << " ");
        xc[6] = (N_sub-jc)*((N_sub-ic-1)* x[0] + (ic+1)* x[3]);
        xc[6] += jc*((N_sub-ic-1)* x[1] + (ic+1)* x[2]);
        xc[6] /=  (N_sub*N_sub);
        yc[6] = (N_sub-jc)*((N_sub-ic-1)* y[0] + (ic+1)* y[3]);
        yc[6] += jc*((N_sub-ic-1)* y[1] + (ic+1)* y[2]);
        yc[6] /=  (N_sub*N_sub);
        //OutPut(xc[3] << " " << yc[3] << endl);
        // coordinates on the edges
        xc[1] = (xc[0]+xc[2])/2;
        yc[1] = (yc[0]+yc[2])/2;
        xc[3] = (xc[0]+xc[6])/2;
        yc[3] = (yc[0]+yc[6])/2;
        xc[5] = (xc[2]+xc[8])/2;
        yc[5] = (yc[2]+yc[8])/2;
        xc[7] = (xc[6]+xc[8])/2;
        yc[7] = (yc[6]+yc[8])/2;
        // coordinates in the barycenter
        xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
        yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;

        // compute local dof at this mesh cell
        // compute dof of the left lower vertex
        // enumeration is lexicographically
        // indizes correspond to indizes of xc, yc
        dof[0] = (2*ic-1)*N_nodes+2*jc-1;
        dof[1] = dof[0] + 1;
        dof[2] = dof[1] + 1;
        dof[3] = dof[0] + N_nodes;
        dof[4] = dof[3] + 1;
        dof[5] = dof[4] + 1;
        dof[6] = dof[3] + N_nodes;
        dof[7] = dof[6] + 1;
        dof[8] = dof[7] + 1;
        // correct for submesh cells with vertices on boundary
        // lower boundary
        if (ic == 0)
        {
          dof[0] = -1;
          dof[1] = -1;
          dof[2] = -1;
        }
        // upper boundary
        if (ic==N_sub-1)
        {
          dof[6] = -1;
          dof[7] = -1;
          dof[8] = -1;
        }
        // left boundary
        if (jc == 0)
        {
          dof[0] = -1;
          dof[3] = -1;
          dof[6] = -1;
        }
        // right boundary
        if (jc==N_sub-1)
        {
          dof[2] = -1;
          dof[5] = -1;
          dof[8] = -1;
        }
        // barycenter
        xs = xc[4];
        ys = yc[4];

        // compute convection field in the center of the submesh cell
        u1->FindGradientLocal(cell,i,xs,ys,u1_values);
        u2->FindGradientLocal(cell,i,xs,ys,u2_values);

        // compute quadrature points
        for (ij=0;ij<4;ij++)
        {
          x_quad[ij] = (xc[1]-xc[0])* gauss2_x[ij] + (xc[3]-xc[0])*gauss2_y[ij] + xc[1]+xc[3];
          x_quad[ij] /= 2;
          y_quad[ij] = (yc[1]-yc[0])* gauss2_x[ij] + (yc[3]-yc[0])*gauss2_y[ij] + yc[1]+yc[3];
          y_quad[ij] /= 2;
        }

        // compute the bilinear function on the mesh cell
        // set matrices for computation of the coefficients
        // of the bilinear function
        for (jj=0;jj<9;jj++)
        {
          a[9*jj] = 1;
          a[9*jj+1] = xc[jj];
          a[9*jj+2] = yc[jj];
          a[9*jj+3] = xc[jj]*yc[jj];
          a[9*jj+4] = xc[jj]*xc[jj];
          a[9*jj+5] = yc[jj]*yc[jj];
          a[9*jj+6] = a[9*jj+4]*yc[jj];
          a[9*jj+7] = a[9*jj+5]*xc[jj];
          a[9*jj+8] = a[9*jj+7]*xc[jj];
          if (dof[jj]!=-1)
          {
            b[jj] =  rhs_sub[dof[jj]];
            b[jj+9] = rhs_sub[dof[jj]+N_nodes2];
          }
          else
            b[jj] = b[jj+9] = 0;
        }

        // solve system for the coefficients of the biquadratic function
        SolveMultipleSystems(a,b,9,9,9,2);

        // loop over the quadrature points
        for(ij=0;ij<4;ij++)
        {
          // compute the value of the test functions in xc[ij], yc[ij]
          TFEDatabase2D::GetRefFromOrig(RefTrans, xc[ij], yc[ij], xi, eta);
          bf->GetDerivatives(D00, xi, eta, test_value);
          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the quad point

          rhs1_x = b[1]+ b[3]*y_quad[ij] + 2*b[4]*x_quad[ij] + 2*b[6]*x_quad[ij]*y_quad[ij]
            + b[7]*y_quad[ij]*y_quad[ij] + 2*b[8]*x_quad[ij]*y_quad[ij]*y_quad[ij];
          rhs1_y = b[2]+ b[3]*x_quad[ij] + 2*b[5]*y_quad[ij] + b[6]*x_quad[ij]*x_quad[ij]
            + 2*b[7]*x_quad[ij]*y_quad[ij] + 2*b[8]*x_quad[ij]*x_quad[ij]*y_quad[ij];
          rhs2_x = b[10]+ b[12]*y_quad[ij] + 2*b[13]*x_quad[ij] + 2*b[15]*x_quad[ij]*y_quad[ij]
            + b[16]*y_quad[ij]*y_quad[ij] + 2*b[17]*x_quad[ij]*y_quad[ij]*y_quad[ij];
          rhs2_y = b[11]+ b[12]*x_quad[ij] + 2*b[14]*y_quad[ij] + b[15]*x_quad[ij]*x_quad[ij]
            + 2*b[16]*x_quad[ij]*y_quad[ij] + 2*b[17]*x_quad[ij]*x_quad[ij]*y_quad[ij];

          val1 = u1_values[0] * rhs1_x + u2_values[0] * rhs1_y;
          val2 = u1_values[0] * rhs2_x + u2_values[0] * rhs2_y;
          if (ij==4)
          {
            val1 *= 2;
            val2 *= 2;
          }
          // update the integrals
          for (j=0;j<N_;j++)
          {
            integral_1[j] += val1 * test_value[j];
            integral_2[j] += val2 * test_value[j];
          }
        }                        // ij
        for (j=0;j<N_;j++)
        {
          integral_1[j] *= area/6;
          integral_2[j] *= area/6;
        }

      }                          // jc
    }                            // ic
    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
    }
    delete test_value;
    delete integral_1;
    delete integral_2;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;
  delete xlocal;
  delete ylocal;
}
#endif

#ifdef __3D__
void VertexCoordinatesSubMeshCell(int N_sub, int ic, int jc, int kc,
double *x, double *y, double *z,
double *xc, double *yc, double *zc)
{
  double val, N_sub3;

  N_sub3 = N_sub*N_sub*N_sub;

  val = kc*(jc*(ic * x[6] + (N_sub-ic)*x[2]));
  val += kc*((N_sub-jc)*(ic * x[5] + (N_sub-ic)*x[1]));
  val += (N_sub-kc)*(jc*(ic * x[7] + (N_sub-ic)*x[3]));
  val += (N_sub-kc)*((N_sub-jc)*(ic * x[4] + (N_sub-ic)*x[0]));
  xc[0] = val/N_sub3;

  val = kc*(jc*(ic * y[6] + (N_sub-ic)*y[2]));
  val += kc*((N_sub-jc)*(ic * y[5] + (N_sub-ic)*y[1]));
  val += (N_sub-kc)*(jc*(ic * y[7] + (N_sub-ic)*y[3]));
  val += (N_sub-kc)*((N_sub-jc)*(ic * y[4] + (N_sub-ic)*y[0]));
  yc[0] = val/N_sub3;

  val = kc*(jc*(ic * z[6] + (N_sub-ic)*z[2]));
  val += kc*((N_sub-jc)*(ic * z[5] + (N_sub-ic)*z[1]));
  val += (N_sub-kc)*(jc*(ic * z[7] + (N_sub-ic)*z[3]));
  val += (N_sub-kc)*((N_sub-jc)*(ic * z[4] + (N_sub-ic)*z[0]));
  zc[0] = val/N_sub3;

  return;
}


void Compute_Q1_Value_Gradient_RFB(double *coeff, double x, double y, double z, double *val)
{
  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + coeff[4] * x * y +  coeff[5] * x * z +  coeff[6] * y * z
    + coeff[7] * x * y *z;
  // x derivative
  val[1] = coeff[1]+ coeff[4]*y + coeff[5]*z + coeff[7]*y*z;
  // y derivative
  val[2] = coeff[2]+ coeff[4]*x + coeff[6]*z + coeff[7]*x*z;
  // z derivative
  val[3] = coeff[3]+ coeff[5]*x + coeff[6]*y + coeff[7]*x*y;
}

void Compute_Q1_Value_Gradient_Center_From_Vertices(double *values, double *val)
{
  // value
  val[0] = values[0] +  values[1] + values[2] + values[3] 
      + values[4] +  values[5] + values[6] + values[7];
  val[0] /= 8.0;
  val[1] = val[2] = val[3] = 0.0;
}


void Compute_Q1_Value_RFB(double *coeff, double x, double y, double z, double *val)
{
  // value
  *val = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + coeff[4] * x * y +  coeff[5] * x * z +  coeff[6] * y * z
    + coeff[7] * x * y *z;
}


void Compute_Q2_Value(double *coeff, double x, double y, double z, double *val)
{
  double xy, xz, xyz, yz;

  xy=x*y;
  xz=x*z;
  xyz=x*y*z;
  yz=y*z;

  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + coeff[4] * xy +  coeff[5] * xz +  coeff[6] * yz
    + coeff[7] * xyz + coeff[8] * x * x + coeff[9] * y * y + coeff[10] * z * z
    + coeff[11] * x * xy + coeff[12] * x * xz + coeff[13] * xy * y
    + coeff[14] * y * yz + coeff[15] * xz * z + coeff[16] * yz * z
    + coeff[17] * xy * xz + coeff[18] * xy * yz
    + coeff[19] * xyz * z + coeff[20] * xy * xy + coeff[21] * xz * xz
    + coeff[22] * yz * yz + coeff[23] * xy * xyz
    + coeff[24] * xyz* xz + coeff[25] * xyz * yz
    + coeff[26] * xyz * xyz;
}


void Compute_Q2_Value_Gradient(double *coeff, double x, double y, double z, double *val)
{
  double xx, xy, xz, yy, xyz, yz, zz;

  xx=x*x;
  xy=x*y;
  xz=x*z;
  yy=y*y;
  xyz=x*y*z;
  yz=y*z;
  zz=z*z;

  // value
  val[0] = coeff[0] +  coeff[1] * x + coeff[2] * y + coeff[3] * z
    + coeff[4] * xy +  coeff[5] * xz +  coeff[6] * yz
    + coeff[7] * xyz + coeff[8] * x * x + coeff[9] * y * y + coeff[10] * z * z
    + coeff[11] * x * xy + coeff[12] * x * xz + coeff[13] * xy * y
    + coeff[14] * y * yz + coeff[15] * xz * z + coeff[16] * yz * z
    + coeff[17] * xy * xz + coeff[18] * xy * yz
    + coeff[19] * xz * yz + coeff[20] * xy * xy + coeff[21] * xz * xz
    + coeff[22] * yz * yz + coeff[23] * xy * xyz
    + coeff[24] * xyz * xz + coeff[25] * xyz * yz
    + coeff[26] * xyz*xyz;
  // x derivative
  val[1] = coeff[1]+ coeff[4]*y + coeff[5]*z + coeff[7]*yz + 2*coeff[8]*x
    + 2*coeff[11]*xy + 2*coeff[12]*xz + coeff[13]*yy + coeff[15]*zz
    + 2*coeff[17]*xyz + coeff[18]*yy*z + coeff[19]*yz*z
    + 2*coeff[20]*xy*y + 2*coeff[21]*xz*z + 2*coeff[23]*xy*yz
    + 2*coeff[24]*xy*zz+ coeff[25]*yy*zz + 2*coeff[26]*xyz*yz;
  // y derivative
  val[2] = coeff[2]+ coeff[4]*x + coeff[6]*z + coeff[7]*x*z + 2*coeff[9]*y
    + coeff[11]*xx + 2*coeff[13]*xy + 2*coeff[14]*yz + coeff[16]*zz
    + coeff[17]*xx*z + 2*coeff[18]*xyz + coeff[19]*xz*z
    + 2*coeff[20]*xx*y + 2*coeff[22]*yz*z + 2*coeff[23]*xx*yz
    + coeff[24]*xx*zz+ 2*coeff[25]*xy*zz + 2*coeff[26]*xyz*xz;
  // z derivative
  val[3] = coeff[3]+ coeff[5]*x + coeff[6]*y + coeff[7]*x*y + 2*coeff[10]*z
    + coeff[12]*xx + coeff[14]*yy + 2*coeff[15]*xz + 2*coeff[16]*yz
    + coeff[17]*xx*y + coeff[18]*xy*y + 2*coeff[19]*xy*z
    + 2*coeff[21]*xx*z + 2*coeff[22]*yy*z + coeff[23]*xx*yy
    + 2*coeff[24]*xx*yz+ 2*coeff[25]*xy*yz + 2*coeff[26]*xyz*xy;
}


void ApproximateRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, CoeffFct3D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, dof[8], N_U, N_;
  int N_sub3, N_nodes, N_nodes3, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[8],y[8], z[8], a[64], b[64], *a_sub, *rhs_sub, test_bary;
  double xc[8], yc[8], zc[8], area, val, xs, ys, zs, u1_values[4];
  double u2_values[4], u3_values[4];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value;
  double rhs1_x, rhs2_x, rhs3_x, rhs1_y, rhs2_y, rhs3_y, rhs1_z, rhs2_z, rhs3_z;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;

  N_sub3 = N_sub * N_sub * N_sub;
  N_nodes = N_sub-1;
  N_nodes3 = N_nodes*N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes3*N_nodes3];
  memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);

  rhs_sub = new double[6 * N_nodes3];
  memset(rhs_sub, 0, 6*N_nodes3*SizeOfDouble);

  coeff = new double[4];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }

    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);
    memset(rhs_sub, 0, 6*N_nodes3*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          //for (jj=0;jj<8;jj++)
          //  OutPut(jj << " " << xc[jj]  << " " << yc[jj]  << " " << zc[jj] << endl);

          // initialize rhs
          memset(b,0,64*SizeOfDouble);

          // set matrices for computation of the coefficients
          // of the trilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            b[9*jj] = 1;
          }
          //for (j=0;j<16;j++)
          //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

          // solve system for the coefficients of the bilinear function
          SolveMultipleSystems(a,b,8,8,8,8);
          /*for (j=0;j<N_V;j++)
          OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
          b[4*j+2] << " " <<  b[4*j+3] << endl);
          */

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }
          //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << " ");
          //OutPut(dof[4] << " " << dof[5] << " " << dof[6] << " " << dof[7] << endl);

          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // compute value of rhs in (xs,ys)
          Coeffs(1, &xs, &ys, &zs, NULL, &coeff);
          //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);

          // update matrix entries
          // the matrix is stored column wise
          // ii -- test function
          for(ii=0;ii<8;ii++)
          {
            if (dof[ii] == -1)
              continue;
            // value of test function in bary center
            test_bary = b[8*ii] +  b[8*ii+1] * xs + b[8*ii+2] * ys + b[8*ii+3] * zs
              + b[8*ii+4] * xs * ys +  b[8*ii+5] * xs * zs +  b[8*ii+6] * ys * zs
              + b[8*ii+7] * xs * ys *zs;
            //OutPut("test bary " << test_bary << endl);

            // compute first rhs
            val = u1_values[0]*u1_values[1] +  u2_values[0]*u1_values[2] +  u3_values[0]*u1_values[3];
            val *= test_bary;
            rhs_sub[dof[ii]] -= area*val;
            // compute second rhs
            val = u1_values[0]*u2_values[1] +  u2_values[0]*u2_values[2]  + u3_values[0]*u2_values[3];
            val *= test_bary;
            rhs_sub[dof[ii]+ N_nodes3] -= area*val;
            // compute third rhs
            val = u1_values[0]*u3_values[1] +  u2_values[0]*u3_values[2]  + u3_values[0]*u3_values[3];
            val *= test_bary;
            rhs_sub[dof[ii]+ 2*N_nodes3] -= area*val;

            // compute third and fourth rhs
            rhs_sub[dof[ii]+ 3*N_nodes3] += area*coeff[1]*test_bary;
            rhs_sub[dof[ii]+ 4*N_nodes3] += area*coeff[2]*test_bary;
            rhs_sub[dof[ii]+ 5*N_nodes3] += area*coeff[3]*test_bary;

            // jj -- ansatz function
            for (jj=0;jj<8;jj++)
            {
              if (dof[jj] == -1)
                continue;
              // entry (jj,ii) exists
              // compute index in array containing the entry

              index = dof[ii]*N_nodes3+dof[jj];
              //OutPut(" ii " << ii << " jj " << jj << " ind " << index )
              // add Laplacian
              // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
              // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
              // approximate integral by midpoint rule
              // x part of the gradient
              val = (b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs)*
                (b[8*ii+1]+ b[8*ii+4]*ys + b[8*ii+5]*zs + b[8*ii+7]*ys*zs);
              // y part of the gradient
              val += (b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs)*
                (b[8*ii+2]+ b[8*ii+4]*xs + b[8*ii+6]*zs + b[8*ii+7]*xs*zs);
              // z part of the gradient
              val += (b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys)*
                (b[8*ii+3]+ b[8*ii+5]*xs + b[8*ii+6]*ys + b[8*ii+7]*xs*ys);
              val *= eps * area;
              a_sub[index] += val;

              // add convective term
              // convection times gradient of ansatz fct.
              val = u1_values[0] * (b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs)
                + u2_values[0] * (b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs)
                + u3_values[0] * (b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys);
              // multiplied with test fct.
              val *= test_bary;
              val *= area;
              a_sub[index] += val;

            }                    // end jj
          }                      // end ii
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*for (jj=0;jj<N_nodes3;jj++)
    {
    for (ii=0;ii<N_nodes3;ii++)
    {
      if (fabs(a_sub[jj+ii*N_nodes3])<1e-10)
    a_sub[jj+ii*N_nodes3] = 0;
      OutPut(a_sub[jj+ii*N_nodes3] << " ");
    }
    OutPut(endl);
    }*/

    // just for comparison with a known shape of solution
    /*for (j=0;j< N_nodes3 ; j++)
      OutPut("rhs " << rhs_sub[j+2*N_nodes3] << endl); */

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    SolveMultipleSystems(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 6);

    // just for comparison with a known shape of solution
    /*    for (j=0;j< N_nodes3 ; j++)
      OutPut(rhs_sub[j] << " " << rhs_sub[j+N_nodes3] << " " <<
      rhs_sub[j+2*N_nodes3] << " " << rhs_sub[j+3*N_nodes3]  << endl);*/
    /*	OutPut(" sol " << xlocal[j] << " " << ylocal[j] << " " <<
      rhs_sub[j+2*N_nodes3] << endl);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    // add contributions from the global function and the rhs
    for (j=0;j< 3 *N_nodes3 ; j++)
      rhs_sub[j] += rhs_sub[j + 3*N_nodes3];

    //OutPut(" rhs " << sqrt(Ddot(3*N_nodes3,rhs_sub,rhs_sub)));

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();
    test_value = new double[N_];
    integral_1 = new double[N_];
    integral_2 = new double[N_];
    integral_3 = new double[N_];
    memset(integral_1,0,N_*SizeOfDouble);
    memset(integral_2,0,N_*SizeOfDouble);
    memset(integral_3,0,N_*SizeOfDouble);

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }

          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // compute the value of the test functions in the bary center
          // of the submesh cell
          TFEDatabase3D::GetRefFromOrig(RefTrans, xs, ys, zs, xi, eta, zeta);
          bf->GetDerivatives(D000, xi, eta, zeta, test_value);
          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the barycenter of the submesh cell
          // compute for this purpose the bilinear function on the mesh cell

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            if (dof[jj]!=-1)
            {
              b[jj] =  rhs_sub[dof[jj]];
              b[jj+8] = rhs_sub[dof[jj]+N_nodes3];
              b[jj+16] = rhs_sub[dof[jj]+2*N_nodes3];
            }
            else
              b[jj] = b[jj+8] = b[jj+16] = 0;
          }

          // solve system for the coefficients of the bilinear function
          SolveMultipleSystems(a,b,8,8,8,3);
          /*for (jj=0;jj<2;jj++)
            OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
            b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/
          // compute gradients of rhs
          rhs1_x = b[1] + b[4]*ys+ b[5]*zs + b[7]*ys*zs;
          rhs1_y = b[2] + b[4]*xs+ b[6]*zs + b[7]*xs*zs;
          rhs1_z = b[3] + b[5]*xs+ b[6]*ys + b[7]*xs*ys;
          rhs2_x = b[9] + b[12]*ys+ b[13]*zs + b[15]*ys*zs;
          rhs2_y = b[10] + b[12]*xs+ b[14]*zs + b[15]*xs*zs;
          rhs2_z = b[11] + b[13]*xs+ b[14]*ys + b[15]*xs*ys;
          rhs3_x = b[17] + b[20]*ys+ b[21]*zs + b[23]*ys*zs;
          rhs3_y = b[18] + b[20]*xs+ b[22]*zs + b[23]*xs*zs;
          rhs3_z = b[19] + b[21]*xs+ b[22]*ys + b[23]*xs*ys;

          // update the integrals
          for (j=0;j<N_;j++)
          {
            val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y + u3_values[0] * rhs1_z)*test_value[j]*area;
            integral_1[j] += val;
            val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y + u3_values[0] * rhs2_z)*test_value[j]*area;
            integral_2[j] += val;
            val = (u1_values[0] * rhs3_x + u2_values[0] * rhs3_y + u3_values[0] * rhs3_z)*test_value[j]*area;
            integral_3[j] += val;
          }
        }                        // kc
      }                          // jc
    }                            // ic
    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
      rhs[gdof+2*N_U] -= integral_3[j];
    }
    delete test_value;
    delete integral_1;
    delete integral_2;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;

}


// =======================================================================
// RFB for time dependent NSE
// using approach of V.Gravemeier
// Q1 fe
// =======================================================================

void ApproximateTimeRFBSolutionQuadNSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, dof[8], N_U, N_;
  int N_sub3, N_nodes, N_nodes3, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof, dof_ii, dof_ii2;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[8],y[8], z[8], a[64], b[64], *a_sub, *rhs_sub, test_bary;
  double xc[8], yc[8], zc[8], area, val, val1, val2, val3, xs, ys, zs, u1_values[4];
  double u2_values[4], u3_values[4], p_values[4], val1_j, val2_j, val3_j;
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value, *test_value_x, *test_value_y, *test_value_z;
  double rhs1, rhs2, rhs3, rhs1_x, rhs2_x, rhs3_x, rhs1_y, rhs2_y, rhs3_y, rhs1_z, rhs2_z, rhs3_z;
  double valU[3], gradU[9], delta, hK, mu_T, tau_K, divU;
  double c_1 = 2, c_2 = 1;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double value_global[3];
  value_global[0] = value_global[1] = value_global[2] = 0;

  OutPut("rfb"<<endl);

  N_sub3 = N_sub * N_sub * N_sub;
  N_nodes = N_sub-1;
  N_nodes3 = N_nodes*N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes3*N_nodes3];
  memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);

  rhs_sub = new double[9 * N_nodes3];
  memset(rhs_sub, 0, 9*N_nodes3*SizeOfDouble);

  coeff = new double[4*N_sub];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }

    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    // compute delta for stabilization

    if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)
      hK = cell->GetShortestEdge();
    else
      hK = cell->GetDiameter();
    // measure of the submesh cells
    hK /= N_sub;
    delta = CharacteristicFilterWidth(hK);

    memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);
    memset(rhs_sub, 0, 3*N_nodes3*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          //for (jj=0;jj<8;jj++)
          //  OutPut(jj << " " << xc[jj]  << " " << yc[jj]  << " " << zc[jj] << endl);

          // initialize rhs
          memset(b,0,64*SizeOfDouble);

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            b[9*jj] = 1;
          }
          //for (j=0;j<16;j++)
          //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

          // solve system for the coefficients of the bilinear function
          //SolveMultipleSystems(a,b,8,8,8,8);
          SolveMultipleSystemsLapack(a,b,8,8,8,8);
          //for (j=0;j<N_V;j++)
          //OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
          //b[4*j+2] << " " <<  b[4*j+3] << endl);
	 
          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }
          //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << " ");
          //OutPut(dof[4] << " " << dof[5] << " " << dof[6] << " " << dof[7] << endl);

          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // compute pressure in the center of the submesh cell
          p->FindGradientLocal(cell,i,xs,ys,zs,p_values);

          // compute value of rhs in (xs,ys)
          Coeffs(1, &xs, &ys, &zs, NULL, &coeff);
          //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);

          // prepare computation of turbulent viscosity
          valU[0]=u1_values[0];
          valU[1]=u2_values[0];
          valU[2]=u3_values[0];
          for(ii=0;ii<3;ii++)
          {
            gradU[3*ii]=u1_values[ii+1];
            gradU[3*ii+1]=u2_values[ii+1];
            gradU[3*ii+2]=u3_values[ii+1];
          }
          // turbulent viscosity
          mu_T = TurbulentViscosity3D(delta,gradU,valU,NULL,NULL,NULL,&zs,-4711);

          // parameter for div-div term
          tau_K = DivDivStab3D(u1_values[0],u2_values[0],u3_values[0],hK,eps);
          divU = u1_values[1] + u2_values[2] + u3_values[3];

          // update matrix entries
          // the matrix is stored column wise
          // ii -- test function
          for(ii=0;ii<8;ii++)
          {
            dof_ii = dof[ii];
            if (dof_ii == -1)
              continue;
            // value of test function in bary center
            test_bary = b[8*ii] +  b[8*ii+1] * xs + b[8*ii+2] * ys + b[8*ii+3] * zs
              + b[8*ii+4] * xs * ys +  b[8*ii+5] * xs * zs +  b[8*ii+6] * ys * zs
              + b[8*ii+7] * xs * ys *zs;
            //OutPut("test bary " << test_bary << endl);

            // values of the derivatives of the test function in bary center
            val1 = b[8*ii+1]+ b[8*ii+4]*ys + b[8*ii+5]*zs + b[8*ii+7]*ys*zs;
            val2 = b[8*ii+2]+ b[8*ii+4]*xs + b[8*ii+6]*zs + b[8*ii+7]*xs*zs;
            val3 = b[8*ii+3]+ b[8*ii+5]*xs + b[8*ii+6]*ys + b[8*ii+7]*xs*ys;

            // compute first rhs
            // convection term
            val = -u1_values[0]*u1_values[1] -  u2_values[0]*u1_values[2] -  u3_values[0]*u1_values[3];
            val *= test_bary;
            // diffusion term
            val -= eps*(u1_values[1]*val1+u1_values[2]*val2+u1_values[3]*val3);
            // div-div term
            val -= tau_K * divU * val1;
            // pressure term
            val += p_values[0]*val1;
            // rhs term
            val += coeff[1]*test_bary;
            rhs_sub[dof_ii] += tau*area*val;

            // compute second rhs
            dof_ii2 = dof_ii + N_nodes3;
            val = -u1_values[0]*u2_values[1] -  u2_values[0]*u2_values[2]  - u3_values[0]*u2_values[3];
            val *= test_bary;
            val -= eps *(u2_values[1]*val1+u2_values[2]*val2+u2_values[3]*val3);
            val -= tau_K * divU * val2;
            val += p_values[0]*val2;
            val += coeff[2]*test_bary;
            rhs_sub[dof_ii2] += tau*area*val;

            // compute third rhs
            dof_ii2 = dof_ii + 2*N_nodes3;
            val = -u1_values[0]*u3_values[1] -  u2_values[0]*u3_values[2]  - u3_values[0]*u3_values[3];
            val *= test_bary;
            val -= eps*(u3_values[1]*val1+u3_values[2]*val2+u3_values[3]*val3);
            val -= tau_K * divU * val3;
            val += p_values[0]*val3;
            val += coeff[3]*test_bary;
            rhs_sub[dof_ii2] += tau*area*val;

            // jj -- ansatz function
            for (jj=0;jj<8;jj++)
            {
              if (dof[jj] == -1)
                continue;
              // entry (jj,ii) exists
              // compute index in array containing the entry

              index = dof_ii*N_nodes3+dof[jj];
              //OutPut(" ii " << ii << " jj " << jj << " ind " << index )
              // add Laplacian
              // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
              // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
              // approximate integral by midpoint rule

              val1_j = b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs;
              val2_j = b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs;
              val3_j = b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys;

              val = val1_j * val1 + val2_j * val2 + val3_j * val3;
              val *= tau*(eps+mu_T) * area;
              a_sub[index] += val;

              // add convective term
              // convection times gradient of ansatz fct.
              val = u1_values[0] * val1_j + u2_values[0] * val2_j + u3_values[0] * val3_j;
              // multiplied with test fct.
              val *= tau*test_bary*area;
              a_sub[index] += val;

              // add time-derivative terms
              val = (b[8*jj] +  b[8*jj+1] * xs + b[8*jj+2] * ys + b[8*jj+3] * zs
                + b[8*jj+4] * xs * ys +  b[8*jj+5] * xs * zs +  b[8*jj+6] * ys * zs
                + b[8*jj+7] * xs * ys *zs) * test_bary;
              val *= area;
              a_sub[index] += val;
            }                    // end jj
          }                      // end ii
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    //SolveMultipleSystems(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 3);
    SolveMultipleSystemsLapack(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 3);
   
    //for (i=0;i<3*N_nodes3;i++)
//	OutPut(i << " " << rhs_sub[i] << endl);
     /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();
    test_value = new double[4*N_];
    test_value_x = test_value + N_;
    test_value_y = test_value_x + N_;
    test_value_z = test_value_y + N_;
    integral_1 = new double[3 * N_];
    memset(integral_1,0,3*N_*SizeOfDouble);
    integral_2 = integral_1 + N_;
    integral_3 = integral_2 + N_;

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }
          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // compute the value of the test functions and its derivatives
          //in the bary center of the submesh cell
          TFEDatabase3D::GetRefFromOrig(RefTrans, xs, ys, zs, xi, eta, zeta);
          bf->GetDerivatives(D000, xi, eta, zeta, test_value);
          bf->GetDerivatives(D100, xi, eta, zeta, test_value_x);
          bf->GetDerivatives(D010, xi, eta, zeta, test_value_y);
          bf->GetDerivatives(D001, xi, eta, zeta, test_value_z);

          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the barycenter of the submesh cell
          // compute for this purpose the bilinear function on the mesh cell

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            if (dof[jj]!=-1)
            {
              b[jj] =  rhs_sub[dof[jj]];
              b[jj+8] = rhs_sub[dof[jj]+N_nodes3];
              b[jj+16] = rhs_sub[dof[jj]+2*N_nodes3];
            }
            else
              b[jj] = b[jj+8] = b[jj+16] = 0;
          }

          // solve system for the coefficients of the bilinear function
          //SolveMultipleSystems(a,b,8,8,8,3);
          SolveMultipleSystemsLapack(a,b,8,8,8,3);
          /*for (jj=0;jj<2;jj++)
            OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
            b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/

          // compute values of tilde_u in the barycenter
          rhs1 = b[0] + b[1]*xs + b[2]*ys + b[3]*zs
            + b[4]*xs*ys + b[5]*xs*zs + b[6]*ys*zs + b[7]*xs*ys*zs;
          rhs2 = b[8] + b[9]*xs + b[10]*ys + b[11]*zs
            + b[12]*xs*ys + b[13]*xs*zs + b[14]*ys*zs + b[15]*xs*ys*zs;
          rhs3 = b[16] + b[17]*xs + b[18]*ys + b[19]*zs
            + b[20]*xs*ys + b[21]*xs*zs + b[22]*ys*zs + b[23]*xs*ys*zs;

          // compute gradients of tilde_u in the barycenter
          rhs1_x = b[1] + b[4]*ys+ b[5]*zs + b[7]*ys*zs;
          rhs1_y = b[2] + b[4]*xs+ b[6]*zs + b[7]*xs*zs;
          rhs1_z = b[3] + b[5]*xs+ b[6]*ys + b[7]*xs*ys;
          rhs2_x = b[9] + b[12]*ys+ b[13]*zs + b[15]*ys*zs;
          rhs2_y = b[10] + b[12]*xs+ b[14]*zs + b[15]*xs*zs;
          rhs2_z = b[11] + b[13]*xs+ b[14]*ys + b[15]*xs*ys;
          rhs3_x = b[17] + b[20]*ys+ b[21]*zs + b[23]*ys*zs;
          rhs3_y = b[18] + b[20]*xs+ b[22]*zs + b[23]*xs*zs;
          rhs3_z = b[19] + b[21]*xs+ b[22]*ys + b[23]*xs*ys;

          // update the integrals
          for (j=0;j<N_;j++)
          {
	      // first component   
            // add convection term
            val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y + u3_values[0] * rhs1_z )
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * rhs1_x + rhs2 * rhs1_y + rhs3 * rhs1_z )
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * u1_values[1] + rhs2 * u1_values[2] + rhs3 * u1_values[3] )
		*test_value[j]*area*theta1*tau;
            //add laplacian
            val += (rhs1_x*test_value_x[j] + rhs1_y*test_value_y[j] + rhs1_z*test_value_z[j]) *area*eps*theta1*tau;
            // add term from time derivative
            val += rhs1*test_value[j]*area;
            integral_1[j] += val;

	      // second component   
            val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y + u3_values[0] * rhs2_z)
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * rhs2_x + rhs2 * rhs2_y + rhs3 * rhs2_z )
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * u2_values[1] + rhs2 * u2_values[2] + rhs3 * u2_values[3] )
		*test_value[j]*area*theta1*tau;
            val += (rhs2_x*test_value_x[j] + rhs2_y*test_value_y[j] + rhs2_z*test_value_z[j]) *area*eps*theta1*tau;
            val += rhs2*test_value[j]*area;
            integral_2[j] += val;

	      // third component   
            val = (u1_values[0] * rhs3_x + u2_values[0] * rhs3_y + u3_values[0] * rhs3_z)
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * rhs3_x + rhs2 * rhs3_y + rhs3 * rhs3_z )
              *test_value[j]*area*theta1*tau;
            val += (rhs1 * u3_values[1] + rhs2 * u3_values[2] + rhs3 * u3_values[3] )
		*test_value[j]*area*theta1*tau;

            val += (rhs3_x*test_value_x[j] + rhs3_y*test_value_y[j] + rhs3_z*test_value_z[j]) *area*eps*theta1*tau;
            val += rhs3*test_value[j]*area;
            integral_3[j] += val;
          }
        }                        // kc
      }                          // jc
    }                            // ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
      rhs[gdof+2*N_U] -= integral_3[j];
      value_global[0] += integral_1[j]*integral_1[j];
      value_global[1] += integral_2[j]*integral_2[j];
      value_global[2] += integral_3[j]*integral_3[j];
    }
    delete test_value;
    delete integral_1;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;

  OutPut("reduced value_global " << sqrt(value_global[0]) << " " << sqrt(value_global[1]) 
	 << " "  << sqrt(value_global[2]) << " " << endl);

}


void ApproximateTimeRFBSolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, dof[27], N_U, N_;
  int N_sub3, N_nodes, N_nodes3, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[27],y[27], z[27], a[729], b[729], *a_sub, *rhs_sub, test_bary;
  double xc[27], yc[27], zc[27], area, val, val1, val2, val3, xs, ys, zs, u1_values[4];
  double u2_values[4], u3_values[4], p_values[4];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value, *test_value_x, *test_value_y, *test_value_z;
  double rhs1, rhs2, rhs3, rhs1_x, rhs2_x, rhs3_x, rhs1_y, rhs2_y, rhs3_y, rhs1_z, rhs2_z, rhs3_z;
  double valU[3], gradU[9], delta, hK, mu_T;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;

  OutPut("rfb"<<endl);

  N_sub3 = N_sub * N_sub * N_sub;
  N_nodes =2* N_sub-1;
  N_nodes3 = N_nodes*N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes3*N_nodes3];
  memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);

  rhs_sub = new double[9 * N_nodes3];
  memset(rhs_sub, 0, 9*N_nodes3*SizeOfDouble);

  coeff = new double[4];

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }

    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    // compute delta for stabilization
    hK = cell->GetDiameter();
    delta = CharacteristicFilterWidth(hK);

    memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);
    memset(rhs_sub, 0, 9*N_nodes3*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell vertices
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);

          // coordinates on the edges
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[16] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[16] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[16] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          //for (jj=0;jj<8;jj++)
          //  OutPut(jj << " " << xc[jj]  << " " << yc[jj]  << " " << zc[jj] << endl);

          // initialize rhs
          memset(b,0,729*SizeOfDouble);

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            a[27*jj] = 1;
            a[27*jj+1] = xc[jj];
            a[27*jj+2] = yc[jj];
            a[27*jj+3] = zc[jj];
            a[27*jj+4] = xc[jj]*yc[jj];
            a[27*jj+5] = xc[jj]*zc[jj];
            a[27*jj+6] = yc[jj]*zc[jj];
            a[27*jj+7] = xc[jj]*yc[jj]*zc[jj];
            a[27*jj+8] = xc[jj]*xc[jj];
            a[27*jj+9] = yc[jj]*yc[jj];
            a[27*jj+10] = zc[jj]*zc[jj];
            a[27*jj+11] = a[27*jj+8]*yc[jj];
            a[27*jj+12] = a[27*jj+8]*zc[jj];
            a[27*jj+13] = xc[jj]*a[27*jj+9];
            a[27*jj+14] = a[27*jj+9]*zc[jj];
            a[27*jj+15] = xc[jj]*a[27*jj+10];
            a[27*jj+16] = yc[jj]*a[27*jj+10];
            a[27*jj+17] = a[27*jj+11]*zc[jj];
            a[27*jj+18] = a[27*jj+13]*zc[jj];
            a[27*jj+19] = xc[jj]*a[27*jj+16];
            a[27*jj+20] = a[27*jj+8]*a[27*jj+9];
            a[27*jj+21] = a[27*jj+8]*a[27*jj+10];
            a[27*jj+22] = a[27*jj+9]*a[27*jj+10];
            a[27*jj+23] = a[27*jj+20]*zc[jj];
            a[27*jj+24] = a[27*jj+21]*yc[jj];
            a[27*jj+25] = xc[jj]*a[27*jj+22];
            a[27*jj+26] = xc[jj]*a[27*jj+25];
            b[28*jj] = 1;
          }
          //for (j=0;j<16;j++)
          //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

          // solve system for the coefficients of the bilinear function

          SolveMultipleSystemsLapack(a,b,27,27,27,27);
          //SolveMultipleSystems(a,b,27,27,27,27);
          /*for (j=0;j<N_V;j++)
          OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
          b[4*j+2] << " " <<  b[4*j+3] << endl);
          */

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (2*ic-1)*N_nodes*N_nodes+ (2*jc-1)*N_nodes+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes*N_nodes;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes*N_nodes;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[2] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;
          }
          //OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << " ");
          //OutPut(dof[4] << " " << dof[5] << " " << dof[6] << " " << dof[7] << endl);

          // barycenter
          xs = xc[13];
          ys = yc[13];
          zs = zc[13];

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          p->FindGradientLocal(cell,i,xs,ys,zs,p_values);

          // compute value of rhs in (xs,ys,zs)
          Coeffs(1, &xs, &ys, &zs, NULL, &coeff);
          //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);

          valU[0]=u1_values[0];
          valU[1]=u2_values[0];
          valU[2]=u3_values[0];
          for(ii=0;ii<3;ii++)
          {
            gradU[3*ii]=u1_values[ii+1];
            gradU[3*ii+1]=u2_values[ii+1];
            gradU[3*ii+2]=u3_values[ii+1];
          }
          // turbulent viscosity
          mu_T = TurbulentViscosity3D(delta,gradU,valU,NULL,NULL,NULL,NULL,-4711);
          // update matrix entries
          // the matrix is stored column wise
          // ii -- test function
          for(ii=0;ii<27;ii++)
          {
            if (dof[ii] == -1)
              continue;

            // value of test function in bary center
            test_bary = b[8*ii] +  b[8*ii+1] * xs + b[8*ii+2] * ys + b[8*ii+3] * zs
              + b[8*ii+4] * xs * ys +  b[8*ii+5] * xs * zs +  b[8*ii+6] * ys * zs
              + b[8*ii+7] * xs * ys *zs + b[8*ii+8] * xs * xs + b[8*ii+9] * ys * ys + b[8*ii+10] * zs * zs
              + b[8*ii+11] * xs * xs * ys + b[8*ii+12] * xs * xs * zs + b[8*ii+13] * xs * ys * ys
              + b[8*ii+14] * ys * ys * zs + b[8*ii+15] * xs * zs * zs + b[8*ii+16] * ys * zs * zs
              + b[8*ii+17] * xs * xs * ys * zs + b[8*ii+18] * xs * ys * ys * zs
              + b[8*ii+19] * xs * ys * zs * zs + b[8*ii+20] * xs * xs * ys * ys + b[8*ii+21] * xs * xs * zs * zs
              + b[8*ii+22] * ys * ys * zs * zs + b[8*ii+23] * xs * xs * ys * ys * zs
              + b[8*ii+24] * xs * xs * ys * zs * zs + b[8*ii+25] * xs * ys * ys * zs * zs
              + b[8*ii+26] * xs * xs * ys * ys * zs * zs;
            //OutPut("test bary " << test_bary << endl);

            // compute first rhs
            val = u1_values[0]*(u1_values[1]+1/tau) +  u2_values[0]*u1_values[2] +  u3_values[0]*u1_values[3];
            val *= test_bary;
            rhs_sub[dof[ii]] -= area*val;
            // values of the derivatives of the test function in bary center
            val1 = b[8*ii+1]+ b[8*ii+4]*ys + b[8*ii+5]*zs + b[8*ii+7]*ys*zs + 2*b[8*ii+8]*xs
              + 2*b[8*ii+11]*xs*ys + 2*b[8*ii+12]*xs*zs + b[8*ii+13]*ys*ys + b[8*ii+15]*zs*zs
              + 2*b[8*ii+17]*xs*ys*zs + b[8*ii+18]*ys*ys*zs + b[8*ii+19]*ys*zs*zs
              + 2*b[8*ii+20]*xs*ys*ys + 2*b[8*ii+21]*xs*zs*zs + 2*b[8*ii+23]*xs*ys*ys*zs
              + 2*b[8*ii+24]*xs*ys*zs*zs+ b[8*ii+25]*ys*ys*zs*zs + 2*b[8*ii+26]*xs*ys*ys*zs*zs;
            val2 = b[8*ii+2]+ b[8*ii+4]*xs + b[8*ii+6]*zs + b[8*ii+7]*xs*zs + 2*b[8*ii+9]*ys
              + b[8*ii+11]*xs*xs + 2*b[8*ii+13]*xs*ys + 2*b[8*ii+14]*ys*zs + b[8*ii+16]*zs*zs
              + b[8*ii+17]*xs*xs*zs + 2*b[8*ii+18]*xs*ys*zs + b[8*ii+19]*xs*zs*zs
              + 2*b[8*ii+20]*xs*xs*ys + 2*b[8*ii+22]*ys*zs*zs + 2*b[8*ii+23]*xs*xs*ys*zs
              + b[8*ii+24]*xs*xs*zs*zs+ 2*b[8*ii+25]*xs*ys*zs*zs + 2*b[8*ii+26]*xs*xs*ys*zs*zs;
            val3 = b[8*ii+3]+ b[8*ii+5]*xs + b[8*ii+6]*ys + b[8*ii+7]*xs*ys + 2*b[8*ii+10]*zs
              + b[8*ii+12]*xs*xs + b[8*ii+14]*ys*ys + 2*b[8*ii+15]*xs*zs + 2*b[8*ii+16]*ys*zs
              + b[8*ii+17]*xs*xs*ys + b[8*ii+18]*xs*ys*ys + 2*b[8*ii+19]*xs*ys*zs
              + 2*b[8*ii+21]*xs*xs*zs + 2*b[8*ii+22]*ys*ys*zs + b[8*ii+23]*xs*xs*ys*ys
              + 2*b[8*ii+24]*xs*xs*ys*zs+ 2*b[8*ii+25]*xs*ys*ys*zs + 2*b[8*ii+26]*xs*xs*ys*ys*zs;
            rhs_sub[dof[ii]] -= (eps+mu_T)*area*(u1_values[1]*val1+u1_values[2]*val2+u1_values[3]*val3);

            // compute second rhs
            val = u1_values[0]*u2_values[1] +  u2_values[0]*(u2_values[2]+1/tau)  + u3_values[0]*u2_values[3];
            val *= test_bary;
            rhs_sub[dof[ii] + N_nodes3] -= area*val;
            rhs_sub[dof[ii] + N_nodes3] -= (eps+mu_T)*area*(u2_values[1]*val1+u2_values[2]*val2+u2_values[3]*val3);

            // compute third rhs
            val = u1_values[0]*u3_values[1] +  u2_values[0]*u3_values[2]  + u3_values[0]*(u3_values[3]+1/tau);
            val *= test_bary;
            rhs_sub[dof[ii]+ 2*N_nodes3] -= area*val;
            rhs_sub[dof[ii]+ 2*N_nodes3] -= (eps+mu_T)*area*(u3_values[1]*val1+u3_values[2]*val2+u3_values[3]*val3);

            // compute fourth, fifth, and sixth rhs
            rhs_sub[dof[ii]+ 3*N_nodes3] += area*(coeff[1]+u1_values[0]/tau)*test_bary;
            rhs_sub[dof[ii]+ 4*N_nodes3] += area*(coeff[2]+u2_values[0]/tau)*test_bary;
            rhs_sub[dof[ii]+ 5*N_nodes3] += area*(coeff[3]+u3_values[0]/tau)*test_bary;

            // compute rhs seven, eight, and nine
            rhs_sub[dof[ii]+ 6*N_nodes3] += area*p_values[0]*val1;
            rhs_sub[dof[ii]+ 7*N_nodes3] += area*p_values[0]*val2;
            rhs_sub[dof[ii]+ 8*N_nodes3] += area*p_values[0]*val3;

            // jj -- ansatz function
            for (jj=0;jj<27;jj++)
            {
              if (dof[jj] == -1)
                continue;
              // entry (jj,ii) exists
              // compute index in array containing the entry

              index = dof[ii]*N_nodes3+dof[jj];
              //OutPut(" ii " << ii << " jj " << jj << " ind " << index )
              // add Laplacian
              // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
              // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
              // approximate integral by midpoint rule
              // x part of the gradient
              val = (b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs + 2*b[8*jj+8]*xs
                + 2*b[8*jj+11]*xs*ys + 2*b[8*jj+12]*xs*zs + b[8*jj+13]*ys*ys + b[8*jj+15]*zs*zs
                + 2*b[8*jj+17]*xs*ys*zs + b[8*jj+18]*ys*ys*zs + b[8*jj+19]*ys*zs*zs
                + 2*b[8*jj+20]*xs*ys*ys + 2*b[8*jj+21]*xs*zs*zs + 2*b[8*jj+23]*xs*ys*ys*zs
                + 2*b[8*jj+24]*xs*ys*zs*zs+ b[8*jj+25]*ys*ys*zs*zs + 2*b[8*jj+26]*xs*ys*ys*zs*zs)*
                (b[8*ii+1]+ b[8*ii+4]*ys + b[8*ii+5]*zs + b[8*ii+7]*ys*zs + 2*b[8*ii+8]*xs
                + 2*b[8*ii+11]*xs*ys + 2*b[8*ii+12]*xs*zs + b[8*ii+13]*ys*ys + b[8*ii+15]*zs*zs
                + 2*b[8*ii+17]*xs*ys*zs + b[8*ii+18]*ys*ys*zs + b[8*ii+19]*ys*zs*zs
                + 2*b[8*ii+20]*xs*ys*ys + 2*b[8*ii+21]*xs*zs*zs + 2*b[8*ii+23]*xs*ys*ys*zs
                + 2*b[8*ii+24]*xs*ys*zs*zs+ b[8*ii+25]*ys*ys*zs*zs + 2*b[8*ii+26]*xs*ys*ys*zs*zs);
              // y part of the gradient
              val +=  (b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs + 2*b[8*jj+9]*ys
                + b[8*jj+11]*xs*xs + 2*b[8*jj+13]*xs*ys + 2*b[8*jj+14]*ys*zs + b[8*jj+16]*zs*zs
                + b[8*jj+17]*xs*xs*zs + 2*b[8*jj+18]*xs*ys*zs + b[8*jj+19]*xs*zs*zs
                + 2*b[8*jj+20]*xs*xs*ys + 2*b[8*jj+22]*ys*zs*zs + 2*b[8*jj+23]*xs*xs*ys*zs
                + b[8*jj+24]*xs*xs*zs*zs+ 2*b[8*jj+25]*xs*ys*zs*zs + 2*b[8*jj+26]*xs*xs*ys*zs*zs)*
                (b[8*ii+2]+ b[8*ii+4]*xs + b[8*ii+6]*zs + b[8*ii+7]*xs*zs + 2*b[8*ii+9]*ys
                + b[8*ii+11]*xs*xs + 2*b[8*ii+13]*xs*ys + 2*b[8*ii+14]*ys*zs + b[8*ii+16]*zs*zs
                + b[8*ii+17]*xs*xs*zs + 2*b[8*ii+18]*xs*ys*zs + b[8*ii+19]*xs*zs*zs
                + 2*b[8*ii+20]*xs*xs*ys + 2*b[8*ii+22]*ys*zs*zs + 2*b[8*ii+23]*xs*xs*ys*zs
                + b[8*ii+24]*xs*xs*zs*zs+ 2*b[8*ii+25]*xs*ys*zs*zs + 2*b[8*ii+26]*xs*xs*ys*zs*zs);
              // z part of the gradient
              val +=  (b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys + 2*b[8*jj+10]*zs
                + b[8*jj+12]*xs*xs + b[8*jj+14]*ys*ys + 2*b[8*jj+15]*xs*zs + 2*b[8*jj+16]*ys*zs
                + b[8*jj+17]*xs*xs*ys + b[8*jj+18]*xs*ys*ys + 2*b[8*jj+19]*xs*ys*zs
                + 2*b[8*jj+21]*xs*xs*zs + 2*b[8*jj+22]*ys*ys*zs + b[8*jj+23]*xs*xs*ys*ys
                + 2*b[8*jj+24]*xs*xs*ys*zs+ 2*b[8*jj+25]*xs*ys*ys*zs + 2*b[8*jj+26]*xs*xs*ys*ys*zs)*
                (b[8*ii+3]+ b[8*ii+5]*xs + b[8*ii+6]*ys + b[8*ii+7]*xs*ys + 2*b[8*ii+10]*zs
                + b[8*ii+12]*xs*xs + b[8*ii+14]*ys*ys + 2*b[8*ii+15]*xs*zs + 2*b[8*ii+16]*ys*zs
                + b[8*ii+17]*xs*xs*ys + b[8*ii+18]*xs*ys*ys + 2*b[8*ii+19]*xs*ys*zs
                + 2*b[8*ii+21]*xs*xs*zs + 2*b[8*ii+22]*ys*ys*zs + b[8*ii+23]*xs*xs*ys*ys
                + 2*b[8*ii+24]*xs*xs*ys*zs+ 2*b[8*ii+25]*xs*ys*ys*zs + 2*b[8*ii+26]*xs*xs*ys*ys*zs);
              val *= (eps+mu_T) * area;
              a_sub[index] += val;

              // add convective term
              // convection times gradient of ansatz fct.
              val = u1_values[0] * (b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs + 2*b[8*jj+8]*xs
                + 2*b[8*jj+11]*xs*ys + 2*b[8*jj+12]*xs*zs + b[8*jj+13]*ys*ys + b[8*jj+15]*zs*zs
                + 2*b[8*jj+17]*xs*ys*zs + b[8*jj+18]*ys*ys*zs + b[8*jj+19]*ys*zs*zs
                + 2*b[8*jj+20]*xs*ys*ys + 2*b[8*jj+21]*xs*zs*zs + 2*b[8*jj+23]*xs*ys*ys*zs
                + 2*b[8*jj+24]*xs*ys*zs*zs+ b[8*jj+25]*ys*ys*zs*zs + 2*b[8*jj+26]*xs*ys*ys*zs*zs);
              val += u2_values[0] * (b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs + 2*b[8*jj+9]*ys
                + b[8*jj+11]*xs*xs + 2*b[8*jj+13]*xs*ys + 2*b[8*jj+14]*ys*zs + b[8*jj+16]*zs*zs
                + b[8*jj+17]*xs*xs*zs + 2*b[8*jj+18]*xs*ys*zs + b[8*jj+19]*xs*zs*zs
                + 2*b[8*jj+20]*xs*xs*ys + 2*b[8*jj+22]*ys*zs*zs + 2*b[8*jj+23]*xs*xs*ys*zs
                + b[8*jj+24]*xs*xs*zs*zs+ 2*b[8*jj+25]*xs*ys*zs*zs + 2*b[8*jj+26]*xs*xs*ys*zs*zs);
              val += u3_values[0] * (b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys + 2*b[8*jj+10]*zs
                + b[8*jj+12]*xs*xs + b[8*jj+14]*ys*ys + 2*b[8*jj+15]*xs*zs + 2*b[8*jj+16]*ys*zs
                + b[8*jj+17]*xs*xs*ys + b[8*jj+18]*xs*ys*ys + 2*b[8*jj+19]*xs*ys*zs
                + 2*b[8*jj+21]*xs*xs*zs + 2*b[8*jj+22]*ys*ys*zs + b[8*jj+23]*xs*xs*ys*ys
                + 2*b[8*jj+24]*xs*xs*ys*zs+ 2*b[8*jj+25]*xs*ys*ys*zs + 2*b[8*jj+26]*xs*xs*ys*ys*zs);
              // multiplied with test fct.
              val *= test_bary;
              val *= area;
              a_sub[index] += val;

              // add time-derivative terms
              val = test_bary * ( b[8*jj] +  b[8*jj+1] * xs + b[8*jj+2] * ys + b[8*jj+3] * zs
                + b[8*jj+4] * xs * ys +  b[8*jj+5] * xs * zs +  b[8*jj+6] * ys * zs
                + b[8*jj+7] * xs * ys *zs + b[8*jj+8] * xs * xs + b[8*jj+9] * ys * ys + b[8*jj+10] * zs * zs
                + b[8*jj+11] * xs * xs * ys + b[8*jj+12] * xs * xs * zs + b[8*jj+13] * xs * ys * ys
                + b[8*jj+14] * ys * ys * zs + b[8*jj+15] * xs * zs * zs + b[8*jj+16] * ys * zs * zs
                + b[8*jj+17] * xs * xs * ys * zs + b[8*jj+18] * xs * ys * ys * zs
                + b[8*jj+19] * xs * ys * zs * zs + b[8*jj+20] * xs * xs * ys * ys + b[8*jj+21] * xs * xs * zs * zs
                + b[8*jj+22] * ys * ys * zs * zs + b[8*jj+23] * xs * xs * ys * ys * zs
                + b[8*jj+24] * xs * xs * ys * zs * zs + b[8*jj+25] * xs * ys * ys * zs * zs
                + b[8*jj+26] * xs * xs * ys * ys * zs * zs);
              val *= area/tau;
              a_sub[index] += val;

            }                    // end jj
          }                      // end ii
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*for (jj=0;jj<N_nodes3;jj++)
    {
    for (ii=0;ii<N_nodes3;ii++)
    {
      if (fabs(a_sub[jj+ii*N_nodes3])<1e-10)
    a_sub[jj+ii*N_nodes3] = 0;
      OutPut(a_sub[jj+ii*N_nodes3] << " ");
    }
    OutPut(endl);
    }*/

    // just for comparison with a known shape of solution
    /*for (j=0;j< N_nodes3 ; j++)
      OutPut("rhs " << rhs_sub[j+2*N_nodes3] << endl); */

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/

    //SolveMultipleSystems(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 9);
    SolveMultipleSystemsLapack(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 9);

    // just for comparison with a known shape of solution
    /*    for (j=0;j< N_nodes3 ; j++)
      OutPut(rhs_sub[j] << " " << rhs_sub[j+N_nodes3] << " " <<
      rhs_sub[j+2*N_nodes3] << " " << rhs_sub[j+3*N_nodes3]  << endl);*/
    /*	OutPut(" sol " << xlocal[j] << " " << ylocal[j] << " " <<
      rhs_sub[j+2*N_nodes3] << endl);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    // add contributions from the global function, the rhs, and the pressure
    for (j=0;j< 3 *N_nodes3 ; j++)
      rhs_sub[j] += rhs_sub[j + 3*N_nodes3]+rhs_sub[j + 6*N_nodes3];

    //OutPut(" rhs " << sqrt(Ddot(3*N_nodes3,rhs_sub,rhs_sub)));

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();
    test_value = new double[N_];
    test_value_x = new double[N_];
    test_value_y = new double[N_];
    test_value_z = new double[N_];
    integral_1 = new double[N_];
    integral_2 = new double[N_];
    integral_3 = new double[N_];
    memset(integral_1,0,N_*SizeOfDouble);
    memset(integral_2,0,N_*SizeOfDouble);
    memset(integral_3,0,N_*SizeOfDouble);

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          // compute coordinates of submesh cell vertices
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);

          // coordinates on the edges
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[15] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[15] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[15] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          dof[0] = (2*ic-1)*N_nodes*N_nodes+ (2*jc-1)*N_nodes+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes*N_nodes;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes*N_nodes;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[2] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;
          }
          // barycenter
          xs = xc[13];
          ys = yc[13];
          zs = zc[13];

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // compute the value of the test functions and its derivatives
          //in the bary center of the submesh cell
          TFEDatabase3D::GetRefFromOrig(RefTrans, xs, ys, zs, xi, eta, zeta);
          bf->GetDerivatives(D000, xi, eta, zeta, test_value);
          bf->GetDerivatives(D100, xi, eta, zeta, test_value_x);
          bf->GetDerivatives(D010, xi, eta, zeta, test_value_y);
          bf->GetDerivatives(D001, xi, eta, zeta, test_value_z);

          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the barycenter of the submesh cell
          // compute for this purpose the bilinear function on the mesh cell

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            a[27*jj] = 1;
            a[27*jj+1] = xc[jj];
            a[27*jj+2] = yc[jj];
            a[27*jj+3] = zc[jj];
            a[27*jj+4] = xc[jj]*yc[jj];
            a[27*jj+5] = xc[jj]*zc[jj];
            a[27*jj+6] = yc[jj]*zc[jj];
            a[27*jj+7] = xc[jj]*yc[jj]*zc[jj];
            a[27*jj+8] = xc[jj]*xc[jj];
            a[27*jj+9] = yc[jj]*yc[jj];
            a[27*jj+10] = zc[jj]*zc[jj];
            a[27*jj+11] = a[27*jj+8]*yc[jj];
            a[27*jj+12] = a[27*jj+8]*zc[jj];
            a[27*jj+13] = xc[jj]*a[27*jj+9];
            a[27*jj+14] = a[27*jj+9]*zc[jj];
            a[27*jj+15] = xc[jj]*a[27*jj+10];
            a[27*jj+16] = yc[jj]*a[27*jj+10];
            a[27*jj+17] = a[27*jj+11]*zc[jj];
            a[27*jj+18] = a[27*jj+13]*zc[jj];
            a[27*jj+19] = xc[jj]*a[27*jj+16];
            a[27*jj+20] = a[27*jj+8]*a[27*jj+9];
            a[27*jj+21] = a[27*jj+8]*a[27*jj+10];
            a[27*jj+22] = a[27*jj+9]*a[27*jj+10];
            a[27*jj+23] = a[27*jj+20]*zc[jj];
            a[27*jj+24] = a[27*jj+21]*yc[jj];
            a[27*jj+25] = xc[jj]*a[27*jj+22];
            a[27*jj+26] = xc[jj]*a[27*jj+25];

            if (dof[jj]!=-1)
            {
              b[jj] =  rhs_sub[dof[jj]];
              b[jj+27] = rhs_sub[dof[jj]+N_nodes3];
              b[jj+54] = rhs_sub[dof[jj]+2*N_nodes3];
            }
            else
              b[jj] = b[jj+27] = b[jj+54] = 0;
          }

          // solve system for the coefficients of the bilinear function

          SolveMultipleSystemsLapack(a,b,27,27,27,3);
          //SolveMultipleSystems(a,b,27,27,27,3);
          /*for (jj=0;jj<2;jj++)
            OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
            b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/

          // compute values of tilde_u in the barycenter
          rhs1 = b[0] +  b[1] * xs + b[2] * ys + b[3] * zs
            + b[4] * xs * ys +  b[5] * xs * zs +  b[6] * ys * zs
            + b[7] * xs * ys *zs + b[8] * xs * xs + b[9] * ys * ys + b[10] * zs * zs
            + b[11] * xs * xs * ys + b[12] * xs * xs * zs + b[13] * xs * ys * ys
            + b[14] * ys * ys * zs + b[15] * xs * zs * zs + b[16] * ys * zs * zs
            + b[17] * xs * xs * ys * zs + b[18] * xs * ys * ys * zs
            + b[19] * xs * ys * zs * zs + b[20] * xs * xs * ys * ys + b[21] * xs * xs * zs * zs
            + b[22] * ys * ys * zs * zs + b[23] * xs * xs * ys * ys * zs
            + b[24] * xs * xs * ys * zs * zs + b[25] * xs * ys * ys * zs * zs
            + b[26] * xs * xs * ys * ys * zs * zs;
          rhs2 = b[27] +  b[28] * xs + b[29] * ys + b[30] * zs
            + b[31] * xs * ys +  b[32] * xs * zs +  b[33] * ys * zs
            + b[34] * xs * ys *zs + b[35] * xs * xs + b[36] * ys * ys + b[37] * zs * zs
            + b[38] * xs * xs * ys + b[39] * xs * xs * zs + b[40] * xs * ys * ys
            + b[41] * ys * ys * zs + b[42] * xs * zs * zs + b[43] * ys * zs * zs
            + b[44] * xs * xs * ys * zs + b[45] * xs * ys * ys * zs
            + b[46] * xs * ys * zs * zs + b[47] * xs * xs * ys * ys + b[48] * xs * xs * zs * zs
            + b[49] * ys * ys * zs * zs + b[50] * xs * xs * ys * ys * zs
            + b[51] * xs * xs * ys * zs * zs + b[52] * xs * ys * ys * zs * zs
            + b[53] * xs * xs * ys * ys * zs * zs;
          rhs3 = b[54] +  b[55] * xs + b[56] * ys + b[57] * zs
            + b[58] * xs * ys +  b[59] * xs * zs +  b[60] * ys * zs
            + b[61] * xs * ys *zs + b[62] * xs * xs + b[63] * ys * ys + b[64] * zs * zs
            + b[65] * xs * xs * ys + b[66] * xs * xs * zs + b[67] * xs * ys * ys
            + b[68] * ys * ys * zs + b[69] * xs * zs * zs + b[70] * ys * zs * zs
            + b[71] * xs * xs * ys * zs + b[72] * xs * ys * ys * zs
            + b[73] * xs * ys * zs * zs + b[74] * xs * xs * ys * ys + b[75] * xs * xs * zs * zs
            + b[76] * ys * ys * zs * zs + b[77] * xs * xs * ys * ys * zs
            + b[78] * xs * xs * ys * zs * zs + b[79] * xs * ys * ys * zs * zs
            + b[80] * xs * xs * ys * ys * zs * zs;

          // compute gradients of tilde_u in the barycenter
          rhs1_x = b[1]+ b[4]*ys + b[5]*zs + b[7]*ys*zs + 2*b[8]*xs
            + 2*b[11]*xs*ys + 2*b[12]*xs*zs + b[13]*ys*ys + b[15]*zs*zs
            + 2*b[17]*xs*ys*zs + b[18]*ys*ys*zs + b[19]*ys*zs*zs
            + 2*b[20]*xs*ys*ys + 2*b[21]*xs*zs*zs + 2*b[23]*xs*ys*ys*zs
            + 2*b[24]*xs*ys*zs*zs+ b[25]*ys*ys*zs*zs + 2*b[26]*xs*ys*ys*zs*zs;
          rhs1_y = b[2]+ b[4]*xs + b[6]*zs + b[7]*xs*zs + 2*b[9]*ys
            + b[11]*xs*xs + 2*b[13]*xs*ys + 2*b[14]*ys*zs + b[16]*zs*zs
            + b[17]*xs*xs*zs + 2*b[18]*xs*ys*zs + b[19]*xs*zs*zs
            + 2*b[20]*xs*xs*ys + 2*b[22]*ys*zs*zs + 2*b[23]*xs*xs*ys*zs
            + b[24]*xs*xs*zs*zs+ 2*b[25]*xs*ys*zs*zs + 2*b[26]*xs*xs*ys*zs*zs;
          rhs1_z = b[3]+ b[5]*xs + b[6]*ys + b[7]*xs*ys + 2*b[10]*zs
            + b[12]*xs*xs + b[14]*ys*ys + 2*b[15]*xs*zs + 2*b[16]*ys*zs
            + b[17]*xs*xs*ys + b[18]*xs*ys*ys + 2*b[19]*xs*ys*zs
            + 2*b[21]*xs*xs*zs + 2*b[22]*ys*ys*zs + b[23]*xs*xs*ys*ys
            + 2*b[24]*xs*xs*ys*zs+ 2*b[25]*xs*ys*ys*zs + 2*b[26]*xs*xs*ys*ys*zs;
          rhs2_x = b[28]+ b[31]*ys + b[32]*zs + b[34]*ys*zs + 2*b[35]*xs
            + 2*b[37]*xs*ys + 2*b[39]*xs*zs + b[40]*ys*ys + b[42]*zs*zs
            + 2*b[44]*xs*ys*zs + b[45]*ys*ys*zs + b[46]*ys*zs*zs
            + 2*b[47]*xs*ys*ys + 2*b[48]*xs*zs*zs + 2*b[50]*xs*ys*ys*zs
            + 2*b[51]*xs*ys*zs*zs+ b[52]*ys*ys*zs*zs + 2*b[53]*xs*ys*ys*zs*zs;
          rhs2_y = b[29]+ b[31]*xs + b[33]*zs + b[34]*xs*zs + 2*b[36]*ys
            + b[38]*xs*xs + 2*b[40]*xs*ys + 2*b[41]*ys*zs + b[43]*zs*zs
            + b[44]*xs*xs*zs + 2*b[45]*xs*ys*zs + b[46]*xs*zs*zs
            + 2*b[47]*xs*xs*ys + 2*b[49]*ys*zs*zs + 2*b[50]*xs*xs*ys*zs
            + b[51]*xs*xs*zs*zs+ 2*b[52]*xs*ys*zs*zs + 2*b[53]*xs*xs*ys*zs*zs;
          rhs2_z = b[30]+ b[32]*xs + b[33]*ys + b[34]*xs*ys + 2*b[37]*zs
            + b[39]*xs*xs + b[41]*ys*ys + 2*b[42]*xs*zs + 2*b[43]*ys*zs
            + b[44]*xs*xs*ys + b[45]*xs*ys*ys + 2*b[46]*xs*ys*zs
            + 2*b[48]*xs*xs*zs + 2*b[49]*ys*ys*zs + b[50]*xs*xs*ys*ys
            + 2*b[51]*xs*xs*ys*zs+ 2*b[52]*xs*ys*ys*zs + 2*b[53]*xs*xs*ys*ys*zs;
          rhs3_x = b[55]+ b[58]*ys + b[59]*zs + b[61]*ys*zs + 2*b[62]*xs
            + 2*b[65]*xs*ys + 2*b[66]*xs*zs + b[67]*ys*ys + b[69]*zs*zs
            + 2*b[71]*xs*ys*zs + b[72]*ys*ys*zs + b[73]*ys*zs*zs
            + 2*b[74]*xs*ys*ys + 2*b[75]*xs*zs*zs + 2*b[77]*xs*ys*ys*zs
            + 2*b[78]*xs*ys*zs*zs+ b[79]*ys*ys*zs*zs + 2*b[80]*xs*ys*ys*zs*zs;
          rhs3_y = b[56]+ b[58]*xs + b[60]*zs + b[61]*xs*zs + 2*b[63]*ys
            + b[65]*xs*xs + 2*b[67]*xs*ys + 2*b[68]*ys*zs + b[70]*zs*zs
            + b[71]*xs*xs*zs + 2*b[72]*xs*ys*zs + b[73]*xs*zs*zs
            + 2*b[74]*xs*xs*ys + 2*b[76]*ys*zs*zs + 2*b[77]*xs*xs*ys*zs
            + b[78]*xs*xs*zs*zs+ 2*b[79]*xs*ys*zs*zs + 2*b[80]*xs*xs*ys*zs*zs;
          rhs3_z = b[57]+ b[59]*xs + b[60]*ys + b[61]*xs*ys + 2*b[64]*zs
            + b[66]*xs*xs + b[68]*ys*ys + 2*b[69]*xs*zs + 2*b[70]*ys*zs
            + b[71]*xs*xs*ys + b[72]*xs*ys*ys + 2*b[73]*xs*ys*zs
            + 2*b[75]*xs*xs*zs + 2*b[76]*ys*ys*zs + b[77]*xs*xs*ys*ys
            + 2*b[78]*xs*xs*ys*zs+ 2*b[79]*xs*ys*ys*zs + 2*b[80]*xs*xs*ys*ys*zs;

          // update the integrals
          for (j=0;j<N_;j++)
          {
            val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y + u3_values[0] * rhs1_z )*test_value[j]*area;
            integral_1[j] += val;

            val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y + u3_values[0] * rhs2_z)*test_value[j]*area;
            integral_2[j] += val;

            val = (u1_values[0] * rhs3_x + u2_values[0] * rhs3_y + u3_values[0] * rhs3_z)*test_value[j]*area;
            integral_3[j] += val;

            //add laplacian
            val = (rhs1_x*test_value_x[j] + rhs1_y*test_value_y[j] + rhs1_z*test_value_z[j]) *area*eps;
            integral_1[j] += val;

            val = (rhs2_x*test_value_x[j] + rhs2_y*test_value_y[j] + rhs2_z*test_value_z[j]) *area*eps;
            integral_2[j] += val;

            val = (rhs3_x*test_value_x[j] + rhs3_y*test_value_y[j] + rhs3_z*test_value_z[j]) *area*eps;
            integral_3[j] += val;

          }
        }                        // kc
      }                          // jc
    }                            // ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= tau*integral_1[j];
      rhs[gdof+N_U] -= tau*integral_2[j];
      rhs[gdof+2*N_U] -= tau*integral_3[j];
    }
    delete test_value;
    delete test_value_x;
    delete test_value_y;
    delete test_value_z;
    delete integral_1;
    delete integral_2;
  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;
}

/***************************************************************************************/
//
// bubble VMS method with solution of coupled system on the mesh cells
//
// velocity is approximated by Q_2 bubbles
// pressure is approximated by Q_1
// on dof of the pressure is set to be zero
//
/***************************************************************************************/

void ApproximateTimeRFB_coupled_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, kk, N_U, N_, iq;
  int dof[27], dofp[8];
  int N_sub3, N_nodes_u, N_nodes_u3, N_vec_u, N_nodes_p, N_nodes_p3, N_all, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof, index1;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[27],y[27], z[27], a[729],ap[64], b[729], bp[64],*a_sub, *rhs_sub;
  double fe_fct_p;
  double xc[27], yc[27], zc[27], xp[8], yp[8], zp[8];
  double area, val, val0;
  double u1_values[4], u2_values[4], u3_values[4], p_values[4];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value, *test_value_x, *test_value_y, *test_value_z;
  double delta, hK, mu_T, temp, temp1;
  double val_test[4], val_ansatz[4], xq, yq, zq, detJK;
  double rhs1[4], rhs2[4], rhs3[4], rhsp, valU[3], gradU[9];
  //double *a0, *a1, *a2, *b0, *b1, *b2;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  //double tau = 1;
  double theta1 = TDatabase::TimeDB->THETA1;

  /*************************************************************************/
  /* data for quadrature rule gauss2                                       */
  /*************************************************************************/
  /*int quad_points = 8;
  double weight[8]={ 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0 };
  double qx[8]={-0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489};
  double qy[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};
  double qz[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};*/
  /*************************************************************************/
  /* data for quadrature rule gauss3                                       */
  /*************************************************************************/
  int quad_points = 27;
  double oneweight[3]={ 0.5555555555555555555555558,
                        0.8888888888888888888888888,
                        0.5555555555555555555555558 };
  double onexi[3]={ -0.7745966692414833770358530,
                     0,
                     0.7745966692414833770358530 };

  double weight[27], qx[27], qy[27],qz[27];
  for(ii=0;ii<3;ii++)
    for(jj=0;jj<3;jj++)
      for(kk=0;kk<3;kk++)
      {
        i=9*ii+3*jj+kk;
        weight[i]=oneweight[ii]*oneweight[jj]*oneweight[kk];
        qx[i]=onexi[kk];
        qy[i]=onexi[jj];
        qz[i]=onexi[ii];
      }
  
  /*************************************************************************/
  /* data for quadrature rule gauss4                                       */
  /*************************************************************************/
  /*int quad_points = 64;
  double oneweight[4]=
  {
    0.652145154862546142626936051,
    0.347854845137453857373063949,
    0.347854845137453857373063949,
    0.652145154862546142626936051
  };

  double onexi[4]=
  {
    0.33998104358485626480266576,
    0.86113631159405257522394649,
    -0.86113631159405257522394649,
    -0.33998104358485626480266576
  };

  double weight[64], qx[64], qy[64],qz[64];
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      for(kk=0;kk<4;kk++)
      {
        i=16*ii+4*jj+kk;
        weight[i]=oneweight[ii]*oneweight[jj]*oneweight[kk];
        qx[i]=onexi[kk];
        qy[i]=onexi[jj];
        qz[i]=onexi[ii];
      }
  */
  OutPut("rfb_coupled"<<endl);

  // compute number of local mesh cells
  N_sub3 = N_sub * N_sub * N_sub;
  // compute number of dof for the local problems
  // do not consider velocity nodes on the boundary of the mesh cell
  // for Q2
  N_nodes_u = 2*N_sub-1;
  N_nodes_u3 = N_nodes_u*N_nodes_u*N_nodes_u;

  // there is no boundary condition for p, for Q1
  N_nodes_p = N_sub+1;
  N_nodes_p3 = N_nodes_p*N_nodes_p*N_nodes_p;
  N_vec_u = 3*N_nodes_u3;
  N_all = N_vec_u+N_nodes_p3;

  //OutPut("u " << N_nodes_u << " u3 " << N_nodes_u3  << " p " << N_nodes_p3 << " all " << N_all << endl);

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is 2*N_sub-1 for the velocity (Q2)
  // and N_sub-1 for the pressure (Q1)

  a_sub = new double[N_all*N_all];
  memset(a_sub, 0, N_all*N_all*SizeOfDouble);

  rhs_sub = new double[N_all];

  coeff = new double[4];

  // write the coordinates of the nodes into the vectors
  // coord: first N_nodes_u3 components correspond to the
  // velocity nodes, and the last N_nodes_p3, to the pressure nodes
  /*double coord_x[N_nodes_u3+N_nodes_p3];
  double coord_y[N_nodes_u3+N_nodes_p3];
  double coord_z[N_nodes_u3+N_nodes_p3];*/

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    //OutPut("cell " << i << endl);
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }

    // The unit cube:
    /* x[0] = 0; y[0] = 0; z[0] = 0;
    x[1] = 1; y[1] = 0; z[1] = 0;
    x[2] = 1; y[2] = 1; z[2] = 0;
    x[3] = 0; y[3] = 1; z[3] = 0;
    x[4] = 0; y[4] = 0; z[4] = 1;
    x[5] = 1; y[5] = 0; z[5] = 1;
    x[6] = 1; y[6] = 1; z[6] = 1;
    x[7] = 0; y[7] = 1; z[7] = 1;*/

    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    detJK = area/8.0;

    // compute delta for stabilization
    /* if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)
      hK = cell->GetShortestEdge();
      else*/
      hK = cell->GetDiameter();
    // measure of the submesh cells
    hK /= N_sub;
    delta = CharacteristicFilterWidth(hK);

    memset(a_sub, 0, N_all*N_all*SizeOfDouble);
    memset(rhs_sub, 0, N_all*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/
    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell vertices
          // for the velocity
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);
          // positions for the pressure dof (Q1)
          xp[0] = xc[0];
          yp[0] = yc[0];
          zp[0] = zc[0];
          //OutPut("press " << xp[0] << " " << yp[0] << " " << zp[0] << endl);
          xp[1] = xc[2];
          yp[1] = yc[2];
          zp[1] = zc[2];

          xp[2] = xc[8];
          yp[2] = yc[8];
          zp[2] = zc[8];

          xp[3] = xc[6];
          yp[3] = yc[6];
          zp[3] = zc[6];

          xp[4] = xc[18];
          yp[4] = yc[18];
          zp[4] = zc[18];

          xp[5] = xc[20];
          yp[5] = yc[20];
          zp[5] = zc[20];

          xp[6] = xc[26];
          yp[6] = yc[26];
          zp[6] = zc[26];

          xp[7] = xc[24];
          yp[7] = yc[24];
          zp[7] = zc[24];

          // coordinates on the edges - for the velocity
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[16] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[16] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[16] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          // initialize rhs
          memset(b,0,729*SizeOfDouble);
          memset(bp,0,64*SizeOfDouble);
          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            index = 27*jj;
            a[index] = 1;
            a[index+1] = xc[jj];
            a[index+2] = yc[jj];
            a[index+3] = zc[jj];
            a[index+4] = xc[jj]*yc[jj];
            a[index+5] = xc[jj]*zc[jj];
            a[index+6] = yc[jj]*zc[jj];
            a[index+7] = xc[jj]*yc[jj]*zc[jj];
            a[index+8] = xc[jj]*xc[jj];
            a[index+9] = yc[jj]*yc[jj];
            a[index+10] = zc[jj]*zc[jj];
            a[index+11] = a[index+8]*yc[jj];
            a[index+12] = a[index+8]*zc[jj];
            a[index+13] = xc[jj]*a[index+9];
            a[index+14] = a[index+9]*zc[jj];
            a[index+15] = xc[jj]*a[index+10];
            a[index+16] = yc[jj]*a[index+10];
            a[index+17] = a[index+11]*zc[jj];
            a[index+18] = a[index+13]*zc[jj];
            a[index+19] = xc[jj]*a[index+16];
            a[index+20] = a[index+8]*a[index+9];
            a[index+21] = a[index+8]*a[index+10];
            a[index+22] = a[index+9]*a[index+10];
            a[index+23] = a[index+20]*zc[jj];
            a[index+24] = a[index+21]*yc[jj];
            a[index+25] = xc[jj]*a[index+22];
            a[index+26] = xc[jj]*a[index+25];
            b[28*jj] = 1;
          }
          //for (jj=0;jj<27*27;jj++) OutPut("i="<<i<<"  ic="<<ic<<"  jc="<<jc<<"  kc="<<kc<<"  a1("<<jj<<")= "<<a[jj]<<endl);
          // for (ii=0;ii<27;ii++)
          // for (jj=0;jj<27;jj++)
          // OutPut("a1("<<ii+1<<","<< jj+1<<")= "<<a[ii*27+jj]<<endl);
          // set matrices for computation of the coefficients
          // of the trilinear function
          for (jj=0;jj<8;jj++)
          {
            index=8*jj;
            ap[index] = 1;
            ap[index+1] = xp[jj];
            ap[index+2] = yp[jj];
            ap[index+3] = zp[jj];
            ap[index+4] = xp[jj]*yp[jj];
            ap[index+5] = xp[jj]*zp[jj];
            ap[index+6] = yp[jj]*zp[jj];
            ap[index+7] = xp[jj]*yp[jj]*zp[jj];
            bp[9*jj] = 1;
          }

          // solve systems for the coefficients

          //SolveMultipleSystems(ap,bp,8,8,8,8);
          //SolveMultipleSystems(a,b,27,27,27,27);

          SolveMultipleSystemsLapack(ap,bp,8,8,8,8);
          SolveMultipleSystemsLapack(a,b,27,27,27,27);
          //         for (jj=0;jj<27;jj++)
          //           OutPut("n0 " <<jj << " " << b[jj] << endl);

          /*for (ii=0;ii<27;ii++)
              for (jj=0;jj<27;jj++)
               {Compute_Q2_Value_Gradient(b+(27*ii),xc[jj],yc[jj],zc[jj],val_test);
                  OutPut(ii << " " << jj << " "<< val_test[0] << endl);
                }
          exit(1);*/

          // compute local velocity dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (2*ic-1)*N_nodes_u*N_nodes_u+ (2*jc-1)*N_nodes_u+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes_u;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes_u;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes_u*N_nodes_u;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes_u;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes_u;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes_u*N_nodes_u;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes_u;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes_u;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          /*OutPut("vor ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          // different than velocity since there are nodes on
          // the boundary of the mesh cell i
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;

          //OutPut("pdof " << dofp[0] << " "<< dofp[1] << " "<< dofp[2] << " "<< dofp[3] << " "
          //	<< dofp[4] << " "<< dofp[5] << " "<< dofp[6] << " "<< dofp[7] << endl);

          // correct (for the velocity) for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
          }
          /*OutPut("nacha ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }

          /*OutPut("nachb ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
          }

          /*OutPut("nachc ");
           for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }

          /*OutPut("nachd ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;
          }

          /*OutPut("nache ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;
          }

          /*OutPut("nachf ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // write the coordinates of the nodes into the vectors
          // coord: first N_nodes_u3 components correspond to the
          // velocity nodes, and the last N_nodes_p3, to the pressure nodes
          /*for (ii=0;ii<27;ii++)
          {
             if (dof[ii] != -1)
             {
               index = dof[ii];
               coord_x[index] = xc[ii];
               coord_y[index] = yc[ii];
               coord_z[index] = zc[ii];
             }
          }
          for (ii=0;ii<8;ii++)
          {
            index = dofp[ii];
            coord_x[N_nodes_u3+index] = xp[ii];
            coord_y[N_nodes_u3+index] = yp[ii];
            coord_z[N_nodes_u3+index] = zp[ii];
          }*/

          // loop over the quadrature points
          for (iq = 0;iq < quad_points; iq++)
          {
            temp1=detJK*weight[iq];
            temp=temp1*tau;

            // quadrature points
            // ONLY FOR PARALLELEPIPED !!!
            xq = ((xc[2]-xc[0])*qx[iq]+(xc[8]-xc[2])*qy[iq]+(xc[26]-xc[8])*qz[iq])/2+xc[13];
            yq = ((yc[2]-yc[0])*qx[iq]+(yc[8]-yc[2])*qy[iq]+(yc[26]-yc[8])*qz[iq])/2+yc[13];
            zq = ((zc[2]-zc[0])*qx[iq]+(zc[8]-zc[2])*qy[iq]+(zc[26]-zc[8])*qz[iq])/2+zc[13];

            // compute convection field in the quadrature point
            u1->FindGradientLocal(cell,i,xq,yq,zq,u1_values);
            u2->FindGradientLocal(cell,i,xq,yq,zq,u2_values);
            u3->FindGradientLocal(cell,i,xq,yq,zq,u3_values);
            p->FindGradientLocal(cell,i,xq,yq,zq,p_values);

            //u1_values[0] = u2_values[0] = u3_values[0]  = 0;
            //p_values[0] = 0;
            // compute value of rhs in quadrature point
            Coeffs(1, &xq, &yq, &zq, NULL, &coeff);
            //OutPut("coeff " << coeff[1] << " " << coeff[2] << endl);
            //coeff[1] = 1; coeff[2] = 0; coeff[3] = 0;
            // prepare computation of turbulent viscosity
            valU[0]=u1_values[0];
            valU[1]=u2_values[0];
            valU[2]=u3_values[0];
            for(ii=0;ii<3;ii++)
            {
              gradU[3*ii]=u1_values[ii+1];
              gradU[3*ii+1]=u2_values[ii+1];
              gradU[3*ii+2]=u3_values[ii+1];
            }
            // turbulent viscosity
            mu_T = TurbulentViscosity3D(delta,gradU,valU,NULL,NULL,NULL,&zq,-4711);
	   
            // update rhs and matrix entries
            // the matrix is stored column wise

            // ii -- test function
            for(ii=0;ii<27;ii++)
            {
              // Dirichlet value
              if (dof[ii] == -1)
                continue;
              index=dof[ii];
              // get values of test function
              Compute_Q2_Value_Gradient(b+(27*ii),xq,yq,zq,val_test);
              //OutPut(val_test[0] << " " << val_test[1] << " " << val_test[2] << " " << val_test[3] << " " << endl);

              // compute first rhs (u1)
              // convective term
              val = u1_values[0]*u1_values[1] +  u2_values[0]*u1_values[2] +  u3_values[0]*u1_values[3];
              val *= val_test[0] * temp;
              rhs_sub[index] -= val;
              // diffusive term
              rhs_sub[index] -= eps*temp*(u1_values[1]*val_test[1]
					  +u1_values[2]*val_test[2]+u1_values[3]*val_test[3]);
              // add (f1,phi_ii)
              //OutPut(detJK<< " " << weight[iq]<< " " <<  coeff[1]<< " " << val_test[0]<< " " << tau << endl);
              rhs_sub[index] += temp*coeff[1]*val_test[0];
              // pressure term
	      rhs_sub[index] += temp*p_values[0]*val_test[1];

              // compute second rhs (u2)
              // convective term
              val = u1_values[0]*u2_values[1] +  u2_values[0]*u2_values[2]  + u3_values[0]*u2_values[3];
              val *= val_test[0] * temp;
              rhs_sub[index + N_nodes_u3] -= val;
              // diffusive term
              rhs_sub[index + N_nodes_u3] -= eps*temp*(u2_values[1]*val_test[1]+
						       u2_values[2]*val_test[2]+u2_values[3]*val_test[3]);
              // add (f2,phi_ii)
              rhs_sub[index + N_nodes_u3] += temp*coeff[2]*val_test[0];
              // pressure term
              rhs_sub[index + N_nodes_u3] += temp*p_values[0]*val_test[2];

              // compute third rhs (u3)
              // convective term
              val = u1_values[0]*u3_values[1] +  u2_values[0]*u3_values[2]  + u3_values[0]*u3_values[3];
              val *= val_test[0]*temp;
              rhs_sub[index+ 2*N_nodes_u3] -= val;
              // diffusive term
              rhs_sub[index+ 2*N_nodes_u3] -= eps*temp*(u3_values[1]*val_test[1]+
							u3_values[2]*val_test[2]+u3_values[3]*val_test[3]);
              // add (f3,phi_ii)
              rhs_sub[index + 2*N_nodes_u3] += temp*coeff[3]*val_test[0];
              // pressure term
              rhs_sub[index + 2*N_nodes_u3] += temp*p_values[0]*val_test[3];

              // jj -- ansatz function
              for (jj=0;jj<27;jj++)
              {
                if (dof[jj] == -1)
                  continue;
                // entry (jj,ii) exists
                // compute index in array containing the entry

                // A11
                index = dof[ii]*N_all+dof[jj];
                // OutPut(" ii " << ii << " jj " << jj << " ind " << index << endl);
                // values of ansatz function in quadrature points
                Compute_Q2_Value_Gradient(b+(27*jj),xq,yq,zq,val_ansatz);
                // add Laplacian
                val = val_ansatz[1]*val_test[1] + val_ansatz[2]*val_test[2] +val_ansatz[3]*val_test[3];
                //val *= eps*temp;
                val *= (eps+mu_T) *temp;
                a_sub[index] += val;

                // add convective term
                // convection times gradient of ansatz fct.
                val = u1_values[0] * val_ansatz[1] + u2_values[0] * val_ansatz[2] + u3_values[0] * val_ansatz[3];
                // multiplied with test fct.
                val *= val_test[0]*temp;
		a_sub[index] += val;

                // add mass term
                // A11
                val = val_ansatz[0] * val_test[0];
                a_sub[index] += val * temp1;

                // copy to the other diagonal blocks
                // A22
                index1 = index + (N_all+1)*N_nodes_u3;
                a_sub[index1] = a_sub[index];

                // A33
                index1 = index + 2*(N_all+1)*N_nodes_u3;
                a_sub[index1] = a_sub[index];
              }                  // end jj

              // right upper pressure matrix
              for (kk=0;kk<8;kk++)
              {
                index = dof[ii]*N_all+N_vec_u+dofp[kk];

                /*if (dof[ii]==13)
                OutPut(dofp[kk] << " " << index << " " << 13*N_all+N_vec_u-1 << " " << 14*N_all << endl);*/

                // value of pressure test function in bary center
                Compute_Q1_Value_RFB(bp+(8*kk), xq, yq, zq, &fe_fct_p);

                // B1

                /*if ((index>=13*N_all+N_vec_u)&&(index<14*N_all)&&((dofp[kk]==12)||(dofp[kk]==12)))
                         OutPut("vor " << index << " a " << a_sub[index] << " u " << dof[ii] << " p "
                  << dofp[kk] << endl);*/

                val =  fe_fct_p * val_test[1];
                a_sub[index] -= val*temp;

                /*if ((index>=13*N_all+N_vec_u)&&(index<14*N_all)&&((dofp[kk]==12)||(dofp[kk]==12)))
                         OutPut("nach " << index << " a " << a_sub[index] << " u " << dof[ii] << " p "
                  << dofp[kk] << endl);*/

                // B2
                index += N_all*N_nodes_u3;
                val = fe_fct_p * val_test[2];
                a_sub[index] -= val*temp;

                // B3
                index += N_all*N_nodes_u3;
                val = fe_fct_p * val_test[3];
                a_sub[index] -= val*temp;
              }                  //end kk (ansatz pressure)
            }                    // end ii (test velocity)

            // left lower pressure block
            // test functions of pressure
            for (kk=0;kk<8;kk++) // test functions of pressure
            {
              // value of pressure test function in quadrature point
              Compute_Q1_Value_RFB(bp+(8*kk), xq, yq, zq, &fe_fct_p);
              // compute forth rhs (continuity eq.)
              val = (u1_values[1] + u2_values[2] + u3_values[3]) * fe_fct_p;
              index=dofp[kk];
              //rhs_sub[index + N_vec_u] -= val*temp;
              // ansatz functions of velocity
              for(jj=0;jj<27;jj++)
              {
                if (dof[jj] == -1)
                  continue;
                index = (N_vec_u+dofp[kk])*N_all+dof[jj];

                // value of pressure test function in quadrature points

                //OutPut(jj << " ind " << index << endl);

                // values of ansatz function in quadrature points
                Compute_Q2_Value_Gradient(b+(27*jj),xq,yq,zq,val_ansatz);
                // B1^T
                val = val_ansatz[1] * fe_fct_p;
                a_sub[index] += val*temp;
                // B2^T
                index += N_nodes_u3;
                val = val_ansatz[2] * fe_fct_p;
                a_sub[index] += val*temp;
                // B3^T
                index += N_nodes_u3;
                val = val_ansatz[3] * fe_fct_p;
                a_sub[index] += val*temp;
              }                  // end jj
            }                    //end kk
          }                      // end iq
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    /* for (ic=0;ic<N_all;ic++)
         for (jc=0;jc<N_all;jc++)
           OutPut("a("<<ic<<","<<jc<<") = " << a_sub[ic*N_all+jc] << endl);*/

    // replace the last line of the matrix on the left hand side
    // with (0, 0, ..., 0, 1), and the last element of the vector
    // on the right hand side with 0

    for(ii=(N_all-1)*N_all;ii<N_all*N_all-1;ii++)
      a_sub[ii]=0;
    a_sub[N_all*N_all-1]=1;
    rhs_sub[N_all-1]=0;

    //	   OutPut("mat: " << sqrt(Ddot(N_all*N_all,a_sub,a_sub)));
    /*   for(ii=0;ii<N_all;ii++)
         OutPut("rhs("<<ii<<") = "<<rhs_sub[ii]<<endl);
    */

    // write out the components of the matrices on the left hand side
    /*
     for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_u3;jj++)
          OutPut("a0("<<ii<<", "<<jj<<")=  "<<a_sub[ii*N_all+jj]<<endl);

    */
    /*for(ii=0;ii<N_all;ii++)
       for(jj=0;jj<N_all;jj++)
       {
    if (fabs(a_sub[ii*N_all+jj] - a_sub[jj*N_all+ii]) > 1e-12)
    {
          OutPut(setw(7) << ii << " " << jj << " " << setprecision(5)<<a_sub[ii*N_all+jj]<<"   ");
          OutPut(setw(7) << setprecision(5)<<a_sub[jj*N_all+ii]<<endl);
    }
       }*/

    /* for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii)*N_all+3*N_nodes_u3+jj])>1e-12)
          OutPut("b0("<<ii<<", "<<jj<<")=  "<<a_sub[(ii)*N_all+N_vec_u+jj]<< endl);

     for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii+N_nodes_u3)*N_all+N_vec_u+jj])>1e-12)
         OutPut("b1("<<ii<<", "<<jj<<")=  "<<a_sub[(ii+N_nodes_u3)*N_all+N_vec_u+jj]<<endl);

     for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii+2*N_nodes_u3)*N_all+N_vec_u+jj])>1e-12)
         OutPut("b2("<<ii<<", "<<jj<<")=  "<<a_sub[(ii+2*N_nodes_u3)*N_all+N_vec_u+jj]<<endl);

    //  OutPut(" mat " << sqrt(Ddot(N_all*N_all,a_sub,a_sub)));

        for(ii=0;ii<3*N_nodes_u3+N_nodes_p3;ii++)
        {
      val = 0;
        for(jj=0;jj<3*N_nodes_u3+N_nodes_p3;jj++)
      val += a_sub[ii*(3*N_nodes_u3+N_nodes_p3)+jj] *sol_sub[jj];
      rhs_sub[ii] = val - rhs_sub[ii];
      if (fabs(rhs_sub[ii])>1e-12)
      OutPut("rhs " << ii << " " << rhs_sub[ii] << endl);
    }
      exit(1);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    //OutPut(" rhs " << rhs_sub[N_all-1]<<"   "<< rhs_sub[N_all-2]<<"   "<<sqrt(Ddot(N_all,rhs_sub,rhs_sub)) << endl);
    //SolveLinearSystem(a_sub, rhs_sub, N_all, N_all);
    SolveLinearSystemTranspose(a_sub, rhs_sub, N_all, N_all);

    // project pressure
    val = 0.0;
    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/
    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          // different than velocity since there are nodes on
          // the boundary of the mesh cell i
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;

	  val0 = 0;
	  // loop over the degrees of freedom
	  for (iq=0;iq<8;iq++)
	  {
	      ii = dofp[iq]+N_vec_u;
	      //OutPut(ii << " ");
	      val0 +=  rhs_sub[ii];
	  }
	  val += val0/8.0;
	}
      }
    }

    val /= N_sub3;

    //OutPut(" new " << val << endl);
    for (ii=0;ii<N_nodes_p3;ii++)
    {
	rhs_sub[ii+N_vec_u] -= val;
    }
    //SolveMultipleSystems(a_sub, rhs_sub, N_all, N_all, 1, 1);
    /*if (i<4)
    {
	OutPut(" sol " << sqrt(Ddot(3*N_nodes_u3,rhs_sub,rhs_sub)/(3*N_nodes_u3)) << endl);
        OutPut("Solution:"<<endl);
        for(ii=3*N_nodes_u3;ii<N_all;ii++)
          OutPut("sol("<<ii<<") = "<<rhs_sub[ii]<<endl);
	  }*/
	//if (i==1)
	// exit(1);
    
    //    OutPut(" sol " << rhs_sub[N_all-1]<<"   "<< rhs_sub[N_all-2]<<"   "<<sqrt(Ddot(3*N_nodes_u3,rhs_sub,rhs_sub)) << endl);
    //    OutPut("a " << i << endl);

    // write out the coordinates of the nodes
    // and the values of the solution (u1, u2, u3, p) in the nodes

    /*OutPut("Node_u  coords                              u1             u2             u3"<<endl);

    for(ii=0;ii<N_nodes_u3;ii++)
      OutPut(ii<<'\t' << "("<<setw(7) << setprecision(5)<<coord_x[ii]<<",   "<<setw(7) << setprecision(5)<<coord_y[ii]<<",   "<<setw(7) << setprecision(5)<<coord_z[ii]<<")"<<setw(15) << setprecision(5)<<
      rhs_sub[ii]<<setw(15) << setprecision(5) <<rhs_sub[ii+N_nodes_u3]<<setw(15) << setprecision(5) <<rhs_sub[ii+2*N_nodes_u3]<<endl);

    OutPut("Node_p  coords                              p"<<endl);
    for(ii=0;ii<N_nodes_p3;ii++)
      OutPut(ii <<'\t' << "("<<setw(7) << setprecision(5)<< coord_x[N_nodes_u3+ii]<<",   "<<setw(7) << setprecision(5)<<coord_y[N_nodes_u3+ii]<<",   "<<setw(7) << setprecision(5)<<coord_z[N_nodes_u3+ii]<<")"<<setw(12) << setprecision(5)<<rhs_sub[ii+3*N_nodes_u3]<<endl);

    // OutPut("N_nodes_u3:  "<<N_nodes_u3<<"   N_nodes_p3:  "<<N_nodes_p3<<endl);

    exit(4711);*/
    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);

    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];

    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);

    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();

    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();

    test_value = new double[4*N_];
    test_value_x = test_value + N_;
    test_value_y = test_value_x + N_;
    test_value_z = test_value_y + N_;
    integral_1 = new double[3 * N_];
    memset(integral_1,0,3*N_*SizeOfDouble);
    integral_2 = integral_1 + N_;
    integral_3 = integral_2 + N_;

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          // compute coordinates of submesh cell vertices
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);

          // for the pressure (Q1)
          xp[0] = xc[0];
          yp[0] = yc[0];
          zp[0] = zc[0];

          xp[1] = xc[2];
          yp[1] = yc[2];
          zp[1] = zc[2];

          xp[2] = xc[8];
          yp[2] = yc[8];
          zp[2] = zc[8];

          xp[3] = xc[6];
          yp[3] = yc[6];
          zp[3] = zc[6];

          xp[4] = xc[18];
          yp[4] = yc[18];
          zp[4] = zc[18];

          xp[5] = xc[20];
          yp[5] = yc[20];
          zp[5] = zc[20];

          xp[6] = xc[26];
          yp[6] = yc[26];
          zp[6] = zc[26];

          xp[7] = xc[24];
          yp[7] = yc[24];
          zp[7] = zc[24];

          // coordinates on the edges - for the velocity
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[16] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[16] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[16] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (2*ic-1)*N_nodes_u*N_nodes_u+ (2*jc-1)*N_nodes_u+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes_u;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes_u;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes_u*N_nodes_u;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes_u;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes_u;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes_u*N_nodes_u;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes_u;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes_u;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;

          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;

          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;

          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;

          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;

          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;

          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;

          }

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            index = 27*jj;
            a[index] = 1;
            a[index+1] = xc[jj];
            a[index+2] = yc[jj];
            a[index+3] = zc[jj];
            a[index+4] = xc[jj]*yc[jj];
            a[index+5] = xc[jj]*zc[jj];
            a[index+6] = yc[jj]*zc[jj];
            a[index+7] = xc[jj]*yc[jj]*zc[jj];
            a[index+8] = xc[jj]*xc[jj];
            a[index+9] = yc[jj]*yc[jj];
            a[index+10] = zc[jj]*zc[jj];
            a[index+11] = a[index+8]*yc[jj];
            a[index+12] = a[index+8]*zc[jj];
            a[index+13] = xc[jj]*a[index+9];
            a[index+14] = a[index+9]*zc[jj];
            a[index+15] = xc[jj]*a[index+10];
            a[index+16] = yc[jj]*a[index+10];
            a[index+17] = a[index+11]*zc[jj];
            a[index+18] = a[index+13]*zc[jj];
            a[index+19] = xc[jj]*a[index+16];
            a[index+20] = a[index+8]*a[index+9];
            a[index+21] = a[index+8]*a[index+10];
            a[index+22] = a[index+9]*a[index+10];
            a[index+23] = a[index+20]*zc[jj];
            a[index+24] = a[index+21]*yc[jj];
            a[index+25] = xc[jj]*a[index+22];
            a[index+26] = xc[jj]*a[index+25];

            if (dof[jj]!=-1)
            {
              index1=dof[jj];
              b[jj] =  rhs_sub[index1];
              b[jj+27] = rhs_sub[index1+N_nodes_u3];
              b[jj+54] = rhs_sub[index1+2*N_nodes_u3];
            }
            else
              b[jj] = b[jj+27] = b[jj+54] = 0;
          }
          //  for (jj=0;jj<27*27;jj++) OutPut("i="<<i<<"  ic="<<ic<<"  jc="<<jc<<"  kc="<<kc<<"  a2("<<jj<<")= "<<a[jj]<<endl);
          // for (ii=0;ii<27;ii++)
          //for (jj=0;jj<27;jj++)
          //OutPut("a2("<<ii+1<<","<< jj+1<<")= "<<a[ii*27+jj]<<endl);
          // solve system for the coefficients of the bilinear function
          //        for (jj=0;jj<27;jj++)
          //           OutPut("b(" <<jj << ")= " << b[jj] << endl);

          SolveMultipleSystemsLapack(a,b,27,27,27,3);

          //         for (jj=0;jj<27;jj++)
          //           OutPut("n " <<jj << " " << b[jj] << endl);

          for (jj=0;jj<8;jj++)
          {
            index=8*jj;
            ap[index] = 1;
            ap[index+1] = xp[jj];
            ap[index+2] = yp[jj];
            ap[index+3] = zp[jj];
            ap[index+4] = xp[jj]*yp[jj];
            ap[index+5] = xp[jj]*zp[jj];
            ap[index+6] = yp[jj]*zp[jj];
            ap[index+7] = xp[jj]*yp[jj]*zp[jj];
            index1=dofp[jj];
            bp[jj] =  rhs_sub[index1+N_vec_u];
          }

          SolveLinearSystemTranspose(ap,bp,8,8);

          for (iq = 0;iq < quad_points; iq++)
          {
            // quadrature points
            // ONLY FOR PARALLELEPIPED !!!
            xq = ((xc[2]-xc[0])*qx[iq]+(xc[8]-xc[2])*qy[iq]+(xc[26]-xc[8])*qz[iq])/2+xc[13];
            yq = ((yc[2]-yc[0])*qx[iq]+(yc[8]-yc[2])*qy[iq]+(yc[26]-yc[8])*qz[iq])/2+yc[13];
            zq = ((zc[2]-zc[0])*qx[iq]+(zc[8]-zc[2])*qy[iq]+(zc[26]-zc[8])*qz[iq])/2+zc[13];

            //OutPut(xq << " " << yq << " " << zq << endl);

            temp=detJK*weight[iq]*tau;
            temp1=temp*theta1;
            // compute convection field in the quadrature point
            u1->FindGradientLocal(cell,i,xq,yq,zq,u1_values);
            u2->FindGradientLocal(cell,i,xq,yq,zq,u2_values);
            u3->FindGradientLocal(cell,i,xq,yq,zq,u3_values);

            // compute the value of the test functions and its derivatives
            //in the quadrature points
            TFEDatabase3D::GetRefFromOrig(RefTrans, xq, yq, zq, xi, eta, zeta);
            bf->GetDerivatives(D000, xi, eta, zeta, test_value);
            bf->GetDerivatives(D100, xi, eta, zeta, test_value_x);
            bf->GetDerivatives(D010, xi, eta, zeta, test_value_y);
            bf->GetDerivatives(D001, xi, eta, zeta, test_value_z);

            //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
            //       test_value[2]  << " " << test_value[3] << endl);

            // compute the derivative of the solution of the RFB problem
            // in the barycenter of the submesh cell
            // compute for this purpose the bilinear function on the mesh cell

            // compute value and gradient of tilde_u in the quadrature points
            Compute_Q2_Value_Gradient(b,xq,yq,zq,rhs1);
            Compute_Q2_Value_Gradient(b+27,xq,yq,zq,rhs2);
            Compute_Q2_Value_Gradient(b+54,xq,yq,zq,rhs3);

            // compute value and gradient of tilde_p in the quadrature points
            Compute_Q1_Value_RFB(bp, xq, yq, zq, &rhsp);

            // update the integrals
            for (j=0;j<N_;j++)
            {
              // add convection term
              val = (u1_values[0] * rhs1[1] + u2_values[0] * rhs1[2] + u3_values[0] * rhs1[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * u1_values[1]  + rhs2[0] * u1_values[2] + rhs3[0] * u1_values[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * rhs1[1] + rhs2[0] * rhs1[2] + rhs3[0] * rhs1[3] )
		  *test_value[j]*temp1;
              //           OutPut(val << " rhs " << rhs1[1] << " " << rhs1[2]<<  " " << rhs1[3]  << endl );
              //add laplacian
              val += (rhs1[1]*test_value_x[j] + rhs1[2]*test_value_y[j] + rhs1[3]*test_value_z[j]) *temp1*eps;
              //            OutPut(val << " ");
              // add term from time derivative
              val += rhs1[0]*test_value[j]*weight[iq]*detJK;
              //           OutPut(val << " ");
              // substract pressure term
              val -= temp*rhsp*test_value_x[j];
              integral_1[j] += val;
              //           OutPut(val << endl);

              val = (u1_values[0] * rhs2[1] + u2_values[0] * rhs2[2] + u3_values[0] * rhs2[3])
                *test_value[j]*temp1;
              val += (rhs1[0] * u2_values[1]  + rhs2[0] * u2_values[2] + rhs3[0] * u2_values[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * rhs2[1] + rhs2[0] * rhs2[2] + rhs3[0] * rhs2[3] )
		  *test_value[j]*temp1;
              //add laplacian
              val += (rhs2[1]*test_value_x[j] + rhs2[2]*test_value_y[j] + rhs2[3]*test_value_z[j]) *temp1*eps;
              // add term from time derivative
              val += rhs2[0]*test_value[j]*weight[iq]*detJK;
              // substract pressure term
              val -= temp*rhsp*test_value_y[j];
              integral_2[j] += val;

              val = (u1_values[0] * rhs3[1] + u2_values[0] * rhs3[2] + u3_values[0] * rhs3[3])
                *test_value[j]*temp1;
              val += (rhs1[0] * u3_values[1]  + rhs2[0] * u3_values[2] + rhs3[0] * u3_values[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * rhs3[1] + rhs2[0] * rhs3[2] + rhs3[0] * rhs3[3] )
		  *test_value[j]*temp1;
              //add laplacian
              val += (rhs3[1]*test_value_x[j] + rhs3[2]*test_value_y[j] + rhs3[3]*test_value_z[j]) *temp1*eps;
              // add term from time derivative
              val += rhs3[0]*test_value[j]*weight[iq]*detJK;
              // substract pressure term
              val -= temp*rhsp*test_value_z[j];
              integral_3[j] += val;

            }                    // end j
          }                      // end iq
        }                        // kc
      }                          // jc
    }                            // ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      rhs[gdof] -= integral_1[j];
      rhs[gdof+N_U] -= integral_2[j];
      rhs[gdof+2*N_U] -= integral_3[j];
    }
    delete test_value;
    delete integral_1;
  }
  delete a_sub;
  delete rhs_sub;
  delete coeff;
  //  for (i=0;i<1000;i++)
  //  OutPut("rhs("<<i<<") = " << rhs[i] << endl);
}

/***************************************************************************************/
//
// bubble VMS method with solution of coupled system on the mesh cells
//
// velocity is approximated by Q_2 bubbles
// pressure is approximated by Q_1
// on dof of the pressure is set to be zero
//
/***************************************************************************************/

void ApproximateTimeRFB_coupled_cn_SolutionQuad_Q2_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *old_small_scales, double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, kk, N_U, N_, iq;
  int dof[27], dofp[8];
  int N_sub3, N_nodes_u, N_nodes_u3, N_vec_u, N_nodes_p, N_nodes_p3, N_all, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof, index1;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[27],y[27], z[27], a[729],ap[64], b[729], bp[64],*a_sub, *rhs_sub;
  double fe_fct_p;
  double xc[27], yc[27], zc[27], xp[8], yp[8], zp[8];
  double area, val, val0;
  double u1_values[4], u2_values[4], u3_values[4], p_values[4];
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value, *test_value_x, *test_value_y, *test_value_z;
  double delta, hK, mu_T, temp, temp1;
  double val_test[4], val_ansatz[4], xq, yq, zq, detJK;
  double rhs1[4], rhs2[4], rhs3[4], rhsp, valU[3], gradU[9];
  double old_small_u1[27], old_small_u2[27], old_small_u3[27], old_small_p[8];
  double coeff_old_small_u1[81], *coeff_old_small_u2, *coeff_old_small_u3;
  double old_small_u1_qp[4], old_small_u2_qp[4], old_small_u3_qp[4], old_small_p_qp;
  //double *a0, *a1, *a2, *b0, *b1, *b2;
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  //double tau = 1;
  double theta1 = TDatabase::TimeDB->THETA1;
  double theta = 0.5, value_global[9];
  value_global[0] = value_global[1] = value_global[2] = 0;
  value_global[3] = value_global[4] = value_global[5] = 0;
  value_global[6] = value_global[7] = value_global[8] = 0;


  /*************************************************************************/
  /* data for quadrature rule gauss2                                       */
  /*************************************************************************/
  /*int quad_points = 8;
  double weight[8]={ 1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0 };
  double qx[8]={-0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489,  0.5773502691896257645091489};
  double qy[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};
  double qz[8]={-0.5773502691896257645091489, -0.5773502691896257645091489,
               -0.5773502691896257645091489, -0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489,
                0.5773502691896257645091489,  0.5773502691896257645091489};*/
  /*************************************************************************/
  /* data for quadrature rule gauss3                                       */
  /*************************************************************************/
  int quad_points = 27;
  double oneweight[3]={ 0.5555555555555555555555558,
                        0.8888888888888888888888888,
                        0.5555555555555555555555558 };
  double onexi[3]={ -0.7745966692414833770358530,
                     0,
                     0.7745966692414833770358530 };

  double weight[27], qx[27], qy[27],qz[27];
  for(ii=0;ii<3;ii++)
    for(jj=0;jj<3;jj++)
      for(kk=0;kk<3;kk++)
      {
        i=9*ii+3*jj+kk;
        weight[i]=oneweight[ii]*oneweight[jj]*oneweight[kk];
        qx[i]=onexi[kk];
        qy[i]=onexi[jj];
        qz[i]=onexi[ii];
      }
  
  /*************************************************************************/
  /* data for quadrature rule gauss4                                       */
  /*************************************************************************/
/*  int quad_points = 64;
  double oneweight[4]=
  {
    0.652145154862546142626936051,
    0.347854845137453857373063949,
    0.347854845137453857373063949,
    0.652145154862546142626936051
  };

  double onexi[4]=
  {
    0.33998104358485626480266576,
    0.86113631159405257522394649,
    -0.86113631159405257522394649,
    -0.33998104358485626480266576
  };

  double weight[64], qx[64], qy[64],qz[64];
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      for(kk=0;kk<4;kk++)
      {
        i=16*ii+4*jj+kk;
        weight[i]=oneweight[ii]*oneweight[jj]*oneweight[kk];
        qx[i]=onexi[kk];
        qy[i]=onexi[jj];
        qz[i]=onexi[ii];
      }
*/
  /*************************************************************************/
  /* data for quadrature rule gauss5                                       */
  /*************************************************************************/

  /*int quad_points = 125;
  double onexi[]={
                0,
                0.538469310105683091036314421,
                0.906179845938663992797626878,
               -0.906179845938663992797626878,
               -0.538469310105683091036314421
             }; 
  double oneweight[]={
                0.568888888888888888888888889,
                0.478628670499366468041291515,
                0.236926885056189087514264041,
                0.236926885056189087514264041,
                0.478628670499366468041291515
             };

  double weight[125], qx[125], qy[125], qz[125];

  int l;
  for(int i=0;i<5;i++)
    for(int j=0;j<5;j++)
      for(int k=0;k<5;k++)
      {
        l=25*i+5*j+k;
        weight[l]=oneweight[i]*oneweight[j]*oneweight[k];
        qx[l]=onexi[k];
        qy[l]=onexi[j];
        qz[l]=onexi[i];
      }

  */
  OutPut("rfb_coupled"<<endl);

  // compute number of local mesh cells
  N_sub3 = N_sub * N_sub * N_sub;
  // compute number of dof for the local problems
  // do not consider velocity nodes on the boundary of the mesh cell
  // for Q2
  N_nodes_u = 2*N_sub-1;
  N_nodes_u3 = N_nodes_u*N_nodes_u*N_nodes_u;

  // there is no boundary condition for p, for Q1
  N_nodes_p = N_sub+1;
  N_nodes_p3 = N_nodes_p*N_nodes_p*N_nodes_p;
  N_vec_u = 3*N_nodes_u3;
  N_all = N_vec_u+N_nodes_p3;

  //OutPut("u " << N_nodes_u << " u3 " << N_nodes_u3  << " p " << N_nodes_p3 << " all " << N_all << endl);

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is 2*N_sub-1 for the velocity (Q2)
  // and N_sub-1 for the pressure (Q1)

  a_sub = new double[N_all*N_all];
  memset(a_sub, 0, N_all*N_all*SizeOfDouble);

  rhs_sub = new double[N_all];

  coeff = new double[4];

  coeff_old_small_u2 = coeff_old_small_u1+27;
  coeff_old_small_u3 = coeff_old_small_u2+27;  
  // write the coordinates of the nodes into the vectors
  // coord: first N_nodes_u3 components correspond to the
  // velocity nodes, and the last N_nodes_p3, to the pressure nodes
  /*double coord_x[N_nodes_u3+N_nodes_p3];
  double coord_y[N_nodes_u3+N_nodes_p3];
  double coord_z[N_nodes_u3+N_nodes_p3];*/

  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
      /*OutPut("cell " << i << " ");
    jj = i*N_all;
    for (ii=0; ii < N_all; ii++)
    {
	OutPut(old_small_scales[jj+ii] << " ");
    }
    OutPut(endl);*/
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //  if (i==379)
      // OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }

    // The unit cube:
    /*
    x[0] = -6; y[0] = 1; z[0] = 0.8;
    x[1] = -4; y[1] = 1; z[1] = 0.8;
    x[2] = -4; y[2] = 1.5; z[2] = 0.8;
    x[3] = -6; y[3] = 1.5; z[3] = 0.8;
    x[4] = -6; y[4] = 1; z[4] = 1;
    x[5] = -4; y[5] = 1; z[5] = 1;
    x[6] = -4; y[6] = 1.5; z[6] = 1;
    x[7] = -6; y[7] = 1.5; z[7] = 1;*/
    /*    
    x[0] = 0; y[0] = 0; z[0] = 0;
    x[1] = 1; y[1] = 0; z[1] = 0;
    x[2] = 1; y[2] = 1; z[2] = 0;
    x[3] = 0; y[3] = 1; z[3] = 0;
    x[4] = 0; y[4] = 0; z[4] = 1;
    x[5] = 1; y[5] = 0; z[5] = 1;
    x[6] = 1; y[6] = 1; z[6] = 1;
    x[7] = 0; y[7] = 1; z[7] = 1;*/
    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    detJK = area/8.0;

    // compute delta for stabilization
    /* if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)
      hK = cell->GetShortestEdge();
      else*/
      hK = cell->GetDiameter();
    // measure of the submesh cells
    hK /= N_sub;
    delta = CharacteristicFilterWidth(hK);

    memset(a_sub, 0, N_all*N_all*SizeOfDouble);
    memset(rhs_sub, 0, N_all*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/
    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell vertices
          // for the velocity
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);
          // positions for the pressure dof (Q1)
          xp[0] = xc[0];
          yp[0] = yc[0];
          zp[0] = zc[0];
          //OutPut("press " << xp[0] << " " << yp[0] << " " << zp[0] << endl);
          xp[1] = xc[2];
          yp[1] = yc[2];
          zp[1] = zc[2];

          xp[2] = xc[8];
          yp[2] = yc[8];
          zp[2] = zc[8];

          xp[3] = xc[6];
          yp[3] = yc[6];
          zp[3] = zc[6];

          xp[4] = xc[18];
          yp[4] = yc[18];
          zp[4] = zc[18];

          xp[5] = xc[20];
          yp[5] = yc[20];
          zp[5] = zc[20];

          xp[6] = xc[26];
          yp[6] = yc[26];
          zp[6] = zc[26];

          xp[7] = xc[24];
          yp[7] = yc[24];
          zp[7] = zc[24];

          // coordinates on the edges - for the velocity
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[16] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[16] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[16] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          // initialize rhs
          memset(b,0,729*SizeOfDouble);
          memset(bp,0,64*SizeOfDouble);
          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            index = 27*jj;
            a[index] = 1;
            a[index+1] = xc[jj];
            a[index+2] = yc[jj];
            a[index+3] = zc[jj];
            a[index+4] = xc[jj]*yc[jj];
            a[index+5] = xc[jj]*zc[jj];
            a[index+6] = yc[jj]*zc[jj];
            a[index+7] = xc[jj]*yc[jj]*zc[jj];
            a[index+8] = xc[jj]*xc[jj];
            a[index+9] = yc[jj]*yc[jj];
            a[index+10] = zc[jj]*zc[jj];
            a[index+11] = a[index+8]*yc[jj];
            a[index+12] = a[index+8]*zc[jj];
            a[index+13] = xc[jj]*a[index+9];
            a[index+14] = a[index+9]*zc[jj];
            a[index+15] = xc[jj]*a[index+10];
            a[index+16] = yc[jj]*a[index+10];
            a[index+17] = a[index+11]*zc[jj];
            a[index+18] = a[index+13]*zc[jj];
            a[index+19] = xc[jj]*a[index+16];
            a[index+20] = a[index+8]*a[index+9];
            a[index+21] = a[index+8]*a[index+10];
            a[index+22] = a[index+9]*a[index+10];
            a[index+23] = a[index+20]*zc[jj];
            a[index+24] = a[index+21]*yc[jj];
            a[index+25] = xc[jj]*a[index+22];
            a[index+26] = xc[jj]*a[index+25];
            b[28*jj] = 1;
          }
          //for (jj=0;jj<27*27;jj++) OutPut("i="<<i<<"  ic="<<ic<<"  jc="<<jc<<"  kc="<<kc<<"  a1("<<jj<<")= "<<a[jj]<<endl);
          // for (ii=0;ii<27;ii++)
          // for (jj=0;jj<27;jj++)
          // OutPut("a1("<<ii+1<<","<< jj+1<<")= "<<a[ii*27+jj]<<endl);
          // set matrices for computation of the coefficients
          // of the trilinear function
          for (jj=0;jj<8;jj++)
          {
            index=8*jj;
            ap[index] = 1;
            ap[index+1] = xp[jj];
            ap[index+2] = yp[jj];
            ap[index+3] = zp[jj];
            ap[index+4] = xp[jj]*yp[jj];
            ap[index+5] = xp[jj]*zp[jj];
            ap[index+6] = yp[jj]*zp[jj];
            ap[index+7] = xp[jj]*yp[jj]*zp[jj];
            bp[9*jj] = 1;
          }

          // solve systems for the coefficients

          //SolveMultipleSystems(ap,bp,8,8,8,8);
          //SolveMultipleSystems(a,b,27,27,27,27);

          SolveMultipleSystemsLapack(ap,bp,8,8,8,8);
          SolveMultipleSystemsLapack(a,b,27,27,27,27);
          //         for (jj=0;jj<27;jj++)
          //           OutPut("n0 " <<jj << " " << b[jj] << endl);

          /*for (ii=0;ii<27;ii++)
              for (jj=0;jj<27;jj++)
               {Compute_Q2_Value_Gradient(b+(27*ii),xc[jj],yc[jj],zc[jj],val_test);
                  OutPut(ii << " " << jj << " "<< val_test[0] << endl);
                }
          exit(1);*/

          // compute local velocity dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (2*ic-1)*N_nodes_u*N_nodes_u+ (2*jc-1)*N_nodes_u+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes_u;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes_u;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes_u*N_nodes_u;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes_u;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes_u;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes_u*N_nodes_u;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes_u;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes_u;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          /*OutPut("vor ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          // different than velocity since there are nodes on
          // the boundary of the mesh cell i
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;
	  
	  jj = i * N_all + N_vec_u;
	  for (ii=0; ii < 8;ii++)
	  {
	      old_small_p[ii] = old_small_scales[jj+ dofp[ii]];
	      //OutPut(old_small_p[ii] << " ");
	  }
	  //OutPut(endl);


          //OutPut("pdof " << dofp[0] << " "<< dofp[1] << " "<< dofp[2] << " "<< dofp[3] << " "
          //	<< dofp[4] << " "<< dofp[5] << " "<< dofp[6] << " "<< dofp[7] << endl);

          // correct (for the velocity) for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
          }
          /*OutPut("nacha ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }

          /*OutPut("nachb ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
          }

          /*OutPut("nachc ");
           for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;
          }

          /*OutPut("nachd ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;
          }

          /*OutPut("nache ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;
          }
	  // position of the mesh cell i in the array of the bubble unknowns
	  //if (i==379)
	  // OutPut("cell " << i);
	  jj = i * N_all;
	  for (ii=0; ii < 27;ii++)
	  {
	      if (dof[ii]!=-1)
	      {
		  old_small_u1[ii] = old_small_scales[jj + dof[ii]];
		  old_small_u2[ii] = old_small_scales[jj + N_nodes_u3 + dof[ii]];
		  old_small_u3[ii] = old_small_scales[jj + 2*N_nodes_u3 + dof[ii]];
		  //OutPut(old_small_u3[ii] << " ");
		  //  if (i==379)
		  //  OutPut(" " << old_small_u3[ii]);
	      }
	      else
	      {
		  old_small_u1[ii] = old_small_u2[ii] = old_small_u3[ii] = 0;
	      }
	  }
	  // if (i==379)
	  //  OutPut(endl);
	  
	  //OutPut(endl);
          /*OutPut("nachf ");
          for (ii=0;ii<27;ii++)
          OutPut(dof[ii] << " ");
          OutPut(endl);*/

          // write the coordinates of the nodes into the vectors
          // coord: first N_nodes_u3 components correspond to the
          // velocity nodes, and the last N_nodes_p3, to the pressure nodes
          /*for (ii=0;ii<27;ii++)
          {
             if (dof[ii] != -1)
             {
               index = dof[ii];
               coord_x[index] = xc[ii];
               coord_y[index] = yc[ii];
               coord_z[index] = zc[ii];
             }
          }
          for (ii=0;ii<8;ii++)
          {
            index = dofp[ii];
            coord_x[N_nodes_u3+index] = xp[ii];
            coord_y[N_nodes_u3+index] = yp[ii];
            coord_z[N_nodes_u3+index] = zp[ii];
          }*/

	   // compute value of old small scales in quadrature point
	   // set matrices for computation of the coefficients
	   // of the bilinear function
	   for (jj=0;jj<27;jj++)
	   {
	       index = 27*jj;
	       a[index] = 1;
	       a[index+1] = xc[jj];
	       a[index+2] = yc[jj];
	       a[index+3] = zc[jj];
	       a[index+4] = xc[jj]*yc[jj];
	       a[index+5] = xc[jj]*zc[jj];
	       a[index+6] = yc[jj]*zc[jj];
	       a[index+7] = xc[jj]*yc[jj]*zc[jj];
	       a[index+8] = xc[jj]*xc[jj];
	       a[index+9] = yc[jj]*yc[jj];
	       a[index+10] = zc[jj]*zc[jj];
	       a[index+11] = a[index+8]*yc[jj];
	       a[index+12] = a[index+8]*zc[jj];
	       a[index+13] = xc[jj]*a[index+9];
	       a[index+14] = a[index+9]*zc[jj];
	       a[index+15] = xc[jj]*a[index+10];
	       a[index+16] = yc[jj]*a[index+10];
	       a[index+17] = a[index+11]*zc[jj];
	       a[index+18] = a[index+13]*zc[jj];
	       a[index+19] = xc[jj]*a[index+16];
	       a[index+20] = a[index+8]*a[index+9];
	       a[index+21] = a[index+8]*a[index+10];
	       a[index+22] = a[index+9]*a[index+10];
	       a[index+23] = a[index+20]*zc[jj];
	       a[index+24] = a[index+21]*yc[jj];
	       a[index+25] = xc[jj]*a[index+22];
	       a[index+26] = xc[jj]*a[index+25];
	       coeff_old_small_u1[jj] = old_small_u1[jj];
	       coeff_old_small_u1[jj+27] = old_small_u2[jj];
	       coeff_old_small_u1[jj+54] = old_small_u3[jj];	    
	   }
	   
	   SolveMultipleSystemsLapack(a,coeff_old_small_u1,27,27,27,3);

          // loop over the quadrature points
          for (iq = 0;iq < quad_points; iq++)
          {
            temp1=detJK*weight[iq];
            temp=temp1*tau;

            // quadrature points
            // ONLY FOR PARALLELEPIPED !!!
            xq = ((xc[2]-xc[0])*qx[iq]+(xc[8]-xc[2])*qy[iq]+(xc[26]-xc[8])*qz[iq])/2+xc[13];
            yq = ((yc[2]-yc[0])*qx[iq]+(yc[8]-yc[2])*qy[iq]+(yc[26]-yc[8])*qz[iq])/2+yc[13];
            zq = ((zc[2]-zc[0])*qx[iq]+(zc[8]-zc[2])*qy[iq]+(zc[26]-zc[8])*qz[iq])/2+zc[13];

            // compute convection field in the quadrature point
            u1->FindGradientLocal(cell,i,xq,yq,zq,u1_values);
            u2->FindGradientLocal(cell,i,xq,yq,zq,u2_values);
            u3->FindGradientLocal(cell,i,xq,yq,zq,u3_values);
            p->FindGradientLocal(cell,i,xq,yq,zq,p_values);
	    
            Compute_Q2_Value_Gradient(coeff_old_small_u1,xq,yq,zq,old_small_u1_qp);
            Compute_Q2_Value_Gradient(coeff_old_small_u2,xq,yq,zq,old_small_u2_qp);
            Compute_Q2_Value_Gradient(coeff_old_small_u3,xq,yq,zq,old_small_u3_qp);
	    
	    //OutPut(old_small_u1_qp[0] << " " << old_small_u2_qp[0] << " "<<
	    //	   old_small_u3_qp[0] << " " << old_small_p_qp << endl);
	    /*
            u1_values[0] = u2_values[0] = u3_values[0]  = 0;
            u1_values[1] = u2_values[1] = u3_values[1]  = 0;
            u1_values[2] = u2_values[2] = u3_values[2]  = 0;
            u1_values[3] = u2_values[3] = u3_values[3]  = 0;
            p_values[0] = 0;*/
            // compute value of rhs in quadrature point
            Coeffs(1, &xq, &yq, &zq, NULL, &coeff);
            //OutPut("coeff " << coeff[1] << " " << coeff[2] << endl);
            //coeff[1] = 1; coeff[2] = 0; coeff[3] = 0;
            // prepare computation of turbulent viscosity
            valU[0]=u1_values[0];
            valU[1]=u2_values[0];
            valU[2]=u3_values[0];
            for(ii=0;ii<3;ii++)
            {
              gradU[3*ii]=u1_values[ii+1];
              gradU[3*ii+1]=u2_values[ii+1];
              gradU[3*ii+2]=u3_values[ii+1];
            }
            // turbulent viscosity
	    valU[0] += old_small_u1[0];	    
	    valU[1] += old_small_u2[0];	    
	    valU[2] += old_small_u3[0];	 
	    gradU[0] += old_small_u1[1];
	    gradU[1] += old_small_u2[1];
	    gradU[2] += old_small_u3[1];
	    gradU[3] += old_small_u1[2];
	    gradU[4] += old_small_u2[2];
	    gradU[5] += old_small_u3[2];
	    gradU[6] += old_small_u1[3];
	    gradU[7] += old_small_u2[3];
	    gradU[8] += old_small_u3[3];
            mu_T = TurbulentViscosity3D(delta,gradU,valU,NULL,NULL,NULL,&zq,-4711);
	    valU[0] -= old_small_u1[0];	    
	    valU[1] -= old_small_u2[0];	    
	    valU[2] -= old_small_u3[0];	    
	    gradU[0] -= old_small_u1[1];
	    gradU[1] -= old_small_u2[1];
	    gradU[2] -= old_small_u3[1];
	    gradU[3] -= old_small_u1[2];
	    gradU[4] -= old_small_u2[2];
	    gradU[5] -= old_small_u3[2];
	    gradU[6] -= old_small_u1[3];
	    gradU[7] -= old_small_u2[3];
	    gradU[8] -= old_small_u3[3];
	    //mu_T = 0;
	    //eps = 1;
            // update rhs and matrix entries
            // the matrix is stored column wise

            // ii -- test function
            for(ii=0;ii<27;ii++)
            {
              // Dirichlet value
              if (dof[ii] == -1)
                continue;
              index=dof[ii];
              // get values of test function
              Compute_Q2_Value_Gradient(b+(27*ii),xq,yq,zq,val_test);
              //OutPut(val_test[0] << " " << val_test[1] << " " << val_test[2] << " " << val_test[3] << " " << endl);

              // compute first rhs (u1)
              // convective term
              val = (u1_values[0]+old_small_u1_qp[0])*(old_small_u1_qp[1])
		  + (u2_values[0]+old_small_u2_qp[0])*(old_small_u1_qp[2]) 
		  + (u3_values[0]+old_small_u3_qp[0])*(old_small_u1_qp[3]);
              val *= val_test[0] * temp*(1-theta);
              rhs_sub[index] -= val;
              val = (u1_values[0]+old_small_u1_qp[0])*(u1_values[1])
		  + (u2_values[0]+old_small_u2_qp[0])*(u1_values[2]) 
		  + (u3_values[0]+old_small_u3_qp[0])*(u1_values[3]);
              val *= val_test[0] * temp;
              rhs_sub[index] -= val;
              // diffusive term
              rhs_sub[index] -= eps*temp*(old_small_u1_qp[1]*val_test[1]
	      		  +old_small_u1_qp[2]*val_test[2]+old_small_u1_qp[3]*val_test[3])*(1-theta);
              rhs_sub[index] -= eps*temp*(u1_values[1]*val_test[1]
	      		  +u1_values[2]*val_test[2]+u1_values[3]*val_test[3]);
              // add (f1,phi_ii)
	      // is is assumed that the rhs of the former time is the same as of the 
	      // present time
              rhs_sub[index] += temp*coeff[1]*val_test[0];
              // pressure term
	      rhs_sub[index] += temp*p_values[0]*val_test[1];
	      // solution from former time
	      //if (TDatabase::TimeDB->CURRENTTIME > 0.2)
	      rhs_sub[index] += temp1 * old_small_u1_qp[0] * val_test[0];
             
	      // compute second rhs (u2)
              // convective term
              val = (u1_values[0]+old_small_u1_qp[0])*(old_small_u2_qp[1])
		  + (u2_values[0]+old_small_u2_qp[0])*(old_small_u2_qp[2])
		  + (u3_values[0]+old_small_u3_qp[0])*(old_small_u2_qp[3]);
	      val *= val_test[0] * temp * (1-theta);
              rhs_sub[index+ N_nodes_u3] -= val;
              val = (u1_values[0]+old_small_u1_qp[0])*(u2_values[1])
		  + (u2_values[0]+old_small_u2_qp[0])*(u2_values[2]) 
		  + (u3_values[0]+old_small_u3_qp[0])*(u2_values[3]);
              val *= val_test[0] * temp;
              rhs_sub[index + N_nodes_u3] -= val;
              // diffusive term
              rhs_sub[index+ N_nodes_u3] -= eps*temp*(old_small_u2_qp[1]*val_test[1]
	      		  +old_small_u2_qp[2]*val_test[2]+old_small_u2_qp[3]*val_test[3])*(1-theta);
              rhs_sub[index+ N_nodes_u3] -= eps*temp*(u2_values[1]*val_test[1]
	      		  +u2_values[2]*val_test[2]+u2_values[3]*val_test[3]);
              // add (f1,phi_ii)
	      // is is assumed that the rhs of the former time is the same as of the 
	      // present time
              rhs_sub[index+ N_nodes_u3] += temp*coeff[2]*val_test[0];
              // pressure term
	      rhs_sub[index+ N_nodes_u3] += temp*p_values[0]*val_test[2];
	      // solution from former time
	      rhs_sub[index+ N_nodes_u3] += temp1 * old_small_u2_qp[0] * val_test[0];

              // compute third rhs (u3)
              // convective term
              val = (u1_values[0]+old_small_u1_qp[0])* (old_small_u3_qp[1])
		  + (u2_values[0]+old_small_u2_qp[0])* (old_small_u3_qp[2])
		  + (u3_values[0]+old_small_u3_qp[0])* (old_small_u3_qp[3]);
              val *= val_test[0] * temp * (1-theta);
              rhs_sub[index+ 2*N_nodes_u3] -= val;
              val = (u1_values[0]+old_small_u1_qp[0])*(u3_values[1])
		  + (u2_values[0]+old_small_u2_qp[0])*(u3_values[2]) 
		  + (u3_values[0]+old_small_u3_qp[0])*(u3_values[3]);
              val *= val_test[0] * temp;
              rhs_sub[index + 2 * N_nodes_u3] -= val;
              // diffusive term
              rhs_sub[index+ 2*N_nodes_u3] -= eps*temp*(old_small_u3_qp[1]*val_test[1]
	     				  +old_small_u3_qp[2]*val_test[2]+old_small_u3_qp[3]*val_test[3])
		  *(1-theta);
              rhs_sub[index+ 2*N_nodes_u3] -= eps*temp*(u3_values[1]*val_test[1]
	      				  + u3_values[2]*val_test[2]+ u3_values[3]*val_test[3]);
              // add (f3,phi_ii)
	      // is is assumed that the rhs of the former time is the same as of the 
	      // present time
              rhs_sub[index+ 2*N_nodes_u3] += temp*coeff[3]*val_test[0];
              // pressure term
	      rhs_sub[index+ 2*N_nodes_u3] += temp*p_values[0]*val_test[3];
	      // solution from former time
	      rhs_sub[index+ 2*N_nodes_u3] += temp1 * old_small_u3_qp[0] * val_test[0];

              // jj -- ansatz function
              for (jj=0;jj<27;jj++)
              {
                if (dof[jj] == -1)
                  continue;
                // entry (jj,ii) exists
                // compute index in array containing the entry

                // A11
                index = dof[ii]*N_all+dof[jj];
                // OutPut(" ii " << ii << " jj " << jj << " ind " << index << endl);
                // values of ansatz function in quadrature points
                Compute_Q2_Value_Gradient(b+(27*jj),xq,yq,zq,val_ansatz);
                // add Laplacian
                val = val_ansatz[1]*val_test[1] + val_ansatz[2]*val_test[2] +val_ansatz[3]*val_test[3];
                //val *= eps*temp;
                val *= (eps+mu_T) *temp * theta;
		a_sub[index] += val;

                // add convective term
                // convection times gradient of ansatz fct.
		val = (u1_values[0]+old_small_u1_qp[0])* val_ansatz[1]
		    + (u2_values[0]+old_small_u2_qp[0])* val_ansatz[2]
		    + (u3_values[0]+old_small_u3_qp[0])* val_ansatz[3];
                // multiplied with test fct.
                val *= val_test[0]*temp * theta;
		a_sub[index] += val;

                // add mass term
                // A11
                val = val_ansatz[0] * val_test[0];
                a_sub[index] += val * temp1;

                // copy to the other diagonal blocks
                // A22
                index1 = index + (N_all+1)*N_nodes_u3;
                a_sub[index1] = a_sub[index];

                // A33
                index1 = index + 2*(N_all+1)*N_nodes_u3;
                a_sub[index1] = a_sub[index];
              }                  // end jj

              // right upper pressure matrix
              for (kk=0;kk<8;kk++)
              {
                index = dof[ii]*N_all+N_vec_u+dofp[kk];

                /*if (dof[ii]==13)
                OutPut(dofp[kk] << " " << index << " " << 13*N_all+N_vec_u-1 << " " << 14*N_all << endl);*/

                // value of pressure test function in bary center
                Compute_Q1_Value_RFB(bp+(8*kk), xq, yq, zq, &fe_fct_p);

                // B1

                /*if ((index>=13*N_all+N_vec_u)&&(index<14*N_all)&&((dofp[kk]==12)||(dofp[kk]==12)))
                         OutPut("vor " << index << " a " << a_sub[index] << " u " << dof[ii] << " p "
                  << dofp[kk] << endl);*/

                val =  fe_fct_p * val_test[1];
                a_sub[index] -= val*temp;

                /*if ((index>=13*N_all+N_vec_u)&&(index<14*N_all)&&((dofp[kk]==12)||(dofp[kk]==12)))
                         OutPut("nach " << index << " a " << a_sub[index] << " u " << dof[ii] << " p "
                  << dofp[kk] << endl);*/

                // B2
                index += N_all*N_nodes_u3;
                val = fe_fct_p * val_test[2];
                a_sub[index] -= val*temp;

                // B3
                index += N_all*N_nodes_u3;
                val = fe_fct_p * val_test[3];
                a_sub[index] -= val*temp;
              }                  //end kk (ansatz pressure)
            }                    // end ii (test velocity)

            // left lower pressure block
            // test functions of pressure
            for (kk=0;kk<8;kk++) // test functions of pressure
            {
              // value of pressure test function in quadrature point
              Compute_Q1_Value_RFB(bp+(8*kk), xq, yq, zq, &fe_fct_p);
              // compute forth rhs (continuity eq.)
              val = (u1_values[1] + u2_values[2] + u3_values[3]) * fe_fct_p;
              index=dofp[kk];
              //rhs_sub[index + N_vec_u] -= val*temp1;
              // ansatz functions of velocity
              for(jj=0;jj<27;jj++)
              {
                if (dof[jj] == -1)
                  continue;
                index = (N_vec_u+dofp[kk])*N_all+dof[jj];

                // value of pressure test function in quadrature points

                //OutPut(jj << " ind " << index << endl);

                // values of ansatz function in quadrature points
                Compute_Q2_Value_Gradient(b+(27*jj),xq,yq,zq,val_ansatz);
                // B1^T
                val = val_ansatz[1] * fe_fct_p;
                a_sub[index] += val*temp1;
                // B2^T
                index += N_nodes_u3;
                val = val_ansatz[2] * fe_fct_p;
                a_sub[index] += val*temp1;
                // B3^T
                index += N_nodes_u3;
                val = val_ansatz[3] * fe_fct_p;
                a_sub[index] += val*temp1;
              }                  // end jj
            }                    //end kk
          }                      // end iq
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    /* for (ic=0;ic<N_all;ic++)
         for (jc=0;jc<N_all;jc++)
           OutPut("a("<<ic<<","<<jc<<") = " << a_sub[ic*N_all+jc] << endl);*/

    // replace the last line of the matrix on the left hand side
    // with (0, 0, ..., 0, 1), and the last element of the vector
    // on the right hand side with 0

    for(ii=(N_all-1)*N_all;ii<N_all*N_all-1;ii++)
      a_sub[ii]=0;
    a_sub[N_all*N_all-1]=1;
    rhs_sub[N_all-1]=0;

    //	   OutPut("mat: " << sqrt(Ddot(N_all*N_all,a_sub,a_sub)));
/*       for(ii=0;ii<N_nodes_u3;ii++)
	   OutPut("rhs("<<ii+1<<") = "<<rhs_sub[ii]<<"; "
		  "rhs("<<N_nodes_u3+ii+1<<") = "<<rhs_sub[N_nodes_u3+ii]<<"; "
		  "rhs("<<2*N_nodes_u3+ii+1<<") = "<<rhs_sub[2+N_nodes_u3+ii]<<";"
		  <<endl);
*/  
/*
    // write out the components of the matrices on the left hand side
    
     for(ii=0;ii<N_all;ii++)
     {
       for(jj=0;jj<N_all;jj++)
	   OutPut("a("<<ii+1<<", "<<jj+1<<")=  "<<a_sub[ii*N_all+jj]<<";"<<endl);
       OutPut(endl);
     }
*/   

    /*for(ii=0;ii<N_all;ii++)
       for(jj=0;jj<N_all;jj++)
       {
    if (fabs(a_sub[ii*N_all+jj] - a_sub[jj*N_all+ii]) > 1e-12)
    {
          OutPut(setw(7) << ii << " " << jj << " " << setprecision(5)<<a_sub[ii*N_all+jj]<<"   ");
          OutPut(setw(7) << setprecision(5)<<a_sub[jj*N_all+ii]<<endl);
    }
       }*/

    /* for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii)*N_all+3*N_nodes_u3+jj])>1e-12)
          OutPut("b0("<<ii<<", "<<jj<<")=  "<<a_sub[(ii)*N_all+N_vec_u+jj]<< endl);

     for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii+N_nodes_u3)*N_all+N_vec_u+jj])>1e-12)
         OutPut("b1("<<ii<<", "<<jj<<")=  "<<a_sub[(ii+N_nodes_u3)*N_all+N_vec_u+jj]<<endl);

     for(ii=0;ii<N_nodes_u3;ii++)
       for(jj=0;jj<N_nodes_p3;jj++)
    if (fabs(a_sub[(ii+2*N_nodes_u3)*N_all+N_vec_u+jj])>1e-12)
         OutPut("b2("<<ii<<", "<<jj<<")=  "<<a_sub[(ii+2*N_nodes_u3)*N_all+N_vec_u+jj]<<endl);

    //   OutPut(" mat " << sqrt(Ddot(N_all*N_all,a_sub,a_sub)));

        for(ii=0;ii<3*N_nodes_u3+N_nodes_p3;ii++)
        {
      val = 0;
        for(jj=0;jj<3*N_nodes_u3+N_nodes_p3;jj++)
      val += a_sub[ii*(3*N_nodes_u3+N_nodes_p3)+jj] *sol_sub[jj];
      rhs_sub[ii] = val - rhs_sub[ii];
      if (fabs(rhs_sub[ii])>1e-12)
      OutPut("rhs " << ii << " " << rhs_sub[ii] << endl);
    }
      exit(1);*/

    /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    //OutPut(" rhs " << rhs_sub[N_all-1]<<"   "<< rhs_sub[N_all-2]<<"   "<<sqrt(Ddot(N_all,rhs_sub,rhs_sub)) << endl);
    //SolveLinearSystemLapack(a_sub, rhs_sub, N_all, N_all);
       value_global[3] += Ddot(N_nodes_u3,rhs_sub,rhs_sub);
       value_global[4] += Ddot(N_nodes_u3,rhs_sub+N_nodes_u3,rhs_sub+N_nodes_u3);
       value_global[5] += Ddot(N_nodes_u3,rhs_sub+2*N_nodes_u3,rhs_sub+2*N_nodes_u3);
    SolveLinearSystemTranspose(a_sub, rhs_sub, N_all, N_all);
       for(ii=0;ii<N_nodes_u3;ii++)
	   /*  OutPut("sol("<<ii+1<<") = "<<rhs_sub[ii]<<"; "
		  "sol("<<N_nodes_u3+ii+1<<") = "<<rhs_sub[N_nodes_u3+ii]<<"; "
		  "sol("<<2*N_nodes_u3+ii+1<<") = "<<rhs_sub[2*N_nodes_u3+ii]<<";"
		  <<endl);*/
       value_global[6] += Ddot(N_nodes_u3,rhs_sub,rhs_sub);
       value_global[7] += Ddot(N_nodes_u3,rhs_sub+N_nodes_u3,rhs_sub+N_nodes_u3);
       value_global[8] += Ddot(N_nodes_u3,rhs_sub+2*N_nodes_u3,rhs_sub+2*N_nodes_u3);


    // project pressure
    val = 0.0;
    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/
    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          // different than velocity since there are nodes on
          // the boundary of the mesh cell i
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;

	  val0 = 0;
	  // loop over the degrees of freedom
	  for (iq=0;iq<8;iq++)
	  {
	      ii = dofp[iq]+N_vec_u;
	      //OutPut(ii << " ");
	      val0 +=  rhs_sub[ii];
	  }
	  val += val0/8.0;
	}
      }
    }

    val /= N_sub3;

    //OutPut(" new " << val << endl);
    for (ii=0;ii<N_nodes_p3;ii++)
    {
	rhs_sub[ii+N_vec_u] -= val;
    }
       
    //SolveMultipleSystems(a_sub, rhs_sub, N_all, N_all, 1, 1);

    /* for(ii=0;ii<N_all;ii++)
	   OutPut("sol("<<ii+1<<") = "<<rhs_sub[ii]<<";"<<endl);
	   exit(1);*/
    /*if (i<4)
    {
	OutPut(" sol " << sqrt(Ddot(3*N_nodes_u3,rhs_sub,rhs_sub)/(3*N_nodes_u3)) << endl);
        OutPut("Solution:"<<endl);
        for(ii=3*N_nodes_u3;ii<N_all;ii++)
          OutPut("sol("<<ii<<") = "<<rhs_sub[ii]<<endl);
	  }*/
	//if (i==1)
	// exit(1);
    
    //    OutPut(" sol " << rhs_sub[N_all-1]<<"   "<< rhs_sub[N_all-2]<<"   "<<sqrt(Ddot(3*N_nodes_u3,rhs_sub,rhs_sub)) << endl);
    //    OutPut("a " << i << endl);

    // write out the coordinates of the nodes
    // and the values of the solution (u1, u2, u3, p) in the nodes

    /*OutPut("Node_u  coords                              u1             u2             u3"<<endl);

    for(ii=0;ii<N_nodes_u3;ii++)
      OutPut(ii<<'\t' << "("<<setw(7) << setprecision(5)<<coord_x[ii]<<",   "<<setw(7) << setprecision(5)<<coord_y[ii]<<",   "<<setw(7) << setprecision(5)<<coord_z[ii]<<")"<<setw(15) << setprecision(5)<<
      rhs_sub[ii]<<setw(15) << setprecision(5) <<rhs_sub[ii+N_nodes_u3]<<setw(15) << setprecision(5) <<rhs_sub[ii+2*N_nodes_u3]<<endl);

    OutPut("Node_p  coords                              p"<<endl);
    for(ii=0;ii<N_nodes_p3;ii++)
      OutPut(ii <<'\t' << "("<<setw(7) << setprecision(5)<< coord_x[N_nodes_u3+ii]<<",   "<<setw(7) << setprecision(5)<<coord_y[N_nodes_u3+ii]<<",   "<<setw(7) << setprecision(5)<<coord_z[N_nodes_u3+ii]<<")"<<setw(12) << setprecision(5)<<rhs_sub[ii+3*N_nodes_u3]<<endl);

    // OutPut("N_nodes_u3:  "<<N_nodes_u3<<"   N_nodes_p3:  "<<N_nodes_p3<<endl);

    exit(4711);*/
    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);

    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];

    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);

    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();

    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();

    test_value = new double[4*N_];
    test_value_x = test_value + N_;
    test_value_y = test_value_x + N_;
    test_value_z = test_value_y + N_;
    integral_1 = new double[3 * N_];
    memset(integral_1,0,3*N_*SizeOfDouble);
    integral_2 = integral_1 + N_;
    integral_3 = integral_2 + N_;

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          // compute coordinates of submesh cell vertices
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[8], &yc[8], &zc[8]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[18], &yc[18], &zc[18]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[20], &yc[20], &zc[20]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[24], &yc[24], &zc[24]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[26], &yc[26], &zc[26]);

          // for the pressure (Q1)
          xp[0] = xc[0];
          yp[0] = yc[0];
          zp[0] = zc[0];

          xp[1] = xc[2];
          yp[1] = yc[2];
          zp[1] = zc[2];

          xp[2] = xc[8];
          yp[2] = yc[8];
          zp[2] = zc[8];

          xp[3] = xc[6];
          yp[3] = yc[6];
          zp[3] = zc[6];

          xp[4] = xc[18];
          yp[4] = yc[18];
          zp[4] = zc[18];

          xp[5] = xc[20];
          yp[5] = yc[20];
          zp[5] = zc[20];

          xp[6] = xc[26];
          yp[6] = yc[26];
          zp[6] = zc[26];

          xp[7] = xc[24];
          yp[7] = yc[24];
          zp[7] = zc[24];

          // coordinates on the edges - for the velocity
          xc[1] = (xc[0]+xc[2])/2;
          yc[1] = (yc[0]+yc[2])/2;
          zc[1] = (zc[0]+zc[2])/2;

          xc[3] = (xc[0]+xc[6])/2;
          yc[3] = (yc[0]+yc[6])/2;
          zc[3] = (zc[0]+zc[6])/2;

          xc[5] = (xc[2]+xc[8])/2;
          yc[5] = (yc[2]+yc[8])/2;
          zc[5] = (zc[2]+zc[8])/2;

          xc[7] = (xc[6]+xc[8])/2;
          yc[7] = (yc[6]+yc[8])/2;
          zc[7] = (zc[6]+zc[8])/2;

          xc[9] = (xc[0]+xc[18])/2;
          yc[9] = (yc[0]+yc[18])/2;
          zc[9] = (zc[0]+zc[18])/2;

          xc[11] = (xc[2]+xc[20])/2;
          yc[11] = (yc[2]+yc[20])/2;
          zc[11] = (zc[2]+zc[20])/2;

          xc[15] = (xc[6]+xc[24])/2;
          yc[15] = (yc[6]+yc[24])/2;
          zc[15] = (zc[6]+zc[24])/2;

          xc[17] = (xc[8]+xc[26])/2;
          yc[17] = (yc[8]+yc[26])/2;
          zc[17] = (zc[8]+zc[26])/2;

          xc[19] = (xc[18]+xc[20])/2;
          yc[19] = (yc[18]+yc[20])/2;
          zc[19] = (zc[18]+zc[20])/2;

          xc[21] = (xc[18]+xc[24])/2;
          yc[21] = (yc[18]+yc[24])/2;
          zc[21] = (zc[18]+zc[24])/2;

          xc[23] = (xc[20]+xc[26])/2;
          yc[23] = (yc[20]+yc[26])/2;
          zc[23] = (zc[20]+zc[26])/2;

          xc[25] = (xc[24]+xc[26])/2;
          yc[25] = (yc[24]+yc[26])/2;
          zc[25] = (zc[24]+zc[26])/2;

          // coordinates in the barycenters
          xc[4] = (xc[0]+xc[2]+xc[6]+xc[8])/4;
          yc[4] = (yc[0]+yc[2]+yc[6]+yc[8])/4;
          zc[4] = (zc[0]+zc[2]+zc[6]+zc[8])/4;

          xc[10] = (xc[0]+xc[2]+xc[18]+xc[20])/4;
          yc[10] = (yc[0]+yc[2]+yc[18]+yc[20])/4;
          zc[10] = (zc[0]+zc[2]+zc[18]+zc[20])/4;

          xc[12] = (xc[0]+xc[6]+xc[18]+xc[24])/4;
          yc[12] = (yc[0]+yc[6]+yc[18]+yc[24])/4;
          zc[12] = (zc[0]+zc[6]+zc[18]+zc[24])/4;

          xc[14] = (xc[2]+xc[8]+xc[20]+xc[26])/4;
          yc[14] = (yc[2]+yc[8]+yc[20]+yc[26])/4;
          zc[14] = (zc[2]+zc[8]+zc[20]+zc[26])/4;

          xc[16] = (xc[6]+xc[8]+xc[24]+xc[26])/4;
          yc[16] = (yc[6]+yc[8]+yc[24]+yc[26])/4;
          zc[16] = (zc[6]+zc[8]+zc[24]+zc[26])/4;

          xc[22] = (xc[18]+xc[20]+xc[24]+xc[26])/4;
          yc[22] = (yc[18]+yc[20]+yc[24]+yc[26])/4;
          zc[22] = (zc[18]+zc[20]+zc[24]+zc[26])/4;

          xc[13] = (xc[3]+xc[5]+xc[21]+xc[23])/4;
          yc[13] = (yc[3]+yc[5]+yc[21]+yc[23])/4;
          zc[13] = (zc[3]+zc[5]+zc[21]+zc[23])/4;

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (2*ic-1)*N_nodes_u*N_nodes_u+ (2*jc-1)*N_nodes_u+2*kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + 1;
          dof[3] = dof[0] + N_nodes_u;
          dof[4] = dof[3] + 1;
          dof[5] = dof[4] + 1;
          dof[6] = dof[3] + N_nodes_u;
          dof[7] = dof[6] + 1;
          dof[8] = dof[7] + 1;
          dof[9] = dof[0] + N_nodes_u*N_nodes_u;
          dof[10] = dof[9] + 1;
          dof[11] = dof[10] + 1;
          dof[12] = dof[9] + N_nodes_u;
          dof[13] = dof[12] + 1;
          dof[14] = dof[13] + 1;
          dof[15] = dof[12] + N_nodes_u;
          dof[16] = dof[15] + 1;
          dof[17] = dof[16] + 1;
          dof[18] = dof[9] + N_nodes_u*N_nodes_u;
          dof[19] = dof[18] + 1;
          dof[20] = dof[19] + 1;
          dof[21] = dof[18] + N_nodes_u;
          dof[22] = dof[21] + 1;
          dof[23] = dof[22] + 1;
          dof[24] = dof[21] + N_nodes_u;
          dof[25] = dof[24] + 1;
          dof[26] = dof[25] + 1;

          // compute local pressure dof at this mesh cell
          // compute dof of the `left lower' vertex
          dofp[0] = ic*N_nodes_p*N_nodes_p+ jc*N_nodes_p+kc;
          dofp[1] = dofp[0] + 1;
          dofp[2] = dofp[1] + N_nodes_p;
          dofp[3] = dofp[0] + N_nodes_p;
          dofp[4] = dofp[0] + N_nodes_p*N_nodes_p;
          dofp[5] = dofp[1] + N_nodes_p*N_nodes_p;
          dofp[6] = dofp[2] + N_nodes_p*N_nodes_p;
          dofp[7] = dofp[3] + N_nodes_p*N_nodes_p;

          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;

          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;
            dof[21] = -1;
            dof[22] = -1;
            dof[23] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;

          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[9] = -1;
            dof[10] = -1;
            dof[11] = -1;
            dof[18] = -1;
            dof[19] = -1;
            dof[20] = -1;

          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[6] = -1;
            dof[7] = -1;
            dof[8] = -1;
            dof[15] = -1;
            dof[16] = -1;
            dof[17] = -1;
            dof[24] = -1;
            dof[25] = -1;
            dof[26] = -1;

          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[9] = -1;
            dof[12] = -1;
            dof[15] = -1;
            dof[18] = -1;
            dof[21] = -1;
            dof[24] = -1;

          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[2] = -1;
            dof[5] = -1;
            dof[8] = -1;
            dof[11] = -1;
            dof[14] = -1;
            dof[17] = -1;
            dof[20] = -1;
            dof[23] = -1;
            dof[26] = -1;

          }

	  // this is now the new solution
	  jj = i * N_all;
	  for (ii=0; ii < 27;ii++)
	  {
	      if (dof[ii]!=-1)
	      {
		  old_small_u1[ii] = old_small_scales[jj + dof[ii]];
		  old_small_u2[ii] = old_small_scales[jj + N_nodes_u3 + dof[ii]];
		  old_small_u3[ii] = old_small_scales[jj + 2*N_nodes_u3 + dof[ii]];
		  //OutPut(old_small_u3[ii] << " ");
	      }
	      else
	      {
		  old_small_u1[ii] = old_small_u2[ii] = old_small_u3[ii] = 0;
	      }
	  }

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<27;jj++)
          {
            index = 27*jj;
            a[index] = 1;
            a[index+1] = xc[jj];
            a[index+2] = yc[jj];
            a[index+3] = zc[jj];
            a[index+4] = xc[jj]*yc[jj];
            a[index+5] = xc[jj]*zc[jj];
            a[index+6] = yc[jj]*zc[jj];
            a[index+7] = xc[jj]*yc[jj]*zc[jj];
            a[index+8] = xc[jj]*xc[jj];
            a[index+9] = yc[jj]*yc[jj];
            a[index+10] = zc[jj]*zc[jj];
            a[index+11] = a[index+8]*yc[jj];
            a[index+12] = a[index+8]*zc[jj];
            a[index+13] = xc[jj]*a[index+9];
            a[index+14] = a[index+9]*zc[jj];
            a[index+15] = xc[jj]*a[index+10];
            a[index+16] = yc[jj]*a[index+10];
            a[index+17] = a[index+11]*zc[jj];
            a[index+18] = a[index+13]*zc[jj];
            a[index+19] = xc[jj]*a[index+16];
            a[index+20] = a[index+8]*a[index+9];
            a[index+21] = a[index+8]*a[index+10];
            a[index+22] = a[index+9]*a[index+10];
            a[index+23] = a[index+20]*zc[jj];
            a[index+24] = a[index+21]*yc[jj];
            a[index+25] = xc[jj]*a[index+22];
            a[index+26] = xc[jj]*a[index+25];

            if (dof[jj]!=-1)
            {
              index1=dof[jj];
              b[jj] =  rhs_sub[index1];
              b[jj+27] = rhs_sub[index1+N_nodes_u3];
              b[jj+54] = rhs_sub[index1+2*N_nodes_u3];
            }
            else
              b[jj] = b[jj+27] = b[jj+54] = 0;
          }
          //  for (jj=0;jj<27*27;jj++) OutPut("i="<<i<<"  ic="<<ic<<"  jc="<<jc<<"  kc="<<kc<<"  a2("<<jj<<")= "<<a[jj]<<endl);
          // for (ii=0;ii<27;ii++)
          //for (jj=0;jj<27;jj++)
          //OutPut("a2("<<ii+1<<","<< jj+1<<")= "<<a[ii*27+jj]<<endl);
          // solve system for the coefficients of the bilinear function
          //        for (jj=0;jj<27;jj++)
          //           OutPut("b(" <<jj << ")= " << b[jj] << endl);

          SolveMultipleSystemsLapack(a,b,27,27,27,3);

          //         for (jj=0;jj<27;jj++)
          //           OutPut("n " <<jj << " " << b[jj] << endl);

          for (jj=0;jj<8;jj++)
          {
            index=8*jj;
            ap[index] = 1;
            ap[index+1] = xp[jj];
            ap[index+2] = yp[jj];
            ap[index+3] = zp[jj];
            ap[index+4] = xp[jj]*yp[jj];
            ap[index+5] = xp[jj]*zp[jj];
            ap[index+6] = yp[jj]*zp[jj];
            ap[index+7] = xp[jj]*yp[jj]*zp[jj];
            index1=dofp[jj];
            bp[jj] =  rhs_sub[index1+N_vec_u];
          }

          SolveLinearSystemTranspose(ap,bp,8,8);

	   // compute value of old small scales in quadrature point
	   // set matrices for computation of the coefficients
	   // of the bilinear function
	   for (jj=0;jj<27;jj++)
	   {
	       index = 27*jj;
	       a[index] = 1;
	       a[index+1] = xc[jj];
	       a[index+2] = yc[jj];
	       a[index+3] = zc[jj];
	       a[index+4] = xc[jj]*yc[jj];
	       a[index+5] = xc[jj]*zc[jj];
	       a[index+6] = yc[jj]*zc[jj];
	       a[index+7] = xc[jj]*yc[jj]*zc[jj];
	       a[index+8] = xc[jj]*xc[jj];
	       a[index+9] = yc[jj]*yc[jj];
	       a[index+10] = zc[jj]*zc[jj];
	       a[index+11] = a[index+8]*yc[jj];
	       a[index+12] = a[index+8]*zc[jj];
	       a[index+13] = xc[jj]*a[index+9];
	       a[index+14] = a[index+9]*zc[jj];
	       a[index+15] = xc[jj]*a[index+10];
	       a[index+16] = yc[jj]*a[index+10];
	       a[index+17] = a[index+11]*zc[jj];
	       a[index+18] = a[index+13]*zc[jj];
	       a[index+19] = xc[jj]*a[index+16];
	       a[index+20] = a[index+8]*a[index+9];
	       a[index+21] = a[index+8]*a[index+10];
	       a[index+22] = a[index+9]*a[index+10];
	       a[index+23] = a[index+20]*zc[jj];
	       a[index+24] = a[index+21]*yc[jj];
	       a[index+25] = xc[jj]*a[index+22];
	       a[index+26] = xc[jj]*a[index+25];
	       coeff_old_small_u1[jj] = old_small_u1[jj];
	       coeff_old_small_u1[jj+27] = old_small_u2[jj];
	       coeff_old_small_u1[jj+54] = old_small_u3[jj];	    
	   }
	   
	   SolveMultipleSystemsLapack(a,coeff_old_small_u1,27,27,27,3);

          for (iq = 0;iq < quad_points; iq++)
          {
            // quadrature points
            // ONLY FOR PARALLELEPIPED !!!
            xq = ((xc[2]-xc[0])*qx[iq]+(xc[8]-xc[2])*qy[iq]+(xc[26]-xc[8])*qz[iq])/2+xc[13];
            yq = ((yc[2]-yc[0])*qx[iq]+(yc[8]-yc[2])*qy[iq]+(yc[26]-yc[8])*qz[iq])/2+yc[13];
            zq = ((zc[2]-zc[0])*qx[iq]+(zc[8]-zc[2])*qy[iq]+(zc[26]-zc[8])*qz[iq])/2+zc[13];

            //OutPut(xq << " " << yq << " " << zq << endl);

            temp=detJK*weight[iq]*tau;
            temp1=temp*theta1;
            // compute convection field in the quadrature point
            u1->FindGradientLocal(cell,i,xq,yq,zq,u1_values);
            u2->FindGradientLocal(cell,i,xq,yq,zq,u2_values);
            u3->FindGradientLocal(cell,i,xq,yq,zq,u3_values);

	    // compute old bubble solution in quadrature point
	    Compute_Q2_Value_Gradient(coeff_old_small_u1,xq,yq,zq,old_small_u1_qp);
            Compute_Q2_Value_Gradient(coeff_old_small_u2,xq,yq,zq,old_small_u2_qp);
            Compute_Q2_Value_Gradient(coeff_old_small_u3,xq,yq,zq,old_small_u3_qp);
	    
            // compute the value of the test functions and its derivatives
            //in the quadrature points
            TFEDatabase3D::GetRefFromOrig(RefTrans, xq, yq, zq, xi, eta, zeta);
            bf->GetDerivatives(D000, xi, eta, zeta, test_value);
            bf->GetDerivatives(D100, xi, eta, zeta, test_value_x);
            bf->GetDerivatives(D010, xi, eta, zeta, test_value_y);
            bf->GetDerivatives(D001, xi, eta, zeta, test_value_z);

            //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
            //       test_value[2]  << " " << test_value[3] << endl);

            // compute the derivative of the solution of the RFB problem
            // in the barycenter of the submesh cell
            // compute for this purpose the bilinear function on the mesh cell

            // compute value and gradient of tilde_u in the quadrature points
            Compute_Q2_Value_Gradient(b,xq,yq,zq,rhs1);
            Compute_Q2_Value_Gradient(b+27,xq,yq,zq,rhs2);
            Compute_Q2_Value_Gradient(b+54,xq,yq,zq,rhs3);

            // compute value and gradient of tilde_p in the quadrature points
            Compute_Q1_Value_RFB(bp, xq, yq, zq, &rhsp);

            // update the integrals
            for (j=0;j<N_;j++)
            {
              // add convection term
              val = (u1_values[0] * rhs1[1] + u2_values[0] * rhs1[2] + u3_values[0] * rhs1[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * rhs1[1] + rhs2[0] * rhs1[2] + rhs3[0] * rhs1[3] )
                *test_value[j]*temp1;
              val += (rhs1[0] * u1_values[1]  + rhs2[0] * u1_values[2] + rhs3[0] * u1_values[3] )
                *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * old_small_u1_qp[1] + old_small_u2_qp[0] * old_small_u1_qp[2] 
		      + old_small_u3_qp[0] * old_small_u1_qp[3] )
		  *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * u1_values[1] + old_small_u2_qp[0] * u1_values[2]
		      + old_small_u3_qp[0] * u1_values[3] )
		      *test_value[j]*temp1;
	      val += (u1_values[0] * old_small_u1_qp[1]  + u2_values[0] * old_small_u1_qp[2] 
		      + u3_values[0] * old_small_u1_qp[3] )
		  *test_value[j]*area*temp1;
              //           OutPut(val << " rhs " << rhs1[1] << " " << rhs1[2]<<  " " << rhs1[3]  << endl );
              //add laplacian
              val += (rhs1[1]*test_value_x[j] + rhs1[2]*test_value_y[j] + rhs1[3]*test_value_z[j]) *temp1*eps;
              //            OutPut(val << " ");
	      val += (old_small_u1_qp[1]*test_value_x[j] + old_small_u1_qp[2]*test_value_y[j] 
		      + old_small_u1_qp[3]*test_value_z[j]) *temp1*eps; 
              // add term from time derivative
              val += rhs1[0]*test_value[j]*weight[iq]*detJK;
              //           OutPut(val << " ");
              val -= old_small_u1_qp[0]*test_value[j]*weight[iq]*detJK;
              // substract pressure term
              val -= temp*rhsp*test_value_x[j];
              integral_1[j] += val;
              //           OutPut(val << endl);

              // add convection term
              val = (u1_values[0] * rhs2[1] + u2_values[0] * rhs2[2] + u3_values[0] * rhs2[3])
		  *test_value[j]*temp1;
              val += (rhs1[0] * rhs2[1] + rhs2[0] * rhs2[2] + rhs3[0] * rhs2[3] )
		  *test_value[j]*temp1;
              val += (rhs1[0] * u2_values[1]  + rhs2[0] * u2_values[2] + rhs3[0] * u2_values[3] )
                *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * old_small_u2_qp[1] + old_small_u2_qp[0] * old_small_u2_qp[2] 
		      + old_small_u3_qp[0] * old_small_u2_qp[3] )
		  *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * u2_values[1] + old_small_u2_qp[0] * u2_values[2]
		      + old_small_u3_qp[0] * u2_values[3] )
		      *test_value[j]*temp1;
	      val += (u1_values[0] * old_small_u2_qp[1]  + u2_values[0] * old_small_u2_qp[2] 
		      + u3_values[0] * old_small_u2_qp[3] )
		  *test_value[j]*area*temp1;	      
              //add laplacian
              val += (rhs2[1]*test_value_x[j] + rhs2[2]*test_value_y[j] + rhs2[3]*test_value_z[j]) *temp1*eps;
	      val += (old_small_u2_qp[1]*test_value_x[j] + old_small_u2_qp[2]*test_value_y[j] 
		      + old_small_u2_qp[3]*test_value_z[j]) *temp1*eps; 
              // add term from time derivative
              val += rhs2[0]*test_value[j]*weight[iq]*detJK;
              val -= old_small_u2_qp[0]*test_value[j]*weight[iq]*detJK;
              // substract pressure term
              val -= temp*rhsp*test_value_y[j];
              integral_2[j] += val;

	      // convective term
              val = (u1_values[0] * rhs3[1] + u2_values[0] * rhs3[2] + u3_values[0] * rhs3[3])
                *test_value[j]*temp1;
              val += (rhs1[0] * rhs3[1] + rhs2[0] * rhs3[2] + rhs3[0] * rhs3[3] )
		  *test_value[j]*temp1;
              val += (rhs1[0] * u3_values[1]  + rhs2[0] * u3_values[2] + rhs3[0] * u3_values[3] )
                *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * old_small_u3_qp[1] + old_small_u2_qp[0] * old_small_u3_qp[2] 
		      + old_small_u3_qp[0] * old_small_u3_qp[3] )
		  *test_value[j]*temp1;
              val += (old_small_u1_qp[0] * u3_values[1] + old_small_u2_qp[0] * u3_values[2]
		      + old_small_u3_qp[0] * u3_values[3] )
		      *test_value[j]*temp1;
	      val += (u1_values[0] * old_small_u3_qp[1]  + u2_values[0] * old_small_u3_qp[2] 
		      + u3_values[0] * old_small_u3_qp[3] )
		  *test_value[j]*area*temp1;
              //add laplacian
              val += (rhs3[1]*test_value_x[j] + rhs3[2]*test_value_y[j] + rhs3[3]*test_value_z[j]) *temp1*eps;
	      val += (old_small_u3_qp[1]*test_value_x[j] + old_small_u3_qp[2]*test_value_y[j] 
		      + old_small_u3_qp[3]*test_value_z[j]) *temp1*eps; 
              // add term from time derivative
              val += rhs3[0]*test_value[j]*weight[iq]*detJK;
              val -= old_small_u3_qp[0]*test_value[j]*weight[iq]*detJK;
              // substract pressure term
              val -= temp*rhsp*test_value_z[j];
              integral_3[j] += val;

            }                    // end j
          }                      // end iq
        }                        // kc
      }                          // jc
    }                            // ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      value_global[0] += integral_1[j]*integral_1[j];
      value_global[1] += integral_2[j]*integral_2[j];
      value_global[2] += integral_3[j]*integral_3[j];
      //if (TDatabase::TimeDB->CURRENTTIME > 0.1)
      {
	  rhs[gdof] -= integral_1[j];
	  rhs[gdof+N_U] -= integral_2[j];
	  rhs[gdof+2*N_U] -= integral_3[j];
      }
    }
    delete test_value;
    delete integral_1;


    // if (i==379)
    //OutPut("cell " << i);
    // save  bubble solution for next time
    jj = i*N_all;
    for (ii=0; ii < N_all; ii++)
    {
	old_small_scales[jj+ii] = rhs_sub[ii];
	if (i==379)
	   OutPut(ii << " " << rhs_sub[ii] << endl);
    }    
  }
  delete a_sub;
  delete rhs_sub;
  delete coeff;
  OutPut("coupled value_global " << sqrt(value_global[0]) << " " << sqrt(value_global[1]) 
	 << " "  << sqrt(value_global[2]) << " " << endl);
  OutPut("coupled value_global_rhs " << sqrt(value_global[3]) << " " << sqrt(value_global[4]) 
	 << " "  << sqrt(value_global[5]) << " " << endl);
  OutPut("coupled value_global_sol " << sqrt(value_global[6]) << " " << sqrt(value_global[7]) 
	 << " "  << sqrt(value_global[8]) << " " << endl);

  //  for (i=0;i<1000;i++)
  //  OutPut("rhs("<<i<<") = " << rhs[i] << endl);
}

void ApproximateTimeRFBSolutionQuad_cn_NSE3D(TCollection *Coll, TFEFunction3D *u1,
TFEFunction3D *u2, TFEFunction3D *u3, TFEFunction3D *p, CoeffFct3D *Coeffs,
double *old_small_scales, double *rhs)
{
  int i, j, N_Cells, N_V, ii, jj, ic, jc, kc, dof[8], N_U, N_;
  int N_sub3, N_nodes, N_nodes3, index, *global_numbers;
  int *begin_index, *N_BaseFunct, *globdof, gdof, dof_ii, dof_ii2;
  int N_sub = TDatabase::ParamDB->RFB_SUBMESH_LAYERS;
  double x[8],y[8], z[8], a[64], b[64], *a_sub, *rhs_sub, test_bary;
  double xc[8], yc[8], zc[8], area, val, val1, val2, val3, xs, ys, zs, u1_values[4];
  double u2_values[4], u3_values[4], p_values[4], val1_j, val2_j, val3_j;
  double eps = 1/TDatabase::ParamDB->RE_NR, *coeff;
  double *integral_1, *integral_2, *integral_3, xi, eta, zeta;
  double *test_value, *test_value_x, *test_value_y, *test_value_z;
  double rhs1, rhs2, rhs3, rhs1_x, rhs2_x, rhs3_x, rhs1_y, rhs2_y, rhs3_y, rhs1_z, rhs2_z, rhs3_z;
  double valU[3], gradU[9], delta, hK, mu_T, tau_K, divU;
  double c_1 = 2, c_2 = 1;
  double old_small_u1[8], old_small_u2[8], old_small_u3[8];
  double coeff_old_small_u1[24], *coeff_old_small_u2, *coeff_old_small_u3;
  double old_small_u1_qp[4], old_small_u2_qp[4], old_small_u3_qp[4];
  TBaseCell *cell;
  TVertex *vertex;
  FE3D CurrentElement;
  TFESpace3D *fespace;
  TFE3D *FE_Obj;
  RefTrans3D RefTrans;
  TBaseFunct3D *bf;
  double tau = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double theta1 = TDatabase::TimeDB->THETA1;
  double val0, theta = 0.5, value_global[3];

  OutPut("rfb"<<endl);
  value_global[0] = value_global[1] = value_global[2] = 0;
  N_sub3 = N_sub * N_sub * N_sub;
  N_nodes = N_sub-1;
  N_nodes3 = N_nodes*N_nodes*N_nodes;

  // allocate RFB matrix
  // N_sub = number of internal mesh cells in one direction
  // number of internal nodes in one direction is N_sub-1;

  a_sub = new double[N_nodes3*N_nodes3];
  memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);

  rhs_sub = new double[9 * N_nodes3];
  memset(rhs_sub, 0, 9*N_nodes3*SizeOfDouble);

  coeff = new double[4];

  coeff_old_small_u2 = coeff_old_small_u1+8;
  coeff_old_small_u3 = coeff_old_small_u2+8;  
  // extract information about the fe space and the global d.o.f.
  fespace = u1->GetFESpace3D();
  N_U = u1->GetLength();
  // array with global numbers of d.o.f.
  global_numbers = fespace->GetGlobalNumbers();
  // array which points to the beginning of the global numbers in
  // global_numbers for each mesh cell
  begin_index = fespace->GetBeginIndex();
  // information on the number of basis functions for the available fe
  N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();

  // number of mesh cells
  N_Cells = Coll->GetN_Cells();

  /*************************************************************************/
  /* loop over the mesh cells of the global grid                           */
  /*************************************************************************/
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_V = cell->GetN_Vertices();
    if (N_V < 8)
    {
      OutPut("RFB stabilization only for hexahedral mesh cells implemented !!!"<<endl);
      exit(4711);
    }
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      vertex = cell->GetVertex(j);
      vertex->GetCoords(x[j], y[j], z[j]);
      //OutPut("coords " << x[j] << " " << y[j] << " " << z[j] << endl);
    }
    // The unit cube:
    /*x[0] = -6; y[0] = 1; z[0] = 0.8;
    x[1] = -4; y[1] = 1; z[1] = 0.8;
    x[2] = -4; y[2] = 1.5; z[2] = 0.8;
    x[3] = -6; y[3] = 1.5; z[3] = 0.8;
    x[4] = -6; y[4] = 1; z[4] = 1;
    x[5] = -4; y[5] = 1; z[5] = 1;
    x[6] = -4; y[6] = 1.5; z[6] = 1;
    x[7] = -6; y[7] = 1.5; z[7] = 1;*/
    /*
    x[0] = 0; y[0] = 0; z[0] = 0;
    x[1] = 1; y[1] = 0; z[1] = 0;
    x[2] = 1; y[2] = 1; z[2] = 0;
    x[3] = 0; y[3] = 1; z[3] = 0;
    x[4] = 0; y[4] = 0; z[4] = 1;
    x[5] = 1; y[5] = 0; z[5] = 1;
    x[6] = 1; y[6] = 1; z[6] = 1;
    x[7] = 0; y[7] = 1; z[7] = 1;
    */
    // area of parallelepiped with the scalar triple product
    // det (P4-P0, P3 - P0, P1 - P0)
    area = (x[4] - x[0])*(y[3] - y[0])*(z[1]-z[0]);
    area += (x[3] - x[0])*(y[1] - y[0])*(z[4]-z[0]);
    area += (x[1] - x[0])*(y[4] - y[0])*(z[3]-z[0]);
    area -= (x[1] - x[0])*(y[3] - y[0])*(z[4]-z[0]);
    area -= (x[3] - x[0])*(y[4] - y[0])*(z[1]-z[0]);
    area -= (x[4] - x[0])*(y[1] - y[0])*(z[3]-z[0]);
    area = fabs(area);
    //OutPut("area " << area);
    // area of submesh cells
    area /= N_sub3 ;
    //OutPut(" " << area << endl);
    // compute delta for stabilization

    if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)
      hK = cell->GetShortestEdge();
    else
      hK = cell->GetDiameter();
    // measure of the submesh cells
    hK /= N_sub;
    delta = CharacteristicFilterWidth(hK);

    memset(a_sub, 0, N_nodes3*N_nodes3*SizeOfDouble);
    memset(rhs_sub, 0, 3*N_nodes3*SizeOfDouble);

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell                                                                  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          //OutPut("cell " << ic*N_sub*N_sub+ jc*N_sub + kc << endl);
          // compute coordinates of submesh cell
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          //for (jj=0;jj<8;jj++)
          //  OutPut(jj << " " << xc[jj]  << " " << yc[jj]  << " " << zc[jj] << endl);

          // initialize rhs
          memset(b,0,64*SizeOfDouble);

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            b[9*jj] = 1;
          }
          //for (j=0;j<16;j++)
          //OutPut("#" <<j << " " <<a[j] << " " << b[j] <<endl);

          // solve system for the coefficients of the bilinear function
          //SolveMultipleSystems(a,b,8,8,8,8);
          SolveMultipleSystemsLapack(a,b,8,8,8,8);
          //for (j=0;j<N_V;j++)
          //OutPut("#" <<j << " " << b[4*j] << " " <<  b[4*j+1]<< " " <<
          //b[4*j+2] << " " <<  b[4*j+3] << endl);
	 
          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }
	  /*if (i==379)
	  {
           OutPut(dof[0] << " " << dof[1] << " " << dof[2] << " " << dof[3] << " ");
          OutPut(dof[4] << " " << dof[5] << " " << dof[6] << " " << dof[7] << endl);
	  }*/
	  jj = i * 3 * N_nodes3;
	  for (ii=0; ii < 8;ii++)
	  {
	      if (dof[ii]!=-1)
	      {
		  old_small_u1[ii] = old_small_scales[jj + dof[ii]];
		  old_small_u2[ii] = old_small_scales[jj + N_nodes3 + dof[ii]];
		  old_small_u3[ii] = old_small_scales[jj + 2*N_nodes3 + dof[ii]];
		  //OutPut(old_small_u3[ii] << " ");
		  //  if (i==379)
		  // OutPut(" " << old_small_u1[ii]);
	      }
	      else
	      {
		  old_small_u1[ii] = old_small_u2[ii] = old_small_u3[ii] = 0;
	      }
	  }
	  // if (i==379)
	  //  OutPut(endl);

          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);
          // compute pressure in the center of the submesh cell
          p->FindGradientLocal(cell,i,xs,ys,zs,p_values);

	  Compute_Q1_Value_Gradient_Center_From_Vertices(old_small_u1,old_small_u1_qp);
	  Compute_Q1_Value_Gradient_Center_From_Vertices(old_small_u2,old_small_u2_qp);
	  Compute_Q1_Value_Gradient_Center_From_Vertices(old_small_u3,old_small_u3_qp);

	  // if (i==379)
	  //  OutPut("center " << old_small_u1_qp[0] << " ");
	  

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            coeff_old_small_u1[jj] = old_small_u1[jj];
            coeff_old_small_u1[jj+8] = old_small_u2[jj];
            coeff_old_small_u1[jj+16] = old_small_u3[jj];	    
          }

          SolveMultipleSystemsLapack(a,coeff_old_small_u1,8,8,8,3);

	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u1,xs,ys,zs,old_small_u1_qp);
	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u2,xs,ys,zs,old_small_u2_qp);
	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u3,xs,ys,zs,old_small_u3_qp);
	  //if (i==379)
	  //  OutPut(old_small_u1_qp[0] << " " << old_small_u1_qp[1] <<  " " << old_small_u3_qp[2]<<endl);
	  
          // compute value of rhs in (xs,ys)
          Coeffs(1, &xs, &ys, &zs, NULL, &coeff);
          //OutPut("rhs " << xs << " " << ys << " " << coeff[1] << " " << coeff[2] << endl);
	  /*u1_values[0] = u2_values[0] = u3_values[0]  = 0;
	  u1_values[1] = u2_values[1] = u3_values[1]  = 0;
	  u1_values[2] = u2_values[2] = u3_values[2]  = 0;
	  u1_values[3] = u2_values[3] = u3_values[3]  = 0;
	  p_values[0] = 0;
	  coeff[1] = 1; coeff[2] = 0; coeff[3] = 0;
	  */
          // prepare computation of turbulent viscosity
          valU[0]=u1_values[0];
          valU[1]=u2_values[0];
          valU[2]=u3_values[0];
          for(ii=0;ii<3;ii++)
          {
            gradU[3*ii]=u1_values[ii+1];
            gradU[3*ii+1]=u2_values[ii+1];
            gradU[3*ii+2]=u3_values[ii+1];
          }
          // turbulent viscosity
	    valU[0] += old_small_u1[0];	    
	    valU[1] += old_small_u2[0];	    
	    valU[2] += old_small_u3[0];	 
	    gradU[0] += old_small_u1[1];
	    gradU[1] += old_small_u2[1];
	    gradU[2] += old_small_u3[1];
	    gradU[3] += old_small_u1[2];
	    gradU[4] += old_small_u2[2];
	    gradU[5] += old_small_u3[2];
	    gradU[6] += old_small_u1[3];
	    gradU[7] += old_small_u2[3];
	    gradU[8] += old_small_u3[3];
	    mu_T = TurbulentViscosity3D(delta,gradU,valU,NULL,NULL,NULL,&zs,-4711);
	    valU[0] -= old_small_u1[0];	    
	    valU[1] -= old_small_u2[0];	    
	    valU[2] -= old_small_u3[0];	    
	    gradU[0] -= old_small_u1[1];
	    gradU[1] -= old_small_u2[1];
	    gradU[2] -= old_small_u3[1];
	    gradU[3] -= old_small_u1[2];
	    gradU[4] -= old_small_u2[2];
	    gradU[5] -= old_small_u3[2];
	    gradU[6] -= old_small_u1[3];
	    gradU[7] -= old_small_u2[3];
	    gradU[8] -= old_small_u3[3];	    
	    
          // parameter for div-div term
          tau_K = DivDivStab3D(u1_values[0],u2_values[0],u3_values[0],hK,eps);
          divU = u1_values[1] + u2_values[2] + u3_values[3];

	  /*mu_T = 0.0;
	  eps = 1.0;
	  tau_K = 0.0;
	  */
          // update matrix entries
          // the matrix is stored column wise
          // ii -- test function
          for(ii=0;ii<8;ii++)
          {
            dof_ii = dof[ii];
            if (dof_ii == -1)
              continue;
            // value of test function in bary center
            test_bary = b[8*ii] +  b[8*ii+1] * xs + b[8*ii+2] * ys + b[8*ii+3] * zs
              + b[8*ii+4] * xs * ys +  b[8*ii+5] * xs * zs +  b[8*ii+6] * ys * zs
              + b[8*ii+7] * xs * ys *zs;
            //OutPut("test bary " << test_bary << endl);

            // values of the derivatives of the test function in bary center
            val1 = b[8*ii+1]+ b[8*ii+4]*ys + b[8*ii+5]*zs + b[8*ii+7]*ys*zs;
            val2 = b[8*ii+2]+ b[8*ii+4]*xs + b[8*ii+6]*zs + b[8*ii+7]*xs*zs;
            val3 = b[8*ii+3]+ b[8*ii+5]*xs + b[8*ii+6]*ys + b[8*ii+7]*xs*ys;

            // compute first rhs
            // convection term
            //val = -u1_values[0]*u1_values[1] -  u2_values[0]*u1_values[2] -  u3_values[0]*u1_values[3];
	    val = (u1_values[0]+old_small_u1_qp[0])*(old_small_u1_qp[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(old_small_u1_qp[2])
		+ (u3_values[0]+old_small_u3_qp[0])*(old_small_u1_qp[3]);
	    val = - val;
            val *= test_bary*(1-theta);
	    val0 = (u1_values[0]+old_small_u1_qp[0])*(u1_values[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(u1_values[2]) 
		+ (u3_values[0]+old_small_u3_qp[0])*(u1_values[3]);
	    val0 = -val0;
	    val0 *= test_bary;
            val += val0;
            // diffusion term
            val -= eps*(u1_values[1]*val1+u1_values[2]*val2+u1_values[3]*val3);
	    val -= eps*(old_small_u1_qp[1]*val1+old_small_u1_qp[2]*val2+old_small_u1_qp[3]*val3)*(1-theta);
            // div-div term
            val -= tau_K * divU * val1;
            // pressure term
            val += p_values[0]*val1;
            // rhs term
            val += coeff[1]*test_bary;
            rhs_sub[dof_ii] += tau*area*val;
	    // mass term, without tau
	    rhs_sub[dof_ii] += area * old_small_u1_qp[0] * test_bary;

            // compute second rhs
            dof_ii2 = dof_ii + N_nodes3;
            //val = -u1_values[0]*u2_values[1] -  u2_values[0]*u2_values[2]  - u3_values[0]*u2_values[3];
	    val = (u1_values[0]+old_small_u1_qp[0])*(old_small_u2_qp[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(old_small_u2_qp[2])
		+ (u3_values[0]+old_small_u3_qp[0])*(old_small_u2_qp[3]);
	    val = - val;
            val *= test_bary*(1-theta);
	    val0 = (u1_values[0]+old_small_u1_qp[0])*(u2_values[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(u2_values[2]) 
		+ (u3_values[0]+old_small_u3_qp[0])*(u2_values[3]);
	    val0 = - val0;
            val0 *= test_bary;
	    val += val0;
            val -= eps *(u2_values[1]*val1+u2_values[2]*val2+u2_values[3]*val3);
	    val -= eps*(old_small_u2_qp[1]*val1+old_small_u2_qp[2]*val2+old_small_u2_qp[3]*val3)*(1-theta);
            val -= tau_K * divU * val2;
            val += p_values[0]*val2;
            val += coeff[2]*test_bary;
            rhs_sub[dof_ii2] += tau*area*val;
	    // mass term, without tau
	    rhs_sub[dof_ii2] += area * old_small_u2_qp[0] * test_bary;

            // compute third rhs
            dof_ii2 = dof_ii + 2*N_nodes3;
            //val = -u1_values[0]*u3_values[1] -  u2_values[0]*u3_values[2]  - u3_values[0]*u3_values[3];
	    val = (u1_values[0]+old_small_u1_qp[0])*(old_small_u3_qp[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(old_small_u3_qp[2])
		+ (u3_values[0]+old_small_u3_qp[0])*(old_small_u3_qp[3]);
	    val = - val;
            val *= test_bary*(1-theta);
            //val = -u1_values[0]*u3_values[1] -  u2_values[0]*u3_values[2]  - u3_values[0]*u3_values[3];
	    val0 = (u1_values[0]+old_small_u1_qp[0])*(u3_values[1])
		+ (u2_values[0]+old_small_u2_qp[0])*(u3_values[2]) 
		+ (u3_values[0]+old_small_u3_qp[0])*(u3_values[3]);
	    val0 = - val0;
            val0 *= test_bary;
	    val += val0;
            val -= eps*(u3_values[1]*val1+u3_values[2]*val2+u3_values[3]*val3);
	    val -= eps*(old_small_u3_qp[1]*val1+old_small_u3_qp[2]*val2+old_small_u3_qp[3]*val3)*(1-theta);
            val -= tau_K * divU * val3;
            val += p_values[0]*val3;
            val += coeff[3]*test_bary;
            rhs_sub[dof_ii2] += tau*area*val;
	    // mass term, without tau 
	    //if (TDatabase::TimeDB->CURRENTTIME > 1.0)
	     rhs_sub[dof_ii2] += area * old_small_u3_qp[0] * test_bary;


            // jj -- ansatz function
            for (jj=0;jj<8;jj++)
            {
              if (dof[jj] == -1)
                continue;
              // entry (jj,ii) exists
              // compute index in array containing the entry

              index = dof_ii*N_nodes3+dof[jj];
              //OutPut(" ii " << ii << " jj " << jj << " ind " << index )
              // add Laplacian
              // gradient of test fct. is (b[4*ii+1]+ b[4*ii+3]*y , b[4*ii+2]+ b[4*ii+3]*x)
              // gradient of ansatz fct. is (b[4*jj+1]+ b[4*jj+3]*y , b[4*jj+2]+ b[4*jj+3]*x)
              // approximate integral by midpoint rule

              val1_j = b[8*jj+1]+ b[8*jj+4]*ys + b[8*jj+5]*zs + b[8*jj+7]*ys*zs;
              val2_j = b[8*jj+2]+ b[8*jj+4]*xs + b[8*jj+6]*zs + b[8*jj+7]*xs*zs;
              val3_j = b[8*jj+3]+ b[8*jj+5]*xs + b[8*jj+6]*ys + b[8*jj+7]*xs*ys;

              val = val1_j * val1 + val2_j * val2 + val3_j * val3;
              val *= theta* tau*(eps+mu_T) * area;
              a_sub[index] += val;

              // add convective term
              // convection times gradient of ansatz fct.
              val = (u1_values[0]+old_small_u1_qp[0]) * val1_j 
		  + (u2_values[0]+old_small_u2_qp[0]) * val2_j 
		  + (u3_values[0]+old_small_u3_qp[0]) * val3_j;
              // multiplied with test fct.
              val *= theta*tau*test_bary*area;
              a_sub[index] += val;

              // add time-derivative terms
              val = (b[8*jj] +  b[8*jj+1] * xs + b[8*jj+2] * ys + b[8*jj+3] * zs
                + b[8*jj+4] * xs * ys +  b[8*jj+5] * xs * zs +  b[8*jj+6] * ys * zs
                + b[8*jj+7] * xs * ys *zs) * test_bary;
              val *= area;
              a_sub[index] += val;
            }                    // end jj
          }                      // end ii
        }                        // end kc
      }                          // end jc
    }                            // end ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* compute solution of the RFB problem                                   */
    /*************************************************************************/
    //SolveMultipleSystems(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 3);
    SolveLinearSystemLapack(a_sub, rhs_sub, N_nodes3, N_nodes3);
    //SolveMultipleSystemsLapack(a_sub, rhs_sub, N_nodes3, N_nodes3, N_nodes3 , 3);
   
    //for (i=0;i<3*N_nodes3;i++)
//	OutPut(i << " " << rhs_sub[i] << endl);
     /*************************************************************************/
    /* the RFB equation is solved on the mesh cell                           */
    /*************************************************************************/

    /*************************************************************************/
    /* compute information on the global d.o.f. connected to the global      */
    /* mesh cell                                                             */
    /*************************************************************************/

    // compute information about the global d.o.f.
    CurrentElement = fespace->GetFE3D(i, cell);
    // number of basis functions (= number of d.o.f.)
    N_ = N_BaseFunct[CurrentElement];
    // get FE object
    FE_Obj = TFEDatabase3D::GetFE3D(CurrentElement);
    // get ID for reference transformation
    RefTrans = FE_Obj->GetRefTransID();
    // set cell for reference transformation
    TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);

    // get base function object
    bf = FE_Obj->GetBaseFunct3D();
    test_value = new double[4*N_];
    test_value_x = test_value + N_;
    test_value_y = test_value_x + N_;
    test_value_z = test_value_y + N_;
    integral_1 = new double[3 * N_];
    memset(integral_1,0,3*N_*SizeOfDouble);
    integral_2 = integral_1 + N_;
    integral_3 = integral_2 + N_;

    // the array which gives the mapping of the local to the global d.o.f.
    globdof = global_numbers+begin_index[i];

    /*************************************************************************/
    /* loop over the mesh cells in the subgrid of the current global mesh    */
    /* cell for evaluating the integrals for the rhs of the global equation  */
    /*************************************************************************/

    for (ic=0;ic<N_sub;ic++)
    {
      for (jc=0;jc<N_sub;jc++)
      {
        for (kc=0;kc<N_sub;kc++)
        {
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc, x, y, z,
            &xc[0], &yc[0], &zc[0]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc, kc+1, x, y, z,
            &xc[1], &yc[1], &zc[1]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc+1, x, y, z,
            &xc[2], &yc[2], &zc[2]);
          VertexCoordinatesSubMeshCell(N_sub, ic, jc+1, kc, x, y, z,
            &xc[3], &yc[3], &zc[3]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc, x, y, z,
            &xc[4], &yc[4], &zc[4]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc, kc+1, x, y, z,
            &xc[5], &yc[5], &zc[5]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc+1, x, y, z,
            &xc[6], &yc[6], &zc[6]);
          VertexCoordinatesSubMeshCell(N_sub, ic+1, jc+1, kc, x, y, z,
            &xc[7], &yc[7], &zc[7]);

          // compute local dof at this mesh cell
          // compute dof of the `left lower' vertex
          dof[0] = (ic-1)*N_nodes*N_nodes+ (jc-1)*N_nodes+kc-1;
          dof[1] = dof[0] + 1;
          dof[2] = dof[1] + N_nodes;
          dof[3] = dof[0] + N_nodes;
          dof[4] = dof[0] + N_nodes*N_nodes;
          dof[5] = dof[1] + N_nodes*N_nodes;
          dof[6] = dof[2] + N_nodes*N_nodes;
          dof[7] = dof[3] + N_nodes*N_nodes;
          // correct for submesh cells with vertices on boundary
          // lower boundary
          if (ic == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[2] = -1;
            dof[3] = -1;
          }
          // upper boundary
          if (ic==N_sub-1)
          {
            dof[4] = -1;
            dof[5] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // front boundary
          if (jc == 0)
          {
            dof[0] = -1;
            dof[1] = -1;
            dof[4] = -1;
            dof[5] = -1;
          }
          // back boundary
          if (jc==N_sub-1)
          {
            dof[2] = -1;
            dof[3] = -1;
            dof[6] = -1;
            dof[7] = -1;
          }
          // left boundary
          if (kc == 0)
          {
            dof[0] = -1;
            dof[3] = -1;
            dof[4] = -1;
            dof[7] = -1;
          }
          // right boundary
          if (kc==N_sub-1)
          {
            dof[1] = -1;
            dof[2] = -1;
            dof[5] = -1;
            dof[6] = -1;
          }

          // barycenter
          xs = xc[0] + xc[1] + xc[2] + xc[3] + xc[4] + xc[5] + xc[6] + xc[7];
          xs /= 8;
          ys = yc[0] + yc[1] + yc[2] + yc[3] + yc[4] + yc[5] + yc[6] + yc[7];
          ys /= 8;
          zs = zc[0] + zc[1] + zc[2] + zc[3] + zc[4] + zc[5] + zc[6] + zc[7];
          zs /= 8;

	  jj = i * 3 * N_nodes3;
	  for (ii=0; ii < 8;ii++)
	  {
	      if (dof[ii]!=-1)
	      {
		  old_small_u1[ii] = old_small_scales[jj + dof[ii]];
		  old_small_u2[ii] = old_small_scales[jj + N_nodes3 + dof[ii]];
		  old_small_u3[ii] = old_small_scales[jj + 2*N_nodes3 + dof[ii]];
		  //OutPut(old_small_u3[ii] << " ");
		  //  if (i==379)
		  //  OutPut(" " << old_small_u3[ii]);
	      }
	      else
	      {
		  old_small_u1[ii] = old_small_u2[ii] = old_small_u3[ii] = 0;
	      }
	  }

          // compute convection field in the center of the submesh cell
          u1->FindGradientLocal(cell,i,xs,ys,zs,u1_values);
          u2->FindGradientLocal(cell,i,xs,ys,zs,u2_values);
          u3->FindGradientLocal(cell,i,xs,ys,zs,u3_values);

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            coeff_old_small_u1[jj] = old_small_u1[jj];
            coeff_old_small_u1[jj+8] = old_small_u2[jj];
            coeff_old_small_u1[jj+16] = old_small_u3[jj];	    
          }

          SolveMultipleSystemsLapack(a,coeff_old_small_u1,8,8,8,3);

	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u1,xs,ys,zs,old_small_u1_qp);
	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u2,xs,ys,zs,old_small_u2_qp);
	  Compute_Q1_Value_Gradient_RFB(coeff_old_small_u3,xs,ys,zs,old_small_u3_qp);

          // compute the value of the test functions and its derivatives
          //in the bary center of the submesh cell
          TFEDatabase3D::GetRefFromOrig(RefTrans, xs, ys, zs, xi, eta, zeta);
          bf->GetDerivatives(D000, xi, eta, zeta, test_value);
          bf->GetDerivatives(D100, xi, eta, zeta, test_value_x);
          bf->GetDerivatives(D010, xi, eta, zeta, test_value_y);
          bf->GetDerivatives(D001, xi, eta, zeta, test_value_z);

          //OutPut("test " << test_value[0] << " " << test_value[1]  << " " <<
          //       test_value[2]  << " " << test_value[3] << endl);

          // compute the derivative of the solution of the RFB problem
          // in the barycenter of the submesh cell
          // compute for this purpose the bilinear function on the mesh cell

          // set matrices for computation of the coefficients
          // of the bilinear function
          for (jj=0;jj<8;jj++)
          {
            a[8*jj] = 1;
            a[8*jj+1] = xc[jj];
            a[8*jj+2] = yc[jj];
            a[8*jj+3] = zc[jj];
            a[8*jj+4] = xc[jj]*yc[jj];
            a[8*jj+5] = xc[jj]*zc[jj];
            a[8*jj+6] = yc[jj]*zc[jj];
            a[8*jj+7] = xc[jj]*yc[jj]*zc[jj];
            if (dof[jj]!=-1)
            {
              b[jj] =  rhs_sub[dof[jj]];
              b[jj+8] = rhs_sub[dof[jj]+N_nodes3];
              b[jj+16] = rhs_sub[dof[jj]+2*N_nodes3];
            }
            else
              b[jj] = b[jj+8] = b[jj+16] = 0;
          }

          // solve system for the coefficients of the bilinear function
          //SolveMultipleSystems(a,b,8,8,8,3);
          SolveMultipleSystemsLapack(a,b,8,8,8,3);
          /*for (jj=0;jj<2;jj++)
            OutPut("#" <<jj << " " << b[4*jj] << " " <<  b[4*jj+1]<< " " <<
            b[4*jj+2] << " " <<  b[4*jj+3] << endl);*/

          // compute values of tilde_u in the barycenter
          rhs1 = b[0] + b[1]*xs + b[2]*ys + b[3]*zs
            + b[4]*xs*ys + b[5]*xs*zs + b[6]*ys*zs + b[7]*xs*ys*zs;
          rhs2 = b[8] + b[9]*xs + b[10]*ys + b[11]*zs
            + b[12]*xs*ys + b[13]*xs*zs + b[14]*ys*zs + b[15]*xs*ys*zs;
          rhs3 = b[16] + b[17]*xs + b[18]*ys + b[19]*zs
            + b[20]*xs*ys + b[21]*xs*zs + b[22]*ys*zs + b[23]*xs*ys*zs;

          // compute gradients of tilde_u in the barycenter
          rhs1_x = b[1] + b[4]*ys+ b[5]*zs + b[7]*ys*zs;
          rhs1_y = b[2] + b[4]*xs+ b[6]*zs + b[7]*xs*zs;
          rhs1_z = b[3] + b[5]*xs+ b[6]*ys + b[7]*xs*ys;
          rhs2_x = b[9] + b[12]*ys+ b[13]*zs + b[15]*ys*zs;
          rhs2_y = b[10] + b[12]*xs+ b[14]*zs + b[15]*xs*zs;
          rhs2_z = b[11] + b[13]*xs+ b[14]*ys + b[15]*xs*ys;
          rhs3_x = b[17] + b[20]*ys+ b[21]*zs + b[23]*ys*zs;
          rhs3_y = b[18] + b[20]*xs+ b[22]*zs + b[23]*xs*zs;
          rhs3_z = b[19] + b[21]*xs+ b[22]*ys + b[23]*xs*ys;

          // update the integrals
          for (j=0;j<N_;j++)
          {
	      /*
            // add convection term
            val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y + u3_values[0] * rhs1_z )
              *test_value[j]*area*theta1*tau;
            //add laplacian
            val += (rhs1_x*test_value_x[j] + rhs1_y*test_value_y[j] + rhs1_z*test_value_z[j]) *area*eps*theta1*tau;
            // add term from time derivative
            //val += rhs1*test_value[j]*area;
            integral_1[j] += val;

            val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y + u3_values[0] * rhs2_z)
              *test_value[j]*area*theta1*tau;
            val += (rhs2_x*test_value_x[j] + rhs2_y*test_value_y[j] + rhs2_z*test_value_z[j]) *area*eps*theta1*tau;
            //val += rhs2*test_value[j]*area;
            integral_2[j] += val;

            val = (u1_values[0] * rhs3_x + u2_values[0] * rhs3_y + u3_values[0] * rhs3_z)
              *test_value[j]*area*theta1*tau;
            val += (rhs3_x*test_value_x[j] + rhs3_y*test_value_y[j] + rhs3_z*test_value_z[j]) *area*eps*theta1*tau;
            //val += rhs3*test_value[j]*area;
            integral_3[j] += val;
	      */
              // add convection term
            val = (u1_values[0] * rhs1_x + u2_values[0] * rhs1_y + u3_values[0] * rhs1_z )
              *test_value[j]*area*theta1*tau;
            val = (rhs1 * u1_values[1] + rhs2 * u1_values[2] + rhs3 * u1_values[3] )
              *test_value[j]*area*theta1*tau;
	    val += (rhs1 * rhs1_x + rhs2 * rhs1_y + rhs3 * rhs1_z )
                *test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * old_small_u1_qp[1] + old_small_u2_qp[0] * old_small_u1_qp[2] 
		      + old_small_u3_qp[0] * old_small_u1_qp[3] )
		*test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * u1_values[1] + old_small_u2_qp[0] * u1_values[2]
		    + old_small_u3_qp[0] * u1_values[3] )
		    *test_value[j]*area*theta1*tau;
            val += (u1_values[0] * old_small_u1_qp[1]  + u2_values[0] * old_small_u1_qp[2] 
		   + u3_values[0] * old_small_u1_qp[3] )
		*test_value[j]*area*theta1*tau;
            //add laplacian
	    val += (rhs1_x*test_value_x[j] + rhs1_y*test_value_y[j] + rhs1_z*test_value_z[j]) 
		*area*eps*theta1*tau;
	    val += (old_small_u1_qp[1]*test_value_x[j] + old_small_u1_qp[2]*test_value_y[j] 
	    	    + old_small_u1_qp[3]*test_value_z[j]) *area*eps*theta1*tau;
            // add term from time derivative
            val += (rhs1-old_small_u1_qp[0])*test_value[j]*area;
	    integral_1[j] += val;

	    // add convection term
            val = (u1_values[0] * rhs2_x + u2_values[0] * rhs2_y + u3_values[0] * rhs2_z )
		*test_value[j]*area*theta1*tau;
	    val += (rhs1 * rhs2_x + rhs2 * rhs2_y + rhs3 * rhs2_z )
                *test_value[j]*area*theta1*tau;
            val = (rhs1 * u2_values[1] + rhs2 * u2_values[2] + rhs3 * u2_values[3] )
              *test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * old_small_u2_qp[1] + old_small_u2_qp[0] * old_small_u2_qp[2] 
		      + old_small_u3_qp[0] * old_small_u2_qp[3] )
		*test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * u2_values[1] + old_small_u2_qp[0] * u2_values[2]
		    + old_small_u3_qp[0] * u2_values[3] )
		    *test_value[j]*area*theta1*tau;
            val += (u1_values[0] * old_small_u2_qp[1]  + u2_values[0] * old_small_u2_qp[2] 
		   + u3_values[0] * old_small_u2_qp[3] )
		*test_value[j]*area*theta1*tau;
	    //add laplacian
	    val += (rhs2_x*test_value_x[j] + rhs2_y*test_value_y[j] + rhs2_z*test_value_z[j]) 
		*area*eps*theta1*tau;
	    val += (old_small_u2_qp[1]*test_value_x[j] + old_small_u2_qp[2]*test_value_y[j] 
	    	    + old_small_u2_qp[3]*test_value_z[j]) *area*eps*theta1*tau;
            // add term from time derivative
            val += (rhs2-old_small_u2_qp[0])*test_value[j]*area;
	    integral_2[j] += val;


	    // add convection term
            val = (u1_values[0] * rhs3_x + u2_values[0] * rhs3_y + u3_values[0] * rhs3_z )
		*test_value[j]*area*theta1*tau;
	    val += (rhs1 * rhs3_x + rhs2 * rhs3_y + rhs3 * rhs3_z )
                *test_value[j]*area*theta1*tau;
            val = (rhs1 * u3_values[1] + rhs2 * u3_values[2] + rhs3 * u3_values[3] )
              *test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * old_small_u3_qp[1] + old_small_u2_qp[0] * old_small_u3_qp[2] 
		      + old_small_u3_qp[0] * old_small_u3_qp[3] )
		*test_value[j]*area*theta1*tau;
	    val += (old_small_u1_qp[0] * u3_values[1] + old_small_u2_qp[0] * u3_values[2]
		    + old_small_u3_qp[0] * u3_values[3] )
		    *test_value[j]*area*theta1*tau;
            val += (u1_values[0] * old_small_u3_qp[1]  + u2_values[0] * old_small_u3_qp[2] 
		   + u3_values[0] * old_small_u3_qp[3] )
		*test_value[j]*area*theta1*tau;
	    //add laplacian
	    val += (rhs3_x*test_value_x[j] + rhs3_y*test_value_y[j] + rhs3_z*test_value_z[j]) 
		*area*eps*theta1*tau;
	    val += (old_small_u3_qp[1]*test_value_x[j] + old_small_u3_qp[2]*test_value_y[j] 
		    + old_small_u3_qp[3]*test_value_z[j]) *area*eps*theta1*tau;
            // add term from time derivative
            val += (rhs3-old_small_u3_qp[0])*test_value[j]*area;
	    integral_3[j] += val;

          }
        }                        // kc
      }                          // jc
    }                            // ic

    /*************************************************************************/
    /* end of loop over the mesh cells in the subgrid of the current         */
    /* global mesh cell                                                      */
    /*************************************************************************/

    /*************************************************************************/
    /* update the global rhs                                                 */
    /*************************************************************************/
    for (j=0;j<N_;j++)
    {
      gdof = globdof[j];
      value_global[0] += integral_1[j]*integral_1[j];
      value_global[1] += integral_2[j]*integral_2[j];
      value_global[2] += integral_3[j]*integral_3[j];
      
      //if (TDatabase::TimeDB->CURRENTTIME > 1.0)
      {
	  rhs[gdof] -= integral_1[j];
	  rhs[gdof+N_U] -= integral_2[j];
	  rhs[gdof+2*N_U] -= integral_3[j];
      }
    }
    delete test_value;
    delete integral_1;

    // save  bubble solution for next time
    //if (i==379)
//	OutPut("cell " << i);
    jj = i*3*N_nodes3;
    for (ii=0; ii < 3*N_nodes3; ii++)
    {
	old_small_scales[jj+ii] = rhs_sub[ii];
//	if (i==379)
	//    OutPut(ii << " " << rhs_sub[ii] << endl);
    }

  }

  delete a_sub;
  delete rhs_sub;
  delete coeff;

  OutPut("value_global " << sqrt(value_global[0]) << " " << sqrt(value_global[1]) 
	 << " "  << sqrt(value_global[2]) << " " << endl);
}

#endif
