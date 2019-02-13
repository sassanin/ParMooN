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
// SSMUM.C
//
// Purpose:     shear slip mesh update method
//
// Authors:     Volker John, 2008/05/22
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
#include <FEM_TVD_FCT.h>
#include <BoundComp.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <BdCircle.h>
#include <SSMUM.h>

#include <stdlib.h>
#include <string.h>
#include <MooNMD_Io.h>
#include <math.h>

/*******************************************************************************/
//
// WriteGridGnu
// writes a gnuplot output for the grid, just for debugging
//
/*******************************************************************************/

int WriteGridGnu(const char *name, TCollection *coll)
{
  int i, ii, j, jj, N_Cells, N_Edges, found, ver_on_strip, common_ver;
  int common_vert0[2], common_vert1[2], not_common_0, not_common_1;
  double x[3], y[3], r;
  const int *TmpEdVer;
  TBaseCell *cell, *cell_n;
  TVertex *ver[3], *ver_n[3];
  TRefDesc *refdesc;

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(8);

  N_Cells = coll->GetN_Cells();

  for(i=0;i<N_Cells;i++)
    //for(i=0;i<30;i++)
  {
    // next cell
    cell = coll->GetCell(i);

    N_Edges = cell->GetN_Edges();
    if (N_Edges != 3)
    {
      OutPut("Only 3 edges are allowed (triangles) !!!" << endl);
      exit(4711);
    }
    for (j=0;j<N_Edges; j++)
    {
      ver[j] = cell->GetVertex(j);
      x[j] = ver[j]->GetX();
      y[j] = ver[j]->GetY();
      dat << x[j] << " " << y[j] << " " << 0 << endl;
    }
    dat << x[0] << " " << y[0] << " " << 0 << endl;
    dat << endl;
  }                                               // end i
  dat.close();

  return 0;
}


/*******************************************************************************/
//
// checks orientation of the vertices
// input: x,y coordinates of the vertices
//
/*******************************************************************************/
double CheckOrientation(double *x, double *y)
{
  double a[2], b[2], val;

  a[0] = x[1] - x[0];
  a[1] = y[1] - y[0];
  b[0] = x[2] - x[0];
  b[1] = y[2] - y[0];
  // compute cross product
  // should be positive for counter clockwise orientation
  val = a[0]*b[1] - b[0]*a[1];
  return val;
}


/*******************************************************************************/
//
// checks if the point (x,y) is in cell
// works only for triangles
// input: cell and x,y coordinates of the vertices
// output 0 - no, 1 yes
//
/*******************************************************************************/

int PointInCell(TBaseCell *cell, double x_coord, double y_coord)
{
  int N_Cells,found,j,N_Vertices;
  double x[3], y[3], v_1x, v_1y, v_0x, v_0y, v_2x, v_2y, det, s, t;
  TVertex *ver[3];

  found=0;
  // number of vertices of cell
  N_Vertices = cell->GetN_Vertices();
  if (N_Vertices != 3)
  {
    OutPut("Only 3 vertices are allowed (triangles) !!!" << endl);
    exit(4711);
  }

  for (j=0;j<N_Vertices; j++)
  {
    ver[j] = cell->GetVertex(j);
    x[j] = ver[j]->GetX();
    y[j] = ver[j]->GetY();
  }

  v_0x=x_coord-x[0];
  v_0y=y_coord-y[0];
  v_1x=x[1]-x[0];
  v_1y=y[1]-y[0];
  v_2x=x[2]-x[0];
  v_2y=y[2]-y[0];
  // check if a point is a convex combination of the three vertices
  det = v_1x*v_2y - v_2x*v_1y;
  s = (v_2y*v_0x - v_2x*v_0y)/det;
  t = (v_1x*v_0y - v_1y*v_0x)/det;

  // check of convex combination
  if ((s>=0)&&(t>=0)&&(s+t)<=1)
    found=1;

  return found;
}


/*******************************************************************************/
//
// SwapEdges
// routine for swapping the edges in the shear slip mesh update method in 2d
//
/*******************************************************************************/

void SwapEdges(const char *name, TCollection *coll,
TFEFunction2D *u1, TFEFunction2D *u2, double *tangential_values_ssl)
{
  int i, ii, j, jj, k, ijk, N_Cells, N_Edges, found, common_ver;
  int common_vert0[2], common_vert1[2], not_common_0, not_common_1;
  int index0, index1, N_Edges_n, N_Edges_nn, me_found;
  int ver_below_strip, ver_above_strip, count;
  const int max_joint_changes = 4096;
  int joint_changes[max_joint_changes][4], N_joint_changes=0, N_joint_changes_1=0;
  double x[3], y[3], r, x_n[3], y_n[3], r_n, x_nn[2], y_nn[2], eps=1e-8, av_rad;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double x_m, y_m, v_1x, v_1y, v_2x,v_2y,s,t,det,v_0x,v_0y;
  double L, tx, ty, val[3], I, I_1, I_2, I_3, I_4[2], r1;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  double inner_radius = TDatabase::ParamDB->SSMUM_INNER_RADIUS;
  double outer_radius = TDatabase::ParamDB->SSMUM_OUTER_RADIUS;
  const int *TmpEdVer, *TmpEdVer_n, *TmpEdVer_nn;
  TBaseCell *cell, *cell_n, *cell_nn;
  TVertex *ver[3], *ver_n[3];
  TRefDesc *refdesc, *refdesc_n, *refdesc_nn;
  TJoint *joint, *joint_n, *joint_nn, *joints[max_joint_changes];

  N_Cells = coll->GetN_Cells();
  // average radius of strip
  av_rad = (inner_radius + outer_radius)/2.0;

  /*  for(i=0;i<N_Cells;i++)
    {
        // next cell
        cell = coll->GetCell(i);

       N_Edges = cell->GetN_Edges();
       if (N_Edges != 3)
       {
     OutPut("Only 3 edges are allowed (triangles) !!!" << endl);
     exit(4711);
       }
  for (j=0;j<N_Edges; j++)
  {
  ver[j] = cell->GetVertex(j);
  x[j] = ver[j]->GetX();
  y[j] = ver[j]->GetY();
  }
  r = CheckOrientation(x,y);
  if (r<0)
  {
  OutPut("Wrong orientation of vertices " << r << endl);
  exit(4711);
  }
  }
  */

  // initialize clipboars
  for(i=0;i<N_Cells;i++)
  {
    cell=coll->GetCell(i);
    cell->SetClipBoard(-i-1);
  }

  // find cells on the shear slip layer
  found = 0;
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);

    N_Edges = cell->GetN_Edges();
    if (N_Edges != 3)
    {
      OutPut("Only 3 edges are allowed (triangles) !!!" << endl);
      exit(4711);
    }
    for (j=0;j<N_Edges; j++)
    {
      ver[j] = cell->GetVertex(j);
      x[j] = ver[j]->GetX();
      y[j] = ver[j]->GetY();
    }
    ver_below_strip = 0;
    ver_above_strip = 0;
    for (j=0;j<N_Edges; j++)
    {
      r = (x[j]-mp_x)*(x[j]-mp_x) + (y[j]-mp_y)*(y[j]-mp_y);
      r = sqrt(r);
      if (r > av_rad)
        ver_above_strip++;
      if (r <= av_rad)
        ver_below_strip++;
      // if ((fabs(r-inner_radius) < 1e-6)||(fabs(r-outer_radius) < 1e-6))
      //    ver_on_strip++;
    }
    if ((ver_below_strip)&&(ver_above_strip))
    {
      if (TDatabase::ParamDB->SC_VERBOSE>1)
      {
        OutPut(found << " " << i << " cell: (" << x[0] <<","<<y[0]<<"), ("
          << x[1] <<","<<y[1]<< "), ("
          << x[2] <<","<<y[2]<< ") " << av_rad << endl);
      }
      cell->SetClipBoard(-cell->GetClipBoard());
      found++;
    }
  }                                               // end i
  //OutPut("found " << found << " inner " << inner_radius << " outer "  << outer_radius << endl);

  count = 0;
  // find the pairs of mesh cells which have a common edge that
  // should be swapped
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    if (cell->GetClipBoard() < 0)
      continue;
    refdesc=cell->GetRefDesc();                   // get refinement descriptor
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
    N_Edges = cell->GetN_Edges();
    if (TDatabase::ParamDB->SC_VERBOSE>1)
      OutPut("cell " << i);
    for (j=0;j<N_Edges; j++)
    {
      ver[j] = cell->GetVertex(j);
      x[j] = ver[j]->GetX();
      y[j] = ver[j]->GetY();
      if (TDatabase::ParamDB->SC_VERBOSE>1)
        OutPut(" ; " << x[j] << " " << y[j]);
                                                  // neighbour cell
      cell_n=cell->GetJoint(j)->GetNeighbour(cell);
      if (TDatabase::ParamDB->SC_VERBOSE>1)
        OutPut(" neigh " << cell_n->GetClipBoard());
    }
    if (TDatabase::ParamDB->SC_VERBOSE>1)
    {
      OutPut(endl);
      for (j=0;j<N_Edges; j++)
      {
        OutPut(" ; " << x[j] << " " << y[j]);
      }
    }
    // check the other mesh cells
    for(ii=0;ii<N_Cells;ii++)
    {
      if (ii==i)
        continue;
      // next cell
      cell_n = coll->GetCell(ii);
      if (cell_n->GetClipBoard() < 0 )
        continue;
      for (j=0;j<N_Edges; j++)
      {
        ver_n[j] = cell_n->GetVertex(j);
      }
      // check if there are two vertices in common with cell i
      common_ver = 0;
      for (j=0;j<N_Edges; j++)
      {
        for (jj=0;jj<N_Edges; jj++)
        {
          if (ver[j] == ver_n[jj])
          {
            common_vert0[common_ver] = j;
            common_vert1[common_ver] = jj;
            common_ver++;
          }
        }
      }
      // found mesh cell with two common vertices
      if (common_ver == 2)
      {
        if (TDatabase::ParamDB->SC_VERBOSE>1)
          OutPut(" common ");
        // compute coordinates of vertices
        for (j=0;j<N_Edges; j++)
        {
          x_n[j] = ver_n[j]->GetX();
          y_n[j] = ver_n[j]->GetY();
          if (TDatabase::ParamDB->SC_VERBOSE>1)
            OutPut(" ; " << x_n[j] << " " << y_n[j]);
        }
        // find the vertices which are not common
        for (j=0;j<N_Edges; j++)
        {
          if ((common_vert0[0]!=j) &&(common_vert0[1]!=j))
          {
            not_common_0 = j;
            break;
          }
        }
        for (j=0;j<N_Edges; j++)
        {
          if ((common_vert1[0]!=j) &&(common_vert1[1]!=j))
          {
            not_common_1 = j;
            break;
          }
        }
        if (TDatabase::ParamDB->SC_VERBOSE>1)
          OutPut(" " << not_common_0 << " " << not_common_1);
        // check if this is the correct pair
        r = (x[not_common_0] - mp_x)*(x[not_common_0] - mp_x)+
          (y[not_common_0] - mp_y)*(y[not_common_0] - mp_y);
        r = sqrt(r);
        // next vertex
        r_n = (x[(not_common_0+1)%3] - mp_x)*(x[(not_common_0+1)%3] - mp_x)+
          (y[(not_common_0+1)%3] - mp_y)*(y[(not_common_0+1)%3] - mp_y);
        r_n = sqrt(r_n);
        if (TDatabase::ParamDB->SC_VERBOSE>1)
          OutPut(" r " << r << " " << r_n << endl);
        // wrong pair
        //if (fabs(r-r_n)>1e-6)
        if (((r<av_rad) && (r_n>=av_rad)) || ((r>=av_rad) && (r_n<av_rad)))
        {
          if (TDatabase::ParamDB->SC_VERBOSE>1)
            OutPut("wrong "  << r << " " << r_n << endl);
          continue;
        }
        if (TDatabase::ParamDB->SC_VERBOSE>1)
        {
          if (TDatabase::ParamDB->SC_VERBOSE>1)
            OutPut("Pair " <<  cell->GetClipBoard()
              << " " << cell_n->GetClipBoard() <<
              " not common " << not_common_0 <<
              " " << not_common_1 << " " << r << " " << r_n << endl);
        }
        // compute line integral for interpolation
        // vertex which belongs only to this mesh cell: not_common_0
        //       coordinates: x[not_common_0], y[not_common_0]
        //       coordinates of other vertices
        //                    x[(not_common_0+1)%3], y[(not_common_0+1)%3]
        //                    x[(not_common_0+2)%3], y[(not_common_0+2)%3]

        // vertex which belongs only to the neighbor mesh cell: not_common_1
        //       coordinates: x_n[not_common_1], y_n[not_common_1]
        // mid point of the non-common vertices
        x_m = (x[not_common_0]+x_n[not_common_1])/2;
        y_m = (y[not_common_0]+y_n[not_common_1])/2;
        switch(TDatabase::ParamDB->SSMUM_INTERPOLATION)
        {
          case 0:
            v_0x = x_m - x[not_common_0];
            v_0y = y_m - y[not_common_0];
            L = sqrt(v_0x*v_0x + v_0y*v_0y);
            // tangential from vertex[not_common_0] to vertex_n[not_common_1]
            tx = v_0x/L;
            ty = v_0y/L;
            // orientation always from interior to outside
            r1 = sqrt((x_n[not_common_1]-mp_x)*(x_n[not_common_1]-mp_x)
              +(y_n[not_common_1]-mp_y)*(y_n[not_common_1]-mp_y));
            if (r>r1)
            {
              tx = -tx;
              ty = -ty;
            }

            // check of convex combination
            if (PointInCell(cell,x_m,y_m))
            {
              // compute velocities and apply Simpson rule
              u1->FindValueLocal(cell,i,x[not_common_0],y[not_common_0],val);
              I_1 = val[0] * tx;
              u2->FindValueLocal(cell,i,x[not_common_0],y[not_common_0],val);
              I_1 += val[0] * ty;
              u1->FindValueLocal(cell,i,x_m,y_m,val);
              I_2 = val[0] * tx;
              //OutPut(val[0] << " " );
              u2->FindValueLocal(cell,i,x_m,y_m,val);
              I_2 += val[0] * ty;
              //OutPut(val[0] << " " << endl);
              u1->FindValueLocal(cell_n,ii,x_n[not_common_1],y_n[not_common_1],val);
              I_3 = val[0] * tx;
              u2->FindValueLocal(cell_n,ii,x_n[not_common_1],y_n[not_common_1],val);
              I_3 += val[0] * ty;
              // factor 2 because L is only from one point to the midpoint
              I= 2*L*(I_1+4*I_2+I_3)/6;
            }
            else
            {
              // compute velocities and apply Simpson rule
              u1->FindValueLocal(cell,i,x[not_common_0],y[not_common_0],val);
              I_1 = val[0] * tx;
              u2->FindValueLocal(cell,i,x[not_common_0],y[not_common_0],val);
              I_1 += val[0] * ty;
              u1->FindValueLocal(cell_n,ii,x_m,y_m,val);
              I_2 = val[0] * tx;
              //OutPut(val[0] << " ");
              u2->FindValueLocal(cell_n,ii,x_m,y_m,val);
              I_2 += val[0] * ty;
              //OutPut(val[0] << " " << endl);
              u1->FindValueLocal(cell_n,ii,x_n[not_common_1],y_n[not_common_1],val);
              I_3 = val[0] * tx;
              u2->FindValueLocal(cell_n,ii,x_n[not_common_1],y_n[not_common_1],val);
              I_3 += val[0] * ty;
              // factor 2 because L is only from one point to the midpoint
              I=2*L*(I_1+4*I_2+I_3)/6;
            }
            OutPut("tangent integral : A " << x[not_common_0] << " " << y[not_common_0]
              << " B " << x_n[not_common_1] << " " << y_n[not_common_1]
              << " M " << x_m << " " << y_m << " " << I << endl);
            //	      OutPut("   " << tx << " " << ty << " " << 2*L << " " << r*r << " "  << r1*r1 << endl);
            tangential_values_ssl[count] = x[not_common_0];
            tangential_values_ssl[count+1] = y[not_common_0];
            tangential_values_ssl[count+2] = I;
            tangential_values_ssl[count+3] = x_n[not_common_1];
            tangential_values_ssl[count+4] = y_n[not_common_1];
            tangential_values_ssl[count+5] = I;
            count += 6;
            break;
            // get the values in the points
          case 1:
            if (PointInCell(cell,x_m,y_m))
            {
              u1->FindValueLocal(cell,i,x_m,y_m,val);
              I_4[0] = val[0];
              u2->FindValueLocal(cell,i,x_m,y_m,val);
              I_4[1] = val[0];
            }
            else
            {
              u1->FindValueLocal(cell_n,ii,x_m,y_m,val);
              I_4[0] = val[0];
              u2->FindValueLocal(cell_n,ii,x_m,y_m,val);
              I_4[1] = val[0];
            }
            OutPut("mittelpunkt : M "<< x_m << " " << y_m << " u1 und u2 in mittelpunkt" << I_4[0] << I_4[1] << endl);

            tangential_values_ssl[count] = x_m;
            tangential_values_ssl[count+1] = y_m;
            tangential_values_ssl[count+2] = I_4[0];
            tangential_values_ssl[count+3] = I_4[1];
            count += 4;
            break;
          case 2:
            // do nothing, compute later linear interpolation
            break;
          default:
            OutPut("SSMUM_INTERPOLATION " <<
              TDatabase::ParamDB->SSMUM_INTERPOLATION <<
              " not implemented !!!" << endl);
            exit(4711);
        }

        // swap vertices
        // compute radius of vertex which is not common on
        // mesh cell i
        r = (x[not_common_0] - mp_x)*(x[not_common_0] - mp_x)+
          (y[not_common_0] - mp_y)*(y[not_common_0] - mp_y);
        r = sqrt(r);
        // find the other vertex with the same radius
        for (j=0;j<N_Edges;j++)
        {
          if  (j==not_common_0)
            continue;
          r_n = (x[j] - mp_x)*(x[j] - mp_x)+
            (y[j] - mp_y)*(y[j] - mp_y);
          r_n = sqrt(r_n);
          // vertex to swap found
          //if (fabs(r-r_n) < 1e-6)
          if (((r>=av_rad)&&(r_n>=av_rad)) || ((r<av_rad)&&(r_n<av_rad)))
          {
            index0 = j;
            if (TDatabase::ParamDB->SC_VERBOSE>1)
              OutPut(r << " "  << not_common_0 << " "  << r_n << " "<< j << endl);
            break;
          }
        }

        // do the same for mesh cell ii
        r_n = (x_n[not_common_1] - mp_x)*(x_n[not_common_1] - mp_x)+
          (y_n[not_common_1] - mp_y)*(y_n[not_common_1] - mp_y);
        r_n = sqrt(r_n);
        for (j=0;j<N_Edges;j++)
        {
          if  (j==not_common_1)
            continue;
          r =   (x_n[j] - mp_x)*(x_n[j] - mp_x)+
            (y_n[j] - mp_y)*(y_n[j] - mp_y);
          r = sqrt(r);
          // vertex to swap found
          //if (fabs(r-r_n) < 1e-6)
          if (((r>=av_rad)&&(r_n>=av_rad)) || ((r<av_rad)&&(r_n<av_rad)))
          {
            index1 = j;
            //cell->SetVertex(j, ver[not_common_1]);
            break;
          }
        }
        if (TDatabase::ParamDB->SC_VERBOSE>1)
        {
          OutPut("index " << index0 << " " << index1 << endl);
        }
        // set new coordinates to vertices
        cell->SetVertex(index0, ver_n[not_common_1]);
        cell_n->SetVertex(index1, ver[not_common_0]);
        if (TDatabase::ParamDB->SC_VERBOSE>1)
        {
          for (j=0;j<N_Edges; j++)
          {
            OutPut(" " << cell->GetVertex(j)->GetX() << " "
              << cell->GetVertex(j)->GetY() << " ") ;
          }
          OutPut(endl);

          for (j=0;j<N_Edges; j++)
          {
            OutPut(" " << cell_n->GetVertex(j)->GetX() << " "
              << cell_n->GetVertex(j)->GetY() << " ") ;
          }
          OutPut(endl);
          OutPut("joints ");
          for (j=0;j<N_Edges; j++)
          {
            OutPut(" " << cell->GetVertex(TmpEdVer[2*j])->GetX() << " "
              << cell->GetVertex(TmpEdVer[2*j])->GetY() << " ") ;
          }
          OutPut(endl);
        }
        // mark cells as done
        cell->SetClipBoard(-cell->GetClipBoard());
        cell_n->SetClipBoard(-cell_n->GetClipBoard());
        break;
      }
    }
  }
  OutPut("count " << count << endl);

  // checking orientation of the vertices
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);

    N_Edges = cell->GetN_Edges();

    for (j=0;j<N_Edges; j++)
    {
      ver[j] = cell->GetVertex(j);
      x[j] = ver[j]->GetX();
      y[j] = ver[j]->GetY();
    }
    r = CheckOrientation(x,y);
    if (r<0)
    {
      OutPut("Wrong orientation of vertices after swapping " << r << endl);
      exit(4711);
    }
  }

  for (ijk=0;ijk<1;ijk++)
  {
    // check and reset neighbour connections
    for(i=0;i<N_Cells;i++)
    {
      // next cell
      cell = coll->GetCell(i);
      refdesc=cell->GetRefDesc();                 // get refinement descriptor
      refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
      N_Edges = cell->GetN_Edges();

      // loop over all edges
      for (j=0;j<N_Edges; j++)
      {
        found = 0;
        joint = cell->GetJoint(j);
        // coordinates of face j
        x[0] = cell->GetVertex(TmpEdVer[2*j])->GetX();
        y[0] = cell->GetVertex(TmpEdVer[2*j])->GetY();
        x[1] = cell->GetVertex(TmpEdVer[2*j+1])->GetX();
        y[1] = cell->GetVertex(TmpEdVer[2*j+1])->GetY();
        // neighbour cell
        cell_n= joint->GetNeighbour(cell);
        if (cell_n == NULL)
          continue;
        if (TDatabase::ParamDB->SC_VERBOSE>1)
        {
          OutPut(i <<  " my neigh " << j << " is " <<  -cell_n->GetClipBoard()-1 <<
            " " << -joint->GetNeighbour(0)->GetClipBoard()-1 <<
            " " << -joint->GetNeighbour(1)->GetClipBoard()-1 << endl);
        }
        // get refinement descriptor of neighbour
        refdesc_n=cell_n->GetRefDesc();
        refdesc_n->GetShapeDesc()->GetEdgeVertex(TmpEdVer_n);
        N_Edges_n = cell_n->GetN_Edges();
        // loop over edges of neighbour
        me_found = 0;
        for (k=0;k<N_Edges_n; k++)
        {
          joint_n = cell_n->GetJoint(k);
          // neighbour cell of neighbour
          cell_nn = joint_n->GetNeighbour(cell_n);
          if (cell_nn == NULL)
            continue;
          if (TDatabase::ParamDB->SC_VERBOSE>1)
            OutPut(" " << k << " : " <<  -cell_nn->GetClipBoard()-1 << " ");
          // this is the original joint
          if (cell_nn == cell)
            // check if the joint is correct
          {
            me_found++;
            x_n[0] = cell_n->GetVertex(TmpEdVer_n[2*k])->GetX();
            y_n[0] = cell_n->GetVertex(TmpEdVer_n[2*k])->GetY();
            x_n[1] = cell_n->GetVertex(TmpEdVer_n[2*k+1])->GetX();
            y_n[1] = cell_n->GetVertex(TmpEdVer_n[2*k+1])->GetY();

            if (TDatabase::ParamDB->SC_VERBOSE>1)
            {
              OutPut(i << " " << j << " " << k << " neigh found " << x[0] << " " << x_n[0]  << " "
                << x[1] << " " << x_n[1] << " "  << y[0] << " " << y_n[0]  << " "
                << y[1] << " " << y_n[1]  << endl);
            }
            // check coordinates of vertices belonging to that face
            if ((fabs(x[0] - x_n[0])<eps)&&(fabs(y[0] - y_n[0])<eps)
              && (fabs(x[1] - x_n[1])<eps)&&(fabs(y[1] - y_n[1])<eps))
              found = 1;
            if ((fabs(x[0] - x_n[1])<eps)&&(fabs(y[0] - y_n[1])<eps)
              && (fabs(x[1] - x_n[0])<eps)&&(fabs(y[1] - y_n[0])<eps))
              found = 1;
            if (found)
              break;                              // loop over k
          }
          // this is the original joint or the reverse way to the
          // original mesh cell is not found
          if (((cell_nn == cell)&&(!found))||((k==N_Edges_n-1)&&(! me_found)))
          {
            if (TDatabase::ParamDB->SC_VERBOSE>1)
              OutPut("cell " << i << " does not posses correct neighbour at joint " << j <<
                " neigh " << -cell_n->GetClipBoard()-1 << " not at joint " << k <<
                " " << me_found << endl);
            // something wrong
            if (ijk)
              exit(4711);
            // loop over all mesh cells to find the correct neighbor
            for (ii=0;ii<N_Cells;ii++)
            {
              if (ii==i)
                continue;
              // get mesh cell
              cell_nn = coll->GetCell(ii);
              // get refinement descriptor
              refdesc_nn=cell->GetRefDesc();
              refdesc_nn->GetShapeDesc()->GetEdgeVertex(TmpEdVer_nn);
              N_Edges_nn = cell_nn->GetN_Edges();
              for (jj=0;jj<N_Edges_nn; jj++)
              {
                found = 0;
                joint_nn = cell_nn->GetJoint(jj);
                // coordinates of face jj
                x_nn[0] = cell_nn->GetVertex(TmpEdVer_nn[2*jj])->GetX();
                y_nn[0] = cell_nn->GetVertex(TmpEdVer_nn[2*jj])->GetY();
                x_nn[1] = cell_nn->GetVertex(TmpEdVer_nn[2*jj+1])->GetX();
                y_nn[1] = cell_nn->GetVertex(TmpEdVer_nn[2*jj+1])->GetY();
                if ((fabs(x[0] - x_nn[0])<eps)&&(fabs(y[0] - y_nn[0])<eps)
                  && (fabs(x[1] - x_nn[1])<eps)&&(fabs(y[1] - y_nn[1])<eps))
                  found = 1;
                if ((fabs(x[0] - x_nn[1])<eps)&&(fabs(y[0] - y_nn[1])<eps)
                  && (fabs(x[1] - x_nn[0])<eps)&&(fabs(y[1] - y_nn[0])<eps))
                  found = 1;
                if (!found)
                  continue;
                // correct face found
                // set the correct joint
                /*cell_nn->SetJoint(jj,joint);
                if (!(cell_nn->GetJoint(jj)==joint))
                    exit(1);
                joint->SetNeighbour(0,cell);
                joint->SetNeighbour(1,cell_nn);
                // replace the wrong joint at the neighbour cell
                cell_n->SetJoint(k,joint_nn);*/
                joint_changes[N_joint_changes][0] = -cell->GetClipBoard()-1;
                joint_changes[N_joint_changes][1] = j;
                joint_changes[N_joint_changes][2] = -cell_nn->GetClipBoard()-1;
                joint_changes[N_joint_changes][3] = jj;
                N_joint_changes++;
                if (N_joint_changes>=max_joint_changes)
                {
                  OutPut("max_joint_changes too small !!!" <<endl);
                  exit(4711);
                }
                if (TDatabase::ParamDB->SC_VERBOSE>1)
                  OutPut(i << " " << j << "  correct neighbour is cell " <<
                    -cell_nn->GetClipBoard()-1 << " joint " << jj << " " <<
                    -cell->GetJoint(j)->GetNeighbour(cell)->GetClipBoard()-1 << endl);

                break;                            // jj
              }
              // neighbour relations corrected, stop loop over mesh cells
              if (found)
                break;                            // ii
            }
          }
        }
        /*if (me_found!=1)
        {
            OutPut("cell " << i << " and neigh " << -cell_n->GetClipBoard()-1
             << " not mutually connected" <<endl);
            exit(4711);
            }*/

        if (TDatabase::ParamDB->SC_VERBOSE>1)
        {
          OutPut(i <<  " my neigh after " << j << " is " <<
            " " << -joint->GetNeighbour(0)->GetClipBoard()-1 <<
            " " << -joint->GetNeighbour(1)->GetClipBoard()-1 << endl);
        }
      }
    }
  }
  //WriteGridGnu(name,coll);
  if (TDatabase::ParamDB->SC_VERBOSE>1)
    OutPut(endl);

  // condense the information to the changes of the joints
  // mark multiple entries
  N_joint_changes_1 = 0;
  for (i=0;i<N_joint_changes;i++)
  {
    if (joint_changes[i][0] == -1)
      continue;
    N_joint_changes_1++;
    for(j=i+1;j<N_joint_changes;j++)
    {
      if ((joint_changes[i][0] == joint_changes[j][2])
        && (joint_changes[i][2] == joint_changes[j][0])
        && (joint_changes[i][1] == joint_changes[j][3])
        && (joint_changes[i][3] == joint_changes[j][1]))
      {
        joint_changes[j][0] = -1;
        break;
      }
    }
  }

  if (TDatabase::ParamDB->SC_VERBOSE>1)
  {
    for (i=0;i<N_joint_changes;i++)
    {
      for(j=0;j<4;j++)
      {
        OutPut(joint_changes[i][j]<< " ");
      }
      OutPut(endl);
    }
  }

  OutPut("joint changes " << N_joint_changes << " " << N_joint_changes_1 << endl);

  // remove multiple entries, concentrate the other entries at the beginning
  // of the array
  for (i=0;i<N_joint_changes_1;i++)
  {
    if (joint_changes[i][0] != -1)
      continue;
    for (j=i+1;j<N_joint_changes;j++)
    {
      if (joint_changes[j][0] != -1)
      {
        joint_changes[i][0] = joint_changes[j][0];
        joint_changes[i][1] = joint_changes[j][1];
        joint_changes[i][2] = joint_changes[j][2];
        joint_changes[i][3] = joint_changes[j][3];
        joint_changes[j][0] = -1;
        break;
      }
    }
  }

  if (TDatabase::ParamDB->SC_VERBOSE>1)
  {
    for (i=0;i<N_joint_changes_1;i++)
    {
      for(j=0;j<4;j++)
      {
        OutPut(joint_changes[i][j]<< " ");
      }
      OutPut(endl);
    }
  }

  ii = 0;
  // loop over all joints which have to be changed
  // save the joints to change in the array joints[]
  for (i=0;i<N_joint_changes_1;i++)
  {
    // get the joint from the joint_changes array
    joint = coll->GetCell(joint_changes[i][0])->GetJoint(joint_changes[i][1]);
    found = 0;
    for (j=0;j<ii;j++)
    {
      if (joints[j] == joint)
      {
        found++;
        break;
      }
    }
    if (!found)
    {
      joints[ii] = joint;
      ii++;
      if (TDatabase::ParamDB->SC_VERBOSE>1)
        OutPut(ii);
    }
    joint = coll->GetCell(joint_changes[i][2])->GetJoint(joint_changes[i][3]);
    found = 0;
    for (j=0;j<ii;j++)
    {
      if (joints[j] == joint)
      {
        found++;
        break;
      }
    }
    if (found)
      continue;
    joints[ii] = joint;
    ii++;
    if (TDatabase::ParamDB->SC_VERBOSE>1)
      OutPut(ii);
  }

  // set new joints
  for (i=0;i<N_joint_changes_1;i++)
  {
    joint = joints[i];
    joint->SetNeighbour(0,coll->GetCell(joint_changes[i][0]));
    joint->SetNeighbour(1,coll->GetCell(joint_changes[i][2]));
    coll->GetCell(joint_changes[i][0])->SetJoint(joint_changes[i][1],joint);
    coll->GetCell(joint_changes[i][2])->SetJoint(joint_changes[i][3],joint);
  }

  // JUST FOR CHECKING
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);
    refdesc=cell->GetRefDesc();                   // get refinement descriptor
    refdesc->GetShapeDesc()->GetEdgeVertex(TmpEdVer);
    N_Edges = cell->GetN_Edges();

    // loop over all edges
    for (j=0;j<N_Edges; j++)
    {
      found = 0;
      joint = cell->GetJoint(j);
      // coordinates of face j
      x[0] = cell->GetVertex(TmpEdVer[2*j])->GetX();
      y[0] = cell->GetVertex(TmpEdVer[2*j])->GetY();
      x[1] = cell->GetVertex(TmpEdVer[2*j+1])->GetX();
      y[1] = cell->GetVertex(TmpEdVer[2*j+1])->GetY();
      // neighbour cell
      cell_n= joint->GetNeighbour(cell);
      if (cell_n == NULL)
        continue;
      if (TDatabase::ParamDB->SC_VERBOSE>1)
      {
        OutPut(i <<  " my neigh " << j << " is " <<  -cell_n->GetClipBoard()-1 <<
          " " << -joint->GetNeighbour(0)->GetClipBoard()-1 <<
          " " << -joint->GetNeighbour(1)->GetClipBoard()-1 << endl);
      }
      // get refinement descriptor of neighbour
      refdesc_n=cell_n->GetRefDesc();
      refdesc_n->GetShapeDesc()->GetEdgeVertex(TmpEdVer_n);
      N_Edges_n = cell_n->GetN_Edges();
      // loop over edges of neighbour
      me_found = 0;
      for (k=0;k<N_Edges_n; k++)
      {
        joint_n = cell_n->GetJoint(k);
        // neighbour cell of neighbour
        cell_nn = joint_n->GetNeighbour(cell_n);
        if (cell_nn == NULL)
          continue;
        // this is the original joint
        if (cell_nn == cell)
          // check if the joint is correct
        {
          me_found++;
          x_n[0] = cell_n->GetVertex(TmpEdVer_n[2*k])->GetX();
          y_n[0] = cell_n->GetVertex(TmpEdVer_n[2*k])->GetY();
          x_n[1] = cell_n->GetVertex(TmpEdVer_n[2*k+1])->GetX();
          y_n[1] = cell_n->GetVertex(TmpEdVer_n[2*k+1])->GetY();

          // check coordinates of vertices belonging to that face
          if ((fabs(x[0] - x_n[0])<eps)&&(fabs(y[0] - y_n[0])<eps)
            && (fabs(x[1] - x_n[1])<eps)&&(fabs(y[1] - y_n[1])<eps))
            found = 1;
          if ((fabs(x[0] - x_n[1])<eps)&&(fabs(y[0] - y_n[1])<eps)
            && (fabs(x[1] - x_n[0])<eps)&&(fabs(y[1] - y_n[0])<eps))
            found = 1;
          if (found)
            break;                                // loop over k
        }
        // this is the original joint or the reverse way to the
        // original mesh cell is not found
        if (((cell_nn == cell)&&(!found))||((k==N_Edges_n-1)&&(! me_found)))
        {
          OutPut("cell " << i << " does not posses correct neighbour at joint " << j <<
            " neigh " << -cell_n->GetClipBoard()-1 << " not at joint " << k <<
            " " << me_found << endl);
          exit(4711);
        }
      }
    }
  }
  // checking orientation of the vertices
  for(i=0;i<N_Cells;i++)
  {
    // next cell
    cell = coll->GetCell(i);

    N_Edges = cell->GetN_Edges();

    for (j=0;j<N_Edges; j++)
    {
      ver[j] = cell->GetVertex(j);
      x[j] = ver[j]->GetX();
      y[j] = ver[j]->GetY();
    }
    r = CheckOrientation(x,y);
    if (r<0)
    {
      OutPut("Wrong orientation of vertices at the end of swapping " << r << endl);
      exit(4711);
    }
  }

}


/*******************************************************************************/
//
// RotateGrid
// controls the grid rotation
//
/*******************************************************************************/

int RotateGrid(const char *name, TCollection *coll, double swap_rotation, double &angle,
double *uoldx, double *uoldy, double *uold1, double *uold2,
TFEFunction2D *u1, TFEFunction2D *u2, double *tangential_values_ssl)
{
  int i, j, jj, k, kk, N_Cells, N_V, swapped=0, count, found;
  int *GlobalNumbers_velo, *BeginIndex_velo, *DOF, index, N_U;
  double x, y, phi, r, phi_plus, t0, t1;
  double val[1], *sol, xp[3], yp[3];
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  TBaseCell *cell, *cell1, *cell_n, *cell_nk;
  TVertex *vertex, *Vertices[4];
  TJoint *joint;
  TBoundComp *BoundComp;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  TFESpace2D *velo_space;
  OutPut("rotate" << endl);

  phi_plus = 2*Pi * TDatabase::TimeDB->TIMESTEPLENGTH * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;

  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  // set clip board of the vertices and cells
  for (j=0;j<N_Cells;j++)
  {
    cell = coll->GetCell(j);
    N_V = cell->GetN_Vertices();
    for (k=0;k<N_V;k++)
    {
      vertex =  cell->GetVertex(k);
      vertex->SetClipBoard(-1);
    }
  }

  // rotate inner mesh cells
  for (j=0;j<N_Cells;j++)
  {
    cell = coll->GetCell(j);
    N_V = cell->GetN_Vertices();
    for (k=0;k<N_V;k++)
    {
      vertex =  cell->GetVertex(k);
      vertex->GetCoords(x, y);
      r = (x - mp_x)*(x - mp_x) + (y - mp_y)*(y - mp_y);
      r = sqrt(r);
      // vertex in outer domain
      if (r>(TDatabase::ParamDB->SSMUM_OUTER_RADIUS+
        TDatabase::ParamDB->SSMUM_INNER_RADIUS)/2.0)
        continue;
      // vertex already moved
      if (vertex->GetClipBoard()==1)
        continue;
      // get argument of the vertex
      vertex->SetClipBoard(1);
      phi = atan2(y-mp_y,x-mp_x);
      //OutPut(x << " " << y << " " << phi << endl);
      phi -= phi_plus;
      x = r * cos(phi) + mp_x;
      y = r * sin(phi) + mp_y;
      vertex->SetCoords(x, y);
    }
  }                                               // end j

  /*    // update parameters on the boundary if necessary
      for (j=0;j<N_Cells;j++)
      {
    cell = coll->GetCell(j);
    N_V = cell->GetN_Vertices();
    for (k=0;k<N_V;k++)
    {
        Vertices[k] = cell->GetVertex(k);
    }
    for (k=0;k<N_V;k++)
    {
  joint = cell->GetJoint(k);
  if(joint->GetType() == BoundaryEdge ||
  joint->GetType() == IsoBoundEdge)
  {
  if(joint->GetType() == BoundaryEdge)
  {
  ((TBoundEdge *)joint)->UpdateParameters(Vertices[k],
  Vertices[(k+1)%N_V]);
  }
  else
  {
  ((TInterfaceJoint *)joint)->UpdateParameters(Vertices[k],
  Vertices[(k+1)%N_V]);
  }

  }
  }
  }
  */
  // test for swapping
  if ((int)(angle/swap_rotation) < (int) ((angle+phi_plus)/swap_rotation))
  {
    OutPut("swapping edges at step " << i << " time " << TDatabase::TimeDB->CURRENTTIME
      << " " << (int)(angle/swap_rotation)
      << " " << (int) ((angle+phi_plus)/swap_rotation) << endl);
    // save the solution on the old mesh in the nodal values of the P2_bubble
    // finite element
    // velocity space
    OutPut("a"<<endl);
    velo_space = u1->GetFESpace2D();
    // number of velocity unknowns
    N_U = u1->GetLength();
    // vector for velocity
    sol = u1->GetValues();
    // information to the degrees of freedom for the velocity
    GlobalNumbers_velo = velo_space->GetGlobalNumbers();
    BeginIndex_velo = velo_space->GetBeginIndex();
    OutPut("b"<<endl);

    count = 0;
    // cells of the old grid
    for (j=0;j<N_Cells;j++)
    {
      cell = coll->GetCell(j);
      // dof belonging to this mesh cell
      DOF =  GlobalNumbers_velo + BeginIndex_velo[j];
      N_V = cell->GetN_Vertices();
      // find coordinates corresponding to the degrees of freedom
      // vertices
      for (k=0;k<N_V;k++)
      {
        vertex =  cell->GetVertex(k);
        vertex->GetCoords(x, y);
        switch(k)
        {
          case 0:
            uoldx[count] = x;
            uoldy[count] = y;
            u1->FindValueLocal(cell,j,x,y,val);
            uold1[count] = val[0];
            u2->FindValueLocal(cell,j,x,y,val);
            uold2[count] = val[0];
            break;
          case 1:
            uoldx[count+2] = x;
            uoldy[count+2] = y;
            u1->FindValueLocal(cell,j,x,y,val);
            uold1[count+2] = val[0];
            u2->FindValueLocal(cell,j,x,y,val);
            uold2[count+2] = val[0];
            break;
          case 2:
            uoldx[count+5] = x;
            uoldy[count+5] = y;
            u1->FindValueLocal(cell,j,x,y,val);
            uold1[count+5] = val[0];
            u2->FindValueLocal(cell,j,x,y,val);
            uold2[count+5] = val[0];
            break;
        }
      }
      // midpoint of edges
      uoldx[count+1] = (uoldx[count] +uoldx[count+2])/2.0;
      uoldy[count+1] = (uoldy[count] +uoldy[count+2])/2.0;
      u1->FindValueLocal(cell,j,uoldx[count+1],uoldy[count+1],val);
      uold1[count+1] = val[0];
      u2->FindValueLocal(cell,j,uoldx[count+1],uoldy[count+1],val);
      uold2[count+1] = val[0];

      uoldx[count+4] = (uoldx[count+2] +uoldx[count+5])/2.0;
      uoldy[count+4] = (uoldy[count+2] +uoldy[count+5])/2.0;
      u1->FindValueLocal(cell,j,uoldx[count+4],uoldy[count+4],val);
      uold1[count+4] = val[0];
      u2->FindValueLocal(cell,j,uoldx[count+4],uoldy[count+4],val);
      uold2[count+4] = val[0];

      uoldx[count+3] = (uoldx[count] +uoldx[count+5])/2.0;
      uoldy[count+3] = (uoldy[count]+uoldy[count+5])/2.0;
      u1->FindValueLocal(cell,j,uoldx[count+3],uoldy[count+3],val);
      uold1[count+3] = val[0];
      u2->FindValueLocal(cell,j,uoldx[count+3],uoldy[count+3],val);
      uold2[count+3] = val[0];

      // bubble part
      uoldx[count+6] = (uoldx[count] +uoldx[count+2] +uoldx[count+5])/3.0;
      uoldy[count+6] = (uoldy[count] +uoldy[count+2] +uoldy[count+5])/3.0;
      index = DOF[6];
      uold1[count+6] = sol[index];
      uold2[count+6] = sol[index+N_U];

      count+=7;
    }

    SwapEdges(name, coll, u1, u2, tangential_values_ssl);
    swapped = 1;
  }

  angle +=  phi_plus;

  // update boundary description
  // ONLY FOR ROTATING FLOW PROBLEM
  /* for (j=0;j<N_Cells;j++)
  {
  cell = coll->GetCell(j);
  N_V = cell->GetN_Vertices();
  for (k=0;k<N_V;k++)
  {
    Vertices[k] = cell->GetVertex(k);
  }
  for (k=0;k<N_V;k++)
  {
    joint = cell->GetJoint(k);
  if(joint->GetType() == BoundaryEdge ||
  joint->GetType() == IsoBoundEdge)
  {
  if(joint->GetType() == BoundaryEdge)
  {
  boundedge = (TBoundEdge *)joint;
  BoundComp = boundedge->GetBoundComp();
  ((TBdCircle*)BoundComp)->SetParams(0,0,1,1,-angle,2*Pi-angle);
  }
  else
  {
  isoboundedge = (TIsoBoundEdge *)joint;
  BoundComp = isoboundedge->GetBoundComp();
  ((TBdCircle*)BoundComp)->SetParams(0,0,1,1,-angle,2*Pi-angle);
  }

  }
  }
  }*/

  // update parameters on the boundary if necessary
  // FOR SLIT PROBLEM
  for (j=0;j<N_Cells;j++)
  {
    cell = coll->GetCell(j);
    N_V = cell->GetN_Vertices();
    for (k=0;k<N_V;k++)
    {
      Vertices[k] = cell->GetVertex(k);
    }
    for (k=0;k<N_V;k++)
    {
      joint = cell->GetJoint(k);
      if(joint->GetType() == BoundaryEdge ||
        joint->GetType() == IsoBoundEdge)
      {
        if(joint->GetType() == BoundaryEdge)
        {
          boundedge = (TBoundEdge *)joint;
          boundedge->GetParameters(t0, t1);
          //OutPut(t0 << " " << t1 << " ");
          //((TBoundEdge *)joint)->UpdateParameters(Vertices[k],
          //Vertices[(k+1)%N_V]);
        }
        else
        {
          boundedge = (TBoundEdge *)joint;
          boundedge->GetParameters(t0, t1);
          //OutPut(t0 << " " << t1 << " ");
          //((TInterfaceJoint *)joint)->UpdateParameters(Vertices[k],
          //					 Vertices[(k+1)%N_V]);
        }

      }
    }
  }

  return swapped;
}                                                 // end i (degree)


/*******************************************************************************/
//
// VelocityInNewPositions
// computes the velicity values in the position after the rotation
//
/*******************************************************************************/

int VelocityAtNewPositions(TCollection *coll,
TFEFunction2D *u1, TFEFunction2D *u2,
double *values)
{
  int i, j, jj, k, kk, N_Cells, N_V, swapped=0, count, found;
  int *GlobalNumbers_velo, *BeginIndex_velo, *DOF, index, N_U;
  double x, y, phi, r, phi_plus, t0, t1;
  double val[1], *sol, xp[3], yp[3];
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  TBaseCell *cell, *cell1, *cell_n, *cell_nk;
  TVertex *vertex, *Vertices[4];
  TJoint *joint;
  TBoundComp *BoundComp;
  TBoundEdge *boundedge;
  TIsoBoundEdge *isoboundedge;
  TFESpace2D *velo_space;

  phi_plus = 2*Pi * TDatabase::TimeDB->TIMESTEPLENGTH * TDatabase::ParamDB->SSMUM_ROT_PER_SECOND;

  // velocity space
  velo_space = u1->GetFESpace2D();
  // number of velocity unknowns
  N_U = u1->GetLength();
  // vector for velocity
  sol = u1->GetValues();
  // information to the degrees of freedom for the velocity
  GlobalNumbers_velo = velo_space->GetGlobalNumbers();
  BeginIndex_velo = velo_space->GetBeginIndex();

  // number of mesh cells
  N_Cells = coll->GetN_Cells();

  // set clip board of the vertices and cells
  for (j=0;j<N_Cells;j++)
  {
    cell = coll->GetCell(j);
    cell->SetClipBoard(j);
  }

  // compute exact valus of velocity in the new coordinates
  // for dof on vertices and mitpoints of edges
  for (j=0;j<N_Cells;j++)
  {
    // mesh cell
    cell = coll->GetCell(j);
    // number of vertices (same as number of edges)
    N_V = cell->GetN_Vertices();
    for (k=0;k<N_V;k++)
    {
      // compute current coordinates
      vertex =  cell->GetVertex(k);
      vertex->GetCoords(xp[k], yp[k]);
    }
    DOF =  GlobalNumbers_velo + BeginIndex_velo[j];
    for (k=0;k<6;k++)
    {
      switch(k)
      {
        case 0:
          x = xp[0];
          y = yp[0];
          break;
        case 1:
          x = (xp[0]+xp[1])/2.0;
          y = (yp[0]+yp[1])/2.0;
          break;
        case 2:
          x = xp[1];
          y = yp[1];
          break;
        case 3:
          x = (xp[0]+xp[2])/2.0;
          y = (yp[0]+yp[2])/2.0;
          break;
        case 4:
          x = (xp[1]+xp[2])/2.0;
          y = (yp[1]+yp[2])/2.0;
          break;
        case 5:
          x = xp[2];
          y = yp[2];
          break;
      }
      // compute coordinates after rotation
      r = (x - mp_x)*(x - mp_x) + (y - mp_y)*(y - mp_y);
      r = sqrt(r);
      // vertex in outer domain
      if (r>(TDatabase::ParamDB->SSMUM_OUTER_RADIUS+
        TDatabase::ParamDB->SSMUM_INNER_RADIUS)/2.0)
        continue;
      index = DOF[k];

      cell1 = NULL;
      // get argument of the vertex
      phi = atan2(y-mp_y,x-mp_x);
      //OutPut(x << " " << y << " " << phi << endl);
      phi -= phi_plus;
      x = r * cos(phi) + mp_x;
      y = r * sin(phi) + mp_y;
      // find mesh cell for new position of the vertex
      // case 1, same mesh cell
      if (PointInCell(cell,x,y))
      {
        cell1 = cell;
      }
      else
      {
        //case 2, check neighbour cell
        for (jj=0;jj<N_V;jj++)
        {
          // neighbour cell
          cell_n=cell->GetJoint(jj)->GetNeighbour(cell);
          // boundary
          if (cell_n ==  NULL)
            continue;
          if (PointInCell(cell_n,x,y))
          {
            cell1 = cell_n;
            break;
          }
        }
        // check the neighbours of the neighbours
        if (cell1 == NULL)
        {
          for (jj=0;jj<N_V;jj++)
          {
            // neighbour cell
            cell_n=cell->GetJoint(jj)->GetNeighbour(cell);
            // boundary
            if (cell_n ==  NULL)
              continue;
            for (kk=0;kk<N_V;kk++)
            {
              cell_nk=cell_n->GetJoint(kk)->GetNeighbour(cell_n);
              // boundary
              if (cell_nk ==  NULL)
                continue;
              if (PointInCell(cell_nk,x,y))
              {
                cell1 = cell_nk;
                break;
              }
            }
            if (cell1 != NULL)
              break;
          }                                       // end jj
          // brute force
          for (jj=0;jj<N_Cells;jj++)
          {
            cell_n =  coll->GetCell(jj);
            if (PointInCell(cell_n,x,y))
            {
              cell1 = cell_n;
              break;
            }
          }
        }
      }

      // some points on the slit boundary do not find a mesh cell on their
      // new position
      if (cell1!=NULL)
      {
        u1->FindValueLocal(cell1,cell1->GetClipBoard(),x,y,val);
        //OutPut("v " << val[0] << " " << sol[index] << endl);
        values[index] = val[0];
        u2->FindValueLocal(cell1,cell1->GetClipBoard(),x,y,val);
        //OutPut(val[0] << " " << sol[index+N_U] << endl);
        values[index+N_U] = val[0];
      }
    }
  }                                               // end j
}


void FillNewVelocity(TCollection *coll,
double *uoldx, double *uoldy, double *uold1, double *uold2,
TFEFunction2D *u1, TFEFunction2D *u2, double *tangential_values_ssl)
{
  int i, j, k, found, N_Cells, N_V, N_U, index, index0, index1, count;
  int found0, found1, found2, found_tan;
  int cell_number;
  int *GlobalNumbers_velo, *BeginIndex_velo, *DOF, *DOF_n;
  double x, y, xx[3], yy[3], eps=1e-8, *sol, val1, val2, x_m,y_m;
  double L,n_x,n_y,T_1,T_2,T_3,I_0,I_1,I_2, t_x, t_y,l,det,P,Q;
  double tangential_value, tangential_value0, r, r_n;
  double u1_unb, u2_unb;
  double mp_x = TDatabase::ParamDB->SSMUM_MP_X;
  double mp_y = TDatabase::ParamDB->SSMUM_MP_Y;
  TBaseCell *cell, *cell_n;
  TVertex *vertex;
  TFESpace2D *velo_space;

  // number of cells
  N_Cells = coll->GetN_Cells();
  // velocity space
  velo_space = u1->GetFESpace2D();
  // number of velocity unknowns
  N_U = u1->GetLength();
  // vector for velocity
  sol = u1->GetValues();
  // information to the degrees of freedom for the velocity
  GlobalNumbers_velo = velo_space->GetGlobalNumbers();
  BeginIndex_velo = velo_space->GetBeginIndex();

  // cells of the new grid
  for (j=0;j<N_Cells;j++)
  {
    // get cell
    cell = coll->GetCell(j);
    // dof belonging to this mesh cell
    DOF =  GlobalNumbers_velo + BeginIndex_velo[j];
    // number of vertices
    N_V = cell->GetN_Vertices();
    // find coordinates corresponding to the degrees of freedom
    // vertices
    for (k=0;k<N_V;k++)
    {
      // get vertex
      vertex =  cell->GetVertex(k);
      // get coordinates of the vertex
      vertex->GetCoords(xx[k], yy[k]);
      // check if vertices are in the list
      found = 0;
      // loop over the information wrt the solution on the old grid
      for (i=0;i<7*N_Cells;i++)
      {
        // coordinates are the same
        if ((fabs(xx[k]-uoldx[i])<eps)&&(fabs(yy[k]-uoldy[i])<eps))
        {
          // set flag
          found = 1;
          // which vertex
          switch(k)
          {
            case 0:
              index = DOF[0];
              sol[index] = uold1[i];              //sol in neue Gitter- erste componente
              sol[index+N_U] = uold2[i];          //sol in neue Gitter- zweite componente
              break;
            case 1:
              index = DOF[2];
              sol[index] = uold1[i];
              sol[index+N_U] = uold2[i];
              break;
            case 2:
              index = DOF[5];
              sol[index] = uold1[i];
              sol[index+N_U] = uold2[i];
              break;
          }
          break;
        }
      }
      if (!found)
      {
        OutPut("Vertex not found!!!"<<endl);
        exit(4711);
      }
    }
    // midpoints of the edges
    x = (xx[0]+xx[1])/2.0;
    y = (yy[0]+yy[1])/2.0;
    found0 = 0;
    // edge mid points which did not change
    for (i=0;i<7*N_Cells;i++)
    {
      if ((fabs(x-uoldx[i])<eps)&&(fabs(y-uoldy[i])<eps))
      {
        found0 = 1;
        index = DOF[1];
        sol[index] = uold1[i];
        sol[index+N_U] = uold2[i];
        break;
      }
    }
    x = (xx[1]+xx[2])/2.0;
    y = (yy[1]+yy[2])/2.0;
    found1 = 0;
    for (i=0;i<7*N_Cells;i++)
    {
      if ((fabs(x-uoldx[i])<eps)&&(fabs(y-uoldy[i])<eps))
      {
        index = DOF[4];
        sol[index] = uold1[i];
        sol[index+N_U] = uold2[i];
        found1 = 1;
        break;
      }
    }
    x = (xx[0]+xx[2])/2.0;
    y = (yy[0]+yy[2])/2.0;
    found2 = 0;
    for (i=0;i<7*N_Cells;i++)
    {
      if ((fabs(x-uoldx[i])<eps)&&(fabs(y-uoldy[i])<eps))
      {
        index = DOF[3];
        sol[index] = uold1[i];
        sol[index+N_U] = uold2[i];
        found2 = 1;
        break;
      }
    }

    // edge mid points which changed
    if ((!found0)||(!found1)||(!found2))
    {
      switch(TDatabase::ParamDB->SSMUM_INTERPOLATION)
      {
        case 0:
        {
          /*sol[DOF[0]]= sol[DOF[2]] = sol[DOF[3]]=sol[DOF[4]] = sol[DOF[5]] = 1;
          sol[DOF[0]+N_U]= sol[DOF[2]+N_U] = sol[DOF[3]+N_U]=sol[DOF[4]+N_U] = sol[DOF[5]+N_U] = 2;
          xx[0] = 0;
          yy[0] = 0;
          xx[1] = 1;
          yy[1] = 0;
          xx[2] = 0;
          yy[2] = 1;
          */

          // divergence constraint
          if (!found0)
          {
            // edge from vertex[0] to vertex[1]
            index = DOF[1];
            sol[index] = sol[index+N_U] = 0;
          }
          if (!found1)
          {
            // edge from vertex[0] to vertex[1]
            index = DOF[4];
            sol[index] = sol[index+N_U] = 0;
          }
          if (!found2)
          {
            // edge from vertex[0] to vertex[1]
            index = DOF[3];
            sol[index] = sol[index+N_U] = 0;
          }

          //  compute outflow with currently available values, use Simpson rule
          // edge from vertex[0] to vertex[1]
          index = DOF[1];
          index0 = DOF[0];
          index1 = DOF[2];
          L=sqrt((xx[0]-xx[1])*(xx[0]-xx[1])+(yy[0]-yy[1])*(yy[0]-yy[1]));
          n_x=yy[1]-yy[0];
          n_y=xx[0]-xx[1];
          n_x /= L;
          n_y /= L;
          T_1=sol[index0]*n_x+sol[index0+N_U]*n_y;
          T_2=4*(sol[index]*n_x+sol[index+N_U]*n_y);
          T_3=sol[index1]*n_x+sol[index1+N_U]*n_y;
          I_0=L*(T_1+T_2+T_3)/6;

          // edge from vertex[1] to vertex[2]
          index = DOF[4];
          index0 = DOF[2];
          index1 = DOF[5];
          L=sqrt((xx[2]-xx[1])*(xx[2]-xx[1])+(yy[2]-yy[1])*(yy[2]-yy[1]));
          n_x=yy[2]-yy[1];
          n_y=xx[1]-xx[2];
          n_x /= L;
          n_y /= L;
          T_1=sol[index0]*n_x+sol[index0+N_U]*n_y;
          T_2=4*(sol[index]*n_x+sol[index+N_U]*n_y);
          T_3=sol[index1]*n_x+sol[index1+N_U]*n_y;
          I_0 += L*(T_1+T_2+T_3)/6;

          // edge from vertex[2] to vertex[0]
          index = DOF[3];
          index0 = DOF[5];
          index1 = DOF[0];
          L=sqrt((xx[2]-xx[0])*(xx[2]-xx[0])+(yy[2]-yy[0])*(yy[2]-yy[0]));
          n_x=yy[0]-yy[2];
          n_y=xx[2]-xx[0];
          n_x /= L;
          n_y /= L;
          T_1=sol[index0]*n_x+sol[index0+N_U]*n_y;
          T_2=4*(sol[index]*n_x+sol[index+N_U]*n_y);
          T_3=sol[index1]*n_x+sol[index1+N_U]*n_y;
          I_0 += L*(T_1+T_2+T_3)/6;

          // edge from vertex[0] to vertex[1]
          if (!found0)
          {
            index = DOF[1];
            // indices of dof in vertices
            index0 = DOF[0];
            index1 = DOF[2];
            // length of the edge
            l=sqrt((xx[0]-xx[1])*(xx[0]-xx[1])+(yy[0]-yy[1])*(yy[0]-yy[1]));
            // tangential from vertex[0] to vertex[1]
            t_x=xx[1]-xx[0];
            t_y=yy[1]-yy[0];
            // norm
            t_x /=l;
            t_y /=l;
            // normal
            n_x = t_y;
            n_y = -t_x;
            // check direction, orientation always from interior to outside
            r = (xx[0]-mp_x)*(xx[0]-mp_x) + (yy[0]-mp_y)*(yy[0]-mp_y);
            r_n = (xx[1]-mp_x)*(xx[1]-mp_x) + (yy[1]-mp_y)*(yy[1]-mp_y);
            if (r_n<r)
            {
              t_x = - t_x;
              t_y = - t_y;
            }
            // find value of integrated tangential velocity which has
            // been stored before swapping the edge
            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              // check coordinates
              if ( fabs(tangential_values_ssl[3*i]-xx[0]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[0]) < eps)
                {
                  // found
                  tangential_value = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }

            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[3*i]-xx[1]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[1]) < eps)
                {
                  tangential_value0 = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
            if (fabs(tangential_value-tangential_value0)>eps)
            {
              OutPut("Wrong tangential value !!! " << fabs(tangential_value-tangential_value0)  << endl);
              exit(4711);
            }
          }

          if (!found1)
          {
            index = DOF[4];
            index0 = DOF[2];
            index1 = DOF[5];
            l=sqrt((xx[2]-xx[1])*(xx[2]-xx[1])+(yy[2]-yy[1])*(yy[2]-yy[1]));
            t_x=xx[2]-xx[1];
            t_y=yy[2]-yy[1];
            t_x /=l;
            t_y /=l;
            n_x=t_y;
            n_y=-t_x;
            // check direction, orientation always from interior to outside
            r = (xx[1]-mp_x)*(xx[1]-mp_x) + (yy[1]-mp_y)*(yy[1]-mp_y);
            r_n = (xx[2]-mp_x)*(xx[2]-mp_x) + (yy[2]-mp_y)*(yy[2]-mp_y);
            if (r_n<r)
            {
              t_x = - t_x;
              t_y = - t_y;
            }

            // find value of integrated tangential velocity which has
            // been stored before swapping the edge
            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[3*i]-xx[1]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[1]) < eps)
                {
                  tangential_value = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }

            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[3*i]-xx[2]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[2]) < eps)
                {
                  tangential_value0 = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
            if (fabs(tangential_value-tangential_value0)>eps)
            {
              OutPut("Wrong tangential value !!! " << fabs(tangential_value-tangential_value0)  << endl);
              exit(4711);
            }
          }
          if (!found2)
          {
            index = DOF[3];
            index0 = DOF[0];
            index1 = DOF[5];
            l=sqrt((xx[0]-xx[2])*(xx[0]-xx[2])+(yy[0]-yy[2])*(yy[0]-yy[2]));
            t_x=xx[0]-xx[2];
            t_y=yy[0]-yy[2];
            t_x /=l;
            t_y /=l;
            n_x=t_y;
            n_y=-t_x;
            // check direction, orientation always from interior to outside
            r = (xx[2]-mp_x)*(xx[2]-mp_x) + (yy[2]-mp_y)*(yy[2]-mp_y);
            r_n = (xx[0]-mp_x)*(xx[0]-mp_x) + (yy[0]-mp_y)*(yy[0]-mp_y);
            if (r_n<r)
            {
              t_x = - t_x;
              t_y = - t_y;
            }

            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[3*i]-xx[0]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[0]) < eps)
                {
                  tangential_value = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }

            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[3*i]-xx[2]) < eps)
              {
                if ( fabs(tangential_values_ssl[3*i+1]-yy[2]) < eps)
                {
                  tangential_value0 = tangential_values_ssl[3*i+2];
                  found_tan++;
                }
              }
            }
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
            if (fabs(tangential_value-tangential_value0)>eps)
            {
              OutPut("Wrong tangential value !!! " << fabs(tangential_value-tangential_value0)  << endl);
              exit(4711);
            }
          }

          // set rhs for 2x2 system for determining the values in the
          // midpoint of the edge
          P = -3*I_0/(2*l);
          Q = 6*tangential_value/l;
          Q -= sol[index0]*t_x + sol[index0+N_U]*t_y;
          Q -= sol[index1]*t_x + sol[index1+N_U]*t_y;
          Q /=4;
          // solve system
          det = n_x*t_y-n_y*t_x;
          //	OutPut((xx[0]+xx[1])/2 << " " << (yy[0]+yy[1])/2 << " :: " <<
          //       sol[index] << " "
          //       << sol[index+N_U] << " :: " << r << " " << r_n << endl);
          sol[index] = (P*t_y-Q*n_y)/det;
          sol[index+N_U] = (Q*n_x-P*t_x)/det;

          //OutPut(j << " " <<
          //	   sol[index] << " " << (sol[index0]+sol[index1])/2 << " : "
          //	   << sol[index+N_U] << " " <<  (sol[index0+N_U]+sol[index1+N_U])/2 << endl);

          /*if (j==1144)
            {
            OutPut(sol[DOF[0]] << " " << sol[DOF[2]]  << " "
            << " " << sol[DOF[3]] << " " << sol[DOF[4]] << " " << sol[DOF[5]]<< endl);
            OutPut(sol[DOF[0]+N_U] << " " << sol[DOF[2]+N_U]  << " "
            << " " << sol[DOF[3]+N_U] << " " << sol[DOF[4]+N_U] << " " << sol[DOF[5]+N_U]<< endl <<endl);
            OutPut(t_x << " " << t_y << " " <<
               l*(sol[DOF[0]]*t_x + sol[DOF[0]+N_U]*t_y
               + 4*sol[DOF[1]]*t_x + 4*sol[DOF[1]+N_U]*t_y
               +sol[DOF[2]]*t_x + sol[DOF[2]+N_U]*t_y)/6 << endl);
          */
          if (0)
            OutPut(t_x << " " << t_y << " " << l << " " << r << " " << r_n << endl);
          /*	       l*(sol[DOF[0]]*t_x + sol[DOF[0]+N_U]*t_y
               + 4*(-0.0611985)*t_x + 4*(0.117472)*t_y
               +sol[DOF[2]]*t_x + sol[DOF[2]+N_U]*t_y)/6 << endl);
          */
          break;
        }
        case 1:
        {
          // edge from vertex[0] to vertex[1]
          if (!found0)
          {
            index = DOF[1];
            x_m = (xx[0]+xx[1])/2;
            y_m = (yy[0]+yy[1])/2;
            // find value of integrated tangential velocity which has
            // been stored before swapping the edge
            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              // check coordinates
              if ( fabs(tangential_values_ssl[4*i]-x_m) < eps)
              {
                if ( fabs(tangential_values_ssl[4*i+1]-y_m) < eps)
                {
                  // found
                  sol[index] = tangential_values_ssl[4*i+2];
                  sol[index+N_U] = tangential_values_ssl[4*i+3];
                  found_tan++;
                }
              }
            }
            OutPut(sol[index] << " " << (sol[DOF[0]]+sol[DOF[2]])/2 << " "
              << sol[index + N_U] << " " << (sol[DOF[0]+N_U]+sol[DOF[2]+N_U])/2 << endl);
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
          }

          if (!found1)
          {
            index = DOF[4];
            x_m = (xx[1]+xx[2])/2;
            y_m = (yy[1]+yy[2])/2;

            // find value of integrated tangential velocity which has
            // been stored before swapping the edge
            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[4*i]-x_m) < eps)
              {
                if ( fabs(tangential_values_ssl[4*i+1]-y_m) < eps)
                {

                  sol[index] = tangential_values_ssl[4*i+2];
                  sol[index+N_U] = tangential_values_ssl[4*i+3];
                  found_tan++;
                }
              }
            }
            OutPut(sol[index] << " " << (sol[DOF[5]]+sol[DOF[2]])/2 << " "
              << sol[index + N_U] << " " << (sol[DOF[5]+N_U]+sol[DOF[2]+N_U])/2 << endl);
            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
          }
          if (!found2)
          {
            index = DOF[3];
            x_m = (xx[2]+xx[0])/2;
            y_m = (yy[2]+yy[0])/2;

            found_tan = 0;
            for (i=0;i<TDatabase::ParamDB->SSMUM_MAX_CELLS_LAYERS;i++)
            {
              if ( fabs(tangential_values_ssl[4*i]-x_m) < eps)
              {
                if ( fabs(tangential_values_ssl[4*i+1]-y_m) < eps)
                {

                  sol[index] = tangential_values_ssl[4*i+2];
                  sol[index+N_U] = tangential_values_ssl[4*i+3];
                  found_tan++;
                }
              }
            }

            OutPut(sol[index] << " " << (sol[DOF[0]]+sol[DOF[5]])/2 << " "
              << sol[index + N_U] << " " << (sol[DOF[0]+N_U]+sol[DOF[5]+N_U])/2 << endl);

            if (found_tan!=1)
            {
              OutPut("Error in finding tangential value !!! " << found_tan << endl);
              exit(4711);
            }
          }
          break;
        }
        case 2:
        {
          if (!found0)
          {
            //OutPut("Midpoint 0 not found!!!"<<endl);
            // linear interpolation
            index = DOF[1];
            index0 = DOF[0];
            index1 = DOF[2];
            sol[index] = (sol[index0]+sol[index1])/2.0;
            sol[index+N_U] = (sol[index0+N_U]+sol[index1+N_U])/2.0;
          }
          if (!found1)
          {
            //OutPut("Midpoint 0 not found!!!"<<endl);
            // linear interpolation
            index = DOF[4];
            index0 = DOF[2];
            index1 = DOF[5];
            sol[index] = (sol[index0]+sol[index1])/2.0;
            sol[index+N_U] = (sol[index0+N_U]+sol[index1+N_U])/2.0;
          }
          if (!found2)
          {
            //OutPut("Midpoint 0 not found!!!"<<endl);
            // linear interpolation
            index = DOF[3];
            index0 = DOF[0];
            index1 = DOF[5];
            sol[index] = (sol[index0]+sol[index1])/2.0;
            sol[index+N_U] = (sol[index0+N_U]+sol[index1+N_U])/2.0;
          }
          break;
        }
      }                                           // edge midpoint not found

      // bubble part
      x = (xx[0]+xx[1]+xx[2])/3.0;
      y = (yy[0]+yy[1]+yy[2])/3.0;
      found = 0;
      for (i=0;i<7*N_Cells;i++)
      {
        if ((fabs(x-uoldx[i])<eps)&&(fabs(y-uoldy[i])<eps))
        {
          index = DOF[6];
          sol[index] = uold1[i];
          sol[index+N_U] = uold2[i];
          found = 1;
          break;
        }
      }
      if (!found)
      {
        index = DOF[6];
        //OutPut(j << " cannot recover bubble part !!! "<<index << endl);
        // linear interpolation
        //index = DOF[3];
        //index0 = DOF[0];
        //index1 = DOF[5];
        sol[index] = -4711;
        //sol[index] = (sol[index0]+sol[index1])/2.0;
        //sol[index+N_U] = (sol[index0+N_U]+sol[index1+N_U])/2.0;
      }

    }
    // missing bubble parts by averaging with neighbors
    // initialize clipboars
    for(i=0;i<N_Cells;i++)
    {
      cell=coll->GetCell(i);
      cell->SetClipBoard(i);
    }
    // cells of the new grid
    for (j=0;j<N_Cells;j++)
    {
      // get cell
      cell = coll->GetCell(j);
      // dof belonging to this mesh cell
      DOF =  GlobalNumbers_velo + BeginIndex_velo[j];
      // bubble index
      index = DOF[6];
      if (sol[index]!=-4711)
        continue;
      // cell in the shear slip layer
      N_V = cell->GetN_Vertices();
      val1 = val2 = 0;
      count = 0;
      for (k=0;k<N_V;k++)
      {
        // get neighbor
        cell_n=cell->GetJoint(k)->GetNeighbour(cell);
        cell_number = cell_n->GetClipBoard();
        DOF_n =  GlobalNumbers_velo + BeginIndex_velo[cell_number];
        index1 = DOF_n[6];
        // neighbor has also no bubble part
        if (sol[index1] == -4711)
          continue;
        count++;
        val1 += sol[index1];
        val2 += sol[index1+N_U];
      }
      sol[index] = val1/count;
      sol[index+N_U] = val2/count;
    }
  }
}

void MakeBubblesDivFree(TCollection *coll,
    TFEFunction2D *u1, TFEFunction2D *u2,
    TFEFunction2D *p,
    TMatrix2D *B1, TMatrix2D *B2)
  {
    int i, ii, j, N_Cells, N_U, N_P, index, index1, row;
    int *GlobalNumbers_velo, *BeginIndex_velo, *DOF_velo;
    int *GlobalNumbers_press, *BeginIndex_press, *DOF_press;
    int *RowPtr, *KCol;
    double *sol, loc_rhs[2], loc_a[2][2], *Entries, det, bub1, bub2;
    TBaseCell *cell;
    TFESpace2D *velo_space, *press_space;

    // number of cells
    N_Cells = coll->GetN_Cells();
    // velocity space
    velo_space = u1->GetFESpace2D();
    // number of velocity unknowns
    N_U = u1->GetLength();
    // vector for velocity
    sol = u1->GetValues();
    // information to the degrees of freedom for the velocity
    GlobalNumbers_velo = velo_space->GetGlobalNumbers();
    BeginIndex_velo = velo_space->GetBeginIndex();
    // pressure space
    press_space = p->GetFESpace2D();
    // number of pressure unknowns
    N_P = p->GetLength();
    // information to the degrees of freedom for the pressure
    GlobalNumbers_press = press_space->GetGlobalNumbers();
    BeginIndex_press = press_space->GetBeginIndex();

    // loop over the cells
    for (j=0;j<N_Cells;j++)
    {
      // get cell
      cell = coll->GetCell(j);
      // dof belonging to this mesh cell
      DOF_velo =  GlobalNumbers_velo + BeginIndex_velo[j];
      DOF_press =  GlobalNumbers_press + BeginIndex_press[j];
      // compute rhs and matrix for local problem which
      // corresponds to the divergence constraint for the bubble
      // functions, the other values at the mesh cell are assumed
      // to be already correct
      // only the pressure test functions DOF_press[1] and DOF_press[2]
      // have an impact on the bubble
      // the bubble function is DOF_velo[6]
      loc_rhs[0] = loc_rhs[1] = 0.0;
      // consider matrix B1
      RowPtr = B1->GetRowPtr();
      KCol = B1->GetKCol();
      Entries = B1->GetEntries();
      // first pressure test function
      row = DOF_press[1];
      // check the corresponding row of B1
      for (i=RowPtr[row];i<RowPtr[row+1];i++)
      {
        index = KCol[i];
        for (ii=0;ii<7;ii++)
        {
          // column corresponds to one of the velo dof
          if (index == DOF_velo[ii])
          {
            // the bubble function
            if (ii==6)
            {
              loc_a[0][0] = Entries[i];
            }
            else
            {
              index1 = DOF_velo[ii];
              loc_rhs[0] -= Entries[i] * sol[index1];
            }
          }
        }
      }                                           // end first pressure test function, B1
      // second pressure test function
      row = DOF_press[2];
      // check the corresponding row of B1
      for (i=RowPtr[row];i<RowPtr[row+1];i++)
      {
        index = KCol[i];
        for (ii=0;ii<7;ii++)
        {
          // column corresponds to one of the velo dof
          if (index == DOF_velo[ii])
          {
            // the bubble function
            if (ii==6)
            {
              loc_a[1][0] = Entries[i];
            }
            else
            {
              index1 = DOF_velo[ii];
              loc_rhs[1] -= Entries[i] * sol[index1];
            }
          }
        }
      }                                           // end second pressure test function, B1

      // consider matrix B2
      RowPtr = B2->GetRowPtr();
      KCol = B2->GetKCol();
      Entries = B2->GetEntries();
      // first pressure test function
      row = DOF_press[1];
      // check the corresponding row of B2
      for (i=RowPtr[row];i<RowPtr[row+1];i++)
      {
        index = KCol[i];
        for (ii=0;ii<7;ii++)
        {
          // column corresponds to one of the velo dof
          if (index == DOF_velo[ii])
          {
            // the bubble function
            if (ii==6)
            {
              loc_a[0][1] = Entries[i];
            }
            else
            {
              index1 = DOF_velo[ii] + N_U;
              loc_rhs[0] -= Entries[i] * sol[index1];
            }
          }
        }
      }                                           // end first pressure test function, B2
      // second pressure test function
      row = DOF_press[2];
      // check the corresponding row of B2
      for (i=RowPtr[row];i<RowPtr[row+1];i++)
      {
        index = KCol[i];
        for (ii=0;ii<7;ii++)
        {
          // column corresponds to one of the velo dof
          if (index == DOF_velo[ii])
          {
            // the bubble function
            if (ii==6)
            {
              loc_a[1][1] = Entries[i];
            }
            else
            {
              index1 = DOF_velo[ii] + N_U;
              loc_rhs[1] -= Entries[i] * sol[index1];
            }
          }
        }
      }                                           // end second pressure test function, B2
      // solve local system
      det = loc_a[0][0] * loc_a[1][1] - loc_a[1][0] * loc_a[0][1];
      bub1 = loc_rhs[0] * loc_a[1][1] - loc_rhs[1] * loc_a[0][1];
      bub1 /= det;
      bub2 = loc_a[0][0] * loc_rhs[1] - loc_a[1][0] * loc_rhs[0] ;
      bub2 /= det;
      //OutPut(sol[DOF_velo[6]] << " " << bub1 << " : " << sol[DOF_velo[6]+N_U] << " " << bub2 << endl);
      sol[DOF_velo[6]] = bub1;
      sol[DOF_velo[6]+N_U] = bub2;
    }
    //OutPut("done"<<endl);
  }

  
#endif // #ifdef __2D__  
