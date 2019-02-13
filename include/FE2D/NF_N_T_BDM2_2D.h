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
   
/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
  used function names in this file must be unique in the entire program. 
  -> use the identifier as prefix
*/

static double NF_N_T_BDM2_2D_Xi[] =  {0};
static double NF_N_T_BDM2_2D_Eta[] = {0};
// NOTE: If you want to use other evaluation points for degress of freedom on
// the edges of a cell, you also have to change basis functions in 
// BF_N_T_BDM2_2D.h
//static double NF_N_T_BDM2_2D_T[] = {-0.5,0,0.5};// equidistant points
//static double NF_N_T_BDM2_2D_T[] = {-0.774596669241483,0,0.774596669241483};//Gauss-points
static double NF_N_T_BDM2_2D_T[] = {-0.866025403784439,0,0.866025403784439};//Tschebyscheff-points

void NF_N_T_BDM2_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
//   static double weights[3] = { 0.5555555555555555555555555555555556,
//                                0.88888888888888888888888888888888889,
//                                0.5555555555555555555555555555555556 };
//   Functionals[0] = ( weights[0]*PointValues[0]
//                     +weights[1]*PointValues[1]
//                     +weights[2]*PointValues[2]) * 0.5;
//   Functionals[1] = ( weights[0]*PointValues[3]
//                     +weights[1]*PointValues[4]
//                     +weights[2]*PointValues[5]) * 0.5;
//   Functionals[2] = ( weights[0]*PointValues[6]
//                     +weights[1]*PointValues[7]
//                     +weights[2]*PointValues[8]) * 0.5;
//   Functionals[3] = ( weights[0]*PointValues[9]
//                     +weights[1]*PointValues[10]
//                     +weights[2]*PointValues[11]) * 0.5;
cout << "Raviart-Thomas elements of order 2 on triangles: "
     << "Nodal functionals are not implemented properly!" << endl;
  for(int i=0; i<13; i++)
    Functionals[i] = PointValues[0];
  
}

void NF_N_T_BDM2_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                            double *PointValues, double *Functionals)
{
// this is needed for setting boundary conditions
  /* the functionals
   * int_Joint v.n q_1  and  int_Joint v.n q_2  and  int_Joint v.n q_3
   * (q_1, q2 and q_3 are three linearly independent polynomials of degree 2)
   * will be multiplied by the length of the Joint (edge). Otherwise one would
   * ensure int_Joint v.n=PointValues[0]. 
   * Example: If you would like to have u.n=1, then without multiplying by 
   *          the edge length l would result in having int_Joint u.n=1 on each
   *          boundary edge. This would mean one gets u.n=1/l on that 
   *          boundary. To avoid this, we introduce the factor l here. 
   * However I am not sure if this causes trouble elsewhere later. 
   * Be carefull!
   *                                            Ulrich Wilbrandt, 11.05.2012
  */
  double l; // length of joint
  double x0,x1,y0,y1;
  #ifdef __2D__
  Cell->GetVertex(Joint)->GetCoords(x0,y0);
  Cell->GetVertex((Joint+1)%3)->GetCoords(x1,y1);// 3=number of edges
  #endif
  l = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1));
  Functionals[0] = PointValues[0]*l;
  Functionals[1] = PointValues[1]*l;
  Functionals[2] = PointValues[2]*l;
}

TNodalFunctional2D *NF_N_T_BDM2_2D_Obj = new TNodalFunctional2D
        (NF_N_T_BDM2_2D, 12, 3, 1, 3, NF_N_T_BDM2_2D_Xi, NF_N_T_BDM2_2D_Eta,
         NF_N_T_BDM2_2D_T, NF_N_T_BDM2_2D_EvalAll, NF_N_T_BDM2_2D_EvalEdge);
