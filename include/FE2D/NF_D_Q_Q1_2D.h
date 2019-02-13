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
   
static double NF_D_Q_Q1_2D_Xi[9] = {
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704,
   -0.77459666924148337704, 0, 0.77459666924148337704 };

static double NF_D_Q_Q1_2D_Eta[9] = {
   -0.77459666924148337704, -0.77459666924148337704, -0.77459666924148337704,
    0, 0, 0,
    0.77459666924148337704, 0.77459666924148337704, 0.77459666924148337704 };

static double *NF_D_Q_Q1_2D_t = NULL;

static double NF_D_Q_Q1_2D_Weight0[9] = {
    0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827,
    0.12345679012345679012, 0.1975308641975308642, 0.12345679012345679012,
    0.077160493827160493827, 0.12345679012345679012, 0.077160493827160493827 };

static double NF_D_Q_Q1_2D_Weight1[9] = {
   -0.17930478454663967061, 0, 0.17930478454663967061,
   -0.28688765527462347298, 0, 0.28688765527462347298,
   -0.17930478454663967061, 0, 0.17930478454663967061 };

static double NF_D_Q_Q1_2D_Weight2[9] = {
   -0.17930478454663967061, -0.28688765527462347298, -0.17930478454663967061,
    0, 0, 0,
    0.17930478454663967061, 0.28688765527462347298, 0.17930478454663967061 };

static double NF_D_Q_Q1_2D_Weight3[9] = {
    0.41666666666666666667, 0, -0.41666666666666666667,
    0, 0, 0,
   -0.41666666666666666667, 0, 0.41666666666666666667 };

void NF_D_Q_Q1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] =  NF_D_Q_Q1_2D_Weight0[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight0[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight0[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight0[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight0[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight0[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight0[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight0[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight0[8]*PointValues[8];

  Functionals[1] =  NF_D_Q_Q1_2D_Weight1[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight1[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight1[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight1[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight1[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight1[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight1[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight1[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight1[8]*PointValues[8];

  Functionals[2] =  NF_D_Q_Q1_2D_Weight2[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight2[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight2[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight2[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight2[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight2[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight2[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight2[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight2[8]*PointValues[8];

  Functionals[3] =  NF_D_Q_Q1_2D_Weight3[0]*PointValues[0]
                   +NF_D_Q_Q1_2D_Weight3[1]*PointValues[1]
                   +NF_D_Q_Q1_2D_Weight3[2]*PointValues[2]
                   +NF_D_Q_Q1_2D_Weight3[3]*PointValues[3]
                   +NF_D_Q_Q1_2D_Weight3[4]*PointValues[4]
                   +NF_D_Q_Q1_2D_Weight3[5]*PointValues[5]
                   +NF_D_Q_Q1_2D_Weight3[6]*PointValues[6]
                   +NF_D_Q_Q1_2D_Weight3[7]*PointValues[7]
                   +NF_D_Q_Q1_2D_Weight3[8]*PointValues[8];
}

void NF_D_Q_Q1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_D_Q_Q1_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_Q1_2D, 4, 0, 9, 0, NF_D_Q_Q1_2D_Xi, NF_D_Q_Q1_2D_Eta,
         NF_D_Q_Q1_2D_t, NF_D_Q_Q1_2D_EvalAll, NULL);
