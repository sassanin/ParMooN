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
*/

static double NF_D_T_SV1_2D_Xi[] = {
                                     7.0/18, 2.0/9, 13.0/18,
                                     7.0/18, 13.0/18, 2.0/9,
                                     2.0/9, 1.0/18, 1.0/18
                                   };
static double NF_D_T_SV1_2D_Eta[] = { 
                                      2.0/9, 1.0/18, 1.0/18,
                                      7.0/18, 2.0/9, 13.0/18,
                                      7.0/18, 13.0/18, 2.0/9
                                    };
static double *NF_D_T_SV1_2D_T = NULL;

void NF_D_T_SV1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = ( 2*PointValues[0]+2*PointValues[1]+2*PointValues[2])/6;
  Functionals[1] = (-4*PointValues[0]-4*PointValues[1]+8*PointValues[2])/6;
  Functionals[2] = ( 8*PointValues[0]-4*PointValues[1]-4*PointValues[2])/6;
  
  Functionals[3] = ( 2*PointValues[3]+2*PointValues[4]+2*PointValues[5])/6;
  Functionals[4] = (-4*PointValues[3]-4*PointValues[4]+8*PointValues[5])/6;
  Functionals[5] = ( 8*PointValues[3]-4*PointValues[4]-4*PointValues[5])/6;

  Functionals[6] = ( 2*PointValues[6]+2*PointValues[7]+2*PointValues[8])/6;
  Functionals[7] = (-4*PointValues[6]-4*PointValues[7]+8*PointValues[8])/6;
  Functionals[8] = ( 8*PointValues[6]-4*PointValues[7]-4*PointValues[8])/6;
}

void NF_D_T_SV1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
}

TNodalFunctional2D *NF_D_T_SV1_2D_Obj = new TNodalFunctional2D
        (NF_D_T_SV1_2D, 9, 0, 9, 0, NF_D_T_SV1_2D_Xi, NF_D_T_SV1_2D_Eta,
         NF_D_T_SV1_2D_T, NF_D_T_SV1_2D_EvalAll, NF_D_T_SV1_2D_EvalEdge);
