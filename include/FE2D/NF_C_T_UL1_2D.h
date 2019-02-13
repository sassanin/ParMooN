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

static double NF_C_T_UL1_2D_Xi[10] = { 0.0, 1.0, 0.0,
                  0.333333333333333333333333333333333,
                  0.797426985353087322398025276169754,
                  0.101286507323456338800987361915123,
                  0.101286507323456338800987361915123,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447, 
                  0.470142064105115089770441209513447 };

static double NF_C_T_UL1_2D_Eta[] = { 0.0, 0.0, 1.0,
                  0.333333333333333333333333333333333, 
                  0.101286507323456338800987361915123,
                  0.797426985353087322398025276169754, 
                  0.101286507323456338800987361915123, 
                  0.470142064105115089770441209513447,
                  0.059715871789769820459117580973106,
                  0.470142064105115089770441209513447 }; 

static double NF_C_T_UL1_2D_W[7] =
                { 0.1125, 
                  0.0629695902724135762978419727500906,
                  0.0629695902724135762978419727500906,
                  0.0629695902724135762978419727500906,
                  0.0661970763942530903688246939165759,
                  0.0661970763942530903688246939165759,
                  0.0661970763942530903688246939165759 };

static double NF_C_T_UL1_2D_T[2] = { -1, 1 };

void NF_C_T_UL1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = 2*(  PointValues[3]*NF_C_T_UL1_2D_W[0]
                      + PointValues[4]*NF_C_T_UL1_2D_W[1]
                      + PointValues[5]*NF_C_T_UL1_2D_W[2]
                      + PointValues[6]*NF_C_T_UL1_2D_W[3]
                      + PointValues[7]*NF_C_T_UL1_2D_W[4]
                      + PointValues[8]*NF_C_T_UL1_2D_W[5]
                      + PointValues[9]*NF_C_T_UL1_2D_W[6] );
}

void NF_C_T_UL1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
}

TNodalFunctional2D *NF_C_T_UL1_2D_Obj = new TNodalFunctional2D
        (NF_C_T_UL1_2D, 4, 2, 10, 2, NF_C_T_UL1_2D_Xi, NF_C_T_UL1_2D_Eta,
         NF_C_T_UL1_2D_T, NF_C_T_UL1_2D_EvalAll, NF_C_T_UL1_2D_EvalEdge);
