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

static double NF_C_Q_Q3_2D_Xi[] = { -1, -0.33333333333333333333,
                                  0.33333333333333333333, 1,
                                  -1, -0.33333333333333333333,
                                  0.33333333333333333333, 1,
                                  -1, -0.33333333333333333333,
                                  0.33333333333333333333, 1,
                                  -1, -0.33333333333333333333,
                                  0.33333333333333333333, 1 };
static double NF_C_Q_Q3_2D_Eta[] = { -1, -1, -1, -1,
                                   -0.33333333333333333333,
                                   -0.33333333333333333333,
                                   -0.33333333333333333333,
                                   -0.33333333333333333333,
                                   0.33333333333333333333,
                                   0.33333333333333333333,
                                   0.33333333333333333333,
                                   0.33333333333333333333,
                                   1, 1, 1, 1 };
static double NF_C_Q_Q3_2D_T[] = { -1, -0.33333333333333333333,
                                 0.33333333333333333333, 1 };

void NF_C_Q_Q3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
/*
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
  Functionals[4] = PointValues[4];
  Functionals[5] = PointValues[5];
  Functionals[6] = PointValues[6];
  Functionals[7] = PointValues[7];
  Functionals[8] = PointValues[8];
  Functionals[9] = PointValues[9];
  Functionals[10] = PointValues[10];
  Functionals[11] = PointValues[11];
  Functionals[12] = PointValues[12];
  Functionals[13] = PointValues[13];
  Functionals[14] = PointValues[14];
  Functionals[15] = PointValues[15];
*/
  memcpy(Functionals, PointValues, 16*SizeOfDouble);
}

void NF_C_Q_Q3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
                           double *Functionals)
{
/*
  Functionals[0] = PointValues[0];
  Functionals[1] = PointValues[1];
  Functionals[2] = PointValues[2];
  Functionals[3] = PointValues[3];
*/
  memcpy(Functionals, PointValues, 4*SizeOfDouble);
}

TNodalFunctional2D *NF_C_Q_Q3_2D_Obj = new TNodalFunctional2D
        (NF_C_Q_Q3_2D, 16, 4, 16, 4, NF_C_Q_Q3_2D_Xi, NF_C_Q_Q3_2D_Eta,
         NF_C_Q_Q3_2D_T, NF_C_Q_Q3_2D_EvalAll, NF_C_Q_Q3_2D_EvalEdge);
