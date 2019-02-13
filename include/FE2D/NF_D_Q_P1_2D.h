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
   
static double NF_D_Q_P1_2D_Xi[9]={ 0.774596669241483, -0.774596669241483,
                                 0.774596669241483, -0.774596669241483,
                                 0.774596669241483, -0.774596669241483,
                                 0,                  0,
                                 0 };
static double NF_D_Q_P1_2D_Eta[9]={  0.774596669241483,  0.774596669241483,
                                  -0.774596669241483, -0.774596669241483,
                                   0,                  0,
                                   0.774596669241483, -0.774596669241483,
                                   0 };

static double *NF_D_Q_P1_2D_t = NULL;

void NF_D_Q_P1_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  static double weights[9]={ 0.308641975308642, 0.308641975308642,
                             0.308641975308642, 0.308641975308642,
                             0.493827160493827, 0.493827160493827,
                             0.493827160493827, 0.493827160493827,
                             0.790123456790123 };
  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2]
                   +weights[3]*PointValues[3]
                   +weights[4]*PointValues[4]
                   +weights[5]*PointValues[5]
                   +weights[6]*PointValues[6]
                   +weights[7]*PointValues[7]
                   +weights[8]*PointValues[8];
  Functionals[1] =  weights[0]*PointValues[0]*NF_D_Q_P1_2D_Xi[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P1_2D_Xi[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P1_2D_Xi[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P1_2D_Xi[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P1_2D_Xi[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P1_2D_Xi[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P1_2D_Xi[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P1_2D_Xi[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P1_2D_Xi[8];
  Functionals[2] =  weights[0]*PointValues[0]*NF_D_Q_P1_2D_Eta[0]
                   +weights[1]*PointValues[1]*NF_D_Q_P1_2D_Eta[1]
                   +weights[2]*PointValues[2]*NF_D_Q_P1_2D_Eta[2]
                   +weights[3]*PointValues[3]*NF_D_Q_P1_2D_Eta[3]
                   +weights[4]*PointValues[4]*NF_D_Q_P1_2D_Eta[4]
                   +weights[5]*PointValues[5]*NF_D_Q_P1_2D_Eta[5]
                   +weights[6]*PointValues[6]*NF_D_Q_P1_2D_Eta[6]
                   +weights[7]*PointValues[7]*NF_D_Q_P1_2D_Eta[7]
                   +weights[8]*PointValues[8]*NF_D_Q_P1_2D_Eta[8];

  Functionals[0] *= 0.25;
  Functionals[1] *= 0.25;
  Functionals[2] *= 0.25;
}

void NF_D_Q_P1_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint, double *PointValues,
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

TNodalFunctional2D *NF_D_Q_P1_2D_Obj = new TNodalFunctional2D
        (NF_D_Q_P1_2D, 3, 0, 9, 0, NF_D_Q_P1_2D_Xi, NF_D_Q_P1_2D_Eta,
         NF_D_Q_P1_2D_t, NF_D_Q_P1_2D_EvalAll, NULL);
