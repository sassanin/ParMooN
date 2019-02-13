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
   
static double NF_N_T_P1MOD_2D_Xi[] = 
        { 0.11270166537925831149, 0.5, 
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149,
          0, 0, 0 };

static double NF_N_T_P1MOD_2D_Eta[] = 
        { 0, 0, 0,
          0.11270166537925831149, 0.5,
          0.88729833462074168851,
          0.88729833462074168851, 0.5,
          0.11270166537925831149 };

static double NF_N_T_P1MOD_2D_T[] = 
        { -0.77459666924148337703585307995647992, 0,
           0.77459666924148337703585307995647992 };

void NF_N_T_P1MOD_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                             double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  static double weights[]={ 0.27777777777777777777777777778,
                            0.44444444444444444444444444444,
                            0.27777777777777777777777777778 };

  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2];
  Functionals[2] =  weights[0]*PointValues[3]
                   +weights[1]*PointValues[4]
                   +weights[2]*PointValues[5];
  Functionals[4] =  weights[0]*PointValues[6]
                   +weights[1]*PointValues[7]
                   +weights[2]*PointValues[8];

  Functionals[1] = 60*( weights[0]*(NF_N_T_P1MOD_2D_Xi[0]-0.5)*PointValues[0]
                       +weights[1]*(NF_N_T_P1MOD_2D_Xi[1]-0.5)*PointValues[1]
                       +weights[2]*(NF_N_T_P1MOD_2D_Xi[2]-0.5)*PointValues[2]);
  Functionals[3] = 60*( weights[0]*(NF_N_T_P1MOD_2D_Eta[3]-0.5)*PointValues[3]
                       +weights[1]*(NF_N_T_P1MOD_2D_Eta[4]-0.5)*PointValues[4]
                       +weights[2]*(NF_N_T_P1MOD_2D_Eta[5]-0.5)*PointValues[5]);
  Functionals[5] = 60*( weights[0]*(0.5-NF_N_T_P1MOD_2D_Eta[6])*PointValues[6]
                       +weights[1]*(0.5-NF_N_T_P1MOD_2D_Eta[7])*PointValues[7]
                       +weights[2]*(0.5-NF_N_T_P1MOD_2D_Eta[8])*PointValues[8]);
  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[1] = -Functionals[1];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[3] = -Functionals[3];
    if(Cell->GetVertex(2) > Cell->GetVertex(0))
      Functionals[5] = -Functionals[5];
  }
  */

  if(Cell)
  {
    OwnNum = Coll->GetIndex(Cell);

    neigh = Cell->GetJoint(0)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[1] = -Functionals[1];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[3] = -Functionals[3];
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 5] = -Functionals[ 5];
    } // endif neigh
  } // endif Cell
}

void NF_N_T_P1MOD_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                              double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  static double weights[3]={ 0.27777777777777777777777777778,
                            0.44444444444444444444444444444,
                            0.27777777777777777777777777778 };

  Functionals[0] =  weights[0]*PointValues[0]
                   +weights[1]*PointValues[1]
                   +weights[2]*PointValues[2];

  Functionals[1] = 60*( weights[0]*NF_N_T_P1MOD_2D_T[0]*PointValues[0]
                       +weights[1]*NF_N_T_P1MOD_2D_T[1]*PointValues[1]
                       +weights[2]*NF_N_T_P1MOD_2D_T[2]*PointValues[2])*0.5;

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%3))
    //   Functionals[1] = -Functionals[1];
    // /*
    neigh = Cell->GetJoint(Joint)->GetNeighbour(Cell);
    if(neigh)
    {
      OwnNum = Coll->GetIndex(Cell);
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[1] = -Functionals[1];
    } // endif neigh
    // */
  }
}

/*
    TNodalFunctional2D(NodalFunctional2D id,
                       int n_allfunctionals, int n_edgefunctionals,
                       int B_T_Pointsall, int B_T_Pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_N_T_P1MOD_2D_Obj = new TNodalFunctional2D
        (NF_N_T_P1MOD_2D, 6, 2, 9, 3, NF_N_T_P1MOD_2D_Xi, NF_N_T_P1MOD_2D_Eta,
         NF_N_T_P1MOD_2D_T, NF_N_T_P1MOD_2D_EvalAll, NF_N_T_P1MOD_2D_EvalEdge);
