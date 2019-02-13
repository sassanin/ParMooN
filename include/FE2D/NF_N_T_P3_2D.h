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
   
static double NF_N_T_P3_2D_Xi[19] = 
        { 
          0.069431844202973712388026755553,
          0.330009478207571867598667120449, 
          0.669990521792428132401332879551,
          0.930568155797026287611973244447, 
          0.930568155797026287611973244447, 
          0.669990521792428132401332879551,
          0.330009478207571867598667120449, 
          0.069431844202973712388026755553,
          0, 0, 0, 0,

          0.3333333333333333, 
          0.7974269853530873, 0.1012865073234563, 
          0.1012865073234563, 
          0.05971587178976982, 0.4701420641051151, 
          0.4701420641051151
        };

static double NF_N_T_P3_2D_Eta[19] = 
        { 
          0, 0, 0, 0,
          0.069431844202973712388026755553,
          0.330009478207571867598667120449, 
          0.669990521792428132401332879551,
          0.930568155797026287611973244447, 
          0.930568155797026287611973244447, 
          0.669990521792428132401332879551,
          0.330009478207571867598667120449, 
          0.069431844202973712388026755553,

          0.3333333333333333, 
          0.1012865073234563, 0.7974269853530873, 
          0.1012865073234563, 
          0.4701420641051151, 0.4701420641051151, 
          0.05971587178976982
        };

static double NF_N_T_P3_2D_T[4] = 
        { 
          -0.861136311594052575223946488893, 
          -0.339981043584856264802665759103, 
           0.339981043584856264802665759103,
           0.861136311594052575223946488893
        };

static double NF_N_T_P3_2D_EdgeWeight0[4] = {
           0.1739274225687269286865319745,
           0.3260725774312730713134680255,
           0.3260725774312730713134680255,
           0.1739274225687269286865319745 };

static double NF_N_T_P3_2D_EdgeWeight1[4] = {
          -0.4493256574676810538434880055,
          -0.3325754854784642078819958520,
           0.3325754854784642078819958520,
           0.4493256574676810538434880055 };

static double NF_N_T_P3_2D_EdgeWeight2[4] = {
          0.532508042018911499194276175,
         -0.532508042018911499194276175,
         -0.532508042018911499194276175,
          0.532508042018911499194276175 };

static double NF_N_T_P3_2D_CellWeight0[7] = {
          0.225,
          0.0382678191976911531218063401284,
          0.0382678191976911531218063401284,
          0.3012819032390987693198956466394,
          0.18673218080230856483484009381144,
          0.0237180967609008544430214777564,
          0.18673218080230856483484009381144,
 };

static double NF_N_T_P3_2D_CellWeight1[7] = {
          0.225,
          0.3012819032390987315381414831913,
          0.0382678191976911153400521766803,
          0.0382678191976911153400521766803,
          0.02371809676090086238667064506676,
          0.1867321808023085727784892611218,
          0.1867321808023085727784892611218 };


static double NF_N_T_P3_2D_CellWeight2[7] = {
          0.225,
          0.0382678191976911153400521766803,
          0.3012819032390987315381414831913,
          0.0382678191976911153400521766803,
          0.1867321808023085727784892611218,
          0.1867321808023085727784892611218,
          0.02371809676090086238667064506676 };

void NF_N_T_P3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_T_P3_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight0[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight0[3]*PointValues[3]);
  Functionals[1] = ( NF_N_T_P3_2D_EdgeWeight0[0]*PointValues[4]
                    +NF_N_T_P3_2D_EdgeWeight0[1]*PointValues[5]
                    +NF_N_T_P3_2D_EdgeWeight0[2]*PointValues[6]
                    +NF_N_T_P3_2D_EdgeWeight0[3]*PointValues[7]);
  Functionals[2] = ( NF_N_T_P3_2D_EdgeWeight0[0]*PointValues[8]
                    +NF_N_T_P3_2D_EdgeWeight0[1]*PointValues[9]
                    +NF_N_T_P3_2D_EdgeWeight0[2]*PointValues[10]
                    +NF_N_T_P3_2D_EdgeWeight0[3]*PointValues[11]);

  Functionals[3] = ( NF_N_T_P3_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight1[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight1[3]*PointValues[3]);
  Functionals[4] = ( NF_N_T_P3_2D_EdgeWeight1[0]*PointValues[4]
                    +NF_N_T_P3_2D_EdgeWeight1[1]*PointValues[5]
                    +NF_N_T_P3_2D_EdgeWeight1[2]*PointValues[6]
                    +NF_N_T_P3_2D_EdgeWeight1[3]*PointValues[7]);
  Functionals[5] = ( NF_N_T_P3_2D_EdgeWeight1[0]*PointValues[8]
                    +NF_N_T_P3_2D_EdgeWeight1[1]*PointValues[9]
                    +NF_N_T_P3_2D_EdgeWeight1[2]*PointValues[10]
                    +NF_N_T_P3_2D_EdgeWeight1[3]*PointValues[11]);

  Functionals[6] = ( NF_N_T_P3_2D_EdgeWeight2[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight2[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight2[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight2[3]*PointValues[3]);
  Functionals[7] = ( NF_N_T_P3_2D_EdgeWeight2[0]*PointValues[4]
                    +NF_N_T_P3_2D_EdgeWeight2[1]*PointValues[5]
                    +NF_N_T_P3_2D_EdgeWeight2[2]*PointValues[6]
                    +NF_N_T_P3_2D_EdgeWeight2[3]*PointValues[7]);
  Functionals[8] = ( NF_N_T_P3_2D_EdgeWeight2[0]*PointValues[8]
                    +NF_N_T_P3_2D_EdgeWeight2[1]*PointValues[9]
                    +NF_N_T_P3_2D_EdgeWeight2[2]*PointValues[10]
                    +NF_N_T_P3_2D_EdgeWeight2[3]*PointValues[11]);

  Functionals[9] =( NF_N_T_P3_2D_CellWeight0[0]*PointValues[12]
                   +NF_N_T_P3_2D_CellWeight0[1]*PointValues[13]
                   +NF_N_T_P3_2D_CellWeight0[2]*PointValues[14]
                   +NF_N_T_P3_2D_CellWeight0[3]*PointValues[15]
                   +NF_N_T_P3_2D_CellWeight0[4]*PointValues[16]
                   +NF_N_T_P3_2D_CellWeight0[5]*PointValues[17]
                   +NF_N_T_P3_2D_CellWeight0[6]*PointValues[18] );
  Functionals[10]=( NF_N_T_P3_2D_CellWeight1[0]*PointValues[12]
                   +NF_N_T_P3_2D_CellWeight1[1]*PointValues[13]
                   +NF_N_T_P3_2D_CellWeight1[2]*PointValues[14]
                   +NF_N_T_P3_2D_CellWeight1[3]*PointValues[15]
                   +NF_N_T_P3_2D_CellWeight1[4]*PointValues[16]
                   +NF_N_T_P3_2D_CellWeight1[5]*PointValues[17]
                   +NF_N_T_P3_2D_CellWeight1[6]*PointValues[18] );
  Functionals[11]=( NF_N_T_P3_2D_CellWeight2[0]*PointValues[12]
                   +NF_N_T_P3_2D_CellWeight2[1]*PointValues[13]
                   +NF_N_T_P3_2D_CellWeight2[2]*PointValues[14]
                   +NF_N_T_P3_2D_CellWeight2[3]*PointValues[15]
                   +NF_N_T_P3_2D_CellWeight2[4]*PointValues[16]
                   +NF_N_T_P3_2D_CellWeight2[5]*PointValues[17]
                   +NF_N_T_P3_2D_CellWeight2[6]*PointValues[18] );

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[3] = -Functionals[3];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[4] = -Functionals[4];
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
        Functionals[3] = -Functionals[3];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 4] = -Functionals[ 4];
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

void NF_N_T_P3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                           double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_T_P3_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight0[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight0[3]*PointValues[3]);
  Functionals[1] = ( NF_N_T_P3_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight1[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight1[3]*PointValues[3]);
  Functionals[2] = ( NF_N_T_P3_2D_EdgeWeight2[0]*PointValues[0]
                    +NF_N_T_P3_2D_EdgeWeight2[1]*PointValues[1]
                    +NF_N_T_P3_2D_EdgeWeight2[2]*PointValues[2]
                    +NF_N_T_P3_2D_EdgeWeight2[3]*PointValues[3]);

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
                       int n_pointsall, int n_pointsedge,
                       double *xi, double *eta, double *t,
                       DoubleFunctVect *evalall,
                       DoubleFunctVect *evaledge);
*/

TNodalFunctional2D *NF_N_T_P3_2D_Obj = new TNodalFunctional2D
        (NF_N_T_P3_2D, 12, 3, 19, 4, NF_N_T_P3_2D_Xi, NF_N_T_P3_2D_Eta,
         NF_N_T_P3_2D_T, NF_N_T_P3_2D_EvalAll, NF_N_T_P3_2D_EvalEdge);
