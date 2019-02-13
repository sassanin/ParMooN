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
   

static double NF_N_Q_Q3_2D_Xi[32] = 
        {
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
           1, 1, 1, 1,
           0.86113631159405257522394649,
           0.33998104358485626480266576,
          -0.33998104358485626480266576,
          -0.86113631159405257522394649,
          -1, -1, -1, -1,

          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649
        };

static double NF_N_Q_Q3_2D_Eta[32] = 
        { -1, -1, -1, -1,
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
           1, 1, 1, 1,
           0.86113631159405257522394649,
           0.33998104358485626480266576,
          -0.33998104358485626480266576,
          -0.86113631159405257522394649,

          -0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.86113631159405257522394649,
          -0.33998104358485626480266576,
          -0.33998104358485626480266576,
          -0.33998104358485626480266576,
          -0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.33998104358485626480266576,
           0.86113631159405257522394649,
           0.86113631159405257522394649,
           0.86113631159405257522394649,
           0.86113631159405257522394649,
        };

static double NF_N_Q_Q3_2D_T[4] = 
        {
               -0.86113631159405257522394649,
               -0.33998104358485626480266576,
                0.33998104358485626480266576,
                0.86113631159405257522394649
        };

static double NF_N_Q_Q3_2D_EdgeWeight0[4] = {
           0.1739274225687269286865319745,
           0.3260725774312730713134680255,
           0.3260725774312730713134680255,
           0.1739274225687269286865319745 };

static double NF_N_Q_Q3_2D_EdgeWeight1[4] = {
          -0.4493256574676810538434880055,
          -0.3325754854784642078819958520,
           0.3325754854784642078819958520,
           0.4493256574676810538434880055 };

static double NF_N_Q_Q3_2D_EdgeWeight2[4] = {
          0.532508042018911499194276175,
         -0.532508042018911499194276175,
         -0.532508042018911499194276175,
          0.532508042018911499194276175 };

static double NF_N_Q_Q3_CellWeight0[16] = {
        .0302507483214005013803030243039, .0567129629629629629629629629461,
        .0567129629629629629629629629461, .0302507483214005013803030243039,
        .0567129629629629629629629629461, .106323325752673572693771049804,
        .106323325752673572693771049804, .0567129629629629629629629629461,
        .0567129629629629629629629629461, .106323325752673572693771049804,
        .106323325752673572693771049804, .0567129629629629629629629629461,
        .0302507483214005013803030243039, .0567129629629629629629629629461,
        .0567129629629629629629629629461, .0302507483214005013803030243039 };

static double NF_N_Q_Q3_CellWeight1[16] = {
       -.0781500534973524151648919963067, -.0578439969988123506087868827606,
        .0578439969988123506087868827606, .0781500534973524151648919963067,
       -.146512775236488111756852006936, -.108443745740419753332211043771,
        .108443745740419753332211043771, .146512775236488111756852006936,
       -.146512775236488111756852006936, -.108443745740419753332211043771,
        .108443745740419753332211043771, .146512775236488111756852006936,
       -.0781500534973524151648919963067, -.0578439969988123506087868827606,
        .0578439969988123506087868827606, .0781500534973524151648919963067 };

static double NF_N_Q_Q3_CellWeight2[16] = {
       -.0781500534973524151648919963067, -.146512775236488111756852006936,
       -.146512775236488111756852006936, -.0781500534973524151648919963067,
       -.0578439969988123506087868827606, -.108443745740419753332211043771,
       -.108443745740419753332211043771, -.0578439969988123506087868827606,
        .0578439969988123506087868827606, .108443745740419753332211043771,
        .108443745740419753332211043771, .0578439969988123506087868827606,
        .0781500534973524151648919963067, .146512775236488111756852006936,
        .146512775236488111756852006936, .0781500534973524151648919963067 };

void NF_N_Q_Q3_2D_EvalAll(TCollection *Coll, TBaseCell *Cell, double *PointValues,
                          double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_Q_Q3_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight0[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight0[3]*PointValues[3]);
  Functionals[1] = ( NF_N_Q_Q3_2D_EdgeWeight0[0]*PointValues[4]
                    +NF_N_Q_Q3_2D_EdgeWeight0[1]*PointValues[5]
                    +NF_N_Q_Q3_2D_EdgeWeight0[2]*PointValues[6]
                    +NF_N_Q_Q3_2D_EdgeWeight0[3]*PointValues[7]);
  Functionals[2] = ( NF_N_Q_Q3_2D_EdgeWeight0[0]*PointValues[8]
                    +NF_N_Q_Q3_2D_EdgeWeight0[1]*PointValues[9]
                    +NF_N_Q_Q3_2D_EdgeWeight0[2]*PointValues[10]
                    +NF_N_Q_Q3_2D_EdgeWeight0[3]*PointValues[11]);
  Functionals[3] = ( NF_N_Q_Q3_2D_EdgeWeight0[0]*PointValues[12]
                    +NF_N_Q_Q3_2D_EdgeWeight0[1]*PointValues[13]
                    +NF_N_Q_Q3_2D_EdgeWeight0[2]*PointValues[14]
                    +NF_N_Q_Q3_2D_EdgeWeight0[3]*PointValues[15]);

  Functionals[4] = ( NF_N_Q_Q3_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight1[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight1[3]*PointValues[3]);
  Functionals[5] = ( NF_N_Q_Q3_2D_EdgeWeight1[0]*PointValues[4]
                    +NF_N_Q_Q3_2D_EdgeWeight1[1]*PointValues[5]
                    +NF_N_Q_Q3_2D_EdgeWeight1[2]*PointValues[6]
                    +NF_N_Q_Q3_2D_EdgeWeight1[3]*PointValues[7]);
  Functionals[6] = ( NF_N_Q_Q3_2D_EdgeWeight1[0]*PointValues[8]
                    +NF_N_Q_Q3_2D_EdgeWeight1[1]*PointValues[9]
                    +NF_N_Q_Q3_2D_EdgeWeight1[2]*PointValues[10]
                    +NF_N_Q_Q3_2D_EdgeWeight1[3]*PointValues[11]);
  Functionals[7] = ( NF_N_Q_Q3_2D_EdgeWeight1[0]*PointValues[12]
                    +NF_N_Q_Q3_2D_EdgeWeight1[1]*PointValues[13]
                    +NF_N_Q_Q3_2D_EdgeWeight1[2]*PointValues[14]
                    +NF_N_Q_Q3_2D_EdgeWeight1[3]*PointValues[15]);

  Functionals[8] = ( NF_N_Q_Q3_2D_EdgeWeight2[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight2[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight2[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight2[3]*PointValues[3]);
  Functionals[9] = ( NF_N_Q_Q3_2D_EdgeWeight2[0]*PointValues[4]
                    +NF_N_Q_Q3_2D_EdgeWeight2[1]*PointValues[5]
                    +NF_N_Q_Q3_2D_EdgeWeight2[2]*PointValues[6]
                    +NF_N_Q_Q3_2D_EdgeWeight2[3]*PointValues[7]);
  Functionals[10]= ( NF_N_Q_Q3_2D_EdgeWeight2[0]*PointValues[8]
                    +NF_N_Q_Q3_2D_EdgeWeight2[1]*PointValues[9]
                    +NF_N_Q_Q3_2D_EdgeWeight2[2]*PointValues[10]
                    +NF_N_Q_Q3_2D_EdgeWeight2[3]*PointValues[11]);
  Functionals[11]= ( NF_N_Q_Q3_2D_EdgeWeight2[0]*PointValues[12]
                    +NF_N_Q_Q3_2D_EdgeWeight2[1]*PointValues[13]
                    +NF_N_Q_Q3_2D_EdgeWeight2[2]*PointValues[14]
                    +NF_N_Q_Q3_2D_EdgeWeight2[3]*PointValues[15]);

  Functionals[12] =( NF_N_Q_Q3_CellWeight0[ 0]*PointValues[16]
                    +NF_N_Q_Q3_CellWeight0[ 1]*PointValues[17]
                    +NF_N_Q_Q3_CellWeight0[ 2]*PointValues[18]
                    +NF_N_Q_Q3_CellWeight0[ 3]*PointValues[19]
                    +NF_N_Q_Q3_CellWeight0[ 4]*PointValues[20]
                    +NF_N_Q_Q3_CellWeight0[ 5]*PointValues[21]
                    +NF_N_Q_Q3_CellWeight0[ 6]*PointValues[22]
                    +NF_N_Q_Q3_CellWeight0[ 7]*PointValues[23]
                    +NF_N_Q_Q3_CellWeight0[ 8]*PointValues[24]
                    +NF_N_Q_Q3_CellWeight0[ 9]*PointValues[25]
                    +NF_N_Q_Q3_CellWeight0[10]*PointValues[26]
                    +NF_N_Q_Q3_CellWeight0[11]*PointValues[27]
                    +NF_N_Q_Q3_CellWeight0[12]*PointValues[28]
                    +NF_N_Q_Q3_CellWeight0[13]*PointValues[29]
                    +NF_N_Q_Q3_CellWeight0[14]*PointValues[30]
                    +NF_N_Q_Q3_CellWeight0[15]*PointValues[31]);
  Functionals[13] =( NF_N_Q_Q3_CellWeight1[ 0]*PointValues[16]
                    +NF_N_Q_Q3_CellWeight1[ 1]*PointValues[17]
                    +NF_N_Q_Q3_CellWeight1[ 2]*PointValues[18]
                    +NF_N_Q_Q3_CellWeight1[ 3]*PointValues[19]
                    +NF_N_Q_Q3_CellWeight1[ 4]*PointValues[20]
                    +NF_N_Q_Q3_CellWeight1[ 5]*PointValues[21]
                    +NF_N_Q_Q3_CellWeight1[ 6]*PointValues[22]
                    +NF_N_Q_Q3_CellWeight1[ 7]*PointValues[23]
                    +NF_N_Q_Q3_CellWeight1[ 8]*PointValues[24]
                    +NF_N_Q_Q3_CellWeight1[ 9]*PointValues[25]
                    +NF_N_Q_Q3_CellWeight1[10]*PointValues[26]
                    +NF_N_Q_Q3_CellWeight1[11]*PointValues[27]
                    +NF_N_Q_Q3_CellWeight1[12]*PointValues[28]
                    +NF_N_Q_Q3_CellWeight1[13]*PointValues[29]
                    +NF_N_Q_Q3_CellWeight1[14]*PointValues[30]
                    +NF_N_Q_Q3_CellWeight1[15]*PointValues[31]);
  Functionals[14] =( NF_N_Q_Q3_CellWeight2[ 0]*PointValues[16]
                    +NF_N_Q_Q3_CellWeight2[ 1]*PointValues[17]
                    +NF_N_Q_Q3_CellWeight2[ 2]*PointValues[18]
                    +NF_N_Q_Q3_CellWeight2[ 3]*PointValues[19]
                    +NF_N_Q_Q3_CellWeight2[ 4]*PointValues[20]
                    +NF_N_Q_Q3_CellWeight2[ 5]*PointValues[21]
                    +NF_N_Q_Q3_CellWeight2[ 6]*PointValues[22]
                    +NF_N_Q_Q3_CellWeight2[ 7]*PointValues[23]
                    +NF_N_Q_Q3_CellWeight2[ 8]*PointValues[24]
                    +NF_N_Q_Q3_CellWeight2[ 9]*PointValues[25]
                    +NF_N_Q_Q3_CellWeight2[10]*PointValues[26]
                    +NF_N_Q_Q3_CellWeight2[11]*PointValues[27]
                    +NF_N_Q_Q3_CellWeight2[12]*PointValues[28]
                    +NF_N_Q_Q3_CellWeight2[13]*PointValues[29]
                    +NF_N_Q_Q3_CellWeight2[14]*PointValues[30]
                    +NF_N_Q_Q3_CellWeight2[15]*PointValues[31]);

  /*
  if(Cell)
  {
    if(Cell->GetVertex(0) > Cell->GetVertex(1))
      Functionals[4] = -Functionals[4];
    if(Cell->GetVertex(1) > Cell->GetVertex(2))
      Functionals[5] = -Functionals[5];
    if(Cell->GetVertex(2) > Cell->GetVertex(3))
      Functionals[6] = -Functionals[6];
    if(Cell->GetVertex(3) > Cell->GetVertex(0))
      Functionals[7] = -Functionals[7];
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
        Functionals[ 4] = -Functionals[ 4];
    } // endif neigh

    neigh = Cell->GetJoint(1)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 5] = -Functionals[ 5];
    } // endif neigh

    neigh = Cell->GetJoint(2)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 6] = -Functionals[ 6];
    } // endif neigh

    neigh = Cell->GetJoint(3)->GetNeighbour(Cell);
    if(neigh)
    {
      NeighNum = Coll->GetIndex(neigh);
      if(NeighNum < OwnNum)
        Functionals[ 7] = -Functionals[ 7];
    } // endif neigh
  } // endif Cell
}

void NF_N_Q_Q3_2D_EvalEdge(TCollection *Coll, TBaseCell *Cell, int Joint,
                           double *PointValues, double *Functionals)
{
  int OwnNum, NeighNum;
  TBaseCell *neigh;

  Functionals[0] = ( NF_N_Q_Q3_2D_EdgeWeight0[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight0[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight0[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight0[3]*PointValues[3]);

  Functionals[1] = ( NF_N_Q_Q3_2D_EdgeWeight1[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight1[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight1[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight1[3]*PointValues[3]);

  Functionals[2] = ( NF_N_Q_Q3_2D_EdgeWeight2[0]*PointValues[0]
                    +NF_N_Q_Q3_2D_EdgeWeight2[1]*PointValues[1]
                    +NF_N_Q_Q3_2D_EdgeWeight2[2]*PointValues[2]
                    +NF_N_Q_Q3_2D_EdgeWeight2[3]*PointValues[3]);

  if(Joint != -1)
  {
    // if(Cell->GetVertex(Joint) > Cell->GetVertex((Joint+1)%4))
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

TNodalFunctional2D *NF_N_Q_Q3_2D_Obj = new TNodalFunctional2D
        (NF_N_Q_Q3_2D, 15, 3, 32, 4, NF_N_Q_Q3_2D_Xi, NF_N_Q_Q3_2D_Eta,
         NF_N_Q_Q3_2D_T, NF_N_Q_Q3_2D_EvalAll, NF_N_Q_Q3_2D_EvalEdge);
