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
// %W% %G%
//
// Class:       TFE3DMapper1Reg
// Superclass:  TFE3DMapper
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//              for 1 regular case
//
// Author:      Gunar Matthies  17.07.2000
//
// History:     start of reimplementation 17.07.2000 (GM)
//
// =======================================================================

#include <FE3DMapper1Reg.h>
#include <FEDatabase3D.h>
#include <HNDesc.h>
#include <HangingNode.h>
#include <MooNMD_Io.h>

/** constructor, filling all data */
TFE3DMapper1Reg::TFE3DMapper1Reg(char *name, char *description,
              int nfine, int ncoarse,
              int n_pairs, int **pairs, 
              int n_noopposite, int **noopposite,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_nodes, int **twistpermutation)
 : TFE3DMapper(name, description, nfine, ncoarse, n_pairs, pairs, 
               n_noopposite, noopposite, 
               n_nodes)
{
  TwistPermutation = twistpermutation;
  CurrentPairs = pairs[0];
  
  if(n_noopposite)
    CurrentNoOp = NoOpposite[0];
  else
    CurrentNoOp = NULL;

  N_Hanging = n_hanging;
  Hanging = hanging;
  HangingTypes = hangingtypes;
  Coupling = coupling;

}

/** map the given local degrees of freedom,
    coarse cell has lower number,
    F0, F1, F2, F3 on finer side, C on coarser side */
void TFE3DMapper1Reg::Map(int *Global,
             int I_KF0, int I_KF1, int I_KF2, int I_KF3,
             int I_KC,
             int *IndicesF0, int *IndicesF1, int *IndicesF2,
             int *IndicesF3, int *IndicesC,
             int TwistIndexF0, int TwistIndexF1, int TwistIndexF2,
             int TwistIndexF3, int TwistIndexC,
             int DirichletBound,
             int &Counter,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers)
{
  static int i, j, v, w, N_;
  static THangingNode *hn;
  int *CurrentTwistPerm;

  // fill in cell local dofs according to twist index

  j=0;
  // cout << "Twist0: " << TwistIndexF0 << endl;
  CurrentTwistPerm = TwistPermutation[TwistIndexF0];
  for(i=0;i<N_DOF0;i++)
  {
    Aux[j]=I_KF0+IndicesF0[CurrentTwistPerm[i]];
    // cout << IndicesF0[i] << " - " << IndicesF0[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "Twist1: " << TwistIndexF1 << endl;
  CurrentTwistPerm = TwistPermutation[TwistIndexF1];
  for(i=0;i<N_DOF0;i++)
  {
    Aux[j]=I_KF1+IndicesF1[CurrentTwistPerm[i]];
    // cout << IndicesF1[i] << " - " << IndicesF1[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "Twist2: " << TwistIndexF2 << endl;
  CurrentTwistPerm = TwistPermutation[TwistIndexF2];
  for(i=0;i<N_DOF0;i++)
  {
    Aux[j]=I_KF2+IndicesF2[CurrentTwistPerm[i]];
    // cout << IndicesF2[i] << " - " << IndicesF2[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "Twist3: " << TwistIndexF3 << endl;
  CurrentTwistPerm = TwistPermutation[TwistIndexF3];
  for(i=0;i<N_DOF0;i++)
  {
    Aux[j]=I_KF3+IndicesF3[CurrentTwistPerm[i]];
    // cout << IndicesF3[i] << " - " << IndicesF3[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  // cout << "MapType: " << TwistIndexC << endl;
  CurrentTwistPerm = TwistPermutation[TwistIndexC];
  for(i=0;i<N_DOF1;i++)
  {
    Aux[j]=I_KC+IndicesC[CurrentTwistPerm[i]];
    // cout << IndicesC[i] << " - " << IndicesC[CurrentTwistPerm[i]];
    // cout << " - " << Aux[j] << endl;
    j++;
  }

  for(i=0;i<N_Pairs;i++)
    if(CurrentPairs[i*2]!=-1)
    {
      // cout << "local: " << CurrentPairs[i*2] << "  ";
      // cout << CurrentPairs[i*2+1] << endl;
      // cout << "global: " << Aux[CurrentPairs[i*2]] << "  ";
      // cout << Aux[CurrentPairs[i*2+1]] << endl;
      MapDOF(Global, Aux[CurrentPairs[i*2]],
             Aux[CurrentPairs[i*2+1]], Counter, vect, numbers);
    }

  for(i=0;i<N_Hanging;i++)
  {
    w=Aux[Hanging[i]];
    while( (v=Global[w]) > -1)
    {
      w=v;
    }

    // cout << i << " type: " << HangingTypes[i] << endl;
    // cout << i << " Global[w]: " << Global[w] << " " << w << endl;

    // create a hanging node only once
    if(Global[w] == HANGINGNODE) continue;

    if(Global[w] != -1)
      if(Global[w] > DirichletBound) continue;

    Global[w]=HANGINGNODE;

    N_=TFEDatabase3D::GetHNDesc3D(HangingTypes[i])->GetN_Nodes();
    hn=new THangingNode(HangingTypes[i], N_, Aux, Coupling[i]);
    vect->AddElement(hn);
    numbers->AddElement(Aux[Hanging[i]]);
    // cout << "HN: " << Aux[Hanging[i]] << endl;
  }

  for(i=0;i<N_NoOpposite;i++)
  {
    w=Aux[CurrentNoOp[i]];
    while( (v=Global[w]) > -1 )
    {
      w=v;
    }

    if( Global[w] == -1 )
    {
      Counter--;
      Global[Aux[CurrentNoOp[i]]]=Counter;
    }
  } // endfor i

}
