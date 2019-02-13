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
// @(#)FE2DMapper1Reg.C        1.1 10/30/98
//
// Class:       TFE2DMapper1Reg
// Superclass:  TFE2DMapper
// Purpose:     find out which of the given local degress of freedom
//              are equivalent to the same global degree of freedom
//              for 1 regular case
//
// Author:      Gunar Matthies  11.11.97
//
// History:     start of reimplementation 03.08.1998 (GM)
//
// =======================================================================

#include <FE2DMapper1Reg.h>
#include <FEDatabase2D.h>
#include <HNDesc.h>
#include <MooNMD_Io.h>

/** constructor, filling all data */
TFE2DMapper1Reg::TFE2DMapper1Reg(char *name, char *description, 
              int n0, int n1, int n2,
              int n_pairs, int *pairs, 
              int n_mid, int *mid,
              int n_hanging, int *hanging,
              HNDesc *hangingtypes, int **coupling,
              int n_farhanging, int *farhanging,
              HNDesc *farhangingtypes, int ****farcoupling,
              int n_noopposite, int *noopposite,
              int n_nodes)
 : TFE2DMapper(name, description, n0, n1, n_pairs, pairs, 
             n_hanging, hanging, hangingtypes, coupling, 
             n_farhanging, farhanging, farhangingtypes, farcoupling, 
             n_noopposite, noopposite, n_nodes)
{
  N_DOF2=n2;

  N_Mid=n_mid;

  Mid=mid;
}

/** map the given local degrees of freedom,
    coarse cell has lower number,
    0 is coarser side 0; 1,2 are on finer side 1 */
void TFE2DMapper1Reg::MapCoarseFine(int *Global, int I_K0, int I_K1, int I_K2,
             int *Indices0, int *Indices1, int *Indices2,
             int LocEdge0, int LocEdge1, int LocEdge2,
             TFEDesc2D *Desc0, TFEDesc2D *Desc1, TFEDesc2D *Desc2,
             int &Counter, int LowerFirstChild,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers)
{
  static int i, v, w, N_;
  static THangingNode *hn;

  for(i=0;i<N_DOF0;i++)
  {
    Aux[i]=I_K0+Indices0[i];
  }

  for(i=0;i<N_DOF1;i++)
  {
    Aux[i+N_DOF0]=I_K1+Indices1[i];
  }

  for(i=0;i<N_DOF2;i++)
  {
    Aux[i+N_DOF0+N_DOF1]=I_K2+Indices2[i];
  }

  for(i=0;i<N_Pairs;i++)
    if(Pairs[i*2]!=-1)
    {
      // cout << Pairs[i*2] << "  " << Pairs[i*2+1] << endl;
      // cout << *array[Pairs[i*2]] << "   " << *array[Pairs[i*2+1]] << endl;
      MapDOF(Global, Aux[Pairs[i*2]], Aux[Pairs[i*2+1]], Counter);
      // cout << Pairs[i*2] << "  " << Pairs[i*2+1] << endl;
    }

  for(i=0;i<N_Mid;i++)
    if(LowerFirstChild)
      MapDOF(Global, Aux[Mid[i*2]], Aux[Mid[i*2+1]], Counter);
    else
      MapDOF(Global, Aux[Mid[i*2+1]], Aux[Mid[i*2]], Counter);

  for(i=0;i<N_Hanging;i++)
  {
    w=Aux[Hanging[i]];
    while( (v=Global[w]) > -1)
    {
      w=v;
    }
    Global[w]=HANGINGNODE;

    N_=TFEDatabase2D::GetHNDesc2D(HangingTypes[i])->GetN_Nodes();
    hn=new THangingNode(HangingTypes[i], N_, Aux, Coupling[i]);
    vect->AddElement(hn);
    numbers->AddElement(Aux[Hanging[i]]);
  }

  for(i=0;i<N_NoOpposite;i++)
  {
    w=Aux[NoOpposite[i]];
    while( (v=Global[w]) > -1 )
    {
      w=v;
    }

    if( Global[w] == -1 )
    {
      Counter--;
      Global[Aux[NoOpposite[i]]]=Counter;
    }
  } // endfor i
}

/** map the given local degrees of freedom,
    coarse cell has bigger number,
    0 is coarser side 0; 1,2 are on finer side 1 */
void TFE2DMapper1Reg::MapFineCoarse(int *Global, int I_K0, int I_K1, int I_K2,
             int *Indices0, int *Indices1, int *Indices2,
             int LocEdge0, int LocEdge1, int LocEdge2,
             TFEDesc2D *Desc0, TFEDesc2D *Desc1, TFEDesc2D *Desc2,
             int &Counter, int LowerFirstChild,
             TVector<THangingNode *> *vect,
             TVector<int> *numbers)
{
  static int i, v, w, N_;
  static THangingNode *hn;

  // cout << "MapFineCoarse" << endl;

  for(i=0;i<N_DOF0;i++)
    Aux[i]=I_K0+Indices0[i];

  for(i=0;i<N_DOF1;i++)
    Aux[i+N_DOF0]=I_K1+Indices1[i];

  for(i=0;i<N_DOF2;i++)
    Aux[i+N_DOF0+N_DOF1]=I_K2+Indices2[i];

  for(i=0;i<N_Pairs;i++)
    if(Pairs[i*2]!=-1)
    {
      // cout << Pairs[i*2] << "  " << Pairs[i*2+1] << endl;
      MapDOF(Global, Aux[Pairs[i*2+1]], Aux[Pairs[i*2]], Counter);
    }

  for(i=0;i<N_Mid;i++)
    if(LowerFirstChild)
      MapDOF(Global, Aux[Mid[i*2]], Aux[Mid[i*2+1]], Counter);
    else
      MapDOF(Global, Aux[Mid[i*2+1]], Aux[Mid[i*2]], Counter);

  for(i=0;i<N_Hanging;i++)
  {
    w=Aux[Hanging[i]];
    while( (v=Global[w]) > -1)
    {
      w=v;
    }
    Global[w]=HANGINGNODE;

    N_=TFEDatabase2D::GetHNDesc2D(HangingTypes[i])->GetN_Nodes();
    hn=new THangingNode(HangingTypes[i], N_, Aux, Coupling[i]);
    vect->AddElement(hn);
    numbers->AddElement(Aux[Hanging[i]]);
  }

  for(i=0;i<N_NoOpposite;i++)
  {
    w=Aux[NoOpposite[i]];
    while( (v=Global[w]) > -1 )
    {
      w=v;
    }

    if( Global[w] == -1 )
    {
      Counter--;
      Global[Aux[NoOpposite[i]]]=Counter;
    }
  } // endfor i
}
