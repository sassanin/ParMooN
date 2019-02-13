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
// @(#)FE2DMapper.C        1.1 10/30/98
//
// Class:       TFE2DMapper
// Purpose:     build connections between
//
// Author:      Gunar Matthies  05.11.97
//
// =======================================================================

#include <FE2DMapper.h>
#include <FEDatabase2D.h>
#include <HNDesc.h>
#include <MooNMD_Io.h>
#include <string.h>
#include <stdlib.h>

/** constructor, filling all data */
TFE2DMapper::TFE2DMapper(char *name, char *description,int n0, int n1,
                     int n_pairs, int *pairs, 
                     int n_hanging, int *hanging,
                     HNDesc *hangingtypes, int **coupling,
                     int n_farhanging, int *farhanging,
                     HNDesc *farhangingtypes, int ****farcoupling,
                     int n_noopposite, int *noopposite,
                     int n_nodes)
{
  Name=strdup(name);

  Description=strdup(description);

  N_DOF0=n0;

  N_DOF1=n1;

  N_Pairs=n_pairs;

  Pairs=pairs;

  N_Hanging=n_hanging;

  Hanging=hanging;

  HangingTypes=hangingtypes;

  Coupling=coupling;

  N_FarHanging=n_farhanging;

  FarHanging=farhanging;

  FarHangingTypes=farhangingtypes;

  FarCoupling=farcoupling;

  N_NoOpposite=n_noopposite;

  NoOpposite=noopposite;

  Aux=new int[n_nodes];
}

void TFE2DMapper::Map(int *Global, int I_K0, int I_K1, 
                    int *Indices0, int *Indices1,
                    int LocEdge0, int LocEdge1,
                    TFEDesc2D *Desc0, TFEDesc2D *Desc1,
                    int &Counter, TVector<THangingNode *> *vect,
                    TVector<int> *numbers)
{
  static int i, v, w, N_;
  static THangingNode *hn;

  // cout << endl;

  for(i=0;i<N_DOF0;i++)
  {
    Aux[i]=I_K0+Indices0[i];
    // cout << Aux[i] << endl;
  }

  for(i=0;i<N_DOF1;i++)
  {
    Aux[i+N_DOF0]=I_K1+Indices1[i];
    // cout << Aux[i+N_DOF0] << endl;
  }

  for(i=0;i<N_Pairs;i++)
    MapDOF(Global, Aux[Pairs[2*i]], Aux[Pairs[2*i+1]], Counter);

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

void TFE2DMapper::MapBound(int *Global, int I_K, int *Indices, 
                         int &BoundCounter)
{
  static int i,j, dof, v, w;

  for(i=0;i<N_DOF0;i++)
  {
    j=Indices[i];
    dof=I_K+j;
    w=dof;

    while( (v=Global[w]) > -1 )
    {
      w=v;
    }

    if(Global[w] == -1) // there are no information
    {
      BoundCounter--;
      Global[dof]=BoundCounter;
    } // endif
    else
    {
      if(Global[w] < BoundCounter)
      {
        // node is mark as a NO boundary node => modify mark
        // cout << "modify mark from " << J_K[j];
        BoundCounter--;
        Global[w]=BoundCounter;
        if(dof!=w) Global[dof]=w; // do not allow recursive linking
        // cout << " to " << J_K[j] << endl;
      } // endif
    } //endelse
  } // endfor i
}

void TFE2DMapper::MapDOF(int *Global, int dof0, int dof1, int &Counter)
{
  int v0, v1, w0, w1;
  int w,e;

  w0=dof0;
  while( (v0=Global[w0]) > -1 )
  {
    w0=v0;
  }

  w1=dof1;
  while( (v1=Global[w1]) > -1 )
  {
    w1=v1;
  }

  // cout << endl;
  // cout << "dof0: " << dof0 << endl;
  // cout << "Global[dof0]: " << Global[dof0] << endl;
  // cout << "dof1: " << dof1 << endl;
  // cout << "Global[dof1]: " << Global[dof1] << endl;
  // cout << "w0: " << w0 << endl;
  // cout << "Global[w0]: " << Global[w0] << endl;
  // cout << "w1: " << w1 << endl;
  // cout << "Global[w1]: " << Global[w1] << endl;

  if( (Global[dof0] == -1) && (Global[dof1] == -1) )
  {
    // there are no information => lower index dof[01] gets right to define
    Counter--;
    if(dof0<dof1)
    {
      Global[dof0]=Counter; 
      Global[dof1]=dof0;
    }
    else
    {
      Global[dof1]=Counter;
      Global[dof0]=dof1;
    }
    // cout << "both" << endl;
  }
  else
    if(Global[dof0] == -1)
    {
      // => dof0=w0
      // cout << "one 0" << endl;
      if(w1>dof0)
      {
        Global[dof0]=Global[w1];
        Global[dof1]=dof0;
        Global[w1]=dof0;
      }
      else
      {
        Global[dof0]=w1;
      }
    }
    else
      if(Global[dof1] == -1)
      {
        // => dof1=w1
        // cout << "one 1" << endl;
        if(w0>dof1)
        {
          Global[dof1]=Global[w0];
          Global[dof0]=dof1;
          Global[w0]=dof1;
        }
        else
        {
          Global[dof1]=w0;
        }
      }
      else
        /*
        {
          // cout << "information on both sides" << endl;

          if( Global[w0] > Global[w1] ) // 0 was first
          {
            // cout << "(*v0) > (*v1)" << endl;
            if( Global[dof0] < -1)
            {
              // => no link from dof0
              if(dof0<w1)
              {
                Global[w1]=dof0;
                Global[dof1]=dof0;  
              }
              else
              {
                Global[w1]=Global[w0];
                Global[w0]=w1;
              }
            }
            else
            {
              Global[dof0]=w0;
              Global[dof1]=w0;
              Global[w1]=w0;
            }
          }
          else
          {
            if( Global[w0] == Global[w1] )
            {
              // cout << "(*v0) == (*v1)" << endl;
              if( Global[dof0] < -1)
              {
                Global[dof1]=w0;
              }
              else
                if( Global[dof1] < -1)
                {
                  Global[dof0]=w1;
                }
                else
                {
                  Global[dof0]=w1;
                  Global[dof1]=w1;
                }
            }
            else
            {
              // cout << "(*v0) < (*v1)" << endl;
              if( Global[dof1] < -1 )
              {
                if(dof1<w0)
                {
                  Global[w0]=dof1;
                  Global[dof0]=dof1;  
                }
                else
                {
                  Global[w0]=Global[w1];
                  Global[w1]=w0;
                }
              }
              else
              {
                Global[dof0]=w1;
                Global[dof1]=w1;
                Global[w0]=w1;
              }
            }
          }

          // cout << "dof's: " << dof0 << "         " << dof1 << endl;
        }
        */
        if(Global[w0]!=Global[w1])
        {
          e=MAX(Global[w0],Global[w1]);
          w=MIN(w0,w1);

          Global[w0]=w;
          Global[w1]=w;
          Global[w]=e;
        }
        else
        {
          if(w0!=w1) cerr << "ungleich" << endl;
        }

      // cout << "dof's: " << dof0 << "         " << dof1 << endl;
  
  if(Global[dof0]>=dof0)
  {
    cerr << "error with dof0" << endl;
    cerr << "dof0: " << dof0 << endl;
  }
  if(Global[w0]>=w0)
  {
    cerr << "error with w0" << endl;
  }
  if(Global[dof1]>=dof1)
  {
    cerr << "error with dof1" << endl;
    cerr << "dof1: " << dof1 << endl;
  }
  if(Global[w1]>=w1)
  {
    cerr << "error with w1" << endl;
  }
  // cout << "Global[dof0]: " << Global[dof0] << endl;
  // cout << "Global[w0]: " << Global[w0] << endl;
  // cout << "Global[dof1]: " << Global[dof1] << endl;
  // cout << "Global[w1]: " << Global[w1] << endl;
}
