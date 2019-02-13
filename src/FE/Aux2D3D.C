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
// Class:       TAux2D3D
// Purpose:     calculate parameters from 3d fefuntion for the use
//              in 2d
//
// Author:      Gunar Matthies (15.05.01)
//
// History:     start of implementation 15.05.01 (Gunar Matthies)
// =======================================================================

#include <Aux2D3D.h>
#include <FEDatabase3D.h>
#include <MooNMD_Io.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>
#include <TetraAffin.h>

/** constructor */
TAux2D3D::TAux2D3D(int *cellnumbers, int *jointnumbers,
                   TFEFunction3D *fefunct, int shift)
{
  CellNumbers = cellnumbers;
  JointNumbers = jointnumbers;
  FEFunct = fefunct;

  FESpace = FEFunct->GetFESpace3D();
  Coll = FESpace->GetCollection();

  Values = FEFunct->GetValues();

  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  Shift = shift;
}

/** calculate gradient for local coordinates (xi,eta) on face
    JointNumbers[num] of cell CellNumbers[num] */
void TAux2D3D::GetGradient(int num, int N_Points, double *t1, double *t2,
                           double **Param)
{
  int i,j,k;
  int CellNumber, JointNumber;
  TBaseCell *Cell;
  FE3D FEID;
  TFE3D *ele;
  TBaseFunct3D *bf;
  int Dim;
  RefTrans3D RefTrans;
  TRefTrans3D *rt;
  double n1, n2, n3;
  double uorig[MaxN_BaseFunctions3D], uxorig[MaxN_BaseFunctions3D];
  double uyorig[MaxN_BaseFunctions3D], uzorig[MaxN_BaseFunctions3D];
  double uref[MaxN_BaseFunctions3D], uxiref[MaxN_BaseFunctions3D];
  double uetaref[MaxN_BaseFunctions3D], uzetaref[MaxN_BaseFunctions3D];
  double xi[MaxN_QuadPoints_2D], eta[MaxN_QuadPoints_2D];
  double zeta[MaxN_QuadPoints_2D];
  double u, ux, uy, uz, val;
  int *Numbers;
  double *param;

  JointNumber = JointNumbers[num];
  CellNumber = CellNumbers[num];
  Cell = Coll->GetCell(CellNumber);
  // OutPut("in GetGradient: " << CellNumber << " " << JointNumber << endl);

  FEID = FESpace->GetFE3D(CellNumber, Cell);
  ele = TFEDatabase3D::GetFE3D(FEID);
  bf = TFEDatabase3D::GetBaseFunct3DFromFE3D(FEID);
  Dim = bf->GetDimension();

  RefTrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FEID);

  switch(RefTrans)
  {
    case HexaAffin:
    case HexaTrilinear:
    case HexaIsoparametric:
      switch(JointNumber)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t2[i];
            eta[i] = t1[i];
            zeta[i] = -1;
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t1[i];
            eta[i] = -1;
            zeta[i] = t2[i];
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = t1[i];
            zeta[i] = t2[i];
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -t1[i];
            eta[i] = 1;
            zeta[i] = t2[i];
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        case 4:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = t2[i];
            zeta[i] = t1[i];
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        case 5:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t1[i];
            eta[i] = t2[i];
            zeta[i] = 1;
            // OutPut(i << " " << t[i] << " " << s[i] << " --- ");
            // OutPut(xi[i] << " " << eta[i] << " " << zeta[i] << endl);
          }
        break;

        default:
          Error("Wrong local joint number" << endl);
      }
    break;

    case TetraAffin:
    case TetraIsoparametric:
    break;
  }

  switch(RefTrans)
  {
    case HexaAffin:
      OutPut("HexaAffin" << endl);
      rt = TFEDatabase3D::GetRefTrans3D(HexaAffin);
      ((THexaAffin *)rt)->SetCell(Cell);
      for(i=0;i<N_Points;i++)
      {
        param = Param[i];

        ((THexaAffin *)rt)->GetOuterNormal(JointNumber, t2[i], t1[i],
                                             n1, n2, n3);
        // OutPut(i << " Normal: " << n1 << " " << n2 << " " << n3 << endl);
        bf->GetDerivatives(D000, xi[i], eta[i], zeta[i], uref);
        bf->GetDerivatives(D100, xi[i], eta[i], zeta[i], uxiref);
        bf->GetDerivatives(D010, xi[i], eta[i], zeta[i], uetaref);
        bf->GetDerivatives(D001, xi[i], eta[i], zeta[i], uzetaref);
        ((THexaAffin *)rt)->GetOrigValues(xi[i], eta[i], zeta[i],
                Dim, uref, uxiref, uetaref, uzetaref,
                uorig, uxorig, uyorig, uzorig);
        u = 0;
        ux = 0;
        uy = 0;
        uz = 0;
        Numbers = GlobalNumbers + BeginIndex[CellNumber];
        for(j=0;j<Dim;j++)
        {
          val = Values[Numbers[j]];
          u  +=  uorig[j]*val;
          ux += uxorig[j]*val;
          uy += uyorig[j]*val;
          uz += uzorig[j]*val;
        } // endfor j

        // OutPut("u: " << u << endl);
        // OutPut("grad(u): " << ux << " " << uy << " " << uz << endl);
        // OutPut(Param[i][0] << " " << Param[i][1] << " " << Param[i][2] << endl);

        // OutPut("PARAM: " << (int)param << endl);
        param[Shift] = u;
        param[Shift+1] = ux;
        param[Shift+2] = uy;
        param[Shift+3] = uz;
        param[Shift+4] = n1;
        param[Shift+5] = n2;
        param[Shift+6] = n3;
      } // endfor i
    break;

    case HexaTrilinear:
      // OutPut("HexaTrilinear" << endl);
      rt = TFEDatabase3D::GetRefTrans3D(HexaTrilinear);
      ((THexaTrilinear *)rt)->SetCell(Cell);
      for(i=0;i<N_Points;i++)
      {
        param = Param[i];

        ((THexaTrilinear *)rt)->GetOuterNormal(JointNumber, t2[i], t1[i],
                                             n1, n2, n3);
        // OutPut(i << " Normal: " << n1 << " " << n2 << " " << n3 << endl);
        bf->GetDerivatives(D000, xi[i], eta[i], zeta[i], uref);
        bf->GetDerivatives(D100, xi[i], eta[i], zeta[i], uxiref);
        bf->GetDerivatives(D010, xi[i], eta[i], zeta[i], uetaref);
        bf->GetDerivatives(D001, xi[i], eta[i], zeta[i], uzetaref);
        ((THexaTrilinear *)rt)->GetOrigValues(xi[i], eta[i], zeta[i],
                Dim, uref, uxiref, uetaref, uzetaref,
                uorig, uxorig, uyorig, uzorig);
        u = 0;
        ux = 0;
        uy = 0;
        uz = 0;
        Numbers = GlobalNumbers + BeginIndex[CellNumber];
        for(j=0;j<Dim;j++)
        {
          val = Values[Numbers[j]];
          u  +=  uorig[j]*val;
          ux += uxorig[j]*val;
          uy += uyorig[j]*val;
          uz += uzorig[j]*val;
        } // endfor j

        // OutPut("u: " << u << endl);
        // OutPut("grad(u): " << ux << " " << uy << " " << uz << endl);
        // OutPut(Param[i][0] << " " << Param[i][1] << " " << Param[i][2] << endl);

        // OutPut("PARAM: " << (int)param << endl);
        param[Shift] = u;
        param[Shift+1] = ux;
        param[Shift+2] = uy;
        param[Shift+3] = uz;
        param[Shift+4] = n1;
        param[Shift+5] = n2;
        param[Shift+6] = n3;
      } // endfor i
    break;

    case TetraAffin:
      OutPut("TetraAffin" << endl);
    break;

    case TetraIsoparametric:
    case HexaIsoparametric:
      Error(__FILE__ << ": Not implemented for iso-elements" << endl);
    break;
  }
} // GetGradient
