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
// @(#)QuadIsoparametric.C        1.7 04/13/00
//
// Class:      TQuadIsoparametric
//
// Purpose:    reference transformations for triangle
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
// 
// =======================================================================

#include <QuadIsoparametric.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <BoundEdge.h>
#include <InterfaceJoint.h>
#include <BoundComp.h>
#include <string.h>
#include <stdlib.h>
#include <IsoInterfaceJoint.h>
#include <IsoBoundEdge.h>
#include <IsoJointEqN.h>

BaseFunct2D TQuadIsoparametric::BaseFunctFromOrder[] = { 
                BF_C_Q_Q0_2D, BF_C_Q_Q1_2D, BF_C_Q_Q2_2D, BF_C_Q_Q3_2D,
                BF_C_Q_Q4_2D, BF_C_Q_Q5_2D, BF_C_Q_Q6_2D, BF_C_Q_Q7_2D,
                BF_C_Q_Q8_2D, BF_C_Q_Q9_2D };

FEDesc2D TQuadIsoparametric::FEDescFromOrder[] = { 
                FE_C_Q_Q0_2D, FE_C_Q_Q1_2D, FE_C_Q_Q2_2D, FE_C_Q_Q3_2D,
                FE_C_Q_Q4_2D, FE_C_Q_Q5_2D, FE_C_Q_Q6_2D, FE_C_Q_Q7_2D,
                FE_C_Q_Q8_2D, FE_C_Q_Q9_2D };

/** constuctor */
TQuadIsoparametric::TQuadIsoparametric()
{
}

/** transfer from reference element to original element */
void TQuadIsoparametric::GetOrigFromRef(double xi, double eta, 
                                        double &X, double &Y)
{
  int N_Points;
  double Xi[1], Eta[1], Xarray[1], Yarray[1], absdetjk[1];

  N_Points = 1;
  Xi[0] = xi; Eta[0] = eta;
  Xarray[0] = X; Yarray[0] = Y;
  GetOrigFromRef(N_Points, Xi, Eta, Xarray, Yarray, absdetjk);
  X = Xarray[0]; Y = Yarray[0];
}

/** transfer a set of points from reference to original element */
void TQuadIsoparametric::GetOrigFromRef(int N_Points, 
                                double *xi, double *eta,
                                double *X, double *Y, double *absdetjk)
{
  int i, j, k;
  double Xi, Eta;
  double dx1, dx2;
  double dy1, dy2;
  double detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;
  
  for(i=0;i<N_Points;i++)
  {
    Xi = xi[i];
    Eta = eta[i];

    dx1 = xc1 + xc3*Eta;
    dx2 = xc2 + xc3*Xi;
    
    dy1 = yc1 + yc3*Eta;
    dy2 = yc2 + yc3*Xi;
    
    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Xi*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Xi*Eta;

    bf = TFEDatabase2D::GetBaseFunct2D(
              BaseFunctFromOrder[ApproximationOrder]);
    bf->GetDerivatives(D00, Xi, Eta, AuxVector);
    bf->GetDerivatives(D10, Xi, Eta, AuxVector+MaxN_BaseFunctions2D);
    bf->GetDerivatives(D01, Xi, Eta, AuxVector+2*MaxN_BaseFunctions2D);

    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];

      dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

      dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
      dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
    }

    detjk = dx1*dy2 - dx2*dy1;
    absdetjk[i] = fabs(detjk);
  }
}

/** transfer from reference element to original element */
void TQuadIsoparametric::GetOrigFromRef(double *ref, double *orig)
{
  int i, j;
  double X, Y, xi, eta;

  xi = ref[0];
  eta = ref[1];

  j = -1;
  for(i=0;i<N_QuadPoints;i++)
  {
    if(fabs(xi-XI[i])<1e-8 && fabs(eta-ETA[i])<1e-8)
    {
      j = i;
      break;
    }
  }

  if(j==-1)
  {
    cout << "error in TQuadIsoparametric::GetOrigFromRef(3)" << endl;
    return;
  }

  X = xc0 + xc1*xi + xc2*eta + xc3*xi*eta;
  Y = yc0 + yc1*xi + yc2*eta + yc3*xi*eta;

  for(i=0;i<N_AuxPoints;i++)
  {
    X += XDistance[i] * FctValues[j][i];
    Y += YDistance[i] * FctValues[j][i];
  }

  orig[0] = X;
  orig[1] = Y;
}

/** transfer from original element to reference element */
void TQuadIsoparametric::GetRefFromOrig(double X, double Y, 
                                        double &xi, double &eta)
{
  cout << "not implemented yet" << endl;
  xi  = 0.0;
  eta = 0.0;
}

/** transfer from original element to reference element */
void TQuadIsoparametric::GetRefFromOrig(double *orig, double *ref)
{
  cout << "not implemented yet" << endl;
  ref[0] = 0.0;
  ref[1] = 0.0;
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadIsoparametric::GetOrigValues(BaseFunct2D BaseFunct,
                               int N_Points, double *xi, double *eta,
                               int N_Functs, QuadFormula2D formula)
{
  cout << "not implemented yet" << endl;
}

/** calculate functions and derivatives from reference element
    to original element, for all given elements */
void TQuadIsoparametric::GetOrigValues(int N_Sets, BaseFunct2D *BaseFuncts,
                               int N_Points, double *xi, double *eta,
                               QuadFormula2D formula,
                               bool *Needs2ndDer)
{
  int i,j,k,N_, start, end;
  double **refvaluesD00, **origvaluesD00;
  double **refvaluesD10, **origvaluesD10;
  double **refvaluesD01, **origvaluesD01;
  double **refvaluesD20, **origvaluesD20;
  double **refvaluesD11, **origvaluesD11;
  double **refvaluesD02, **origvaluesD02;
  double *refD00, *origD00;
  double *refD10, *origD10;
  double *refD01, *origD01;
  double *refD20, *origD20;
  double *refD11, *origD11;
  double *refD02, *origD02;
  double r20, r11, r02, o20, o11, o02;
  double *aux;
  double GeoData[3][3];
  double Eye[3][3];
  BaseFunct2D BaseFunct;
  int N_Functs;
  bool SecondDer;
  double rec_detjk, a11, a12, a21, a22;

  SecondDer = FALSE;
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();

    refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D00);
    if(refvaluesD00==NULL)
    {
      TFEDatabase2D::GetBaseFunct2D(BaseFunct)->MakeRefElementData(formula);
      refvaluesD00=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                                    formula, D00);
    }
  
    origvaluesD00=TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);
    if(origvaluesD00==NULL)
    {
      origvaluesD00 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD00[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D00, origvaluesD00);
    }
  
    for(j=0;j<N_Points;j++)
    {
      refD00 = refvaluesD00[j];
      origD00 = origvaluesD00[j];

      memcpy(origD00, refD00, N_Functs*SizeOfDouble);
    } // endfor j

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);
    if(origvaluesD10==NULL)
    {
      origvaluesD10 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD10[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D10, origvaluesD10);
    }
  
    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);
    if(origvaluesD01==NULL)
    {
      origvaluesD01 = new double* [MaxN_QuadPoints_2D];
      aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
      for(j=0;j<MaxN_QuadPoints_2D;j++)
        origvaluesD01[j] = aux+j*MaxN_BaseFunctions2D;
      TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D01, origvaluesD01);
    }
  
    if(Needs2ndDer[i])
    {
      SecondDer = TRUE;

      refvaluesD20=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D20);
      origvaluesD20=TFEDatabase2D::GetOrigElementValues(BaseFunct, D20);
      if(origvaluesD20==NULL)
      {
        origvaluesD20 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD20[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D20, origvaluesD20);
      }
    
      refvaluesD11=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D11);
      origvaluesD11=TFEDatabase2D::GetOrigElementValues(BaseFunct, D11);
      if(origvaluesD11==NULL)
      {
        origvaluesD11 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD11[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D11, origvaluesD11);
      }
    
      refvaluesD02=TFEDatabase2D::GetRefElementValues(BaseFunct, 
                                        formula, D02);
      origvaluesD02=TFEDatabase2D::GetOrigElementValues(BaseFunct, D02);
      if(origvaluesD02==NULL)
      {
        origvaluesD02 = new double* [MaxN_QuadPoints_2D];
        aux = new double [MaxN_QuadPoints_2D*MaxN_BaseFunctions2D];
        for(j=0;j<MaxN_QuadPoints_2D;j++)
          origvaluesD02[j] = aux+j*MaxN_BaseFunctions2D;
        TFEDatabase2D::RegisterOrigElementValues(BaseFunct, D02, origvaluesD02);
      }
    } // endfor Needs2ndDer[i]
  } // endfor i
  
  // D10 and D01
  for(i=0;i<N_Sets;i++)
  {
    BaseFunct=BaseFuncts[i];
    N_Functs = TFEDatabase2D::GetBaseFunct2D(BaseFunct)->GetDimension();

    refvaluesD10=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D10);
    origvaluesD10=TFEDatabase2D::GetOrigElementValues(BaseFunct, D10);

    refvaluesD01=TFEDatabase2D::GetRefElementValues(BaseFunct, formula, D01);
    origvaluesD01=TFEDatabase2D::GetOrigElementValues(BaseFunct, D01);

    for(j=0;j<N_Points;j++)
    {
      refD10 = refvaluesD10[j];
      origD10 = origvaluesD10[j];
  
      refD01 = refvaluesD01[j];
      origD01 = origvaluesD01[j];
  
      a11 = xc1 + xc3 * eta[j];
      a21 = xc2 + xc3 * xi[j];
      a12 = yc1 + yc3 * eta[j];
      a22 = yc2 + yc3 * xi[j];
    
      for(k=0;k<N_AuxPoints;k++)
      {
        a11 += XDistance[k] * XiDerValues[j][k];
        a21 += XDistance[k] * EtaDerValues[j][k];
    
        a12 += YDistance[k] * XiDerValues[j][k];
        a22 += YDistance[k] * EtaDerValues[j][k];
      }
    
      rec_detjk = 1 / (a11*a22 - a12*a21);

      for(k=0;k<N_Functs;k++)
      {
        origD10[k]=( a22*refD10[k]-a12*refD01[k]) * rec_detjk;
        origD01[k]=(-a21*refD10[k]+a11*refD01[k]) * rec_detjk;
      } // endfor k
    } // endfor j
  } // endfor i

  // leave if no second derivatives are needed
  if(!SecondDer) return;
  
  // implement calculation of second derivatives later

}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadIsoparametric::GetOrigValues(double xi, double eta, 
                int N_BaseFunct,
                double *uref, double *uxiref, double *uetaref,
                double *uorig, double *uxorig, double *uyorig)
{
  int i, j, k;
  double dx1, dx2;
  double dy1, dy2;
  double rec_detjk;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;


  dx1 = xc1 + xc3*eta;
  dx2 = xc2 + xc3*xi;
    
  dy1 = yc1 + yc3*eta;
  dy2 = yc2 + yc3*xi;
    
  bf = TFEDatabase2D::GetBaseFunct2D(
            BaseFunctFromOrder[ApproximationOrder]);
  bf->GetDerivatives(D10, xi, eta, AuxVector+MaxN_BaseFunctions2D);
  bf->GetDerivatives(D01, xi, eta, AuxVector+2*MaxN_BaseFunctions2D);

  for(k=0;k<N_AuxPoints;k++)
  {
    j = IntAux[k];

    dx1 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dx2 += XDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];

    dy1 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D];
    dy2 += YDistance[k] * AuxVector[j+MaxN_BaseFunctions2D*2];
  }

  rec_detjk = 1/(dx1*dy2 - dx2*dy1);

  // D00
  for(i=0;i<N_BaseFunct;i++)
    uorig[i] = uref[i];

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( dy2*uxiref[i]-dy1*uetaref[i]) * rec_detjk;
    uyorig[i]=(-dx2*uxiref[i]+dx1*uetaref[i]) * rec_detjk;
  } // endfor i
}

/** calculate functions and derivatives from reference element
    to original element */
void TQuadIsoparametric::GetOrigValues(int joint, double zeta,
            int N_BaseFunct,
            double *uref, double *uxiref, double *uetaref,
            double *uorig, double *uxorig, double *uyorig)
{
  int i, j, k;
  double a11, a12, a21, a22, rec_detjk;
  double xi, eta;
  TBaseFunct2D *bf;
  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);
  double valxi[MaxN_BaseFunctions2D];
  double valeta[MaxN_BaseFunctions2D];
  double *values[1];

  switch(joint)
  {
    case 0:
      xi  = zeta; eta = -1;
    break;

    case 1:
      xi = 1; eta = zeta;
    break;

    case 2:
      xi  = -zeta; eta = 1;
    break;

    case 3:
      xi = -1 ; eta = -zeta;
    break;
  }

  a11 = xc1 + xc3*eta;
  a21 = xc2 + xc3*xi;
  a12 = yc1 + yc3*eta;
  a22 = yc2 + yc3*xi;

  values[0] = valxi;
  bf->GetValues(1, &zeta, joint, D10, values);
  values[0] = valeta;
  bf->GetValues(1, &zeta, joint, D01, values);

  // add correction due to isoparamatric boundary approximation
  for(i=0;i<N_AuxPoints;i++)
  {
    k = IntAux[i];
    a11 += XDistance[i] * valxi[k];
    a21 += XDistance[i] * valeta[k];

    a12 += YDistance[i] * valxi[k];
    a22 += YDistance[i] * valeta[k];
  } // endfor i

  rec_detjk = 1 / (a11*a22 - a12*a21);

  // D10 and D01
  for(i=0;i<N_BaseFunct;i++)
  {
    uxorig[i]=( a22*uxiref[i]-a12*uetaref[i]) * rec_detjk;
    uyorig[i]=(-a21*uxiref[i]+a11*uetaref[i]) * rec_detjk;
    uorig[i] = uref[i];
  } // endfor i
}

void TQuadIsoparametric::SetCell(TBaseCell *cell)
{
  int i, j, k, N_;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp2D *comp;
  TInterfaceJoint *interface;
  JointType type;
  double t0, t1, t, dt;
  double xa, ya, xe, ye, xm, ym, xp, yp, dx, dy;
  int compid;
  BoundTypes bdtype;
  TBaseFunct2D *bf;
  TFEDesc2D *fedesc;
  int *JointDOF;
  int N_Vertices;
  TVertex **Vertices;
#ifdef __3D__
  double z[4], zp;
#endif
  int CurvedJoint;

  N_AuxPoints = 0;
  CurvedJoint = -1;

  Cell = cell;

  TFEDatabase2D::GetQuadFormula2D(QuadFormula)
        ->GetFormulaData(N_QuadPoints, W, XI, ETA);

  // cout << endl;
  for(i=0;i<4;i++)
  {
#ifdef __3D__
    Cell->GetVertex(i)->GetCoords(x[i], y[i], z[i]);
#else
    Cell->GetVertex(i)->GetCoords(x[i], y[i]);
#endif
    // cout << setw(20) << x[i] << setw(20) << y[i] << endl;
  }

  xc0=( x[0] + x[1] + x[2] + x[3]) * 0.25;
  xc1=(-x[0] + x[1] + x[2] - x[3]) * 0.25;
  xc2=(-x[0] - x[1] + x[2] + x[3]) * 0.25;
  xc3=( x[0] - x[1] + x[2] - x[3]) * 0.25;

  yc0=( y[0] + y[1] + y[2] + y[3]) * 0.25;
  yc1=(-y[0] + y[1] + y[2] - y[3]) * 0.25;
  yc2=(-y[0] - y[1] + y[2] + y[3]) * 0.25;
  yc3=( y[0] - y[1] + y[2] - y[3]) * 0.25;

  for(i=0;i<4;i++)
  {
    // check whether the joint i is curved
    joint = Cell->GetJoint(i);
    type = joint->GetType();
    if(type == BoundaryEdge || type == InterfaceJoint)
    {
      // cout << "joint: " << i << endl;
      xa = x[i];
      ya = y[i];
      xe = x[(i+1) % 4];
      ye = y[(i+1) % 4];

      if(type == BoundaryEdge)
      {
        boundedge = (TBoundEdge *)(joint);
        comp = boundedge->GetBoundComp();
        boundedge->GetParameters(t0, t1);
      }
      else
      {
        interface = (TInterfaceJoint *)(joint);
        comp = interface->GetBoundComp();
        if(Cell == interface->GetNeighbour(0))
          interface->GetParameters(t0, t1); // correct order
        else
          interface->GetParameters(t1, t0); // backward order

      }

      fedesc = TFEDatabase2D::GetFEDesc2D(
                      FEDescFromOrder[ApproximationOrder]);
      bf = TFEDatabase2D::GetBaseFunct2D(
                      BaseFunctFromOrder[ApproximationOrder]);

      JointDOF = fedesc->GetJointDOF(i);

      dt = (t1-t0)/ApproximationOrder;
      dx = (xe-xa)/ApproximationOrder;
      dy = (ye-ya)/ApproximationOrder;
      for(j=1;j<ApproximationOrder;j++)
      {
        t = t0 + j * dt;
        comp->GetXYofT(t, xp, yp);
        xm = xa + j * dx;
        ym = ya + j * dy;
 
        // cout << "(" << xm << ", " << ym << ")" << endl;
        // cout << "(" << xp << ", " << yp << ")" << endl;

        if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8)
        {
          // curved boundary
          XDistance[N_AuxPoints] = xp - xm;
          YDistance[N_AuxPoints] = yp - ym;
          IntAux[N_AuxPoints] = JointDOF[j];
          // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
          if(N_AuxPoints == 0)
            CurvedJoint = i;

          if(CurvedJoint != i)
          {
            OutPut("There is only one curved joint allowed" << endl);
            exit(-1);
          }

          N_AuxPoints++;
        } // endif 
      } // endfor
    } // endif
    else
    {
      if(type == IsoInterfaceJoint || type == IsoBoundEdge ||
         type == IsoJointEqN)
      {
        switch(type)
        {
          case IsoInterfaceJoint:
            N_Vertices = ((TIsoInterfaceJoint *)joint)->GetN_Vertices();
            Vertices = ((TIsoInterfaceJoint *)joint)->GetVertices();
          break;

          case IsoBoundEdge:
            N_Vertices = ((TIsoBoundEdge *)joint)->GetN_Vertices();
            Vertices = ((TIsoBoundEdge *)joint)->GetVertices();
          break;

          case IsoJointEqN:
            N_Vertices = ((TIsoJointEqN *)joint)->GetN_Vertices();
            Vertices = ((TIsoJointEqN *)joint)->GetVertices();
          break;
	  default:
	    cout<< "Not a 2D BD type "<<endl;
	   break;
        } // endswitch

        // cout << "N_AuxVertices: " << N_Vertices << endl;
        ApproximationOrder = N_Vertices+1;

        fedesc = TFEDatabase2D::GetFEDesc2D(
                        FEDescFromOrder[ApproximationOrder]);
        bf = TFEDatabase2D::GetBaseFunct2D(
                        BaseFunctFromOrder[ApproximationOrder]);

        JointDOF = fedesc->GetJointDOF(i);

        xa = x[i];
        ya = y[i];
        xe = x[(i+1) % 4];
        ye = y[(i+1) % 4];

        dt = (t1-t0)/ApproximationOrder;
        dx = (xe-xa)/ApproximationOrder;
        dy = (ye-ya)/ApproximationOrder;
        for(j=1;j<ApproximationOrder;j++)
        {
          t = t0 + j * dt;
#ifdef __3D__
          Vertices[j-1]->GetCoords(xp, yp, zp);
#else
          Vertices[j-1]->GetCoords(xp, yp);
#endif
          xm = xa + j * dx;
          ym = ya + j * dy;
  
          // cout << "(" << xm << ", " << ym << ")" << endl;
          // cout << "(" << xp << ", " << yp << ")" << endl;
  
          if(fabs(xp-xm) > 1e-8 || fabs(yp-ym) > 1e-8)
          {
            // curved boundary
            XDistance[N_AuxPoints] = xp - xm;
            YDistance[N_AuxPoints] = yp - ym;
            IntAux[N_AuxPoints] = JointDOF[j];
            // cout << setw(3) << N_AuxPoints << setw(3) << JointDOF[j] << endl;
            if(N_AuxPoints == 0)
              CurvedJoint = i;
  
            if(CurvedJoint != i)
            {
              OutPut("Only one curved joint per quadrilateral is allowed" << endl);
              exit(-1);
            }
  
            N_AuxPoints++;
          } // endif 
        } // endfor
      } // endif
    } // endelse
  } // endfor

  if(N_AuxPoints)
  {
    // OutPut("ApproximationOrder: " << ApproximationOrder << endl);
    if(ApproximationOrder == 2)
    {
      XDistance[N_AuxPoints] = XDistance[0]*0.5;
      YDistance[N_AuxPoints] = YDistance[0]*0.5;
      IntAux[N_AuxPoints] = 4;
      N_AuxPoints++;
    } // ApproximationOrder == 2

    if(ApproximationOrder == 3)
    {
      XDistance[N_AuxPoints] = XDistance[0]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/3;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/3;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/3;
      N_AuxPoints++;

      N_AuxPoints -= 4;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 5;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
        break;
      } // endswitch
    } // ApproximationOrder == 3

    if(ApproximationOrder == 4)
    {
      XDistance[N_AuxPoints] = XDistance[0]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*3.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*2.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/4;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*1.0/4;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/4;
      N_AuxPoints++;

      N_AuxPoints -= 9;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 11;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 6;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 17;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 12;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 18;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
        break;
      } // endswitch
    } // ApproximationOrder == 4

    if(ApproximationOrder == 5)
    {
      XDistance[N_AuxPoints] = XDistance[0]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*4.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*4.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*3.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*3.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*2.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*2.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[0]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[0]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[1]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[1]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[2]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[2]*1.0/5;
      N_AuxPoints++;

      XDistance[N_AuxPoints] = XDistance[3]*1.0/5;
      YDistance[N_AuxPoints] = YDistance[3]*1.0/5;
      N_AuxPoints++;

      N_AuxPoints -= 16;

      switch(CurvedJoint)
      {
        case 0:
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
        break;

        case 1:
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
        break;

        case 2:
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
        break;

        case 3:
          IntAux[N_AuxPoints] = 25;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 19;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 13;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 7;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 26;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 20;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 14;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 8;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 27;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 21;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 15;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 9;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 28;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 22;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 16;
          N_AuxPoints++;
          IntAux[N_AuxPoints] = 10;
          N_AuxPoints++;
        break;
      } // switch
    } // ApproximationOrder == 5
  }

  if(N_AuxPoints)
  {
    for(i=0;i<N_QuadPoints;i++)
    {
      bf->GetDerivatives(D00, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        FctValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D10, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        XiDerValues[i][j] = DoubleAux[IntAux[j]]; 

      bf->GetDerivatives(D01, XI[i], ETA[i], DoubleAux);
      for(j=0;j<N_AuxPoints;j++)
        EtaDerValues[i][j] = DoubleAux[IntAux[j]]; 

    } // endfor i
  } // endif
} // TQuadIsoparametric::SetCell

/** return outer normal vector */
void TQuadIsoparametric::GetOuterNormal(int j, double zeta,
                                        double &n1, double &n2)
{
  double len;
  double t1, t2;

  GetTangent(j, zeta, t1, t2);
  len = sqrt(t1*t1+t2*t2);

  n1 =  t2/len;
  n2 = -t1/len;
} // end GetOuterNormal

/** return tangent */
void TQuadIsoparametric::GetTangent(int j, double zeta,
                                        double &t1, double &t2)
{
  double valx[MaxN_BaseFunctions2D];
  double valy[MaxN_BaseFunctions2D];
  double *values[1];
  double xip, etap;
  int i, k;
  double a11, a12, a21, a22, rec_detjk;
  double xi, eta;
  TBaseFunct2D *bf;
  bf = TFEDatabase2D::GetBaseFunct2D(BaseFunctFromOrder[ApproximationOrder]);

  switch(j)
  {
    case 0:
      xi  = zeta; eta = -1;
      xip = 1; etap = 0;
    break;

    case 1:
      xi = 1; eta = zeta;
      xip = 0; etap = 1;
    break;

    case 2:
      xi  = -zeta; eta = 1;
      xip = -1; etap = 0;
    break;

    case 3:
      xi = -1 ; eta = -zeta;
      xip = 0; etap = -1;
    break;
  }

  // cout << xi << " " << eta << endl;

  a11 = xc1 + xc3*eta;
  a21 = xc2 + xc3*xi;
  a12 = yc1 + yc3*eta;
  a22 = yc2 + yc3*xi;

  values[0] = valx;
  bf->GetValues(1, &zeta, j, D10, values);
  values[0] = valy;
  bf->GetValues(1, &zeta, j, D01, values);

  // add correction due to isoparamatric boundary approximation
  for(i=0;i<N_AuxPoints;i++)
  {
    k = IntAux[i];
    a11 += XDistance[i] * valx[k];
    a21 += XDistance[i] * valy[k];

    a12 += YDistance[i] * valx[k];
    a22 += YDistance[i] * valy[k];
  } // endfor i

  // rec_detjk = 1 / (a11*a22 - a12*a21);

  t1 = xip*a11 + etap*a21;
  t2 = xip*a12 + etap*a22;
} // end GetTangent

/** return volume of cell according to isoparametric boundary */
double TQuadIsoparametric::GetVolume()
{
  double locvol;

  int i, j, k;
  double Xi, Eta, a11, a12, a21, a22;
  double absdetjk;

  locvol = 0;
  for(i=0;i<N_QuadPoints;i++)
  {
    Xi = XI[i];
    Eta = ETA[i];

    a11 = xc1 + xc3*Eta;
    a21 = xc2 + xc3*Xi;
    a12 = yc1 + yc3*Eta;
    a22 = yc2 + yc3*Xi;

    for(k=0;k<N_AuxPoints;k++)
    {
      a11 += XDistance[k] * XiDerValues[i][k];
      a21 += XDistance[k] * EtaDerValues[i][k];

      a12 += YDistance[k] * XiDerValues[i][k];
      a22 += YDistance[k] * EtaDerValues[i][k];
    }
    absdetjk = fabs(a11*a22 - a12*a21);

    locvol += absdetjk*W[i];
  } // endfor i
  
  return locvol;
}


/** transfer a  set of boundary points from reference to original element */
void TQuadIsoparametric::GetOrigBoundFromRef( int joint, int N_Points, 
                                              double *zeta, double *X, double *Y)
{
  int i, j, k;
  double Xi, Eta;
  double AuxVector[3*MaxN_BaseFunctions2D];
  TBaseFunct2D *bf;


 for(i=0;i<N_Points;i++)
 {
  switch(joint)
  {
    case 0:
      Xi  = zeta[i]; Eta = -1;
    break;

    case 1:
      Xi = 1; Eta = zeta[i];
    break;

    case 2:
      Xi  = -zeta[i]; Eta = 1;
    break;

    case 3:
      Xi = -1 ; Eta = -zeta[i];
    break;
  }


//    cout<< joint << "  xc0 " << xc0 <<  " Xc1  "<< xc1 << " Xc2  "<< xc2 << endl;

    X[i] = xc0 + xc1*Xi + xc2*Eta + xc3*Xi*Eta;
    Y[i] = yc0 + yc1*Xi + yc2*Eta + yc3*Xi*Eta;

//    cout<< joint << " X[i] " << X[i] <<  " Xi  "<< Xi << endl;

    bf = TFEDatabase2D::GetBaseFunct2D(
              BaseFunctFromOrder[ApproximationOrder]);
    bf->GetDerivatives(D00, Xi, Eta, AuxVector);


    for(k=0;k<N_AuxPoints;k++)
    {
      j = IntAux[k];
      X[i] += XDistance[k] * AuxVector[j];
      Y[i] += YDistance[k] * AuxVector[j];
    }
  }

}
