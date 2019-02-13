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
   
#ifdef __2D__

#include <FreeSurface2D.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <FEVectFunct2D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <InterfaceJoint.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
   
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i)  Diognal = k;
	sum1 -= Entries11[k] * sol[col] 
	      +Entries12[k] * sol[col+N_DOF];
	sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 1000) break;
  } // end while
// OutPut("Grid Solver: Number iteration "<<iter<<endl);
// exit(0);

}


void Solver_3dia(int N_Splines, double *a, double *b, double *c, double *rhs, double *sol)
{
  double *alpha, *beta, *y;
  int i, N;

  N = N_Splines+1;
  alpha = new double[N]; beta = new double[N]; y = new double[N];

  alpha[0] = a[0]; y[0] = rhs[0];
  for(i=1;i<N;i++)
  {
    beta[i] = b[i]/alpha[i-1];
    alpha[i] = a[i]-beta[i]*c[i-1];
    y[i] = rhs[i]-beta[i]*y[i-1];
  }

  sol[N-1] = y[N-1]/alpha[N-1];
  for(i=N-2;i>=0;i--)
    sol[i] = (y[i]-c[i]*sol[i+1])/alpha[i];

  delete [] alpha; delete [] beta; delete [] y;
}




void MoveGrid_imping(double **Entries, double *Sol, double *d, double *Rhs,
                  int *KCol, int *RowPtr,
                  TFEVectFunct2D *GridPos,
                  TFEVectFunct2D *Velocity, double dt,
                  TFEVectFunct2D *NewGridPos, 
                  TVertex ***MovBoundVert, int *N_MovVert,
                  TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                  bool &reparam, int &N_ReParam, TFEVectFunct2D *Stress_FEVectFunc)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength, polydegree;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;

  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y, IsoX, IsoY, Ay;
  double x0, x1, y0, y1, h_tot, res, oldres;
  double *gridvelo, *Nx, *Ny, *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, x_max, x_min, temp, eps=1e-6;  
  
  TIsoBoundEdge *isojoint;
  TMGLevel2D *Level;  
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;  
  QuadFormula2D QuadFormula;
  BoundCond Cond0, Cond1;  
 
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   d = new double[ 2*GridLength];

  if(TDatabase::ParamDB->P5 > 0)
  {  
   Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
   memset(Nx, 0, 2*N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, 2*N_BoundaryNodes*SizeOfDouble);
  }
  //cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
 // cout << GridLength << " " << N_DOF << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;

    // Outward normal no need if we move boundary with velocity
  if(TDatabase::ParamDB->P5 > 0)
  {
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
//             ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
//             ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
}

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();


      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) && (ValuesY[l] != 0) )
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
 
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }                
          }
        }
 //    Due to spline approximation solid boundary end vertices may take negative y value
         if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-12 ) NewValuesX[l] = 0.0;
      } // endfor j
     } // endif
   } // endfor i

/*

   for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
 cout << i <<"  ---  " <<NewValuesX[i] << "  ---  " << NewValuesY[i]<<endl;

exit(0);
*/

   MovBoundVert[0][0]->GetCoords(x, Ay);   
    
   NewGridPos->DataToGrid();
    
   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point
   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0];
   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);
 
  
   // axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   
   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
//    N--;
   
   for(i=1;i<N;i++)
    {
//      MovBoundVert[1][i]->GetCoords(x, y);
//      cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl;      
     y = ((double)(N-i))*h_tot;
     MovBoundVert[1][i]->SetCoords(x, y);   
    }       
 

  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================     
   if(reparam)  
   {   
    // update the iso points and then reparmetrize the free surf
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     if(m == k+2)
      {
       JointDOF = FEDesc->GetJointDOF(j);
       DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
       for(l=0;l<k;l++)
        {
         m = DOF[JointDOF[l+1]];
         Vertices[l]->GetCoords(IsoX, IsoY);
         if(TDatabase::ParamDB->P5 > 0)
          {
           un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
            
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(IsoX, IsoY);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[2];i++)     
//   if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8)
// {
//    cout << "test ReParam_axial3D_U start " << reparam << endl;
// //   exit(0);
// }  
    if(TDatabase::ParamDB->TENSOR_TYPE == 0)
    { 
    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, TRUE);   
    }
     else if(TDatabase::ParamDB->TENSOR_TYPE == 1 || TDatabase::ParamDB->TENSOR_TYPE == 2)
    {
     ReParam_axial3D_UAndStress(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, Stress_FEVectFunc, TRUE); 
    }
//    if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8)
// {
//    cout << "test ReParam_axial3D_U end " << reparam << endl;
//   exit(0);
// }  
    OutPut("ReParam CURRENT TIME: ");
    OutPut(TDatabase::TimeDB->CURRENTTIME << endl);   
   }  //if (reparam)  
   
   NewGridPos->GridToData();
   GridPos->DataToGrid();

   memset(Rhs, 0, 2*GridLength*SizeOfDouble);
   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);


  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);

// for(i=GridLength-N_BoundaryNodes;i<GridLength;i++)
//  cout << i << "  ---  "<<NewValuesX[i] << "  ---  " << NewValuesX[i+GridLength] << endl;
//cout << "test " << endl;
  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]); 
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
 
  if (!reparam)  
   {  
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " "
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {  
              IsoX += dt*ValuesVX[m];
              IsoY += dt*ValuesVY[m];
            }
           if(IsoY<=0) IsoY = 1e-5;
//            if(fabs(IsoX)<1e-12) IsoX = 0.;
           Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j  
   } // if(reparam)
  } // endfor i

   if(reparam)
    {
     N_ReParam++; 
     reparam = FALSE;
    }    
    
//   delete [] d;
  
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny;}
  
} // MoveGrid_imping




void FreeSurf_axial3D_new(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                          double *rhs1, double *rhs2,
                          BoundCondFunct2D *BoundaryCondition,
                          double dt, double *Ucl, TFEFunction2D *Surfact, double *param)
{
  int i, j, k, l, DOF_R, DOF_L, m;
  int *KCol, *RowPtr, *JointDOF, N_DOF;
  int N_LinePoints;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  int count=0, count1=0, count2=0;
  int N_BaseFunct, *N_BaseFuncts;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;  
  int comp, N_U, test_L=1, test_R=1;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;  
    
  double r2, r;  
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz, tmp;
  double val, theta, factor1, factor2, angle;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;  
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  double t0, t1, n0, n1, normn, line_wgt;
  double T_Marangoni, Gamma, ngrad_Gamma, GammaE1, GammaE2;    
  double D = TDatabase::ParamDB->REACTOR_P12 / TDatabase::ParamDB->REACTOR_P16 ;

  
  TBaseCell *cell;
  TFEDesc2D *FeDesc;
  BaseFunct2D *BaseFuncts;
  TCollection *Coll;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  BoundCond Cond0, Cond1;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;  
  
#ifdef __SURFACT__ 
  TFESpace2D *surfactantspace;  
  FE2D TFEId;
  TFEDesc2D *TFeDesc;  
  
  double **Turef, **Tuxiref, **Tuetaref;  
  double Tuorig[MaxN_BaseFunctions2D], Tuxorig[MaxN_BaseFunctions2D];
  double Tuyorig[MaxN_BaseFunctions2D];
  double T_val[3], *S_Values;
  
  
  int  *TGlobalNumbers, *TBeginIndex, local_dof;
  int TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;

// surfactant elasticity E
  double E = TDatabase::ParamDB->REACTOR_P10;
//Equation of state, 0 linear, 1 non-linear
  int EOS = int(TDatabase::ParamDB->REACTOR_P11);
//\Gamma_1/Gamma_\infty

//   double CHAR_L = TDatabase::ParamDB->CHAR_L0;
  
  surfactantspace = Surfact->GetFESpace2D();
  S_Values=Surfact->GetValues();
  TGlobalNumbers = surfactantspace->GetGlobalNumbers();
  TBeginIndex = surfactantspace->GetBeginIndex();  
  
//   Gamma_Max = 0;
//   N_Surf = Surfact->GetLength();

//   for(i=0;i<N_Surf;i++)
//    if(Gamma_Max<S_Values[i]) Gamma_Max=S_Values[i];
// 
//    OutPut("Gamma_Max " << Gamma_Max<<endl);  
  
#endif   
  
 
  
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double Re = TDatabase::ParamDB->RE_NR;
  double We = TDatabase::ParamDB->WB_NR, U;
  double Ca = We/Re, D_Angle;
  double beta = TDatabase::ParamDB->FRICTION_CONSTANT;

 	 double EQ_Angle = TDatabase::ParamDB->EQ_CONTACT_ANGLE;
  
  EQ_Angle = (3.141592654/180)*EQ_Angle;
   
 for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      
      if(joint->GetType() == IsoBoundEdge)
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);
        
        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
     } // endfor j
    
    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
      FEId = fespace->GetFE2D(i, cell);
       
#ifdef __SURFACT__ 
     
      TFEId = surfactantspace->GetFE2D(i, cell);
      
#endif
       
      for(j=0;j<N_IsoJoints;j++)
      {
//      cout << "Cell " << i << " has free surface." << endl;
        IJoint = JointNumbers[j];
        // cout << "joint number: " << IJoint << endl;
        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);
        //   if(y0==0||y1==0)
        //   cout<< " y0= " <<y0<<" y1= "<<y1<<"  x0= "<<x0<<"  x1= "<<x1<<endl;
        //   cout<< " N_LinePoints= " <<N_LinePoints<<endl;
      

	
// 	 Patterened surface
// 	double innerlength = TDatabase::ParamDB->P13;
// 	
//         if(x1<innerlength)
// 	  EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->EQ_CONTACT_ANGLE;
//         else if(x1<(2*innerlength))
// 	  EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->P12;
// 	 else if(x1<(3*innerlength))
// 	  EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->EQ_CONTACT_ANGLE;
// 	else if(x1<(4*innerlength))
// 	 EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->P12;
// 	else if(x1<(5*innerlength))
// 	 EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->EQ_CONTACT_ANGLE;
// 	else 
// 	EQ_Angle = (3.141592654/180)*TDatabase::ParamDB->P12;
// 	
	
     // entries for wetting DOF
      if(y0==0) // right wett point edge (bottom)
       {
        FeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        JointDOF = FeDesc->GetJointDOF(IJoint);
        N_DOF = FeDesc->GetN_JointDOF();
        for(m=0;m<N_DOF;m++)
         {
          DOF_R =  GlobalNumbers[BeginIndex[i]+JointDOF[m]];
          fespace->GetDOFPosition(DOF_R, x, y);
          if(y==0) // right wett point
          {
           U = Ucl[DOF_R];
          switch((int)TDatabase::ParamDB->CONTACT_ANGLE_TYPE)
          {
           case 0:
              D_Angle = EQ_Angle;
           break;

           case 1:
            if(U>0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->AD_CONTACT_ANGLE; } // advanving angle
            else if(U<0.)
	     { D_Angle = (3.141592654/180.)*TDatabase::ParamDB->RE_CONTACT_ANGLE; }// receding angle
            else
             {
              D_Angle = EQ_Angle;
//                        + tanh(50.*U)*(TDatabase::ParamDB->AD_CONTACT_ANGLE
//                                        - TDatabase::ParamDB->RE_CONTACT_ANGLE);
             }
             
             
           break;

           case 2:
//   Hocking's expression
              D_Angle = pow(EQ_Angle, (double)3.0) 
                           + 9.0 * Ca * (fabs(U))* log(beta);

              D_Angle = pow(fabs(D_Angle), 1./3.);
           break;

           case 3:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Jiang et al (1979) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*tanh( 4.96*pow(Ca,0.702) )   );
           break;
           case 4:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Bracke et al (1989) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.*pow(Ca,0.5) )   );
           break;
           case 5:
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
//   Berg et al (1992) expression
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 2.24*pow(Ca,0.54) )   );
           break;

           case 6:
//  Berg et al (1992) expression
              Ca *= fabs(U); // capillary number w.r.t contact line velocity
              D_Angle = acos( cos(EQ_Angle) - (cos(EQ_Angle) + 1.)*( 4.47*pow(Ca,0.42) )   );
           break;

          }
// OutPut("  x= "<< x <<"  y= "<< y << " U " << U<<  " D_Angle: " << (180./Pi)*D_Angle<< endl);
// exit(0);
          param[0] = x;
          param[1] = y;
          param[2] = U;
          r_axial = x;       // r value in the axial symmetric integral
          rhs1[DOF_R] +=  r_axial*((cos(D_Angle))/We);   break;
         }
        }
       }

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
   
#ifdef __SURFACT__ 
       TFEDatabase2D::GetBaseFunct2DFromFE2D(TFEId)->MakeRefElementData(LineQuadFormula);
       TFeDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(TFEId);
       TN_BaseFunct = N_BaseFuncts[TFEId];
       TJointDOF = TFeDesc->GetJointDOF(IJoint);
       TN_DOF_Local = TFeDesc->GetN_JointDOF();
       TDOF = TGlobalNumbers + TBeginIndex[i];
#endif

      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(IJoint, N_LinePoints, zeta, X_B, Y_B);

        break;
      } // endswitch


       uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], LineQuadFormula, IJoint);
       uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D10);
       uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], LineQuadFormula, IJoint, D01);

#ifdef __SURFACT__ 
       Turef = TFEDatabase2D::GetJointValues2D(BaseFuncts[TFEId], LineQuadFormula, IJoint);
       Tuxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D10);
       Tuetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[TFEId], LineQuadFormula, IJoint, D01);       
#endif
  
  
       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
#ifdef __SURFACT__                         
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
#endif 
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
                        
#ifdef __SURFACT__                         
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        TN_BaseFunct, Turef[k], Tuxiref[k], Tuetaref[k],
                        Tuorig, Tuxorig, Tuyorig);
#endif                         
            break;
          } // endswitch

          // modify matrices
         F_K->GetTangent(IJoint, zeta[k], t0, t1);  // old line
         r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
         normn = sqrt(t0*t0+t1*t1);
         n0 =  t1/normn;
         n1 = -t0/normn;
       

#ifdef __SURFACT__   
 
       T_val[0] = 0.; T_val[1] = 0.; T_val[2] = 0.;
       
       for(l=0;l<TN_DOF_Local;l++)
        {
          // assumed that the velo space and 2D surfactant space are same fe space
          local_dof   = TJointDOF[l];
          m = TDOF[local_dof];

//           if(S_Values[m]<0)
// 	  {
// 	   OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<S_Values[m]<<endl);
// 	    S_Values[m] = 0.;
// 	  }
	  
	  
          val = S_Values[m];


          T_val[0] += val*Tuorig[local_dof];  // Surfactant C
          T_val[1] += val*Tuxorig[local_dof];  // C_x
          T_val[2] += val*Tuyorig[local_dof];  // C_y
        } // for(l=0;l<TN_

        if(T_val[0]<0. )
	 {
          OutPut(i<< "x : "<<X_B[k]<< " y: " << Y_B[k] <<"  Surfactant exceeds the reference value, T_val= " <<T_val[0]<<endl);
//        for(l=0;l<TN_DOF_Local;l++)
//         {
//           // assumed that the velo space and 2D surfactant space are same fe space
//           local_dof   = TJointDOF[l];
//           m = TDOF[local_dof];
//           val = S_Values[m];
// 	  
// 	   OutPut(l << " sval  "<<val<<endl);
// 	
// 	   
// 	}
          //numerical correction
          T_val[0]=0.; 
         }  

       //calculate the Marangoni effect
          if(EOS==0)
          {
           Gamma =(1. + E*(D - T_val[0]) ); 
           T_Marangoni = normn*E*(t0*T_val[1] + t1*T_val[2])/We;           
           }
          else
          {
           Gamma =(1. + E*log(1. - T_val[0]) );
   
          if(T_val[0]!=1)
           {T_Marangoni = normn*E*( t0*T_val[1] + t1*T_val[2]) / ( We*(1. - T_val[0]) );}
          else
           {T_Marangoni = 0.; }
     
           //see SolubleSurf JCP paper
           if(Gamma<0.1)
            {
             Gamma = 0.1;
             T_Marangoni = 0.;
            }
          }              
#else
    D = 0.;
    Gamma = 1;
    T_Marangoni = 0.;    
#endif       
          // Multiply with time step dt in the main program not here
          r = normn/We;
          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];

           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;

#ifdef __SURFACT__   
            d1 *= Gamma;
            d2 *= Gamma;

            ngrad_Gamma = n0*T_val[1] + n1*T_val[2];
            GammaE1 = (T_val[1] - ngrad_Gamma*n0)*uorig[l];
            GammaE2 = (T_val[2] - ngrad_Gamma*n1)*uorig[l];
            
            if(EOS==0)
             {
              GammaE1 *=E;
              GammaE2 *=E;
             }
            else
             {
              GammaE1 *=(E/(1.- T_val[0]));
              GammaE2 *=(E/(1.- T_val[0]));
             }
#else
              GammaE1 = 0.;  
              GammaE2 = 0.;  
#endif

// rhs1
            val = r_axial*( (1.-n0*n0)*(d1 - GammaE1) - n0*n1*(d2 -GammaE2) );
            val +=(Gamma*uorig[l]); // due to axialsymmetric
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;

           //Marangoni convection
#ifdef __SURFACT__  
            val = r_axial*t0*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs1[TestDOF] -= val;
#endif
   
// rhs2  
            val =  r_axial*( -n1*n0*(d1 - GammaE1) + (1.-n1*n1)*(d2 - GammaE2) );
            val *= LineWeights[k]*r;
            rhs2[TestDOF] -= val;
 
            // Marangoni convection
#ifdef __SURFACT__  
            val = r_axial*t1*uorig[l];
            val *= LineWeights[k]*T_Marangoni;
            rhs2[TestDOF] -= val;
#endif
            index2 = RowPtr[TestDOF+1];
//               cout << TestDOF  << " RhsA  " << rhs1[TestDOF] << " RhsB  " << rhs2[TestDOF] << endl;
	    
            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              ngrad_ansatz= n0*uxorig[m] + n1*uyorig[m];
              e1 = uxorig[m] - ngrad_ansatz*n0;
              e2 = uyorig[m] - ngrad_ansatz*n1;


              val =(d1 - GammaE1)*e1 + (d2 -GammaE2)*e2 + Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val =(d1 - GammaE1)*e1 + (d2 -GammaE2)*e2;
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j

    } // end (N_IsoJoints > 0)
  } // endfor i
  
 } //FreeSurf_axial3D_new(TSquareMatr



void ReParam_axial3D_U(int N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, 
                       TFEVectFunct2D *Velocity,  TFEFunction2D *Surfactant, 
                       TFEFunction2D *SurfSurfactant, bool UpdateU)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
  int *SurfactGlobalNumbers, *SurfactBeginIndex, *SJointDOF, *SDOF, SN_DOF_Joint, *Surf_DOF;
  int *SurfSurfactGlobalNumbers, *SurfSurfactBeginIndex, *SurfSJointDOF, *SurfSDOF;
  int *SurfSurf_DOF, SurfSN_DOF_Joint;
  
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *Mx, *My,*Mu1, *Mu2, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, *u1_spl, *u2_spl;
  double *ValuesUX, *ValuesUY, *Surfact, *surf_spl, *srhs, *Msurf, surf;
  double *SurfSurfact, *Surfsurf_spl, *Surfsrhs, *SurfMsurf, Surfsurf;
   
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace, *SurfactSpace, *SurfSurfactSpace;
  FE2D FEId, SFEId, SurfSFEId;
  TFE2D *ele, *Sele, *SurfSele;
  TFEDesc2D *FeDesc, *SFeDesc, *SurfSFeDesc;
  TCollection *coll;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

 if (abs(VSP) == 201)
    {ORDER = 2;}
  else if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+1 + N_E*(ORDER-1);

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V];
  
  
  VelocitySpace = Velocity->GetFESpace2D();  
  coll = VelocitySpace->GetCollection();  
  
  if(UpdateU)
  {
   u1rhs = new double[N_Splines+1];
   u2rhs = new double[N_Splines+1];
   u1_spl = new double[N_Splines+1];
   u2_spl = new double[N_Splines+1];
   Mu1 = new double[N_Splines+1];
   Mu2 = new double[N_Splines+1];  
   FEParams = new double [4*2*N_Splines]; // u1, u2, surfact, surfsurfact
   U_DOF = new int[N_V];

   VeloBeginIndex = VelocitySpace->GetBeginIndex();
   VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
   ValuesUX = Velocity->GetValues();
   ValuesUY = ValuesUX + Velocity->GetLength();
   
#ifdef __SURFACT__      
   Surf_DOF = new int[N_V];
   surf_spl = new double[N_V];  
   srhs = new double[N_Splines+1];   
   Msurf = new double[N_Splines+1];
 
   SurfSurf_DOF = new int[N_V];
   Surfsurf_spl = new double[N_V];  
   Surfsrhs = new double[N_Splines+1];   
   SurfMsurf = new double[N_Splines+1];  
   
   
   SurfactSpace = Surfactant->GetFESpace2D();
   SurfactBeginIndex = SurfactSpace->GetBeginIndex();
   SurfactGlobalNumbers = SurfactSpace->GetGlobalNumbers();
   Surfact = Surfactant->GetValues(); 
   
   SurfSurfactSpace = SurfSurfactant->GetFESpace2D();
   SurfSurfactBeginIndex = SurfSurfactSpace->GetBeginIndex();
   SurfSurfactGlobalNumbers = SurfSurfactSpace->GetGlobalNumbers();
   SurfSurfact = SurfSurfactant->GetValues();    
#endif   
// SurfSurfactant   
  }
// if(fabs(TDatabase::TimeDB->CURRENTTIME - 0.00220711)<1e-8 && UpdateU)
// {
//    cout << "test movegrid "  << endl;
//   exit(0);
// }   
   m = 0;
   m1 = 0;
   for(i=0;i<N_E;i++) // i<N_E
   {
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge "<<endl;
      exit(0);
     }

// for velocity
   if(UpdateU)
   {
    FEId = VelocitySpace->GetFE2D(CellNo[i], Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    FeDesc = ele->GetFEDesc2D();   // fe descriptor
    JointDOF = FeDesc->GetJointDOF(EdgeNo[i]);
    N_DOF_Joint = FeDesc->GetN_JointDOF();
    DOF = VeloGlobalNumbers + VeloBeginIndex[CellNo[i]];
   
#ifdef __SURFACT__
    // for surfactant
    SFEId = SurfactSpace->GetFE2D(CellNo[i], Me);
    Sele = TFEDatabase2D::GetFE2D(SFEId);
    SFeDesc = Sele->GetFEDesc2D();   // fe descriptor
    SJointDOF = SFeDesc->GetJointDOF(EdgeNo[i]);
    SN_DOF_Joint = SFeDesc->GetN_JointDOF();
    SDOF = SurfactGlobalNumbers + SurfactBeginIndex[CellNo[i]];    
    
    // for Surfsurfactant
    SurfSFEId = SurfSurfactSpace->GetFE2D(CellNo[i], Me);
    SurfSele = TFEDatabase2D::GetFE2D(SurfSFEId);
    SurfSFeDesc = SurfSele->GetFEDesc2D();   // fe descriptor
    SurfSJointDOF = SurfSFeDesc->GetJointDOF(EdgeNo[i]);
    SurfSN_DOF_Joint = SurfSFeDesc->GetN_JointDOF();
    SurfSDOF = SurfSurfactGlobalNumbers + SurfSurfactBeginIndex[CellNo[i]];       
#endif
    
    if((N_DOF_Joint-1)!=ORDER)
     {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " (N_DOF_Joint-1) " << N_DOF_Joint-1 << " ORDER " << ORDER <<endl;
      exit(0);
     }

    if(i !=N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
     N_DOF_Joint--; // assumed that velocity and surfactant having same no. of dof on edge

     // //   cout << " CellNo[i] " << CellNo[i] << endl;
     for (i3=0;i3<N_DOF_Joint;i3++)
       {
         U_DOF[m1] = DOF[JointDOF[i3]]; // needed for later update
         u1_spl[m1] = ValuesUX[DOF[JointDOF[i3]]];
         u2_spl[m1] = ValuesUY[DOF[JointDOF[i3]]];
   
#ifdef __SURFACT__
         Surf_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
         surf_spl[m1] = Surfact[SDOF[SJointDOF[i3]]];
 
         SurfSurf_DOF[m1] = SurfSDOF[SurfSJointDOF[i3]]; // needed for later update
         Surfsurf_spl[m1] = SurfSurfact[SurfSDOF[SurfSJointDOF[i3]]];	 
#endif 
         m1++;
       }
     } //  if(UpdateU)       
              
    } // for(i=0;i<N_E
  
  
//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

 
  if(m+1!=m1 && UpdateU)
    {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " m " << m << " m1 " << m1 <<endl;
      exit(0);
     }

//  for(i=0;i<N_V;i++)
//    OutPut("OldX: "<< i <<' '<<x[i] <<' '<< y[i] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[0] <<' '<< y[0] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-2] <<' '<< y[N_V-2] <<endl);
//  OutPut("OldX: "<< i <<' '<<x[N_V-1] <<' '<< y[N_V-1] <<endl);
// cout << "Surfact[Surf_DOF[0]] " <<Surfact[Surf_DOF[0]] << " Surfact[Surf_DOF[m1]] " << Surfact[Surf_DOF[m1-1]] <<endl;
  
  h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;  
    b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
    c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, Mx);

  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 4] = Mx[i]*h[i];
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 6] = -Mx[i+1];
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;

   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  
 if(UpdateU)
 {
  // ===============================================================
  // u1 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u1_spl[i-1];
     u1 = u1_spl[i];
     u2 = u1_spl[i+1];

     u1rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u1rhs[0] = u1rhs[1];
   u1rhs[N_Splines] = u1rhs[N_Splines-1];


  
  Solver_3dia(N_Splines, a, b, c, u1rhs, Mu1);     
      
  // ===============================================================
  // u2 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u2_spl[i-1];
     u1 = u2_spl[i];
     u2 = u2_spl[i+1];

     u2rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u2rhs[0] = u2rhs[1];
   u2rhs[N_Splines] = u2rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u2rhs, Mu2);
   
#ifdef __SURFACT__  
// ===============================================================
// surfactant
  for(i=1;i<N_Splines;i++)
   {
     u0 = surf_spl[i-1];
     u1 = surf_spl[i];
     u2 = surf_spl[i+1];

     srhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   srhs[0] = srhs[1];
   srhs[N_Splines] = srhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, srhs, Msurf);
// ===============================================================
// surfactant
  for(i=1;i<N_Splines;i++)
   {
     u0 = Surfsurf_spl[i-1];
     u1 = Surfsurf_spl[i];
     u2 = Surfsurf_spl[i+1];

     Surfsrhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Surfsrhs[0] = Surfsrhs[1];
   Surfsrhs[N_Splines] = Surfsrhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Surfsrhs, SurfMsurf);
  // ===============================================================
#endif  
  
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*8;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2  ] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
   
#ifdef __SURFACT__
    FEParams[ISpline + 4  ] = -Msurf[i]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
    FEParams[ISpline + 5] = Msurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((surf_spl[i+1]-surf_spl[i])/h[i+1]-h[i+1]/6.*(Msurf[i+1]-Msurf[i]))*h[i+1];
  
    FEParams[ISpline + 6  ] = -SurfMsurf[i]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1];
    FEParams[ISpline + 7] = SurfMsurf[i+1]*h[i+1]*h[i+1]/2. +
                          ((Surfsurf_spl[i+1]-Surfsurf_spl[i])/h[i+1]-h[i+1]/6.*(SurfMsurf[i+1]-SurfMsurf[i]))*h[i+1]; 
#endif
  
  }
  // ===================================================================
 } // if(UpdateU)

   teta = 1.0/N_Splines;
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
// int iso = 0;
   for(j=0;j<N_E;j++)
    {
   
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       USpline = (i-1)*8;
       FeDof   = i-1;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
         break;
        }
      }

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

    if(Y < 0 || fabs(Y)<1e-8) Y = 0.0; // no penetration on solid boundary
    cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);

//         if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)

//  if(UpdateU)
//        cout<<"NewX :"<<' '<< iso++ <<' '<<X<< ' '<< Y<< endl;
  
//     OutPut("NewX:"<<' '<< m <<' '<<X<<' '<< Y<<endl);
    m++;

// =========================================================================
 if(UpdateU)
 {
  // for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;
   
#ifdef __SURFACT__
    surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;


    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;      
#endif    
    if(j!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
   
#ifdef __SURFACT__      
      Surfact[Surf_DOF[m1]] = surf;
      SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf; 
#endif         
     }
    m1++;
 }
// ====================================================================
// interpolation for isopoints
// no need if reparam is only for grid velo calculation
// ====================================================================
    Joint = cell[j]->GetJoint(EdgeNo[j]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
        
  if(UpdateU)
   {    
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {
       T = double(m)*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         USpline = (i-1)*8;
         FeDof   = i-1;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
           // further T must be from [0;1] on a subspline
           //          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
           //          cout<< ' ' << T <<endl;
          break;
         }
       }
  
     phi1 = (2.*T*T - 3.*T)*T + 1.;
     phi2 = (-2.*T + 3.)*T*T;
     phi3 = (T*T - 2.*T + 1.)*T;
     phi4 = (T - 1)*T*T;

     X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
         Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
     Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
         Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

     IsoVertices[i3]->SetCoords(X, Y);

//       if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)
//        cout<<"NewX iso:"<<' '<< iso++ <<' '<<X<<' '<< Y<<endl;
//  
     m++;

  // ====================================================================
  // for fe values
  // ====================================================================

    u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
    u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    ValuesUX[U_DOF[m1]] = u0;
    ValuesUY[U_DOF[m1]] = u1;
#ifdef __SURFACT__
   surf = surf_spl[FeDof]*phi1 + surf_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;

    Surfsurf = Surfsurf_spl[FeDof]*phi1 + Surfsurf_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;

    Surfact[Surf_DOF[m1]] = surf;    
    SurfSurfact[SurfSurf_DOF[m1]] = Surfsurf; 
#endif        
    m1++;
// ====================================================================
    }  // for(i3=0;i3<k
    }   // if(k==ORDER-1)    

  }  // if(UpdateU) 
  else if (k==ORDER-1)
  {m += k; }
  
 }  //  for(j=0;j<N_E

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
  
 if(UpdateU)
  { 
   delete [] u1rhs;  delete [] u2rhs; delete [] u1_spl;
   delete [] u2_spl; delete [] Mu1; delete [] Mu2;
   delete [] U_DOF;  delete [] FEParams;  
#ifdef __SURFACT__   
   delete [] Surf_DOF; delete [] surf_spl;  
   delete [] srhs; delete [] Msurf;
   delete [] SurfSurf_DOF; delete [] Surfsurf_spl;  
   delete [] Surfsrhs; delete [] SurfMsurf;  
#endif   
  }
  
}


void ReParam_axial3D_UAndStress(int N_E, TBaseCell **cell, int *EdgeNo,  int *CellNo, 
                       TFEVectFunct2D *Velocity,  TFEVectFunct2D *Stress_FEVectFunc, 
                       bool UpdateU)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, m1, k, i3, USpline, FeDof;
  int *VeloBeginIndex, *VeloGlobalNumbers, *JointDOF, *DOF, N_DOF_Joint, *U_DOF;
  int *StressGlobalNumbers, *StressBeginIndex, *SJointDOF, *SDOF, SN_DOF_Joint, *Stress_DOF;
  
  double *h, *t, u0, u1, u2;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *u1rhs, *u2rhs, *Mx, *My,*Mu1, *Mu2, *Params, *Param9, *FEParams;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, *u1_spl, *u2_spl;
  double *ValuesUX, *ValuesUY, *Tauxx_spl, *Tauxy_spl, *Tauyy_spl, *Tauxxrhs, *Tauxyrhs, *Tauyyrhs, tau11, tau12, tau22;
  double *MTauxx, *MTauxy, *MTauyy, *ValuesTauxx, *ValuesTauxy, *ValuesTauyy;
   
  TIsoBoundEdge *isojoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;
  TFESpace2D *VelocitySpace, *StressSpace;
  FE2D FEId, SFEId;
  TFE2D *ele, *Sele;
  TFEDesc2D *FeDesc, *SFeDesc;
  TCollection *coll;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

 if (abs(VSP) == 201)
    {ORDER = 2;}
  else if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+1 + N_E*(ORDER-1);

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V];
  
  
  VelocitySpace = Velocity->GetFESpace2D();  
  coll = VelocitySpace->GetCollection();  
  
  if(UpdateU)
  {
   u1rhs = new double[N_Splines+1];
   u2rhs = new double[N_Splines+1];
   u1_spl = new double[N_Splines+1];
   u2_spl = new double[N_Splines+1];
   Mu1 = new double[N_Splines+1];
   Mu2 = new double[N_Splines+1];  
   FEParams = new double [5*2*N_Splines]; // 5 fe functions, u1, u2, tauxx, tauxy, tauyy
   U_DOF = new int[N_V];

   VeloBeginIndex = VelocitySpace->GetBeginIndex();
   VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
   ValuesUX = Velocity->GetValues();
   ValuesUY = ValuesUX + Velocity->GetLength();
   
   Tauxxrhs = new double[N_Splines+1];
   Tauxyrhs = new double[N_Splines+1];
   Tauyyrhs = new double[N_Splines+1];
   
   Tauxx_spl = new double[N_Splines+1];
   Tauxy_spl = new double[N_Splines+1];
   Tauyy_spl = new double[N_Splines+1];

   MTauxx = new double[N_Splines+1];
   MTauxy = new double[N_Splines+1];
   MTauyy = new double[N_Splines+1];
   
   Stress_DOF = new int[N_V];  
   
   StressSpace = Stress_FEVectFunc->GetFESpace2D();
   StressBeginIndex = StressSpace->GetBeginIndex();
   StressGlobalNumbers = StressSpace->GetGlobalNumbers();
   
   ValuesTauxx = Stress_FEVectFunc->GetValues();
   ValuesTauxy =  ValuesTauxx + Stress_FEVectFunc->GetLength();
   ValuesTauyy =  ValuesTauxx + 2*Stress_FEVectFunc->GetLength(); 
 
  }

   m = 0;
   m1 = 0;
   for(i=0;i<N_E;i++) // i<N_E
   {
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge "<<endl;
      exit(0);
     }

   if(UpdateU)
   {
    FEId = VelocitySpace->GetFE2D(CellNo[i], Me);
    ele = TFEDatabase2D::GetFE2D(FEId);
    FeDesc = ele->GetFEDesc2D();   // fe descriptor
    JointDOF = FeDesc->GetJointDOF(EdgeNo[i]);
    N_DOF_Joint = FeDesc->GetN_JointDOF();
    DOF = VeloGlobalNumbers + VeloBeginIndex[CellNo[i]];
   
    // for viscoelastic stress
    SFEId = StressSpace->GetFE2D(CellNo[i], Me);
    Sele = TFEDatabase2D::GetFE2D(SFEId);
    SFeDesc = Sele->GetFEDesc2D();   // fe descriptor
    SJointDOF = SFeDesc->GetJointDOF(EdgeNo[i]);
    SN_DOF_Joint = SFeDesc->GetN_JointDOF();
    SDOF = StressGlobalNumbers + StressBeginIndex[CellNo[i]];    

    
    if((N_DOF_Joint-1)!=ORDER)
     {
      // only second order conforming elements implimented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " (N_DOF_Joint-1) " << N_DOF_Joint-1 << " ORDER " << ORDER <<endl;
      exit(0);
     }

    if(i !=N_E-1)// -1 due to end dof will be the start dof of the next edge except on last edge
    {  N_DOF_Joint--; // assumed that velocity and viscoelastic stress having same no. of dof on edge
      SN_DOF_Joint--;
    }
    
     // //   cout << " CellNo[i] " << CellNo[i] << endl;
     
      if(N_DOF_Joint != SN_DOF_Joint)
    {
     cout<<"Both stress and velocity space need to be of same order !!!! \n";
     exit(4711);
    }
     else
    {
     for (i3=0;i3<N_DOF_Joint;i3++)
       {
         U_DOF[m1] = DOF[JointDOF[i3]]; // needed for later update
         u1_spl[m1] = ValuesUX[DOF[JointDOF[i3]]];
         u2_spl[m1] = ValuesUY[DOF[JointDOF[i3]]];
   
         Stress_DOF[m1] = SDOF[SJointDOF[i3]]; // needed for later update
         Tauxx_spl[m1] = ValuesTauxx[SDOF[SJointDOF[i3]]];
         Tauxy_spl[m1] = ValuesTauxy[SDOF[SJointDOF[i3]]];
         Tauyy_spl[m1] = ValuesTauyy[SDOF[SJointDOF[i3]]];
 
         m1++;
       }
    }
   } //  if(UpdateU)       
              
  } // for(i=0;i<N_E
  
  
//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

 
  if(m+1!=m1 && UpdateU)
    {
      // only second order conforming elements implemented
      cout<< " No match in no. velo-nodal functunals and vertices on the free surface edge  "<<endl;
      cout<< " m " << m << " m1 " << m1 <<endl;
      exit(0);
     }
     
  h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;  
    b[i] = h[i]/(h[i]+h[i+1]); // \mu_i in PhD thesis
    c[i] = h[i+1]/(h[i]+h[i+1]); // \lambda_i in PhD thesis
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, Mx);

  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);

  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 4] = Mx[i]*h[i];
    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

//     Params[ISpline + 6] = -Mx[i+1];
    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;

   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  
 if(UpdateU)
 {
  // ===============================================================
  // u1 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u1_spl[i-1];
     u1 = u1_spl[i];
     u2 = u1_spl[i+1];

     u1rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u1rhs[0] = u1rhs[1];
   u1rhs[N_Splines] = u1rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u1rhs, Mu1);     
      
  // ===============================================================
  // u2 component
  // ===============================================================
  for(i=1;i<N_Splines;i++)
   {
     u0 = u2_spl[i-1];
     u1 = u2_spl[i];
     u2 = u2_spl[i+1];

     u2rhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   u2rhs[0] = u2rhs[1];
   u2rhs[N_Splines] = u2rhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, u2rhs, Mu2);
   

  // Tauxx component
  for(i=1;i<N_Splines;i++)
   {
     u0 = Tauxx_spl[i-1];
     u1 = Tauxx_spl[i];
     u2 = Tauxx_spl[i+1];

     Tauxxrhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Tauxxrhs[0] = Tauxxrhs[1];
   Tauxxrhs[N_Splines] = Tauxxrhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Tauxxrhs, MTauxx);
  
    // Tauxy component
  for(i=1;i<N_Splines;i++)
   {
     u0 = Tauxy_spl[i-1];
     u1 = Tauxy_spl[i];
     u2 = Tauxy_spl[i+1];

     Tauxyrhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Tauxyrhs[0] = Tauxyrhs[1];
   Tauxyrhs[N_Splines] = Tauxyrhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Tauxyrhs, MTauxy);
  
      // Tauyy component
  for(i=1;i<N_Splines;i++)
   {
     u0 = Tauyy_spl[i-1];
     u1 = Tauyy_spl[i];
     u2 = Tauyy_spl[i+1];

     Tauyyrhs[i] = 6./(h[i]+h[i+1])*((u2-u1)/h[i+1]-(u1-u0)/h[i]);
    }

   Tauyyrhs[0] = Tauyyrhs[1];
   Tauyyrhs[N_Splines] = Tauyyrhs[N_Splines-1];

  Solver_3dia(N_Splines, a, b, c, Tauyyrhs, MTauyy);

  
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;

    FEParams[ISpline  ] = -Mu1[i]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];
    FEParams[ISpline + 1] = Mu1[i+1]*h[i+1]*h[i+1]/2. +
                          ((u1_spl[i+1]-u1_spl[i])/h[i+1]-h[i+1]/6.*(Mu1[i+1]-Mu1[i]))*h[i+1];


    FEParams[ISpline + 2  ] = -Mu2[i]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    FEParams[ISpline + 3] = Mu2[i+1]*h[i+1]*h[i+1]/2. +
                          ((u2_spl[i+1]-u2_spl[i])/h[i+1]-h[i+1]/6.*(Mu2[i+1]-Mu2[i]))*h[i+1];
    
    FEParams[ISpline + 4  ] = -MTauxx[i]*h[i+1]*h[i+1]/2. +
                          ((Tauxx_spl[i+1]-Tauxx_spl[i])/h[i+1]-h[i+1]/6.*(MTauxx[i+1]-MTauxx[i]))*h[i+1];
    FEParams[ISpline + 5] = MTauxx[i+1]*h[i+1]*h[i+1]/2. +
                          ((Tauxx_spl[i+1]-Tauxx_spl[i])/h[i+1]-h[i+1]/6.*(MTauxx[i+1]-MTauxx[i]))*h[i+1];
			  
    FEParams[ISpline + 6  ] = -MTauxy[i]*h[i+1]*h[i+1]/2. +
                          ((Tauxy_spl[i+1]-Tauxy_spl[i])/h[i+1]-h[i+1]/6.*(MTauxy[i+1]-MTauxy[i]))*h[i+1];
    FEParams[ISpline + 7] = MTauxy[i+1]*h[i+1]*h[i+1]/2. +
                          ((Tauxy_spl[i+1]-Tauxy_spl[i])/h[i+1]-h[i+1]/6.*(MTauxy[i+1]-MTauxy[i]))*h[i+1];
			  
    FEParams[ISpline + 8  ] = -MTauyy[i]*h[i+1]*h[i+1]/2. +
                          ((Tauyy_spl[i+1]-Tauyy_spl[i])/h[i+1]-h[i+1]/6.*(MTauyy[i+1]-MTauyy[i]))*h[i+1];
    FEParams[ISpline + 9] = MTauyy[i+1]*h[i+1]*h[i+1]/2. +
                          ((Tauyy_spl[i+1]-Tauyy_spl[i])/h[i+1]-h[i+1]/6.*(MTauyy[i+1]-MTauyy[i]))*h[i+1];
  }
  // ===================================================================
 } // if(UpdateU)

   teta = 1.0/N_Splines;
   T = 0;

   Param9[0] = 0;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   m1 = 0;
// int iso = 0;
   for(j=0;j<N_E;j++)
    {
   
     T = double(m)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       USpline = (i-1)*10;
       FeDof   = i-1;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
         break;
        }
      }

   phi1 = (2.*T*T - 3.*T)*T + 1.;
   phi2 = (-2.*T + 3.)*T*T;
   phi3 = (T*T - 2.*T + 1.)*T;
   phi4 = (T - 1)*T*T;

   X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
       Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
   Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
       Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

    if(Y < 0 || fabs(Y)<1e-8) Y = 0.0; // no penetration on solid boundary
    cell[j]->GetVertex(EdgeNo[j])->SetCoords(X, Y);


    m++;

// =========================================================================
 if(UpdateU)
 {
  // for fe values

     u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
     u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;
   
     tau11 = Tauxx_spl[FeDof]*phi1 + Tauxx_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;
	      
     tau12 = Tauxy_spl[FeDof]*phi1 + Tauxy_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;
	      
     tau22 = Tauyy_spl[FeDof]*phi1 + Tauyy_spl[FeDof+1]*phi2 +
              FEParams[USpline+8]*phi3 + FEParams[USpline + 9]*phi4;
	      
    if(j!=0) // endpoints no need to set
     {
      ValuesUX[U_DOF[m1]] = u0;
      ValuesUY[U_DOF[m1]] = u1;
        
      ValuesTauxx[Stress_DOF[m1]] = tau11;
      ValuesTauxy[Stress_DOF[m1]] = tau12;
      ValuesTauyy[Stress_DOF[m1]] = tau22;
        
     }
    m1++;
 }
// ====================================================================
// interpolation for isopoints
// no need if reparam is only for grid velo calculation
// ====================================================================
    Joint = cell[j]->GetJoint(EdgeNo[j]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
        
  if(UpdateU)
   {    
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {
       T = double(m)*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         USpline = (i-1)*10;
         FeDof   = i-1;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
           // further T must be from [0;1] on a subspline
           //          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
           //          cout<< ' ' << T <<endl;
          break;
         }
       }
  
     phi1 = (2.*T*T - 3.*T)*T + 1.;
     phi2 = (-2.*T + 3.)*T*T;
     phi3 = (T*T - 2.*T + 1.)*T;
     phi4 = (T - 1)*T*T;

     X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
         Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
     Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
         Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

     IsoVertices[i3]->SetCoords(X, Y);

//       if(fabs(dx0-X)>1e-4 || fabs(dy0-Y)>1e-4)
//        cout<<"NewX iso:"<<' '<< iso++ <<' '<<X<<' '<< Y<<endl;
//  
     m++;

  // ====================================================================
  // for fe values
  // ====================================================================

    u0 = u1_spl[FeDof]*phi1 + u1_spl[FeDof+1]*phi2 +
              FEParams[USpline]*phi3 + FEParams[USpline + 1]*phi4;
    u1 = u2_spl[FeDof]*phi1 + u2_spl[FeDof+1]*phi2 +
              FEParams[USpline+2]*phi3 + FEParams[USpline + 3]*phi4;

    ValuesUX[U_DOF[m1]] = u0;
    ValuesUY[U_DOF[m1]] = u1;

     tau11 = Tauxx_spl[FeDof]*phi1 + Tauxx_spl[FeDof+1]*phi2 +
              FEParams[USpline+4]*phi3 + FEParams[USpline + 5]*phi4;
	      
     tau12 = Tauxy_spl[FeDof]*phi1 + Tauxy_spl[FeDof+1]*phi2 +
              FEParams[USpline+6]*phi3 + FEParams[USpline + 7]*phi4;
	      
     tau22 = Tauyy_spl[FeDof]*phi1 + Tauyy_spl[FeDof+1]*phi2 +
              FEParams[USpline+8]*phi3 + FEParams[USpline + 9]*phi4;
	      
      ValuesTauxx[Stress_DOF[m1]] = tau11;
      ValuesTauxy[Stress_DOF[m1]] = tau12;
      ValuesTauyy[Stress_DOF[m1]] = tau22;

     
    m1++;
// ====================================================================
    }  // for(i3=0;i3<k
    }   // if(k==ORDER-1)    

  }  // if(UpdateU) 
  else if (k==ORDER-1)
  {m += k; }
  
 }  //  for(j=0;j<N_E

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
  
 if(UpdateU)
  { 
   delete [] u1rhs;  delete [] u2rhs; delete [] u1_spl;
   delete [] u2_spl; delete [] Mu1; delete [] Mu2;
   delete [] U_DOF;  delete [] FEParams; 
   
   delete []  Tauxxrhs;
   delete []  Tauxyrhs;
   delete []  Tauyyrhs;
   delete [] Stress_DOF;    delete []  Tauxx_spl;
   delete []  Tauxy_spl;
   delete []  Tauyy_spl;
     delete []  MTauxx;
   delete []  MTauxy;
   delete []  MTauyy;  
   
  }
  
}


void GridVelo_imping(double **Entries, double *Sol, double *d, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, 
                     TVertex ***MovBoundVert, int *N_MovVert,
                     TBaseCell **Free_Cells, int **IsoCellEdgeNos,
                     bool &reparam, TFEVectFunct2D *RefGridPos, TFEVectFunct2D *Stress_FEVectFunc)
{
  int i,j,k,l,m,comp, N;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels, *DOF, *JointDOF, GridLength;
  int N_BoundaryNodes, N_LinePoints, IIso, N_Inner, N_;
  
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *RefValueX, *RefValueY;
  double *ValuesVX, *ValuesVY, *NewValuesX, *NewValuesY;
  double s, t, x, y, h_tot, x0, x1, y0, y1;
  double res, oldres, *gridvelo, *Nx, *Ny, Ay;
  double *LineWeights, *zeta;
  double normalx, normaly, tangenx, tangeny, nx, ny, tx, ty;
  double un, hE, t0,t1, temp2, eps=1e-6;
  double h, hmin, hmax, hlimit, *IsoX, *IsoY;
   
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  TIsoBoundEdge *isojoint;  
  TMGLevel2D *Level;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  BoundTypes bdtype;
  TBoundEdge *BoundEdge;
  TBoundComp2D *BoundComp;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  boolean OnBoundary;
  TJoint *joint;
  TVertex **Vertices;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BoundCond Cond0, Cond1;
    
  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridPos->GridToData();
  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();

  if(TDatabase::ParamDB->P5 > 0)
  {
   Nx = new double[N_BoundaryNodes];
   Ny = new double[N_BoundaryNodes];
   memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
   memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);
  }

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  RefValueX = RefGridPos->GetValues();
  RefValueY = RefValueX + GridLength;   
  
  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
//  cout << "N_Cells: " <<N_Cells<< endl;
  // determine outer normal vectors

  IIso = N_BoundaryNodes;
  // Outward normal no need if we move boundary with velocity
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }
 }
  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = FALSE;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) || 
           (cell->GetJoint(j)->GetType() == InterfaceJoint)  || 
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {
        OnBoundary = TRUE;
        joint = cell->GetJoint(j);
       }
    } // endfor j

    if(OnBoundary)
    {

      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[i];

      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];

        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];

      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if((TDatabase::ParamDB->P5 > 0) )
          { 
           un = VX[j]*Nx[k] + VY[j]*Ny[k];
            
           if(ValuesX[l] == 0 )
            { NewValuesX[l] = ValuesX[l]; }
           else
            { NewValuesX[l] = ValuesX[l] + dt*un*Nx[k]; }

            if(ValuesY[l] == 0) 
             { NewValuesY[l] = ValuesY[l];  }
            else    
            { NewValuesY[l] = ValuesY[l] + dt*un*Ny[k]; }
          }
          else
          {
	    if(ValuesX[l] == 0 )
	    { NewValuesX[l] = ValuesX[l]; }
	    else
	    { NewValuesX[l] = ValuesX[l] + dt*VX[j]; }

            if(ValuesY[l] == 0) 
	     { NewValuesY[l] = ValuesY[l];  }
            else    
	     { NewValuesY[l] = ValuesY[l] + dt*VY[j];  }
	    
	  }
       } //  if(k>=0)
 //    Due to spline approximation solid boundary end vertices may take negative y value
        if(NewValuesY[l]<0.0 ) NewValuesY[l] = 0.0;
        if( fabs(NewValuesX[l]) < 1e-10 ) NewValuesX[l] = 0.0;
      } // endfor j
    } // endif
  } // endfor i

// cout << " dt " << dt <<endl;
/*
   for(i=0;i<GridLength;i++)
 cout << i <<"  ---  " <<ValuesX[i] << "  ---  " << NewValuesX[i] << endl;
*/
// exit(0);

   MovBoundVert[0][0]->GetCoords(x, Ay);   
   AuxGridPos->DataToGrid();
   
  //======================================================================  
  //  Reparametrization of free surface - Begin
  //======================================================================  
  if(!reparam)
  {
   h_tot = 0;
   hmin = 100;
   hmax = 0.0;
   
   for(k=0;k<N_MovVert[2];k++)
    {
     MovBoundVert[2][k]->GetCoords(x1, y1);

     if(k==N_MovVert[2]-1)
     { MovBoundVert[1][0]->GetCoords(x, y);}
     else
     { MovBoundVert[2][k+1]->GetCoords(x, y); }

     h = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
     h_tot +=h;
     if (h < hmin) hmin = h;
     if (h > hmax) hmax = h;
    } // for(k=0;k<N_MovVert[2];k++) 
   
    h_tot /= (double)N_MovVert[2];
    hlimit =  0.8*h_tot;   
   }
   
   if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   { 
 
    //before reparam update iso points, since we use interpolated cubic spline
    //which pass through iso points also
    IsoX = new double[IIso];
    IsoY = new double[IIso];
    
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
      j = IsoCellEdgeNos[1][i];
      joint = cell->GetJoint(j);
      isojoint = (TIsoBoundEdge *)joint;
      k = isojoint->GetN_Vertices();
      Vertices = isojoint->GetVertices();
      FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
      FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
      m = FEDesc->GetN_JointDOF();
      if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[IsoCellEdgeNos[0][i]];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              x  = IsoX[IIso] + dt*un*Nx[IIso+N_BoundaryNodes];
              y  = IsoY[IIso] + dt*un*Ny[IIso+N_BoundaryNodes];
            }
            else
            {
             x  = IsoX[IIso] + dt*ValuesVX[m];
             y  = IsoY[IIso] + dt*ValuesVY[m];
            }
            
           if(y<=0) y = 1e-5;
//            if(fabs(x)<1e-12) x = 0.;
           
           Vertices[l]->SetCoords(x, y);
           IIso++;
          } // endfor l
        }
       else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        } //  if(m == k+2)     
    }// for(i=0;i<N_MovVert[2];i++)

    if(TDatabase::ParamDB->TENSOR_TYPE == 0)
    { 
    ReParam_axial3D_U(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, NULL, NULL, FALSE);  
    }
    else if(TDatabase::ParamDB->TENSOR_TYPE == 1 || TDatabase::ParamDB->TENSOR_TYPE == 2)
    {
     ReParam_axial3D_UAndStress(N_MovVert[2], Free_Cells,  IsoCellEdgeNos[1], IsoCellEdgeNos[0],
                      Velocity, Stress_FEVectFunc, FALSE); 
    }
    reparam = TRUE;   
    RefGridPos->GridToData();    
    Daxpy(2*GridLength, -1, ValuesX, RefValueX); // now reparamdisp in RefValueX

    //back to orig mesh (no reparam movement in calculation of free surf w)
    AuxGridPos->DataToGrid(); 
    
   //restore iso points
    IIso = 0;
    for(i=0;i<N_MovVert[2];i++)
    {
     cell = Free_Cells[i];
     j = IsoCellEdgeNos[1][i];      
     joint = cell->GetJoint(j);
     isojoint = (TIsoBoundEdge *)joint;
     k = isojoint->GetN_Vertices();
     Vertices = isojoint->GetVertices();
     FEId = VelocitySpace->GetFE2D(IsoCellEdgeNos[0][i], cell);
     FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
     m = FEDesc->GetN_JointDOF();
     for(l=0;l<k;l++)
      {
        Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
        IIso++;
      }
    }// for(i=0;i<N_Cells;i++)    
    
    delete [] IsoX; delete [] IsoY; 
   } // if ( ((hmin < hlimit) || (hmax > 3.*h_tot/2.)) ||  reparam )
   //======================================================================       

   MovBoundVert[2][0]->GetCoords(x, y); // right wetting point

   y=0.;
   h_tot = x;
   h_tot /= (double)N_MovVert[0]; 

   for(i=1;i<N_MovVert[0];i++)
     MovBoundVert[0][i]->SetCoords(h_tot*(double)i, y);
 

// axial boundary
   MovBoundVert[1][0]->GetCoords(x, y);
   
   N=N_MovVert[1];      
   h_tot = (y-Ay)/(double)N;   
//    N--;
   
    for(i=1;i<N;i++)
    {
//      MovBoundVert[1][i]->GetCoords(x, y);
     // cout<< " y " << y <<" new y " << y +((double)(N-i))*h_tot<<endl;      
     y = ((double)(N-i))*h_tot;
     MovBoundVert[1][i]->SetCoords(x, y);   
    }      
   
//    exit(0);
   
//    x=0.;
//    h_tot = -y;
//    h_tot /= (double)N_MovVert[1];
//    for(i=1;i<N_MovVert[1];i++)
//     MovBoundVert[1][i]->SetCoords(x,  y +  h_tot*(double)i );
//      
   AuxGridPos->GridToData();   
   GridPos->DataToGrid();
   
   
//       MovBoundVert[1][0]->GetCoords(x, y);
//     cout << " x " << x <<" y " << y <<endl;  
//     exit(0); 
    
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
	 

//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< Rhs[i] << "  ---  " << Rhs[i+GridLength] << endl;
//     
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);

  //======================================================================  
  //  Reparametrization of free surface 
  //======================================================================      
  if(reparam)
  {
   memset(Rhs, 0, 2*GridLength*SizeOfDouble);    
   memcpy(Rhs + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
 
  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), RefValueX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes), RefValueX+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble); 
 
  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);
  
  Dscal(2*GridLength, 1./dt, Sol);
   
  //only inner mesh velo, since U will be interpolated in free surf reparm     
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol, gridvelo);
  Daxpy(GridLength-N_BoundaryNodes, 1., Sol+GridLength, gridvelo+GridLength);  
  
//     cout<< "reparam : "<<endl;
//     exit(0);
  }  
  
//     for(i=0;i<GridLength;i++)
//      cout<< i <<"  ---  "<< gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

//   delete [] d;
   
  if(TDatabase::ParamDB->P5 > 0)
  { delete [] Nx;  delete [] Ny; }
  
  
} // GridVelo_imping





// ====================================================================
// modify matrices due to integrals on free surface
// ====================================================================
void FreeSurfInt(TSquareMatrix2D *A11, TSquareMatrix2D *A12,
                 TSquareMatrix2D *A21, TSquareMatrix2D *A22,
                 double *rhs1, double *rhs2,
                 BoundCondFunct2D *BoundaryCondition,
                 double dt, double factor)
{
  int i,j,k,l,m;
  TBaseCell *cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp;
  double t0, t1, n0, n1, normn;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1;
  int N_BaseFunct, *N_BaseFuncts;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  double val;

  double Ca = TDatabase::ParamDB->P9;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA12 = A12->GetEntries();
  ValuesA21 = A21->GetEntries();
  ValuesA22 = A22->GetEntries();

  for(i=0;i<N_Cells;i++)
  {
    // cout << endl << "CELL number: " << i << endl;
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge) 
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
    } // endfor j

    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
      // cout << "Cell " << i << " has free surface." << endl;
      FEId = fespace->GetFE2D(i, cell);
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
        break;
      } // endswitch

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];

      for(j=0;j<N_IsoJoints;j++)
      {
        IJoint = JointNumbers[j];
        // cout << "joint number: " << IJoint << endl;
        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);

        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;
          // cout << "zeta: " << zeta[k] << endl;
          // cout << "k= " << k << "  tangent: " << t0 << " " << t1 << endl;
          // cout << "length: " << sqrt(t0*t0+t1*t1) << endl;
          uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint);
          uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint, D10);
          uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint, D01);
          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;
          } // endswitch

          // modify matrices
          r2 = 1/(t0*t0+t1*t1);
          r = dt * sqrt(t0*t0+t1*t1)/Ca;
          for(l=0;l<N_BaseFunct;l++)
          {
            TestDOF = DOF[l];

            // updating rhs
            val = -r2* t0 * (uxorig[l]*t0+uyorig[l]*t1);
            val *= LineWeights[k]*r;
            // cout << "Rhs1: " << TestDOF << " " << val;
            rhs1[TestDOF] += val;
            // cout << " = " << rhs1[TestDOF] << endl;

            val = -r2* t1 * (uxorig[l]*t0+uyorig[l]*t1);
            val *= LineWeights[k]*r;
            // cout << "Rhs2: " << TestDOF << " " << val;
            rhs2[TestDOF] += val;
            // cout << " = " << rhs2[TestDOF] << endl;

            val = sqrt(t0*t0+t1*t1) * LineWeights[k] *
                        n1*n1 * factor * uorig[l]*dt;
            // cout << "l= " << l << " bf: " << uorig[l] << endl;
            rhs1[TestDOF] += val*n0;
            rhs2[TestDOF] += val*n1;

            index2 = RowPtr[TestDOF+1];
            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              val = r2*dt*( 
                        (uxorig[m]*t0 + uyorig[m]*t1)  *
                        (uxorig[l]*t0 + uyorig[l]*t1) );
              val *= LineWeights[k]*r;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val = r2*dt*( 
                        (uxorig[m]*t0 + uyorig[m]*t1)  *
                        (uxorig[l]*t0 + uyorig[l]*t1) );
              val *= LineWeights[k]*r;
              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j
    } // end (N_IsoJoints > 0)
    else
    {
      // cout << "Cell " << i << " has NO free surface." << endl;
    }
  } // endfor i
}



 // ====================================================================
// modify matrices due to integrals on the interface for two-phase flows
// axialsymmetric case
// ====================================================================
void Interface_2PhaseSurfAxial3D(TSquareMatrix2D *A11, TSquareMatrix2D *A22,
                     double *rhs1, double *rhs2,
                     BoundCondFunct2D *BoundaryCondition,
                     double dt)
{
  int i, j, k, l, DOF_R, DOF_L, m;
  TBaseCell *cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges, Semi_implicit=0;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp, N_U, test_L=1, test_R=1, ORDER;
  double t0, t1, n0, n1, normn, line_wgt;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  FE2D FEId, TFEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1,tx,ty,mod_t, x, y;
  int N_BaseFunct, *N_BaseFuncts, TN_BaseFunct, *TJointDOF, TN_DOF_Local, *TDOF;
  double **uref, **uxiref, **uetaref;
  double **Turef, **Tuxiref, **Tuetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r, T_val[3], T_Marangoni, *S_Values;
  int *KCol, *RowPtr, *JointDOF, N_DOF, N_Surf;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2, Phase_No;
  double val, theta, factor1, factor2, angle, Gamma;
  int count=0, count1=0, count2=0;
  double  X_B[100], Y_B[100], r_axial, d1, d2, e1, e2, ngrad_test, ngrad_ansatz;
  int  local_dof;
  double Gamma_Max, Gamma_Infty;

  TFEDesc2D *FeDesc, *TFeDesc;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA22 = A22->GetEntries();

  double We = TDatabase::ParamDB->WB_NR;

 for(i=0;i<N_Cells;i++)
  {
//      cout << endl << "CELL number: " << i << endl;
   cell = Coll->GetCell(i);
   Phase_No = cell->GetRegionID();
   if(Phase_No==0)
    {
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
     {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoInterfaceJoint || joint->GetType() == InterfaceJoint)
      {
       FEId = fespace->GetFE2D(i, cell);
       DOF = GlobalNumbers + BeginIndex[i];
       N_BaseFunct = N_BaseFuncts[FEId];
       ele = TFEDatabase2D::GetFE2D(FEId);
       RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
       ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEId);
       l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
//        l = 3;
       LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
       qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
       qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
       TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)->MakeRefElementData(LineQuadFormula);
         
       switch(RefElement)
       {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
          ((TQuadIsoparametric *)F_K)->GetOrigBoundFromRef(j, N_LinePoints, zeta, X_B, Y_B);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
          ((TTriaIsoparametric *)F_K)->GetOrigBoundFromRef(j, N_LinePoints, zeta, X_B, Y_B);
        break;
       } // endswitch


     uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId],
                        LineQuadFormula, j);
     uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, j, D10);
     uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId],
                        LineQuadFormula, j, D01);

       for(k=0;k<N_LinePoints;k++)
        {
         switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
             break;
            case BFUnitTriangle:

              ((TTriaIsoparametric *)F_K)->GetOrigValues(j,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
	      break;
          } // endswitch

          // modify matrices
         F_K->GetTangent(j, zeta[k], t0, t1);  // old line
         if(TDatabase::ParamDB->Axial3D==1)
	 {
	   r_axial = fabs(X_B[k]);   // r value in the axial symmetric integral
	 if(X_B[k]<=0)
          {
           cell->GetVertex(j)->GetCoords(x0, y0);       
           cell->GetVertex((j+1) % N_Edges)->GetCoords(x1, y1);    
           cout <<"X_B[k] negative in freesurface int, change Quad rule " <<  X_B[k] <<endl;
           cout << j << " " << k <<" X  " << x0<<" Y  " << y0<<" X  " << x1<<" Y  " << y1 <<endl;
           cout << "If the quad formula is correct, then the time step might be large !!!! " <<endl; 
           exit(0);
          }
	 }
	 else if(TDatabase::ParamDB->Axial3D==0)
	 {
	   r_axial = 1.0;
	 }
          
         normn = sqrt(t0*t0+t1*t1);
	 n0 =  t1/normn;
         n1 = -t0/normn;
       
    Gamma = 1;   
           
          r = normn/We;           
          for(l=0;l<N_BaseFunct;l++)
          {
           TestDOF = DOF[l];
           
           // updating rhs
            ngrad_test= n0*uxorig[l] + n1*uyorig[l];
            d1 = uxorig[l] - ngrad_test*n0;
            d2 = uyorig[l] - ngrad_test*n1;
            
	                 
// rhs1
            val = r_axial*( (1.-n0*n0)*d1 - n0*n1*d2);  
	    if(TDatabase::ParamDB->Axial3D==1)
	    {
            val +=(Gamma*uorig[l]); // due to axialsymmetric
	    }
            val *= LineWeights[k]*r;
            rhs1[TestDOF] -= val;
// rhs2
            val =  r_axial*( -n1*n0*d1 + (1.-n1*n1)*d2 );    
            val *= LineWeights[k]*r;
            rhs2[TestDOF] -= val;

            index2 = RowPtr[TestDOF+1];

            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              ngrad_ansatz= n0*uxorig[m] + n1*uyorig[m];
              e1 = uxorig[m] - ngrad_ansatz*n0;
              e2 = uyorig[m] - ngrad_ansatz*n1;

              val =d1*e1 + d2*e2 ;
	      if(TDatabase::ParamDB->Axial3D==1)
	     {
              val+= Gamma*(uorig[l]*uorig[m]/(r_axial*r_axial));
	     }
              val *= dt*LineWeights[k]*r*r_axial;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val =d1*e1 + d2*e2;
              val *= dt*LineWeights[k]*r*r_axial;

              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;

            } // endfor m
          } // endfor l
        } // endfor k
       } // 
     } // endfor j
     }    // if(Phase_No==0)
   } // endfor i
 }




// ====================================================================
// determine grid velocity in whole domain
// ====================================================================
void GetGridVelocity(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     double *Nx, double *Ny,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  double *Rhs, *Sol;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  int IIso;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_;
  double un;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
  for(i=0;i<N_Cells;i++)
  {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k 

        DOF = GridGlobalNumbers + GridBeginIndex[i];

/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/
// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

/* 
        // not needed
        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            IIso++;
          } // endfor
        } // endif
*/
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
            NewValuesX[l] = ValuesX[l] + dt*VX[j];
            NewValuesY[l] = ValuesY[l] + dt*VY[j];
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  N_Levels = GridMG->GetN_Levels();
  Level = GridMG->GetLevel(N_Levels-1);
  Rhs = Level->GetRhs();
  Sol = Level->GetSolution();
  
  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesX, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesX, Sol, GridLength*SizeOfDouble);

  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesY, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesY, Sol, GridLength*SizeOfDouble);

  // put solution into grid position
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;
      
      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
  } // endfor i

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
  Dscal(2*GridLength, 1/dt, gridvelo);

  // for(i=0;i<GridLength;i++)
  //   cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

} // GetGridVelocity

// ====================================================================
// determine new grid position
// ====================================================================
void MoveGrid(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
              double *IsoX, double *IsoY,
              double *Nx, double *Ny,
              TFEVectFunct2D *Velocity, double dt,
              TFEVectFunct2D *NewGridPos, 
              double *NewIsoX, double *NewIsoY,
              TFEVectFunct2D *GridVelocity)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  double *Rhs, *Sol;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TIsoBoundEdge *isojoint;
  TVertex **Vertices;
  int IIso;
  double *gridvelo;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_;
  double un;
  int polydegree;
  QuadFormula2D QuadFormula;

  if(TDatabase::ParamDB->P5>0)
  {
    OutPut("Boundary update using (u.n)n" << endl);
  }
  else
  {
    OutPut("Boundary update using u" << endl);
  }

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - N_Inner;
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
            ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
            ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k 

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
            NewValuesX[l] = ValuesX[l] + dt*VX[j];
            NewValuesY[l] = ValuesY[l] + dt*VY[j];
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  N_Levels = GridMG->GetN_Levels();
  Level = GridMG->GetLevel(N_Levels-1);
  Rhs = Level->GetRhs();
  Sol = Level->GetSolution();
  
  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesX, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesX, Sol, GridLength*SizeOfDouble);

  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesY, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE
    oldres = res;
    j++;
  }

  memcpy(NewValuesY, Sol, GridLength*SizeOfDouble);

  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;
      
      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch

    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes] 
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX[IIso] += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY[IIso] += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " " 
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {
              IsoX[IIso] += dt*ValuesVX[m];
              IsoY[IIso] += dt*ValuesVY[m];
            }
            Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
            IIso++;
          } // endfor l
        }
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j
  } // endfor i

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
  Dscal(2*GridLength, 1/dt, gridvelo);

} // MoveGrid


void GetGridVelocity(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, int *Velo_CellNo)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints, Phase_No;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty;
  int IIso, Velo_N_Cells;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_i;
  double un, hE;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells = Velo_Coll->GetN_Cells();


  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
//   cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   cout << GridLength << " N_BoundaryNodes:" << endl;
// 
// exit(0);

  d = new double[2*GridLength];
  Nx = new double[N_BoundaryNodes];
  Ny = new double[N_BoundaryNodes];

  memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
  memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);

//   cout << GridLength << " N_BoundaryNodes:" << endl;
  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();


//   for(i=0;i<N_Cells;i++)
//     {
// cout << "cell : " << i << " number " << Velo_CellNo[i] <<endl;
// }

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
//     cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();


    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
            || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
             || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
      {
        // cout << "joint: " << j << endl;
        Phase_No = cell->GetRegionID();
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
//        
        Velo_i = cell->GetGlobalCellNo();

        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
  if(Phase_No==1)  // outer phase has inward normal
   {
    Nx[i] /= -t;
    Ny[i] /= -t;
    }
   else
    {
    Nx[i] /= t;
    Ny[i] /= t;
    }

  }
 }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {

    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
         || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//           cout<<"Cell : "<<i<<"\n";
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();

      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
         if(TDatabase::ParamDB->P5 > 0)
         {
          un = VX[j]*Nx[k] + VY[j]*Ny[k];
          cout << l <<"  ---  "<< VX[j] << "  ---  " << VY[j] << endl;
          cout << "update normal direction is not good in rising bubble" << endl;
          NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
          NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
         }
         else
         {
//           if(fabs(ValuesX[l])<1e-8 && fabs(VX[j])>1e-16)   
//            cout << " x " << ValuesX[l] <<  " y " << ValuesY[l] << " VX " << VX[j] << endl;

          if(ValuesX[l]!=0)
           NewValuesX[l] = ValuesX[l] + dt*VX[j]; //axial bd no need x update
          NewValuesY[l] = ValuesY[l] + dt*VY[j];
         }
        }
      } // endfor j
    } // endif
  } // endfor i


  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);


//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

  delete [] d;
  delete [] Nx;
  delete [] Ny;

//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

} // GetGridVelocity

void GetGridVelo_outer(double **Entries, double *Sol, double *Rhs,
                       int *KCol, int *RowPtr,
                       TFEVectFunct2D *GridPos,
                       TFEVectFunct2D *AuxGridPos,
                       TFEVectFunct2D *Velocity, double dt,
                       TFEVectFunct2D *GridVelocity, int *Velo_CellNo)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints, Phase_No;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty;
  int IIso, Velo_N_Cells;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_i;
  double un, hE;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells = Velo_Coll->GetN_Cells();


  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
//   cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   cout << GridLength << " N_BoundaryNodes:" << endl;


  d = new double[2*GridLength];
  Nx = new double[N_BoundaryNodes];
  Ny = new double[N_BoundaryNodes];

  memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
  memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);

//   cout << GridLength << " N_BoundaryNodes:" << endl;
  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();


//   for(i=0;i<N_Cells;i++)
//     {
// cout << "cell : " << i << " number " << Velo_CellNo[i] <<endl;
// }

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
//     cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();


    for(j=0;j<N_Edges;j++)
    {
     if( !(cell->GetJoint(j)->InnerJoint()) 
        || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
        || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
      {
        // cout << "joint: " << j << endl;
        Phase_No = cell->GetRegionID();
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
  if(Phase_No==1)  // outer phase has inward normal
   {
    Nx[i] /= -t;
    Ny[i] /= -t;
    }
   else
    {
    Nx[i] /= t;
    Ny[i] /= t;
    }

  }
 }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
        || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
        || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();
      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
         if(TDatabase::ParamDB->P5 > 0)
         {
          un = VX[j]*Nx[k] + VY[j]*Ny[k];
          cout << l <<"  ---  "<< VX[j] << "  ---  " << VY[j] << endl;
          cout << "update normal direction is not good in raising bubble" << endl;
          NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
          NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
         }
         else
         {
/*          if(fabs(ValuesX[l])<1e-8 && fabs(VX[j])>1e-16)   
           cout << "GetGridVelo_outer x " << ValuesX[l] <<  " y " << ValuesY[l] << " VX " << VX[j] << endl;*/  
          if(ValuesX[l]!=0) 
           NewValuesX[l] = ValuesX[l] + dt*VX[j];
          NewValuesY[l] = ValuesY[l] + dt*VY[j];
         }
        }
      } // endfor j
    } // endif
  } // endfor i


  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs,
                    KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1/dt, gridvelo);


//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

  delete [] d;
  delete [] Nx;
  delete [] Ny;

//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

} // GetGridVelocity




void MoveGrid_2Phase(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *NewGridPos,
                     TFEVectFunct2D *Velocity, double dt,  
                     int *Velo_CellNo, int isoupdate)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int *DOF, *JointDOF;
  TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TIsoInterfaceJoint *isojoint;
  TVertex **Vertices;
  int IIso;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_N_Cells;
  double un;
  int polydegree, Velo_i;
  QuadFormula2D QuadFormula;
  double IsoX, IsoY;


  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells =  Velo_Coll->GetN_Cells();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - N_Inner;
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
  d = new double[ 2*GridLength];

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);
  Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
  Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5>0)
  {
   for(i=0;i<N_Cells;i++)
    {
     cell = Coll->GetCell(i);
     N_Edges = cell->GetN_Edges();

     for(j=0;j<N_Edges;j++)
      {
       if( !(cell->GetJoint(j)->InnerJoint())
            || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell == Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
            ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
            ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge ||
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint))
        {
          FEId = VelocitySpace->GetFE2D(Velo_CellNo[i], cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
} // if(TDatabase::ParamDB->P5 > 0)


//          if(TY[0]>=0 && TY[1]>=0 )
//             cout<< " x0: "<< setw(10)<< TX[0]<< setw(10)<<" x1: " << TX[1]<< setw(10)<<
// 	                 "y1 - y0: "<< setw(10) << TX[1] - TX[0] <<endl;



  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
         || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();
      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
           if(ValuesX[l]!=0)
             NewValuesX[l] = ValuesX[l] + dt*VX[j];
           NewValuesY[l] = ValuesY[l] + dt*VY[j]; 
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs,
                    KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);


/*
  for(i=0;i<GridLength;i++)
   cout << i << "  ---  "<<Sol[i] << "  ---  " << Sol[i+GridLength] << endl;
*/
  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch


  if(isoupdate)
   {
    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge || joint->GetType() == IsoInterfaceJoint)
      {
        isojoint = (TIsoInterfaceJoint *)joint;
        k = isojoint->GetN_Vertices();
//         cout <<" k " << k << endl;
        Vertices = isojoint->GetVertices();

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[Velo_i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
    
           if(TDatabase::ParamDB->P5>0)
            {
             un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                 + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
             IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
             IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
           else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
       
//            if(IsoX<=0.01 && IsoY<0.5)
//             {
// //              cout << IsoY << "  IsoX less than zero !!!!!!!!!!"  << IsoX << endl;
//               cell->GetVertex(j)->GetCoords(x0, y0);
//               cell->GetVertex((j+1) % N_Edges)->GetCoords(x1, y1);           
//  
//              cout << "X: " << x0  << " " <<  IsoX  << " " <<x1 << " "  << (x0+x1)/2.<< " "  << IsoX - (x0+x1)/2.  << endl; 
//              cout << "Y: " << y0  << " " <<  IsoY  << " " <<y1 << " "  <<  (y0+y1)/2. << " "  <<  IsoY - (y0+y1)/2. << endl; 
// 	     
//              cout << "U1: " << ValuesVX[DOF[JointDOF[l]]]  << " " <<  ValuesVX[DOF[JointDOF[l+1]]]  << " " << ValuesVX[DOF[JointDOF[l+2]]]  << " " << endl; 
//              cout << "U2: " << ValuesVY[DOF[JointDOF[l]]]  << " " <<  ValuesVY[DOF[JointDOF[l+1]]]  << " " << ValuesVY[DOF[JointDOF[l+2]]]  << endl;    
//             }
       
            Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;

          } // endfor l
        }
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j
  }
  } // endfor i

 // gridvelo = GridVelocity->GetValues();
//  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
//  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
//  Dscal(2*GridLength, 1/dt, gridvelo);
  delete [] d;
  delete [] Nx;
  delete [] Ny;
} // MoveGrid

// ====================================================================
// Get the inner angles of the cells in whole domain
// ====================================================================
void Getcellangle(TFESpace2D *Space, double *MinMaxAngle)
{
 int i,j,k,l, N_Cells, N_Edges;
 int found,  N_LinePoints;

 double TX[4], TY[4], hE[4], Theta, tx, ty, Test, MQI=0.;
 TBaseCell *cell;
 FE2D FEId;
 BF2DRefElements RefElement;
 TRefTrans2D *F_K;
 RefTrans2D RefTrans;
 TCollection *Cells;

  MinMaxAngle[0] = 180;  // Min_Angel = 180
  MinMaxAngle[1] = 0;  // Max_Angel = 0
  Cells = Space->GetCollection();
  N_Cells = Cells->GetN_Cells();
     
//      TX      = new double[4];  // Max no edges in 2d
//      TY      = new double[4];  // Max no edges in 2d
//      hE      = new double[4];  // Max no edges in 2d

  for(i=0;i<N_Cells;i++)
   {
     cell    = Cells->GetCell(i);
     N_Edges = cell->GetN_Edges();

     FEId = Space->GetFE2D(i, cell);
     RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

     switch(RefElement)
        {
         case BFUnitTriangle:

            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin *)F_K)->SetCell(cell);

          break;

          case BFUnitSquare:

            RefTrans = QuadAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadAffin *)F_K)->SetCell(cell);

          break;

          default:
            Error("only triangles and quadrilaterals are allowed" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
          } // endswitch

     for(j=0;j<N_Edges;j++)
      {
        F_K->GetTangent(j, 0, tx, ty);
        TX[j] = tx;
        TY[j] = ty;
        hE[j] = sqrt(tx*tx+ty*ty);

   // cout <<"cell : " <<i << "  j= " << j << ": " <<TX[j]<< "------ " << TY[j] << endl;
       } // endfor j

//      Test = 0;
      k = N_Edges -1;
      for(j=0;j<N_Edges;j++)
      {
       if(hE[j]==0.0 || hE[k]== 0.0 )
        Theta = 0.0;
       else
        Theta = acos(-(TX[j]*TX[k]+TY[j]*TY[k])/(hE[j]*hE[k]))*(180/3.141592654);

       k = j;
//        Test +=Theta;
       if(MinMaxAngle[0]>Theta) MinMaxAngle[0] = Theta;
       if(MinMaxAngle[1]<Theta) MinMaxAngle[1] = Theta;
//        cout <<"cell : " <<i << "  j= " << j << ": " << " Theta : " << Theta << endl;
//  *****************************************************
//  Grid test

      MQI += (60. - Theta)*(60. - Theta);
//  *****************************************************

     }
//       cout <<"cell : " <<i <<  " sum of 3 angels : " << Test << endl;
     //  cout<<endl;

   } // endfor i

   MQI /=double(3*N_Cells);
   MQI = sqrt(MQI);

// OutPut("Mesh Quality Indicator: "<< MQI<< endl);
//    delete [] TX;
//    delete [] TY;
//    delete [] hE;
 //cout<< " Min_Angel: "<< MinMaxAngle[0]<< "  Max_Angel : "<<MinMaxAngle[1]<< endl;
// exit(0);
}



double Volume(TFESpace2D *FESpace)
{
  TCollection *Coll;
  double vol, locvol;
  TBaseCell *cell;
  int i, j, N_Cells, N_Edges;
  FE2D FEId;
  RefTrans2D RefTrans;
  TRefTrans2D *rt;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  int polydegree;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  vol = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEId);
    N_Edges = cell->GetN_Edges(); 

    IsIsoparametric = false;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = true;
      }
    } // endif 

    if(IsIsoparametric)
    {
      switch(N_Edges)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    rt = TFEDatabase2D::GetRefTrans2D(RefTrans);
    switch(RefTrans)
    {
      case TriaAffin:
        ((TTriaAffin *)rt)->SetCell(cell);
        locvol = ((TTriaAffin *)rt)->GetVolume();
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        locvol = ((TTriaIsoparametric *)rt)->GetVolume();
      break;

      case QuadAffin:
        ((TQuadAffin *)rt)->SetCell(cell);
        locvol = ((TQuadAffin *)rt)->GetVolume();
      break;

      case QuadBilinear:
        ((TQuadBilinear *)rt)->SetCell(cell);
        locvol = ((TQuadBilinear *)rt)->GetVolume();
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        locvol = ((TQuadIsoparametric *)rt)->GetVolume();
      break;
    }

    vol += locvol;
  } // endfor i

  return vol;
}


void ReParametrize_pts(int &N_Edges, TBaseCell **cell, int *EdgeNo, double h_min, double **FreePts)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, k, i3, N_E = N_Edges;
  double *h, *t;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *Mx, *My, *Params, *Param9;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, Isoteta;
  TIsoBoundEdge *isojoint;
  TIsoInterfaceJoint *isoIntjoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

    if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+ N_E*(ORDER-1) + 1; // add end point for spline

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V];

   m = 0;
   for(i=0;i<N_E;i++)
   {
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge in Remesc2D.C "<<endl;
      exit(0);
     }
   }

//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

 h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;
    b[i] = h[i]/(h[i]+h[i+1]);
    c[i] = h[i+1]/(h[i]+h[i+1]);
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);
  
  Solver_3dia(N_Splines, a, b, c, rhs, Mx);


  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);
  
  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;
  
   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }

  /*
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    //cout<<i+1<< " subspline"<<endl;
    cout<<setprecision(15);
    cout<<"  "<<Params[ISpline    ]<<'\t'<<Params[ISpline + 1]<<endl;
    cout<<"  "<<Params[ISpline + 2]<<'\t'<<Params[ISpline + 3]<<endl;
    cout<<"  "<<Params[ISpline + 4]<<'\t'<<Params[ISpline + 5]<<endl;
    cout<<"  "<<Params[ISpline + 6]<<'\t'<<Params[ISpline + 7]<<endl;
    cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  */

   delete [] FreePts[0];
   delete [] FreePts[1];

   N_E = int(t[N_Splines]/h_min);
  // minimum points will not be less than the initial number of points (used in stbubble problem)
  //    if(N_E<int(TDatabase::ParamDB->P6))
 //         N_E = int(TDatabase::ParamDB->P6);

   N_Edges = N_E;
   FreePts[0] = new double[N_E+1];
   FreePts[1] = new double[N_E+1];

   teta = 1.0/N_E;

   T = 0.;
   Param9[0] = 0.;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   for(j=0;j<=N_E;j++)
    {
     T = double(j)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
//           cout<< ISpline << ' ' << T;
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//           cout<< ' ' << T;
         break;
        }
      }
  
    phi1 = (2.*T*T - 3.*T)*T + 1.;
    phi2 = (-2.*T + 3.)*T*T;
    phi3 = (T*T - 2.*T + 1.)*T;
    phi4 = (T - 1)*T*T;

    X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
        Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
    Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
        Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

    if(fabs(X)<1.e-12)
      X = 0.;
    
    FreePts[0][j] = X;
    FreePts[1][j] = Y;
//        cout<< j  << " X " << X << " Y " << Y<<endl; 
    if(j<N_E)
     {
     // need to set it for surface interpolation maping
     Joint = cell[j]->GetJoint(EdgeNo[j]);
     isoIntjoint = (TIsoInterfaceJoint *)Joint;
     k = isojoint->GetN_Vertices();
     IsoVertices = isojoint->GetVertices();
     Isoteta = 1./(double)(k+1);

     for(i3=0;i3<k;i3++)
      {   
       T = ( (double)(j) + ((double)(i3+1.))*Isoteta )*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
      // further T must be from [0;1] on a subspline
//          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//          cout<< ' ' << T <<endl;
          break;
         }
	}
      phi1 = (2.*T*T - 3.*T)*T + 1.;
      phi2 = (-2.*T + 3.)*T*T;
      phi3 = (T*T - 2.*T + 1.)*T;
      phi4 = (T - 1)*T*T;

      X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
	  Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
      Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
	  Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

      IsoVertices[i3]->SetCoords(X, Y);
      //     cout<< j << "ISO X " << X << " ISO Y " << Y<<endl;
     } //   for(i3=0;i3<k;i3
    } // if(j<N_E)

   }  // for j
//   end vertex of the freeboundary
 

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
}




void IntUn(TFEVectFunct2D *u, double *Nx, double *Ny)
{
  int i,j,k,l,m,n;
  TBaseCell *cell;
  TCollection *coll;
  int N_Cells, N_Edges, N_DOF;
  int N_LinePoints;
  double *zeta, *LineWeights;
  double *u1, *u2;
  TFESpace2D *USpace;
  TJoint *joint;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  double **JointValues, *Values;
  double x0, x1, y0, y1, s, t;
  double un, int_un, loc_un;
  double n1, n2, t1, t2;
  double u1loc, u2loc;
  int *DOF;
  int *BeginIndex, *GlobalNumbers;

  USpace = u->GetFESpace2D();
  coll = USpace->GetCollection();
  N_Cells = coll->GetN_Cells();
  N_DOF = u->GetLength();
  u1 = u->GetValues();
  u2 = u1+N_DOF;
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();

  int_un = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        loc_un = 0;

        FEId = USpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
        TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                    ->MakeRefElementData(LineQuadFormula);
        BF = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId);
        N_DOF = TFEDatabase2D::GetN_BaseFunctFromFE2D(FEId);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch RefElement

        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

        JointValues = TFEDatabase2D::GetJointValues2D(BF, LineQuadFormula, j);
        for(k=0;k<N_LinePoints;k++)
        {
          un = 0;
          Values = JointValues[k];
          F_K->GetOuterNormal(j, zeta[k], n1, n2);
          u1loc = u2loc = 0;
          for(l=0;l<N_DOF;l++)
          {
            // cout << "l= " << l << ": " << Values[l];
            // cout << "zeta[k]: " << zeta[k] << endl;
            m = DOF[l];
            u1loc += u1[m]*Values[l];
            u2loc += u2[m]*Values[l];
          } // endfor l
          un += u1loc*n1 + u2loc*n2;
          F_K->GetTangent(j, zeta[k], t1, t2);
          loc_un += sqrt(t1*t1+t2*t2)*LineWeights[k]*un;
        } // endfor k
        int_un += loc_un;
      } // endif
    } // endfor j
  } // endfor i

  OutPut("Int(u.n) = " << int_un << endl);
  // exit(-1);
}

#endif
