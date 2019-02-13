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
// Class:      TBaseFunct3D
//
// Purpose:    represents the set of base functions for a finite element
//             in three dimensions
//
// Author:     Gunar Matthies
//
// History:    19.11.99 start implementation
// 
// =======================================================================

#include <Constants.h>
#include <BaseFunct3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <JointEqN.h>

/** constructor, fill in all information */
TBaseFunct3D::TBaseFunct3D(int dimension, BaseFunct3D basefunct,
                           BF3DRefElements refelement, 
                           DoubleFunct3D* functions, 
                           DoubleFunct3D* derivativesXi,
                           DoubleFunct3D* derivativesEta,
                           DoubleFunct3D* derivativesZeta,
                           DoubleFunct3D* derivativesXiXi,
                           DoubleFunct3D* derivativesXiEta,
                           DoubleFunct3D* derivativesXiZeta,
                           DoubleFunct3D* derivativesEtaEta,
                           DoubleFunct3D* derivativesEtaZeta,
                           DoubleFunct3D* derivativesZetaZeta,
                           int polynomialdegree,
                           int accuracy,
                           int n_bf2change,
                           int ***bf2change,
                           int baseVectDim)
{
  Dimension=dimension;

  RefElement = refelement;

  BaseFunct = basefunct;

  Functions[D000]=functions;
  Functions[D100]=derivativesXi;
  Functions[D010]=derivativesEta;
  Functions[D001]=derivativesZeta;
  Functions[D200]=derivativesXiXi;
  Functions[D110]=derivativesXiEta;
  Functions[D101]=derivativesXiZeta;
  Functions[D020]=derivativesEtaEta;
  Functions[D011]=derivativesEtaZeta;
  Functions[D002]=derivativesZetaZeta;

  changable = false;

  PolynomialDegree = polynomialdegree;
  Accuracy = accuracy;

  N_BF2Change = n_bf2change;
  BF2Change = bf2change;
  BaseVectDim = baseVectDim;
}

/** constructor without filling data structure */
TBaseFunct3D::TBaseFunct3D(int dimension)
{
  Dimension = dimension;

  changable = true;
}

/** return the values for derivative MultiIndex at all
    quadrature points */
void TBaseFunct3D::GetDerivatives(MultiIndex3D MultiIndex, 
                        TQuadFormula3D *formula, double **values)
{
  int i, N_;
  int j;
  double *w, *xi, *eta, *zeta;

  formula->GetFormulaData(N_, w, xi, eta, zeta);

  for(i=0;i<N_;i++)
    GetDerivatives(MultiIndex, xi[i], eta[i], zeta[i], values[i]);
}

/** return values on joint */
void TBaseFunct3D::GetValues(int N_Points, double *t, double *s, 
                             int joint, double **Values)
{
  int i;
  double xi[MaxN_QuadPoints_2D], eta[MaxN_QuadPoints_2D];
  double zeta[MaxN_QuadPoints_2D];

  switch(RefElement)
  {
    case BFUnitTetrahedron:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = s[i];
            zeta[i] = 0;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = 0;
            zeta[i] = t[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = 1-t[i]-s[i];
            zeta[i] = s[i];
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 0;
            eta[i] = t[i];
            zeta[i] = s[i];
          }
          break;
      }
      break;
    case BFUnitHexahedron:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = s[i];
            zeta[i] = -1;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = -1;
            zeta[i] = t[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = s[i];
            zeta[i] = t[i];
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -s[i];
            eta[i] = 1;
            zeta[i] = t[i];
          }
          break;
        case 4:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = t[i];
            zeta[i] = s[i];
          }
          break;
        case 5:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = t[i];
            zeta[i] = 1;
          }
          break;      } // endswitch
      break;
  } // endswitch

  for(i=0;i<N_Points;i++)
    GetDerivatives(D000, xi[i], eta[i], zeta[i], Values[i]);
}

/** return values of derivative index on joint */
void TBaseFunct3D::GetValues(int N_Points, double *t, double *s, int joint, 
                             MultiIndex3D index, double **Values)
{
  int i;
  double xi[MaxN_QuadPoints_2D], eta[MaxN_QuadPoints_2D];
  double zeta[MaxN_QuadPoints_2D];

  switch(RefElement)
  {
    case BFUnitTetrahedron:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = s[i];
            zeta[i] = 0;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = 0;
            zeta[i] = t[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = 1-t[i]-s[i];
            zeta[i] = s[i];
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 0;
            eta[i] = t[i];
            zeta[i] = s[i];
          }
          break;
      }
      break;
    case BFUnitHexahedron:
      switch(joint)
      {
        case 0:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = t[i];
            eta[i] = s[i];
            zeta[i] = -1;
          }
          break;
        case 1:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = -1;
            zeta[i] = t[i];
          }
          break;
        case 2:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = 1;
            eta[i] = s[i];
            zeta[i] = t[i];
          }
          break;
        case 3:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -s[i];
            eta[i] = 1;
            zeta[i] = t[i];
          }
          break;
        case 4:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = -1;
            eta[i] = t[i];
            zeta[i] = s[i];
          }
          break;
        case 5:
          for(i=0;i<N_Points;i++)
          {
            xi[i] = s[i];
            eta[i] = t[i];
            zeta[i] = 1;
          }
          break;      
      } // endswitch
      break;
  } // endswitch

  for(i=0;i<N_Points;i++)
    GetDerivatives(index, xi[i], eta[i], zeta[i], Values[i]);
}

/** set function for derivative MultiIndex */
void TBaseFunct3D::SetFunction(MultiIndex3D MultiIndex, 
                               DoubleFunct3D* function)
{
  if(changable)
    Functions[MultiIndex] = function;
}

/** make data on reference element */
void TBaseFunct3D::MakeRefElementData(QuadFormula2D FaceQF)
{
  int i, j, N_Joints;
  double **Values, *AllValues;
  TQuadFormula2D *qf2;
  double *w, *t, *s;
  int N_Points;

  qf2 = TFEDatabase3D::GetQuadFormula2D(FaceQF);
  qf2->GetFormulaData(N_Points, w, t, s);

  switch(RefElement)
  {
    case BFUnitTetrahedron:
      N_Joints = 4;
      break;
    case BFUnitHexahedron:
      N_Joints = 6;
      break;
    default:
      N_Joints = 0;
  }

  // joint values
  for(i=0;i<N_Joints;i++)
  {
    Values=TFEDatabase3D::GetJointValues3D(BaseFunct, FaceQF, i);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension;
      GetValues(N_Points, t, s, i, Values);
      TFEDatabase3D::RegisterJointValues3D(BaseFunct, FaceQF, 
                       i, Values);
    }

    Values=TFEDatabase3D::GetJointDerivatives3D(BaseFunct, FaceQF, i, D100);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension;
      GetValues(N_Points, t, s, i, D100, Values);
      TFEDatabase3D::RegisterJointDerivatives3D(BaseFunct, FaceQF, 
                       i, D100, Values);
    } // endif

    Values=TFEDatabase3D::GetJointDerivatives3D(BaseFunct, FaceQF, i, D010);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension;
      GetValues(N_Points, t, s, i, D010, Values);
      TFEDatabase3D::RegisterJointDerivatives3D(BaseFunct, FaceQF, 
                       i, D010, Values);
    } // endif

    Values=TFEDatabase3D::GetJointDerivatives3D(BaseFunct, FaceQF, i, D001);
    if(Values == NULL)
    {
      // data not generated yet
      Values = new double* [N_Points];
      AllValues = new double [N_Points*Dimension];
      for(j=0;j<N_Points;j++)
        Values[j] = AllValues+j*Dimension;
      GetValues(N_Points, t, s, i, D001, Values);
      TFEDatabase3D::RegisterJointDerivatives3D(BaseFunct, FaceQF, 
                       i, D001, Values);
    } // endif
  } // endfor
}

/** make data on reference element */
void TBaseFunct3D::MakeRefElementData(QuadFormula3D QuadFormula)
{
  int j;
  double **Values, *AllValues;
  TQuadFormula3D *qf2;

  qf2 = TFEDatabase3D::GetQuadFormula3D(QuadFormula);
  int n_points = qf2->GetN_QuadPoints();

  // D000
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D000);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension;
    GetDerivatives(D000, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D000, Values);
  }

  // D100
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D100);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension;
    GetDerivatives(D100, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D100, Values);
  }

 // D010
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D010);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension;
    GetDerivatives(D010, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D010, Values);
  }

  // D001
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D001);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension;
    GetDerivatives(D001, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                          D001, Values);
  }
  
  // D200
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D200);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D200, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D200, Values);
  }

  // D110
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D110);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D110, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D110, Values);
  }

  // D101
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D101);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D101, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D101, Values);
  }
  
  // D020
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D020);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D020, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D020, Values);
  }
  
  // D011
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D011);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D011, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D011, Values);
  }
  
  // D002
  Values=TFEDatabase3D::GetRefElementValues(BaseFunct, QuadFormula, D002);
  if( Values==NULL)
  {
    Values = new double* [n_points];
    AllValues = new double [n_points*Dimension*BaseVectDim];
    for(j=0;j<n_points;j++)
      Values[j] = AllValues+j*Dimension*BaseVectDim;
    GetDerivatives(D002, qf2, Values);
    TFEDatabase3D::RegisterRefElementValues(BaseFunct, QuadFormula, 
                                            D002, Values);
  }
}

TGridCell *TBaseFunct3D::GenerateRefElement()
{
  TGridCell *Cell;
  TVertex *v[8];
  TJointEqN *b[6];
  int i;

  switch(RefElement)
  {
    case BFUnitHexahedron:
      Cell = new TGridCell(TDatabase::RefDescDB[Hexahedron], 0);
      v[0] = new TVertex(-1.0, -1.0, -1.0);
      v[1] = new TVertex( 1.0, -1.0, -1.0);
      v[2] = new TVertex( 1.0,  1.0, -1.0);
      v[3] = new TVertex(-1.0,  1.0, -1.0);
      v[4] = new TVertex(-1.0, -1.0,  1.0);
      v[5] = new TVertex( 1.0, -1.0,  1.0);
      v[6] = new TVertex( 1.0,  1.0,  1.0);
      v[7] = new TVertex(-1.0,  1.0,  1.0);
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<8;i++)
        Cell->SetVertex(i, v[i]);

      for(i=0;i<6;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetJoint(i, b[i]);
      }

    break;

    case BFUnitTetrahedron:
      Cell = new TGridCell(TDatabase::RefDescDB[Tetrahedron], 0);
      v[0] = new TVertex(0.0, 0.0, 0.0);
      v[1] = new TVertex(1.0, 0.0, 0.0);
      v[2] = new TVertex(0.0, 1.0, 0.0);
      v[3] = new TVertex(0.0, 0.0, 1.0);
    
      Cell->SetClipBoard(Refinement);

      for(i=0;i<4;i++)
      {
        b[i] = new TJointEqN(Cell);
        Cell->SetVertex(i, v[i]);
        Cell->SetJoint(i, b[i]);
      } // endfor i
    break;

    default:
      cout << "unknown reference element!" << endl;
  }

  return Cell;
}

/** change basis functions on cell if needed */
void TBaseFunct3D::ChangeBF(TCollection *Coll, TBaseCell *Cell, double *Values)
{
  int i, j, k, maptype, *JointArray;
  TJoint *joint;
  TBaseCell *neigh;

  if(BF2Change == NULL) return;

  switch(RefElement)
  {
    case BFUnitTetrahedron:
    break;

    case BFUnitHexahedron:
      for(i=0;i<6;i++)
      {
        joint = Cell->GetJoint(i);
        maptype = joint->GetMapType();
        if(maptype == 1 || maptype == 2)
        {
          neigh = (TGridCell *)(joint->GetNeighbour(Cell));
          if(neigh != NULL && neigh-Cell < 0)
          {
            JointArray = BF2Change[0][i];
            for(j=0;j<N_BF2Change;j++)
            {
              k = JointArray[j];
              Values[k] = -Values[k];
            }
          }
        }
        if(maptype == 2 || maptype == 3)
        {
          neigh = (TGridCell *)(joint->GetNeighbour(Cell));
          if(neigh != NULL && neigh-Cell < 0)
          {
            JointArray = BF2Change[1][i];
            for(j=0;j<N_BF2Change;j++)
            {
              k = JointArray[j];
              Values[k] = -Values[k];
            }
          }
        }
      } // endfor i
    break;
  } // end switch
}

// #ifdef _HYBRID
// /** change basis functions on cell if needed */
// void TBaseFunct3D::ChangeBF(TCollection *Coll, TBaseCell *Cell, double *Values, int &maptype, TJoint *joint)
// {
//   int i, j, k, *JointArray;
//   TBaseCell *neigh;
// 
//   if(BF2Change == NULL) return;
// 
//   switch(RefElement)
//   {
//     case BFUnitTetrahedron:
//     break;
// 
//     case BFUnitHexahedron:
//       for(i=0;i<6;i++)
//       {
//         joint = Cell->GetJoint(i);
//         maptype = joint->GetMapType();
//         if(maptype == 1 || maptype == 2)
//         {
//           neigh = (TGridCell *)(joint->GetNeighbour(Cell));
//           if(neigh != NULL && neigh-Cell < 0)
//           {
//             JointArray = BF2Change[0][i];
//             for(j=0;j<N_BF2Change;j++)
//             {
//               k = JointArray[j];
//               Values[k] = -Values[k];
//             }
//           }
//         }
//         if(maptype == 2 || maptype == 3)
//         {
//           neigh = (TGridCell *)(joint->GetNeighbour(Cell));
//           if(neigh != NULL && neigh-Cell < 0)
//           {
//             JointArray = BF2Change[1][i];
//             for(j=0;j<N_BF2Change;j++)
//             {
//               k = JointArray[j];
//               Values[k] = -Values[k];
//             }
//           }
//         }
//       } // endfor i
//     break;
//   } // end switch
// }
// #endif

/** change basis functions on cell in all points if needed */
void TBaseFunct3D::ChangeBF(TCollection *Coll, TBaseCell *Cell, int N_Points, double **Values)
{
  int i, j, k, l, maptype, *JointArray;
  TJoint *joint;
  TBaseCell *neigh;
  double *Array;

  if(BF2Change == NULL) return;

  switch(RefElement)
  {
    case BFUnitTetrahedron:
    break;

    case BFUnitHexahedron:
      for(i=0;i<6;i++)
      {
        joint = Cell->GetJoint(i);
        maptype = joint->GetMapType();
        if(maptype == 1 || maptype == 2)
        {
          neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
          if(neigh != NULL && neigh-Cell < 0)
          {
            JointArray = BF2Change[0][i];
            for(l=0;l<N_Points;l++)
            {
              Array = Values[l];
              for(j=0;j<N_BF2Change;j++)
              {
                k = JointArray[j];
                Array[k] = -Array[k];
              }
            }
          }
        }
        if(maptype == 2 || maptype == 3)
        {
          neigh = (TBaseCell *)(joint->GetNeighbour(Cell));
          if(neigh != NULL && neigh-Cell < 0)
          {
            JointArray = BF2Change[1][i];
            for(l=0;l<N_Points;l++)
            {
              Array = Values[l];
              for(j=0;j<N_BF2Change;j++)
              {
                k = JointArray[j];
                Array[k] = -Array[k];
              }
            }
          }
        }
      } // endfor i
    break;
  } // end switch
}
