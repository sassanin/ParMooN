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
// @(#)FEFunction1D.C
// 
// Class:       TFEFunction1D
// Purpose:     a function from a finite element space in 1D
// 
// Author:      Sashikumaar Ganesan (17.05.2007)
//
// History:     start of implementation 17.05.2007
//
// =======================================================================

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEFunction1D.h>
#include <string.h>
#include <AllRefTrans.h>
#include <MooNMD_Io.h>

#include <BaseFunct1D.h>
#include <NodalFunctional1D.h>
#include <stdlib.h>
/** constructor with vector initialization */
TFEFunction1D::TFEFunction1D(TFESpace1D *fespace1D, char *name, 
                             char *description, double *values, int length)
{
  FESpace1D=fespace1D;

  Name=strdup(name);

  Description=strdup(description);

  Values=values;

  Length=length;
}

/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(DoubleFunct2D *Exact)
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double X[MaxN_PointsForNodal1D], Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;


  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  memset(Values, 0, SizeOfDouble*N_DOFs);
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);

    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, X, 
#ifdef __2D__
                                       Y,
#endif
                                       AbsDetjk);

    for(j=0;j<N_Points;j++)
    {
      FctVal[0] = double(cell->GetRegionID());
      Exact(X[j], Y[j], FctVal);
      PointValues[j] = FctVal[0];
    }

    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];


    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];

  } // for(i=0;i<N_Cells

}


// Interpolate given 1D + 1d function, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(int ConstCoord, double x, DoubleFunct2D *Exact)
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();


  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Y, AbsDetjk);

    for(j=0;j<N_Points;j++)
     {
      if(ConstCoord==0)
       {
        Exact(x, Y[j], FctVal);
        //cout << " x " << x << " y " <<Y[j] << endl;
       }
      else
       {
        Exact(Y[j], x, FctVal);
       }

      PointValues[j] = FctVal[0];
     }


    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];


    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];

   } // for(i=0; i<N_Cells; i++)

}

// Interpolate given 2D + 1d function, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::Interpolate(double x, double y, DoubleFunct3D *Exact)
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();


  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);

    for(j=0;j<N_Points;j++)
     {
      Exact(x, y, Z[j], FctVal);
      //cout << " x " << x << " y " <<Z[j] << endl;
      PointValues[j] = FctVal[0];
     }

    nf->GetAllFunctionals(PointValues, FunctionalValues);
    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];

   } // for(i=0; i<N_Cells; i++)

}


// Interpolate given 2D/3D + 1d function's at nodal interpolation points, PBE
/** calculate the interpolation of an exact function */
void TFEFunction1D::InterpolateNodalPts(int N_Coord, double *Coords, DoubleFunctND *Exact, double *val)
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points, disp, N_GlobNodalPts;
  int *DOF, *IndexArray, *NodalPtIndex;

  double *xi, *eta;
  double Z[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4], x, y, z;

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  x = Coords[0];
  y = Coords[1];
#ifdef __3D__  
  z = Coords[2];
#endif  
  N_Coord++; // this 1D Coord
  
  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();

  NodalPtIndex = FESpace1D->GetIntlPtIndexOfPts();
  N_GlobNodalPts = FESpace1D->GetN_RootNodalPts();
  IndexArray = new int[N_GlobNodalPts];
  memset(IndexArray, 0, SizeOfInt*N_GlobNodalPts);
  memset(val, 0, SizeOfDouble*N_GlobNodalPts);
  
  disp = 0;
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);
    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();
    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Z, AbsDetjk);

    for(j=0;j<N_Points;j++)
     {
      //       Exact(x, y, Z[j], FctVal);
      Coords[N_Coord-1] = Z[j];
      Exact(N_Coord, Coords, FctVal);
      k = NodalPtIndex[disp + j];
      val[k] += FctVal[0];
      IndexArray[k]++;
     }

    disp +=N_Points;
   } // for(i=0; i<N_Cells; i++)


   for(i=0; i<N_GlobNodalPts; i++)
    {
     if(IndexArray[i]==0) 
      {
       cout << "Error in TFEFunction1D::InterpolateNodalPts: "<< IndexArray[i] << endl;
       exit(0);
      }
     val[i] /= (double)IndexArray[i];
    }

}




/** convert current grid to vector-values FE function */
void TFEFunction1D::GridToData()
{
  int i,j,k,l;
  int N_Cells;
  int N_DOFs, N_LocalDOFs;
  int *BeginIndex, *GlobalNumbers;
  int N_Points;
  int *DOF;

  double *xi, *eta;
  double Y[MaxN_PointsForNodal1D];
  double AbsDetjk[MaxN_PointsForNodal1D];
  double PointValues[MaxN_PointsForNodal1D];
  double FunctionalValues[MaxN_PointsForNodal1D];
  double FctVal[4];

  TBaseCell *cell;
  TCollection *Coll;
  FE1D FEId;
  TFE1D *Element;
  TFE1D *FE_Obj;
  TNodalFunctional1D *nf;
  TRefTrans1D *rt;
  TBaseFunct1D *bf;
  RefTrans1D RefTrans, *RefTransArray;

  Coll = FESpace1D->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace1D->GetBeginIndex();
  GlobalNumbers = FESpace1D->GetGlobalNumbers();
  N_DOFs = FESpace1D->GetN_DegreesOfFreedom();


  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    FEId = FESpace1D->GetFE1D(i, cell);
    Element = TFEDatabase2D::GetFE1D(FEId);

    FE_Obj = TFEDatabase2D::GetFE1D(FEId);
    bf = FE_Obj->GetBaseFunct1D();
    nf = Element->GetNodalFunctional1D();
    nf->GetPointsForAll(N_Points, xi, eta);
    N_LocalDOFs = Element->GetN_DOF();

    rt = TFEDatabase2D::GetRefTrans1D(LineAffin);
    ((TLineAffin *)rt)->SetCell(cell);
    ((TLineAffin *)rt)->GetOrigFromRef(N_Points, xi, Y, AbsDetjk);

    for(j=0;j<N_Points;j++)
      PointValues[j] = Y[j];

    nf->GetAllFunctionals(PointValues, FunctionalValues);

    DOF = GlobalNumbers+BeginIndex[i];

    for(j=0;j<N_LocalDOFs;j++)
      Values[DOF[j]] = FunctionalValues[j];

   } // for(i=0; i<N_Cells; i++)

}


TFEFunction1D::~TFEFunction1D()
{
  if(Name) delete Name;
  if(Description) delete Description;
}
