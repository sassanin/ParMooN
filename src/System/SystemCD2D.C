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
   
/** ************************************************************************ 
* @brief     source file for TSystemCD2D
* @author    Sashikumaar Ganesan, 
* @date      08.08.14
* @History 
 ************************************************************************  */
#ifdef __2D__
#include <Database.h>
#include <SystemCD2D.h>
#include <SquareStructure2D.h>
#include <DiscreteForm2D.h>
#include <Assemble2D.h>
#include <AuxParam2D.h>
#include <LocalProjection.h>
#include <DirectSolver.h>
#include <Solver.h>
#include <stdlib.h>
#include <string.h>
// #include <sstream>
// #include <MooNMD_Io.h>


#include <MacroCell.h>
#include <BoundEdge.h>
#include <IsoBoundEdge.h>
#include <gridgen.h>
// #include <Remesh2D.h>
#include <IsoInterfaceJoint.h>
#include <BdLine.h>
#include <BdCircle.h>
#include <GridCell.h>

#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>


TSystemCD2D::TSystemCD2D(TFESpace2D *fespace, int disctype, int solver)
{
  //store the FEspace
  FeSpace = fespace;
  
  // N_DOF of FeSpace
  N_DOF = FeSpace->GetN_DegreesOfFreedom();
  N_Active =  FeSpace->GetActiveBound(); 
  
  //set the discretization type
  Disctype = disctype;
  
  //set the solver type
  SOLVER = solver;
  
  // build matrices
  // first build matrix structure
  sqstructure = new TSquareStructure2D(fespace);
  sqstructure->Sort();  // sort column numbers: numbers are in increasing order

  /** A is the stiffness/system mat for stationary problem   */
  sqmatrixA = new TSquareMatrix2D(sqstructure);  
  N_Matrices = 1;

}

TSystemCD2D::~TSystemCD2D()
{
//   delete sqstructure;
//   delete sqmatrixA;
}
  
  
void TSystemCD2D::Init(CoeffFct2D *BilinearCoeffs, BoundCondFunct2D *BoundCond, BoundValueFunct2D *BoundValue)
{
  BoundaryConditions[0] =  BoundCond;
  BoundaryValues[0] = BoundValue;
    
  TDiscreteForm2D *DiscreteFormUpwind;  
  TDiscreteForm2D *DiscreteFormGalerkin;
  TDiscreteForm2D *DiscreteFormSDFEM;
  TDiscreteForm2D *DiscreteFormGLS;  
  TDiscreteForm2D *DiscreteFormHeatLine;
  
  if(Disctype==HEATLINE)
   {
    InitializeDiscreteForms_HeatLine(DiscreteFormHeatLine, BilinearCoeffs);  
   }
  else
   {
    InitializeDiscreteForms_Stationary(DiscreteFormUpwind, DiscreteFormGalerkin, DiscreteFormSDFEM, DiscreteFormGLS,
                                     BilinearCoeffs);
   }
    switch(Disctype)
     {
      case GALERKIN:
      case LOCAL_PROJECTION:
           DiscreteFormARhs = DiscreteFormGalerkin;
      break;

      case SUPG:
           DiscreteFormARhs = DiscreteFormSDFEM;
      break;

      case UPWIND:
           DiscreteFormARhs = DiscreteFormUpwind;
      break;      
      
      case GLS:
           DiscreteFormARhs = DiscreteFormGLS;
      break;
      
      case HEATLINE:
           DiscreteFormARhs = DiscreteFormHeatLine;
      break;

      
      default:
            OutPut("Unknown DISCTYPE" << endl);
            exit(4711);;
     }  
     
     
} // TSystemCD2D::Init


void TSystemCD2D::Assemble(TAuxParam2D *aux, double *sol, double *rhs)
{
  int N_DirichletDof;
  double *RHSs[1];
 

  TFESpace2D *fesp[1], *ferhs[1];
   
    N_DirichletDof = N_DOF - N_Active;
    
    RHSs[0] = rhs;
    memset(rhs, 0, N_DOF*SizeOfDouble);
  
    fesp[0] = FeSpace;
    ferhs[0] = FeSpace;
    
    // initialize matrices
    SQMATRICES[0] = sqmatrixA;
    SQMATRICES[0]->Reset(); 
    
    
    if(aux==NULL)
     { aux = new TAuxParam2D(1, 0, 0, 0, fesp, NULL, NULL, NULL, NULL, 0, NULL); }
    
    // assemble
    Assemble2D(1, fesp,
               N_Matrices, SQMATRICES,
               0, NULL,
               1, RHSs, ferhs,
               DiscreteFormARhs,
               BoundaryConditions,
               BoundaryValues,
               aux);
 
     delete aux;
     
     // apply local projection stabilization method
     if(Disctype==LOCAL_PROJECTION && TDatabase::ParamDB->LP_FULL_GRADIENT>0)
      {
       if(TDatabase::ParamDB->LP_FULL_GRADIENT==1)
        { 
         UltraLocalProjection(sqmatrixA, FALSE);
        }
       else
        {
         OutPut("Check! LP_FULL_GRADIENT needs to be one to use LOCAL_PROJECTION" << endl);
         exit(4711);;
        }
      }
      
      
    // set rhs for Dirichlet nodes
    memcpy(sol+N_Active, rhs+N_Active, N_DirichletDof*SizeOfDouble);     
     
} // void TSystemCD2D::Assemble(T


void TSystemCD2D::Solve(double *sol, double *rhs)
{
    switch(SOLVER)
     {
      case AMG_SOLVE:
         Solver(sqmatrixA, rhs, sol);
      break;

      case GMG:
        cout << "GMG solver not yet implemented " <<endl;
      break;

      case DIRECT:
        DirectSolver(sqmatrixA, rhs, sol);
      break;      
 
      default:
            OutPut("Unknown Solver" << endl);
            exit(4711);;
     }    
  }

void TSystemCD2D::GetMassAndArea(TFEFunction2D *fefunction, double *parameters)
 {
  int i,j,k,l, polydegree, N_QFPoints, ORDER;
  int N_Cells, N_Joints, N_Vertices;
  int *BeginIndex, *GlobalNumbers, *DOF, N_BF;

  double *U, Mult, r_axial, val, mass, volume, Concentration;
  double *weights, *xi, *eta;
  double values[MaxN_QuadPoints_2D][MaxN_BaseFunctions2D];
  double AbsDetjk[MaxN_QuadPoints_2D], X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D];

  TJoint *joint;
  TFESpace2D *FeSpace;
  TBaseCell *cell;
  TCollection *coll;
  JointType jointtype;
  BoundTypes bdtype;
  RefTrans2D RefTrans;
  boolean IsIsoparametric;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  FE2D FEid;
  TBaseFunct2D *bf;
  TRefTrans2D *F_K;

  FeSpace = fefunction->GetFESpace2D();
  BeginIndex = FeSpace->GetBeginIndex();
  GlobalNumbers = FeSpace->GetGlobalNumbers();
  U = fefunction->GetValues();

  coll = FeSpace->GetCollection();
  N_Cells = coll->GetN_Cells();

  mass = 0.;
  volume = 0.;
  Concentration = 0.;
  for(i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    FEid = FeSpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEid);
    N_Joints = cell->GetN_Joints();
    IsIsoparametric = FALSE;
    if(TDatabase::ParamDB->USE_ISOPARAMETRIC)
     {
     for(j=0;j<N_Joints;j++)
      {
       joint = cell->GetJoint(j);
       jointtype = joint->GetType();
       if(jointtype == BoundaryEdge)
        {
         bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
         if(bdtype != Line)  IsIsoparametric = TRUE;
        }
       if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = TRUE;
        }
        if(jointtype == IsoInterfaceJoint || jointtype == IsoBoundEdge)
         IsIsoparametric = TRUE;

      } // for(j=0;j<
     } // if(TDatabase::ParamDB->USE_ISOPARAMETRIC)

   if(IsIsoparametric)
    {
      switch(N_Joints)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
    ORDER = TFEDatabase2D::GetAccuracyFromFE2D(FEid);
    switch(RefTrans)
    {
      case TriaAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaAffin *)F_K)->SetCell(cell);
//         locvol = ((TTriaAffin *)rt)->GetVolume();
        ((TTriaAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(9);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(ORDER);
        ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TTriaIsoparametric *)F_K)->GetVolume();
        ((TTriaIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadAffin:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadAffin *)F_K)->SetCell(cell);
//         locvol = ((TQuadAffin *)rt)->GetVolume();
        ((TQuadAffin *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadBilinear:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadBilinear *)F_K)->SetCell(cell);
        ((TQuadBilinear *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEid);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
	qf2 = TFEDatabase2D::GetQuadFormula2D(QuadFormula);
        qf2->GetFormulaData(N_QFPoints, weights, xi, eta);
        ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)F_K)->SetCell(cell);
//         locvol = ((TQuadIsoparametric *)rt)->GetVolume();
        ((TQuadIsoparametric *)F_K)->GetOrigFromRef(N_QFPoints, xi, eta,
                                         X, Y, AbsDetjk);
      break;
    }

    // find basis functions on cell i
    bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FEid);
    N_BF = bf->GetDimension();
    DOF = GlobalNumbers + BeginIndex[i];

    for(k=0;k<N_QFPoints;k++)
     {
      bf->GetDerivatives(D00, xi[k], eta[k], values[k]);
      
      if(TDatabase::ParamDB->Axial3D==0)
       { r_axial = 1.; }
      else
      { 
       if(TDatabase::ParamDB->Axial3DAxis==1)
        {
         r_axial = 2.*Pi*fabs(Y[k]); // r in axial3D (X: symmetric) problems (Sashikumaar Ganesan)      
        }
       else
        {
         r_axial = 2.*Pi*fabs(X[k]); // r in axial3D (Y: symmetric) problems (Sashikumaar Ganesan)      
        }
      if(X[k]<0)
       {
        cout <<"X[k] negative in Get_KE change Quad rule " <<  X[k] <<endl;
//         exit(0);
       }        
      }
      


      Mult = r_axial*weights[k]*AbsDetjk[k];
      val = 0.;
      for(l=0;l<N_BF;l++)
       {
        j = DOF[l];
        val += U[j]*values[k][l];
       }

     mass += val*Mult;
     volume += Mult;
    } //  for(k=0;k<N_QFPoints;
   } //  for(i=0;i<N_Cells;i++)

   parameters[0] = mass;
   parameters[1] = volume;
   Concentration = mass/volume;
   parameters[2] = Concentration;


   OutPut( "Time, C Mass " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[0]<< " "<<endl);
//    OutPut( "Time, C Surface area " <<TDatabase::TimeDB->CURRENTTIME<< " " <<parameters[1]<< " "<<endl);
// exit(0);
 }


#endif // #ifdef __2D__
