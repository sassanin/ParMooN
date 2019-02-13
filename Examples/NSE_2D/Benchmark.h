// Navier-Stokes problem, Benchmark problem
// 
// u(x,y) = unknown
// p(x,y) = unknown

void ExampleFile()
{
  OutPut("Example: Benchmark.h" << endl) ;
}

#define __BENCH__

#include <Joint.h>
#include <BoundEdge.h>
#include <BoundComp.h>
#include <FE2D.h>
#include <FEDesc2D.h>
// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
 switch(i)
  {
    case 0: 
       cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
       TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION=1;
//        cond = DIRICHLET;
     break;      
    case 2: 
//        cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
//        TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION=1;
       cond = DIRICHLET;
     break;
     
    case 1:
      cond = NEUMANN;
//        cond = DIRICHLET;
      TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;      
      break;
    
    case 3:       
    case 4:         
      cond = DIRICHLET;
      break;
    default: cout << "wrong boundary part number: " << i << endl;
  }   
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value=0;
            break;
    case 2: value = 0;
            break;
    case 3:
//           value=1.2*Param*(1-Param) ;
          value= Param;
//             cout<<Param<<" "<<value<<endl;
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number: " << BdComp << endl;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y,
               double **parameters, double **coeffs)
{
  double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}


/** calculate characteristic values */
void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,
             TFEFunction2D *pfct, TFEFunction2D *u1old,
             TFEFunction2D *u2old,
             double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Edges,comp,ed_nr;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D];
  double Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  int N_LocalUsedElements;
  FE2D LocalUsedElements[2], CurrentElement;
  int *DOF;
  double **OrigFEValues, *Orig;
  bool SecondDer[2] = { FALSE, FALSE };
  double *u1, *u2, *p;
  TFESpace2D *USpace, *PSpace;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  int *N_BaseFunct, N_Cells;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double value, value1, value2, value3;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValues2[MaxN_BaseFunctions2D];
  double FEFunctValues3[MaxN_BaseFunctions2D];
  int N_DerivativesU = 3;
  double *Derivatives[MaxN_BaseFunctions2D];
  MultiIndex2D NeededDerivatives[3] = { D00, D10, D01 };
  TFEFunction2D *vfct;
  double *v, nu = 1/TDatabase::ParamDB->RE_NR;
  double *Der, *aux;
  TJoint *joint;
  TBoundEdge *boundedge;
  TBoundComp *BoundComp;
  TFE2D *eleCell;
  FE2D FEEle;
  TFEDesc2D *FEDesc;
  int N_DOF_Circ, *DOF_Circ;
  char VString[] = "v";

  u1 = u1fct->GetValues();
  u2 = u2fct->GetValues();
  p = pfct->GetValues();

  USpace = u1fct->GetFESpace2D();
  PSpace = pfct->GetFESpace2D();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  aux = new double [MaxN_QuadPoints_2D*10];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*10;

  N_ = u1fct->GetLength();
  v = new double[N_];
  memset(v,0,N_*SizeOfDouble);
  vfct = new TFEFunction2D(USpace, VString, VString, v, N_);

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
 
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges=cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryEdge)||
          (joint->GetType() == IsoBoundEdge)) // boundary edge 
      {
        
        boundedge=(TBoundEdge *)joint;  
        BoundComp=boundedge->GetBoundComp();  // get boundary component
        comp=BoundComp->GetID();              // boundary id 
        if (comp==4) 
          {
            FEEle = USpace->GetFE2D(i,cell);   // finite element of cell
            eleCell =  TFEDatabase2D::GetFE2D(FEEle); 
            FEDesc = eleCell->GetFEDesc2D();   // fe descriptor
            N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
            DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
            DOF = UGlobalNumbers + UBeginIndex[i]; // pointer to global dofs
            for (k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
              v[DOF[DOF_Circ[k]]] = 1;
          }
      }      
    }
  }
  
  cd = 0;
  cl = 0;

// ########################################################################
// loop over all cells
// ########################################################################
  Coll = USpace->GetCollection(); // all spaces use same Coll
  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);

    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    N_LocalUsedElements = 2;
    LocalUsedElements[0] = USpace->GetFE2D(i, cell);
    LocalUsedElements[1] = PSpace->GetFE2D(i, cell);

    // ####################################################################
    // calculate values on original element
    // ####################################################################
    TFEDatabase2D::GetOrig(N_LocalUsedElements, LocalUsedElements, Coll, cell, SecondDer, N_Points, xi, eta, weights, X, Y, AbsDetjk);

    // calculate all needed values of p 
    CurrentElement = LocalUsedElements[1];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = PGlobalNumbers + PBeginIndex[i];
    for(l=0;l<N_;l++)
      FEFunctValues[l] = p[DOF[l]];

    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
      value = 0;
      for(l=0;l<N_;l++)
        value += FEFunctValues[l] * Orig[l];

      Derivatives[j][0] = value;
    }

    // calculate all needed values of u1, u2 
    CurrentElement = LocalUsedElements[0];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = UGlobalNumbers + UBeginIndex[i];
    for(l=0;l<N_;l++)
    {
      FEFunctValues1[l] = u1[DOF[l]];
      FEFunctValues2[l] = u2[DOF[l]];
      FEFunctValues3[l] = v[DOF[l]];
    }

    for(k=0;k<N_DerivativesU;k++)
    {
      OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct,
                                      NeededDerivatives[k]);
      for(j=0;j<N_Points;j++)
      {
        Orig = OrigFEValues[j];
        value1 = 0;
        value2 = 0;
        value3 = 0;
        for(l=0;l<N_;l++)
        {
          value1 += FEFunctValues1[l] * Orig[l];
          value2 += FEFunctValues2[l] * Orig[l];
          value3 += FEFunctValues3[l] * Orig[l];
        } // endfor l
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+4] = value2;
        Derivatives[j][k+7] = value3;
      } // endfor j
    } // endfor k

    // calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      
      // nu * (u1_x*v_x, u1_y*v_y), v= (v,0)
      value1  = nu*(Der[2]*Der[8]+Der[3]*Der[9]);
      // (u1 * u1_x + u2* u1_y) * (1,0)
      value1 += (Der[1]*Der[2]+Der[4]*Der[3])*Der[7];
      // pressure times divergence of test function (1,0)
      value1 -= Der[0]*Der[8];

      value2  = nu*(Der[5]*Der[8]+Der[6]*Der[9]);
      value2 += (Der[1]*Der[5]+Der[4]*Der[6])*Der[7];
      value2 -= Der[0]*Der[9];

      cd += AbsDetjk[j]*weights[j] * value1;
      cl += AbsDetjk[j]*weights[j] * value2;
    }

  } // endfor i

  cd *= -500;
  cl *= -500;

  delete Derivatives[0];
  delete vfct;
  delete v;
}
