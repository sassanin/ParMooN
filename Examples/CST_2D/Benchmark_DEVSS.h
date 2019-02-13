// Navier-Stokes coupled with conformation stress tensor problem, Benchmark problem


void ExampleFile()
{
  OutPut("Example: Benchmark_FlowPastCylinder_DEVSS" << endl) ;
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
      cond = DIRICHLET;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value = 0;
            break;
    case 1: value=1.2*Param*(1.0-Param);
            break;
    case 2: value = 0;
            break;
    case 3: value=1.2*Param*(1.0-Param) ;
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
  double eps = 1/TDatabase::ParamDB->RE_NR, beta = TDatabase::ParamDB->P3;
  int i;
  double *coeff, nondim;
 
  if (TDatabase::ParamDB->TENSOR_TYPE == 1)
  nondim = beta*eps;
  else if (TDatabase::ParamDB->TENSOR_TYPE == 2)
  nondim = eps;
  
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
  }
}

// ========================================================================
// exact solution
// ========================================================================
void ExactS1(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactS2(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactS3(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_CST(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;

}

void S1BoundValue(int BdComp, double Param, double &value)
{
  double Wei = TDatabase::ParamDB->WEI_NR;
  switch(BdComp)
  {
    case 0: value=1.0;
            break;
    case 1: value= (pow(Wei*1.2*(1.0-(2.0*Param))/0.41,2) * 2.0 )+1.0;
            break;
    case 2: value=1.0;
            break;
    case 3: value=(pow(Wei*1.2*(1.0-(2.0*(1.0-Param)))/0.41,2)* 2.0 )+1.0;
            break;
    case 4: value=1.0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void S2BoundValue(int BdComp, double Param, double &value)
{
  double Wei = TDatabase::ParamDB->WEI_NR;
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=Wei*1.2*(1.0-(2.0*Param))/0.41;
            break;
    case 2: value=0;
            break;
    case 3: value=Wei*1.2*(1.0-(2.0*(1.0-Param)))/0.41;
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}



void S3BoundValue(int BdComp, double Param, double &value)
{
 value=1.0;
  if(BdComp>4) cout << "wrong boundary part number: " << BdComp << endl;
}

void LinCoeffs_CST(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  double nu=1./TDatabase::ParamDB->WEI_NR;

  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
  
    coeff[0] = nu;
    coeff[1] = nu;   //f1
    coeff[2] = 0.0;  //f2
    coeff[3] = nu;   //f3
    

  }
}


// ========================================================================
// exact solution
// ========================================================================
void ExactD1(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactD2(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

void ExactD3(double x, double y, double *values)
{
    values[0] = 0;
    values[1] = 0;
    values[2] = 0;
    values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition_DFT(int i, double t, BoundCond &cond)
{
  cond = DIRICHLET;

}

void D1BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}

void D2BoundValue(int BdComp, double Param, double &value)
{
  switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0.6*(1-2*Param)/0.41;
            break;
    case 2: value=0;
            break;
    case 3: value=0.6*(1-2*(1-Param))/0.41;
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}



void D3BoundValue(int BdComp, double Param, double &value)
{
switch(BdComp)
  {
    case 0: value=0;
            break;
    case 1: value=0;
            break;
    case 2: value=0;
            break;
    case 3: value=0;
            break;
    case 4: value=0;
            break;
    default: cout << "wrong boundary part number" << endl;
            break;
  }
  return;
}


void LinCoeffs_DFT(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y, *param;
  double u1x, u1y, u2x, u2y;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];

//     param = parameters[i];
//       
//     u1x = param[2];
//     u1y = param[4];
//     u2x =  param[3];
//     u2y = param[5];

    coeff[0] = 1;
//     coeff[1] = u1x;
//     coeff[2] = (u1y+u2x)*0.5;
//     coeff[3] = u2y;
    
  }
}


/** calculate characteristic values */
void GetCdCl(TFEFunction2D *u1fct, TFEFunction2D *u2fct,
             TFEFunction2D *pfct, TFEFunction2D *tau1fct,
             TFEFunction2D *tau2fct, TFEFunction2D *tau3fct,
             double &cd, double &cl)
{
  int i,j,k,l, N_;
  int N_Points,N_Edges,comp,ed_nr;
  double *weights, *xi, *eta;
  double X[MaxN_QuadPoints_2D];
  double Y[MaxN_QuadPoints_2D];
  double AbsDetjk[MaxN_QuadPoints_2D];
  int N_LocalUsedElements;
  FE2D LocalUsedElements[3], CurrentElement;
  int *DOF;
  double **OrigFEValues, *Orig;
  bool SecondDer[3] = { FALSE, FALSE, FALSE };
  double *u1, *u2, *p, *tau1, *tau2, *tau3;
  TFESpace2D *USpace, *PSpace, *TauSpace;
  int *UGlobalNumbers, *UBeginIndex;
  int *PGlobalNumbers, *PBeginIndex;
  int *TauGlobalNumbers, *TauBeginIndex;
  int *N_BaseFunct, N_Cells;
  BaseFunct2D BaseFunct, *BaseFuncts;
  TCollection *Coll;
  TBaseCell *cell;
  double value, value1, value2, value3, value4, value5, value6;
  double FEFunctValues[MaxN_BaseFunctions2D];
  double FEFunctValues1[MaxN_BaseFunctions2D];
  double FEFunctValues2[MaxN_BaseFunctions2D];
  double FEFunctValues3[MaxN_BaseFunctions2D];
  double FEFunctValues4[MaxN_BaseFunctions2D];
  double FEFunctValues5[MaxN_BaseFunctions2D];
  double FEFunctValues6[MaxN_BaseFunctions2D];
  
  
  int N_DerivativesU = 3;
  double *Derivatives[MaxN_BaseFunctions2D];
  MultiIndex2D NeededDerivatives[3] = { D00, D10, D01 };
  TFEFunction2D *vfct;
  double *v, nu = 1/TDatabase::ParamDB->RE_NR;
  double wei = 1/TDatabase::ParamDB->WEI_NR, beta = TDatabase::ParamDB->P3;
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
  tau1 = tau1fct->GetValues();
  tau2 = tau2fct->GetValues();
  tau3 = tau3fct->GetValues();

  USpace = u1fct->GetFESpace2D();
  PSpace = pfct->GetFESpace2D();
  TauSpace = tau1fct->GetFESpace2D();

  UGlobalNumbers = USpace->GetGlobalNumbers();
  UBeginIndex = USpace->GetBeginIndex();

  PGlobalNumbers = PSpace->GetGlobalNumbers();
  PBeginIndex = PSpace->GetBeginIndex();
  
  TauGlobalNumbers = TauSpace->GetGlobalNumbers();
  TauBeginIndex = TauSpace->GetBeginIndex();

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  aux = new double [MaxN_QuadPoints_2D*13];
  for(j=0;j<MaxN_QuadPoints_2D;j++)
    Derivatives[j] = aux + j*13;

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
    N_LocalUsedElements = 3;
    LocalUsedElements[0] = USpace->GetFE2D(i, cell);
    LocalUsedElements[1] = PSpace->GetFE2D(i, cell);
    LocalUsedElements[2] = TauSpace->GetFE2D(i, cell);

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

    
    // calculate all needed values of tau1, tau2, tau3 
    CurrentElement = LocalUsedElements[2];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];

    DOF = TauGlobalNumbers + TauBeginIndex[i];
    for(l=0;l<N_;l++)
    {  
      FEFunctValues4[l] = tau1[DOF[l]];
      FEFunctValues5[l] = tau2[DOF[l]];
      FEFunctValues6[l] = tau3[DOF[l]];
      
    }
    OrigFEValues = TFEDatabase2D::GetOrigElementValues(BaseFunct, D00);

    for(j=0;j<N_Points;j++)
    {
      Orig = OrigFEValues[j];
        value4 = 0;
        value5 = 0;
        value6 = 0;
        for(l=0;l<N_;l++)
        {
          value4 += FEFunctValues4[l] * Orig[l];
          value5 += FEFunctValues5[l] * Orig[l];
          value6 += FEFunctValues6[l] * Orig[l];
        } // endfor l
      Derivatives[j][10] = value4;
      Derivatives[j][11] = value5;
      Derivatives[j][12] = value6;
    }
    
    
    
    // calculation
    for(j=0;j<N_Points;j++)
    {
      Der = Derivatives[j];
      
      // beta*nu * (u1_x*v_x, u1_y*v_y), v= (v,0)
      value1  = beta*nu*(Der[2]*Der[8]+Der[3]*Der[9]);
      // (u1 * u1_x + u2* u1_y) * (1,0)
      value1 += (Der[1]*Der[2]+Der[4]*Der[3])*Der[7];
      // pressure times divergence of test function (1,0)
      value1 -= Der[0]*Der[8];
      // (1-beta)*nu*wei* (tau1*v_x + tau2*v_y)
      value1 += (1.0-beta)*nu*wei*(Der[10]*Der[8] + Der[11]*Der[9]);

      value2  = beta*nu*(Der[5]*Der[8]+Der[6]*Der[9]);
      value2 += (Der[1]*Der[5]+Der[4]*Der[6])*Der[7];
      value2 -= Der[0]*Der[9];
      value2 += (1.0-beta)*nu*wei*(Der[11]*Der[8] + Der[12]*Der[9]);

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
