// Navier-Stokes problem, 3D Channel flow
// 
// u(x,y) = unknown
// p(x,y) = unknown
#include<BoundFace.h>

void ExampleFile()
{
#ifdef _MPI
 int rank;
 MPI_Comm_rank(TDatabase::ParamDB->Comm, &rank);

 if(rank==TDatabase::ParamDB->Par_P0)
#endif
  { 
   OutPut("Example: Aerofoil_3D_ALE.h " << endl);
  }

}

void InitialU1(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialU3(double x, double y, double z, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y,  double z, double *values)
{
  values[0] = 0;
}

// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU2(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactU3(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void ExactP(double x, double y,  double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// ========================================================================
// Description for moving grid
// ========================================================================
void GridU1(double x, double y,double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void GridU2(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

void GridU3(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}


void GridBoundCondition(int BdComp, double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void GridBoundValue(int BdComp, double x, double y, double z, double &value)
{
  value = 0;
}

void GridCoeffs(int n_points, double *x, double *y,double *z,
        double **parameters, double **coeffs)
{
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    coeff[0] = 1;

    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;

    coeff[4] = 0;
  }
}
//=============================================================================


void ModifyCoords(double x, double y, double z, 
		  double &X, double &Y, double &Z, double t)
{

 X = x;
 Y = y;
 Z = z;
  }
  
  void ModifyBdCoords(int N_MovVert, TVertex **MovBoundVert, 
		      int N_Movfaces, int * Movfaces, 
		      TBaseCell ** MovCells, double t)
{
 int i, j, k, l, m;
   
 double x, y, z, Fact=0.06;
 double disp, theta, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 

 const int *TmpFV, *TmpLen;
 int MaxLen,pt,jid;
 TBaseCell* cell;
 TBoundFace* Bdface;
 TBoundComp3D* Bdcomp;
 TVertex* vert;
   for(i=0;i<N_MovVert;i++)
    {
     MovBoundVert[i]->GetCoords(x, y, z);     
     disp = 0.05*sin(3.93*t);
//      disp = 5*sin(10*t);
     y += disp;      
     MovBoundVert[i]->SetCoords(x, y, z);    
    }
     
   
    for(i=0;i<N_Movfaces;i++)
    { 
      
    
      cell = MovCells[i];
       
      jid = Movfaces[i];
       
  

        
       cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      
     
      Bdface = (TBoundFace*)cell->GetJoint(jid);
      
      Bdcomp = Bdface->GetBoundComp();
      
      
      pt = TmpLen[jid];
      for(int ii=0;ii<pt;ii++)
      {  
	vert = cell->GetVertex(TmpFV[ jid*MaxLen + ii]);
	
	Bdcomp->GetTSofXYZ(vert->GetX(), vert->GetY(), vert->GetZ(), 
		 T[ii+1], S[ii+1]);
	 
	
      }
	
      Bdface->SetParameters(T,S);
      
    }
     
}// ModifyBdCoords



void  GetMovingBoundData(TCollection *coll, 
			 int &N_MovVert, 
			 TVertex ** &MovBoundVert,
                         int &N_Movfaces , 
			 int * &Movfaces,
			 TBaseCell ** &MovCells)
{
 int i, j, k, l, m0, m1, m2, N_Cells, comp;
 int temp;
 const int *TmpFV, *TmpLen;
 int MaxLen, N_Points;
 int bdid,pt;
 
 double  x, y, z, x1, y1, z1;
  
 TBaseCell *Me, *temp_cell;
 TJoint *Joint;
 TBoundComp *BoundComp;  
 TBoundFace* Bdface;
 TVertex *temp_vert, **Vertices;

  N_Cells = coll->GetN_Cells();
  N_MovVert = 0;
  N_Movfaces = 0;     
       for(j=0;j<N_Cells;j++)
      {
	 Me = coll->GetCell(j);
	 k = Me->GetN_Vertices();
	 for(l=0;l<k;l++)
        {
	  temp_vert = Me->GetVertex(l);
	  temp_vert->SetClipBoard(-1);
	}
      }
  
     for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
	Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
	
        k = Me->GetN_Joints();
        
	for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
         if(Joint->GetType() == BoundaryFace)
	 {
	  Bdface = (TBoundFace*)Joint;
	  BoundComp = Bdface->GetBoundComp();
	  bdid = BoundComp->GetID();
	  if(bdid == 6 || bdid == 7 )
	  {
	    pt = TmpLen[l];
	    N_Movfaces ++; // no need of clips as they belong to atmost 1 cell
	    for(int ii =0;ii<pt;ii++)
	    {
	      temp_vert = Me->GetVertex(TmpFV[ l*MaxLen + ii]);
	      if(temp_vert->GetClipBoard()==-1)
	      {
		temp_vert->SetClipBoard(1);
		N_MovVert++; 
	      }	      
	    }
	  }
	 }
	}// endfor l
      }
      
      MovBoundVert = new TVertex*[N_MovVert];
      Movfaces = new int[N_Movfaces];
      MovCells = new TBaseCell*[N_Movfaces];
      
      cout << "N_MovVert" << N_MovVert  << " N_Movfaces " << N_Movfaces << endl;
      
      N_MovVert = 0;
      N_Movfaces = 0;
      
      for(j=0;j<N_Cells;j++)
      {
        Me = coll->GetCell(j);
	Me->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
        k = Me->GetN_Joints();
        for(l=0;l<k;l++)
        {
         Joint = Me->GetJoint(l);
         if(Joint->GetType() == BoundaryFace)
	 {
	  Bdface = (TBoundFace*)Joint;
	  BoundComp = Bdface->GetBoundComp();
	  bdid = BoundComp->GetID();
	  if(bdid == 6 || bdid == 7)
	  {
	    pt = TmpLen[l];
	    Movfaces[N_Movfaces] = l;
	    if(l>4)
	    {
	      cout << "garbage " << endl;
	    }
	      
	    MovCells[N_Movfaces] = Me;
	   
	    for(int ii =0;ii<pt;ii++)
	    {
	      temp_vert = Me->GetVertex(TmpFV[ l*MaxLen + ii]);
	      if(temp_vert->GetClipBoard()==1)
	      {
		temp_vert->SetClipBoard(-1);
		MovBoundVert[N_MovVert] = temp_vert;
		N_MovVert++; 
	      }	      
	    }	    
	    N_Movfaces++;
	  }
	 }
	}// endfor l
      }
      
       cout << "N_MovVert" << N_MovVert  << " N_Movfaces " << N_Movfaces << endl;
       
}// GetMovingBoundData



// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{ 
   TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
 
     switch(CompID)
  {
      // 6 and 7 are the top and the bottom surface of the aerofoil    
    case 0: case 1: case 2: case 3:case 5:case 6:case 7:
      cond = DIRICHLET;
       break;
      
    case 4:  
     cond = NEUMANN;
      break;
      
    default: 
      cout << "wrong boundary part number" << endl;
            break;     
   
  } 
 
 
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{


   switch(CompID)
  {
            
    case 0: 
            value=0;
            break;
    case 1: 
            value=0;
            break;
    case 2: 
            value=((6-y)*(6+y)*z*(6-z))/324.00;// 12 by 12  by 6 domain sd7003_001.mesh
            break;
    case 3: 
            value=0;
            break;
    case 4: 
            value=0;
            break;
    case 5: 
            value=0;
            break;
    case 6: 
            value=0;
            break;
    case 7: 
            value=0;
            break;
    default: 
      cout << "wrong boundary part number" << endl;
            break;
  }  

 
}

// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *x, double *y, double *z,
               double **parameters, double **coeffs)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    coeff[0] = eps;
    coeff[1] = 0; // f1
    coeff[2] = 0; // f2
    coeff[3] = 0; // f3
  }
}
