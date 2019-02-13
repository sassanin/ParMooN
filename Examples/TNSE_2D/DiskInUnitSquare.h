// Navier-Stokes problem, DiskInUnitSquare FDM
//

extern "C"
{
  #include <gridgen.h>    
  void triangulate(char*, struct triangulateio*,
                   struct triangulateio*, struct triangulateio*);
}

void ExampleFile()
{
  OutPut("Example: DiskInUnitSquare.h" << endl) ;
}
// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialU2(double x, double y, double *values)
{
  values[0] = 0;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


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
   case 1:     
   case 2:    
   case 3: 
     cond = DIRICHLET;
     break; 
   case 4: 
       cond = NEUMANN;
   break;    
   default: cout << "wrong boundary part number" << endl;
       break;
  }   
     
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  double t = TDatabase::TimeDB->CURRENTTIME;

  switch(BdComp)
  {
    case 0: 
            value=0;
            break;
    case 1: 
            value=0;
            break;
    case 2:  
            value=1; // top moving side velocity
            break;
    case 3: 
            value=0;
            break;
    case 4: 
            value=0;
            break;   
    default: cout << "wrong boundary part number" << endl;
            break;
  }
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  int i;
  double *coeff, x, y; 
  static double eps=1/TDatabase::ParamDB->RE_NR;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = eps;

    coeff[1] = 0;  // f1
    coeff[2] = 0;  // f2
  }
}


// kind of boundary condition (for FE space needed)
void LevelSetBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = DIRICHLET;
}




// value of boundary condition
void LevelSetBoundValue(int BdComp, double Param, double &value)
{
  value = 0.;
}

// ========================================================================
// initial solution
// ========================================================================
void InitialLevelSet(double x, double y, double *values)
{
  values[0] = 0.25 - sqrt(x*x + y*y);
}

// kind of boundary condition (for FE space needed)
void LagrangeBoundCondition(int BdComp, double t, BoundCond &cond)
{
  cond = NEUMANN;
}

// value of boundary condition
void LagrangeBoundValue(int BdComp, double Param, double &value)
{
  value = 0.;
}


// // ========================================================================
// // coefficients for Stokes form: A, B1, B2, f1, f2
// // ========================================================================
// void LagrangeLinCoeffs(int n_points, double *X, double *Y,
//                double **parameters, double **coeffs)
// {
//   int i;
//   double *coeff, x, y; 
//   static double eps=1/TDatabase::ParamDB->RE_NR;
// 
//   for(i=0;i<n_points;i++)
//   {
//     coeff = coeffs[i];
//     x = X[i];
//     y = Y[i];
// 
//     coeff[0] = eps;
// 
//     coeff[1] = 0;  // f1
//     coeff[2] = 0;  // f2
//   }
// }









// **************************************************************************
// Triangular Mesh Generation
// **************************************************************************

void  TriaReMeshGen(TDomain *&Domain)
{
  int j, ID, k, N_G, *PartMarker, *PointNeighb, maxEpV=0;
  int a, b, len1, len2, Neighb_tmp, BDpart;
  int i, temp, N_Cells, N_RootCells, CurrVertex, N_Joints, N_Vertices;
  int N_Interface_Vert, N_Verti, N_Hori, N_SlipBound_Vert, N_BDVertices;
  int CurrComp, In_Index, *Triangles, Neib[2], CurrNeib;
  
  double deviation, hi, x0, y0, x, y, phi1, phi2;
  double T_a, T_b, C_x, C_y, s, theta;  
  double dt, area = TDatabase::ParamDB->Area;
  double *Coordinates, *Hole_List, left, right, bottom, top;  
  
  double Xi[4] = {-2., 2., 2.,-2.};
  double Yi[4] = {-5.,-5., 2., 2.}; 

  TBaseCell **CellTree, *cell;
  TBoundPart *BoundPart;
  TJoint *Joint;
  TCollection *coll;
  TVertex **VertexDel, **NewVertices;  
  TBdLine *UpdateBound[12];
  TBdCircle *UpdateIntface;

  
  boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
  
  struct triangulateio In, Out;
  std::ostringstream opts;
  opts << " ";
  
  Out.pointlist = NULL;
  Out.pointattributelist = NULL;
  Out.pointmarkerlist = NULL;
  Out.trianglelist = NULL;
  Out.triangleattributelist = NULL;
  Out.trianglearealist = NULL;
  Out.neighborlist = NULL;
  Out.segmentlist = NULL;
  Out.segmentmarkerlist = NULL;
  Out.holelist = NULL;
  Out.regionlist = NULL;
  Out.edgelist = NULL;
  Out.edgemarkerlist = NULL;
  Out.normlist = NULL;

  opts.seekp(std::ios::beg); 
  

  BoundPart = Domain->GetBdPart(0);
  UpdateBound[0]  = (TBdLine*)BoundPart->GetBdComp(0);
  UpdateBound[1]  = (TBdLine*)BoundPart->GetBdComp(1);
  UpdateBound[2]  = (TBdLine*)BoundPart->GetBdComp(2);
  UpdateBound[3]  = (TBdLine*)BoundPart->GetBdComp(3);
  BoundPart = Domain->GetBdPart(1);
  UpdateIntface = (TBdCircle*)BoundPart->GetBdComp(0);
  
  
//OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);

  opts<<'p'; // Constrained Delaunay Triangulation:
           // initial values - only points defined on the boundary of the domain;
           // triangulation near boundary may variate from Delaunay criterion
  opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
              // Quality mesh generation with no angles smaller than 20 degrees;
  opts<<"a"<< area; // Imposes a maximum triangle area.
  opts<<'e'; // Outputs a list of edges of the triangulation
  opts<<'z'; // Numbers if items starting from 0
//   opts<<"VVVV"; // Gives detailed information about what Triangle is doing
  opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
//   opts<<'Y'; // Supress adding vertices on boundary edges
  opts<<ends;
  
  N_Interface_Vert = 50;    //Freesurf except end point
  N_Hori  = 15;      // number of horrizontal BD vertices
  N_Verti = 50;       // number of vertical BD vertices
  N_SlipBound_Vert = 2*N_Hori + 2*N_Verti;  

  N_BDVertices = N_Interface_Vert+N_SlipBound_Vert;
  In.numberofpoints = N_BDVertices;
  In.pointlist = new double[2*In.numberofpoints];
  In.pointmarkerlist = new int[In.numberofpoints];
  In.numberofpointattributes = 0;

  In.numberofsegments = In.numberofpoints;
  In.segmentlist = new int[2*In.numberofsegments];
  In.segmentmarkerlist = new int[In.numberofsegments]; 
  In.numberofregions = 0;
  In.regionlist = NULL;
  
  In.numberofholes = 0;
//   Hole_List = NULL;
//   Hole_List[0] = 0.;
//   Hole_List[1] = 0.;
  In.holelist = NULL;
    
  In_Index = 0;
  CurrComp = 1;  
  
  hi = (Xi[1] - Xi[0])/(double)N_Hori;
  x0 = Xi[0];
  y0 = Yi[0];
  x  = Xi[0];  
  
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Hori;i++) // without last point
   {
    x = x0 + (double)i*hi;
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
   }
   
  CurrComp++;  
//    cout<<endl; 
  hi = (Yi[2] - Yi[1])/(double)N_Verti;
  x0 = Xi[1];
  y0 = Yi[1];
  y  = Yi[1];
  // points and segments on the horizontal boundary (marker=1)
  for(i=0;i<N_Verti;i++) // without last point
   {
    y = y0 + (double)i*hi; 
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;

   }
   
  CurrComp++;  
//   cout<<endl;

  hi = (Xi[3] - Xi[2])/(double)N_Hori;
  x0 = Xi[2];
  y0 = Yi[2];
  x  = Xi[2];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Hori;i++) // without last point
   {
    x = x0 + (double)i*hi;
    In.pointlist[2*In_Index] = x;
    In.pointlist[2*In_Index+1] = y0;
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;

   }
   
  CurrComp++;
//   cout<<endl;

  hi = (Yi[0] - Yi[3])/(double)N_Verti;
  x0 = Xi[3];
  y0 = Yi[3];
  y  = Yi[3];
  // points and segments on the horizontal boundary (marker=1)
 for(i=0;i<N_Verti;i++) // without last point
   {
    y = y0 + (double)i*hi;     
    In.pointlist[2*In_Index] = x0;
    In.pointlist[2*In_Index+1] = y;
//     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
    In.pointmarkerlist[In_Index] = CurrComp;
    In.segmentlist[2*In_Index] = In_Index;
    In.segmentlist[2*In_Index+1] = In_Index+1;
    In.segmentmarkerlist[In_Index] = CurrComp;
    In_Index++;
 
   }
  CurrComp++;  
  
  In.segmentlist[2*(In_Index-1)+1] = 0;
  temp=In_Index;
//    cout<<end/*l*/;
  
//   T_a = TDatabase::ParamDB->P1; // x axis value in ellipse
//   T_b = TDatabase::ParamDB->P2; // y axis value in ellipse !!! set in modifyGausscoord funct also
  T_a = 0.5;
  T_b = 0.5;
//   deviation = fabs( T_a - T_b);  


   C_x = 0.0;  // center of the inner phase
   C_y = 0.0, // center of the inner phase
   phi1 = 0.000000E+0000; // end degree value of interface
   phi2 = 2.*Pi;
 
//   C_x = 0.0;  // center of the inner phase
//   C_y = 0.0, // center of the inner phase
  s = (phi2- phi1)/(double)N_Interface_Vert; 
 
  // points and segments on the interface (marker=2)
//   theta = 0.;
//   t0 = theta;
//   cout << " s " << s << endl;
   for(i=0;i<N_Interface_Vert;i++)
    {
      theta = phi1 + (double)i*s;       
//      cout<<" theta : "<< theta <<endl;
      
//       if(fabs(I_FaceX[i]) < 1e-10 ) I_FaceX[i] = 0.0;
      In.pointlist[2*In_Index] =   T_a*cos(theta);;
      In.pointlist[2*In_Index+1] =  T_b*sin(theta);

//       if(i==0) 
//        {
//         FreeBD_X = I_FaceX[i];   FreeBD_Y = I_FaceY[i]; 
//        cout << " sorting " << FreeBD_X << ' ' << FreeBD_Y<<endl;
//        }

//       if(i==1)
//        {
//         h_interface = sqrt((I_FaceX[i-1]-I_FaceX[i])*(I_FaceX[i-1]-I_FaceX[i]) +
//                             (I_FaceY[i-1]-I_FaceY[i])*(I_FaceY[i-1]-I_FaceY[i]) );
//         OutPut("h_interface " <<h_interface << endl);
//        }
//       cout<<(180./Pi)*theta<< " x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;

      In.pointmarkerlist[In_Index] = CurrComp;
      In.segmentlist[2*In_Index] = In_Index;
      In.segmentlist[2*In_Index+1] = In_Index+1;
      In.segmentmarkerlist[In_Index] = CurrComp;
      In_Index++;
    }

  In.segmentlist[2*(In_Index-1)+1] = temp;  
   

 
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) {
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;}
    
    
  // call triangle
  triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);
  

 
  Domain->GetTreeInfo(CellTree, N_RootCells);
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  
  // remove all existing vertices and joints
  VertexDel = new TVertex*[3*N_RootCells];
  
   CurrVertex = 0;

  for(i=0;i<N_Cells;i++)
    {
      cell = coll->GetCell(i);
      N_Joints = cell->GetN_Joints();
      N_Vertices = cell->GetN_Vertices();
      for(j=0;j<N_Joints;j++)
        {
         if(CurrVertex==0)
          {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
           }
          else
           {
            ID = 0;
            for(k=0;k<CurrVertex;k++)
            if(VertexDel[k]==cell->GetVertex(j))
             {
              ID = 1; break;
             }
            if(ID!=1)
             {
              VertexDel[CurrVertex] = cell->GetVertex(j);
              CurrVertex++;
             }
           } // else if(CurrVertex==0)

          ID = 0;
          for(k=0;k<CurrVertex;k++)
          if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
           {
            ID = 1; break;
           }
            if(ID!=1)
           {
            VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
            CurrVertex++;
           }
        } // for j
    } // for i
  for(i=0;i<CurrVertex;i++)
    delete VertexDel[i];

  delete []VertexDel;
  OutPut(CurrVertex<<" vertices were deleted"<<endl);

 // remove all existing cells and joints
  for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
  OutPut(N_RootCells<<" cells were deleted"<<endl);
   delete CellTree;
   delete coll;

       
//  exit(0);
   
 // Solid Bound startx, starty, x length and y length
  UpdateBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
  UpdateBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
  UpdateBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
  UpdateBound[3]->SetParams(Xi[3], Yi[3], Xi[0]-Xi[3],Yi[0]-Yi[3]);
  

// Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
  UpdateIntface->SetParams(C_x, C_y, T_a, T_b, phi1, phi2);   
   
   
  N_RootCells = Out.numberoftriangles; 

  // allocate auxillary fields
  Coordinates = Out.pointlist;
  Triangles = Out.trianglelist;
  PartMarker = new int[Out.numberofpoints];

  // generate new vertices
  N_G = Out.numberofpoints;
  NewVertices = new TVertex*[N_G];

  for (i=0;i<N_G;i++)
     NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);

      // set bounding box
  left = bottom = 1e8;
  right = top = -1e8;

   for(i=0;i<In.numberofpoints;i++)
    {
      if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
      if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
      if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
      if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
    }  
 
//   OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);

  Domain->SetBoundBox(right-left,top-bottom);
  Domain->SetBoundBoxstart(left,bottom);
  
 // generate cells
  CellTree = new TBaseCell*[N_RootCells];

  for (i=0;i<N_RootCells;i++)
  {
    CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);

    CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
    CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
    CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);

      ((TMacroCell *) CellTree[i])->SetSubGridID(0);
  }

  Domain->SetTreeInfo(CellTree, N_RootCells);

  // initialize iterators
  TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
  TDatabase::IteratorDB[It_LE]->SetParam(Domain);
  TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
  TDatabase::IteratorDB[It_Between]->SetParam(Domain);
  TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);   
  
  // search neighbours
  N_G = Out.numberofpoints;
  PointNeighb = new int[N_G];

  memset(PointNeighb, 0, N_G *SizeOfInt);

  for (i=0;i<3*N_RootCells;i++)
    PointNeighb[Triangles[i]]++;

  maxEpV=0;

  for (i=0;i<N_G;i++)
    if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
  delete [] PointNeighb;

  PointNeighb = new int[++maxEpV * N_G];

  memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);
      
   // first colomn contains the number of following elements
   // for every point at first column we set the number of neighbour points
   // at further columns we set the index of corresponding cells
  for(i=0;i<3*N_RootCells;i++)
  {
    j = Triangles[i]*maxEpV;
    PointNeighb[j]++;
    PointNeighb[j + PointNeighb[j]] = i / 3;
  }
 
  // generate new edges
  N_G = Out.numberofedges;
  for (i=0;i<N_G;i++)
  {
    a = Out.edgelist[2*i];
    b = Out.edgelist[2*i+1];
    Neib[0] = -1;
    Neib[1] = -1;
    CurrNeib = 0;

    len1 = PointNeighb[a*maxEpV];
    len2 = PointNeighb[b*maxEpV];

  // find indexes of cells containing the current edge
   for (j=1;j<=len1;j++)
    {
      Neighb_tmp = PointNeighb[a*maxEpV + j];
       for (k=1;k<=len2;k++)
        if (Neighb_tmp == PointNeighb[b*maxEpV + k])
        {
          Neib[CurrNeib++] = Neighb_tmp;
          break;
        }
      if (CurrNeib == 2) break;
    }

//          cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
 
 
    
    if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
    {
      CurrComp = Out.edgemarkerlist[i] - 1;
      if (CurrComp >= 100000) CurrComp -= 100000;
      
      BDpart=Domain->GetBdPartID(CurrComp);
      CurrComp= Domain->GetLocalBdCompID(CurrComp);

//      cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
 
      if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
          Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->GetTofXY(
            NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
       {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
        //  exit(0);
       }
      
//  exit(0);   
// ===============================================================================       
 //  due to the hole the orientation of the circle is colck-wise  from out side
// the parameter of the starting edge (in 0, 2pi) is for e.g (0, 0.9) (colck-wise) 
// while refining to get the mid point we should get  the mid point parameter as 0.95 
// but we will get only (0+0.9)/2 =0.45 (wrong mid point). So we change.
//  Note only if the orientation is colck-wise  !!!!!!!
// ===============================================================================  

     if(BDpart==1 && CurrComp==0  && fabs(T_a)==0 ) T_a=1;      
       
       
//       cout<<BDpart << " BDpart CurrComp "<< CurrComp <<endl;
     
     
      if (CurrNeib == 2)    // 2 cells contain the current edge
       {
        if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->IsFreeBoundary())
	 {
          Joint = new TIsoInterfaceJoint(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
	 }
        else
	 {
          Joint = new TInterfaceJoint(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp),
                  T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
	 }
       }
      else
       {
        if(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp)->IsFreeBoundary())
	 {
          Joint = new TIsoBoundEdge(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp), T_a, T_b);
// 	   cout<<" FreeBoundary"<<endl;
	 }
        else
         { 
          Joint = new TBoundEdge(Domain->GetBdPart(BDpart)->GetBdComp(CurrComp), T_a, T_b);
         }
       }
    }
    else // inner edge
    {
    if (CurrNeib != 2)
        cerr << "Error!!!!!!!! not enough neighbours!" << endl;

    Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
    }

    // find the local index for the point 'a' on the cell
    for (j=0;j<3;j++)
      if (Triangles[3*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<3;k++)
      if (Triangles[3*Neib[0]+k] == b) break;

   k = k*10 + j;

    switch (k)
    {
      case  1:
      case 10:
        j = 0;
        break;
      case 12:
      case 21:
        j = 1;
        break;
      case  2:
      case 20:
        j = 2;
        break;
    }

  CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1)
    {
      // find the local index for the point 'a' on the cell
      for (j=0;j<3;j++)
        if (Triangles[3*Neib[1]+j] == a) break;

      // find the local index for the point 'b' on the cell
      for (k=0;k<3;k++)
        if (Triangles[3*Neib[1]+k] == b) break;

      k = k*10 + j;

      switch (k) // j will contain the local index for the current
      {
        case  1:
        case 10:
          j = 0;
          break;
        case 12:
        case 21:
          j = 1;
          break;
        case  2:
        case 20:
          j = 2;
          break;
      }

      CellTree[Neib[1]]->SetJoint(j, Joint);
    }

  if (Joint->GetType() == InterfaceJoint ||
        Joint->GetType() == IsoInterfaceJoint)
      ((TInterfaceJoint *) Joint)->CheckOrientation();
  }
 

  delete [] NewVertices;
  delete [] PointNeighb;
  delete [] In.pointlist;
  delete [] In.pointmarkerlist;
  delete [] In.segmentlist;
  delete [] In.segmentmarkerlist;
  
  if(Out.pointlist!=NULL) {
    free(Out.pointlist); Out.pointlist = NULL;}
  if(Out.pointattributelist!=NULL) { 
    free(Out.pointattributelist); Out.pointattributelist = NULL;}
  if(Out.pointmarkerlist!=NULL) {
    free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
  if(Out.trianglelist!=NULL) {
    free(Out.trianglelist); Out.trianglelist = NULL;}
  if(Out.triangleattributelist!=NULL) {
    free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
  if(Out.trianglearealist!=NULL) {
    free(Out.trianglearealist); Out.trianglearealist = NULL;}
  if(Out.neighborlist!=NULL) {
    free(Out.neighborlist); Out.neighborlist = NULL;}
  if(Out.segmentlist!=NULL) {
    free(Out.segmentlist); Out.segmentlist = NULL;}
  if(Out.segmentmarkerlist!=NULL) {
    free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
  if(Out.holelist!=NULL) {
    free(Out.holelist); Out.holelist = NULL;}
  if(Out.regionlist!=NULL) {
    free(Out.regionlist); Out.regionlist = NULL;}
  if(Out.edgelist!=NULL) {
    free(Out.edgelist); Out.edgelist = NULL;}
  if(Out.edgemarkerlist!=NULL) {
    free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
  if(Out.normlist!=NULL) {
    free(Out.normlist); Out.normlist = NULL;} 
  
//======================================================================
// Triangular for grid generation --end
//======================================================================


//   cout<< "tetgen" <<endl;
//   exit(0);
  
 
} // TriaReMeshGen


// void  TriaReMeshGen(TDomain *&Domain, int &N_LineBDComp, int &N_Cells_P1, int &N_Cells_P2)
// {
//  int rank, size;
// 
// 
//   MPI_Comm_rank(comm, &rank);
//   MPI_Comm_size(comm, &size);
// 
//   TBoundComp *BoundComp;
// //   TBdLine *UpdateSlipBound[4];
//   TBdCircle *UpdateIntface;
//   TBaseCell **CellTree, *cell;
//   TBoundPart *BoundPart;
//   TJoint *Joint;
//   TCollection *coll;
//   TVertex **VertexDel, **NewVertices;
// 
//   double dt, area = TDatabase::ParamDB->Area;
//   double mode_diff = TDatabase::ParamDB->P4;
//   double r, phi, t, T_a, T_b, x, y, theta;
//   double *Coordinates, temp;
//   double left, right, bottom, top;
//   double y_begin, y_end, dx, deviation;
//   double hi, x0, y0, C_x, C_y, phi1, phi2;
//   double Xi[4] = {0., 1., 1., 0.};
//   double Yi[4] = {0., 0., 5., 5.};
//   double *I_FaceX, *I_FaceY, h_interface;
// 
//   int i, j, k, ID, In_Index, N_Cells, N_G, *PointNeighb, maxEpV=0;
//   int a, b, len1, len2, Neighb_tmp, mode = int (TDatabase::ParamDB->P4);
//   int N_Interf_Vertices, CurrComp, CurrVertex, N_Joints, N_Vertices;
//   int N_RootCells, *PartMarker, *Triangles, Neib[2], CurrNeib;
//   int N_SlipBound_Vert, N_Hori, N_Verti, temp_segment;
//   int N_Interface_Vert, N_Old_Face_Vert, N_P;
// 
//   boolean AllowEdgeRef = (boolean) TDatabase::ParamDB->MESHGEN_ALLOW_EDGE_REF;
//  // cout<< " MESHGEN_ALLOW_EDGE_REF" <<AllowEdgeRef<<endl;
// 
//   struct triangulateio In, Out;
//   std::ostringstream opts;
//   opts << " ";
// 
//   deviation = fabs(TDatabase::ParamDB->P1 -  TDatabase::ParamDB->P2);
//   BoundPart = Domain->GetBdPart(0);
// 
//   N_LineBDComp = 4;
//   UpdateSlipBound = new TBdLine *[4];
// 
//   UpdateSlipBound[0]  = (TBdLine*)BoundPart->GetBdComp(0);
//   UpdateSlipBound[1]  = (TBdLine*)BoundPart->GetBdComp(1);
//   UpdateSlipBound[2]  = (TBdLine*)BoundPart->GetBdComp(2);
//   UpdateSlipBound[3]  = (TBdLine*)BoundPart->GetBdComp(3);
// 
//   UpdateIntface = (TBdCircle*)BoundPart->GetBdComp(4);
// 
//   Out.pointlist = NULL;
//   Out.pointattributelist = NULL;
//   Out.pointmarkerlist = NULL;
//   Out.trianglelist = NULL;
//   Out.triangleattributelist = NULL;
//   Out.trianglearealist = NULL;
//   Out.neighborlist = NULL;
//   Out.segmentlist = NULL;
//   Out.segmentmarkerlist = NULL;
//   Out.holelist = NULL;
//   Out.regionlist = NULL;
//   Out.edgelist = NULL;
//   Out.edgemarkerlist = NULL;
//   Out.normlist = NULL;
// 
//   opts.seekp(std::ios::beg);
// 
// 
// 
// //OutPut("MESHGEN_REF_QUALIT " << TDatabase::ParamDB->MESHGEN_REF_QUALITY << endl);
// 
//  opts<<'p'; // Constrained Delaunay Triangulation:
//            // initial values - only points defined on the boundary of the domain;
//            // triangulation near boundary may variate from Delaunay criterion
//  opts<<"q"<<  TDatabase::ParamDB->MESHGEN_REF_QUALITY;
//               // Quality mesh generation with no angles smaller than 20 degrees;
// 
//   opts<<"a"<< area; // Imposes a maximum triangle area.
//   opts<<'e'; // Outputs a list of edges of the triangulation
//   opts<<'z'; // Numbers if items starting from 0
//   //opts<<"VVVV"; // Gives detailed information about what Triangle is doing
//   opts<<'Q'; // Supress all explanation of what Triangle is doing, unless an error occurs
// //   opts<<'Y'; // Supress adding vertices on boundary edges
//   opts<<ends;
// 
//   N_Interface_Vert = int (TDatabase::ParamDB->P6);    //Freesurf except end point
//   N_Hori  = 20;      // number of horoyontal vertices
//   N_Verti = 60;       // number of horoyontal vertices
//   N_SlipBound_Vert = 2*N_Hori + 2*N_Verti;
// 
//   N_Interf_Vertices = N_Interface_Vert+N_SlipBound_Vert;
//   In.numberofpoints = N_Interf_Vertices;
//   In.pointlist = new double[2*In.numberofpoints];
//   In.pointmarkerlist = new int[In.numberofpoints];
//   In.numberofpointattributes = 0;
// 
//   In.numberofsegments = In.numberofpoints;
//   In.segmentlist = new int[2*In.numberofsegments];
//   In.segmentmarkerlist = new int[In.numberofsegments];
//   In.numberofholes = 0;
//   In.holelist = NULL;
//   In.numberofregions = 0;
//   In.regionlist = NULL;
// 
//   In_Index = 0;
//   CurrComp = 1;
// 
//   hi = (Xi[1] - Xi[0])/N_Hori;
//   x0 = Xi[0];
//   y0 = Yi[0];
//   x  = Xi[0];
//   // points and segments on the horizontal boundary (marker=1)
//  for(i=0;i<N_Hori;i++) // without last point
//    {
//     In.pointlist[2*In_Index] = x;
//     In.pointlist[2*In_Index+1] = y0;
//     //     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
//     In.pointmarkerlist[In_Index] = CurrComp;
//     In.segmentlist[2*In_Index] = In_Index;
//     In.segmentlist[2*In_Index+1] = In_Index+1;
//     In.segmentmarkerlist[In_Index] = CurrComp;
//     In_Index++;
//     x = x0 + (i+1)*hi;
//    }
//   CurrComp++;
// 
//   hi = (Yi[2] - Yi[1])/N_Verti;
//   x0 = Xi[1];
//   y0 = Yi[1];
//   y  = Yi[1];
//   // points and segments on the horizontal boundary (marker=1)
//  for(i=0;i<N_Verti;i++) // without last point
//    {
//     In.pointlist[2*In_Index] = x0;
//     In.pointlist[2*In_Index+1] = y;
//     //     cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
//     In.pointmarkerlist[In_Index] = CurrComp;
//     In.segmentlist[2*In_Index] = In_Index;
//     In.segmentlist[2*In_Index+1] = In_Index+1;
//     In.segmentmarkerlist[In_Index] = CurrComp;
//     In_Index++;
//     y = y0 + (i+1)*hi;
//    }
//   CurrComp++;
// 
// 
//   hi = (Xi[3] - Xi[2])/N_Hori;
//   x0 = Xi[2];
//   y0 = Yi[2];
//   x  = Xi[2];
//   // points and segments on the horizontal boundary (marker=1)
//  for(i=0;i<N_Hori;i++) // without last point
//    {
//     In.pointlist[2*In_Index] = x;
//     In.pointlist[2*In_Index+1] = y0;
// // cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
//     In.pointmarkerlist[In_Index] = CurrComp;
//     In.segmentlist[2*In_Index] = In_Index;
//     In.segmentlist[2*In_Index+1] = In_Index+1;
//     In.segmentmarkerlist[In_Index] = CurrComp;
//     In_Index++;
//     x = x0 + (i+1)*hi;
//    }
//   CurrComp++;
// 
// 
//   hi = (Yi[0] - Yi[3])/N_Verti;
//   x0 = Xi[3];
//   y0 = Yi[3];
//   y  = Yi[3];
//   // points and segments on the horizontal boundary (marker=1)
//  for(i=0;i<N_Verti;i++) // without last point
//    {
//     In.pointlist[2*In_Index] = x0;
//     In.pointlist[2*In_Index+1] = y;
//     // cout<<" x : "<< In.pointlist[2*In_Index] << " y : "<< In.pointlist[2*In_Index+1] <<endl;
//     In.pointmarkerlist[In_Index] = CurrComp;
//     In.segmentlist[2*In_Index] = In_Index;
//     In.segmentlist[2*In_Index+1] = In_Index+1;
//     In.segmentmarkerlist[In_Index] = CurrComp;
//     In_Index++;
//     y = y0 + (i+1)*hi;
//    }
//   CurrComp++;
// 
//   In.segmentlist[2*(In_Index-1)+1] = 0;
//   temp_segment=In_Index;
// 
// 
//  T_a = 0.25; // x axis value in ellipse
//  T_b = 0.25; // y axis value in ellipse !!! set in modifyGausscoord funct also
// 
//  C_x = .50;  // center of the inner phase
//  C_y = .50, // center of the inner phase
// 
//  phi1 = 0.000000E+0000; // end degree value of interface
//  phi2 = 6.28318530717958647688E+0000;
//  t = (phi2-phi1)/N_Interface_Vert;
// 
//   N_Old_Face_Vert = N_Interface_Vert;
//   I_FaceX = new double[N_Interface_Vert]; // for spline construction
//   I_FaceY = new double[N_Interface_Vert]; // for spline construction
// 
// 
//  // points and segments on the interface (marker=2)
//   theta = phi1;
//   double t0 = theta;
//   for(i=0;i<N_Interface_Vert;i++)
//     {
// //      cout<<" theta : "<< theta <<endl;
//       I_FaceX[i] = C_x + T_a*cos(theta);
//       I_FaceY[i] = C_y + T_b*sin(theta);
//       
//       if(fabs(I_FaceX[i]) < 1e-10 ) I_FaceX[i] = 0.0;
//       In.pointlist[2*In_Index] = I_FaceX[i];
//       In.pointlist[2*In_Index+1] = I_FaceY[i];
// 
//       if(i==0) 
//        {
//         FreeBD_X = I_FaceX[i];   FreeBD_Y = I_FaceY[i]; 
//        cout << " sorting " << FreeBD_X << ' ' << FreeBD_Y<<endl;
//        }
// 
//       if(i==1)
//        {
//         h_interface = sqrt((I_FaceX[i-1]-I_FaceX[i])*(I_FaceX[i-1]-I_FaceX[i]) +
//                             (I_FaceY[i-1]-I_FaceY[i])*(I_FaceY[i-1]-I_FaceY[i]) );
//         cout << "h_interface " <<h_interface << endl;
//        }
// //       cout<<(180./Pi)*theta<< " x : "<< I_FaceX[i] << " y : "<< I_FaceY[i] <<endl;
// 
//       In.pointmarkerlist[In_Index] = CurrComp;
//       In.segmentlist[2*In_Index] = In_Index;
//       In.segmentlist[2*In_Index+1] = In_Index+1;
//       In.segmentmarkerlist[In_Index] = CurrComp;
//       In_Index++;
//       theta = t0 + double(i+1)*t; 
//     }
// 
//   In.segmentlist[2*(In_Index-1)+1] = temp_segment;
// 
// // exit(0);
//   if(Out.pointlist!=NULL) {
//     free(Out.pointlist); Out.pointlist = NULL;}
//   if(Out.pointattributelist!=NULL) {
//     free(Out.pointattributelist); Out.pointattributelist = NULL;}
//   if(Out.pointmarkerlist!=NULL) {
//     free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
//   if(Out.trianglelist!=NULL) {
//     free(Out.trianglelist); Out.trianglelist = NULL;}
//   if(Out.triangleattributelist!=NULL) {
//     free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
//   if(Out.trianglearealist!=NULL) {
//     free(Out.trianglearealist); Out.trianglearealist = NULL;}
//   if(Out.neighborlist!=NULL) {
//     free(Out.neighborlist); Out.neighborlist = NULL;}
//   if(Out.segmentlist!=NULL) {
//     free(Out.segmentlist); Out.segmentlist = NULL;}
//   if(Out.segmentmarkerlist!=NULL) {
//     free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
//   if(Out.holelist!=NULL) {
//     free(Out.holelist); Out.holelist = NULL;}
//   if(Out.regionlist!=NULL) {
//     free(Out.regionlist); Out.regionlist = NULL;}
//   if(Out.edgelist!=NULL) {
//     free(Out.edgelist); Out.edgelist = NULL;}
//   if(Out.edgemarkerlist!=NULL) {
//     free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
//   if(Out.normlist!=NULL) {
//     free(Out.normlist); Out.normlist = NULL;}
// 
// 
// if(rank==0)
//  {
// /*
//   for(i=0;i<In.numberofpoints;i++)
//     OutPut(i<<' '<<In.pointmarkerlist[i]<<' '<<
// 	   In.pointlist[2*i]<<' '<<In.pointlist[2*i+1]<<endl);
// cout<<endl;
// //exit(0);
// */
//   triangulate((char*)opts.str().c_str(), &In, &Out, (struct triangulateio *)NULL);
// 
// /*
//   for(i=0;i<Out.numberofpoints;i++)
//      OutPut(i<<' '<<Out.pointmarkerlist[i]<<' '<<
//       Out.pointlist[2*i]<<' '<<Out.pointlist[2*i+1]<<endl);
//   */
//  }
// 
//   MPI_Bcast(&Out.numberofpoints, 1, MPI_INT, 0, comm);
//   MPI_Bcast(&Out.numberoftriangles, 1, MPI_INT, 0, comm);
//   MPI_Bcast(&Out.numberofedges, 1, MPI_INT, 0, comm);
// 
// 
//   if(rank!=0)
//    {
//     Out.pointlist = new double[2*Out.numberofpoints];
//     Out.trianglelist = new int[3*Out.numberoftriangles];
//     Out.edgelist = new int[2*Out.numberofedges];
//     Out.edgemarkerlist = new int[Out.numberofedges];
//    }
// 
//    MPI_Bcast(Out.pointlist, 2*Out.numberofpoints, MPI_DOUBLE, 0, comm);
//    MPI_Bcast(Out.trianglelist, 3*Out.numberoftriangles, MPI_INT, 0, comm);
//    MPI_Bcast(Out.edgelist, 2*Out.numberofedges, MPI_INT, 0, comm);
//    MPI_Bcast(Out.edgemarkerlist, Out.numberofedges, MPI_INT, 0, comm);
// 
// //         cout << "Triangulation end\n" << endl;
// //         MPI_Finalize();
// //         exit(0);
// 
//   Domain->GetTreeInfo(CellTree,N_RootCells);
//   coll = Domain->GetCollection(It_Finest, 0);
//   N_Cells = coll->GetN_Cells();
// 
//   // remove all existing vertices and joints
//   VertexDel = new TVertex*[3*N_RootCells];
//   // DelCell = new TGridCell*[N_Cells];
// 
//   CurrVertex = 0;
// 
//   for(i=0;i<N_Cells;i++)
//     {
//       cell = coll->GetCell(i);
//       N_Joints = cell->GetN_Joints();
//       N_Vertices = cell->GetN_Vertices();
//       for(j=0;j<N_Joints;j++)
//         {
//          if(CurrVertex==0)
//           {
//               VertexDel[CurrVertex] = cell->GetVertex(j);
//               CurrVertex++;
//            }
//           else
//            {
//             ID = 0;
//             for(k=0;k<CurrVertex;k++)
//             if(VertexDel[k]==cell->GetVertex(j))
//              {
//               ID = 1; break;
//              }
//             if(ID!=1)
//              {
//               VertexDel[CurrVertex] = cell->GetVertex(j);
//               CurrVertex++;
//              }
//            } // else if(CurrVertex==0)
// 
//           ID = 0;
//           for(k=0;k<CurrVertex;k++)
//           if(VertexDel[k]==cell->GetVertex((j+1)%N_Vertices))
//            {
//             ID = 1; break;
//            }
//             if(ID!=1)
//            {
//             VertexDel[CurrVertex] = cell->GetVertex((j+1)%N_Vertices);
//             CurrVertex++;
//            }
//         } // for j
//     } // for i
//   for(i=0;i<CurrVertex;i++)
//     delete VertexDel[i];
// 
//   delete []VertexDel;
//   OutPut(CurrVertex<<" vertices were deleted"<<endl);
// 
//  // remove all existing cells and joints
// //   for(i=0;i<N_RootCells;i++)
// //     delete (TGridCell*)CellTree[i];
// //   OutPut(N_RootCells<<" cells were deleted"<<endl);
// //    delete CellTree;
// //    delete coll;
// // 
// //   N_RootCells = Out.numberoftriangles;
// // 
// //   // allocate auxillary fields
// //   Coordinates = Out.pointlist;
// //   Triangles = Out.trianglelist;
// //   PartMarker = new int[Out.numberofpoints];
// // 
// //   // generate new vertices
// //   N_G = Out.numberofpoints;
// //   NewVertices = new TVertex*[N_G];
// // 
// //   for (i=0;i<N_G;i++)
// //      NewVertices[i] = new TVertex(Coordinates[2*i], Coordinates[2*i+1]);
// // 
//       // set bounding box
//   left = bottom = 1e8;
//   right = top = -1e8;
// // 
// //    for(i=0;i<In.numberofpoints;i++)
// //     {
// //       if(left>In.pointlist[2*i]) left = In.pointlist[2*i];
// //       if(right<In.pointlist[2*i]) right = In.pointlist[2*i];
// //       if(top<In.pointlist[2*i+1]) top = In.pointlist[2*i+1];
// //       if(bottom>In.pointlist[2*i+1]) bottom = In.pointlist[2*i+1];
// //     }
// // 
// //   OutPut("left: "<<left<<" right: "<<right<<" top: "<<top<<" bottom: "<<bottom<<endl);
// // 
// //   Domain->SetBoundBox(right-left,top-bottom);
// //   Domain->SetBoundBoxstart(left,bottom);
// // 
// //  // Solid Bound startx, starty, x length and y length
// //   UpdateSlipBound[0]->SetParams(Xi[0], Yi[0], Xi[1]-Xi[0],Yi[1]-Yi[0]);
// //   UpdateSlipBound[1]->SetParams(Xi[1], Yi[1], Xi[2]-Xi[1],Yi[2]-Yi[1]);
// //   UpdateSlipBound[2]->SetParams(Xi[2], Yi[2], Xi[3]-Xi[2],Yi[3]-Yi[2]);
// //   UpdateSlipBound[3]->SetParams(Xi[3], Yi[3], Xi[0]-Xi[3],Yi[0]-Yi[3]);
// // 
// // // Free boundary xmid, ymid, radius_a, radius_b, start angle, end angle
// //   UpdateIntface->SetParams(C_x, C_y, T_a, T_b, phi1, phi2);
// // 
// // //  cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[2*N_SlipBound_Vert]<< endl;
// //  // cout << " In.pointlist[N_SlipBound_Vert] : " << In.pointlist[0]<< endl;
// //  // generate cells
// //    CellTree = new TBaseCell*[N_RootCells];
// // 
// //   for (i=0;i<N_RootCells;i++)
// //   {
// //     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Triangle], 0);
// // 
// //     CellTree[i]->SetVertex(0, NewVertices[Out.trianglelist[3*i    ]]);
// //     CellTree[i]->SetVertex(1, NewVertices[Out.trianglelist[3*i + 1]]);
// //     CellTree[i]->SetVertex(2, NewVertices[Out.trianglelist[3*i + 2]]);
// // 
// //       ((TMacroCell *) CellTree[i])->SetSubGridID(0);
// //   }
// // 
// //   Domain->SetTreeInfo(CellTree, N_RootCells);
// //  
// //   // initialize iterators
// //   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
// //   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
// //   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
// //   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
// //   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);
// // 
// //   // search neighbours
// //   N_G = Out.numberofpoints;
// //   PointNeighb = new int[N_G];
// // 
// //   memset(PointNeighb, 0, N_G *SizeOfInt);
// // 
// //   for (i=0;i<3*N_RootCells;i++)
// //     PointNeighb[Triangles[i]]++;
// // 
// //   maxEpV=0;
// // 
// //   for (i=0;i<N_G;i++)
// //     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
// //   delete [] PointNeighb;
// // 
// //   PointNeighb = new int[++maxEpV * N_G];
// // 
// //   memset(PointNeighb, 0, maxEpV * N_G *SizeOfInt);
// // 
// //    // first colomn contains the number of following elements
// //    // for every point at first column we set the number of neighbour points
// //    // at further columns we set the index of corresponding cells
// //   for(i=0;i<3*N_RootCells;i++)
// //   {
// //     j = Triangles[i]*maxEpV;
// //     PointNeighb[j]++;
// //     PointNeighb[j + PointNeighb[j]] = i / 3;
// //   }
// // 
// //   // generate new edges
// //   N_G = Out.numberofedges;
// //   for (i=0;i<N_G;i++)
// //   {
// //     a = Out.edgelist[2*i];
// //     b = Out.edgelist[2*i+1];
// //     Neib[0] = -1;
// //     Neib[1] = -1;
// //     CurrNeib = 0;
// // 
// //     len1 = PointNeighb[a*maxEpV];
// //     len2 = PointNeighb[b*maxEpV];
// // 
// //   // find indexes of cells containing the current edge
// //    for (j=1;j<=len1;j++)
// //     {
// //       Neighb_tmp = PointNeighb[a*maxEpV + j];
// //        for (k=1;k<=len2;k++)
// //         if (Neighb_tmp == PointNeighb[b*maxEpV + k])
// //         {
// //           Neib[CurrNeib++] = Neighb_tmp;
// //           break;
// //         }
// //       if (CurrNeib == 2) break;
// //     }
// // 
// //     if (Out.edgemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
// //     {
// //       CurrComp = Out.edgemarkerlist[i] - 1;
// //       if (CurrComp >= 100000) CurrComp -= 100000;
// // 
// // 
// //       if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
// //             NewVertices[a]->GetX(), NewVertices[a]->GetY(), T_a) ||
// //           Domain->GetBdPart(0)->GetBdComp(CurrComp)->GetTofXY(
// //             NewVertices[b]->GetX(), NewVertices[b]->GetY(), T_b))
// //        {
// //           cerr<<"Error: could not set parameter values"<<endl;
// //           OutPut(NewVertices[a]<<endl);
// //           OutPut(NewVertices[b]<<endl);
// //         //  exit(0);
// //        }
// // 
// //       if (CurrNeib == 2)    // 2 cells contain the current edge
// //         if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
// //           Joint = new TIsoInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
// //                   T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
// //         else
// //           Joint = new TInterfaceJoint(Domain->GetBdPart(0)->GetBdComp(CurrComp),
// //                   T_a, T_b, CellTree[Neib[0]], CellTree[Neib[1]]);
// //       else
// //         if(Domain->GetBdPart(0)->GetBdComp(CurrComp)->IsFreeBoundary())
// //           Joint = new TIsoBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
// //         else
// //           Joint = new TBoundEdge(Domain->GetBdPart(0)->GetBdComp(CurrComp), T_a, T_b);
// //     }
// //     else // inner edge
// //     {
// //     if (CurrNeib != 2)
// //         cerr << "Error!!!!!!!! not enough neighbours!" << endl;
// // 
// //     Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
// //     }
// // 
// //     // find the local index for the point 'a' on the cell
// //     for (j=0;j<3;j++)
// //       if (Triangles[3*Neib[0]+j] == a) break;
// // 
// //     // find the local index for the point 'b' on the cell
// //     for (k=0;k<3;k++)
// //       if (Triangles[3*Neib[0]+k] == b) break;
// // 
// //    k = k*10 + j;
// // 
// //     switch (k)
// //     {
// //       case  1:
// //       case 10:
// //         j = 0;
// //         break;
// //       case 12:
// //       case 21:
// //         j = 1;
// //         break;
// //       case  2:
// //       case 20:
// //         j = 2;
// //         break;
// //     }
// // 
// //   CellTree[Neib[0]]->SetJoint(j, Joint);
// // 
// //    if (Neib[1] != -1)
// //     {
// //       // find the local index for the point 'a' on the cell
// //       for (j=0;j<3;j++)
// //         if (Triangles[3*Neib[1]+j] == a) break;
// // 
// //       // find the local index for the point 'b' on the cell
// //       for (k=0;k<3;k++)
// //         if (Triangles[3*Neib[1]+k] == b) break;
// // 
// //       k = k*10 + j;
// // 
// //       switch (k) // j will contain the local index for the current
// //       {
// //         case  1:
// //         case 10:
// //           j = 0;
// //           break;
// //         case 12:
// //         case 21:
// //           j = 1;
// //           break;
// //         case  2:
// //         case 20:
// //           j = 2;
// //           break;
// //       }
// // 
// //       CellTree[Neib[1]]->SetJoint(j, Joint);
// //     }
// // 
// //   if (Joint->GetType() == InterfaceJoint ||
// //         Joint->GetType() == IsoInterfaceJoint)
// //       ((TInterfaceJoint *) Joint)->CheckOrientation();
// //   }
// // 
// // 
// //   delete [] NewVertices;
// //   delete [] PointNeighb;
// //   delete [] In.pointlist;
// //   delete [] In.pointmarkerlist;
// //   delete [] In.segmentlist;
// //   delete [] In.segmentmarkerlist;
// //   
// //   if(Out.pointlist!=NULL) {
// //     free(Out.pointlist); Out.pointlist = NULL;}
// //   if(Out.pointattributelist!=NULL) { 
// //     free(Out.pointattributelist); Out.pointattributelist = NULL;}
// //   if(Out.pointmarkerlist!=NULL) {
// //     free(Out.pointmarkerlist); Out.pointmarkerlist = NULL;}
// //   if(Out.trianglelist!=NULL) {
// //     free(Out.trianglelist); Out.trianglelist = NULL;}
// //   if(Out.triangleattributelist!=NULL) {
// //     free(Out.triangleattributelist); Out.triangleattributelist = NULL;}
// //   if(Out.trianglearealist!=NULL) {
// //     free(Out.trianglearealist); Out.trianglearealist = NULL;}
// //   if(Out.neighborlist!=NULL) {
// //     free(Out.neighborlist); Out.neighborlist = NULL;}
// //   if(Out.segmentlist!=NULL) {
// //     free(Out.segmentlist); Out.segmentlist = NULL;}
// //   if(Out.segmentmarkerlist!=NULL) {
// //     free(Out.segmentmarkerlist); Out.segmentmarkerlist = NULL;}
// //   if(Out.holelist!=NULL) {
// //     free(Out.holelist); Out.holelist = NULL;}
// //   if(Out.regionlist!=NULL) {
// //     free(Out.regionlist); Out.regionlist = NULL;}
// //   if(Out.edgelist!=NULL) {
// //     free(Out.edgelist); Out.edgelist = NULL;}
// //   if(Out.edgemarkerlist!=NULL) {
// //     free(Out.edgemarkerlist); Out.edgemarkerlist = NULL;}
// //   if(Out.normlist!=NULL) {
// //     free(Out.normlist); Out.normlist = NULL;} 
//   
// // //======================================================================
// // // Triangular for grid generation --end
// // //  split the doamin in to interior and exterior -- begin
// // //=====================================================================
// // int sp_No, sp_no0, sp_no1;
// // double x_short, y_short,  temp0, temp1;
// // double x_mid, y_mid, hmin, tx, ty, nx, ny, Pos_Indicator;
// // 
// //     coll=Domain->GetCollection(It_Finest, 0);
// //     N_Cells = coll->GetN_Cells();
// //     N_Cells_P1 = 0;
// //     N_Cells_P2 = 0;
// // 
// //   for(i=0;i<N_Cells;i++)
// //    {
// //     cell = coll->GetCell(i);
// // 
// //      x_mid = 0.;  y_mid = 0.;
// // 
// //     for(j=0;j<3;j++)
// //      {
// //       cell->GetVertex(j)->GetCoords(x, y);
// //       x_mid +=x;
// //       y_mid +=y;
// //      } // for j
// // 
// //       x_mid /=3.;
// //       y_mid /=3.;
// //       hmin = 1000.0;
// // // find a point P on  the interface for the i'th cell mid point X
// // // such that ||P - X|| is minimum
// //       for(k=0;k<N_Old_Face_Vert;k++)
// //        {
// //         temp = sqrt( (I_FaceX[k]-x_mid)*(I_FaceX[k]-x_mid) + (I_FaceY[k]-y_mid)*(I_FaceY[k]-y_mid) );
// //         if(temp<hmin)
// //          {
// //           hmin=temp;  T_a = I_FaceX[k];  T_b = I_FaceY[k];
// //           x_short = T_a;  y_short = T_b;
// //           sp_No = k;
// //          }
// //        } // for k
// // 
// // // find next shortest point 
// // // i.e other shortest point of a spline containing (T_a,T_b)
// // // previous point
// //    if(sp_No==0)
// //     sp_no0 = N_Old_Face_Vert-1;
// //    else
// //     sp_no0 = sp_No-1;
// // 
// //     temp0 = sqrt( (I_FaceX[sp_no0]-x_mid)*(I_FaceX[sp_no0]-x_mid)
// //              + (I_FaceY[sp_no0]-y_mid)*(I_FaceY[sp_no0]-y_mid) );
// // 
// // // next point 
// //    if(sp_No==N_Old_Face_Vert-1)
// //     sp_no1 = 0;
// //    else
// //     sp_no1 = sp_No+1;
// // 
// //     temp1 = sqrt( (I_FaceX[sp_no1]-x_mid)*(I_FaceX[sp_no1]-x_mid) 
// //                 + (I_FaceY[sp_no1]-y_mid)*(I_FaceY[sp_no1]-y_mid) );
// // 
// //    if( temp0 < temp1)
// //      {
// //      T_a = I_FaceX[sp_no0];  T_b = I_FaceY[sp_no0];
// //      C_x = I_FaceX[sp_No];  C_y = I_FaceY[sp_No];
// //      }
// //    else 
// //      {
// //      T_a = I_FaceX[sp_No];  T_b = I_FaceY[sp_No];
// //      C_x = I_FaceX[sp_no1];  C_y = I_FaceY[sp_no1];
// //      }
// // 
// //      tx = x_mid - x_short; //x distance between the point and the shortest distance point on interface
// //      ty = y_mid - y_short; //y distance between the point and the shortest distance point on interface
// // 
// //      nx = -(C_y - T_b); // (-) normal at (T_a,T_b) pointing into the inner domain
// //      ny = -(T_a - C_x); // (-) normal at (T_a,T_b) pointing into the inner domain
// // 
// //      Pos_Indicator = tx*nx + ty*ny;
// // 
// //    if(Pos_Indicator > 0.)
// //     {
// //      //  cell is in inner domain
// //      cell->SetPhase_ID(0);
// //      N_Cells_P1++;
// // 
// //     }
// //    else if(Pos_Indicator < 0.)
// //     {
// //      //  cell is in outer domain
// //      cell->SetPhase_ID(1);
// //      N_Cells_P2++;
// //     }
// //    else
// //     {
// //      cout<< " error in identifying phase check remesh2d.c" << endl;
// //      exit(0);
// //     }
// //   }// for i
// // 
// // 
// //   delete [] I_FaceX;
// //   delete [] I_FaceY;
// // 
// //  
// 
// 
// }
// 
// 
// 




