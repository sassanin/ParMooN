// ==========================================================================
// instationary problem
// ==========================================================================
//===========================================================================
// example file
// =========================================================================
// exact solution in unit cube
void ExampleFile()
{
#define __ROBINBC__ 
  OutPut("Example: import3Dmesh.h" << endl);
}

void GridBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

void SolidBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}

void FluidBoundCondition(double x, double y, double z, BoundCond &cond)
{
  cond = DIRICHLET;
}
// exact solution
void Exact(double x, double y, double z, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
  values[4] = 0;
}

// initial conditon
void InitialCondition(double x, double y, double z, double *values)
{
 values[0] = TDatabase::ParamDB->P0;
 values[0] = x;
}

// kind of boundary condition
void BoundCondition(double x, double y, double z, BoundCond &cond)
{
 

   if(TDatabase::ParamDB->P13==0)
    { //Lside 
     if(x<225 && y>76 )  
       { cond = DIRICHLET; }
     else
       { cond =  NEUMANN;  }
    }
   else
   {
    if(x< 240 && y<-72. && z>191)  
     { cond = DIRICHLET; }
    else
     { cond =  NEUMANN;  }
   }   
}

// value of boundary condition
void BoundValue(double x, double y, double z, double &value)
{
  value = 0;
}

void BilinearCoeffs(int n_points, double *X, double *Y, double *Z, double **parameters, double **coeffs)
{
  double eps=0;
  int i, Region_ID;
  double *coeff;
  double x, y, z, hK;
 
  
  double rhsfact, alpha,  char_L, beam_r, DimlessBeam_r;
  double xp=0., yp=0., zp=0., SourceCoord, Sourceradius;
    
   if(TDatabase::ParamDB->P13==0)
    { //   center point is (291.977  -1.8733333333333333 159.89 ), old, wrong
      //   2.922130000000000e+02 -1.8733333333333333e+00 1.583800000000000e+02 1       
      xp=2.922130000000000e+02;
      zp=1.583800000000000e+02; 
    }   
   else
    { //  side  300.565 69.66666666666666 157.038 
      xp= 300.565;
      zp=157.038; 
    }  
 
    
  beam_r = TDatabase::ParamDB->P11;
  char_L = TDatabase::ParamDB->P12;

  hK = coeffs[0][2];  
  Region_ID = coeffs[0][1];
  DimlessBeam_r =  beam_r/char_L;
  

  // cout << "Region_ID = " << Region_ID << endl;
   
   if(Region_ID==0) //0-Fat/adipost
    {
     eps = TDatabase::ParamDB->P1/(1100.*4483*char_L*char_L); ;   
     alpha = TDatabase::ParamDB->P7;     
     rhsfact = (alpha*TDatabase::ParamDB->P4)/(1100.0*4483*(22./7.)*beam_r*beam_r);   
    }
   else  
    {
     eps = TDatabase::ParamDB->P2/(1040.0*3500.*char_L*char_L); 
     alpha = TDatabase::ParamDB->P8;        
     rhsfact = alpha*TDatabase::ParamDB->P4/(1040.0*3500.*(22./7.)*beam_r*beam_r);
    }

//      //2-CSF  
//      eps = TDatabase::ParamDB->P2/(997.0*3710*char_L*char_L); 
//      alpha = TDatabase::ParamDB->P8;        
//      rhsfact = alpha*TDatabase::ParamDB->P4/(997.0*3710*(22./7.)*beam_r*beam_r);
//      
  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];

    x = X[i];
    y = Y[i];
    z = Z[i];

    // diffusion
    coeff[0] = eps;
    // convection in x direction
    coeff[1] = 0;
    // convection in y direction
    coeff[2] = 0;
    // convection in z direction
    coeff[3] = 0;
    // reaction
    coeff[4] = 0;
     // rhs

     if(TDatabase::ParamDB->P13==0)
      { //Lside  
       Sourceradius = sqrt( (x-xp)*(x-xp) + (z-zp)*(z-zp));
       
       if(y<0.) y=0.;
       SourceCoord = -y;
      }   
     else
      { //  side  
       Sourceradius = sqrt( (x-xp)*(x-xp) + (z-zp)*(z-zp));
       
       if(y>0.) y=0.;       
       SourceCoord = y;
      } 
      
    
    if(Sourceradius<=DimlessBeam_r)
     {           
      coeff[5] = rhsfact*exp(alpha*SourceCoord*char_L);// f  
      
//       cout << "Sourceradius : " << Sourceradius << endl;
      if(TDatabase::ParamDB->P14< hK) TDatabase::ParamDB->P14=hK;        
     }
    else
     {coeff[5] = 0; }// f
     
    coeff[6] = 0;
    
  }
}

void ReadAdditionalNModes(char *NODE, tetgenio &InAddPts)
{
  int i, N_Vertices, index, N_Atribude, BDMarker;
  char line[100];  
 std::ifstream dat(NODE);

  if (!dat)
  {
    cerr << "cannot open '" << NODE << "' for input" << endl;
    exit(-1);
  }      
  
  dat.getline (line, 99); // first line is command line    
  dat >> N_Vertices >> N_Atribude >> BDMarker;  
  dat.getline (line, 99);  
 
  InAddPts.numberofpoints = N_Vertices;
  InAddPts.pointlist = new double[3*N_Vertices];  

  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> index >> InAddPts.pointlist[3*i] >> InAddPts.pointlist[3*i+1] >> InAddPts.pointlist[3*i+2]; 
     cout<< i << " vert X: " <<InAddPts.pointlist[3*i] << " vert Y: " <<InAddPts.pointlist[3*i+1] <<" vert Z: " <<InAddPts.pointlist[3*i+2] <<endl;
 
    }
 
   dat.close();

// cout << "N_Vertices " <<N_Vertices <<endl;
// 
// exit(0);

} // ReadAdditionalNModes

void ReadMeditMeshWithCells(char *SMESH, tetgenio &In)
{  
  int i, j, k, dimension, N_FVert, N_Faces, N_Vertices, N_CVert, N_Cells;
  int BDComp_Min=10000000;
  int *plist;
  
  double x;
  
  char line[100];
  tetgenio::facet *F;
  tetgenio::polygon *P;    

 std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
    exit(-1);
  }     
   
  // check the dimension
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Dimension"))  ||  (!strcmp(line, "dimension")) ||  (!strcmp(line, "DIMENSION")))
    {
     dat.getline (line, 99);
     dat >> dimension;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }
  
  if(dimension!=3)
   {
    cerr << "dimension: " << dimension << endl;
    exit(-1);
   }
  
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Vertices")) ||  (!strcmp(line, "vertices"))   ||  (!strcmp(line, "VERTICES"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_Vertices;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  } 
  
  In.numberofpoints = N_Vertices;
  In.pointlist = new double[3*N_Vertices];  

  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> x >> In.pointlist[3*i+1] >> In.pointlist[3*i+2];
     
     In.pointlist[3*i] = x - 3.606596999999999830777142; //center point
     
//      if(i==1063 || i==1061 || i==1062 || i==162 || i==1175 || i==517 || i==1064 || i==518 ) // i==1063  as center
     if(i==1062 || i==916 || i==341 || i==914 || i==162 || i==1063 || i==1061) // i==1062  as center
      cout<< i << " vert X: " <<In.pointlist[3*i] << " vert Y: " <<In.pointlist[3*i+1] <<" vert Z: " <<In.pointlist[3*i+2] <<endl;   
    }
    
  // find the N_Triangles
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Triangles")) ||  (!strcmp(line, "triangles"))   ||  (!strcmp(line, "TRIANGLES"))   ) 
    {
     N_FVert = 3;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
    else if ( (!strcmp(line, "Quadrilaterals")) ||  (!strcmp(line, "quadrilaterals"))   ||  (!strcmp(line, "QUADRILATERALS"))   ) 
    {
     N_FVert = 4;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
   
    // read until end of line
    dat.getline (line, 99);   
  } 
 
    In.numberoffacets = N_Faces;
    In.facetlist = new tetgenio::facet[In.numberoffacets];
    In.facetmarkerlist = new int[In.numberoffacets];
    
    for(i=0;i<N_Faces; i++)
     {
      dat.getline (line, 99);       
      
      F = &In.facetlist[i];
      tetgenio::init(F);      
      F->numberofpolygons = 1;
      F->polygonlist = new tetgenio::polygon[F->numberofpolygons];
      F->numberofholes = 0;
      F->holelist = NULL;
      P = &F->polygonlist[0];
      tetgenio::init(P);
      P->numberofvertices = N_FVert;
      P->vertexlist = new int[P->numberofvertices];
      
      for(j=0;j<N_FVert;j++)
      {
        dat >> k;
        P->vertexlist[j] = k-1;  // c numbering 
      }
      
      dat >>  In.facetmarkerlist[i];        
      
      if(BDComp_Min > In.facetmarkerlist[i])
        BDComp_Min=In.facetmarkerlist[i];
      
      // cout << i<<  " P->vertexlist[j]:  " <<  In.facetmarkerlist[i]  << endl;  
     } //   for(i=0;i<N_Faces; i++)
   
    BDComp_Min--;
   
    for(i=0;i<N_Faces; i++)
      In.facetmarkerlist[i] -= BDComp_Min;
         
    
    /**load cells */
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Tetrahedron")) ||  (!strcmp(line, "tetrahedron"))   ||  (!strcmp(line, "TETRAHEDRON"))   ) 
    {
     N_CVert = 4;
     dat.getline (line, 99);
     dat >> N_Cells;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }     
 
     In.numberoftetrahedra = N_Cells;
     In.numberofcorners = N_CVert;
     In.numberoftetrahedronattributes = 1;
     In.tetrahedronlist = new int[N_Cells * 4];
     In.tetrahedronattributelist = new REAL[N_Cells];   
     In.tetrahedronvolumelist = new REAL[N_Cells];  
     
     for(i=0; i<N_Cells; i++) 
      {  
       dat.getline (line, 99);       
       plist = &(In.tetrahedronlist[i * 4]);

       for(j=0;j<N_CVert;j++)
        {
         dat >> k;
         plist[j] = k -1;       
        }//  for(j=0;j<N_CV
        
        
       dat >> In.tetrahedronattributelist[i];   
//        In.tetrahedronvolumelist [i] = 100.;
      } // for(i=0; i<N_
     
     
     
//      cout << i<<  "  N_Cells :  " <<  N_Cells << endl;     
//      exit(0);
     
       dat.close();
//    exit(0);    
} // ReadMeditMeshWithCells 

void ReadMeditMesh(char *SMESH, tetgenio &In)
{  
  int i, j, k, dimension, N_FVert, N_Faces, N_Vertices;
  int BDComp_Min=10000000;
  
  char line[100];
  tetgenio::facet *F;
  tetgenio::polygon *P;    

 std::ifstream dat(SMESH);

  if (!dat)
  {
    cerr << "cannot open '" << SMESH << "' for input" << endl;
    exit(-1);
  }     
   
  // check the dimension
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Dimension"))  ||  (!strcmp(line, "dimension")) ||  (!strcmp(line, "DIMENSION")))
    {
     dat.getline (line, 99);
     dat >> dimension;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }
  
  if(dimension!=3)
   {
    cerr << "dimension: " << dimension << endl;
    exit(-1);
   }
  
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Vertices")) ||  (!strcmp(line, "vertices"))   ||  (!strcmp(line, "VERTICES"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_Vertices;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  } 
  
  In.numberofpoints = N_Vertices;
  In.pointlist = new double[3*N_Vertices];  

  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> In.pointlist[3*i] >> In.pointlist[3*i+1] >> In.pointlist[3*i+2];
     //cout<< i << " vert X: " <<In.pointlist[3*i] << " vert Y: " <<In.pointlist[3*i+1] <<endl;   
    }
    
  // find the N_Triangles
  // find the N_Vertices
  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Triangles")) ||  (!strcmp(line, "triangles"))   ||  (!strcmp(line, "TRIANGLES"))   ) 
    {
     N_FVert = 3;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
    else if ( (!strcmp(line, "Quadrilaterals")) ||  (!strcmp(line, "quadrilaterals"))   ||  (!strcmp(line, "QUADRILATERALS"))   ) 
    {
     N_FVert = 4;
     dat.getline (line, 99);
     dat >> N_Faces;
     break;
    }    
   
    // read until end of line
    dat.getline (line, 99);   
  } 
   
    In.numberoffacets = N_Faces;
    In.facetlist = new tetgenio::facet[In.numberoffacets];
    In.facetmarkerlist = new int[In.numberoffacets];
    
    for(i=0;i<N_Faces; i++)
     {
      dat.getline (line, 99);       
      
      F = &In.facetlist[i];
      tetgenio::init(F);      
      F->numberofpolygons = 1;
      F->polygonlist = new tetgenio::polygon[F->numberofpolygons];
      F->numberofholes = 0;
      F->holelist = NULL;
      P = &F->polygonlist[0];
      tetgenio::init(P);
      P->numberofvertices = N_FVert;
      P->vertexlist = new int[P->numberofvertices];
      
      for(j=0;j<N_FVert;j++)
      {
        dat >> k;
        P->vertexlist[j] = k-1;  // c numbering 
      }
      
      dat >>  In.facetmarkerlist[i];        
      
      if(BDComp_Min > In.facetmarkerlist[i])
        BDComp_Min=In.facetmarkerlist[i];
      
      // cout << i<<  " P->vertexlist[j]:  " <<  In.facetmarkerlist[i]  << endl;  
     } //   for(i=0;i<N_Faces; i++)
   
    BDComp_Min--;
   
    for(i=0;i<N_Faces; i++)
      In.facetmarkerlist[i] -= BDComp_Min;
         
//         cout << i<<  " P->vertexlist[j]:  " <<  BDComp_Min << endl;     
     
     
       dat.close();
//    exit(0);    
} // ReadMeditMesh

void DeleteDomain(TDomain *&Domain)
{
  int i, j, k, N_Cells, N_RootCells;
  int CurrVertex, N_Faces, N_Vertices, ID;  
  TBaseCell **CellTree,  **SurfCellTree, *cell;
  TGridCell **DelCell;
  TVertex **VertexDel;  
  TCollection *coll;  
  
    Domain->GetTreeInfo(CellTree,N_RootCells);
    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();    

    VertexDel = new TVertex*[8*N_RootCells];
    CurrVertex = 0;   
    
   for(i=0;i<N_Cells;i++)
     {
       cell = coll->GetCell(i);
       N_Faces = cell->GetN_Faces();
       N_Vertices = cell->GetN_Vertices();
       for(j=0;j<N_Faces;j++)
         {
          if(CurrVertex==0)
           {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
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
               VertexDel[CurrVertex++] = cell->GetVertex(j);
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
             VertexDel[CurrVertex++] = cell->GetVertex((j+1)%N_Vertices);
            }
 
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+2)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+2)%N_Vertices);
            }
         if(N_Faces==6) // If old cell is hexahedrol
          {  
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+4)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+4)%N_Vertices);
            }        
        }
      } // for j
    } // for i
    
    for(i=0;i<CurrVertex;i++)
     delete VertexDel[i];
     delete [] VertexDel;
     OutPut(CurrVertex<<" vertices were deleted"<<endl);

  // remove all existing cells and joints
   for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
    delete [] CellTree;
    OutPut(N_RootCells<<" cells were deleted"<<endl);  
}

TDomain* TetrameshCreate()
{
  int i, j, k, l, dimension, N_Vertices;
  int N_Faces, N_RootCells;
  int v1, v2, v3, v4, CellMarker, RefLevel=0, *Tetrahedrals, *PointNeighb, maxEpV, *Tetrahedrals_loc;
  int a, b, c, Neib[2], Neighb_tmp, CurrNeib, len1, len2, len3, CurrComp=0, N_Points ;
  int N_Cells, MaxLen, jj, CompID=0;  
  const int *EdgeVertex, *TmpFV, *TmpLen;  
  double X, Y, Z, DispX, DispY, DispZ;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ;   
  char *MESH, line[100];
  MESH = TDatabase::ParamDB->SMESHFILE;
  TDomain *Domain= new TDomain;
  TVertex **NewVertices;
  TBaseCell **CellTree, *cell;  
  TJoint *Joint;  
  TBoundComp3D *bdcomp;  
  TShapeDesc *ShapeDesc;
  TCollection *coll;
  TBoundPart *someboundpart;
  std::ifstream dat(MESH);

  if (!dat)
  {
    cerr << "cannot open '" << MESH << "' for input" << endl;
    exit(-1);
  }

  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Dimension"))  ||  (!strcmp(line, "dimension")) ||  (!strcmp(line, "DIMENSION")))
    {
     dat.getline (line, 99);
     dat >> dimension;
     break;
    }    

    dat.getline (line, 99);   
  }

  if(dimension!=3)
   {
    cerr << "dimension: " << dimension << endl;
    exit(-1);
   }

  while (!dat.eof())
  {
    dat >> line;    
    if ( (!strcmp(line, "Vertices")) ||  (!strcmp(line, "vertices"))   ||  (!strcmp(line, "VERTICES"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_Vertices;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }   
   // generate new vertices
   cout <<"N_Vertices "<<N_Vertices<<endl;
   NewVertices = new TVertex*[N_Vertices];

   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> X >> Y >> Z;      
     NewVertices[i] = new TVertex(X, Y, Z);

      if (X > Xmax) Xmax = X;
      if (X < Xmin) Xmin = X;
      if (Y > Ymax) Ymax = Y;
      if (Y < Ymin) Ymin = Y;
      if (Z > Zmax) Zmax = Z;
      if (Z < Zmin) Zmin = Z;  
    } 
   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;

   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);   

  while (!dat.eof())
  {
    dat >> line;
    
    if ( (!strcmp(line, "Tetrahedra")) ||  (!strcmp(line, "tetrahedra"))   ||  (!strcmp(line, "TETRAHEDRA"))   ) 
    {
     dat.getline (line, 99);
     dat >> N_RootCells;
     break;
    }    
    // read until end of line
    dat.getline (line, 99);   
  }   

  // generate new cells
   cout << "number of root cells "<<N_RootCells <<endl;   
   CellTree = new TBaseCell*[N_RootCells];
   Tetrahedrals = new int[4*N_RootCells];
  
   for (i=0;i<N_RootCells;i++)
   {
     dat.getline (line, 99);
     dat >> v1 >> v2 >> v3 >> v4 >> CellMarker;  
     Tetrahedrals[4*i    ] = v1-1;
     Tetrahedrals[4*i + 1] = v2-1;
     Tetrahedrals[4*i + 2] = v3-1;
     Tetrahedrals[4*i + 3] = v4-1;   
     
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
     CellTree[i]->SetRegionID(CellMarker);
     CellTree[i]->SetAsLayerCell(1);      
     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);
     CellTree[i]->SetClipBoard(i);
     
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);
   }
   dat.close();
   
   Domain->SetTreeInfo(CellTree, N_RootCells);
   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   // search neighbours
   PointNeighb = new int[N_Vertices];     
   memset(PointNeighb, 0, N_Vertices*SizeOfInt);     
  
     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;
     
   maxEpV = 0;
   for (i=0;i<N_Vertices;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;        
   cout<<"maximum edges per vertex "<< maxEpV<<endl;
   PointNeighb = new int[++maxEpV * N_Vertices];
   memset(PointNeighb, 0, maxEpV*N_Vertices*SizeOfInt);

   // every vertex contains "maxEpV" columns
   // for every vertex at first colomn contains the number of cells containing this vertex
   // at further columns we set the index of corresponding cells containing this vertex
   // cout<<"maxEpV*N_Vertices "<<maxEpV*N_Vertices<<endl;
   for(i=0;i<4*N_RootCells;i++)
    {      
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
/** now generate faces */ //===================================================== 
   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//     for(i=0;i<N_Cells;i++)
//     {
//       cell = coll->GetCell(i);
//       ShapeDesc= cell->GetShapeDesc();   
//       ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
//       N_Faces = cell->GetN_Faces(); 
//          cout<<"N_Faces "<<N_Faces<<endl;
//       for(j=0;j<N_Faces;j++)
// 	cout<<"TmpLen["<<j<<"] is "<<TmpLen[j]<<endl;      
//     }
// exit(0);   
   for(i=0;i<N_Cells;i++)
    {
      cell = coll->GetCell(i);
      ShapeDesc= cell->GetShapeDesc();   
      ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      /** which vertices are on each face(local indexing) */      
      /** number of  vertices on each face */
      /** maximum number of vertices per face */
      N_Faces = cell->GetN_Faces();      
      Tetrahedrals_loc = Tetrahedrals+4*i;
//      cout <<"chek 9 "<< i<<endl;
       for(jj=0;jj<N_Faces;jj++)
         if(cell->GetJoint(jj) == NULL)
          {
           N_Points = TmpLen[jj];   
           if(N_Points!=3)
            {     
             cerr << "Only Tria faces are allowed!!! N_FVert: " << N_Points << endl;
             exit(-1);     
            }
           a = Tetrahedrals_loc[TmpFV[jj*MaxLen]];
           b = Tetrahedrals_loc[TmpFV[jj*MaxLen + 1]];
           c = Tetrahedrals_loc[TmpFV[jj*MaxLen + 2]];	   
// 	   cout<<"a "<<a<<endl;
// 	   cout<<"b "<<b<<endl;
// 	   cout<<"c "<<c<<endl;
// 	   cout<<"TmpFV1 "<<TmpFV[jj*MaxLen]<<endl;
// 	   cout<<"TmpFV2 "<<TmpFV[jj*MaxLen+1]<<endl;
// 	   cout<<"TmpFV3 "<<TmpFV[jj*MaxLen+2]<<endl;
// 	   	   cout<<"TmpFV1 "<<TmpFV[jj*MaxLen]<<endl;
// 	   cout<<"TmpFV2 "<<TmpFV[(jj+1)*MaxLen+1]<<endl;
// 	   cout<<"TmpFV3 "<<TmpFV[(jj+1)*MaxLen+2]<<endl;
// 	   	   cout<<"TmpFV1 "<<TmpFV[(jj+1)*MaxLen]<<endl;
// 	   cout<<"TmpFV2 "<<TmpFV[(jj+2)*MaxLen+1]<<endl;
// 	   cout<<"TmpFV3 "<<TmpFV[(jj+2)*MaxLen+2]<<endl;
// 	   	   cout<<"TmpFV1 "<<TmpFV[(jj+2)*MaxLen]<<endl;
// 	   cout<<"TmpFV2 "<<TmpFV[(jj+3)*MaxLen+1]<<endl;
// 	   cout<<"TmpFV3 "<<TmpFV[(jj+3)*MaxLen+2]<<endl;
// 	   exit(0);	   
// 	      a = TmpFV[jj*MaxLen];
//            b = TmpFV[jj*MaxLen + 1];
//            c = TmpFV[jj*MaxLen + 2];	   
           Neib[0] = -1;
           Neib[1] = -1;
           CurrNeib = 0;
           len1 = PointNeighb[a*maxEpV];
           len2 = PointNeighb[b*maxEpV];
           len3 = PointNeighb[c*maxEpV]; 
// 	   cout <<"chek 10 "<< i<<endl;
           // find the index of the cells containing current face with node(vertex) indices a,b,c
           for (j=1;j<=len1;j++)
            {
             Neighb_tmp = PointNeighb[a*maxEpV + j];
             for (k=1;k<=len2;k++)
              {
               if(Neighb_tmp == PointNeighb[b*maxEpV + k])
                {
                 for(l=1;l<=len3;l++)
                  if(Neighb_tmp == PointNeighb[c*maxEpV + l])
                   {
                    Neib[CurrNeib++] = Neighb_tmp;
                    break;
                   }
                 }
               }
             if (CurrNeib == 2) break;
            } // for (j=1;j<=len1;j++)    
//             cout <<"chek 11 "<< i<<endl;
        if(CurrNeib==1)
         { 
// 	   cout <<"chek 11.10 "<<endl;
	   someboundpart = Domain->GetBdPart(0);
//           if(someboundpart==NULL)	   
// 	   cout <<"chek 11.12 NULL"<<endl;
// 	   fflush(stdout);  
// 	   cout <<"CurrComp "<<CurrComp<<endl;
// 	  scanf("%d",&CurrComp);
	   bdcomp = new TBdPlane(CompID++);
// 	  bdcomp = someboundpart->GetBdComp(CurrComp);
// 	   cout <<"chek 11.15 "<<endl;
// 	   cout <<"chek 11.20 "<<endl;
// 	   fflush(stdout);
// 	  scanf("%d",&CurrComp);
// 	   exit(0);
// 	  double vertx, verty, vertz, t1, s1;
// 	  vertx = NewVertices[a]->GetX();
// 	  verty = NewVertices[a]->GetY();
// 	  vertz = NewVertices[a]->GetZ();
// 	   cout <<"chek 11.11 "<<NewVertices[a]->GetX()<<" "<<NewVertices[a]->GetY()<<" "<<NewVertices[a]->GetZ()<<endl;
// 	  cout <<"chek 11.11 "<<endl;
// 	   bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(), NewVertices[a]->GetZ(), T[1], S[1]);
// 	  bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(), NewVertices[a]->GetZ(), T[1], S[1]);
// 	  cout <<T[1]<<" "<<S[1]<<endl;
// 	  cout<<"b is "<<b<<endl;
// 	  cout <<NewVertices[239]->GetX()<<" "<<NewVertices[239]->GetY()<<endl;
// 	  bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(), NewVertices[b]->GetZ(), T[2], S[2]);
// 	  bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(), NewVertices[c]->GetZ(), T[3], S[3]);
//            bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);  
// 	   exit(0);
// 	   cout<<"chek 11.12"<<endl;
// 	   cout <<"chek 11.25 "<< i<<endl;
//            cout <<"chek 11.35 "<<NewVertices[a]->GetX()<<" "<<NewVertices[a]->GetY()<<" "<<NewVertices[a]->GetZ()<<endl;
// 	   cout <<"bdcomp "<<bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(), NewVertices[a]->GetZ(), T[1], S[1])<<endl;
// 	   exit(0);
           if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(), NewVertices[a]->GetZ(), T[1], S[1]) 
	                                            ||
              bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(), NewVertices[b]->GetZ(), T[2], S[2])
	                                            ||
              bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(), NewVertices[c]->GetZ(), T[3], S[3]) )
             {
	      cout <<"chek 11.5 "<< i<<endl;
              cerr<<"Error: could not set parameter values"<<endl;
              OutPut(NewVertices[a]<<endl);
              OutPut(NewVertices[b]<<endl);
              OutPut(NewVertices[c]<<endl);
              exit(0);
             }
           Joint = new TBoundFace(bdcomp, T, S);
// 	   cout <<"chek 11.75 "<< i<<endl;
         }           
        else
        {
         if (CurrNeib != 2)
         cerr << "Error !!!!!!!! not enough neighbours!" << endl;
         Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
        } 
//         cout <<"chek 12 "<< i<<endl;
        // First element containing the current face
        // find the local index for the point 'a' on the cell        
        for (j=0;j<4;j++)
         if (Tetrahedrals[4*Neib[0]+j] == a) break;
        // find the local index for the point 'b' on the cell
        for (k=0;k<4;k++)
         if (Tetrahedrals[4*Neib[0]+k] == b) break;
        // find the local index for the point 'c' on the cell
        for (l=0;l<4;l++)
         if (Tetrahedrals[4*Neib[0]+l] == c) break;
        l = l*100 + k*10 + j;       

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

       default:
        Error("Unable to set the face !!!!!!!!!!!!" << endl);
        exit(0);
      } // switch (l)  
      cout <<"chek 13 "<< i<<endl;
     CellTree[Neib[0]]->SetJoint(j, Joint);

     // second cell containing the current face 
    if (Neib[1] != -1) // second element containing the current face
     {  
       cout <<"chek 13.5 "<< i<<endl;
     // find the local index for the point 'a' on the cell
     for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

     // find the local index for the point 'b' on the cell
     for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

      // find the local index for the point 'c' on the cell
     for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;
     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
       }
       CellTree[Neib[1]]->SetJoint(j, Joint);   
     } // if (Neib[1] != -1)      
   
      if (Joint->GetType() == InterfaceJoint3D || Joint->GetType() == IsoInterfaceJoint3D)
        {
         ((TInterfaceJoint3D*)Joint)->SetMapType();
         ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
        }
       else 
        if (Joint->GetType() == JointEqN)
           Joint->SetMapType();       
       } //  if(cell->GetJoint(j)
      }  // for(i=0;i<N_Cells;i++)       
delete [] PointNeighb;
 cout<< "tetrameshcreatesuccessful " <<"\n";
 return Domain;
}

void TetrameshGen(TDomain *&Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, *N_FVerts, N_Faces, *Facets;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PointNeighb, dimension;
  int *Facelist, *Facemarkerlist;
 
  double *Coordinates, N_x, N_y, N_z;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out;
   
  TBaseCell **CellTree,  **SurfCellTree;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices; 
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;  

  char *SMESH, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;

#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;
    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
#endif  

  /** only the root call the mesh generator and broadcast the mesh info to all */
#ifdef _MPI    
  if(rank==0)
#endif  
   {
   
//     In.initialize();
//     Out.initialize();
   
    std::ostringstream opts;
    opts << " "; 
  
    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
//     opts<<'r'; // -r  Reconstructs a previously generated mesh     
//     opts<<"q"<<1.49; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<0.1; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
//     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
//     opts<<'V';  //verbose mode
    opts<<ends;    
  
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);  
  
//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);
    
   } // if(rank==0)
     
    
 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);

    Domain->GetTreeInfo(CellTree,N_RootCells);
    coll = Domain->GetCollection(It_Finest, 0);
    N_Cells = coll->GetN_Cells();

//     cout<<"N_RootCells: "<<N_RootCells<<endl;
   // remove all existing vertices and joints
    VertexDel = new TVertex*[8*N_RootCells];
    CurrVertex = 0;

   for(i=0;i<N_Cells;i++)
     {
       cell = coll->GetCell(i);
       N_Faces = cell->GetN_Faces();
       N_Vertices = cell->GetN_Vertices();
       for(j=0;j<N_Faces;j++)
         {
          if(CurrVertex==0)
           {
               VertexDel[CurrVertex++] = cell->GetVertex(j);
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
               VertexDel[CurrVertex++] = cell->GetVertex(j);
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
             VertexDel[CurrVertex++] = cell->GetVertex((j+1)%N_Vertices);
            }
 
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+2)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+2)%N_Vertices);
            }
         if(N_Faces==6) // If old cell is hexahedrol
          {  
           ID = 0;
           for(k=0;k<CurrVertex;k++)
           if(VertexDel[k]==cell->GetVertex((j+4)%N_Vertices))
            {
             ID = 1; break;
            }
             if(ID!=1)
            {
             VertexDel[CurrVertex++] = cell->GetVertex((j+4)%N_Vertices);
            }        
        }
      } // for j
    } // for i

    for(i=0;i<CurrVertex;i++)
     delete VertexDel[i];
     delete [] VertexDel;
     OutPut(CurrVertex<<" vertices were deleted"<<endl);

  // remove all existing cells and joints

   for(i=0;i<N_RootCells;i++)
    delete (TGridCell*)CellTree[i];
    delete [] CellTree;
    OutPut(N_RootCells<<" cells were deleted"<<endl);

   
#ifdef _MPI     
     if(rank==0) 
#endif  
     {
       N_RootCells = Out.numberoftetrahedra;
       N_G = Out.numberofpoints;   

       // allocate auxillary fields
       Coordinates = Out.pointlist;  
       Tetrahedrals = Out.tetrahedronlist;       
     }

     
     
#ifdef _MPI
   MPI_Bcast(&N_RootCells, 1, MPI_INT, 0, Comm);
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);  
   
   
  if(rank!=0)
   {
    Coordinates = new double [3*N_G];
    Tetrahedrals = new int[4*N_RootCells];
   }
   
   MPI_Bcast(Coordinates, 3*N_G, MPI_DOUBLE, 0, Comm);  
   MPI_Bcast(Tetrahedrals, 4*N_RootCells, MPI_INT, 0, Comm);     
#endif      
   
  // generate new vertices
   NewVertices = new TVertex*[N_G];

   for (i=0;i<N_G;i++)
    {
      NewVertices[i] = new TVertex(Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);

      // set bounding box
      if (Coordinates[3*i] > Xmax) Xmax = Coordinates[3*i];
      if (Coordinates[3*i] < Xmin) Xmin = Coordinates[3*i];
      if (Coordinates[3*i+1] > Ymax) Ymax = Coordinates[3*i+1];
      if (Coordinates[3*i+1] < Ymin) Ymin = Coordinates[3*i+1];
      if (Coordinates[3*i+2] > Zmax) Zmax = Coordinates[3*i+2];
      if (Coordinates[3*i+2] < Zmin) Zmin = Coordinates[3*i+2];
   }
   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;

   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);   
//        cout<<Xmin <<"  "<<Ymin <<"  "<<Zmin<<endl;
//        cout<<Xmax <<"  "<<Ymax <<"  "<<Zmax<<endl;

   CellTree = new TBaseCell*[N_RootCells];
 //  output of each tetraheron vertex indices (four vertices for each)
//   for (i=0;i<N_RootCells;i++)
//        cout<< Tetrahedrals[4*i]<<"  "<<Tetrahedrals[4*i + 1]<<"  "
//          <<Tetrahedrals[4*i + 2]<<"  "<<Tetrahedrals[4*i + 3]<<endl;


   for (i=0;i<N_RootCells;i++)
   {
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],RefLevel);

     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);

    CellTree[i]->SetClipBoard(i);
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
   PointNeighb = new int[N_G];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

#ifdef _MPI    
     if(rank==0)   
#endif   
   cout<<"maxEpV "<< maxEpV<<endl;

   PointNeighb = new int[++maxEpV * N_G];

   memset(PointNeighb, 0, maxEpV*N_G*SizeOfInt);

    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of corresponding cells containing this vertex
    // cout<<"maxEpV*N_G "<<maxEpV*N_G<<endl;

   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
 //  output of PointNeighb columns for each point
//   for (i=0;i<N_G;i++)
//    {
//     for (j=0;j<maxEpV;j++)
//     cout<<"  "<< PointNeighb[i*maxEpV+j];
//     cout<<endl;
//    }

   // generate new faces 
  
#ifdef _MPI     
     if(rank==0) 
#endif  
     {    
      N_G =  Out.numberoftrifaces;
      Facelist = Out.trifacelist;
      Facemarkerlist = Out.trifacemarkerlist;
     }
      
#ifdef _MPI
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);  

  if(rank!=0)
   {
    Facelist = new int [3*N_G];
    Facemarkerlist = new int[N_G];
   }   
   
   
   MPI_Bcast(Facelist, 3*N_G, MPI_INT, 0, Comm); 
   MPI_Bcast(Facemarkerlist, N_G, MPI_INT, 0, Comm); 
  
  if(rank==0) 
#endif      
   cout<<"numberoftrifaces "<<N_G<<endl;
    
   for (i=0;i<N_G;i++)
   {
     a = Facelist[3*i];
     b = Facelist[3*i+1];
     c = Facelist[3*i+2];

//      cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;

     Neib[0] = -1;
     Neib[1] = -1;
     CurrNeib = 0;

     len1 = PointNeighb[a*maxEpV];
     len2 = PointNeighb[b*maxEpV];
     len3 = PointNeighb[c*maxEpV];

   // find the index of the cells containing current face with point indices a,b,c
    for (j=1;j<=len1;j++)
     {
       Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if (Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            for (l=1;l<=len3;l++)
             if (Neighb_tmp == PointNeighb[c*maxEpV + l])
             {
              Neib[CurrNeib++] = Neighb_tmp;
              break;
             }
           }
          }
       if (CurrNeib == 2) break;
     }
 //   cout<<"CurrNeib " << CurrNeib <<endl;
// cout<<"Facemarkerlist[i] : "<<Facemarkerlist[i]<<endl;
     if (Facemarkerlist[i]) // 0 for inner edges and Boundcomp+1 for Boundedge respect
      {

       CurrComp = Facemarkerlist[i] - 1;
//        cout<<"Boundary face CurrComp: "<<CurrComp<<endl;
// exit(0);
       CurrComp = 0; // not yet implemented fully
       
       
       
       bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);
       

       if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
                             NewVertices[a]->GetY(), T[1], S[1]) ||
          bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
                             NewVertices[b]->GetY(), T[2], S[2]) ||
          bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
                             NewVertices[c]->GetY(), T[3], S[3])    )
         {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          OutPut(NewVertices[c]<<endl);
          exit(0);
         }
 
       if (CurrNeib == 2)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
        }
       else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
// //       cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough neighbours!" << endl;

       Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);

      }

    // First element containing the current face
    // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[0]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[0]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[0]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
     }
      CellTree[Neib[0]]->SetJoint(j, Joint);

   if (Neib[1] != -1) // second element containing the current face
    {
          // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
      }
      CellTree[Neib[1]]->SetJoint(j, Joint);
     }

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();

  }

#ifdef _MPI    
  if(rank==0)    
#endif
  {
//    In.deinitialize();
//    Out.deinitialize();
  }
  

 
 
#ifdef _MPI    
  if(rank!=0) 
   {
    delete [] Tetrahedrals ;  
    delete [] Coordinates;      
    delete [] Facelist;
    delete [] Facemarkerlist;
   }
#endif  
  
  delete [] PointNeighb;

// #ifdef _MPI    
//      if(rank==3)
//     printf(" TetrameshGen  %d \n",   N_G );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//      

}

/** same as TetrameshGen but without using the face list info from tetgen */
void TetraGen(TDomain *&Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, *N_FVerts, N_Faces, *Facets;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PointNeighb, dimension;
  int jj, N_Points, MaxLen, *Tetrahedrals_loc;
  const int *EdgeVertex, *TmpFV, *TmpLen;

  double *Coordinates, N_x, N_y, N_z;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out;
   
  TBaseCell **CellTree,  **SurfCellTree;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices; 
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;  
  TShapeDesc *ShapeDesc;
   
  char *SMESH, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;

#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size); 
#endif    

  /** only the root call the mesh generator and broadcast the mesh info to all */
#ifdef _MPI    
  if(rank==0)
#endif  
   {   
//     In.initialize();
//     Out.initialize();   
    std::ostringstream opts;
    opts << " "; 
  
    opts.seekp(std::ios::beg);
    opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
//     opts<<'r'; // -r  Reconstructs a previously generated mesh     
    opts<<"q"<<1.20; // quality mesh generation(q) with a specified quality bound
//     opts<<"a"<<0.1; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
//     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
//     opts<<'V';  //verbose mode
    opts<<ends;    
  
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);
  
//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);

 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);

 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);    
    
   } // if(rank==0)
     
  /** delete all cells, vertices in the old domain */
  DeleteDomain(Domain);
   
#ifdef _MPI     
     if(rank==0) 
#endif  
     {
       N_RootCells = Out.numberoftetrahedra;
       N_G = Out.numberofpoints;
       // allocate auxillary fields
       Coordinates = Out.pointlist;  
       Tetrahedrals = Out.tetrahedronlist;       
     }
 
#ifdef _MPI
   MPI_Bcast(&N_RootCells, 1, MPI_INT, 0, Comm);
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);     
  if(rank!=0)
   {
    Coordinates = new double [3*N_G];
    Tetrahedrals = new int[4*N_RootCells];
   }   
   MPI_Bcast(Coordinates, 3*N_G, MPI_DOUBLE, 0, Comm);  
   MPI_Bcast(Tetrahedrals, 4*N_RootCells, MPI_INT, 0, Comm);     
#endif      
   
  // generate new vertices
   NewVertices = new TVertex*[N_G];

   for (i=0;i<N_G;i++)
    {
      NewVertices[i] = new TVertex(Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);
      // set bounding box
      if (Coordinates[3*i] > Xmax) Xmax = Coordinates[3*i];
      if (Coordinates[3*i] < Xmin) Xmin = Coordinates[3*i];
      if (Coordinates[3*i+1] > Ymax) Ymax = Coordinates[3*i+1];
      if (Coordinates[3*i+1] < Ymin) Ymin = Coordinates[3*i+1];
      if (Coordinates[3*i+2] > Zmax) Zmax = Coordinates[3*i+2];
      if (Coordinates[3*i+2] < Zmin) Zmin = Coordinates[3*i+2];
    }

   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;

   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);   
//        cout<<Xmin <<"  "<<Ymin <<"  "<<Zmin<<endl;
//        cout<<Xmax <<"  "<<Ymax <<"  "<<Zmax<<endl;
   CellTree = new TBaseCell*[N_RootCells];
 //  output of each tetraheron vertex indices (four vertices for each)
//   for (i=0;i<N_RootCells;i++)
//        cout<< Tetrahedrals[4*i]<<"  "<<Tetrahedrals[4*i + 1]<<"  "
//          <<Tetrahedrals[4*i + 2]<<"  "<<Tetrahedrals[4*i + 3]<<endl;


   for (i=0;i<N_RootCells;i++)
    {
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);

     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);
     CellTree[i]->SetClipBoard(i);
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
   PointNeighb = new int[N_G];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

#ifdef _MPI    
     if(rank==0)   
#endif   
   cout<<"maxEpV "<< maxEpV<<endl;
   PointNeighb = new int[++maxEpV * N_G];
   memset(PointNeighb, 0, maxEpV*N_G*SizeOfInt);
    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of corresponding cells containing this vertex
    // cout<<"maxEpV*N_G "<<maxEpV*N_G<<endl;

   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }
   // generate new faces    
   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();       
   N_G = 0;
   for (i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    ShapeDesc= cell->GetShapeDesc();   
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);      
    N_Faces = cell->GetN_Faces();
    Tetrahedrals_loc = Tetrahedrals+4*i;
           
    for(jj=0;jj<N_Faces;jj++)
     if(cell->GetJoint(jj) == NULL)
      {
       N_G++;    
       N_Points = TmpLen[jj];
   
       if(N_Points!=3)
        {     
         cerr << "Only Tria faces are allowed!!! N_FVert: " << N_Points << endl;
         exit(-1);     
        }     
     
       //printf(" TetrameshGen  Null joint  \n" );             
       a = Tetrahedrals_loc[TmpFV[jj*MaxLen]];
       b = Tetrahedrals_loc[TmpFV[jj*MaxLen + 1]];
       c = Tetrahedrals_loc[TmpFV[jj*MaxLen + 2]];
//      cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;
       Neib[0] = -1;
       Neib[1] = -1;
       CurrNeib = 0;

       len1 = PointNeighb[a*maxEpV];
       len2 = PointNeighb[b*maxEpV];
       len3 = PointNeighb[c*maxEpV];

       // find the index of the cells containing current face with point indices a,b,c
       for (j=1;j<=len1;j++)
       {
        Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if (Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            for (l=1;l<=len3;l++)
             if (Neighb_tmp == PointNeighb[c*maxEpV + l])
             {
              Neib[CurrNeib++] = Neighb_tmp;
              break;
             }
            }
           }
         if (CurrNeib == 2) break;
        }// for (j=1;j<=len1;j++)        
      //   cout<<"CurrNeib " << CurrNeib <<endl;
      // cout<<"Facemarkerlist[i] : "<<Facemarkerlist[i]<<endl;
      if (CurrNeib == 1) // 0 for inner edges and Boundcomp+1 for Boundedge respect
       {
        CurrComp = 0; // not yet implemented fully     
       
        bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);       
 
        if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
                             NewVertices[a]->GetY(), T[1], S[1]) ||
          bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
                             NewVertices[b]->GetY(), T[2], S[2]) ||
          bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
                             NewVertices[c]->GetY(), T[3], S[3])    )
         {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          OutPut(NewVertices[c]<<endl);
          exit(0);
         }
 
       if (CurrNeib == 2)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
        }
       else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
       // cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough neighbours!" << endl;

       Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
      }

     // First element containing the current face
     // find the local index for the point 'a' on the cell
     for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[0]+j] == a) break;

     // find the local index for the point 'b' on the cell
     for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[0]+k] == b) break;

      // find the local index for the point 'c' on the cell
     for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[0]+l] == c) break;   

      l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

      switch (l) // j will contain the local index for the current face
       {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

       default:
        Error("Unable to set the face !!!!!!!!!!!!" << endl);
        exit(0);
      }
      CellTree[Neib[0]]->SetJoint(j, Joint);

    if (Neib[1] != -1) // second element containing the current face
     {
          // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;  
//      cout<<""<< l <<endl;
     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
      }
      CellTree[Neib[1]]->SetJoint(j, Joint);
     } // if (Neib[1] != -1)

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();
   } // if(cell->GetJoint(jj) == NULL)
  } //    for (i=0;i<N_Cells;i++)
  
#ifdef _MPI    
    if(rank==0) 
#endif      
   cout<<"numberoftrifaces after "<<N_G<<endl;
  
#ifdef _MPI    
  if(rank==0)    
#endif
  {
//    In.deinitialize();
//    Out.deinitialize();
  } 
#ifdef _MPI    
  if(rank!=0) 
   {
    delete [] Tetrahedrals ;  
    delete [] Coordinates;      
//     delete [] Facelist;
//     delete [] Facemarkerlist;
   }
#endif   
  delete [] PointNeighb;
// #ifdef _MPI    
//      if(rank==3)
//     printf(" TetrameshGen  %d \n",   N_G );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//      
}

/** input includes tetrahedron */
void TetraGenWithInputCells(TDomain *&Domain)
{
 //======================================================================
 // Tetgen for grid generation begin
 //======================================================================
  int i, j, k, l, N_Coord, *N_FVerts, N_Faces, *Facets;
  int N, N_RootCells, N_Cells, CurrVertex, N_Vertices, ID, N_G, RefLevel=0;
  int CurrNeib, len1, len2, len3, maxEpV = 0, a, b, c, Neib[2], Neighb_tmp, CurrComp;
  int *Tetrahedrals, *PointNeighb, dimension;
  int jj, N_Points, MaxLen, *Tetrahedrals_loc;
  const int *EdgeVertex, *TmpFV, *TmpLen;
  double *TetRegionlist;

  double *Coordinates, N_x, N_y, N_z;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out, InAddPts;
   
  TBaseCell **CellTree,  **SurfCellTree;
  TGridCell **DelCell;
  TVertex **VertexDel, **NewVertices, **NewSurfVertices; 
  TBoundPart *BoundPart;
  TBdPlane **UpdateFaceParams;
  TJoint *Joint;
  TBoundComp3D *bdcomp; 
  TCollection *coll, *SurfColl;
  TBaseCell *cell;
  TBdSphere *UpdateParam;  
  TShapeDesc *ShapeDesc;
   
  char *SMESH, *NODE, line[100];
  SMESH = TDatabase::ParamDB->SMESHFILE;
  NODE  = TDatabase::ParamDB->GRAPEBASENAME;
  
#ifdef _MPI  
  int rank, size;  
  MPI_Comm Comm = TDatabase::ParamDB->Comm;    
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size); 
#endif    

  /** only the root call the mesh generator and broadcast the mesh info to all */
#ifdef _MPI    
  if(rank==0)
#endif  
   {   
//     In.initialize();
//     Out.initialize();   
    std::ostringstream opts;
    opts << " "; 
  
    opts.seekp(std::ios::beg);
//     opts<<'p'; // Tetrahedralize the PLC. Switches are chosen to read a PLC (p)
    opts<<'r'; // -r  Reconstructs a previously generated mesh     
//     opts<<'q'<<1.2; // quality mesh generation(q) with a specified quality bound
//     opts<<'a'<<100; // maximum volume constraint
//     opts<<'a'; // maximum volume constraint
//     opts<<'i'; // Inserts a list of additional points into mesh.
    opts<<'z'; // numbers all output items starting from zero
//     opts<<'d'; // Detect intersections of PLC facets.
    opts<<'f'; // Outputs all  faces (including non-boundary) 
    opts<<'e'; // Outputs a list of edges of the triangulation
//     opts<<'I'; // Suppresses mesh iteration numbers.
    opts<<'C'; // Checks the consistency of the final mesh.
//     opts<<'Q'; // Quiet: No terminal output except errors.
//     opts<<'g'; // Outputs mesh to .mesh file for viewing by Medit
    opts<<'Y'; // Suppresses boundary facets/segments splitting
//     opts<<'A'; // output region attribute list
//     opts<<'V';  //verbose mode
//     opts<<'V';
//     opts<<'V';
//     opts<<'V';
//     opts<<'V';    
    opts<<ends;
    
     In.initialize();
     Out.initialize();   
//      InAddPts.initialize();     
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMeshWithCells(SMESH, In);
//   ReadAdditionalNModes(NODE, InAddPts);   
//      In.load_medit(SMESH, 1);     

//      InAddPts.load_node(NODE);
//      cout<< " InAddPts.numberofpoints " << InAddPts.numberofpoints <<endl; 
//      exit(0);
//     for(i=0;i<In.numberofpoints;i++)
//       OutPut(i<<" (x, y, z) =  "<<
//        In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);
//      exit(0);
 // Calling  tetrahedralize function of 3dtetgen mesh generator
        tetrahedralize((char*)opts.str().c_str(), &In, &Out);
//     tetrahedralize((char*)opts.str().c_str(), &In, &Out, &InAddPts);
//     exit(0);
 //   output: coordinates of all vertices
//  for(i=0;i<Out.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);    
    
   } // if(rank==0)
     
  /** delete all cells, vertices in the old domain */
  DeleteDomain(Domain);
   
#ifdef _MPI     
     if(rank==0) 
#endif  
     {
       N_RootCells = Out.numberoftetrahedra;
       N_G = Out.numberofpoints;   
       // allocate auxillary fields
       Coordinates = Out.pointlist;  
       Tetrahedrals = Out.tetrahedronlist;          
       TetRegionlist = Out.tetrahedronattributelist;
     }

#ifdef _MPI
   MPI_Bcast(&N_RootCells, 1, MPI_INT, 0, Comm);
   MPI_Bcast(&N_G, 1, MPI_INT, 0, Comm);     
   
  if(rank!=0)
   {
    Coordinates = new double [3*N_G];
    Tetrahedrals = new int[4*N_RootCells];
    TetRegionlist = new double [N_RootCells];    
   }
   
   MPI_Bcast(Coordinates, 3*N_G, MPI_DOUBLE, 0, Comm);  
   MPI_Bcast(Tetrahedrals, 4*N_RootCells, MPI_INT, 0, Comm);
   MPI_Bcast(TetRegionlist, N_RootCells, MPI_DOUBLE, 0, Comm);    
#endif      
   
  // generate new vertices
   NewVertices = new TVertex*[N_G];

   for (i=0;i<N_G;i++)
    {
      NewVertices[i] = new TVertex(Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);

      // set bounding box
      if (Coordinates[3*i] > Xmax) Xmax = Coordinates[3*i];
      if (Coordinates[3*i] < Xmin) Xmin = Coordinates[3*i];
      if (Coordinates[3*i+1] > Ymax) Ymax = Coordinates[3*i+1];
      if (Coordinates[3*i+1] < Ymin) Ymin = Coordinates[3*i+1];
      if (Coordinates[3*i+2] > Zmax) Zmax = Coordinates[3*i+2];
      if (Coordinates[3*i+2] < Zmin) Zmin = Coordinates[3*i+2];
    }

   // set bounding box
    StartX = Xmin;
    StartY = Ymin;
    StartZ = Zmin;
    BoundX = Xmax - Xmin;
    BoundY = Ymax - Ymin;
    BoundZ = Zmax - Zmin;
   Domain->SetBoundBox(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);   
//        cout<<Xmin <<"  "<<Ymin <<"  "<<Zmin<<endl;
//        cout<<Xmax <<"  "<<Ymax <<"  "<<Zmax<<endl;
   CellTree = new TBaseCell*[N_RootCells];
 //  output of each tetraheron vertex indices (four vertices for each)
//   for (i=0;i<N_RootCells;i++)
//        cout<< Tetrahedrals[4*i]<<"  "<<Tetrahedrals[4*i + 1]<<"  "
//          <<Tetrahedrals[4*i + 2]<<"  "<<Tetrahedrals[4*i + 3]<<endl;


   for (i=0;i<N_RootCells;i++)
    {
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron], RefLevel);
     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);
     CellTree[i]->SetClipBoard(i);
      ((TMacroCell *) CellTree[i])->SetSubGridID(0);      
//        cout << i<<  "  TetRegionlist :  " <<  TetRegionlist[i] << endl;     
       CellTree[i]->SetRegionID((int)TetRegionlist[i]);
    }
//      exit(0);
   Domain->SetTreeInfo(CellTree, N_RootCells);
   // initialize iterators
   TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
   TDatabase::IteratorDB[It_LE]->SetParam(Domain);
   TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
   TDatabase::IteratorDB[It_Between]->SetParam(Domain);
   TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);   
  
   // search neighbours
   PointNeighb = new int[N_G];
#ifdef _MPI    
     if(rank==0)   
#endif       
   cout<<"numberofpoints "<<N_G<<endl;
   memset(PointNeighb, 0, N_G *SizeOfInt);

     for (i=0;i<4*N_RootCells;i++)
     PointNeighb[Tetrahedrals[i]]++;

   for (i=0;i<N_G;i++)
     if (PointNeighb[i] > maxEpV) maxEpV = PointNeighb[i];
   delete [] PointNeighb;

#ifdef _MPI    
     if(rank==0)   
#endif   
   cout<<"maxEpV "<< maxEpV<<endl;

   PointNeighb = new int[++maxEpV * N_G];

   memset(PointNeighb, 0, maxEpV*N_G*SizeOfInt);
    // every vertex contains "maxEpV" columns
    // for every vertex at first colomn contains the number of cells containing this vertex
    // at further columns we set the index of corresponding cells containing this vertex
    // cout<<"maxEpV*N_G "<<maxEpV*N_G<<endl;

   for(i=0;i<4*N_RootCells;i++)
    {
     j = Tetrahedrals[i]*maxEpV;
     PointNeighb[j]++;
     //cout<<"j + PointNeighb[j] " << j <<endl;
     PointNeighb[j + PointNeighb[j]] = i / 4;
    }

   // generate new faces    
   coll=Domain->GetCollection(It_Finest, 0);
   N_Cells = coll->GetN_Cells();       
   N_G = 0;
   for (i=0;i<N_Cells;i++)
   {
    cell = coll->GetCell(i);
    ShapeDesc= cell->GetShapeDesc();   
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);      
    N_Faces = cell->GetN_Faces();
    Tetrahedrals_loc = Tetrahedrals+4*i;
           
    for(jj=0;jj<N_Faces;jj++)
     if(cell->GetJoint(jj) == NULL)
      {
       N_G++;    
       N_Points = TmpLen[jj];
   
       if(N_Points!=3)
        {     
         cerr << "Only Tria faces are allowed!!! N_FVert: " << N_Points << endl;
         exit(-1);     
        }      
     
       //printf(" TetrameshGen  Null joint  \n" );             
       a = Tetrahedrals_loc[TmpFV[jj*MaxLen]];
       b = Tetrahedrals_loc[TmpFV[jj*MaxLen + 1]];
       c = Tetrahedrals_loc[TmpFV[jj*MaxLen + 2]];

//      cout<<"  "<< a<<"  "<< b<<"  "<< c<<endl;

       Neib[0] = -1;
       Neib[1] = -1;
       CurrNeib = 0;

       len1 = PointNeighb[a*maxEpV];
       len2 = PointNeighb[b*maxEpV];
       len3 = PointNeighb[c*maxEpV];

       // find the index of the cells containing current face with point indices a,b,c
       for (j=1;j<=len1;j++)
       {
        Neighb_tmp = PointNeighb[a*maxEpV + j];
        for (k=1;k<=len2;k++)
         {
          if (Neighb_tmp == PointNeighb[b*maxEpV + k])
           {
            for (l=1;l<=len3;l++)
             if (Neighb_tmp == PointNeighb[c*maxEpV + l])
             {
              Neib[CurrNeib++] = Neighb_tmp;
              break;
             }
            }
           }
         if (CurrNeib == 2) break;
        }// for (j=1;j<=len1;j++)        
      //   cout<<"CurrNeib " << CurrNeib <<endl;
      // cout<<"Facemarkerlist[i] : "<<Facemarkerlist[i]<<endl;
      if (CurrNeib == 1) // 0 for inner edges and Boundcomp+1 for Boundedge respect
       {
        CurrComp = 0; // not yet implemented fully      
       
        bdcomp = Domain->GetBdPart(0)->GetBdComp(CurrComp);       
 
        if(bdcomp->GetTSofXYZ(NewVertices[a]->GetX(), NewVertices[a]->GetY(),
                             NewVertices[a]->GetY(), T[1], S[1]) ||
          bdcomp->GetTSofXYZ(NewVertices[b]->GetX(), NewVertices[b]->GetY(),
                             NewVertices[b]->GetY(), T[2], S[2]) ||
          bdcomp->GetTSofXYZ(NewVertices[c]->GetX(), NewVertices[c]->GetY(),
                             NewVertices[c]->GetY(), T[3], S[3])    )
         {
          cerr<<"Error: could not set parameter values"<<endl;
          OutPut(NewVertices[a]<<endl);
          OutPut(NewVertices[b]<<endl);
          OutPut(NewVertices[c]<<endl);
          exit(0);
         }
 
       if (CurrNeib == 2)
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoInterfaceJoint3D(bdcomp, T, S,
                       CellTree[Neib[0]], CellTree[Neib[1]] );
          else
           Joint = new TInterfaceJoint3D(bdcomp, T, S, 
                       CellTree[Neib[0]], CellTree[Neib[1]] );
        }
       else
        {
         if(bdcomp->IsFreeBoundary())
           Joint = new TIsoBoundFace(bdcomp, T, S);
         else
           Joint = new TBoundFace(bdcomp, T, S);
        }
      }
     else
      {
       // cout<<"Inner face"<<endl;
       if (CurrNeib != 2)
       cerr << "Error !!!!!!!! not enough neighbours!" << endl;

       Joint = new TJointEqN(CellTree[Neib[0]], CellTree[Neib[1]]);
      }

     // First element containing the current face
     // find the local index for the point 'a' on the cell
     for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[0]+j] == a) break;

     // find the local index for the point 'b' on the cell
     for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[0]+k] == b) break;

      // find the local index for the point 'c' on the cell
     for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[0]+l] == c) break;   

      l = l*100 + k*10 + j;  
//      cout<<""<< l <<endl;

      switch (l) // j will contain the local index for the current face
       {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

       default:
        Error("Unable to set the face !!!!!!!!!!!!" << endl);
        exit(0);
      }
      CellTree[Neib[0]]->SetJoint(j, Joint);

    if (Neib[1] != -1) // second element containing the current face
     {
          // find the local index for the point 'a' on the cell
    for (j=0;j<4;j++)
      if (Tetrahedrals[4*Neib[1]+j] == a) break;

    // find the local index for the point 'b' on the cell
    for (k=0;k<4;k++)
      if (Tetrahedrals[4*Neib[1]+k] == b) break;

       // find the local index for the point 'c' on the cell
    for (l=0;l<4;l++)
      if (Tetrahedrals[4*Neib[1]+l] == c) break;   

     l = l*100 + k*10 + j;  

//      cout<<""<< l <<endl;

     switch (l) // j will contain the local index for the current face
      {
        case 210: case 21: case 102:
        case 120: case 12: case 201:
          j = 0;
          break;  
        case 310: case 31: case 103:
        case 130: case 13: case 301:
          j = 1;
          break;  
        case 321: case 132: case 213:
        case 231: case 123: case 312:
          j = 2;
          break;  
        case 230: case 23: case 302:
        case 320: case 32: case 203:
          j = 3;
          break; 

      default:
       Error("Unable to set the face !!!!!!!!!!!!" << endl);
       exit(0);
      }
      CellTree[Neib[1]]->SetJoint(j, Joint);
     } // if (Neib[1] != -1)

  if (Joint->GetType() == InterfaceJoint3D ||
      Joint->GetType() == IsoInterfaceJoint3D)
      {
        ((TInterfaceJoint3D*)Joint)->SetMapType();
        ((TInterfaceJoint3D*)(Joint))->CheckOrientation();
      }
      else 
       if (Joint->GetType() == JointEqN)
           Joint->SetMapType();
   } // if(cell->GetJoint(jj) == NULL)
  } //    for (i=0;i<N_Cells;i++)
  
#ifdef _MPI    
    if(rank==0) 
#endif      
   cout<<"numberoftrifaces after "<<N_G<<endl;
  
#ifdef _MPI    
  if(rank==0)    
#endif
  {
//    In.deinitialize();
//    Out.deinitialize();
  }
 
#ifdef _MPI    
  if(rank!=0) 
   {
    delete [] Tetrahedrals ;  
    delete [] Coordinates;      
    delete [] TetRegionlist;
//     delete [] Facemarkerlist;
   }
#endif  
  
  delete [] PointNeighb;

// #ifdef _MPI    
//      if(rank==3)
//     printf(" TetrameshGen  %d \n",   N_G );  
//     MPI_Finalize(); 
// #endif    
//   exit(0);
//      
}
