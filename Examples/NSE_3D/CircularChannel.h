// channel with circular cross section



#include <InterfaceJoint3D.h>
#include <IsoInterfaceJoint3D.h>
#include <IsoBoundFace.h>
#include <MacroCell.h>
#include <BdSphere.h>
#include <tetgen.h>


void ExampleFile()
{
  OutPut("Example: CircularChannel.h" << endl) ;
  
  #define __Cylinder__   
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
 
  values[0] = 1-(x*x+y*y);
  values[1] = -2*x;
  values[2] = -2*y;
  values[3] = 0;
  values[4] = -4;
}

void ExactP(double x, double y,  double z, double *values)
{
  static double eps = 1/TDatabase::ParamDB->RE_NR;

  values[0] = 4*eps*(10-z);
  values[1] = 0;
  values[2] = 0;
  values[3] = -4*eps;
  values[4] = 0;
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int CompID, double x, double y, double z, BoundCond &cond)
{
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  if (fabs(z-10)<1e-6)
  {
    cond = NEUMANN;
    TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE = 0;
  }
  else
    cond  = DIRICHLET;
}


// value of boundary condition
void U2BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U1BoundValue(int CompID, double x, double y, double z, double &value)
{
  value = 0;
}

// value of boundary condition
void U3BoundValue(int CompID, double x, double y, double z, double &value)
{
   double r2;
   if (fabs(z)<1e-8)
   {
      r2 = x*x+y*y;
      value =  1-r2;
   }
   else
   {
      value = 0;
   }
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y, double *Z,
               double **parameters, double **coeffs)
{
  static double eps = 1./TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y, z;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
      
    coeff[0] = eps;
    coeff[1] = 0;
    coeff[2] = 0;
    coeff[3] = 0;
  }
}

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
  In.pointmtrlist = new double[N_Vertices];  
  In.numberofpointmtrs = 1;
  
  //read from file
   for(i=0;i<N_Vertices; i++)
    {
     dat.getline (line, 99);
     dat >> In.pointlist[3*i] >> In.pointlist[3*i+1] >> In.pointlist[3*i+2];
     In.pointmtrlist[i] = 10.;
     //cout<< i << " vert X: " <<In.pointlist[3*i] << " vert Y: " <<In.pointlist[3*i+1] <<endl;   
    }
 
//       In.pointmtrlist[0] = 0.75;
 
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

  double X, Y, Z, DispX, DispY, DispZ;  
  double *Coordinates, N_x, N_y, N_z;   
  double *Vertices;
  double Xmin = 1e10, Xmax = -1e10, Ymin = 1e10, Ymax = -1e10;  
  double Zmin = 1e10, Zmax = -1e10, T[4]={0,0,0,0}, S[4]={0,0,0,0}; 
  double StartX, StartY, StartZ, BoundX, BoundY, BoundZ; 
  
  tetgenio In, Out, Out_New;
   
  TBaseCell **CellTree;
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
    opts<<"q"<<1.5; // quality mesh generation(q) with a specified quality bound
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
//     opts<<'Y'; // Suppresses boundary facets/segments splitting
//     opts<<'m'; //  
   
    
    opts<<'V';  //verbose mode
    opts<<ends;
   
//   cout << " SMESHFILE is " << SMESH << endl;
  /** load the medit mesh file into Tetgen */  
  ReadMeditMesh(SMESH, In);
 
//    for(i=0;i<In.numberofpoints;i++)
//   OutPut(" (x, y, z) =  "<< In.pointlist[3*i]<<' '<<In.pointlist[3*i+1]<<' '<<In.pointlist[3*i+2]<<endl);
// exit(0);
 // Calling  tetrahedralize function of 3dtetgen mesh generator
    tetrahedralize((char*)opts.str().c_str(), &In, &Out);

   } // if(rank==0)
  
//   MPI_Finalize();
//   exit(0);
//   
    
 //   output: coordinates of all vertices
 for(i=0;i<Out.numberofpoints;i++)
  OutPut(" (x, y, z) =  "<< Out.pointlist[3*i]<<' '<<Out.pointlist[3*i+1]<<' '<<Out.pointlist[3*i+2]<<endl);
exit(0);

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
     NewVertices[i] = new TVertex( Coordinates[3*i], Coordinates[3*i+1], Coordinates[3*i+2]);
     
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
     CellTree[i] = new TMacroCell(TDatabase::RefDescDB[Tetrahedron],
                                    RefLevel);

     CellTree[i]->SetVertex(0, NewVertices[Tetrahedrals[4*i    ]]);
     CellTree[i]->SetVertex(1, NewVertices[Tetrahedrals[4*i + 1]]);
     CellTree[i]->SetVertex(2, NewVertices[Tetrahedrals[4*i + 2]]);
     CellTree[i]->SetVertex(3, NewVertices[Tetrahedrals[4*i + 3]]);

    CellTree[i]->SetClipBoard(i);
     ((TMacroCell *) CellTree[i])->SetSubGridID(0);
     
     //default all cells will be treated as skull
     CellTree[i]->SetRegionID(1);
     CellTree[i]->SetAsLayerCell(1);     
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



