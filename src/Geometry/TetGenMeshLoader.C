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
// @(#)Domain.C        1.23 05/05/00
//
// Class:       TTetGenMeshLoader
// Purpose:     creates domain with tetgen from a .smesh file
//
// Author:      Andreas Hahn 26.04.2010
//
// History:
//
// =======================================================================


#include <BoundPart.h>
#include <BoundComp3D.h>
#include <BdPlane.h>
#include <BdCylinder.h>
#include <Vertex.h>
#include <MacroCell.h>
#include <Database.h>
#include <JointEqN.h>
#include <BoundFace.h>
#include <MooNMD_Io.h>
#include <TetGenMeshLoader.h>

#include <sstream>
#include <string>

// CTOR
TTetGenMeshLoader::TTetGenMeshLoader()
{
  mFileName      = NULL;
  mBoundFile     = NULL;
  mAllBoundComps = NULL;
  mAllVertices   = NULL;
  mAllJoints     = NULL;
  mTrifaceHash   = NULL;
  plc     	= false;
  reconstruct   = false;
  insertpoints  = false;
}

TTetGenMeshLoader::TTetGenMeshLoader(const char *FileName)
{
  int n, status;

  // count string length
  n=0;
  while (FileName[n] != 0) {++n;}

  mFileName = new char [n+1];
  mBoundFile = new char [n+5];

  // copy filename
  for (int i=0;i<n;++i)
  {
    mFileName[i] = FileName[i];
    mBoundFile[i] = FileName[i];
  }
  mFileName[n] = 0;

  mAllBoundComps = NULL;
  mAllVertices   = NULL;
  mAllJoints     = NULL;
  mTrifaceHash   = NULL;

  plc     	= false;
  reconstruct   = false;
  insertpoints  = false;
}


// DTOR
TTetGenMeshLoader::~TTetGenMeshLoader()
{
  if ( mFileName )  delete [] mFileName;
  if ( mBoundFile )   delete [] mBoundFile;
  if ( mAllBoundComps ) delete [] mAllBoundComps;
  if ( mAllVertices )   delete [] mAllVertices;
  if ( mAllJoints )   delete [] mAllJoints;

  if ( mTrifaceHash )
  {
    for (int i=0;i<3*mTetOut.numberofpoints;++i)
    {
      if ( mTrifaceHash[i] != NULL )
        delete [] mTrifaceHash[i];
    }

    delete [] mTrifaceHash;
  }
}


// Methods

int TTetGenMeshLoader::BuildBoundary(TBoundPart **&BdParts, int &N_BoundParts,
int &N_BoundComps, int *&StartBdCompID,
int *&Interfaces)
{
  //   if ( mFileName == NULL)
  //   {
  //     cerr << "No file name given !" << endl;
  //     return -1;
  //   }
  //
  //   // read smesh file
  //   mTetIn.load_poly(mFileName);

  // count number of boundary components
  N_BoundComps = 0;
  for (int i=0;i<mTetIn.numberoffacets;++i)
  {
    if ( N_BoundComps < mTetIn.facetmarkerlist[i] ) N_BoundComps = mTetIn.facetmarkerlist[i];
  }

  cout << "N_BoundParts: " << (N_BoundParts = mTetIn.numberofholes + 1) << endl;
  cout << "N_BoundComps: " << N_BoundComps << endl;

  if  (N_BoundParts != 1)
  {
    cerr << "Domains with holes not supportet !" << endl;
    return -1;
  }

  // alloc arrays for parts and components
  BdParts        = new TBoundPart* [N_BoundParts];
  Interfaces     = new int [N_BoundParts];
  mAllBoundComps = new TBoundComp3D* [N_BoundComps];
  StartBdCompID  = new int [N_BoundParts+1];

  ReadBoundFile(N_BoundComps);

  SetBdPlaneParams(N_BoundComps);

  // set components to part
  BdParts[0] = new TBoundPart(N_BoundComps);

  for (int i=0;i<N_BoundComps;++i)
  {
    BdParts[0]->SetBdComp(i, mAllBoundComps[i]);
  }

  StartBdCompID[0] = 0;
  StartBdCompID[1] = N_BoundComps;

  Interfaces[0] = 1;

  return 0;
}


int TTetGenMeshLoader::AllocRootCells(TBaseCell **&CellTree, int &N_RootCells)
{
  TMacroCell *Cell;
  TVertex *Vertex;
  int index, orient;

  N_RootCells = mTetOut.numberoftetrahedra;

  CellTree = new TBaseCell* [N_RootCells];

//   OutPut("Number of region: " << mTetIn.numberofregions << endl);
  OutPut("# of tetrahedron attributes: " << mTetOut.numberoftetrahedronattributes << endl);
  
  for (int i=0;i<N_RootCells;++i)
  {
    Cell = new TMacroCell (TDatabase::RefDescDB[Tetrahedron], 0);
    CellTree[i] = Cell;

    if ( mTetOut.numberofcorners != 4 )
    {
      cout << "Wrong number of corners !" << endl;
      return -1;
    }

    index = 4*i;

    for (int j=0;j<4;++j)
    {
      //       cout << "j: " << j << endl;
      Vertex = mAllVertices[mTetOut.tetrahedronlist[index+j]];

      Cell->SetVertex(j, Vertex);
    }
    
    if ( mTetOut.numberoftetrahedronattributes > 0 )
    {
      Cell->SetRegionID((int) mTetOut.tetrahedronattributelist[i]);
//       OutPut("tetrahedron attribute: " << Cell->GetRegionID() << endl);
    }
  }

  OutPut("number of tetrahedra: " << mTetOut.numberoftetrahedra << endl);

  return 0;
}


int TTetGenMeshLoader::AllocVertices(double &StartX, double &StartY, double &StartZ,
double &BoundX, double &BoundY, double &BoundZ)
{
  double x, y, z;
  double xmin, xmax, ymin, ymax, zmin, zmax;

  mAllVertices = new TVertex* [mTetOut.numberofpoints];

  for (int i=0;i<mTetOut.numberofpoints;++i)
  {
    x = mTetOut.pointlist[3*i  ];
    y = mTetOut.pointlist[3*i+1];
    z = mTetOut.pointlist[3*i+2];

    if ( i > 0 )
    {
      if ( xmin > x ) xmin = x;
      if ( xmax < x ) xmax = x;

      if ( ymin > y ) ymin = y;
      if ( ymax < y ) ymax = y;

      if ( zmin > z ) zmin = z;
      if ( zmax < z ) zmax = z;
    }
    else
    {
      xmin = xmax = x;
      ymin = ymax = y;
      zmin = zmin = z;
    }

    mAllVertices[i] = new TVertex (x, y, z);
  }

  StartX = xmin;
  BoundX = xmax - xmin;

  StartY = ymin;
  BoundY = ymax - ymin;

  StartZ = zmin;
  BoundZ = zmax - zmin;

  OutPut("number of vertices: " << mTetOut.numberofpoints << endl);

  return 0;
}


int TTetGenMeshLoader::SetBdPlaneParams(int N_BoundComps)
{
  TBoundComp3D *BoundComp;

  for (int i=0;i<N_BoundComps;++i)
  {
    BoundComp = mAllBoundComps[i];

    if ( BoundComp->GetType() == Plane )
    {
      double p[3], a[3], b[3], n[3];
      int tetVert1, tetVert2, tetVert3;
      tetgenio::facet *tetFacet = NULL;

      // find facet belonging to boundcomp
      for (int j=0;j<mTetIn.numberoffacets;++j)
      {
        if ( mTetIn.facetmarkerlist[j] == i+1 )
        {
          tetFacet = mTetIn.facetlist + j;
          break;
        }
      }

      if ( tetFacet == NULL )
      {
        cerr << "could not find facet!" << endl;
        return -1;
      }

      // check vertex count
      if ( tetFacet->polygonlist->numberofvertices < 3 )
      {
        cout << "Facet has less than 3 vertices !" << endl;
        return -1;
      }

      // get first three vertices
      tetVert1 = tetFacet->polygonlist->vertexlist[0];
      tetVert2 = tetFacet->polygonlist->vertexlist[1];
      tetVert3 = tetFacet->polygonlist->vertexlist[2];

      // calc vectors for plane
      for (int k=0;k<3;++k)
      {
        p[k] = mTetIn.pointlist[3*tetVert1+k];
        a[k] = mTetIn.pointlist[3*tetVert2+k] - p[k];
        b[k] = mTetIn.pointlist[3*tetVert3+k] - p[k];
      }

      Normalize3D(a);
      Normalize3D(b);
      Cross3D(b, a, n);

      ((TBdPlane*) BoundComp)->SetParams(p[0], p[1], p[2], a[0], a[1], a[2],
        n[0], n[1], n[2]);
    }
  }

  return 0;
}


void TTetGenMeshLoader::Normalize3D(double *vec)
{
  double abs;

  abs = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);

  vec[0] /= abs;
  vec[1] /= abs;
  vec[2] /= abs;
}


void TTetGenMeshLoader::Cross3D(double *vec1, double *vec2, double *res)
{
  res[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  res[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  res[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  Normalize3D(res);
}


int TTetGenMeshLoader::ReadBoundComp(std::ifstream *in, int &CurrBoundComp, int &Range)
{
  char c;

  *in >> CurrBoundComp;
  *in >> c;

  // check if range is given
  if ( c == '-' )
  {
    *in >> Range;
  }
  else
  {
    in->unget();
    Range = CurrBoundComp;
  }

  //   cout << CurrBoundComp << " - " << Range << "  ";

  --CurrBoundComp;

  return 0;
}


int TTetGenMeshLoader::ReadBoundType(std::ifstream *in, char *buff)
{
  *in >> buff;

  //   cout << buff << endl;

  return 0;
}


int TTetGenMeshLoader::ReadBoundParams(std::ifstream *in, char *type,
int CurrBoundComp, int Range)
{
  if ( strcmp(type, "Plane") == 0 )
  {
    for (int i=CurrBoundComp;i<Range;++i)
      mAllBoundComps[i] = new TBdPlane (i);
  }
  else if ( strcmp(type, "Cylinder") == 0 )
  {
    double radius, px, py, pz, ax, ay, az, nx, ny, nz;

    ReadBoundParamsCylinder(in, radius, px, py, pz, ax, ay, az,
      nx, ny, nz);

    for (int i=CurrBoundComp;i<Range;++i)
    {
      mAllBoundComps[i] = new TBdCylinder (i);
      ((TBdCylinder*)mAllBoundComps[i])->SetParams(radius, px, py, pz,
        ax, ay, az, nx, ny, nz);
    }
  }
  else
  {
    cerr << "Boundary type \"" << type << "\" not yet supported !" << endl;
    return -1;
  }
    return(0);
}


int TTetGenMeshLoader::ExtractValue(const char *buff, int increment, double &val)
{
  std::istringstream str(buff+increment);

  str >> val;
  //   cout << val << endl;

  return 0;
}


int TTetGenMeshLoader::ExtractVector(const char *buff, int increment,
double &val1, double &val2, double &val3)
{
  char buff1[128], c;
  const char *p;

  p = buff + increment + 1;

  // remove brakets and komma
  for (int i=0; (buff1[i]=0, c=p[i]) != ')'; ++i)
  {
    if ( c != ',' )
      buff1[i] = c;

    else
      buff1[i] = ' ';
  }

  //   cout << buff1 << endl;

  std::istringstream str(buff1);

  str >> val1 >> val2 >> val3;

  return 0;
}


int TTetGenMeshLoader::ReadBoundParamsCylinder(std::ifstream *in, double &radius,
double &px, double &py, double &pz,
double &ax, double &ay, double &az,
double &nx, double &ny, double &nz)
{
  int eqsign;
  char buff[128];

  // four parameters
  for (int i=0;i<4;++i)
  {
    *in >> buff;

    // find equality sign
    for (int j=0; j<128;++j)
    {
      if ( buff[j] == '=' )
      {
        eqsign = j;
        break;
      }
    }

    // switch params
    if ( strncmp(buff, "radius", 6) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractValue(buff, eqsign+1, radius);
      //       cout << "radius: " << radius << endl;
    }
    else if ( strncmp(buff, "zero", 4) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractVector(buff, eqsign+1, px, py, pz);
      //       cout << "px: " << px << " , py: " << py << " , pz: " << pz << endl;
    }
    else if ( strncmp(buff, "axis", 4) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractVector(buff, eqsign+1, ax, ay, az);
      //       cout << "ax: " << ax << " , ay: " << ay << " , az: " << az << endl;
    }
    else if ( strncmp(buff, "normal", 6) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractVector(buff, eqsign+1, nx, ny, nz);
      //       cout << "nx: " << nx << " , ny: " << ny << " , nz: " << nz << endl;
    }
    else
    {
      cerr << "parameter unexpected for cylinder !" << endl;
      return -1;
    }
  }

  return 0;
}

int TTetGenMeshLoader::ReadBoundParamsSphere(std::ifstream *in, double &radius, 
					     double &mx, double &my, double &mz)
{
  int eqsign;
  char buff[128];

  // four parameters
  for (int i=0;i<4;++i)
  {
    *in >> buff;

    // find equality sign
    for (int j=0; j<128;++j)
    {
      if ( buff[j] == '=' )
      {
        eqsign = j;
        break;
      }
    }

    // switch params
    if ( strncmp(buff, "radius", 6) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractValue(buff, eqsign+1, radius);
      //       cout << "radius: " << radius << endl;
    }
    else if ( strncmp(buff, "midpoint", 8) == 0 )
    {
      //       cout << endl << buff << endl;
      ExtractVector(buff, eqsign+1, mx, my, mz);
      //       cout << "px: " << px << " , py: " << py << " , pz: " << pz << endl;
    }
    else
    {
      cerr << "parameter unexpected for cylinder !" << endl;
      return -1;
    }
  }

  return 0;
}

int TTetGenMeshLoader::ReadBoundFile(int N_BoundComps)
{
  char postfix[] = ".bound";
  int CurrBoundComp, Range;

  char buff[128];

  strcat (mBoundFile, postfix);
  //   cout << mBoundFile << endl;

  // open bound file
  std::ifstream in (mBoundFile);
  if ( !in.is_open() )
  {
    cout << "Warning: Could not find \"" << mBoundFile << "\"! " << endl;
    cout << " -- All boundary components are assumed to be planes!" << endl;

    for (int i=0;i<N_BoundComps;++i)
    {
      mAllBoundComps[i] = new TBdPlane (i);
    }

    return 0;
  }

  // read bound file
  while ( !in.eof() )
  {
    ReadBoundComp(&in, CurrBoundComp, Range);
    ReadBoundType(&in, buff);
    ReadBoundParams(&in, buff, CurrBoundComp, Range);
  }

  in.close();
  
  return(0);
}


void TTetGenMeshLoader::MakeOptionString(std::string &opts)
{
  std::ostringstream osq, osa;

  osq << TDatabase::ParamDB->TETGEN_QUALITY;
  osa << TDatabase::ParamDB->TETGEN_VOLUMEN;

  if ( plc )
  {
    opts += "p";                                  // Tetrahedralizes a piecewise linear comple.
  }

  if ( reconstruct )
  {
    opts += "r";                                  // reconstruct mesh
  }

  if ( TDatabase::ParamDB->TETGEN_QUALITY >= 1.0 /*&& plc*/ )
  {
    opts += "q";
    opts += osq.str();
  }

  if ( TDatabase::ParamDB->TETGEN_VOLUMEN > 0.0 /*&& plc*/ )
  {
    opts += "a";
    opts += osa.str();
  }
  
  if ( insertpoints )
  {
    opts += "i";
  }

  opts += "z";                                    // Numbers all output items starting from zero.
  opts += "f";                                    // Outputs faces (including non-boundary faces).

  if ( TDatabase::ParamDB->TETGEN_STEINER == 1 && plc )
  {
    opts += "Y";                                  // Suppresses boundary facets/segments splitting.
  }

  opts += "C";                                    // Checks the consistency of the final mesh.
  //   opts += "n"; // Outputs tetrahedra neighbors.
//   opts += "Q";                                 // Quiet
//   opts += "V"; 					  // Verbose: Detailed information, more terminal output.
  opts += "M";                                    // Does not merge coplanar facets.
  
  opts += "A"; // assigns region attributes
}


int TTetGenMeshLoader::AllocJoints(TBaseCell **CellTree)
{
  int left, right;
  int bdcomp;
  double param1[4], param2[4], x, y, z;

  TBoundComp3D *BoundComp;

  mAllJoints = new TJoint* [mTetOut.numberoftrifaces];

  for (int i=0;i<mTetOut.numberoftrifaces;++i)
  {
    if ( mTetOut.trifacemarkerlist[i] == 0 )      // inner joints
    {
      //       GetTetrasOfTriface(i, left, right);
      left = mTetOut.adjtetlist[2*i];
      right = mTetOut.adjtetlist[2*i+1];

      mAllJoints[i] = new TJointEqN (CellTree[left],
        CellTree[right]);
    }
    else                                          // boundary joints
    {
      //       cout << "boundary joint" << endl;

      bdcomp = mTetOut.trifacemarkerlist[i] - 1;

      BoundComp = mAllBoundComps[bdcomp];

      //       JointIsFaceOfTet(i, left, right);

      //       cout << "\t" << left << "   " << right << endl;

      //       JointIsFaceOfTet(i, left, right);

      //       for (int j=0; j<3; ++j)
      //       {
      // 	x = mTetOut.pointlist[3*mTetOut.trifacelist[3*i+j]   ];
      // 	y = mTetOut.pointlist[3*mTetOut.trifacelist[3*i+j] +1];
      // 	z = mTetOut.pointlist[3*mTetOut.trifacelist[3*i+j] +2];
      //
      // 	BoundComp->GetTSofXYZ(x, y, z, param1[j], param2[j]);
      //       }

      mAllJoints[i] = new TBoundFace (BoundComp/*, param1, param2*/);
      //       mAllJoints[i]->SetNeighbour (0, CellTree[left] );
    }
  }

  OutPut("number of joints: " << mTetOut.numberoftrifaces << endl);

  return 0;
}


int TTetGenMeshLoader::Load_msh(const char *filename)
{
  char buff[256];
  int N_, k, l;
  REAL x, y, z;

  std::ifstream dat(filename);

  if (!dat)
  {
    cerr << "could not open file \"" << filename << "\"" << endl;
    return 1;
  }

  dat.getline(buff, 255);
  dat.getline(buff, 255);
  dat.getline(buff, 255);
  dat.getline(buff, 255);

  dat >> N_;
  cout << "numberofpoints: " << N_ << endl;

  mTetIn.mesh_dim = 3;
  mTetIn.numberofpoints = N_;
  mTetIn.numberofpointattributes = 0;
  mTetIn.numberofpointmtrs = 0;

  mTetIn.pointlist = new REAL [3*N_];
  //   mTetIn.pointmarkerlist = new int [N_];

  for (int i=0;i<N_;++i)
  {
    dat >> k >> x >> y >> z;
    mTetIn.pointlist[3*i  ] = x;
    mTetIn.pointlist[3*i+1] = y;
    mTetIn.pointlist[3*i+2] = z;
    //     cout << k << " " << x << " " << y << " " << z << endl;
  }

  dat.getline(buff, 255);
  dat.getline(buff, 255);
  dat.getline(buff, 255);

  dat >> N_;
  cout << "number of elements: " << N_ << endl;

  mTetIn.numberoftetrahedra = 0;
  mTetIn.numberofcorners = 4;
  mTetIn.numberoftetrahedronattributes = 0;

  mTetIn.tetrahedronlist = new int [4*N_];

  int elm_number, elm_type, n_tags;
  for (int i=0;i<N_;++i)
  {
    dat >> elm_number >> elm_type >> n_tags;
    if ( elm_type == 4 )
    {
      for (int j=0;j<n_tags;++j)
      {
        dat >> k;
      }

      k = mTetIn.numberoftetrahedra++;
      dat >> l;
      mTetIn.tetrahedronlist[4*k  ] = l-1;
      dat >> l;
      mTetIn.tetrahedronlist[4*k+1] = l-1;
      dat >> l;
      mTetIn.tetrahedronlist[4*k+2] = l-1;
      dat >> l;
      mTetIn.tetrahedronlist[4*k+3] = l-1;
    }

    dat.getline(buff, 255);
  }

  cout << "numberoftetrahedra: " << mTetIn.numberoftetrahedra << endl;

  dat.close();

  return 0;
}


int TTetGenMeshLoader::Generate(TDomain *Domain)
{
  char *suffix;

  suffix = strrchr(mFileName, '.');

  if ( !strcmp(suffix, ".smesh") )
  {
    *suffix = 0;

    suffix = strrchr(mBoundFile, '.');
    *suffix = 0;

    plc = true;

    // read smesh file
    mTetIn.load_poly(mFileName);

    if ( mTetIn.facetmarkerlist != NULL )
    {
      BuildBoundary(Domain->BdParts, Domain->N_BoundParts,
        Domain->N_BoundComps, Domain->StartBdCompID, Domain->Interfaces);

      Tetgen();

      BuildMooNMDMesh(Domain->CellTree, Domain->N_RootCells,
        Domain->StartX, Domain->StartY, Domain->StartZ,
        Domain->BoundX, Domain->BoundY, Domain->BoundZ);
    }
    else
    {
      OutPut("NO MARKERS" << endl);

//       MakeBoundaryLayer_smesh();
      
      Tetgen();

      BuildBoundary2(Domain->BdParts, Domain->N_BoundParts,
        Domain->N_BoundComps, Domain->StartBdCompID, Domain->Interfaces);

      BuildMooNMDMesh (Domain->CellTree, Domain->N_RootCells,
        Domain->StartX, Domain->StartY, Domain->StartZ,
        Domain->BoundX, Domain->BoundY, Domain->BoundZ);
    }

    return 1;
  }

  if ( !strcmp(suffix, ".msh") )
  {
    OutPut("gmsh" << endl);

    reconstruct = true;

    Load_msh(mFileName);
    
    MakeBoundaryLayer_msh();

    Tetgen();

    BuildBoundary2(Domain->BdParts, Domain->N_BoundParts,
      Domain->N_BoundComps, Domain->StartBdCompID, Domain->Interfaces);

    BuildMooNMDMesh (Domain->CellTree, Domain->N_RootCells,
      Domain->StartX, Domain->StartY, Domain->StartZ,
      Domain->BoundX, Domain->BoundY, Domain->BoundZ);

    return 1;
  }

  return 0;
}


int  TTetGenMeshLoader::FindTriface(int a, int b, int c)
{
  int hash = a+b+c;
  int count ;
  int *trifacelist = mTetOut.trifacelist;
  int triface, vertex, found;

  assert(mTrifaceHash[hash] != NULL);
  count = mTrifaceHash[hash][0];

  for (int i=1;i<=count;++i)
  {
    triface = mTrifaceHash[hash][i];
    found = 0;

    for (int j=0;j<3;++j)
    {
      vertex = trifacelist[3*triface+j];

      if ( a == vertex ||
        b == vertex ||
        c == vertex   )
        ++found;
    }

    if ( found == 3)
      return triface;
  }

  return -1;
}


int TTetGenMeshLoader::CreateAdjacency()
{
  int n_trifaces = mTetOut.numberoftrifaces;
  int n_tetrahedra = mTetOut.numberoftetrahedra;
  int triface;
  int N_BoundComp = n_trifaces;

  int FaceVertex[][3] =
  {
    {
      0, 1, 2
    }
    ,
    {
      0, 3, 1
    }
    ,
    {
      2, 1, 3
    }
    ,
    {
      0, 2, 3
    }
  };

  if ( mTetOut.adjtetlist ) delete [] mTetOut.adjtetlist;

  mTetOut.adjtetlist = new int [2*n_trifaces];
  for (int i=0;i<2*n_trifaces;++i)
    mTetOut.adjtetlist[i] = -1;

  for (int i=0;i<n_tetrahedra;++i)
  {
    for (int face=0;face<4;++face)
    {
      triface = FindTriface(mTetOut.tetrahedronlist[4*i+FaceVertex[face][0]],
        mTetOut.tetrahedronlist[4*i+FaceVertex[face][1]],
        mTetOut.tetrahedronlist[4*i+FaceVertex[face][2]]);

      assert ( triface != -1 );

      if (mTetOut.adjtetlist[2*triface] == -1)
        mTetOut.adjtetlist[2*triface] = i;
      else
      {
        assert(mTetOut.adjtetlist[2*triface+1] == -1);
        mTetOut.adjtetlist[2*triface+1] = i;
        --N_BoundComp;
      }
    }
  }

  return N_BoundComp;
}


void TTetGenMeshLoader::HashTrifaces()
{
  int left, right, N_BoundComp=0;
  int a, b, c, hash;
  int n_trifaces = mTetOut.numberoftrifaces;
  int n_points = mTetOut.numberofpoints;

  mTrifaceHash = new int* [3*n_points];
  memset(mTrifaceHash, 0, 3*n_points*sizeof(mTrifaceHash[0]));

  int *BucketCount = new int [3*n_points];
  memset(BucketCount, 0, 3*n_points*sizeof(BucketCount[0]));

  for (int i=0;i<mTetOut.numberoftrifaces;++i)
  {
    a = mTetOut.trifacelist[3*i  ];
    b = mTetOut.trifacelist[3*i+1];
    c = mTetOut.trifacelist[3*i+2];

    hash = a+b+c;
    BucketCount[hash]++;
  }

  for (int i=0;i<mTetOut.numberoftrifaces;++i)
  {
    a = mTetOut.trifacelist[3*i  ];
    b = mTetOut.trifacelist[3*i+1];
    c = mTetOut.trifacelist[3*i+2];

    hash = a+b+c;

    if ( mTrifaceHash[hash] == NULL )
    {
      mTrifaceHash[hash] = new int [BucketCount[hash]+1];
      mTrifaceHash[hash][0] = 1;
      mTrifaceHash[hash][1] = i;
    }
    else
    {
      int pos = ++mTrifaceHash[hash][0];
      mTrifaceHash[hash][pos] = i;
    }
  }

  //   for (int i=0;i<n_points;++i)
  //     cout << BucketCount[i] << " ";

  //   cout << endl;

  delete [] BucketCount;
}


int TTetGenMeshLoader::BuildBoundary2 (TBoundPart **&BdParts, int &N_BoundParts,
int &N_BoundComps, int *&StartBdCompID, int *&Interfaces)
{
  int n_trifaces = mTetOut.numberoftrifaces;
  int *trifaces = mTetOut.trifacelist;
  int counter=0;
  double p[3], a[3], b[3], n[3];

  HashTrifaces();
  N_BoundComps = CreateAdjacency();
  OutPut ("N_BoundComps: " << N_BoundComps << endl);

  N_BoundParts = 1;
  BdParts = new TBoundPart* [N_BoundParts];
  StartBdCompID = new int [N_BoundParts+1];

  StartBdCompID[0] = 0;
  StartBdCompID[1] = N_BoundComps;

  mAllBoundComps = new TBoundComp3D* [N_BoundComps];

  if ( mTetOut.trifacemarkerlist == NULL )
    mTetOut.trifacemarkerlist = new int [n_trifaces];

  for (int i=0;i<n_trifaces;++i)
  {
    if ( mTetOut.adjtetlist[2*i+1] == -1 ||
      mTetOut.adjtetlist[2*i  ] == -1   )
    {
      mAllBoundComps[counter] = new TBdPlane(counter);

      p[0] = mTetOut.pointlist[3*trifaces[3*i  ]  ];
      a[0] = mTetOut.pointlist[3*trifaces[3*i+1]  ] - p[0];
      b[0] = mTetOut.pointlist[3*trifaces[3*i+2]  ] - p[0];
      p[1] = mTetOut.pointlist[3*trifaces[3*i  ]+1];
      a[1] = mTetOut.pointlist[3*trifaces[3*i+1]+1] - p[1];
      b[1] = mTetOut.pointlist[3*trifaces[3*i+2]+1] - p[1];
      p[2] = mTetOut.pointlist[3*trifaces[3*i  ]+2];
      a[2] = mTetOut.pointlist[3*trifaces[3*i+1]+2] - p[2];
      b[2] = mTetOut.pointlist[3*trifaces[3*i+2]+2] - p[2];

      Normalize3D (a);
      Normalize3D (b);
      Cross3D (a, b, n);

      ((TBdPlane*) mAllBoundComps[counter])->SetParams(p[0], p[1], p[2],
        a[0], a[1], a[2],
        n[0], n[1], n[2]);

      ++counter;

      mTetOut.trifacemarkerlist[i] = counter;
    }
    else
      mTetOut.trifacemarkerlist[i] = 0;
  }
  //   assert( counter == N_BoundComps );

  return 0;
}


int TTetGenMeshLoader::DistributeJoints(TBaseCell **CellTree)
{
  int trifaces[4];
  const int *TmpFV, *TmpLen;
  int MaxLen, triface;
  TShapeDesc *ShapeDesc;

  for (int i=0;i<mTetOut.numberoftetrahedra;++i)
  {
    ShapeDesc = CellTree[i]->GetShapeDesc();
    ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);

    for (int j=0;j<4;++j)
    {
      triface = FindTriface(mTetOut.tetrahedronlist[4*i+TmpFV[j*MaxLen  ]],
        mTetOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+1]],
        mTetOut.tetrahedronlist[4*i+TmpFV[j*MaxLen+2]]);

      assert ( triface != -1 );

      CellTree[i]->SetJoint(j, mAllJoints[triface]);

      // correct params
      if ( mAllJoints[triface]->GetType() == BoundaryFace )
      {
        // 	cout << "Correct params" << endl;

        TBoundFace *BoundFace = (TBoundFace*) mAllJoints[triface];
        TBoundComp3D *BoundComp = BoundFace->GetBoundComp();
        double x, y, z, t, s;
        double param1[4], param2[4];

        for (int k=0;k<TmpLen[j];++k)
        {
          CellTree[i]->GetVertex(TmpFV[MaxLen*j+k])->GetCoords(x,y,z);
          BoundComp->GetTSofXYZ(x,y,z, t, s);

          param1[k] = t;
          param2[k] = s;
        }

        BoundFace->SetParameters(param1, param2);
      }                                           // end if if ( mAllJoints[trifaces[j]]->GetType() == BoundaryFace )
    }
  }

  // set map type
  for (int i=0;i<mTetOut.numberoftrifaces;++i)
  {
    mAllJoints[i]->SetMapType();
  }
  return(0);  
  
}


int TTetGenMeshLoader::BuildMooNMDMesh(TBaseCell **&CellTree, int &N_RootCells,
double &StartX, double &StartY, double &StartZ,
double &BoundX, double &BoundY, double &BoundZ)
{
  AllocVertices(StartX, StartY, StartZ, BoundX, BoundY, BoundZ);

  AllocRootCells(CellTree, N_RootCells);

  AllocJoints(CellTree);

  DistributeJoints(CellTree);

  return 0;
}


int TTetGenMeshLoader::Tetgen()
{
  std::string opts;

  MakeOptionString(opts);

#ifdef __TETGEN_14X__
  if (mTetBeh.parse_commandline((char*) opts.c_str()))
  {
    if ( insertpoints )
    {
      cout << "adding points" << endl;
      tetrahedralize(&mTetBeh, &mTetIn, &mTetOut, &mTetAddIn);
    }
    else    
      tetrahedralize(&mTetBeh, &mTetIn, &mTetOut);      
  }
  else
  {
    cerr << "parse_commandline failed !" << endl;
    return 0;
  }
#else
  tetrahedralize((char*) opts.c_str(), &mTetIn, &mTetOut);
#endif

  OutPut("TetGen - Mesh generated" << endl);

  HashTrifaces();
  CreateAdjacency();

  return 1;
}

void TTetGenMeshLoader::FindBoundaryPoints_msh(int *pointlist, double *normallist)
{ 
  int N_Tetrahedra = mTetIn.numberoftetrahedra;
  int N_Points = mTetIn.numberofpoints;
  int a, b, c, d;
  int r, s, t;
  int i, j, k, l;
  
  int *facemarker = new int [4*N_Tetrahedra]; // 1 inner ; 0 outer
  
  for (i=0;i<4*N_Tetrahedra;++i)
    facemarker[i] = 0;
  
  for (i=0;i<N_Points;++i)
    pointlist[i] = 0;
  
  for (i=0;i<3*N_Points;++i)
    normallist[i] = 0.0;
  
  int FaceVertex[][3] =
  {
    {
      0, 1, 2
    }
    ,
    {
      0, 3, 1
    }
    ,
    {
      2, 1, 3
    }
    ,
    {
      0, 2, 3
    }
  };
  
  for (i=0;i<N_Tetrahedra;++i)
  {
    for (j=0;j<4;++j)
    {
      a = mTetIn.tetrahedronlist[4*i+FaceVertex[j][0]];
      b = mTetIn.tetrahedronlist[4*i+FaceVertex[j][1]];
      c = mTetIn.tetrahedronlist[4*i+FaceVertex[j][2]];
      
      for (k=i+1;k<N_Tetrahedra;++k)
      {
	for (l=0;l<4;++l)
	{
	  r = mTetIn.tetrahedronlist[4*k+FaceVertex[l][0]];
	  s = mTetIn.tetrahedronlist[4*k+FaceVertex[l][1]];
	  t = mTetIn.tetrahedronlist[4*k+FaceVertex[l][2]];
	  
	  if ( ( a == r || a == s || a == t ) &&
	       ( b == r || b == s || b == t ) &&
	       ( c == r || c == s || c == t )    )
	  {
	    facemarker[4*i+j] = 1;
	    facemarker[4*k+l] = 1;
	  }
	}
      }
    }
  }
  
  double vec1[3], vec2[3], norm[3];
  
  for (i=0;i<N_Tetrahedra;++i)
  {
    double barycenter[3];
    
    a = mTetIn.tetrahedronlist[4*i  ];
    b = mTetIn.tetrahedronlist[4*i+1];
    c = mTetIn.tetrahedronlist[4*i+2];
    d = mTetIn.tetrahedronlist[4*i+3];
    
    
    
    barycenter[0] = 0.25*(mTetIn.pointlist[3*a  ] +
			  mTetIn.pointlist[3*b  ] +
			  mTetIn.pointlist[3*c  ] +
			  mTetIn.pointlist[3*d  ]  );
			  
    barycenter[1] = 0.25*(mTetIn.pointlist[3*a+1] +
			  mTetIn.pointlist[3*b+1] +
			  mTetIn.pointlist[3*c+1] +
			  mTetIn.pointlist[3*d+1]  );
			 
    barycenter[2] = 0.25*(mTetIn.pointlist[3*a+2] +
			  mTetIn.pointlist[3*b+2] +
			  mTetIn.pointlist[3*c+2] +
			  mTetIn.pointlist[3*d+2]  );
    
    for (j=0;j<4;++j)
    {
      if ( facemarker[4*i+j] == 0 )
      {
	double facebarycenter[3];
	double inside[3], scalar;
	
	a = mTetIn.tetrahedronlist[4*i+FaceVertex[j][0]];
	b = mTetIn.tetrahedronlist[4*i+FaceVertex[j][1]];
	c = mTetIn.tetrahedronlist[4*i+FaceVertex[j][2]];
	
	facebarycenter[0] = 1.0/3.0*(mTetIn.pointlist[3*a  ] +
				     mTetIn.pointlist[3*b  ] +
				     mTetIn.pointlist[3*c  ]  );
			  
	facebarycenter[1] = 1.0/3.0*(mTetIn.pointlist[3*a+1] +
				     mTetIn.pointlist[3*b+1] +
				     mTetIn.pointlist[3*c+1]  );
				     
	facebarycenter[2] = 1.0/3.0*(mTetIn.pointlist[3*a+2] +
				     mTetIn.pointlist[3*b+2] +
				     mTetIn.pointlist[3*c+2]  );
	
	pointlist[a] += 1;
	pointlist[b] += 1;
	pointlist[c] += 1;
	
	vec1[0] = mTetIn.pointlist[3*b  ] - mTetIn.pointlist[3*a  ];
	vec1[1] = mTetIn.pointlist[3*b+1] - mTetIn.pointlist[3*a+1];
	vec1[2] = mTetIn.pointlist[3*b+2] - mTetIn.pointlist[3*a+2];
	
	vec2[0] = mTetIn.pointlist[3*c  ] - mTetIn.pointlist[3*a  ];
	vec2[1] = mTetIn.pointlist[3*c+1] - mTetIn.pointlist[3*a+1];
	vec2[2] = mTetIn.pointlist[3*c+2] - mTetIn.pointlist[3*a+2];
	
	Cross3D(vec1, vec2, norm);
	
	inside[0] = barycenter[0] - facebarycenter[0];
	inside[1] = barycenter[1] - facebarycenter[1];
	inside[2] = barycenter[2] - facebarycenter[2];
	
	scalar = norm[0]*inside[0] + norm[1]*inside[1] + norm[2]*inside[2];
	if ( scalar < 0 )
	{
	  norm[0] = - norm[0];
	  norm[1] = - norm[1];
	  norm[2] = - norm[2];
	}
	
	normallist[3*a  ] += norm[0];
	normallist[3*a+1] += norm[1];
	normallist[3*a+2] += norm[2];
	
	normallist[3*b  ] += norm[0];
	normallist[3*b+1] += norm[1];
	normallist[3*b+2] += norm[2];
	
	normallist[3*c  ] += norm[0];
	normallist[3*c+1] += norm[1];
	normallist[3*c+2] += norm[2];
      }
    }
  }
  
  for (i=0;i<N_Points;++i)
  {
    if ( pointlist[i] > 0 )
    {
      normallist[3*i  ] /= pointlist[i];
      normallist[3*i+1] /= pointlist[i];
      normallist[3*i+2] /= pointlist[i];
      
//       double bla = mTetIn.pointlist[3*i  ]*normallist[3*i  ] + 
// 		   mTetIn.pointlist[3*i+1]*normallist[3*i+1] +
// 		   mTetIn.pointlist[3*i+2]*normallist[3*i+2];
// 		   
//       double r =   sqrt(mTetIn.pointlist[3*i  ]*mTetIn.pointlist[3*i  ] + 
// 		        mTetIn.pointlist[3*i+1]*mTetIn.pointlist[3*i+1] +
// 			mTetIn.pointlist[3*i+2]*mTetIn.pointlist[3*i+2]);
// 	      
//       cout << normallist[3*i  ] << " ";
//       cout << normallist[3*i+1] << " ";
//       cout << normallist[3*i+2] << " " << bla << " " << r << endl;
    }
  }
  
  delete [] facemarker;
}

void TTetGenMeshLoader::MakeBoundaryLayer_msh()
{
  int N_Points = mTetIn.numberofpoints;
  int *pointlist;
  double *normallist, delta = 0.05;
  int N_NewPoints=0;
  int count;
  
  pointlist = new int [N_Points];
  normallist = new double [3*N_Points];
  
  FindBoundaryPoints_msh(pointlist, normallist);
  
  for (int i=0;i<N_Points;++i)
  {
    if ( pointlist[i] > 0 )
    {
      ++N_NewPoints;
    }
  }
  
  mTetAddIn.pointlist = new double [3*N_NewPoints];
  mTetAddIn.pointmarkerlist = new int [N_NewPoints];
  mTetAddIn.numberofpoints = N_NewPoints;
  mTetAddIn.numberofpointattributes = 0;
  mTetAddIn.numberofpointmtrs = 0;
  
  count = 0;
  for (int i=0;i<N_Points;++i)
  {
    if ( pointlist[i] > 0 )
    {
      mTetAddIn.pointlist[3*count  ] = mTetIn.pointlist[3*i  ] + delta*normallist[3*i  ];
      mTetAddIn.pointlist[3*count+1] = mTetIn.pointlist[3*i+1] + delta*normallist[3*i+1];
      mTetAddIn.pointlist[3*count+2] = mTetIn.pointlist[3*i+2] + delta*normallist[3*i+2];
      
      mTetAddIn.pointmarkerlist[count] = 0;
      
      ++count;
    }
  }
  
  insertpoints = true;
  
  delete [] pointlist;
  delete [] normallist;
  
}

void TTetGenMeshLoader::MakeBoundaryLayer_smesh()
{
  int N_Points = mTetIn.numberoffacets;
  
  mTetAddIn.pointlist = new REAL [3*N_Points];
  mTetAddIn.pointmarkerlist = new int [3*N_Points];
  mTetAddIn.numberofpoints = N_Points;
  mTetAddIn.numberofpointattributes = 0;
  mTetAddIn.numberofpointmtrs = 0;
  
//   cout << mTetIn.numberofpoints << endl;
//   for ( int i=0;i<mTetIn.numberofpoints;++i)
//   {
//     double norm;
//     
//     norm = sqrt(mTetIn.pointlist[3*i]*mTetIn.pointlist[3*i] +
// 		mTetIn.pointlist[3*i+1]*mTetIn.pointlist[3*i+1] +
// 		mTetIn.pointlist[3*i+2]*mTetIn.pointlist[3*i+2]);
// 		
//     cout << norm << endl;
//   }
  
  for (int i=0;i<N_Points;++i)
  {
    int a, b, c;
    REAL p1[3], p2[3], n[3];
    REAL barycenter;
    tetgenio::polygon poly = mTetIn.facetlist[i].polygonlist[0];
    
//     cout << poly.numberofvertices << endl;
    
    a = poly.vertexlist[0] - mTetIn.firstnumber;
    b = poly.vertexlist[1] - mTetIn.firstnumber;
    c = poly.vertexlist[2] - mTetIn.firstnumber;
    
    p1[0] = mTetIn.pointlist[3*b  ] - mTetIn.pointlist[3*a  ];
    p1[1] = mTetIn.pointlist[3*b+1] - mTetIn.pointlist[3*a+1];
    p1[2] = mTetIn.pointlist[3*b+2] - mTetIn.pointlist[3*a+2];
    
    p2[0] = mTetIn.pointlist[3*c  ] - mTetIn.pointlist[3*a  ];
    p2[1] = mTetIn.pointlist[3*c+1] - mTetIn.pointlist[3*a+1];
    p2[2] = mTetIn.pointlist[3*c+2] - mTetIn.pointlist[3*a+2];
    
    Cross3D(p1, p2, n);
    
//     cout << n[0] << " , " << n[1] << " , " << n[2] << endl;
//     cout << mTetIn.pointlist[3*a  ] << " , " << mTetIn.pointlist[3*a+1] << " , " << mTetIn.pointlist[3*a+2] << " | ";
//     cout << n[0]*mTetIn.pointlist[3*a  ] + n[1]*mTetIn.pointlist[3*a+1] + n[2]*mTetIn.pointlist[3*a+2]<< endl;
    for (int j=0;j<3;++j)
    {
      barycenter = 1.0/3.0 * (mTetIn.pointlist[3*a+j] + mTetIn.pointlist[3*b+j] + mTetIn.pointlist[3*c+j]);
      mTetAddIn.pointlist[3*i+j] = barycenter - 0.05 * n[j];
    }
  }
  
  insertpoints = true;
}

int TTetGenMeshLoader::Generate(int N_Points, double *Points, int N_Facets, int *Facets,
				int N_Regions, double *Regions, TDummyDomain *Domain)
{
//   tetgenio::facet *Facet;
//   
//   mTetIn.numberofpoints = N_Points;
//   mTetIn.numberofpointattributes = 0;
//   mTetIn.numberofpointmtrs = 0;
//   mTetIn.pointlist = Points;
//   
//   mTetIn.facetlist = new tetgenio::facet [N_Facets];
//   mTetIn.numberoffacets = N_Facets;
//   for (int i=0;i<N_Facets;++i)
//   {
//     Facet = mTetIn.facetlist+i;
//     
//     mTetIn.init(Facet);        
//     
//     Facet->numberofpolygons = 1;
//     Facet->polygonlist = new tetgenio::polygon [1];
//     
//     mTetIn.init(Facet->polygonlist);
//     
//     Facet->polygonlist->numberofvertices = 3;
//     
//     Facet->polygonlist->vertexlist = new int [3];
//     
//     Facet->polygonlist->vertexlist[0] = Facets[3*i+0];
//     Facet->polygonlist->vertexlist[1] = Facets[3*i+1];
//     Facet->polygonlist->vertexlist[2] = Facets[3*i+2];
//   }
  
//   mTetIn.numberofregions = N_Regions;
//   mTetIn.regionlist = Regions;
  
  plc = true;
  
  mTetIn.load_poly(TDatabase::ParamDB->SMESHFILE);
  
  Tetgen();
  
  BuildBoundary2(Domain->BdParts, Domain->N_BoundParts,
		 Domain->N_BoundComps, Domain->StartBdCompID, Domain->Interfaces);

  BuildMooNMDMesh (Domain->CellTree, Domain->N_RootCells,
        Domain->StartX, Domain->StartY, Domain->StartZ,
        Domain->BoundX, Domain->BoundY, Domain->BoundZ);
  
  mTetIn.pointlist = NULL;
  mTetIn.regionlist = NULL;
	
  return 0;
}
