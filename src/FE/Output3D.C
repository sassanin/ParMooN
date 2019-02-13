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
// Class:       TOutput3D
// Purpose:     store given data and realize output
//
// Author:      Gunar Matthies (24.07.2000)
//
// History:     start of implementation 24.07.2000 (Gunar Matthies)
//              added parvtk, tecplot, mesh data storage 12.03.12 (Sashikumaar)
// =======================================================================

#include <FEDatabase3D.h>
#include <Output3D.h>
#include <BaseCell.h>
#include <Joint.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <TECIO.h>

#include <TetraAffin.h>
#include <HexaAffin.h>
#include <HexaTrilinear.h>

#ifdef __COMPAQ__
  // include byte reordering routines for Compaq only
  #include <Utilities.h>
#endif

#include <sstream>
// #include <malloc.h>
#include <dirent.h> 
#include <unistd.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

/** constructor maximum number of these things */
TOutput3D::TOutput3D(int maxn_fespaces, int maxn_scalar,
                     int maxn_vect, int maxn_parameters, TDomain *domain, TCollection *coll, const char *name)
{
  MaxN_FESpaces=maxn_fespaces;
  N_FESpaces=0;

  MaxN_ScalarVar=maxn_scalar;
  N_ScalarVar=0;

  MaxN_VectorVar=maxn_vect;
  N_VectorVar=0;

  MaxN_Parameters=maxn_parameters;
  N_Parameters=0;

  FESpaceArray=new TFESpace3D*[MaxN_FESpaces];

  FEFunctionArray=new TFEFunction3D*[MaxN_ScalarVar];

  FEVectFunctArray=new TFEVectFunct3D*[MaxN_VectorVar];

  ParameterValues=new double[MaxN_Parameters];

  ParameterDescription=new const char*[MaxN_Parameters];

  Coll = coll;

  Domain = domain;

  Data = NULL;

  if (name) Name = strdup(name);
  else Name = strdup("none");
 
}

/** add a FESpace into this output object (internal use) */
int TOutput3D::AddFESpace(TFESpace3D *fespace)
{
  TFESpace3D **NewStorage;
  int i;

  if(N_FESpaces==0)
  { 
    // no fespace in this output object    
    Coll=fespace->GetCollection();
  }

  // check whether fespace is already stored
  for(i=0;i<N_FESpaces;i++)
    if(FESpaceArray[i]==fespace) return i+1;

  // this space is based on a different Collection
  if(fespace->GetCollection()!=Coll) return 0;

  if(MaxN_FESpaces<=N_FESpaces)
  {
    // enlarge storage
    NewStorage=new TFESpace3D*[MaxN_FESpaces+5];
    memcpy(NewStorage, FESpaceArray, sizeof(TFESpace3D *)*MaxN_FESpaces);
    MaxN_FESpaces +=5;
    delete FESpaceArray;
    FESpaceArray=NewStorage;
  }

  // store space on the next place
  FESpaceArray[N_FESpaces]=fespace;
  N_FESpaces++;

  return N_FESpaces;
}

/** add a FEFunction into this output object */
int TOutput3D::AddFEFunction(TFEFunction3D *fefunction)
{
  TFEFunction3D **NewStorage;

  if(!(AddFESpace(fefunction->GetFESpace3D()))) return 0;

  if(MaxN_ScalarVar<=N_ScalarVar)
  {
    // enlarge storage
    NewStorage=new TFEFunction3D*[MaxN_ScalarVar+5];
    memcpy(NewStorage, FEFunctionArray, MaxN_ScalarVar*sizeof(TFEFunction3D *));
    MaxN_ScalarVar +=5;
    delete FEFunctionArray;
    FEFunctionArray=NewStorage;
  }

  // store function on the next place
  FEFunctionArray[N_ScalarVar]=fefunction;
  N_ScalarVar++;

  return N_ScalarVar;
}

/** add a FEVectFunct into this output object */
int TOutput3D::AddFEVectFunct(TFEVectFunct3D *fevectfunct)
{
  TFEVectFunct3D **NewStorage;

  if(!(AddFESpace(fevectfunct->GetFESpace3D()))) return 0;

  if(MaxN_VectorVar<=N_VectorVar)
  {
    // enlarge storage
    NewStorage=new TFEVectFunct3D*[MaxN_VectorVar+5];
    memcpy(NewStorage, FEVectFunctArray, 
                MaxN_VectorVar*sizeof(TFEVectFunct3D *));
    MaxN_VectorVar +=5;
    delete FEVectFunctArray;
    FEVectFunctArray=NewStorage;
  }

  // store function on the next place
  FEVectFunctArray[N_VectorVar]=fevectfunct;
  N_VectorVar++;

  return N_VectorVar;
}

/** add parameter into this output object */
int TOutput3D::AddParameter(double value, const char *descr)
{
  double *NewValues;
  const char **NewDescr;

  if(MaxN_Parameters<=N_Parameters)
  {
    // enlarge storage
    NewValues=new double[MaxN_Parameters+5];
    NewDescr=new const char *[MaxN_Parameters+5];
    memcpy(NewValues, ParameterValues, sizeof(double)*MaxN_Parameters);
    memcpy(NewDescr, ParameterDescription, sizeof(const char *)*MaxN_Parameters);
    MaxN_Parameters +=5;
    delete ParameterValues;
    delete ParameterDescription;
    ParameterValues=NewValues;
    ParameterDescription=NewDescr;
  }

  // value and description on next place
  ParameterValues[N_Parameters]=value;
  ParameterDescription[N_Parameters]=strdup(descr);
  N_Parameters++;

  return N_Parameters;
}

static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr=new int[2*len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete rr;

}

static int GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }
  return m;
}

/** write stored data into a .vtk-file*/
int TOutput3D::Write(std::string basename, int i, double t)
{
  std::ostringstream os;
  os.seekp(std::ios::beg);
  os << basename << i << ".vtk"<<ends;
  
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    cout << " Output3D:: writing " << os.str() << endl;
    this->WriteVtk(os.str().c_str());
  }
  
  return 0;
  // put more here if needed
}


/** write stored data into a grape file */
int TOutput3D::WriteGrape(const char *name)
{
  int i,j,k,l,m,p,N_;
  int headerlength;
  int *header, first[3];
  int format, version, N_Elements, N_Vertices;
  int N_LocVertices, *N_DOF, *N_Comp, *FESpaceNumber;
  int N_LocFaces, *Neighbours;
  TBaseCell *cell, *neigh;
  TJoint *joint;
  TVertex **Vertices, *Last, *Current;
  int *NumberVertex;
  double *ParamPtr;
  TFESpace3D *fespace;
  TFEFunction3D *fefunction;
  TFEVectFunct3D *fevectfunct;
  double *Coords;
  int *VertexNumbers;
  int *BaseFuncts;
  FE3D FE_ID;
  BaseFunct3D BaseFunct;
  const char TimeString[] = "Time";

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  if(N_Parameters == 0)
    AddParameter(TDatabase::TimeDB->CURRENTTIME, TimeString);

  format=1;
  version=1; 
  headerlength=50+20*(N_FESpaces+N_ScalarVar+N_VectorVar+N_Parameters);
  // cout << "Header length: " << headerlength << endl;
  header=new int[headerlength];
  first[0]=format;
  first[1]=version;
  first[2]=headerlength;
  ParameterValues[0] = TDatabase::TimeDB->CURRENTTIME;

  N_Elements = Coll->GetN_Cells();
  // cout << "number of elements: " << N_Elements << endl;

  for(i=0;i<N_Elements;i++)
  {
    cell=Coll->GetCell(i);
    k=cell->GetN_Joints();
    for(j=0;j<k;j++)
    {  
      neigh=cell->GetJoint(j)->GetNeighbour(cell); 
      if(neigh) neigh->SetClipBoard(-1);
    }       
    cell->SetClipBoard(-1); 
  } // endfor i             

  N_LocVertices = 0;
  N_LocFaces = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
    N_LocFaces += cell->GetN_Faces();
    cell->SetClipBoard(i);
  }

  Vertices = new TVertex*[N_LocVertices];
  Neighbours = new int[N_LocFaces];
  N_ = 0;
  m = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }

    k=cell->GetN_Faces();
    for(j=0;j<k;j++)
    {
      joint = cell->GetJoint(j);
      neigh = joint->GetNeighbour(cell);

      if(neigh)
      {
        if((l=neigh->GetClipBoard()) == -1)
          Neighbours[m] = -3;
        else
          Neighbours[m] = l;
      }
      else
      {
        if(joint->GetType() == BoundaryFace)
          Neighbours[m] = -1;
        else
          Neighbours[m] = -2;
      }
      m++;
    }
  }

  // sort the Vertices array
  Sort(Vertices, N_);

  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }

  // fill in geometry information
  header[0]=3; // dimension
  header[1]=N_Elements;
  header[2]=N_Vertices;
  header[3]=N_LocVertices;
  header[4]=N_LocFaces;
  header[5]=N_FESpaces;
  header[6]=N_ScalarVar;
  header[7]=N_VectorVar;
  header[8]=N_Parameters;
  // header[9] is unused
#ifdef __COMPAQ__
  SwapIntArray(header, 10);
#endif
  strncpy((char *)(header+10),"Output generated by 'MooN_MD' program",160);

  // cout << "header[0]= " << header[0] << endl;
  // cout << "header[1]= " << header[1] << endl;
  // cout << "header[2]= " << header[2] << endl;
  // cout << "header[3]= " << header[3] << endl;
  // cout << "header[4]= " << header[4] << endl;
  // cout << "header[5]= " << header[5] << endl;
  // cout << "header[6]= " << header[6] << endl;
  // cout << "header[7]= " << header[7] << endl;

  // putting in parameter information
  N_=50;
  for(i=0;i<N_Parameters;i++)
  {    
    ParamPtr=(double *)(header+N_);
    *ParamPtr = ParameterValues[i];
#ifdef __COMPAQ__
    SwapDoubleArray(ParamPtr, 1);
#endif
    OutPut("time " << ParameterValues[i] << " ");
    memcpy(header+N_,ParamPtr,SizeOfDouble);
    N_ += 2;
    strncpy((char *)(header+N_),ParameterDescription[i],72);
    N_ += 18;
  }

  // putting in fespace information
  N_DOF=new int[N_FESpaces];
  for(i=0;i<N_FESpaces;i++)
  {
    fespace=FESpaceArray[i];
    N_DOF[i]=fespace->GetN_DegreesOfFreedom();
    header[N_]=N_DOF[i];
    // cout << "N_DOF " << header[N_] << endl;
    N_++;
    j=(fespace->GetBeginIndex())[N_Elements];
    header[N_]=j;
    // cout << "N_LocDOF " << header[N_] << endl;
    N_++;
#ifdef __COMPAQ__
    SwapIntArray(header+N_-2, 2);
#endif
    strncpy((char *)(header+N_),fespace->GetName(),72);
    // cout << (char *)(header+N_) << endl;
    N_ += 18;
  }

  // putting in fefunction information
  N_Comp=new int[N_ScalarVar+N_VectorVar];
  FESpaceNumber=new int[N_ScalarVar+N_VectorVar];
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comp[i]=1;
    header[N_]=N_Comp[i];
    N_++;

    // find fespace for this fefunction
    fespace=FEFunctionArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
    header[N_]=j;
    N_++;

#ifdef __COMPAQ__
    SwapIntArray(header+N_-2, 2);
#endif
    strncpy((char *)(header+N_),FEFunctionArray[i]->GetName(),72);
    // cout << (char *)(header+N_) << endl;
    N_ += 18;
  }

  // putting in fevectfunct information
  k=N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comp[k]=FEVectFunctArray[i]->GetN_Components();
    header[N_]=N_Comp[k];
    N_++;

    // find fespace for this fefunction
    fespace=FEVectFunctArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
    header[N_]=j;
    N_++;

#ifdef __COMPAQ__
    SwapIntArray(header+N_-2, 2);
#endif
    strncpy((char *)(header+N_),FEVectFunctArray[i]->GetName(),72);
    // cout << (char *)(header+N_) << endl;
    N_ += 18;
  }

  // cout << "N_= " << N_ << endl;

#ifdef __COMPAQ__
  SwapIntArray(first, 3);
#endif
  dat.write((char *)first, sizeof(int)*3);
  dat.write((char *)header, sizeof(int)*headerlength);
  delete header;

  Coords=new double[3*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1],Coords[N_+2]);
      k++;
      N_ += 3;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

#ifdef __COMPAQ__
  SwapDoubleArray(Coords, 3*N_Vertices);
#endif
  dat.write((char *)Coords,sizeof(double)*3*N_Vertices);
  delete Coords;

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      // cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

#ifdef __COMPAQ__
  SwapIntArray(VertexNumbers, N_LocVertices);
  SwapIntArray(Neighbours, N_LocFaces);
#endif
  dat.write((char *)VertexNumbers,sizeof(int)*N_LocVertices);
  dat.write((char *)Neighbours,sizeof(int)*N_LocFaces);
  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;
  delete Neighbours;

  BaseFuncts=new int[N_Elements];
  for(p=0;p<N_FESpaces;p++)
  {
    fespace=FESpaceArray[p];

    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);

      FE_ID = fespace->GetFE3D(i, cell);
      BaseFunct = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D_ID();
      BaseFuncts[i] = BaseFunct;
    } // endfor i

#ifdef __COMPAQ__
    SwapIntArray(BaseFuncts, N_Elements);
#endif
    dat.write((char *)BaseFuncts,sizeof(int)*N_Elements);
    m=(fespace->GetBeginIndex())[N_Elements];
#ifdef __COMPAQ__
    SwapIntArray(fespace->GetGlobalNumbers(), m);
#endif
    dat.write((char *)fespace->GetGlobalNumbers(),sizeof(int)*m);
#ifdef __COMPAQ__
    // putting bytes in right order back
    SwapIntArray(fespace->GetGlobalNumbers(), m);
#endif

  } // endfor p
  delete BaseFuncts;

  for(p=0;p<N_ScalarVar;p++)
  {
    fefunction=FEFunctionArray[p];
    m=fefunction->GetLength();
    // cout << m << "  " << N_DOF[FESpaceNumber[p]] << endl;
#ifdef __COMPAQ__
    SwapDoubleArray(fefunction->GetValues(), m);
#endif
    dat.write((char *)fefunction->GetValues(),sizeof(double)*m);
#ifdef __COMPAQ__
    // putting bytes in right order back
    SwapDoubleArray(fefunction->GetValues(), m);
#endif
  }

  for(p=0;p<N_VectorVar;p++)
  {
    fevectfunct=FEVectFunctArray[p];
    m=fevectfunct->GetLength() * N_Comp[p+N_ScalarVar];
#ifdef __COMPAQ__
    SwapDoubleArray(fevectfunct->GetValues(), m);
#endif
    dat.write((char *)fevectfunct->GetValues(),sizeof(double)*m);
#ifdef __COMPAQ__
    // putting bytes in right order back
    SwapDoubleArray(fevectfunct->GetValues(), m);
#endif
  }

  dat.close();

  OutPut("wrote output into file: " << name << endl);

  delete N_DOF;
  delete N_Comp;
  delete FESpaceNumber;

  return 0;
}

/** write stored data into a tecplot file */
int TOutput3D::WriteTecplot(const char *name)
{
  int i,j,k,l,m,n, N_;
  int N_Elements, N_Vertices, N_LocVertices;
  int N_LocDOF;
  int *VertexNumbers, *NumberVertex;
  int *IntArray;
  double *Coords, *DoubleArray;
  Shapes ShapeType, ShapeType2;
  TBaseCell *cell;
  TVertex **Vertices, *Current, *Last;
  FE3D FE_ID;
  TBaseFunct3D *bf;
  TFESpace3D *fespace;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double BFValues[MaxN_BaseFunctions3D], *Coeffs;
  double xi, eta, zeta, value;
  int MaxN_Comp, N_Comp, Length;

  static double HexaCoords[] =
    { -1, -1, -1,  1, -1, -1,  1, 1, -1,  -1, 1, -1 ,
      -1, -1,  1,  1, -1,  1,  1, 1,  1,  -1, 1,  1 };

  static double TetraCoords[] =
    { 0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1 };

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  N_Elements = Coll->GetN_Cells();
  ShapeType = Coll->GetCell(0)->GetType();
  if(ShapeType == Brick) ShapeType = Hexahedron;
  N_LocVertices = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    ShapeType2 = cell->GetType();
    if(ShapeType2 == Brick) ShapeType2 = Hexahedron;
    if(ShapeType2 != ShapeType)
    {
      Error("All cells have to have the same shape!" << endl);
      dat.close();
      return -1;
    }

    N_LocVertices += cell->GetN_Vertices();
  } // endfor i

  Vertices = new TVertex*[N_LocVertices];
  N_ = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_] = cell->GetVertex(j);
      N_++;
    } // endfor k
  } // endfor i

  // sort the Vertices array
  Sort(Vertices, N_);

  Last = NULL;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  Coords = new double[3*N_Vertices];
  VertexNumbers = new int[N_LocVertices];
  NumberVertex = new int[N_LocVertices];

  dat << "Variables = \"X\", \"Y\", \"Z\"";

  for(i=0;i<N_ScalarVar;i++)
  {
    dat << ", \"" << FEFunctionArray[i]->GetName() << "\"";
  }
  MaxN_Comp = 1;

  for(i=0;i<N_VectorVar;i++)
  {
    k = FEVectFunctArray[i]->GetN_Components();
    for(j=0;j<k;j++)
    {
      dat << ", \"" << FEVectFunctArray[i]->GetName();
      dat << "_" << j+1 << "\"";
    }
    if(k>MaxN_Comp) MaxN_Comp = k;
  }
  dat << endl;
 
  
  DoubleArray = new double[MaxN_Comp*N_Vertices];
  cout << "max: " << MaxN_Comp*N_Vertices << endl;
  IntArray = new int[N_Vertices];

  dat << "ZONE N= " << N_Vertices << ", E= " << N_Elements;
  dat << ", F=FEBLOCK, ";

  switch(ShapeType)
  {
    case Tetrahedron:
      dat << "ET=TETRAHEDRON" << endl;
      break;
    case Hexahedron:
    case Brick:
      dat << "ET=BRICK" << endl;
      break;
    default:
      
    break;  
  }

  cout << N_Vertices << endl;
  cout << N_LocVertices << endl;

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1],Coords[N_+2]);
      k++;
      N_ += 3;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  // dat.write((char *)Coords,sizeof(double)*3*N_Vertices);
  for(i=0;i<N_Vertices;i++)
    dat << setw(12) << Coords[3*i] << endl;
  dat << endl;
  for(i=0;i<N_Vertices;i++)
    dat << setw(12) << Coords[3*i+1] << endl;
  dat << endl;
  for(i=0;i<N_Vertices;i++)
    dat << setw(12) << Coords[3*i+2] << endl;
  dat << endl;
  delete Coords;
  exit(0);
  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      // cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

  // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace3D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
    
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetType())
        {
          case Tetrahedron: 
            xi = TetraCoords[3*j];
            eta = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
          break;

          case Hexahedron: 
          case Brick: 
            xi = HexaCoords[3*j];
            eta = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
          break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        value = 0;
        for(l=0;l<N_LocDOF;l++)
          value += BFValues[l] * Coeffs[DOF[l]];
        DoubleArray[VertexNumbers[m]] += value;
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    for(i=0;i<N_Vertices;i++)
    {
      DoubleArray[i] /= IntArray[i];
      dat << DoubleArray[i] << endl;
    }
    dat << endl;

  } // endfor k

  // write vector-valued variables into file
  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace3D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    cout << "n: " << N_Vertices*N_Comp << endl;
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
    
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetType())
        {
          case Tetrahedron: 
            xi = TetraCoords[3*j];
            eta = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
          break;

          case Hexahedron: 
          case Brick: 
            xi = HexaCoords[3*j];
            eta = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
          break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        for(n=0;n<N_Comp;n++)
        {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l] + n*Length];
          DoubleArray[VertexNumbers[m] + n*N_Vertices] += value;
        }
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    l = 0;
    for(j=0;j<N_Comp;j++)
    {
      for(i=0;i<N_Vertices;l++,i++)
      {
        // DoubleArray[l] /= IntArray[i];
        dat << DoubleArray[l] << endl;
      }
      dat << endl;
    } // endfor l
    dat << endl;
  } // endfor k

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      dat << setw(5) << VertexNumbers[m]+1;
      m++;
    } // endfor j
    dat << endl;
  } //endfor i

  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;
  delete DoubleArray;
  delete IntArray;

  dat.close();
  
  return 0;
}

/** write stored data into a GMV file */
int TOutput3D::WriteGMV(const char *name)
{
  int i,j,k,l,m,n, N_;
  int N_Elements, N_Vertices, N_LocVertices;
  int N_LocDOF;
  int *VertexNumbers, *NumberVertex;
  int *IntArray;
  double *Coords, *DoubleArray;
  Shapes ShapeType, ShapeType2;
  TBaseCell *cell;
  TVertex **Vertices, *Current, *Last;
  FE3D FE_ID;
  TBaseFunct3D *bf;
  TFESpace3D *fespace;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double BFValues[MaxN_BaseFunctions3D], *Coeffs;
  double xi, eta, zeta, value;
  int MaxN_Comp, N_Comp, Length;
  int Number;
  int MaxSubGridID;

  char  matx[8];

  char gmvinput[] = "gmvinput";
  char ieeei4r8[] = "ieeei4r8";
  char nodes[]    = "nodes   ";
  char cells[]    = "cells   ";
  char hex[]      = "hex     ";
  char tetra[]    = "tetra   ";
  char velocity[] = "velocity";
  char variable[] = "variable";
  char endvars[]  = "endvars ";
  char endgmv[]   = "endgmv  ";
  char material[] = "material";

  static double HexaCoords[] =
    { -1, -1, -1,  1, -1, -1,  1, 1, -1,  -1, 1, -1 ,
      -1, -1,  1,  1, -1,  1,  1, 1,  1,  -1, 1,  1 };

  static double TetraCoords[] =
    { 0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1 };

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  MaxSubGridID = -1;
  N_Elements = Coll->GetN_Cells();
  N_LocVertices = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    if( (Number = cell->GetSubGridID()) > MaxSubGridID )
      MaxSubGridID = Number;
    N_LocVertices += cell->GetN_Vertices();
  } // endfor i

  Vertices = new TVertex*[N_LocVertices];
  N_ = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_] = cell->GetVertex(j);
      N_++;
    } // endfor k
  } // endfor i

  // sort the Vertices array
  Sort(Vertices, N_);

  Last = NULL;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  Coords = new double[3*N_Vertices];
  VertexNumbers = new int[N_LocVertices];
  NumberVertex = new int[N_LocVertices];
  
  dat.write(gmvinput, 8);
  dat.write(ieeei4r8, 8);
  dat.write(nodes, 8);

  Number = N_Vertices;
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char*)(&Number), SizeOfInt);

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+N_Vertices],
                             Coords[N_+2*N_Vertices]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

#ifdef __COMPAQ__
  SwapDoubleArray(Coords, 3*N_Vertices);
#endif
  dat.write((char *)Coords,sizeof(double)*3*N_Vertices);

  dat.write(cells, 8);
  Number = N_Elements;
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char*)(&Number), SizeOfInt);

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l]+1;
      // cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
    switch(k)
    {
      case 4:
        dat.write(tetra, 8);
        Number = 4;
#ifdef __COMPAQ__
         SwapIntArray(&Number, 1);
#endif
        dat.write((char *)(&Number), SizeOfInt);
      break;

      case 8:
        dat.write(hex, 8);
        Number = 8;
#ifdef __COMPAQ__
         SwapIntArray(&Number, 1);
#endif
        dat.write((char *)(&Number), SizeOfInt);
      break;

      default:
        Error("Unknown cell type! " << endl);
        exit(-1);
    }
#ifdef __COMPAQ__
    SwapIntArray(&VertexNumbers[m-k], k);
#endif
    dat.write((char *)(&VertexNumbers[m-k]), k*SizeOfInt);
#ifdef __COMPAQ__
    // swap back
    SwapIntArray(&VertexNumbers[m-k], k);
#endif
  } //endfor i

  MaxSubGridID++;
  dat.write(material, 8);
  Number = MaxSubGridID;
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char*)(&Number), SizeOfInt);
  Number = 0; // cell data
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char*)(&Number), SizeOfInt);
  for(i=0;i<MaxSubGridID;i++)
  {
    sprintf(matx,"mat%i",i);
    dat.write(matx, 8);
  }
  
  for(i=0;i<N_Elements;i++)
  {
    Number = Coll->GetCell(i)->GetSubGridID()+1;
#ifdef __COMPAQ__
    SwapIntArray(&Number, 1);
#endif
    dat.write((char*)(&Number), SizeOfInt);
  }

  DoubleArray = new double[3*N_Vertices];
  IntArray = new int[N_Vertices];

  if(N_VectorVar)
  {
    dat.write(velocity, 8);
    Number = 1;
#ifdef __COMPAQ__
    SwapIntArray(&Number, 1);
#endif
    dat.write((char*)(&Number), SizeOfInt);
    // write vector-valued variables into file
    for(k=0;k<N_VectorVar;k++)
    {
      fespace = FEVectFunctArray[k]->GetFESpace3D();
      N_Comp = FEVectFunctArray[k]->GetN_Components();
      Length = FEVectFunctArray[k]->GetLength();
      Coeffs = FEVectFunctArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
  
      memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
      memset(IntArray, 0, SizeOfInt*N_Vertices);
      m = 0;
      
      for(i=0;i<N_Elements;i++)
      {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();
  
        // find FE data for this element
        FE_ID = fespace->GetFE3D(i, cell);
        bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
        {
          switch(cell->GetType())
          {
            case Tetrahedron: 
              xi = TetraCoords[3*j];
              eta = TetraCoords[3*j+1];
              zeta = TetraCoords[3*j+2];
            break;
  
            case Hexahedron: 
            case Brick: 
              xi = HexaCoords[3*j];
              eta = HexaCoords[3*j+1];
              zeta = HexaCoords[3*j+2];
            break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
          }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);
          for(n=0;n<N_Comp;n++)
          {
            value = 0;
            for(l=0;l<N_LocDOF;l++)
              value += BFValues[l] * Coeffs[DOF[l] + n*Length];
            DoubleArray[VertexNumbers[m]-1 + n*N_Vertices] += value;
          }
          IntArray[VertexNumbers[m]-1]++;
          m++;
        } // endfor j
      } // endfor i
  
      l = 0;
      for(j=0;j<N_Comp;j++)
      {
        for(i=0;i<N_Vertices;l++,i++)
          DoubleArray[l] /= IntArray[i];
      } // endfor l
#ifdef __COMPAQ__
      SwapDoubleArray(DoubleArray, N_Vertices*N_Comp);
#endif
      dat.write((char*)DoubleArray, SizeOfDouble*N_Vertices*N_Comp);
#ifdef __COMPAQ__
      // swap back
      SwapDoubleArray(DoubleArray, N_Vertices*N_Comp);
#endif
    } // endfor k
  }

  dat.write(variable, 8);
  // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace3D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
      
    strncpy(matx, FEFunctionArray[k]->GetName(), 8);
    dat.write(matx, 8);
    Number = 1;
#ifdef __COMPAQ__
    SwapIntArray(&Number, 1);
#endif
    dat.write((char*)(&Number), SizeOfInt);
    
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetType())
        {
          case Tetrahedron: 
            xi = TetraCoords[3*j];
            eta = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
          break;

          case Hexahedron: 
          case Brick: 
            xi = HexaCoords[3*j];
            eta = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
          break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        value = 0;
        for(l=0;l<N_LocDOF;l++)
          value += BFValues[l] * Coeffs[DOF[l]];
        DoubleArray[VertexNumbers[m]-1] += value;
        IntArray[VertexNumbers[m]-1]++;
        m++;
      } // endfor j
    } // endfor i

    for(i=0;i<N_Vertices;i++)
      DoubleArray[i] /= IntArray[i];

#ifdef __COMPAQ__
    SwapDoubleArray(DoubleArray, N_Vertices); 
#endif
    dat.write((char*)DoubleArray, N_Vertices*SizeOfDouble);
#ifdef __COMPAQ__
    // swap back
    SwapDoubleArray(DoubleArray, N_Vertices); 
#endif

  } // endfor k
  dat.write(endvars, 8);

  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;
  delete DoubleArray;
  delete IntArray;
  delete Coords;

  dat.write(endgmv, 8);

  dat.close();
  
  OutPut("wrote output into file: " << name << endl);

  return 0;
}

/** write stored data into an amira file */
int TOutput3D::WriteAmira(const char *name)
{
  int i,j,k,l,m,n, N_;
  int N_Elements, N_Vertices, N_LocVertices;
  int N_LocDOF;
  int *VertexNumbers, *NumberVertex;
  int *IntArray;
  double *Coords, *DoubleArray;
  Shapes ShapeType, ShapeType2;
  TBaseCell *cell;
  TVertex **Vertices, *Current, *Last;
  FE3D FE_ID;
  TBaseFunct3D *bf;
  TFESpace3D *fespace;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double BFValues[MaxN_BaseFunctions3D], *Coeffs;
  double xi, eta, zeta, value;
  int MaxN_Comp, N_Comp, Length;
  int Number;
  int MaxSubGridID;
  int VariableCounter, VariableStart = 3;
  char *VarName;

  char AmiraMesh[] = "# AmiraMesh 3D BINARY 1.0";
  char DefineNodes[] = "define Nodes ";
  char DefineHexa[] = "define Hexahedra ";
  char Nodes[] = "Nodes { double[3] Coordinates } @1";
  char Hexahedra[] = "Hexahedra { int[8] Nodes } @2";
  char ScalarData[] = "Nodes { double Data } @";
  char VectorData[] = "Nodes { double[3] Data } @";
  char ScalarField[] = "Field { double ";
  char VectorField[] = "Field { double[3] ";
  char Linear[] = " } Linear(@";
  char Closed[] = ")";
  char At1[] = "@1";
  char At2[] = "@2";

  // format (x_i, y_i, z_i)
  static double HexaCoords[] =
    { -1, -1, -1,  1, -1, -1,  1, 1, -1,  -1, 1, -1 ,
      -1, -1,  1,  1, -1,  1,  1, 1,  1,  -1, 1,  1 };

  // format (x_i, y_i, z_i)
  static double TetraCoords[] =
    { 0, 0, 0,   1, 0, 0,   0, 1, 0,   0, 0, 1 };

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  for(k=0;k<N_VectorVar;k++)
  {
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    if(N_Comp != 3)
    {
      Error("WriteAmira can handle right now only three-component vector fields." << endl);
      dat.close();
      return -1;
    }
  }

  MaxSubGridID = -1;
  N_Elements = Coll->GetN_Cells();
  N_LocVertices = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    if( (Number = cell->GetSubGridID()) > MaxSubGridID )
      MaxSubGridID = Number;
    k = cell->GetN_Vertices();
    if(k!=8)
    {
      Error("WriteAmira can handle right now only hexadral meshes!" << endl);
      dat.close();
      return -1;
    }
    N_LocVertices += k;
  } // endfor i

  Vertices = new TVertex*[N_LocVertices];
  N_ = 0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k = cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_] = cell->GetVertex(j);
      N_++;
    } // endfor k
  } // endfor i

  // sort the Vertices array
  Sort(Vertices, N_);

  Last = NULL;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  Coords = new double[3*N_Vertices];
  VertexNumbers = new int[N_LocVertices];
  NumberVertex = new int[N_LocVertices];

  dat << AmiraMesh << endl;
  dat << endl;
  
  // for coordinates
  dat << DefineNodes << N_Vertices << endl;
  // for hexahedra
  dat << DefineHexa << N_Elements << endl;
  
  // for each scalar variable
  for(i=0;i<N_ScalarVar;i++)
    dat << DefineNodes << N_Vertices << endl;

  // for each vector-values variable
  for(i=0;i<N_VectorVar;i++)
    dat << DefineNodes << N_Vertices << endl;
  
  dat << endl;

  // Material = SubGridID at this place
  // Material { { ... } }
  
  // for coordinates
  dat << Nodes << endl;
  // for hexahedra
  dat << Hexahedra << endl;

  // association hexahedron <-> material goes here

  VariableCounter = VariableStart;
  // scalar data
  for(i=0;i<N_ScalarVar;i++)
  {
    dat << ScalarData << VariableCounter << endl;
    VariableCounter++;
  }

  // vector-valued data
  for(i=0;i<N_VectorVar;i++)
  {
    dat << VectorData << VariableCounter << endl;
    VariableCounter++;
  }
  dat << endl;

  VariableCounter = VariableStart;
  // scalar data
  for(i=0;i<N_ScalarVar;i++)
  {
    VarName = FEFunctionArray[i]->GetName();
    dat << ScalarField << VarName << Linear;
    dat << VariableCounter << Closed << endl;
    VariableCounter++;
  }

  // vector-valued data
  for(i=0;i<N_VectorVar;i++)
  {
    VarName = FEVectFunctArray[i]->GetName();
    dat << VectorField << VarName << Linear;
    dat << VariableCounter << Closed << endl;
    VariableCounter++;
  }
  dat << endl;

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[3*N_],Coords[3*N_+1], Coords[3*N_+2]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

#ifdef __COMPAQ__
  SwapDoubleArray(Coords, 3*N_Vertices);
#endif

  dat << At1 << endl;
  dat.write((char *)Coords,sizeof(double)*3*N_Vertices);
  dat << endl;

  m=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l]+1;
      // cout << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

#ifdef __COMPAQ__
  SwapIntArray(VertexNumbers, N_LocVertices);
#endif
  dat << endl << At2 << endl;
  dat.write((char *)(VertexNumbers), N_LocVertices*SizeOfInt);
  dat << endl;
#ifdef __COMPAQ__
  // swap back
  SwapIntArray(VertexNumbers, N_LocVertices);
#endif

  // write material data here

  DoubleArray = new double[3*N_Vertices];
  IntArray = new int[N_Vertices];

  VariableCounter = VariableStart;
  // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace3D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetType())
        {
          case Tetrahedron: 
            xi = TetraCoords[3*j];
            eta = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
          break;

          case Hexahedron: 
          case Brick: 
            xi = HexaCoords[3*j];
            eta = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
          break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        value = 0;
        for(l=0;l<N_LocDOF;l++)
          value += BFValues[l] * Coeffs[DOF[l]];
        DoubleArray[VertexNumbers[m]-1] += value;
        IntArray[VertexNumbers[m]-1]++;
        m++;
      } // endfor j
    } // endfor i

    for(i=0;i<N_Vertices;i++)
      DoubleArray[i] /= IntArray[i];

#ifdef __COMPAQ__
    SwapDoubleArray(DoubleArray, N_Vertices); 
#endif
    dat << endl << "@" << VariableCounter << endl;
    VariableCounter++;
    dat.write((char*)DoubleArray, N_Vertices*SizeOfDouble);
    dat << endl;
#ifdef __COMPAQ__
    // swap back
    SwapDoubleArray(DoubleArray, N_Vertices); 
#endif

  } // endfor k

  // write vector-valued variables into file
  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace3D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetType())
        {
          case Tetrahedron: 
            xi = TetraCoords[3*j];
            eta = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
          break;

          case Hexahedron: 
          case Brick: 
            xi = HexaCoords[3*j];
            eta = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
          break;
	  default:
	     cout<<"Handles only Tetra or Hexa " <<endl;
	  break;  
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        for(n=0;n<N_Comp;n++)
        {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l] + n*Length];
          DoubleArray[N_Comp*(VertexNumbers[m]-1) + n] += value;
        }
        IntArray[VertexNumbers[m]-1]++;
        m++;
      } // endfor j
    } // endfor i

    l = 0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
        DoubleArray[l] /= IntArray[i];
        l++;
      }
    } // endfor l
#ifdef __COMPAQ__
    SwapDoubleArray(DoubleArray, N_Vertices*N_Comp);
#endif
    dat << endl << "@" << VariableCounter << endl;
    VariableCounter++;
    dat.write((char*)DoubleArray, SizeOfDouble*N_Vertices*N_Comp);
    dat << endl;
#ifdef __COMPAQ__
    // swap back
    SwapDoubleArray(DoubleArray, N_Vertices*N_Comp);
#endif
  } // endfor k
  dat << endl;

  dat.close();
  
  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;
  delete DoubleArray;
  delete IntArray;
  delete Coords;

  OutPut("wrote output into file: " << name << endl);

  return 0;
}
/** write stored data into a VTK file */
/** start of implementation by Piotr Skrzypacz 24.03.04 */

int TOutput3D::WriteVtk(const char *name)
{
  int i,j,k,l,m,n,p;
  int N_Cells, N_Vertices, N_CellVertices, N_Comps, MaxN_VerticesPerCell;
  // *N_DOF, *N_Comp, *FESpaceNumber;
  int  N_, N_Elements, N_LocVertices;
  TVertex **Vertices, *Last, *Current;
  TFEFunction3D *fefunction;
  TFEVectFunct3D *fevectfunct;
  int N_BaseFunct, N_DOF;
  TVertex *vertex;
  TBaseCell *cell;
  double x[4], y[4], z[4];
  int *FESpaceNumber;
  TFESpace3D *fespace;
  char *Comment;
  double xi, eta, zeta;
  double value;
  double *Coeffs;
  int N_LocDOF;
  int Length, N_Comp;
  double t;
  
  TBaseFunct3D *bf;
  BaseFunct3D BaseFunct;
  FE3D FE_ID;
  double BFValues[3*MaxN_BaseFunctions3D]; // 3 for vector valued basis functions
  double *FEValues;
  int *GlobalNumbers, *BeginIndex, *Index, *DOF;
  double *Coords;
  int *VertexNumbers;
  int *NumberVertex;
  double LocValues[MaxN_BaseFunctions3D];
  double s;
  char str[15];
  int *BaseFuncts;
  int *IntArray;
  double *DoubleArray;
  // 12 - hexaedron
  // 10 - tetra
  const int CELL_TYPES_HEXA=12;
  const int CELL_TYPES_TETRA=10;

 // format (x_i, y_i, z_i)
  static double HexaCoords[] =
    { -1, -1, -1,
       1, -1, -1,
       1,  1, -1,
      -1,  1, -1,
      -1, -1,  1,  
       1, -1,  1,
       1,  1,  1,
      -1,  1,  1 
    };

  // format (x_i, y_i, z_i)
  static double TetraCoords[] =
    { 0, 0, 0,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1 };


  MaxN_VerticesPerCell = 8; // 3D case

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(4);

  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];

  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace3D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  // determine data for vtk file

  N_Elements=Coll->GetN_Cells();
//   cout << "N_Elements: " <<  N_Elements << endl;
  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  Vertices=new TVertex*[N_LocVertices];
  N_=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }
  
  // check for discontinuous scalar variables. In such a case write a new file
  // for this variable. ParaView will really display a discontinuous function
  // instead of projecting it onto P1/Q1 space.
  // However there are some drawbacks. 
  for(i = 0; i < N_ScalarVar; i++)
  {
    j = FEFunctionArray[i]->GetFESpace3D()->IsDGSpace();
    if(j==1)
    {
      // draw all discontinuous functions not just this i-th one.
      WriteVtkDiscontinuous(name,N_LocVertices,Vertices);
      break;
    }
  }
  
//   cout << "N_" << N_ << endl;
  Sort(Vertices, N_);
  //Sort(cell_types, N_);
  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  //cout << "N_Vertices: " << N_Vertices << endl;
  Coords=new double[3*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[3*N_],Coords[3*N_+1], Coords[3*N_+2]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }
  
  m=0;
  for(i=0;i<N_Elements;i++)
    {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      // cout << (int)(Current) << endl;
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      //cout << "Vertex Number: " << VertexNumbers[m] << endl;
      m++;
    } // endfor j
  } //endfor i

 


  // one additional column for absolute values of velocity
  N_Comps++;

  // to check
  //cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps << endl;
  //cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << endl;
  //cout << "N_Comps" << N_Comps << endl;
  
  dat << std::scientific;
  dat.precision(6);
  dat << "# vtk DataFile Version 4.2" << endl;
  dat << "file created by ParMooN" << endl;


  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl << endl;
  dat << "POINTS " << N_Vertices << " double" << endl;
  N_=0;
  //cout << "N_LocVertices: " << N_LocVertices << endl;
  for(i=0;i<N_Vertices;i++)
  {  
    dat << Coords[N_] << " " <<  Coords[N_+1] << " " << Coords[N_+2] << endl;
    N_ +=3;
  }
  dat << endl;
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  l=0;
  for(i=0;i<N_Elements;i++)
  {
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(j=0;j<N_CellVertices;j++)
    {
      dat << VertexNumbers[l] << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  dat << "CELL_TYPES " << N_Elements << endl;
  for(i=0;i<N_Elements;i++)
  {  
    N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
    case 4: dat << 10 << " ";
      break;
    case 8: dat << 12 << " ";
      break; 
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;  

  // function values
  
  DoubleArray = new double[3*N_Vertices];
  IntArray = new int[N_Vertices];

   // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace3D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;
    
    // will be set to true for vector valued basis functions (for example 
    // Raviart-Thomas or Brezzi-Douglas-Marini)
    bool VectOutput = false;
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      int BaseVectDim = bf->GetBaseVectDim();
      if(BaseVectDim == 3) 
        VectOutput = true;
      else if(BaseVectDim != 1)
        ErrMsg("unkown number of basis function components, assume 1.");
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
	  //  Tetrahedron
          case 4: 
	    xi   = TetraCoords[3*j];
            eta  = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
	    break;
	  // Hexahedron
	  case 8: 
            xi   = HexaCoords[3*j];
            eta  = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
	    break;
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        if(!VectOutput)
        {
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
        }
        else // VectOutput
        {
          // transform values using the Piola transform
          RefTrans3D RefTrans = TFEDatabase3D::GetRefTrans3D_IDFromFE3D(FE_ID);
          TRefTrans3D *F_K = TFEDatabase3D::GetRefTrans3D(RefTrans);
          TFEDatabase3D::SetCellForRefTrans(cell, RefTrans);
          double *BFValuesOrig = new double[3*N_LocDOF];
          switch(RefTrans)
          {
            case TetraAffin:
            case HexaAffin:
              F_K->PiolaMapOrigFromRef(N_LocDOF, BFValues, BFValuesOrig);
              break;
            case HexaTrilinear:
              ErrMsg("Piola transform for trilinear reference map not yet " << 
                     "implemented");
              break;
            default:
              ErrMsg("unknown reference transformation");
              exit(0);
              break;
          }
          
          double value_x = 0, value_y = 0, value_z = 0;
          for( l = 0; l < N_LocDOF; l++)
          {
            int face = TFEDatabase3D::GetFE3D(FE_ID)->GetFEDesc3D()
                ->GetJointOfThisDOF(l);
            int nsign = 1;
            if(face != -1)
              nsign = cell->GetNormalOrientation(face);
            value_x += BFValuesOrig[l             ] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[l +   N_LocDOF] * Coeffs[DOF[l]]*nsign;
            value_z += BFValuesOrig[l + 2*N_LocDOF] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[BaseVectDim*VertexNumbers[m] + 0] += value_x;
          DoubleArray[BaseVectDim*VertexNumbers[m] + 1] += value_y;
          DoubleArray[BaseVectDim*VertexNumbers[m] + 2] += value_z;
          
          delete [] BFValuesOrig;
        }
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    if(!VectOutput)
    {
      // non conforming
      for(i=0;i<N_Vertices;i++)
        DoubleArray[i] /= IntArray[i];
    }
    else // VectOutput
    {
      for(i = 0; i < N_Vertices; i++)
      {
        DoubleArray[3*i    ] /= IntArray[i];
        DoubleArray[3*i + 1] /= IntArray[i];
        DoubleArray[3*i + 2] /= IntArray[i];
      }
    }

    if(!VectOutput)
    {
      dat << "SCALARS " << FEFunctionArray[k]->GetName();
      dat << " double"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j] << endl;
      dat << endl;
      dat << endl;
    }
    else
    {
      // vector output, we don't write each component individually, 
      dat << "VECTORS " << FEFunctionArray[k]->GetName() << " double\n";
      for(i = 0; i < N_Vertices; i++)
      {
        for(j = 0; j < 3; j++)
        {
          dat << DoubleArray[3 * i + j] << " ";
        }
        dat << endl;
      }
      dat << endl;
      // reset
      VectOutput = false;
    }
  } // endfor k

 
  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace3D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(IntArray, 0, SizeOfInt*N_Vertices);
    m = 0;

    
    //for(k=0;k<FEVectFunctArray[i]->GetN_Components();k++)
    //{
      
    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE3D(i, cell);
      bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 4: 
	    xi   = TetraCoords[3*j];
            eta  = TetraCoords[3*j+1];
            zeta = TetraCoords[3*j+2];
	    break;

	  case 8: 
	    xi   = HexaCoords[3*j];
            eta  = HexaCoords[3*j+1];
            zeta = HexaCoords[3*j+2];
	    break;
        }
        bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
       
	for(n=0;n<N_Comp;n++)
        {
	  value = 0;
	  for(l=0;l<N_LocDOF;l++)
	    value += BFValues[l] * Coeffs[DOF[l]+n*Length];
	  DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
	} 
        IntArray[VertexNumbers[m]]++;
        m++;
      } // endfor j
    } // endfor i

    // midle value

    l = 0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
        DoubleArray[l] /= IntArray[i];
	l++;
      }
    } // endfor l
    /*
    for(i=0;i<2*N_Vertices;i++)
    {
      cout << "Do[" << i << "]" << DoubleArray[i] << endl;
    }
    */
    
    for(j=0;j<N_Comp;j++)
    {
      dat << "SCALARS " << FEVectFunctArray[k]->GetName() << j;
      dat << " double"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(i=0;i<N_Vertices;i++)
      {
	dat << DoubleArray[i*N_Comp+j] << endl;
      }
      dat << endl << endl;
    }
    
    dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
    dat << " double"<< endl;
    dat << "LOOKUP_TABLE " << "default" << endl;
    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      t=0;
      for(j=0;j<N_Comp;j++)
      {	
       t+=DoubleArray[l]*DoubleArray[l];
       l++;
      }
      dat << sqrt(t)<< endl;
    }
    dat << endl << endl;

    dat << "VECTORS " << FEVectFunctArray[k]->GetName();
    dat << " double"<< endl;
    
    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
	dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << endl;
    }
    dat << endl;
  } // endfor k

     
    // write the RegionID
    dat << endl;
    dat <<  "CELL_DATA "<< N_Elements<<endl;
    dat <<  "SCALARS RegionID int  1"<<endl;
    dat <<  "LOOKUP_TABLE default  "<<endl;
#ifdef _MPI    
    for(i=0;i<N_Elements;i++)     
      dat << (Coll->GetCell(i))->GetSubDomainNo()<<endl;
#endif    
     for(i=0;i<N_Elements;i++)
      dat << (Coll->GetCell(i))->GetRegionID() <<endl;      
    
  dat.close();
  
  delete [] IntArray;
  delete [] DoubleArray;
  delete [] NumberVertex;
  delete [] VertexNumbers;
  delete [] Vertices;
  delete [] Coords;
  delete [] FESpaceNumber;

  OutPut("wrote output into vtk file: " << name << endl);
  return 0;
}


/*
  Ulrich Wilbrandt, December 2013.

  writes an extra vtk-file for (scalar) discontinuous functions. ParaView can
  then display discontinuous data. The "POINTS" in the resulting vtk-file are 
  the vertices of the mesh, but every vertex is put here as many times as there 
  are cells this vertex belongs to. That means if a vertex belongs to 4 cells 
  in the mesh it will appear 4 times in the list "POINTS" in the resulting 
  vtk-file. The input "TVertex **Vertices" should be in this pattern. It can be
  generated by

  N_Elements=Coll->GetN_Cells();
  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
  }
  Vertices=new TVertex*[N_LocVertices];
  N_=0;
  
  for(i=0;i<N_Elements;i++)
  {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
    }
  }

  as is done in WriteVtk(cont char *name).
  WARNING: This destroys the topology of the mesh. Some filters in ParaView 
  might work incorrectly. Warp by scalar works though.
  The way this is done here is not very elegant, but I couldn't find a better
  solution.
  Also note that in each element the function is projected onto P1/Q1.
*/
void TOutput3D::WriteVtkDiscontinuous(const char *fileName,
int N_LocVertices, TVertex **Vertices)
{
  char Disc[80];             // the new file name
  strcpy(Disc,fileName);     // copy the file name ...
  strcat(Disc,"_disc.vtk");  // ... add a string to the output file name
  
  std::ofstream dat(Disc);
  if (!dat)
  {
    Error("cannot open file for output\n");
    exit(-1);
  }
  
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  dat << "# vtk DataFile Version 4.2" << endl;
  dat << "file created by MooNMD." << endl;

  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << N_LocVertices << " float" << endl;

  for(int i = 0; i < N_LocVertices; i++)
  {
    double x,y,z;
    Vertices[i]->GetCoords(x,y,z);
    dat << x << " " <<  y << " " << z << endl;
  }
  dat << endl;
  int N_Elements = Coll->GetN_Cells();
  
  // writing which vertices belong to which cells, here it is ignored that
  // a vertex might belong to multiple cells.
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  int l = 0;
  for(int i = 0; i < N_Elements; i++)
  {
    TBaseCell* current_cell = Coll->GetCell(i);
    int N_CellVertices = current_cell->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(int j = 0; j < N_CellVertices; j++)
    {
      dat << l << " ";
      l++;
    }
    dat << endl;
  }
  dat << endl;
  
  // the cell types tell paraview if this is a tetrahedron or a hexahedron
  // (export of other types is not supported here)
  dat << "CELL_TYPES " << N_Elements << endl;
  for(int i = 0; i < N_Elements; i++)
  {
    int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
    switch(N_CellVertices)
    {
    case 4: dat << 10 << " ";
      break;
    case 8: dat << 12 << " ";
      break; 
    }
  }
  dat << endl << endl;
  
  // write the function values, only for scalar functions, which includes
  // vector valued basis functions (such as Raviart-Thomas), because these 
  // are handled like scalar basis functions
  dat << "POINT_DATA " << N_LocVertices << endl;
  for(int space_number = 0; space_number < N_ScalarVar; space_number++)
  {
    TFEFunction3D* fefunction = FEFunctionArray[space_number];
    if(fefunction->GetFESpace3D()->IsDGSpace() != 1)
      continue;
    // this is a discontinuous space
    int BaseVectDim = TFEDatabase3D::GetFE3D(
      fefunction->GetFESpace3D()->GetFE3D(0,Coll->GetCell(0)))
      ->GetBaseFunct3D()->GetBaseVectDim();   // ugly, but we need to know this
    // scalar valued basis functions (normal case)
    if (BaseVectDim == 1)
    {
      dat << endl << endl;
      dat << "SCALARS " << fefunction->GetName();
      dat << " double" << endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      double *function_value = new double[1];
      for(int i = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x,y,z;
          current_cell->GetVertex(j)->GetCoords(x,y,z);
          fefunction->FindValueLocal(current_cell,i,x,y,z,function_value);
          dat << function_value[0] << endl;
        }
      }
      delete [] function_value;
    }
    // vector valued basis functions (e.g. Raviart-Thomas)
    else if(BaseVectDim == 3)
    {
      // find values for all components
      double* function_value = new double[3];   // 3==BaseVectDim
      // store function values at all vertices (three components)
      double **allValues = new double*[3];       // 3==BaseVectDim
      for(int l = 0; l < BaseVectDim; l++)
        allValues[l] = new double[N_LocVertices];
      for(int i = 0, k = 0; i < N_Elements; i++)
      {
        TBaseCell* current_cell = Coll->GetCell(i);
        int N_CellVertices = current_cell->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          double x,y,z;
          current_cell->GetVertex(j)->GetCoords(x, y, z);
          // FindValueLocal includes the necessary sign changes due to global
          // normals (for Raviart-Thomas elements)
          fefunction->FindValueLocal(current_cell,i,x,y,z,function_value);
          //for(int l = 0; l < BaseVectDim; l++)
          //  allValues[l][k] = function_value[l];
          allValues[0][k] = function_value[0];
          allValues[1][k] = function_value[1];
          allValues[2][k] = function_value[2];
          k++;
        }
      }
      delete [] function_value;
      // write the function values to the vtk-file
      dat << endl << endl;
      dat << "VECTORS " << fefunction->GetName();
      dat << " double"<< endl;
      for(int i = 0, k = 0; i < N_Elements; i++)
      {
        int N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
        for(int j = 0; j < N_CellVertices; j++)
        {
          dat << allValues[0][k] << "\t" << allValues[1][k] 
              << "\t" << allValues[2][k] << endl;
          k++;
        }
      }
      
      for(int l = 0; l < BaseVectDim; l++)
        delete [] allValues[l];
      delete [] allValues;
    }
    else
      ErrMsg("TOutput3D::WriteVtkDiscontinuous: Basis functions of dimension "
             << BaseVectDim << " are not supported");
    
    
  }
  dat.close();
}



/** write stored PARALLEL data into a pvtu and vtu files (XML files for paraview) (Sashikumaar Ganesan) */

int TOutput3D::Write_ParVTK(
#ifdef _MPI
                                MPI_Comm comm,
#endif
                               int img, char *subID)
{
  int i, j, k,l,m,n,p, rank, size, N_, N_Elements, N_LocVertices, N_BaseFunct, N_DOF;
  int AnsatzSpace, N_Cells, N_Vertices, N_CellVertices, N_Comps, MaxN_VerticesPerCell;
  int *FESpaceNumber, N_LocDOF, Length, N_Comp, *GlobalNumbers, *BeginIndex, *Index, *DOF;
  int *VertexNumbers, *NumberVertex, *BaseFuncts, begin, ID;

  double xi, eta, zeta, t, value, *Coeffs, *WArray, *DoubleArray;
  double BFValues[MaxN_BaseFunctions3D], LocValues[MaxN_BaseFunctions3D];
  double *FEValues, *Coords, s, z;
  static double HexaCoords[] = { -1, -1, -1, 1, -1, -1, 1,  1, -1, -1,  1, -1,
                                 -1, -1,  1, 1, -1,  1, 1,  1,  1, -1,  1,  1  };
  static double TetraCoords[] = { 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  char *VtkBaseName, Dquot;
  char *Comment;
  const char vtudir[] = "VTU";
  time_t rawtime;
  struct tm * timeinfo;

  TVertex **Vertices, *Last, *Current;
  TFEFunction3D *fefunction;
  TFEVectFunct3D *fevectfunct;
  TVertex *vertex;
  TBaseCell *cell;
  TFESpace3D *fespace;
  TBaseFunct3D *bf;
  BaseFunct3D BaseFunct;
  FE3D FE_ID;

  Dquot = 34; //  see ASCII Chart
  VtkBaseName = TDatabase::ParamDB->BASENAME;
  char *output_directory = TDatabase::ParamDB->OUTPUTDIR;
  AnsatzSpace = int(TDatabase::ParamDB->ANSATZ_ORDER);
#ifdef _MPI
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
#else
  rank = 0;
  size =1;
#endif

  if(rank==0)
   {
//     remove(vtudir);
    mkdir(vtudir, 0777); // create the folder to store SubDomain vtu files
   }

#ifdef _MPI
  MPI_Barrier (TDatabase::ParamDB->Comm);
#endif

  std::ostringstream os;
  os << " ";

  time ( &rawtime );
  timeinfo = localtime (&rawtime );


// write the master pvtu file
  if(rank==0)
   {
   if( TDatabase::ParamDB->SC_VERBOSE > 0 )
    OutPut("writing output into "<< output_directory << "/" << VtkBaseName 
           <<subID << "*." <<img<< " xml vtk file"<< endl);
    os.seekp(std::ios::beg);
    os << output_directory << "/" << VtkBaseName << subID;
    if(img<10) os << ".0000"<<img<<".pvtu" << ends;
    else if(img<100) os << ".000"<<img<<".pvtu" << ends;
    else if(img<1000) os << ".00"<<img<<".pvtu" << ends;
    else if(img<10000) os << ".0"<<img<<".pvtu" << ends;
    else  os << "."<<img<<".pvtu" << ends;
    std::ofstream dat(os.str().c_str());
    if (!dat)
     {
      cerr << "cannot open file for output" << endl;
#ifdef _MPI
      MPI_Abort(MPI_COMM_WORLD, 0);
#else
     exit(0);
#endif
     }

    dat <<  "<?xml version="<<Dquot<<"1.0"<<Dquot<<"?>"<< endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: Master file for parallel vtk data" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "Problem Current Time " <<TDatabase::TimeDB->CURRENTTIME << endl;
    dat <<  "  -->" << endl;

    dat << endl;


    dat <<  "<VTKFile type="<<Dquot<<"PUnstructuredGrid"<<Dquot<<" version="<<Dquot<<"0.1"
             <<Dquot<<" byte_order="<<Dquot<<"LittleEndian"<<Dquot<<">"<<endl;
    dat <<  "<PUnstructuredGrid GhostLevel="<<Dquot<<0<<Dquot<<">"<<endl;
    dat << endl;

    dat <<  " <PPoints>"<<endl;
    dat <<  "   <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
        <<"Position"<<Dquot<<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<"/>"<<endl;
    dat <<   "</PPoints>"<<endl;
    dat << endl;

    dat <<   "<PCells>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"connectivity"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"offsets"<<Dquot
        <<"      NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"UInt8"<<Dquot<<" Name="<<Dquot<<"types"<<Dquot
        <<"        NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "</PCells>"<<endl;
    dat << endl;

    dat <<   "<PPointData Vectors="<<Dquot<<"Vectors"<<Dquot<<" "<<"Scalars="<<Dquot<<"Scalars"<<Dquot<<">"<<endl;
    for(i=0;i<N_VectorVar;i++)
    dat <<   "  <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<FEVectFunctArray[i]->GetName()<<Dquot
        <<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
      dat << endl;

    for(i=0;i<N_VectorVar;i++)
    {
    for(j=0;j<FEVectFunctArray[i]->GetN_Components();j++)
     {
      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[i]->GetName()<<j<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
      dat << endl;
     }
    }

    for(i=0;i<N_ScalarVar;i++)
    dat <<   "  <PDataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<FEFunctionArray[i]->GetName()<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<"/>"<<endl;
    dat <<   "</PPointData>"<<endl;
    dat << endl;

    dat <<   "<PCellData Scalars="<<Dquot<<"SubDomainAndRegionID"<<Dquot<<">"<<endl;
#ifdef _MPI    
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<"   Name="<<Dquot<<"SubDomain"<<Dquot
        <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
#endif
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<"   Name="<<Dquot<<"RegionID"<<Dquot
        <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;

    dat <<   "</PCellData>"<<endl;
    dat << endl; 

 
    // root not take part in computation
//     begin = 1;
    begin = 0; // root take part in computation
    
    for(i=begin;i<size;i++)
     {
      ID = i; // root not take part in computation

      if(img<10)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".0000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<100)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".000"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<1000)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".00"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else if(img<10000)
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".0"<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
      else
       {
        if(i<10)        dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<100)  dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else if(i<1000) dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
        else            dat <<   "  <Piece Source="<<Dquot<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<"."<<img<<".vtu"<<Dquot<<"/>"<<endl;
       }
     }
    dat << endl;

    dat <<  "</PUnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;
    dat.close();
   } // if(rank==0
   
  // root take part in computation 
//   else
   {

    // determine data for vtu file of each processor
    MaxN_VerticesPerCell = 8; // 3D case
    FESpaceNumber = new int[N_ScalarVar+N_VectorVar];
    //cout << "N_ScalarVar: " <<  N_ScalarVar << endl;
    N_Comps = 0;
    for(i=0;i<N_ScalarVar;i++)
    {
     N_Comps++;
     fespace=FEFunctionArray[i]->GetFESpace3D();
     j=0;
     while(FESpaceArray[j]!=fespace) j++;
     FESpaceNumber[i]=j;
    }

    k = N_ScalarVar;
    for(i=0;i<N_VectorVar;i++,k++)
    {
     N_Comps += FEVectFunctArray[i]->GetN_Components();
     fespace=FEVectFunctArray[i]->GetFESpace3D();
     j=0;
     while(FESpaceArray[j]!=fespace) j++;
     FESpaceNumber[k]=j;
    }


#ifdef _MPI
  N_Elements=Coll->GetN_OwnCells();
#else
  N_Elements=Coll->GetN_Cells();
#endif


  N_LocVertices=0;
  for(i=0;i<N_Elements;i++)
   {
    cell = Coll->GetCell(i);
    N_LocVertices += cell->GetN_Vertices();
   }

  if(N_LocVertices)
   Vertices=new TVertex*[N_LocVertices];
  N_=0;
  for(i=0;i<N_Elements;i++)
   {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
     {
      Vertices[N_]=cell->GetVertex(j);
      N_++;
     }
   }


  if(N_)
   Sort(Vertices, N_);

  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }

  if(N_LocVertices)
   {
    Coords=new double[3*N_Vertices];
    VertexNumbers=new int[N_LocVertices];
    NumberVertex=new int[N_LocVertices];
   }
   Last=NULL;
   N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      Vertices[i]->GetCoords(Coords[3*N_],Coords[3*N_+1], Coords[3*N_+2]);
      k++;
      N_++;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  m=0;
  for(i=0;i<N_Elements;i++)
    {
    cell = Coll->GetCell(i);
    k=cell->GetN_Vertices();
    for(j=0;j<k;j++)
    {
      Current=cell->GetVertex(j);
      l=GetIndex(Vertices, N_LocVertices, Current);
      VertexNumbers[m]=NumberVertex[l];
      m++;
    } // endfor j
  } //endfor i


  // additional column for absolute values of velocity
  N_Comps++;

  // to check
  //cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps << endl;
  //cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << endl;
  //cout << "N_Comps" << N_Comps << endl;

  if(N_LocVertices)
   {
    DoubleArray = new double[3*N_Vertices];
    WArray = new double[N_Vertices];
   }

  ID = rank;
    os.seekp(std::ios::beg);
      if(img<10)
       {
        if(rank<10)        os<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".0000"<<img<<".vtu" <<ends;
        else if(rank<100)  os<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".0000"<<img<<".vtu" <<ends;
        else if(rank<1000) os<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".0000"<<img<<".vtu" <<ends;
        else            os<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".0000"<<img<<".vtu" <<ends;
       }
      else if(img<100)
       {
        if(rank<10)        os<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".000"<<img<<".vtu" <<ends;
        else if(rank<100)  os<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".000"<<img<<".vtu" <<ends;
        else if(rank<1000) os<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".000"<<img<<".vtu" <<ends;
        else            os<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".000"<<img<<".vtu" <<ends;
       }
      else if(img<1000)
       {
        if(rank<10)        os<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".00"<<img<<".vtu" <<ends;
        else if(rank<100)  os<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".00"<<img<<".vtu" <<ends;
        else if(rank<1000) os<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".00"<<img<<".vtu" <<ends;
        else            os<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".00"<<img<<".vtu" <<ends;
       }
      else if(img<10000)
       {
        if(rank<10)        os<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<".0"<<img<<".vtu" <<ends;
        else if(rank<100)  os<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<".0"<<img<<".vtu" <<ends;
        else if(rank<1000) os<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<".0"<<img<<".vtu" <<ends;
        else            os<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<".0"<<img<<".vtu" <<ends;
       }
      else
       {
        if(rank<10)        os<<"VTU/"<<VtkBaseName<<subID<<".000"<<ID<<"."<<img<<".vtu" <<ends;
        else if(rank<100)  os<<"VTU/"<<VtkBaseName<<subID<<".00" <<ID<<"."<<img<<".vtu" <<ends;
        else if(rank<1000) os<<"VTU/"<<VtkBaseName<<subID<<".0"  <<ID<<"."<<img<<".vtu" <<ends;
        else            os<<"VTU/"<<VtkBaseName<<subID<<"."   <<ID<<"."<<img<<".vtu" <<ends;
       }

    std::ofstream dat(os.str().c_str());
    if (!dat)
     {
      cerr << "cannot open file for output" << endl;
#ifdef _MPI
      MPI_Abort(MPI_COMM_WORLD, 0);
#else
     exit(0);
#endif
     }
    dat << setprecision(8);

    dat <<  "<?xml version="<<Dquot<<"1.0"<<Dquot<<"?>" << endl;
    dat << endl;
    dat <<  "<!--" << endl;
    dat <<  "      Title: SubDomain data for master ptvu file" << endl;
    dat <<  "    Program: ParMooN " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
    dat <<  "  -->" << endl;
    dat << endl;

    dat <<  "<VTKFile type="<<Dquot<<"UnstructuredGrid"<<Dquot<<" version="<<Dquot<<"0.1"
             <<Dquot<<" byte_order="<<Dquot<<"LittleEndian"<<Dquot<<">"<<endl;
    dat <<  "<UnstructuredGrid>"<<endl;
    dat << endl;

    dat <<  "<Piece NumberOfPoints="<<Dquot<<N_Vertices<<Dquot<<" NumberOfCells="<<Dquot<<N_Elements<<Dquot<<">"<<endl;
    dat <<  "<Points>"<<endl;
    dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<"Position"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    N_=0;
    for(i=0;i<N_Vertices;i++)
    {
     dat <<  "   " << Coords[N_] << " " <<  Coords[N_+1] << " " << Coords[N_+2] << endl;
     N_ +=3;
    }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Points>"<<endl;


    dat <<  "<Cells>"<<endl;
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"connectivity"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    l=0;
    for(i=0;i<N_Elements;i++)
     {
      N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
      dat << "       ";
      for(j=0;j<N_CellVertices;j++)
       {
        dat << VertexNumbers[l] << " ";
        l++;
       }
      dat << endl;
      }
    dat <<  "  </DataArray>"<<endl;
    dat << endl;

    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"offsets"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=1;i<=N_Elements;i++)
      {
        N_CellVertices=Coll->GetCell(i-1)->GetN_Vertices();
        dat <<  i*N_CellVertices<<"  " ;
      }
    dat <<  "   </DataArray>"<<endl;
    dat << endl;

    dat <<  "  <DataArray type="<<Dquot<<"UInt8"<<Dquot<<"  Name="<<Dquot<<"types"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=0;i<N_Elements;i++)
     {
      N_CellVertices=Coll->GetCell(i)->GetN_Vertices();
      switch(N_CellVertices)
       {
        case 4: dat << 10 << " ";
           break;
        case 8: dat << 12 << " ";
           break; 
    }
     }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Cells>"<<endl;
    dat <<endl;
    dat <<  "<PointData Vectors="<<Dquot<<"Velocity"<<Dquot<<" "<<"Scalars="<<Dquot<<"Scalars"<<Dquot<<">"<<endl;
    dat << endl;


  // write vector variables into file
    if(N_LocVertices)
    for(k=0;k<N_VectorVar;k++)
     {
      fespace = FEVectFunctArray[k]->GetFESpace3D();
      N_Comp = FEVectFunctArray[k]->GetN_Components();
      Length = FEVectFunctArray[k]->GetLength();
      Coeffs = FEVectFunctArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
      memset(WArray, 0, SizeOfDouble*N_Vertices);
      m = 0;

      for(i=0;i<N_Elements;i++)
       {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        FE_ID = fespace->GetFE3D(i, cell);
        bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
         {
          switch(cell->GetN_Vertices())
           {
            case 4: 
	      xi   = TetraCoords[3*j];
              eta  = TetraCoords[3*j+1];
              zeta = TetraCoords[3*j+2];
            break;

	    case 8: 
	      xi   = HexaCoords[3*j];
              eta  = HexaCoords[3*j+1];
              zeta = HexaCoords[3*j+2];
	    break;
           }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);

	  for(n=0;n<N_Comp;n++)
           {
	    value = 0;
	    for(l=0;l<N_LocDOF;l++)
	      value += BFValues[l] * Coeffs[DOF[l]+n*Length];
	    DoubleArray[N_Comp*VertexNumbers[m] + n] += value;
	   } 
          WArray[VertexNumbers[m]] +=1.;
          m++;
         } // endfor j
        } // endfor i


       // midle value
       l = 0;
       for(i=0;i<N_Vertices;i++)
        {
         for(j=0;j<N_Comp;j++)
          {
           if(WArray[i]>1.)
            DoubleArray[l] /= WArray[i];
	   l++;
          }
         } // endfor l


//        for(i=0;i<2*N_Vertices;i++)
//         cout << "Do[" << i << "]" << DoubleArray[i] << endl;


      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[k]->GetName()<<Dquot<<" NumberOfComponents="<<Dquot
          <<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
      l=0;
      l=0;
      for(i=0;i<N_Vertices;i++)
       {
        for(j=0;j<N_Comp;j++)
	 dat << DoubleArray[N_Comp*i+j] << " ";
       }
      dat << endl;

      dat <<  "  </DataArray>"<<endl;
      dat << endl;

      for(j=0;j<N_Comp;j++)
       {
        dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[k]->GetName()<<j<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
        for(i=0;i<N_Vertices;i++)
	 dat << DoubleArray[i*N_Comp+j] << " ";

        dat << endl;
        dat <<  "  </DataArray>"<<endl;
        dat << endl;
      }

     } // for(k=0;k<N_Vec

    // write scalar variables into file
    if(N_LocVertices)
    for(k=0;k<N_ScalarVar;k++)
     {
      fespace = FEFunctionArray[k]->GetFESpace3D();
      Coeffs = FEFunctionArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble*N_Vertices);
      memset(WArray, 0, SizeOfDouble*N_Vertices);
      m = 0;

      for(i=0;i<N_Elements;i++)
       {
        cell = Coll->GetCell(i);
        N_ = cell->GetN_Vertices();

        // find FE data for this element
        FE_ID = fespace->GetFE3D(i, cell);
        bf = TFEDatabase3D::GetFE3D(FE_ID)->GetBaseFunct3D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
         {
          switch(cell->GetN_Vertices())
           {
	    //  Tetrahedron
            case 4: 
	      xi   = TetraCoords[3*j];
              eta  = TetraCoords[3*j+1];
              zeta = TetraCoords[3*j+2];
	    break;
	    // Hexahedron
	    case 8: 
              xi   = HexaCoords[3*j];
              eta  = HexaCoords[3*j+1];
              zeta = HexaCoords[3*j+2];
	    break;
           }
          bf->GetDerivatives(D000, xi, eta, zeta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);
          value = 0;
          for(l=0;l<N_LocDOF;l++)
            value += BFValues[l] * Coeffs[DOF[l]];
          DoubleArray[VertexNumbers[m]] += value;
          WArray[VertexNumbers[m]] +=1.;
          m++;
        } // endfor j
       } // endfor i

      // non conforming
      for(i=0;i<N_Vertices;i++)
       {
        if(WArray[i]>1.)
         DoubleArray[i] /= WArray[i];
       }

      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEFunctionArray[k]->GetName()<<Dquot<<" NumberOfComponents="<<Dquot
          <<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j]<< " ";
       dat << endl;
      dat <<  "  </DataArray>"<<endl;

       dat << endl;

     }// for(k=0;k<N_ScalarVar

     
//     dat <<  "</PointData>"<<endl;
//     dat << endl;
// 
//     dat <<  "<CellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
//     dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"Region"<<Dquot
//         <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
//     for(i=0;i<N_Elements;i++)
//       dat << (Coll->GetCell(i))->GetRegionID()   << " ";
//     dat <<  "  </DataArray>"<<endl;
//     dat <<  "</CellData>"<<endl;
// 
//  
     
    dat <<  "</PointData>"<<endl;
    dat << endl;

    dat <<  "<CellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
#ifdef _MPI    
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"SubDomain"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=0;i<N_Elements;i++)     
      dat << (Coll->GetCell(i))->GetSubDomainNo() << " ";
    dat <<  "  </DataArray>"<<endl;
#endif    
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"RegionID"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
     for(i=0;i<N_Elements;i++)
      dat << (Coll->GetCell(i))->GetRegionID() << " ";       
    dat <<  "  </DataArray>"<<endl;    

    dat <<  "</CellData>"<<endl;

    dat <<  "</Piece>"<<endl;
    dat <<  "</UnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;

    dat.close();



    delete [] FESpaceNumber;

  if(N_LocVertices)
   {
    delete [] NumberVertex;
    delete [] VertexNumbers;
    delete [] Vertices;
    delete [] DoubleArray;
    delete [] WArray;
    delete [] Coords;
   }
  } //else root
//  if(rank==1)


  return 0;
} //TOutput3D::Write_ParVTK

TOutput3D::~TOutput3D()
{
  int i;

  delete [] FESpaceArray;
  delete [] FEFunctionArray;
  delete [] FEVectFunctArray;
  delete [] ParameterValues;

  for(i=0;i<N_Parameters;i++)
    free((void*) ParameterDescription[i]);

  delete [] ParameterDescription;

  if (Data) delete Data;

  free(Name);
}

TOutput3D::TOutputData::~TOutputData()
{
  if (Nodes) delete [] Nodes;
  if (ConList) delete [] ConList;
  if (FEFuncValues) 
  {
    delete [] FEFuncValues[0];
    delete [] FEFuncValues;
  }
}

void TOutput3D::ComputeOutputData()
{
  TBaseCell *Cell;
  TVertex *Vert;
  int N_Cells, N_, counter;
  
  if (Data) delete Data;
  Data = new TOutputData();

  // reset clipboard
  N_Cells = Coll->GetN_Cells();
  for(int i=0; i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();

    switch(N_)
    {
      case 4:
	Data->Type = TOutputData::TETRAHEDRON;
	break;
      case 8:
	Data->Type = TOutputData::BRICK;
	break;
    }

    for(int j=0;j<N_;j++)
    {
      Cell->GetVertex(j)->SetClipBoard(-1);
    }
  }

  // count vertices
  counter = 0;
  for(int i=0; i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();

    for(int j=0;j<N_;j++)
    {
      if( (Vert = Cell->GetVertex(j))->GetClipBoard() == -1 )
      {
	Vert->SetClipBoard(counter++);
      } 
    }
  }
  Data->N_Nodes = counter;

  // find all vertices and put into an array, create conn list
  Data->Nodes = new TVertex* [Data->N_Nodes];
  Data->ConList = new int [N_*N_Cells];
  for(int i=0;i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();
    for(int j=0;j<N_;j++)
    {
      counter = (Vert = Cell->GetVertex(j))->GetClipBoard();
      Data->Nodes[counter] = Vert;
      Data->ConList[i*N_+j] = counter+1;
    }
  }
    
  ComputeFEValues();

}

void TOutput3D::ComputeFEValues()
{
  TBaseCell *Cell;
  TFESpace3D *fespace;
  TBaseFunct3D *BaseFunc;
  FE3D FeID;
  double x, y, z, *xi, *eta, *zeta;
  double *FECoeffs;
  double BaseFuncValues[MaxN_BaseFunctions3D];
  int N_ , N_Cells, N_BF, counter;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int N_Comps, Length, index, N_VectorTotal;

  double xi_tetra[]   = {0, 1, 0, 0};
  double eta_tetra[]  = {0, 0, 1, 0};
  double zeta_tetra[] = {0, 0, 0, 1};
  double xi_brick[]   = {-1,  1,  1, -1, -1,  1,  1, -1};
  double eta_brick[]  = {-1, -1,  1,  1, -1, -1,  1,  1};
  double zeta_brick[] = {-1, -1, -1, -1,  1,  1,  1,  1};

  N_VectorTotal = 4*N_VectorVar;

  N_ = 3 + N_ScalarVar + N_VectorTotal;
  Data->N_Data = Data->N_Nodes * N_;

  N_Cells = Coll->GetN_Cells();

  switch(Data->Type)
  {
    case TOutputData::TETRAHEDRON:
      xi = xi_tetra;
      eta = eta_tetra;
      zeta = zeta_tetra;
    case TOutputData::BRICK:
      xi = xi_brick;
      eta = eta_brick;
      zeta = zeta_brick;
      break;
  }

  double *tmp = new double [Data->N_Data];
  Data->FEFuncValues = new double* [N_];
  for(int i=0;i<N_;i++)
    Data->FEFuncValues[i] = tmp + i*Data->N_Nodes;

  // coords
  for(int i=0;i<Data->N_Nodes;i++)
  {
    Data->Nodes[i]->GetCoords(x,y,z);

    Data->FEFuncValues[0][i] = x;
    Data->FEFuncValues[1][i] = y;
    Data->FEFuncValues[2][i] = z;
  }

  // find fe values (scalar)
  for(int i=0;i<N_ScalarVar;i++)
  {
    counter = 0;
    fespace = FEFunctionArray[i]->GetFESpace3D();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    FECoeffs = FEFunctionArray[i]->GetValues();

    for(int j=0;j<N_Cells;j++)
    {
      Cell = Coll->GetCell(j);
      FeID = fespace->GetFE3D(j, Cell);
      N_ = Cell->GetN_Vertices();
      DOF = GlobalNumbers + BeginIndex[j];
      BaseFunc = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
      N_BF = BaseFunc->GetDimension();

      for(int k=0;k<N_;k++)
      {
	double val = 0;
	
	if ( Cell->GetVertex(k)->GetClipBoard() == counter )
	{
	  BaseFunc->GetDerivatives(D000, xi[k], eta[k], zeta[k],
				   BaseFuncValues);
	  for(int l=0;l<N_BF;l++)
	  {
	    val += BaseFuncValues[l] * FECoeffs[DOF[l]];
	  }
	  Data->FEFuncValues[i+3][counter++] = val;
	}   
      }
    } // end for j (Cells)
  } // end for i (ScalarVar)

  // find fe values (vector)
  index = 3 + N_ScalarVar;
  for(int i=0;i<N_VectorVar;i++)
  {
    counter = 0;
    fespace = FEVectFunctArray[i]->GetFESpace3D();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    N_Comps = FEVectFunctArray[i]->GetN_Components();
    Length = FEVectFunctArray[i]->GetLength();
    FECoeffs = FEVectFunctArray[i]->GetValues();

    for(int j=0;j<N_Cells;j++)
    {	
      Cell = Coll->GetCell(j);
      FeID = fespace->GetFE3D(j, Cell);
      N_ = Cell->GetN_Vertices();
      DOF = GlobalNumbers + BeginIndex[j];
      BaseFunc = TFEDatabase3D::GetBaseFunct3DFromFE3D(FeID);
      N_BF = BaseFunc->GetDimension();

      for(int k=0;k<N_;k++)
      {
	if ( Cell->GetVertex(k)->GetClipBoard() == counter )
	{
	  double abs = 0;

	  BaseFunc->GetDerivatives(D000, xi[k], eta[k], zeta[k],
				   BaseFuncValues);

	  for(int p=0;p<N_Comps;p++)
	  {
	    double val = 0;
	    for(int l=0;l<N_BF;l++)
	    {
	      val += BaseFuncValues[l] * FECoeffs[p*Length+DOF[l]];
	    }
	    abs += val*val;
	    Data->FEFuncValues[index+p][counter] = val;
	  } // end for p (Components)
	  Data->FEFuncValues[index+N_Comps][counter] = sqrt(abs);
	  counter++;
	}
	
      } // end for k (Vertices)
    } // end for j (Cells)
    index += N_Comps + 1;
  } // end for i (VectorVar)
}

int TOutput3D::WriteBinaryPlt(const char *filename)
{
  int ret, N_Cells, N_Comps;
  int Zero=0, One=1;
  int PassiveVarList[3+N_ScalarVar+4*N_VectorVar];

  char Title[] = "tecplot file created by MooNMD";
  char ScratchDir[] = ".";
  char *ZoneTitle = Name;
  char *Variables;

  cout << "writing tecplot file ...";

  N_Cells = Coll->GetN_Cells();

  memset(PassiveVarList, 0, sizeof(PassiveVarList));

  ComputeOutputData();

  // variables names
  std::ostringstream os;
  std::string string;
  os << "X, Y, Z";
  for(int i=0;i<N_ScalarVar;i++)
  {
    os << ", " << FEFunctionArray[i]->GetName() ;
  }
  for(int i=0;i<N_VectorVar;i++)
  {
    N_Comps = FEVectFunctArray[i]->GetN_Components();
    for(int j=0;j<N_Comps;j++)
    {
      os << ", " << FEVectFunctArray[i]->GetName() << "_" << j;
    }
    os << ", |" << FEVectFunctArray[i]->GetName() << "|";
  }
  os << ends << std::flush;
  string = os.str();

  Variables = (char*) string.c_str();
  
  ret = TECINI112(Title, Variables, (char*) filename, ScratchDir, &Zero,
		  &Zero, &One);
  if ( ret == -1 )
  {
    cerr << "TECINI112() failed !" << endl;
    return -1;
  }

  ret = TECZNE112(ZoneTitle, (int*) &Data->Type, &Data->N_Nodes, &N_Cells,
		  &Zero, &Zero, &Zero, &Zero,
		  &TDatabase::TimeDB->CURRENTTIME, &Zero, &Zero, &One, 
		  &Zero, &Zero, &Zero, &Zero, &Zero, PassiveVarList,
		  NULL, NULL, &Zero); 
  if ( ret == -1 )
  {
    cerr << "TECZNE112() failed !" << endl;
    TECEND();
    return -1;
  }

  ret = TECDAT112(&Data->N_Data, Data->FEFuncValues[0], &One);
  if ( ret == -1 )
  {
    cerr << "TECDAT112() failed !" << endl;
    TECEND112();
    return -1;
  }

  ret = TECNOD112(Data->ConList);
  if ( ret == -1 )
  {
    cerr << "TECNOD112() failed !" << endl;
    TECEND112();
    return -1;
  }

  TECEND112();

  cout << " done " << endl;

  return 0;
}
