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
// @(#)Output2D.C        1.7 05/05/00
//
// Class:       TOutput2D
// Purpose:     store given data and realize output
//
// Author:      Gunar Matthies (21.08.1998)
//
// History:     start of implementation 21.08.1998 (Gunar Matthies)
//        :     parallel vtk output 26.09.09 (Sashikumaar Ganesan)
// =======================================================================

#include <FESpace2D.h>
#include <FEFunction2D.h>
#include <FEVectFunct2D.h>
#include <FEDatabase2D.h>
#include <Output2D.h>
#include <BaseCell.h>
#include <Joint.h>
#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __COMPAQ__
// include byte reordering routines for Compaq only
#include <Utilities.h>
#endif

#  include <sstream>
// #  include <malloc.h>
#  include <dirent.h> 
#  include <unistd.h>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#  include <sys/stat.h>
#  include <sys/types.h>
#include <math.h>

#include <QuadAffin.h>
#include <TriaAffin.h>
#include <QuadBilinear.h>

// tecplot
#include <TECIO.h>

/** constructor maximum number of these things */
TOutput2D::TOutput2D(int maxn_fespaces, int maxn_scalar, int maxn_vect,
int maxn_parameters, TDomain *domain)
{
  MaxN_FESpaces=maxn_fespaces;
  N_FESpaces=0;

  MaxN_ScalarVar=maxn_scalar;
  N_ScalarVar=0;

  MaxN_VectorVar=maxn_vect;
  N_VectorVar=0;

  MaxN_Parameters=maxn_parameters;
  N_Parameters=0;

  FESpaceArray=new TFESpace2D*[MaxN_FESpaces];

  FEFunctionArray=new TFEFunction2D*[MaxN_ScalarVar];

  FEVectFunctArray=new TFEVectFunct2D*[MaxN_VectorVar];

  ParameterValues=new double[MaxN_Parameters];

  ParameterDescription=new const char*[MaxN_Parameters];

  Coll = NULL;

  Domain = domain;
  
  Data = NULL;
}


/** add a FESpace into this output object (internal use) */
int TOutput2D::AddFESpace(TFESpace2D *fespace)
{
  TFESpace2D **NewStorage;
  int i;

  if(N_FESpaces==0)
  {
    // no fespace in this output object
    Coll=fespace->GetCollection();
  }

  // check whether fespace is already stored
  for(i=0;i<N_FESpaces;i++)
    if(FESpaceArray[i]==fespace) return i+1;

  // this space is based on a different collection
  if(fespace->GetCollection()!=Coll) return 0;

  if(MaxN_FESpaces<=N_FESpaces)
  {
    // enlarge storage
    NewStorage=new TFESpace2D*[MaxN_FESpaces+5];
    memcpy(NewStorage, FESpaceArray, sizeof(TFESpace2D *)*MaxN_FESpaces);
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
int TOutput2D::AddFEFunction(TFEFunction2D *fefunction)
{
  TFEFunction2D **NewStorage;

  if(!(AddFESpace(((TFEFunction2D *)fefunction)->GetFESpace2D()))) return 0;

  if(MaxN_ScalarVar<=N_ScalarVar)
  {
    // enlarge storage
    NewStorage=new TFEFunction2D*[MaxN_ScalarVar+5];
    memcpy(NewStorage, FEFunctionArray, MaxN_ScalarVar*sizeof(TFEFunction2D *));
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
int TOutput2D::AddFEVectFunct(TFEVectFunct2D *fevectfunct)
{
  TFEVectFunct2D **NewStorage;

  if(!(AddFESpace(fevectfunct->GetFESpace2D()))) return 0;

  if(MaxN_VectorVar<=N_VectorVar)
  {
    // enlarge storage
    NewStorage=new TFEVectFunct2D*[MaxN_VectorVar+5];
    memcpy(NewStorage, FEVectFunctArray,
      MaxN_VectorVar*sizeof(TFEVectFunct2D *));
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
int TOutput2D::AddParameter(double value, const char *descr)
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
  int i, j, *rr, len;
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

  delete [] rr;
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
int TOutput2D::Write(std::string basename, int i, double _current_time)
{
  std::ostringstream os;

  os.seekp(std::ios::beg);
  os << basename << i << ".vtk"<<ends;
  os.seekp(std::ios::beg);
  os << basename << i << ".vtk"<<ends;
  if(TDatabase::ParamDB->WRITE_VTK)
  {
    cout << " Output2D:: writing " << os.str() << endl;
    WriteVtk(os.str().c_str());
  }
  
  os.seekp(std::ios::beg);
  os << basename << i << ".gnu"<<ends;
  if(TDatabase::ParamDB->WRITE_GNU) 
  {
    cout << " Output2D:: writing " << os.str() << endl;
    WriteGnuplot(os.str().c_str());
  }
  // add more here if needed
  return 0;
}

/** write stored data into a grape file */
int TOutput2D::WriteGrape(const char *name)
{
  int i,j,k,l,m,p,N_;
  int headerlength;
  int *header, first[3];
  int format, version, N_Elements, N_Vertices;
  int N_LocVertices, *N_DOF, *N_Comp, *FESpaceNumber;
  TBaseCell *cell;
  TVertex **Vertices, *Last, *Current;
  int *NumberVertex;
  double *ParamPtr;
  TFESpace2D *fespace;
  TFEFunction2D *fefunction;
  TFEVectFunct2D *fevectfunct;
  double *Coords;
  int *VertexNumbers;
  int *BaseFuncts;
  FE2D FE_ID;
  BaseFunct2D BaseFunct;
  #ifdef __3D__
  double z;
  #endif

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }

  format=1;
  version=5;                     // for new FE concept
  headerlength=48+20*(N_FESpaces+N_ScalarVar+N_VectorVar+N_Parameters);
  // cout << "Header length: " << headerlength << endl;
  header=new int[headerlength];
  first[0]=format;
  first[1]=version;
  first[2]=headerlength;
  ParameterValues[0] = TDatabase::TimeDB->CURRENTTIME;

  N_Elements=Coll->GetN_Cells();
  // cout << "number of elements: " << N_Elements << endl;

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
  header[0]=2;                   // dimension
  header[1]=N_Elements;
  header[2]=N_Vertices;
  header[3]=N_LocVertices;
  header[4]=N_FESpaces;
  header[5]=N_ScalarVar;
  header[6]=N_VectorVar;
  header[7]=N_Parameters;


  #ifdef __COMPAQ__
  SwapIntArray(header, 8);
  #endif



  strncpy((char *)(header+8),"Output generated by 'MooN_MD' program",160);

  // cout << "header[0]= " << header[0] << endl;
  // cout << "header[1]= " << header[1] << endl;
  // cout << "header[2]= " << header[2] << endl;
  // cout << "header[3]= " << header[3] << endl;
  // cout << "header[4]= " << header[4] << endl;
  // cout << "header[5]= " << header[5] << endl;
  // cout << "header[6]= " << header[6] << endl;
  // cout << "header[7]= " << header[7] << endl;

  // putting in parameter information
  N_=48;
  for(i=0;i<N_Parameters;i++)
  {
    ParamPtr=(double *)(header+N_);
    *ParamPtr=ParameterValues[i];


    #ifdef __COMPAQ__
    SwapDoubleArray(ParamPtr, 1);
    #endif


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
    fespace=FEFunctionArray[i]->GetFESpace2D();
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
    fespace=FEVectFunctArray[i]->GetFESpace2D();
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

  Coords=new double[2*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];

  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
      #ifdef __3D__
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1], z);
      #else
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1]);
      #endif
      // cout << N_/2 << "  ";
      // cout << Coords[N_] << "       " << Coords[N_+1] << endl;
      k++;
      N_ += 2;
      Last=Current;
    }
    NumberVertex[i]=k;
  }

  #ifdef __COMPAQ__
  SwapDoubleArray(Coords, 2*N_Vertices);
  #endif

  dat.write((char *)Coords,sizeof(double)*2*N_Vertices);
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
    }                            // endfor j
  }                              //endfor i



  #ifdef __COMPAQ__
  SwapIntArray(VertexNumbers, N_LocVertices);
  #endif



  dat.write((char *)VertexNumbers,sizeof(int)*N_LocVertices);
  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;

  BaseFuncts=new int[N_Elements];
  for(p=0;p<N_FESpaces;p++)
  {
    fespace=FESpaceArray[p];

    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);

      FE_ID = fespace->GetFE2D(i, cell);
      BaseFunct = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D_ID();
      BaseFuncts[i] = BaseFunct;
    }                            // endfor i


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



  }                              // endfor p
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

  cout << endl;
  cout << "wrote output into file: " << name << endl;

  delete N_DOF;
  delete N_Comp;
  delete FESpaceNumber;

  return 0;
}


/** write stored data into a gunplot file */
int TOutput2D::WriteGnuplot(const char *name)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Vertices, N_Comps, MaxN_VerticesPerCell;
  int N_BaseFunct, N_DOF;
  TBaseCell *cell;
  double x[4], y[4];
  double *Data;
  int *FESpaceNumber;
  TFESpace2D *fespace;
  char *Comment;
  double xi[4], eta[4];
  TBaseFunct2D *bf;
  FE2D FE_ID;
  double BFValues[4][MaxN_BaseFunctions2D];
  double *FEValues;
  int *GlobalNumbers, *BeginIndex, *Index;
  double LocValues[MaxN_BaseFunctions2D];
  double s;
  char str[15];
  #ifdef __3D__
  double z;
  #endif

  MaxN_VerticesPerCell = 4;      // 2D case

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(8);

  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];

  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  // one additional column for absolute values of velocity
  N_Comps++;

  Data = new double[MaxN_VerticesPerCell*N_Comps];

  dat << "# gnuplot data generated by MooNMD" << endl;
  dat << "#" << setw(11) << "x" << setw(12) << "y";
  for(i=0;i<N_ScalarVar;i++)
    dat << setw(12) << FEFunctionArray[i]->GetName();

  for(i=0;i<N_VectorVar;i++)
  {
    Comment = FEVectFunctArray[i]->GetName();
    for(k=0;k<FEVectFunctArray[i]->GetN_Components();k++)
      dat << setw(10) << Comment << "-" << setw(1) << k+1;
  }
  if(N_VectorVar>0)
  {
    str[0] = '|'; str[1] = '\0';
    Comment = FEVectFunctArray[0]->GetName();
    strncat(str, Comment, 10);
    strcat(str,"|");
    dat << setw(12) << str;
  }

  dat << endl;
  dat << "#" << endl;

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Vertices = cell->GetN_Vertices();
    switch(N_Vertices)
    {
      case 3:
        xi[0] = 0; xi[1] = 1; xi[2] = 0;
        eta[0] = 0; eta[1] = 0; eta[2] = 1;
        break;

      case 4:
        xi[0] = -1; xi[1] = 1; xi[2] = 1; xi[3] = -1;
        eta[0] = -1; eta[1] = -1; eta[2] = 1; eta[3] = 1;
        break;
    }                            // endswitch

    for(j=0;j<N_Vertices;j++)
    {
      #ifdef __3D__
      cell->GetVertex(j)->GetCoords(x[j],y[j], z);
      #else
      cell->GetVertex(j)->GetCoords(x[j],y[j]);
      #endif
    }

    for(j=0;j<N_ScalarVar;j++)
    {
      fespace = FEFunctionArray[j]->GetFESpace2D();
      FEValues = FEFunctionArray[j]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
      Index = GlobalNumbers+BeginIndex[i];
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
      N_BaseFunct = bf->GetDimension();
      for(k=0;k<N_Vertices;k++)
      {
        bf->GetDerivatives(D00, xi[k], eta[k], BFValues[k]);
        bf->ChangeBF(Coll, cell, BFValues[k]);
      }
      for(k=0;k<N_BaseFunct;k++)
        LocValues[k] = FEValues[Index[k]];

      for(k=0;k<N_Vertices;k++)
      {
        s = 0;
        for(l=0;l<N_BaseFunct;l++)
          s += LocValues[l]*BFValues[k][l];
        Data[j*MaxN_VerticesPerCell+k] = s;
      }                          // endfor k
    }                            // endfor j

    n = N_ScalarVar;
    for(j=0;j<N_VectorVar;j++)
    {
      fespace = FEVectFunctArray[j]->GetFESpace2D();
      N_DOF = fespace->GetN_DegreesOfFreedom();
      FEValues = FEVectFunctArray[j]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
      Index = GlobalNumbers+BeginIndex[i];
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
      N_BaseFunct = bf->GetDimension();
      for(k=0;k<N_Vertices;k++)
      {
        bf->GetDerivatives(D00, xi[k], eta[k], BFValues[k]);
        bf->ChangeBF(Coll, cell, BFValues[k]);
      }

      for(m=0;m<FEVectFunctArray[j]->GetN_Components();m++)
      {
        for(k=0;k<N_BaseFunct;k++)
          LocValues[k] = FEValues[m*N_DOF+Index[k]];

        for(k=0;k<N_Vertices;k++)
        {
          s = 0;
          for(l=0;l<N_BaseFunct;l++)
            s += LocValues[l]*BFValues[k][l];
          Data[n*MaxN_VerticesPerCell+k] = s;
        }                        // endfor k
        n++;
      }                          // endfor m
      if(j==0)
        for(k=0;k<N_Vertices;k++)
          Data[k+MaxN_VerticesPerCell*(N_Comps-1)] =
            sqrt(  Data[k+MaxN_VerticesPerCell*N_ScalarVar]
            * Data[k+MaxN_VerticesPerCell*N_ScalarVar]
            + Data[k+MaxN_VerticesPerCell*(N_ScalarVar+1)]
            * Data[k+MaxN_VerticesPerCell*(N_ScalarVar+1)]);
    }                            // endfor j

    for(j=0;j<N_Vertices;j++)
    {
      dat << setw(12) << x[j] << setw(12) << y[j];
      for(k=0;k<N_Comps;k++)
        dat << setw(12) << Data[k*MaxN_VerticesPerCell+j];
      dat << endl;
    }                            // endfor j
    dat << setw(12) << x[0] << setw(12) << y[0];
    for(k=0;k<N_Comps;k++)
      dat << setw(12) << Data[k*MaxN_VerticesPerCell+0];
    dat << endl << endl << endl; // two blank lines
  }                              // endfor i

  dat.close();

  delete Data;

  cout << endl;
  cout << "wrote output into gnuplot file: " << name << endl;

  return 0;
}



/** write stored data into a VTK file */
/** start of implementation by Piotr Skrzypacz 24.03.04 */
/** modified by Sashi  26.09.09 **/

int TOutput2D::WriteVtk(const char *name)
{
  int i,j,k,l,m,n;
  int N_Vertices, N_CellVertices, N_Comps;
  int  N_, N_Elements, N_LocVertices;
  int *FESpaceNumber;
  int N_LocDOF, Length, N_Comp;
  int *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers, *NumberVertex;

  double xi, eta, value, *Coeffs, t;
  double BFValues[MaxN_BaseFunctions2D];
  double *Coords, *WArray, *DoubleArray;
  /** this is needed to handle vector valued basis functions such as 
   * Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM) element */
  double BFValuesOrig[MaxN_BaseFunctions2D];
  bool VectOutput = false;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  int BaseVectDim;
  double value_y;
  int edge, nsign;
  
  double QuadCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1};
  double TriaCoords[] = { 0, 0, 1, 0,  0, 1};

  char *variable;
  char var1[2] = {'r', 'z'};
  char var2[3] = {'1', '2','3'};  
  
  if(TDatabase::ParamDB->Axial3D)
   { variable = var1;}
  else
   { variable = var2;}
  
  time_t rawtime;
  struct tm * timeinfo;
  
  TVertex **Vertices, *Last, *Current;
  TBaseCell *cell;
  TFESpace2D *fespace;
  TBaseFunct2D *bf;
  FE2D FE_ID;

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(9);

  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];
   //cout << "N_ScalarVar: " <<  N_ScalarVar << endl;
  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  // determine data for vtk file

  N_Elements=Coll->GetN_Cells();
  // 
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
  // instead of projecting it onto P1/Q1 space. However there are some
  // drawbacks.
  for(i=0;i<N_ScalarVar;i++)
  {
    j = FEFunctionArray[i]->GetFESpace2D()->IsDGSpace();
    if(j==1)
    {
      WriteVtkDiscontinuous(name, N_LocVertices, Vertices);
      break;
    }
  }
  
  //cout << "N_" << N_ << endl;
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
  Coords=new double[2*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
#ifdef __3D__
      double z;
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1], z);
#else
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1]);
#endif
      //cout << N_/2 << "  ";
      //cout << Coords[N_] << "       " << Coords[N_+1] << endl;
      k++;
      N_ += 2;
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
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  
//     dat << "File created by MooNMD"<<  "Date & Time: " <<asctime (timeinfo) << " Computational Time < " <<TDatabase::TimeDB->CURRENTTIME <<" >" 
// #ifdef __PSD__    
//   << " L < " << TDatabase::ParamDB->REACTOR_P29 
// #endif  
//   <<" >" <<endl;
  
  
  dat << "# vtk DataFile Version 4.2" << endl;
  dat <<"File created by MooNMD, "<< "Computational Time : " << TDatabase::TimeDB->CURRENTTIME <<", " 
#ifdef __PSD__    
  << " L < " << TDatabase::ParamDB->REACTOR_P29 <<" >"
#endif    
  << "Date & Time: " <<asctime (timeinfo) <<endl;


  dat << "ASCII" << endl;
  dat << "DATASET UNSTRUCTURED_GRID" << endl;
  dat << "POINTS " << N_Vertices << " float" << endl;
  N_=0;
  //cout << "N_LocVertices: " << N_LocVertices << endl;
  for(i=0;i<N_Vertices;i++)
  {
    dat << Coords[N_] << " " <<  Coords[N_+1] << " ";
    dat << double(0) << endl;
    N_ +=2;
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
    case 4: dat << 9 << " ";
      break;
    case 3: dat << 5 << " ";
      break;
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_Vertices << endl;

  // function values
  DoubleArray = new double[2*N_Vertices];
  WArray = new double[N_Vertices];

   // write scalar variables into file
  for(k=0;k<N_ScalarVar;k++)
  {
    fespace = FEFunctionArray[k]->GetFESpace2D();
    Coeffs = FEFunctionArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

    // get dimension of basis functions
    cell = Coll->GetCell(0);
    FE_ID = fespace->GetFE2D(i, cell);
    bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
    BaseVectDim = bf->GetBaseVectDim();
    N_Comp = BaseVectDim;

    memset(DoubleArray, 0, SizeOfDouble*2*N_Vertices);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;
    
    // set to true if basis functions are vectors
    if (BaseVectDim>1)  
      VectOutput = true;
    else 
      VectOutput = false;

    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      
      RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FE_ID);
      //BF2DRefElements RefElement = TFEDatabase2D::GetRefElementFromFE2D(FE_ID);
      F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
      switch(RefTrans)
      {
        case TriaAffin:
          ((TTriaAffin*)F_K)->SetCell(cell);
          break;
        case QuadAffin:
          ((TQuadAffin*)F_K)->SetCell(cell);
          break;
        case QuadBilinear:
          ((TQuadBilinear*)F_K)->SetCell(cell);
          break;
        default: 
          cout << "Output2D():: no such reference transformation allowed" 
               << endl;
          break;
      }                                           // endswitch

//       for(j=0; j<N_LocDOF; j++)
// 	if(Coeffs[DOF[l]]>1)
//        cout <<i<< " SolPbe " <<Coeffs[DOF[l]] <<endl;

      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
	    break;

	  case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
	    break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);
        value = 0;
        //for(l=0;l<N_LocDOF;l++)
        //  value += BFValues[l] * Coeffs[DOF[l]];
  
        //DoubleArray[VertexNumbers[m]] += value;
        value_y = 0;
        if(!VectOutput)
        {                                         // standard
          for(l=0;l<N_LocDOF;l++)
          {
            value += BFValues[l] * Coeffs[DOF[l]];
          }
          DoubleArray[VertexNumbers[m]] += value;
        }
        else // VectOutput
        {
          // apply Piola transform 
          switch(RefTrans)
          {
            case TriaAffin:
            case QuadAffin:
              F_K->PiolaMapOrigFromRef(N_LocDOF,BFValues,BFValuesOrig);
              break;
            case QuadBilinear: 
              // for non affine reference transformations one needs to know the 
              // point, because the determinant is not constant in this case.
              ((TQuadBilinear*)F_K)->PiolaMapOrigFromRefNotAffine(
                                        N_LocDOF,BFValues,BFValuesOrig,xi,eta);
              break;
            default:
              cout << "Output2D():: no such reference transformation allowed" 
               << endl;
              break;
          }
          for(l=0;l<N_LocDOF;l++)
          {
            // change sign of basis functions according to global normal
            edge=TFEDatabase2D::GetFE2D(FE_ID)->GetFEDesc2D()->GetJointOfThisDOF(l);
            if (edge != -1)
            {
              nsign = cell->GetNormalOrientation(edge);
            }
            else
              nsign=1;
            value += BFValuesOrig[l] * Coeffs[DOF[l]]*nsign;
            value_y += BFValuesOrig[N_LocDOF+l] * Coeffs[DOF[l]]*nsign;
          }
          DoubleArray[N_Comp*VertexNumbers[m] + 0] += value;
          DoubleArray[N_Comp*VertexNumbers[m] + 1] += value_y;
        }
        
        WArray[VertexNumbers[m]] +=1.;
        m++;

      } // endfor j
    } // endfor i

    // non conforming
    if (!VectOutput)
    {
      for(i=0;i<N_Vertices;i++)
      {
        if(WArray[i]!=0.)
        {
          DoubleArray[i] /= WArray[i];
        }
      }
    }
    else // VectOutput
    {
      l = 0;
      for(i=0;i<N_Vertices;i++)
      {
        for(j=0;j<N_Comp;j++)
        {
          if(WArray[i]!=0.)
            DoubleArray[l] /= WArray[i];
          l++;
        }
      }                                           // endfor i
    }


//       for(j=0; j<N_Vertices; j++)
//      if(DoubleArray[j]>1)
//        cout <<j << " SolPbe " <<DoubleArray[j] <<endl;
     
     
    if (!VectOutput)
    {
      // standard output writing
      dat << "SCALARS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(j=0;j<N_Vertices;j++)
        dat << DoubleArray[j] << endl;
      dat << endl;
      dat << endl;
    }
    else // VectOutput
    {
      // scalar components
      for(j=0;j<N_Comp;j++)
      {
        dat << "SCALARS " << FEFunctionArray[k]->GetName() << j;
        dat << " float"<< endl;
        dat << "LOOKUP_TABLE " << "default" << endl;
        for(i=0;i<N_Vertices;i++)
        {
          dat << DoubleArray[i*N_Comp+j] << endl;
        }
        dat << endl << endl;
      }

      // absolute value
      dat << "SCALARS " << "|" << FEFunctionArray[k]->GetName() << "|";
      dat << " float"<< endl;
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

      dat << "VECTORS " << FEFunctionArray[k]->GetName();
      dat << " float"<< endl;

      l=0;
      for(i=0;i<N_Vertices;i++)
      {
        for(j=0;j<N_Comp;j++)
        {
          dat << DoubleArray[N_Comp*i+j] << " ";
        }
        dat << double(0) << " " << endl;
      }
      dat << endl;
      // reset for next iteration in loop over all scalar variables
      VectOutput = false;
    }
  } // endfor k (loop over scalar variables


  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace2D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

   // cout << "N_Comp  " << N_Comp << endl;
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
            break;
          case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
            break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

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

    // mean value

    l = 0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
        if(WArray[i]!=0.)
          DoubleArray[l] /= WArray[i];
        l++;
      }
    } // endfor l

    for(j=0;j<N_Comp;j++)
    {
      dat << "SCALARS " << FEVectFunctArray[k]->GetName() <<variable[j];
      dat << " float"<< endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      for(i=0;i<N_Vertices;i++)
      {
        dat << DoubleArray[i*N_Comp+j] << endl;
      }
      dat << endl << endl;
    }

    dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
    dat << " float"<< endl;
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
    dat << " float"<< endl;

    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
        dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << double(0) << " " << endl;
    }
    dat << endl;
  } // endfor k

  dat << endl;

//   delete [] FESpaceNumber;
//   delete [] NumberVertex;
//   delete [] VertexNumbers;
//   delete [] Vertices;
//   delete [] DoubleArray;
//   delete [] WArray;
//   delete [] Coords;

  dat.close();
//   if( TDatabase::ParamDB->SC_VERBOSE > 0)
    OutPut("wrote output into vtk file: " << name << endl);
  return 0;
}


/*
  Ulrich Wilbrandt, November 2011.

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
void TOutput2D::WriteVtkDiscontinuous(const char *fileName,
int N_LocVertices, TVertex **Vertices)
{
  double x,y;                // coordinates of a vertex
  int N_Elements;            // number of all elements in this mesh
  int N_CellVertices;        // number of all vertices in this cell
  int i,j,k,l;               // loop variables
  TFEFunction2D *fefunction; // this function, which is discontinuous
  TBaseCell *current_cell;   
  double *function_value;    // value of function at a particular vertex
  double **allValues;        // in case of vector valued basis functions, store
                             // all values of all components
  int space_number;          
  char Disc[80];             // the new file name
  strcpy(Disc,fileName);     // copy the file name ...
  strcat(Disc,"_disc.vtk");  // ... add a string to the output file name

  N_Elements = Coll->GetN_Cells();

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

  for(i=0;i<N_LocVertices;i++)
  {
    #ifdef __3D__
    //Vertices[i]->GetCoords(x,y,z);
    OutPut("WriteVtkDiscontinuous has to be checked for 3D"<<endl);
    exit(4711);
    #else
    Vertices[i]->GetCoords(x,y);
    #endif
    dat << x << " " <<  y << " " << double(0) << endl;
    // in 2D: set last entry 0
  }
  dat << endl;
  dat << "CELLS " << N_Elements << " " <<  N_Elements+N_LocVertices << endl;
  l=0;
  for(i=0;i<N_Elements;i++)
  {
    current_cell = Coll->GetCell(i);
    N_CellVertices = current_cell->GetN_Vertices();
    dat <<  N_CellVertices << " ";
    for(j=0;j<N_CellVertices;j++)
    {
      dat << l << " ";
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
      case 4: dat << 9 << " ";
      break;
      case 3: dat << 5 << " ";
      break;
    }
  }
  dat << endl << endl;
  dat << "POINT_DATA " << N_LocVertices << endl;
  for(space_number=0; space_number<N_ScalarVar; space_number++)
  {
    fefunction = FEFunctionArray[space_number];
    if(fefunction->GetFESpace2D()->IsDGSpace() != 1)
      continue;
    // this is a discontinuous space

    int BaseVectDim =
      TFEDatabase2D::GetFE2D(
      fefunction->GetFESpace2D()->GetFE2D(0,Coll->GetCell(0)))
      ->GetBaseFunct2D()->GetBaseVectDim();   // ugly, but we need to know this
    if (BaseVectDim==1)
    {
      dat << endl << endl;
      dat << "SCALARS " << fefunction->GetName();
      dat << " double" << endl;
      dat << "LOOKUP_TABLE " << "default" << endl;
      function_value = new double[3];  // only need to allocate this first entry
      for(i=0;i<N_Elements;i++)
      {
        current_cell = Coll->GetCell(i);
        N_CellVertices=current_cell->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
          #ifdef __2D__
          current_cell->GetVertex(j)->GetCoords(x,y);
          #endif
          fefunction->FindValueLocal(current_cell,i,x,y,function_value);
          dat << function_value[0] << endl;
        }
      }
    }
    else if(BaseVectDim==2)
    {
      // find values for all components
      function_value = new double[2];   // 2==BaseVectDim
      allValues = new double*[2];       // 2==BaseVectDim
      for(l=0; l<BaseVectDim; l++)
        allValues[l] = new double[N_LocVertices];
      k=0;
      for(i=0;i<N_Elements;i++)
      {
        current_cell = Coll->GetCell(i);
        N_CellVertices=current_cell->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
          #ifdef __2D__
          current_cell->GetVertex(j)->GetCoords(x,y);
          #endif
          // FindValueLocal includes the necessary sign changes due to global
          // normals (for Raviart-Thomas elements)
          fefunction->FindValueLocal(current_cell,i,x,y,function_value);
          //for(l=0; l<BaseVectDim; l++)
          //  allValues[l][k] = function_value[l];
          allValues[0][k] = function_value[0];
          allValues[1][k] = function_value[1];
          k++;
        }
      }
      delete [] function_value;
      // write the function values to the vtk-file
      for(l=0; l<BaseVectDim; l++)
      {
        dat << endl << endl;
        dat << "SCALARS " << fefunction->GetName() << l;
        dat << " double" << endl;
        dat << "LOOKUP_TABLE " << "default" << endl;
        k=0;
        for(i=0;i<N_Elements;i++)
        {
          N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
          for(j=0;j<N_CellVertices;j++)
          {
            dat << allValues[l][k] << endl;
            k++;
          }
        }
      }
      dat << endl << endl;
      dat << "VECTORS " << fefunction->GetName();
      dat << " double"<< endl;
      k=0;
      for(i=0;i<N_Elements;i++)
      {
        N_CellVertices = Coll->GetCell(i)->GetN_Vertices();
        for(j=0;j<N_CellVertices;j++)
        {
          dat << allValues[0][k] << "\t" << allValues[1][k] // 2D 
              << "\t" << double(0) << endl;
          k++;
        }
      }
      for(l=0; l<BaseVectDim; l++)
        delete [] allValues[l];
      delete [] allValues;
    }
    else
      OutPut("TOutput2D::WriteVtkDiscontinuous: Basis functions of dimension "
        << BaseVectDim << " are not supported." << endl);
  }
  dat.close();
}




/******************************************************************************/
/*                                                                            */
/** write stored data into a MATLAB file                                      */
/** this works only for regular grids on rectangular domains !!!              */
/** this works only for one scalar field !!!                                  */
/*                                                                            */
/******************************************************************************/
#ifdef __2D__


int TOutput2D::WriteMatlab(const char *name)
{
  int i,j,k,n,nx,ny ,counter;
  int N_Cells, N_Vertices, N_Comps, MaxN_VerticesPerCell;
  TBaseCell *cell;
  double *Data_sorted;
  double max_x, max_y, min_x, min_y, maxdist_x, maxdist_y;
  int *FESpaceNumber;
  TFESpace2D *fespace;
  double val[5], *X, *Y, *Val;
  int cordx, cordy;
  TFEFunction2D *u;
  #ifdef __3D__
  double z;
  #endif

  MaxN_VerticesPerCell = 4;      // 2D case

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(8);

  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];

  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  OutPut("WriteMatlab: " << N_ScalarVar << " scalar fields, "
	 << N_VectorVar << " vector fields; " <<
	 " WRITING ONLY FIRST SCALAR FIELD !!!"<< endl);
  //int N_Unknowns = fespace->GetN_DegreesOfFreedom();
  N_Cells = Coll->GetN_Cells();
  N_Vertices = Coll->GetCell(0)->GetN_Vertices();
  if (N_Vertices == 3)
  {
      // triangles
      i = N_Cells/2;
      nx = ny = (int) (sqrt(i)+1e-6);
  }
  else
  {
      // quads 
      nx = ny = (int) (sqrt(N_Cells)+1e-6);
  }
  OutPut("WriteMatlab: grid dimension " << nx + 1 << endl);
  X = new double[N_Cells*MaxN_VerticesPerCell];
  Y = new double[N_Cells*MaxN_VerticesPerCell];
  Val = new double[N_Cells*MaxN_VerticesPerCell];
  Data_sorted= new double[(nx+1)*(ny+1)];
  u = FEFunctionArray[0];

  counter = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Vertices = cell->GetN_Vertices();
    for (j=0; j< N_Vertices; j++)
    {
	cell->GetVertex(j)->GetCoords(X[counter],Y[counter]);
	u->FindGradientLocal(cell, i, X[counter], Y[counter], val);
	Val[counter] = val[0];
	counter++;
    }
  }  
 // for(i=0;i<counter;i++)
     // dat << setw(12) << X[i] << setw(12) << Y[i] << setw(12) << Val[i] << endl;

 //initialisation
 max_x =-1000;
 max_y =-1000;
 min_x = 1000;
 min_y = 1000; 

// find the size of the field
for(n=0;n<counter;n++)
   { 
    if(X[n]>max_x) max_x=X[n];
    if(Y[n]>max_y) max_y=Y[n];
    if(X[n]<min_x) min_x=X[n];
    if(Y[n]<min_y) min_y=Y[n];
   }

 maxdist_x = max_x - min_x;
 maxdist_y = max_y - min_y;


 
for(n=0;n<counter;n++)
  { //x-coordinate of the sorted array 
    cordx=(int) (((X[n]-min_x)*nx)/maxdist_x + 1e-6);
    //x-coordinate of the sorted array
    cordy=(int) (((Y[n]-min_y)*ny)/maxdist_y + 1e-6);
    //write the Data in the array
    Data_sorted[cordx+nx*cordy] = Val[n];
  }


// write the data into the file
//first line: vector x
for(i=0;i<=nx;i++)
   {
    dat << setw(12) << (min_x+i/(maxdist_x*nx));
   }
dat << endl;

//second line: vector y
for(j=0;j<=ny;j++)
   {
    dat << setw(12) << (min_y+j/(maxdist_y*ny)) ;
   }
dat << endl;

//sorted array, transposed for matlab
//dim nx+1 times nx+1  
for(j=0;j<=ny;j++)
    { 
    for(i=0;i<=nx;i++)
       {
        dat << setw(12) <<  Data_sorted[i+nx*j];
       }  
    dat << endl; 
    }

  dat.close();
  
  delete X;
  delete Y;
  delete Val;
  delete Data_sorted;
  delete FESpaceNumber;

  cout << endl;
  cout << "wrote output into matlab file: " << name << endl;
  
  return 0;
}

#endif

/** write stored data into a MATLAB file */
int TOutput2D::WriteMatlabOld(const char *name)
{
  int i,j,k,l,m,n;
  int N_Cells, N_Vertices, N_Comps, MaxN_VerticesPerCell;
  int N_BaseFunct, N_DOF;
  TBaseCell *cell;
  double x[4], y[4];
  double *Data;
  int *FESpaceNumber;
  TFESpace2D *fespace;
  char *Comment;
  double xi[4], eta[4];
  TBaseFunct2D *bf;
  FE2D FE_ID;
  double BFValues[4][MaxN_BaseFunctions2D];
  double *FEValues;
  int *GlobalNumbers, *BeginIndex, *Index;
  double LocValues[MaxN_BaseFunctions2D];
  double s;
  char str[15];
  #ifdef __3D__
  double z;
  #endif

  MaxN_VerticesPerCell = 4;      // 2D case

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
    fespace=FEFunctionArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[k]=j;
  }

  // one additional column for absolute values of velocity
  N_Comps++;

  Data = new double[MaxN_VerticesPerCell*N_Comps];
  //double *Data_aux = new double[MaxN_VerticesPerCell*N_Comps];
  //double *Data_sorted = new double[MaxN_VerticesPerCell*N_Comps];

  dat << "# Matlab data generated by MooNMD" << endl;
  dat << "# useful only for rectangular domains" << endl;
  dat << "#" << setw(11) << "x" << setw(12) << "y";
  for(i=0;i<N_ScalarVar;i++)
    dat << setw(12) << FEFunctionArray[i]->GetName();

  for(i=0;i<N_VectorVar;i++)
  {
    Comment = FEVectFunctArray[i]->GetName();
    for(k=0;k<FEVectFunctArray[i]->GetN_Components();k++)
      dat << setw(10) << Comment << "-" << setw(1) << k+1;
  }
  if(N_VectorVar>0)
  {
    str[0] = '|'; str[1] = '\0';
    Comment = FEVectFunctArray[0]->GetName();
    strncat(str, Comment, 10);
    strcat(str,"|");
    dat << setw(12) << str;
  }

  dat << endl;
  dat << "#" << endl;

  N_Cells = Coll->GetN_Cells();
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Vertices = cell->GetN_Vertices();
    switch(N_Vertices)
    {
      case 3:
        xi[0] = 0; xi[1] = 1; xi[2] = 0;
        eta[0] = 0; eta[1] = 0; eta[2] = 1;
        break;

      case 4:
        xi[0] = -1; xi[1] = 1; xi[2] = 1; xi[3] = -1;
        eta[0] = -1; eta[1] = -1; eta[2] = 1; eta[3] = 1;
        break;
    }                            // endswitch

    for(j=0;j<N_Vertices;j++)
    {
      #ifdef __3D__
      cell->GetVertex(j)->GetCoords(x[j],y[j], z);
      #else
      cell->GetVertex(j)->GetCoords(x[j],y[j]);
      #endif
    }

    //loop over scalar values
    for(j=0;j<N_ScalarVar;j++)
    {
      fespace = FEFunctionArray[j]->GetFESpace2D();
      FEValues = FEFunctionArray[j]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
      Index = GlobalNumbers+BeginIndex[i];
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
      N_BaseFunct = bf->GetDimension();
      for(k=0;k<N_Vertices;k++)
      {
        bf->GetDerivatives(D00, xi[k], eta[k], BFValues[k]);
        bf->ChangeBF(Coll, cell, BFValues[k]);
      }
      for(k=0;k<N_BaseFunct;k++)
        LocValues[k] = FEValues[Index[k]];

      for(k=0;k<N_Vertices;k++)
      {
        s = 0;
        for(l=0;l<N_BaseFunct;l++)
          s += LocValues[l]*BFValues[k][l];
        Data[j*MaxN_VerticesPerCell+k] = s;
        Data[j*MaxN_VerticesPerCell+k] = s;
      }                          // endfor k
    }                            // endfor j

    n = N_ScalarVar;
    // loop over vector valued variables
    for(j=0;j<N_VectorVar;j++)
    {
      fespace = FEVectFunctArray[j]->GetFESpace2D();
      N_DOF = fespace->GetN_DegreesOfFreedom();
      FEValues = FEVectFunctArray[j]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
      Index = GlobalNumbers+BeginIndex[i];
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
      N_BaseFunct = bf->GetDimension();
      for(k=0;k<N_Vertices;k++)
      {
        bf->GetDerivatives(D00, xi[k], eta[k], BFValues[k]);
        bf->ChangeBF(Coll, cell, BFValues[k]);
      }

      for(m=0;m<FEVectFunctArray[j]->GetN_Components();m++)
      {
        for(k=0;k<N_BaseFunct;k++)
          LocValues[k] = FEValues[m*N_DOF+Index[k]];

        for(k=0;k<N_Vertices;k++)
        {
          s = 0;
          for(l=0;l<N_BaseFunct;l++)
            s += LocValues[l]*BFValues[k][l];
          Data[n*MaxN_VerticesPerCell+k] = s;
        }                        // endfor k
        n++;
      }                          // endfor m
      if(j==0)
        for(k=0;k<N_Vertices;k++)
          Data[k+MaxN_VerticesPerCell*(N_Comps-1)] =
            sqrt(  Data[k+MaxN_VerticesPerCell*N_ScalarVar]
            * Data[k+MaxN_VerticesPerCell*N_ScalarVar]
            + Data[k+MaxN_VerticesPerCell*(N_ScalarVar+1)]
            * Data[k+MaxN_VerticesPerCell*(N_ScalarVar+1)]);
    }                            // endfor j

    //dat << "#" << N_Vertices << endl;
    /*
    for(j=0;j<N_Vertices;j++)
    {
      Datax[j] << endl;
    }
    */
    for(j=0;j<N_Vertices;j++)
    {
      dat << setw(12) << x[j] << setw(12) << y[j];
      for(k=0;k<N_Comps;k++)
        dat << setw(12) << Data[k*MaxN_VerticesPerCell+j];
      dat << endl;
    }                            // endfor j
    dat << endl << endl;         // two blank lines
  }                              // endfor i

  dat.close();

  delete Data;

  cout << endl;
  cout << "wrote output into matlab file: " << name << endl;

  return 0;
}


TOutput2D::~TOutput2D()
{
  if(Data) delete Data;
  
  delete [] FESpaceArray;
  delete [] FEFunctionArray;
  delete [] FEVectFunctArray;
  delete [] ParameterValues;
  delete [] ParameterDescription; 
}


/** write streamline function into a gunplot file */
int TOutput2D::WriteGNU_iso(const char *name, int scalar)
{
  int i, k, idx, nextx, nexty, LocVertex;
  int N_Cells, N_Verts, N_BaseFct;
  double startx, starty, boundx, boundy, X, Y;
  bool found = false;
  TBaseCell *CellX, *CellY;
  TJoint *JointX, *JointY;
  TFESpace2D *FESpace;
  double *FEValues, value, LocValues[4];
  double BFValues[4][MaxN_BaseFunctions2D];
  int *GlobalNumbers, *BeginIndex, *Index;
  TBaseFunct2D *BaseFct;
  FE2D FE_ID;
  double xi[] =
  {
    -1,  1,  1, -1
  }
  , eta[] =
  {
    -1, -1,  1,  1
  };
  #ifdef __3D__
  double z;
  #endif

  std::ofstream dat(name);
  if (!dat)
  {
    cerr << "Error in WriteGNU_psi: cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::scientific);
  dat << setprecision(4);

  // check streamline function
  cout << "write scalar function " << FEFunctionArray[scalar]->GetName()
    << " into " << name << " (for Gnuplot)" << endl;

  if (scalar >= N_ScalarVar)
  {
    cerr << "Error in WriteGNU_iso: wrong function" << endl;
    return -1;
  }

  FESpace = FEFunctionArray[scalar]->GetFESpace2D();
  FEValues = FEFunctionArray[scalar]->GetValues();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  BeginIndex = FESpace->GetBeginIndex();

  //put collection index into clipboard
  N_Cells = Coll->GetN_Cells();
  for (i=0;i<N_Cells;i++)
    Coll->GetCell(i)->SetClipBoard(i);

  // get lower left cell
  Domain->GetBoundBox(startx, starty, boundx, boundy);
  if (boundx < 0) startx += boundx;
  if (boundy < 0) starty += boundy;

  for (i=0;i<N_Cells;i++)
  {
    CellX = Coll->GetCell(i);
    N_Verts = CellX->GetN_Vertices();
    if (N_Verts != 4)
    {
      cerr << "Error in WriteGNU_psi: it only works for quadrangles"
        << endl;
      return -1;
    }

    for (nextx=0;nextx<4;nextx++)
    {
      #ifdef __3D__
      CellX->GetVertex(nextx)->GetCoords(X, Y, z);
      #else
      CellX->GetVertex(nextx)->GetCoords(X, Y);
      #endif
      if (ABS(X - startx) < 1e-10)
        if (ABS(Y - starty) < 1e-10)
      {
        found = true;
        break;
      }
    }

    if (found) break;
  }

  if (!found)
  {
    cerr << "Error in WriteGNU_psi: start cell not found" << endl;
    return -1;
  }

  dat << "# gnuplot data generated by MooNMD" << endl;
  dat << "#" << setw(11) << "x" << setw(12) << "y" << setw(12) <<
    "stream" << endl;

  nextx = ++nextx % 4;
  nexty = -1;
  do
  {
    if (nexty != -1)
    {
      if (JointX->GetType() == BoundaryEdge) break;

      CellX = JointX->GetNeighbour(CellX);

      for (nextx=0;nextx<4;nextx++)
        if (CellX->GetJoint(nextx) == JointX) break;

      nextx = (nextx + 2) % 4;
      nexty = -1;
    }

    CellY = CellX;

    LocVertex = (nextx + 3) % 4;
    idx = CellY->GetClipBoard();
    Index = GlobalNumbers + BeginIndex[idx];
    FE_ID = FESpace->GetFE2D(idx, CellY);
    BaseFct = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
    N_BaseFct = BaseFct->GetDimension();
    BaseFct->GetDerivatives(D00, xi[LocVertex], eta[LocVertex],
      BFValues[LocVertex]);
    BaseFct->ChangeBF(Coll, CellY, BFValues[LocVertex]);
    for(k=0;k<N_BaseFct;k++)
      LocValues[k] = FEValues[Index[k]];

    value = 0;
    for(k=0;k<N_BaseFct;k++)
      value += LocValues[k]*BFValues[LocVertex][k];

    #ifdef __3D__
    CellY->GetVertex(LocVertex)->GetCoords(X, Y, z);
    #else
    CellY->GetVertex(LocVertex)->GetCoords(X, Y);
    #endif
    dat << setw(12) << X << setw(12) << Y << setw(12) << value << endl;

    do
    {
      if (nexty != -1)
      {
        if (!(CellY = JointY->GetNeighbour(CellY))) break;

        for (nexty=0;nexty<4;nexty++)
          if (CellY->GetJoint(nexty) == JointY) break;

        nexty = (nexty + 2) % 4;
      }
      else
        nexty = (nextx + 1) % 4;

      LocVertex = (nexty + 1) % 4;
      idx = CellY->GetClipBoard();
      Index = GlobalNumbers + BeginIndex[idx];
      FE_ID = FESpace->GetFE2D(idx, CellY);
      BaseFct = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
      N_BaseFct = BaseFct->GetDimension();
      BaseFct->GetDerivatives(D00, xi[LocVertex], eta[LocVertex],
        BFValues[LocVertex]);
      BaseFct->ChangeBF(Coll, CellY, BFValues[LocVertex]);
      for(k=0;k<N_BaseFct;k++)
        LocValues[k] = FEValues[Index[k]];

      value = 0;
      for(k=0;k<N_BaseFct;k++)
        value += LocValues[k]*BFValues[LocVertex][k];

#ifdef __3D__
      CellY->GetVertex(LocVertex)->GetCoords(X, Y, z);
#else
      CellY->GetVertex(LocVertex)->GetCoords(X, Y);
#endif
      dat << setw(12) << X << setw(12) << Y << setw(12) << value << endl;

    } while ( (JointY = CellY->GetJoint(nexty) ) );

    dat << endl;

  } while ( (JointX = CellX->GetJoint(nextx)));

  CellY = CellX;

  LocVertex = nextx;
  idx = CellY->GetClipBoard();
  Index = GlobalNumbers + BeginIndex[idx];
  FE_ID = FESpace->GetFE2D(idx, CellY);
  BaseFct = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
  N_BaseFct = BaseFct->GetDimension();
  BaseFct->GetDerivatives(D00, xi[LocVertex], eta[LocVertex],
    BFValues[LocVertex]);
  BaseFct->ChangeBF(Coll, CellY, BFValues[LocVertex]);
  for(k=0;k<N_BaseFct;k++)
    LocValues[k] = FEValues[Index[k]];

  value = 0;
  for(k=0;k<N_BaseFct;k++)
    value += LocValues[k]*BFValues[LocVertex][k];

  #ifdef __3D__
  CellY->GetVertex(LocVertex)->GetCoords(X, Y, z);
  #else
  CellY->GetVertex(LocVertex)->GetCoords(X, Y);
  #endif
  dat << setw(12) << X << setw(12) << Y << setw(12) << value << endl;

  nexty = -1;

  do
  {
    if (nexty != -1)
    {
      if (!(CellY = JointY->GetNeighbour(CellY))) break;

      for (nexty=0;nexty<4;nexty++)
        if (CellY->GetJoint(nexty) == JointY) break;

      nexty = (nexty + 2) % 4;
    }
    else
      nexty = (nextx + 1) % 4;

    LocVertex = nexty;
    idx = CellY->GetClipBoard();
    Index = GlobalNumbers + BeginIndex[idx];
    FE_ID = FESpace->GetFE2D(idx, CellY);
    BaseFct = TFEDatabase2D::GetBaseFunct2DFromFE2D(FE_ID);
    N_BaseFct = BaseFct->GetDimension();
    BaseFct->GetDerivatives(D00, xi[LocVertex], eta[LocVertex],
      BFValues[LocVertex]);
    BaseFct->ChangeBF(Coll, CellY, BFValues[LocVertex]);
    for(k=0;k<N_BaseFct;k++)
      LocValues[k] = FEValues[Index[k]];

    value = 0;
    for(k=0;k<N_BaseFct;k++)
      value += LocValues[k]*BFValues[LocVertex][k];

    #ifdef __3D__
    CellY->GetVertex(LocVertex)->GetCoords(X, Y, z);
    #else
    CellY->GetVertex(LocVertex)->GetCoords(X, Y);
    #endif
    dat << setw(12) << X << setw(12) << Y << setw(12) << value << endl;

  } while ((JointY = CellY->GetJoint(nexty)) );

  dat.close();

  return 0;
}

#ifdef __2D__
/* WRITE STORED DATA INTO A GMV FILE **/

int TOutput2D::WriteGMV(const char *name)
{
  int i,j,k,l,m,n, N_;
  int N_Elements, N_Vertices, N_LocVertices;
  int N_LocDOF;
  int *VertexNumbers, *NumberVertex;
  int *IntArray;
  double *Coords, *DoubleArray;
  TBaseCell *cell;
  TVertex **Vertices, *Current, *Last;
  FE2D FE_ID;
  TBaseFunct2D *bf;
  TFESpace2D *fespace;
  int *GlobalNumbers, *BeginIndex, *DOF;
  double BFValues[MaxN_BaseFunctions2D], *Coeffs;
  double xi, eta, value;
  int N_Comp, Length;
  int Number;
  int MaxSubGridID;

  char matx[8];

  char gmvinput[] = "gmvinput";
  char ieeei4r8[] = "ieeei4r8";
  char nodes[]    = "nodes   ";
  char cells[]    = "cells   ";
  char quad[]	  = "quad    ";					// anstatt hexa
  char tria[]	  = "tria    ";					// anstatt tetra
  char velocity[] = "velocity";
  char variable[] = "variable";
  char endvars[]  = "endvars ";
  char endgmv[]   = "endgmv  ";
  char material[] = "material";

  static double QuadCoords[] =					// statische Variablen,
    { -1, -1,  1, -1,  1, 1,  -1, 1};				// Alle Richtungen im 2D Koordinatensystem
								// Macht es Sinn, den Namen von static Datein zu ändern?
  static double TriaCoords[] =					// Alternativ die dritte Komponente konstant 0 ?
      { 0, 0,  1, 0,  0, 1};

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
  }								// endfor i
  OutPut("N_LocVertices " << N_LocVertices << endl);
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
    }								// endfor k
  }								// endfor i
  Sort(Vertices, N_);						// sort the Vertices array

  Last = NULL;
  N_Vertices = 0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  Coords = new double[3*N_Vertices];				// 2 anstatt 3
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
	Vertices[i]->GetCoords(Coords[N_],Coords[N_+N_Vertices]);
	//OutPut(Coords[N_] << " " <<  Coords[N_+N_Vertices] << " " );
	Coords[N_+2*N_Vertices] = 0.0;
	k++;
	N_++;
	Last=Current;
    }
    NumberVertex[i]=k;
  }

#ifdef __COMPAQ__
  SwapDoubleArray(Coords, 2*N_Vertices);			// anstatt 3 (?)
#endif


  dat.write((char *)Coords,sizeof(double)*3*N_Vertices);	// 3 muss bleiben, sonst "invalid pointer" beim einlesen

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
    }								// endfor j
    switch(k)
    {
      case 3:
	dat.write(tria, 8);					// tria statt tetra
	Number = 3;


#ifdef __COMPAQ__
	 SwapIntArray(&Number, 1);
#endif


	dat.write((char *)(&Number), SizeOfInt);
      break;

      case 4:
	dat.write(quad, 8);					// quad statt hex
	Number = 4;


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
    SwapIntArray(&VertexNumbers[m-k], k);			// swap back (?)
#endif


  }								// endfor i
  MaxSubGridID++;
  dat.write(material, 8);
  Number = MaxSubGridID;
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char *)(&Number), SizeOfInt);
  Number = 0;							// cell data
#ifdef __COMPAQ__
  SwapIntArray(&Number, 1);
#endif
  dat.write((char *)(&Number), SizeOfInt);
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

  DoubleArray = new double[3*N_Vertices];			// 2 statt 3
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
      fespace = FEVectFunctArray[k]->GetFESpace2D();
      //      N_Comp = FEVectFunctArray[k]->GetN_Components();
      N_Comp =  3;
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
        FE_ID = fespace->GetFE2D(i, cell);
        bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
        DOF = GlobalNumbers+BeginIndex[i];
        N_LocDOF = bf->GetDimension();
        for(j=0;j<N_;j++)
        {
          switch(cell->GetType())
          {
            case Triangle: 
              xi = TriaCoords[2*j];
              eta = TriaCoords[2*j+1];
            break;
  
	    case Quadrangle:
              case Parallelogram:
	      case Rectangle:
              xi = QuadCoords[2*j];
              eta = QuadCoords[2*j+1];
            break;
            default:
              ErrMsg("unsupported cell type");
              exit(1);
          }
          bf->GetDerivatives(D00, xi, eta, BFValues);
          bf->ChangeBF(Coll, cell, BFValues);
          for(n=0;n<N_Comp-1;n++)
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

  k=0;

  dat.write(variable, 8);
  for(k=0;k<N_ScalarVar;k++)					//write scalar variables into file
  {
    fespace = FEFunctionArray[k]->GetFESpace2D();
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
      FE_ID = fespace->GetFE2D(i, cell);			// find FE data for this element
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
	switch(cell->GetType())
	{
	  case Triangle:					// statt Tetrahedron (wie oben)
	    xi = TriaCoords[2*j];				// wegen Veränderung der Deklaration!
	    eta = TriaCoords[2*j+1];				// aber nur wenn Namensänderung das Programm nicht zerstört.
	  break;
	    case  Quadrangle:
	      case Parallelogram:
	      case Rectangle:	  
	    xi = QuadCoords[2*j];
	    eta = QuadCoords[2*j+1];
	  break;
        default:
          ErrMsg("unsupported cell type");
          exit(1);
	}
        bf->GetDerivatives(D00, xi, eta, BFValues);
        bf->ChangeBF(Coll, cell, BFValues);
        value = 0;
        for(l=0;l<N_LocDOF;l++)
          value += BFValues[l] * Coeffs[DOF[l]];
        DoubleArray[VertexNumbers[m]-1] += value;
        IntArray[VertexNumbers[m]-1]++;
        m++;
      }								// endfor j
    }								// endfor i

    for(i=0;i<N_Vertices;i++)
      DoubleArray[i] /= IntArray[i];

#ifdef __COMPAQ__
    SwapDoubleArray(DoubleArray, N_Vertices);
#endif
    dat.write((char*)DoubleArray, N_Vertices*SizeOfDouble);
#ifdef __COMPAQ__
    SwapDoubleArray(DoubleArray, N_Vertices);			// swap back
#endif
  }								// endfor k
  dat.write(endvars, 8);

  delete NumberVertex;
  delete VertexNumbers;
  delete Vertices;
  delete DoubleArray;
  delete IntArray;
  delete Coords;

  dat.write(endgmv, 8);
  dat.close();

  if(TDatabase::ParamDB->SC_VERBOSE > 0)
    OutPut("wrote output into file: " << name << endl);

  return 0;
}

/** write stored PARALLEL data into a pvtu and vtu files (XML files for paraview) */

int TOutput2D::Write_ParVTK(
#ifdef _MPI
                                MPI_Comm comm,
#endif
                               int img, char *subID)
{
  int i, j, k,l,m,n, rank, size, N_, N_Elements, N_LocVertices;
  int N_Vertices, N_CellVertices, N_Comps;
  int *FESpaceNumber, N_LocDOF, Length, N_Comp, *GlobalNumbers, *BeginIndex, *DOF;
  int *VertexNumbers, *NumberVertex, begin, ID;

  double xi, eta, value, *Coeffs, *WArray, *DoubleArray;
//   double *BubbleArray;
  double BFValues[MaxN_BaseFunctions2D];
  double *Coords;
  double QuadCoords[] = { -1, -1, 1, -1, 1, 1, -1, 1};
  double TriaCoords[] = { 0, 0, 1, 0, 0, 1};

  char *VtkBaseName, Dquot;
  const char vtudir[] = "VTU";
  time_t rawtime;
  struct tm * timeinfo;

  TVertex **Vertices, *Last, *Current;
  TBaseCell *cell;
  TFESpace2D *fespace;
  TBaseFunct2D *bf;
  FE2D FE_ID;

  Dquot = 34; //  see ASCII Chart
  VtkBaseName = TDatabase::ParamDB->BASENAME;
  char *output_directory = TDatabase::ParamDB->OUTPUTDIR;
  // int AnsatzSpace = int(TDatabase::ParamDB->ANSATZ_ORDER);
#ifdef _MPI
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
#else
  rank = 0;
  size =1;
#endif


  if(!rank)
    mkdir(vtudir, 0777); // create the folder for each SubDomain vtu files

#ifdef _MPI
  MPI_Barrier (MPI_COMM_WORLD );
#endif

  std::ostringstream os;
  os << " ";

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

// write the master pvtu file
  if(rank==0)
   {
//    if( TDatabase::ParamDB->SC_VERBOSE > 0 )
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
    dat <<  "    Program: MooNMD " << endl;
    dat <<  "    Version: v1.0.0 " << endl;
    dat <<  "Date & Time: " <<asctime (timeinfo) << endl;
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

    dat <<   "<PCellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
    dat <<   "  <PDataArray type="<<Dquot<<"Int32"<<Dquot<<"   Name="<<Dquot<<"SubDomain"<<Dquot
        <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
    dat <<   "</PCellData>"<<endl;
    dat << endl;

//     dat <<   "<PCellData Scalars="<<Dquot<<"Bubble"<<Dquot<<">"<<endl;
//     dat <<   "  <PDataArray type="<<Dquot<<"Float32"<<Dquot<<"   Name="<<Dquot<<"Bubble"<<Dquot
//         <<"  NumberOfComponents="<<Dquot<<"1"<<Dquot<<"/>"<<endl;
//     dat <<   "</PCellData>"<<endl;
//     dat << endl;

    begin = 0;
    if(!TDatabase::ParamDB->Par_P1) // if root not take part in computation
      begin = 1;


    for(i=begin;i<size;i++)
     {
      if(TDatabase::ParamDB->Par_P1) // root take part in computation
       {
        ID = i+1;
       }
      else
       {
        ID = i;
       }

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
  } // if(!rank)

  ID = rank;

  if(rank==0 && !TDatabase::ParamDB->Par_P1)
    ID = -1;

 if(rank==ID)
  {
  // determine data for vtu file of each processor
  //int MaxN_VerticesPerCell = 4; // 2D case
  FESpaceNumber = new int[N_ScalarVar+N_VectorVar];
   //cout << "N_ScalarVar: " <<  N_ScalarVar << endl;
  N_Comps = 0;
  for(i=0;i<N_ScalarVar;i++)
  {
    N_Comps++;
    fespace=FEFunctionArray[i]->GetFESpace2D();
    j=0;
    while(FESpaceArray[j]!=fespace) j++;
    FESpaceNumber[i]=j;
  }

  k = N_ScalarVar;
  for(i=0;i<N_VectorVar;i++,k++)
  {
    N_Comps += FEVectFunctArray[i]->GetN_Components();
    fespace=FEVectFunctArray[i]->GetFESpace2D();
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
  Sort(Vertices, N_);

  Last=NULL;
  N_Vertices=0;
  for(i=0;i<N_LocVertices;i++)
    if((Current=Vertices[i])!=Last)
    {
      N_Vertices++;
      Last=Current;
    }
  Coords=new double[2*N_Vertices];
  VertexNumbers=new int[N_LocVertices];
  NumberVertex=new int[N_LocVertices];
  Last=NULL;
  N_=0; k=-1;
  for(i=0;i<N_LocVertices;i++)
  {
    if((Current=Vertices[i])!=Last)
    {
#ifdef __3D__
      double z;
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1], z);
#else
      Vertices[i]->GetCoords(Coords[N_],Coords[N_+1]);
//       if(rank==0) cout<< Coords[N_] << " " << Coords[N_]<<endl;
#endif
      k++;
      N_ += 2;
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
      m++;
    } // endfor j
  } //endfor i

  // additional column for absolute values of velocity
  N_Comps++;


//   cout << "MaxN_VerticesPerCell*N_Comps" << MaxN_VerticesPerCell*N_Comps << endl;
//   cout << "MaxN_VerticesPerCell" << MaxN_VerticesPerCell << endl;
//   cout << "N_Comps" << N_Comps << endl;

  // function values

  DoubleArray = new double[2*N_Vertices];
  WArray = new double[N_Vertices];
//   BubbleArray = new double[N_Elements];

    if(TDatabase::ParamDB->Par_P1) // root take part in computation
     {
      ID = rank+1;
     }
    else
     {
      ID = rank;
     }

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
       dat << "          " << Coords[N_]  << " " <<  Coords[N_+1]  << " ";
       dat << double(0) << endl;
       N_ +=2;
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
       case 4: dat << 9 << " ";
         break;
       case 3: dat << 5 << " ";
        break;
       }
     }
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</Cells>"<<endl;
    dat <<endl;
    dat <<  "<PointData Vectors="<<Dquot<<"Velocity"<<Dquot<<" "<<"Scalars="<<Dquot<<"Scalars"<<Dquot<<">"<<endl;
    dat << endl;
  // write vector variables into file


  for(k=0;k<N_VectorVar;k++)
  {
    fespace = FEVectFunctArray[k]->GetFESpace2D();
    N_Comp = FEVectFunctArray[k]->GetN_Components();
    Length = FEVectFunctArray[k]->GetLength();
    Coeffs = FEVectFunctArray[k]->GetValues();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();

   // cout << "N_Comp  " << N_Comp << endl;
    memset(DoubleArray, 0, SizeOfDouble*N_Vertices*N_Comp);
    memset(WArray, 0, SizeOfDouble*N_Vertices);
    m = 0;

    for(i=0;i<N_Elements;i++)
     {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
	    break;

	  case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
	    break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);

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
      	if(WArray[i]!=0.)
         DoubleArray[l] /= WArray[i];
	l++;
      }
    } // endfor l

//     for(i=0;i<2*N_Vertices;i++)
//     {
//       cout << "Do[" << i << "]" << DoubleArray[i] << endl;
//     }




//     dat << "SCALARS " << "|" << FEVectFunctArray[k]->GetName() << "|";
//     dat << " float"<< endl;
//     dat << "LOOKUP_TABLE " << "default" << endl;
//     l=0;
//     for(i=0;i<N_Vertices;i++)
//     {
//       t=0;
//       for(j=0;j<N_Comp;j++)
//       {
//        t+=DoubleArray[l]*DoubleArray[l];
//        l++;
//       }
//       dat << sqrt(t)<< endl;
//     }
//     dat << endl << endl;



      dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot
          <<FEVectFunctArray[k]->GetName()<<Dquot<<" NumberOfComponents="<<Dquot
          <<"3"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    l=0;
    for(i=0;i<N_Vertices;i++)
    {
      for(j=0;j<N_Comp;j++)
      {
	dat << DoubleArray[N_Comp*i+j] << " ";
      }
      dat << double(0) << " ";
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

  } // endfor k

  // write scalar variables into file
   for(k=0;k<N_ScalarVar;k++)
    {
      fespace = FEFunctionArray[k]->GetFESpace2D();
      Coeffs = FEFunctionArray[k]->GetValues();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();

      memset(DoubleArray, 0, SizeOfDouble*2*N_Vertices);
      memset(WArray, 0, SizeOfDouble*N_Vertices);
//       memset(BubbleArray, 0, SizeOfDouble*N_Elements);

      m = 0;

   for(i=0;i<N_Elements;i++)
    {
      cell = Coll->GetCell(i);
      N_ = cell->GetN_Vertices();

      // find FE data for this element
      FE_ID = fespace->GetFE2D(i, cell);
      bf = TFEDatabase2D::GetFE2D(FE_ID)->GetBaseFunct2D();
      DOF = GlobalNumbers+BeginIndex[i];
      N_LocDOF = bf->GetDimension();
//       switch(AnsatzSpace) // no bubble part
//        {
//          case 100:
//            N_LocDOF--; 
//          break; 
// 
//          case 201:
//            N_LocDOF -=3;
//          break;
//        }


      for(j=0;j<N_;j++)
      {
        switch(cell->GetN_Vertices())
        {
          case 3:
            xi = TriaCoords[2*j];
            eta = TriaCoords[2*j+1];
	    break;

	  case 4:
            xi = QuadCoords[2*j];
            eta = QuadCoords[2*j+1];
	    break;
        }
        bf->GetDerivatives(D00, xi, eta, BFValues);
        value = 0;

        for(l=0;l<N_LocDOF;l++)
          value += BFValues[l] * Coeffs[DOF[l]];
        DoubleArray[VertexNumbers[m]] += value;
        WArray[VertexNumbers[m]] +=1.;
        m++;
     } // endfor j

// // for bubble function
//         switch(cell->GetN_Vertices())
//          {
//           case 3:
//             xi = 1./3.;
//             eta = 1./3.;
// 	    break;
// 
// 	  case 4:
//             xi = 0.;
//             eta = 0.;
// 	    break;
//          }
// 
//       bf->GetDerivatives(D00, xi, eta, BFValues);
//       value = 0;
// 
//       switch(AnsatzSpace)
//        {
//          case 100:
//            for(l=3;l<4;l++)
//              value += BFValues[l] * Coeffs[DOF[l]];
//          break; 
// 
//          case 201:
//            for(l=6;l<9;l++)
//              value += BFValues[l] * Coeffs[DOF[l]];
//          break;
//        }
//       BubbleArray[i] = value;

    } // endfor i

    // non conforming
    for(i=0;i<N_Vertices;i++)
     {
     if(WArray[i]!=0)
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
     } // endfor k


    dat <<  "</PointData>"<<endl;
    dat << endl;

    dat <<  "<CellData Scalars="<<Dquot<<"SubDomain"<<Dquot<<">"<<endl;
    dat <<  "  <DataArray type="<<Dquot<<"Int32"<<Dquot<<" Name="<<Dquot<<"SubDomain"<<Dquot
        <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
    for(i=0;i<N_Elements;i++)
      dat << rank << " ";
    dat <<  "  </DataArray>"<<endl;
    dat <<  "</CellData>"<<endl;


//     dat <<  "<CellData Scalars="<<Dquot<<"Bubble"<<Dquot<<">"<<endl;
//     dat <<  "  <DataArray type="<<Dquot<<"Float32"<<Dquot<<" Name="<<Dquot<<"Bubble"<<Dquot
//         <<" NumberOfComponents="<<Dquot<<"1"<<Dquot<<" format="<<Dquot<<"ascii"<<Dquot<<">"<<endl;
//     for(i=0;i<N_Elements;i++)
//       dat << BubbleArray[i] << " ";
//     dat <<  "  </DataArray>"<<endl;
//     dat <<  "</CellData>"<<endl;

    dat <<  "</Piece>"<<endl;
    dat <<  "</UnstructuredGrid>"<<endl;
    dat <<  "</VTKFile>"<<endl;

    dat.close();

  delete [] FESpaceNumber;
  delete [] NumberVertex;
  delete [] VertexNumbers;
  delete [] Vertices;
  delete [] DoubleArray;
  delete [] WArray;
  delete [] Coords;

  }

//   img++;
  return 0;
}
#endif

/** write ascii file for tecplot **/
int TOutput2D::WriteAsciiPlt(const char *filename)
{
  TBaseCell *Cell;
  int N_Cells, N_;

  char tria[] = "ZONETYPE=FETRIANGLE";
  char quad [] = "ZONETYPE=FEQUADRILATERAL";
  char *type;

  cout << "Write plt file \"" << filename << "\" ...";
// exit(0);

  ComputeOutputData();

  N_Cells = Coll->GetN_Cells();

  // open ascii file
  std::ofstream dat(filename);
  if (!dat)
  {
    cerr << "cannot open file for output" << endl;
    return -1;
  }
  dat.setf(std::ios::fixed);
  dat << setprecision(4);

  // write header
  if (Data->Type == TOutputData::TRIA)
    type = tria;
  else
    type = quad;

  dat << "TITLE = \"plt file written by MooNMD\"" << endl;
  dat << "VARIABLES = \"X\", \"Y\"";
  for(int i=0;i<N_ScalarVar;i++)
    dat << ", \"" <<  FEFunctionArray[i]->GetName() << "\"";
  dat << endl;
  dat << "ZONE " << "DATAPACKING=POINT, ";
  dat << "NODES=" << Data->N_Nodes << ", ELEMENTS=" << N_Cells;
  dat << ", " << type << endl;

  // write points and point values
  for(int i=0;i<Data->N_Nodes;i++)
  {
    for(int j=0;j<N_ScalarVar+2;j++)
      dat << Data->FEFuncValues[j][i] << "\t";

    dat << endl;
  }

  // write elements
  for(int i=0;i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();
    
    for(int j=0;j<N_;j++)
    {
      dat << Cell->GetVertex(j)->GetClipBoard() + 1 << " ";
    }
    dat << endl;
  }

  // auxiliary data record
  dat << endl;
  for(int i=0;i<N_Parameters;i++)
  {
    dat << "DATASETAUXDATA ";
    dat << ParameterDescription[i] << " = ";
    dat << "\"" << ParameterValues << "\"" << endl;
  }

  dat.close();

  cout << " done" << endl;

  return 0;
}

/** write binary file for tecplot **/
int TOutput2D::WriteBinaryPlt(const char *filename)
{    
  int ret, N_Cells, N_Comps;
  int FileType = 0; // full
  int Debug = 0; // no debug
  int Precision = 1; // double precision
  int Zero = 0, One = 1;
  int PassiveVarList[2+N_ScalarVar+3*N_VectorVar];

  char Title[] = "tecplot file created by MooNMD";
  char ScratchDir[] = ".";
  const char *Variables;
  char ZoneTitle[] = "test";

  cout << "writing binary plt file \"" << filename << "\" ... ";

  std::ostringstream os;
  std::string string;

  N_Cells = Coll->GetN_Cells();

  memset(PassiveVarList, 0, sizeof(PassiveVarList));

  // variable names
  os << "X,Y";
  for(int i=0;i<N_ScalarVar;i++) // scalar
  {
    os << ",";
    os << FEFunctionArray[i]->GetName();
  }
  for(int i=0;i<N_VectorVar;i++) // vector
  {
    N_Comps = FEVectFunctArray[i]->GetN_Components();
    for(int j=0;j<N_Comps;j++)
    {
      os << ",";
      os << FEVectFunctArray[i]->GetName();
      os << "_";
      os << j;
    }
    os << ",|";
    os << FEVectFunctArray[i]->GetName();
    os << "|";
  }
  os << ends << std::flush;
  string = os.str();
  Variables = string.c_str();
//   cout << endl << Variables << endl;
//   exit(0);

  ret = TECINI112(Title, (char*) Variables, (char*) filename,
		  ScratchDir, &FileType, &Debug, &Precision);  
  if( ret == -1 ) 
  {
    cerr << "TECINI112 failed!" << endl;
    return -1;
  }

  ComputeOutputData();
  
  ret = TECZNE112(ZoneTitle, (int*) &Data->Type, &Data->N_Nodes, &N_Cells,
		  &Zero, &Zero, &Zero, &Zero,
		  &TDatabase::TimeDB->CURRENTTIME, &Zero, &Zero, &One, &Zero, &Zero, &Zero, &Zero, &Zero,
		  PassiveVarList, NULL, NULL, &Zero);
  if ( ret == -1 ) 
  {
    cerr << "TECZNE112 failed !" << endl;
    TECEND();
    return -1;
  }

  ret = TECDAT112(&Data->N_Data, Data->FEFuncValues[0], &One);
  if ( ret == -1 ) 
  {
    cerr << "TECDAT112 failed !" << endl;
    TECEND();
    return -1;
  }

  ret = TECNOD112(Data->ConList);
  if ( ret == -1 )
  {
    cerr << "TECNOD112 failed !" << endl;
    TECEND();
    return -1;
  }

  TECEND();

  cout << "done"  << endl;

  return 0;
}

void TOutput2D::ComputeOutputData()
{
  TBaseCell *Cell;
  TVertex   *Vert;
  int N_Cells, N_, counter;

  if (Data) delete Data;
  Data = new TOutputData();

  // reset clipboard of vertices
  N_Cells = Coll->GetN_Cells();
  for(int i=0;i<N_Cells;i++)
  {
    Cell = Coll->GetCell(i);
    N_ = Cell->GetN_Vertices();
    for(int j=0;j<N_;j++)
    {
      Cell->GetVertex(j)->SetClipBoard(-1);
    }
  }

  // count vertices
  counter = 0;
  for(int i=0;i<N_Cells;i++)
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

  // find all vertices and put into an array, create conlist
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

void TOutput2D::ComputeFEValues()
{
  TFESpace2D *fespace;
  TBaseFunct2D *BaseFunc;
  FE2D FeID;
  TBaseCell *Cell;
  TVertex *Vert;
  int N_VectorTotal, N_Comps, Length, index;
  int *GlobalNumbers, *BeginIndex, *DOF, N_Cells, N_, N_BF;
  int counter;
  double *FECoeffs;
  double BaseFuncValues[MaxN_BaseFunctions2D];
  
  double xi_tria[] = {0, 1, 0};
  double eta_tria[] = {0, 0, 1};
  double xi_quad[] = {-1, 1, 1,-1};
  double eta_quad[] = {-1,-1, 1, 1};
  double *xi, *eta, x, y;

  N_Cells = Coll->GetN_Cells();

  // count all vector components
  N_VectorTotal = 0;
  for(int i=0;i<N_VectorVar;i++)
  {
    N_VectorTotal += FEVectFunctArray[i]->GetN_Components()+1; // one for abs
  }

  // storage for FE values
  N_ = N_ScalarVar+2+N_VectorTotal;
  double *tmp = new double [Data->N_Nodes*N_];
  Data->FEFuncValues = new double* [N_];
  Data->N_Data = Data->N_Nodes*N_;

  for(int i=0;i<N_;i++)
    Data->FEFuncValues[i] = tmp + i*Data->N_Nodes;

  // find fe values (scalar)

  for(int i=0;i<N_ScalarVar;i++)
  {
      counter = 0;
      fespace = FEFunctionArray[i]->GetFESpace2D();
      GlobalNumbers = fespace->GetGlobalNumbers();
      BeginIndex = fespace->GetBeginIndex();
      FECoeffs = FEFunctionArray[i]->GetValues();

      for(int j=0;j<N_Cells;j++)
      {
	Cell = Coll->GetCell(j);
	FeID = fespace->GetFE2D(j, Cell);
	N_ = Cell->GetN_Vertices();
	DOF = GlobalNumbers + BeginIndex[j];
	BaseFunc = TFEDatabase2D::GetBaseFunct2DFromFE2D(FeID);
	N_BF = BaseFunc->GetDimension();

	switch (N_)
	{
	  case 3:
	    xi = xi_tria;
	    eta = eta_tria;
	    Data->Type = TOutputData::TRIA;
	    break;
      
	  case 4:
	    xi = xi_quad;
	    eta = eta_quad;
	    Data->Type = TOutputData::QUAD;
	    break;
	}

	for(int k=0;k<N_;k++)
	{
	  double val = 0;

	  Vert = Cell->GetVertex(k);
	  if ( Vert->GetClipBoard() == counter ) // skip redundant nodes
	  {
	    BaseFunc->GetDerivatives(D00, xi[k], eta[k], BaseFuncValues);

	    for(int l=0;l<N_BF;l++)
	    {
	      val += BaseFuncValues[l] * FECoeffs[DOF[l]];
	    }	  

	    Data->FEFuncValues[i+2][counter] = val;
	    counter++;
	  }
	}
      }     
  }

  // find fe values (vector)
  index = N_ScalarVar + 2;
  for(int i=0;i<N_VectorVar;i++)
  {
    counter = 0;
    fespace = FEVectFunctArray[i]->GetFESpace2D();
    GlobalNumbers = fespace->GetGlobalNumbers();
    BeginIndex = fespace->GetBeginIndex();
    N_Comps = FEVectFunctArray[i]->GetN_Components();
    Length = FEVectFunctArray[i]->GetLength();
    FECoeffs = FEVectFunctArray[i]->GetValues();

    for(int j=0;j<N_Cells;j++)
      {
	Cell = Coll->GetCell(j);
	FeID = fespace->GetFE2D(j, Cell);
	N_ = Cell->GetN_Vertices();
	DOF = GlobalNumbers + BeginIndex[j];
	BaseFunc = TFEDatabase2D::GetBaseFunct2DFromFE2D(FeID);
	N_BF = BaseFunc->GetDimension();

	switch (N_)
	{
	  case 3:
	    xi = xi_tria;
	    eta = eta_tria;
	    Data->Type = TOutputData::TRIA;
	    break;
      
	  case 4:
	    xi = xi_quad;
	    eta = eta_quad;
	    Data->Type = TOutputData::QUAD;
	    break;
	}

	for(int k=0;k<N_;k++)
	{
	  Vert = Cell->GetVertex(k);
	  if ( Vert->GetClipBoard() == counter ) // skip redundant nodes
	  {
	    double abs = 0;

	    BaseFunc->GetDerivatives(D00, xi[k], eta[k], BaseFuncValues);

	    for(int p=0;p<N_Comps;p++)
	    {
	      double val = 0;
	      for(int l=0;l<N_BF;l++)
	      {
		val += BaseFuncValues[l] * FECoeffs[p*Length+DOF[l]];
	      }	  
	      abs += val*val;
	      Data->FEFuncValues[index+p][counter] = val;	      
	    }
	    Data->FEFuncValues[index+N_Comps][counter] = sqrt(abs);
	    counter++;
	  }
	} // end for k
      } // end for j (cells) 

    index += N_Comps+1;
  } // end for i (FEVectFunct)

  // coordinates
  for(int i=0;i<Data->N_Nodes;i++)
  {
#ifdef __2D__
    Data->Nodes[i]->GetCoords(x,y);
#endif

    Data->FEFuncValues[0][i] = x;
    Data->FEFuncValues[1][i] = y;
  }
}

TOutput2D::TOutputData::~TOutputData()
{
  if (Nodes) delete [] Nodes;
  if (ConList) delete [] ConList;

  if (FEFuncValues) 
  {
    delete [] FEFuncValues[0];
    delete [] FEFuncValues;
  }
}
